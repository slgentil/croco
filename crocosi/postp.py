"""
Module postp
 - Run class holds metadata about a CROCO run
"""

import os, fnmatch
import os.path as path
import numpy as np
from numpy import ma
import pandas as pd
import xarray as xr
from pyamg import ruge_stuben_solver, solve       
from glob import glob

second2day = 1./86400.
grav = 9.81

# ideal dask size
_chunk_size_threshold = 4000*3000

import crocosi.gridop as gop

def tofloat(x):
    """
    Python doesn't recognise the fortran d/D for floating point numbers...
    """
    return float(x.replace('d','e').replace('D','e'))

class Run(object):
    """
    Class Run contains several xarray classes for each netCDF file type
    (ave, his, etc.), a grid object, and online output diagnostics (e.g. energy, ...).
    """
    def __init__(self,
                 dirname, 
                 prefix='',
                 suffix='',
                 outputs=[],
                 read_zarr=True,
                 tdir_max=0,
                 chunks={},
                 open_kwargs={},
                 persist=False,
                 grid_params={},
                 grid_periodicity=False,
                 grid_regular=True,
                 verbose=0):
        """ Run object that gathers output and grid information

        Parameters
        ----------
        dirname: str
            Path to base directory where model output lies
        prefix: str, optional
            Prefix to all netcdf output files (e.g. 'file_')
        outputs: list, str, optional
            List of outputs to load, file names should look like: 
            [prefix+nc+'*.nc' for nc in outputs]
            Default is to load no outputs, i.e. outputs is empty
            If 'all', all available outputs will be loaded
        read_zarr: boolean, optional
            Turn zarr file reading on/off. Default is True.
        tdir_max: int, optional
            Maximum run iteration loaded, default is 0
        chunks: dict, optional
            Chunks that will be either passed to file opener (zarr or netcdf) either
            prescribed in a rechunking operation
            must be a dict of dict with outputs type as parent dict keys and dimensions as child dict keys
        open_kwargs: dict, optional
            Keyword arguments passed to data opener (open_dataset, open_zarr)
        persist: boolean, optional
            Persists data into memory, you need to have enought memory to do that obviously
            Default if False
        grid_params: dict, optional
            Relevant grid parameters: y0, beta, yr_beta
        grid_periodicity: boolean, list, optional
            Passed to xgcm:
            Whether the grid is periodic (i.e. "wrap-around"). If a list is
            specified (e.g. ``['xi', 'eta']``), the axis names in the list will be
            be periodic and any other axes founds will be assumed non-periodic.
        verbose: int, optional
            Prints different levels of information. 0 (minimal info) by default
            3 levels (0,1,2) at the moment
        """
        self.dirname = os.path.expanduser(dirname)
        self.verbose = verbose
        self.prefix = prefix
        self.suffix = suffix
        self.zarr_dir = path.join(self.dirname, 'zarr')
        self.tdir_max = tdir_max # limits the number of t directories
        #
        self._grid_params = grid_params
        self.grid_periodicity = grid_periodicity
        self.grid_regular = grid_regular
        #
        self._explore_tree(outputs)   # Find files that we know how to handle
        self._read_input_params()  # Scan croco.in for parameters
        self._read_output_stats()  # Scan output.mpi for parameters and stats
        #
        self._sort_out_chunks(chunks)
        self._open_data_files(persist, read_zarr, **open_kwargs)   # Open the NetCDF files as scDatasets
        self._read_grid()    # Read the horizontal/vertical grid

    def __repr__(self):
        return ("Run: "+self.dirname+"\n"
                "  output keys: "+" / ".join(s for s in self.outputs)+"\n"
        )
        
    #def _repr_html_(self):
    #    return ("<b>"+self.dirname+"</b>"
    #            "<p>"+"".join(self.ds[s].__repr__() for s in self.outputs)+"</p>"
    #    )
        
    #def __del__(self):
    #    """ Close any files linked to the datasets
    #    """
    #    for s in self.outputs:
    #        self.ds[s].close()

    def __getitem__(self, key):
        """ Load data set by providing key
        """
        # grid is for backward compatibility, should be in self.outputs
        _gettable_attrs = ['grid', 'xgrid']
        # assert key in self.outputs
        if key in self.outputs:
            return self.ds[key]
        elif key in _gettable_attrs:
            return getattr(self, key)
        elif key in self.params_output.keys():
            return self.params_output[key]
        elif key in self.params_input.keys():
            return self.params_input[key]

    def __setitem__(self, key, item):
        """
        Load data set by providing key
        """
        if key in self.outputs:
            self.ds[key] = item
        elif key in self.params_input.keys():
            self.params_input[key] = item
        elif key in self.params_output.keys():
            self.params_output[key] = item
    
    def __iter__(self):
        """ this allows looping over datasets
        """
        for item in self.outputs:
            yield item, self.ds[item]

    def _explore_tree(self, outputs):
        if self.verbose>0:
            print("Analysing directory " + self.dirname)

        # Find the list of segments (t1, t2, ...), a.k.a. chains
        _segs = []
        _segs_avail = []
        for item in os.listdir(self.dirname):
            if path.isdir(path.join(self.dirname, item)):
                if item != "t0" and item[0] == "t" and item[1:].isdigit():
                    i = int(item[1:])
                    _segs_avail.append(i)
                    # Drop the 't' prefix so we can sort by increasing integer
                    if (self.tdir_max>0 and i<=self.tdir_max) or (self.tdir_max==0):
                        _segs.append(i)
        # Sort and restore the 't'
        self.segs_avail = ['t' + str(x) for x in sorted(_segs_avail)]
        self.segs = ['t' + str(x) for x in sorted(_segs)]
        if self.verbose>0:
            print("Found " + str(len(self.segs)) + " segments")
        #
        # find all avaible outputs
        self._find_available_outputs()
        if path.isfile(path.join(self.dirname,'t1/{}grid{}.nc'.format(self.prefix,self.suffix))):
            _outputs = ['{}grid{}'.format(self.prefix,self.suffix)]
        else:
            # for backward compatibility
            _outputs = ['his']
        self.grid_name = _outputs[0]
        if isinstance(outputs, str) and outputs=='all':
            _outputs += self.outputs_avail
        elif isinstance(outputs, list):
            _outputs += outputs
        # remove duplicated values while preserving order:
        self.outputs = list(dict.fromkeys(_outputs).keys())
        #
        # Now loop over segments in sequential order
        self.log_files = []
        self.nc_files = {s: [] for s in self.outputs}
        for segname in self.segs:
            if path.isfile(path.join(self.dirname, segname, "output.mpi")):
                self.log_files.append((segname, "output.mpi"))
            for key in self.outputs:
                # We use intermediate lists for each segment so we can sort them
                filename = []
                for cfile in os.listdir(path.join(self.dirname, segname)):
                    if fnmatch.fnmatchcase(cfile, self.prefix+"*.nc"):
                        if fnmatch.fnmatchcase(cfile, self.prefix+key+"*.nc"):
                            filename.append((segname, cfile))
                # Append sorted intermediate lists, such that the entire list is
                # in the correct order for accumulation with scDataset
                self.nc_files[key].extend(sorted(filename))
        #
        if self.verbose>1:
            for key in self.outputs:
                print("Found " + str(len(self.nc_files[key])) + " " + key + " files")

    def _find_available_outputs(self):
        """ Find all available outputs
        Note that outputs must not have an underscore (e.g. 'my_output...' not allowed)
        """
        outputs = set()
        for file_path in glob(path.join(self.dirname, self.segs[0], self.prefix+'*.nc')):
            file = file_path.split('/')[-1].split('.')[0]
            if 'rst' not in file and 'vmodes' not in file:
                outputs.add(file
                            .replace(self.prefix,'')
                            .split('_')[0]
                           )
        self.outputs_avail = list(outputs)
                                
    def _open_data_files(self, persist, read_zarr, **kwargs):
        """ 
        Open data files (preferentially zarr archives, otherwise ncfiles)
        """
        # Open datasets found in list self.outputs
        if self.verbose>0:
            print('Opening datasets: '+' / '.join(self.outputs))
        self.ds = {}
        for key in self.outputs:
            # check exsistence of zarr archives
            zarr_archive = path.join(self.zarr_dir, key+'.zarr')
            if path.isdir(zarr_archive) and read_zarr:
                ds = xr.open_zarr(zarr_archive, **kwargs)
                if 'time_counter' in ds.dims:
                    ds=ds.rename({'time_counter':'time'})
                # rechunk
                if self.chunks[key]:
                    ds = ds.chunk(self.chunks[key])
            else:
                ds = self._open_netcdf(key, **kwargs)
            #
            if persist:
                ds = ds.persist()
            #
            self.ds[key] = ds
            #
            if self.verbose>1:
                print("  {} - {:.1f} GB".format(key, ds.nbytes/1e9))

    def _open_netcdf(self, key, **kwargs):
        # Helper function to synthesise inputs and call xarray
        assert len(self.nc_files[key]) > 0, \
                'No output with key {} to be found'.format(key)
        ncset = self.nc_files[key]
        tdir = [x[0] for x in ncset]
        files = [path.join(self.dirname, x[0], x[1]) for x in ncset]
        data_files = []
        for f, td in zip(files, tdir):
            if 'grid' in f:
                _blacklist = ["x_rho", "y_rho", "x_psi","y_psi"]
                try:
                    _blacklist.extend(kwargs['drop_variables'])
                except Exception:
                    pass
                ds = xr.open_dataset(f, drop_variables= _blacklist)
                ds = self._adjust_grid(ds)
                return ds
            else:
                data_files.append(f)
        #
        def _preprocess_nc_file(ds):
            fname = ds.encoding['source']
            td = fname.split('/')[-2]
            if 'time_counter' in ds:
                t0 = pd.Timestamp(ds.time_counter.time_origin).to_datetime64()
                _timevars = (t for t in ['time_counter', 'time_center', 'time_instant'] 
                             if t in ds)
                for tvar in _timevars:
                    t = ds[tvar]
                    t = ((t-t0)/np.timedelta64(1, 's') + self.t0[td])*second2day
                    ds[tvar] = t
            return ds
        # load one file to figure out dimensions
        ds = xr.open_dataset(data_files[0])
        _chunks = {d:c for d, c in self.chunks[key].items() if d in ds.dims}
        #
        open_kwargs = {'concat_dim': 'time_counter',
                       'combine': 'nested',
                       'coords': 'minimal',
                       'parallel': False,
                       'compat': 'override'}
        open_kwargs.update(**kwargs)
        ds = xr.open_mfdataset(data_files,
                               preprocess=_preprocess_nc_file,
                               chunks=_chunks,
                               **open_kwargs
                              )
        #
        ds = ds.rename_dims({'time_counter':'time'})
        ds = ds.set_coords([c for c in ds.data_vars if 'time' in c])
        ds = self._adjust_grid(ds)
        return ds
    
    def _sort_out_chunks(self, chunks):
        """ sort out chunks for all outputs
        """
        if 'default' in chunks:
            _default_chunks = chunks['default']
        else:
            _default_chunks = {}
        # initiate all time chunks to default value
        self.chunks = {nc:_default_chunks for nc in self.outputs}
        # overwrite with prescribed values
        for key in chunks:
            if key in self.outputs and key != 'default':
                self.chunks[key] = chunks[key]
                
    def _read_input_params(self):
        """
        Short function to find parameters from a croco run directory.
        Currently we only examine croco.in.
        """
        # Read croco.in to extract parameters
        romsfile=path.join(self.dirname, self.segs[0], "croco.in")
        if self.verbose>0:
            print('Search for parameters in croco.in :')
        if path.isfile(romsfile):
            f = open(romsfile)
            pline=[] #previous line
            params = {}
            for line in iter(f):
                if 'time_stepping:' in pline:
                    params['dt']=tofloat(line.split()[1])
                elif 'S-coord:' in pline:
                    tmp = [tofloat(x) for x in line.split()]
                    params['theta_s'], params['theta_b'] = tmp[0], tmp[1]
                    params['Hc'] = tmp[2]
                elif 'rho0:' in pline:
                    params['rho0']=tofloat(line.split()[0])
                elif 'bottom_drag:' in pline:
                    tmp = [tofloat(x) for x in line.split()]
                    params['rdrg'] = tmp[0]
                    params['rdrg2'] = tmp[1]
                pline=line
            f.close()
            self.params_input = params
        else:
            print("File not found: "+romsfile)
            self.params_input = None
        #
        if self.params_input and self.verbose>1:
            for key, value in self.params_input.items():
                print('  {} = {}'.format(key, value))

    def _read_output_stats(self):
        # Now read each output.mpi and get the energy diagnostics
        self.t0 = dict()
        self.params_output = dict()
        n=0
        nbstats=None
        if self.verbose>0:
            print('Parameters detected in output.mpi :')
        for ii, cfile in enumerate(self.log_files):
            f = open(path.join(self.dirname, cfile[0], cfile[1]))
            search = False
            firstline = True
            skipfirstline = ii>0
            pline=[] # previous line
            for line in iter(f):
                if ii==0 and 'hmax' in pline and 'grdmin' in pline:
                    self.H=tofloat(line.split()[1])
                    if self.verbose>1:
                        print("  H = " + str(self.H) + " m")
                if 'MAIN' in line and 'started' in line:
                    search = True # Enable conversions
                if 'MAIN' in line and 'DONE' in line:
                    break # We're done with this file
                if 'OSI:' in line:
                    if ('(') in line: # parameter is a vector
                        varname = line.split()[1].strip()[:-2]
                        data=np.asfarray(line.split()[2:-1],float)
                        self.params_output[varname]=data
                    else:
                        varname = line.split()[1].strip()[:-1]
                        self.params_output[varname]=float(line.split()[2])
                if search and ii==0 and 'STEP' in line:
                    # Found header; save titles and create empty storage array
                    statnames = line.split()
                    nbstats = len(statnames)
                    statdata  = np.empty([5000,nbstats])
                    if self.verbose>1:
                        print("  Found " + str(len(statnames)) + " columns in output.mpi:")
                        for _s in statnames:
                            print('    {}'.format(_s))
                elif search and len(line.split())==nbstats and not 'STEP' in line:
                    if firstline:
                        # record the model starting offset
                        idx = [i for i, x in enumerate(statnames) if x=="time[DAYS]"][0]
                        tmpdata = [float(x) for x in line.split()]
                        self.t0[cfile[0]] = tmpdata[idx]*86400 # convert days to seconds
                        firstline=False
                    if skipfirstline:
                        # We skip the first line of data for the second file
                        # onwards because it is duplicated at the end of the
                        # previous file.
                        skipfirstline = False
                    else:
                        # Expand the array if needed
                        if n >= statdata.shape[0]:
                            statdata.resize((statdata.shape[0]*2, len(statnames)))
                        # Add these values to the arrays
                        try:
                            statdata[n,:] = [float(x) for x in line.split()]
                            n=n+1
                        except:
                            pass
                pline=line
            f.close()
        # Truncate unused end of the array
        statdata.resize((n, len(statnames)))
        self.stats = pd.DataFrame(statdata, columns=statnames).set_index('time[DAYS]')

    def _adjust_grid(self, ds):
        # relevant to regular/analytical grid for now
        #
        ds = ds.reset_coords([c for c in ds.coords if 'nav' in c])
        ds = ds.squeeze()
        # rename redundant dimensions
        _dims = (d for d in ['x', 'y', 'x_v', 'y_u', 'x_w', 'y_w'] if d in ds.dims)
        for d in _dims:
            ds = ds.rename({d: d[0]+'_rho'})
        #
        # rename variables nav_lon and nav_lat without any suffix
        _coords = (d for d in ['nav_lon', 'nav_lat'] if d in ds.data_vars.keys())
        for d in _coords:
            ds = ds.rename({d: d+'_rho'})
            
        _coords = [d for d in [d for d in ds.data_vars.keys()] if "nav_" in d]
        if self.grid_regular:
        # slice nav variables
            for v in _coords:
                xy = [d for d in ds[v].dims if d[0]=='y']
                if 'nav_lat' in v:
                    xy = [d for d in ds[v].dims if d[0]=='x']
                if len(xy)>0:
                    ds[v] = ds[v].isel({xy[0]: 0}).squeeze()
            # change nav variables in coordinates        
            ds = ds.set_coords(_coords)        
            # rename coordinates
            eta_suff={}
            for c in ds.coords:
                new_c = c.replace('nav_lat','y').replace('nav_lon','x')
                ds = ds.rename({c:new_c})
                # reset names and units
                ds[new_c] = (ds[new_c].assign_attrs(units='m', 
                                                    standard_name=new_c,
                                                    long_name=new_c)
                            )
        else:
            # change nav variables in coordinates        
            ds = ds.set_coords(_coords)        
            
        # fills in grid parameters, f, f0, beta
        _f0, _beta, _yrbeta = None, None, 0
        if 'f0' in self._grid_params:
            _f0 = self._grid_params['f0']
        elif 'f0' in self.params_output:
            _f0 = self.params_output['f0']
        #
        if 'beta' in self._grid_params:
            _beta = self._grid_params['beta']
        elif 'beta' in self.params_output:
            _beta = self.params_output['beta']
        if 'yrbeta' in self._grid_params:
            _yrbeta = self._grid_params['yrbeta']
        elif "yrbeta" in self.params_output:
            _yrbeta = self.params_output['yrbeta']
        #
        if "grid" in self.ds:
            ds = ds.assign_coords(f=self.ds['grid'].f) 
        elif _f0 and _beta:
            y_coords = [c for c in ds.coords if c[0]=='y']
            for c in y_coords:
                ds = ds.assign_coords(**{'f_'+c.split('_')[1]: \
                                            _beta*(ds[c]-_yrbeta)+_f0})
            ds = ds.assign_coords(f=ds.f_rho)
        return ds

    def _read_grid(self, check=False):
        """ !!! to document !!!
        """
        from xgcm import Grid

        def s_coordinate(sc, theta_s, theta_b):
            '''
            Allows use of theta_b > 0 (July 2009)
            '''
            if theta_s > 0.:
                csrf = ((1 - np.cosh(theta_s * sc)) /
                        (np.cosh(theta_s) - 1))
            else:
                csrf = -sc ** 2
            sc1 = csrf + 1
            if theta_b > 0.:
                Cs = ((np.exp(theta_b * sc1) - 1) /
                      (np.exp(theta_b) - 1) - 1)
            else:
                Cs = csrf
            return Cs

        # Store grid sizes
    
#         if 'grid' in self.outputs:
#             ds = self.ds['grid']
#         elif 'his' in self.outputs: # backward compatibility
#             ds = self.ds['his']
        ds = self.ds[self.grid_name]
        self.L = ds.sizes['x_rho']
        self.M = ds.sizes['y_rho']
        self.Lm = self.L - 1
        self.Mm = self.M - 1
        self.N  = ds.sizes['s_rho']
        self.Np = self.N + 1
        
        if self.grid_regular:
            # rename dimensions to be consistent with future xgrid object
            if 'x_psi' not in ds: # backward compatibility
                ds = ds.assign_coords(
                        x_psi = .5*(ds['x_rho'].shift(x_rho=-1)+ds['x_rho'])[:-1]
                            .drop('x_rho')
                            .rename({'x_rho': 'x_psi'}) )
            if 'y_psi' not in ds: # backward compatibility
                ds = ds.assign_coords(
                        y_psi = .5*(ds['y_rho'].shift(y_rho=-1)+ds['y_rho'])[:-1]
                            .drop('y_rho')
                            .rename({'y_rho': 'y_psi'}) )
            #
            if 'x_u' not in ds:
                ds = ds.rename({'x_psi': 'x_u'})
            if 'y_v' not in ds:
                ds = ds.rename({'y_psi': 'y_v'})

        # add S-coordinate stretching curves at RHO-points in dataset if not in
        if 'sc_r' not in list(ds.data_vars):
            sc = ((np.arange(1, self.N + 1, dtype=np.float64)) - self.N - 0.5) / self.N
            ds["sc_r"]=(['s_rho'],  sc)
        if 'sc_w' not in list(ds.data_vars):
            sc = (np.arange(self.N + 1, dtype=np.float64) - self.N) / self.N
            ds["sc_w"]=(['s_w'],  sc)
        if 'Cs_r' not in list(ds.data_vars):
            cs = s_coordinate(ds["sc_r"].values, self.params_input['theta_s'], self.params_input['theta_b'])
            ds["Cs_r"]=(['s_rho'],  cs)
        if 'Cs_w' not in list(ds.data_vars):
            cs = s_coordinate(ds["sc_w"].values, self.params_input['theta_s'], self.params_input['theta_b'])
            ds["Cs_w"]=(['s_w'],  cs)

        # Add topography in dataset if not in file
        if 'h' not in list(ds.data_vars):
            ds['h']=(['y_rho','x_rho'],  self.H*np.ones((self.M,self.L)))

        if self.verbose>0:
            print("Grid size: (L ,M, N) = (" + str(self.L) + ", " + str(self.M) + ", " + str(self.N) + ")")

        # Create xgcm grid
        # axis
        coords={'xi': {'center':'x_rho', 'inner':'x_u'}, 
                'eta': {'center':'y_rho', 'inner':'y_v'}, 
                's': {'center':'s_rho', 'outer':'s_w'}}
        # add metrics terms
        if self.grid_regular:
            ds, metrics = _compute_metrics(ds)
        else:
            ds, coords, metrics = self._compute_metrics_curvilinear()
        # generate xgcm grid
        self.xgrid = Grid(ds, 
                          periodic=self.grid_periodicity,
                          coords=coords, 
                          metrics=metrics,
                          boundary='extend')
        
        if 'grid' not in self.outputs: # backward compatibility
            self.grid = ds
        else:
            self.ds['grid'] = ds
    
    ### store data to zarr archives
    def _is_zarr_archive(self, key):
        """ Utils, test existence of a zarr archive
        """
        zarr_archive = path.join(self.zarr_dir, key+'.zarr')
        return path.isdir(zarr_archive)

    def store_zarr(self, 
                   chunks={}, 
                   auto_rechunk=True, 
                   tdir_max_overwrite=False,
                   return_ds=False,
                   **kwargs):
        """ writes data in zarr archives
        
        Parameters
        ----------
        chunks: dict, optional
            Dictionary with output keys as keys and dimension chunk sizes 
            as values
        auto_rechunk: boolean, optional
            Activate automatic rechunking which will ensure chunks are larger than
            _chunk_size_threshold (see postp.py). Default is True.
        tdir_max_overwrite: boolean, optional
            Allows creating zarr archives that do not consider the full run.
        return_ds: boolean, optional
            Skip writing to disk and return a dictionary of datasets, along with
            options for file writing
            (debug feature that should disapear)
        **kwargs:
            Passed to xr.to_zarr method
        """
        assert tdir_max_overwrite or len(self.segs)==len(self.segs_avail), \
            'Cannot store zarr archives if the full dataset is not loaded.' \
            +'Reload without limiting tdir_max.'
        D = {}
        for key in self.outputs:
            ds = self[key]
            ds = _move_singletons_as_attrs(ds)
            #
            if key in chunks:
                ds = ds.chunk(chunks[key])
            elif auto_rechunk:
                ds = _auto_rechunk(ds)
            # loop around vars, coords to check chunk sizes
            for k, da in ds.items(): # data_vars
                _check_chunks_sizes(da)
            for k, da in ds.coords.items():
                _check_chunks_sizes(da)
            # datasets are first collected in a list
            D[key] = ds
        D_out={}
        for key, ds in D.items():
            zarr_archive = path.join(self.zarr_dir, key+'.zarr')
            # need to get rid of several variables whose (time) units crashes to_zarr
            _to_delete = ['time_centered', 'time_centered_bounds', 
                          'time_counter_bounds', 'time_instant_bounds']
            for v in _to_delete:
                if v in ds:
                    del ds[v]
            # fix encoding inplace for nonchunked data
            _fix_nochunk_encoding(ds)
            #
            if return_ds:
                D_out[zarr_archive] = ds
            else:
                ds.to_zarr(zarr_archive, **kwargs)
                print('- {} stored'.format(key))
                print_zarr_archive_info(zarr_archive)
        if return_ds:
            return D_out, kwargs
            
    def delete_nc(self, outputs=None, test=True):
        """ Delete netcdf files if zarr corresponding zarr archives have
        been produced
        
        Parameters
        ----------
        outputs: list of str, optional
            Output keys to delete, all keys available otherwise
        """
        if outputs is None:
            outputs = self.outputs_avail
        for key in outputs:
            # check there is a corresponding zarr archive
            assert self._is_zarr_archive(key), \
                    '{} output has not been converted to zarr'.format(key)
            files = [path.join(self.dirname, x[0], x[1]) for x in self.nc_files[key]]
            for f in files:
                if path.isfile(f):
                    if not test:
                        os.remove(f)
                    print('{} deleted'.format(f))
    
    ### store/load diagnostics
    def store_diagnostic(self, name, data, 
                         overwrite=False,
                         file_format=None,
                         directory='diagnostics/',
                         auto_rechunk=True,
                         **kwargs
                        ):
        """ Write diagnostics to disk
        
        Parameters
        ----------
        name: str
            Name of a diagnostics to store on disk
        data: xr.Dataset, xr.DataArray (other should be implemented)
            Data to be stored
        overwrite: boolean, optional
            Overwrite an existing diagnostic. Default is False
        file_format: str, optional
            Storage file format (supported at the moment: zarr, netcdf)
        directory: str, optional
            Directory where diagnostics will be stored (absolute or relative to output directory).
            Default is 'diagnostics/'
        auto_rechunk: boolean, optional
            Automatically rechunk diagnostics such as o ensure they are not too small.
            Default is True.
        **kwargs:
            Any keyword arguments that will be passed to the file writer
        """
        # create diagnostics dir if not present
        _dir = _check_diagnostic_directory(directory, self.dirname, 
                                           create=True)
        if auto_rechunk:
            data = _auto_rechunk(data)
        #
        if isinstance(data, xr.DataArray):
            self.store_diagnostic(name, data.to_dataset(),
                                  overwrite=overwrite,
                                  file_format=file_format,
                                  directory=directory,
                                  auto_rechunk=True,                                  
                                  **kwargs
                                 )
        elif isinstance(data, xr.Dataset):
            success=False
            if file_format is None or file_format.lower() in ['zarr', '.zarr']:
                _file = path.join(_dir, name+'.zarr')
                write_kwargs = dict(kwargs)
                if overwrite:
                    write_kwargs.update({'mode': 'w'})
                data = _move_singletons_as_attrs(data)
                data = _reset_chunk_encoding(data)
                data.to_zarr(_file, **write_kwargs)
                success=True
            elif file_format.lower() in ['nc', 'netcdf']:
                _file = path.join(_dir, name+'.nc')
                write_kwargs = dict(kwargs)
                if overwrite:
                    write_kwargs.update({'mode': 'w'})
                data.to_netcdf(_file, **write_kwargs)
                success=True
            if success:
                print('data stored in {}'.format(_file))

    def load_diagnostic(self, name, 
                        directory='diagnostics/', 
                        **kwargs):
        """ Load diagnostics from disk
        
        Parameters
        ----------
        name: str, list
            Name of a diagnostics or list of names of diagnostics to load
        directory: str, optional
            Directory where diagnostics will be stored (absolute or relative to output directory).
            Default is 'diagnostics/'
        **kwargs:
            Any keyword arguments that will be passed to the file reader
        """
        _dir = _check_diagnostic_directory(directory, self.dirname)
        # find the diagnostic file
        _file = glob(path.join(_dir,name+'.*'))
        assert len(_file)==1, 'More that one diagnostic file {}'.format(_file)
        _file = _file[0]
        # get extension
        _extension = _file.split('.')[-1]
        if _extension=='zarr':
            return xr.open_zarr(_file, **kwargs)
        elif _extension=='nc':
            return xr.open_dataset(_file, **kwargs)
        else:
            raise NotImplementedError('{} extension not implemented yet'
                                      .format(_extension))
        
    ### wrappers (gop, xgcm grid)
    
    # horizontal grid moving
    def x2u(self, v):
        return gop.x2u(v, self.xgrid)
    def x2rho(self, v):
        return gop.x2rho(v, self.xgrid)
    def x2v(self, v):
        return gop.x2v(v, self.xgrid)
    def x2psi(self, v):
        return gop.x2psi(v, self.xgrid)
    # Vertical grid moving
    def x2w(self, v):
        return gop.x2w(v, self.xgrid)

    # xgcm functions
    def diff(self, *args, **kwargs):
        return self.xgrid.diff(*args, **kwargs)

    def interp(self, *args, **kwargs):
        return self.xgrid.interp(*args, **kwargs)

    def derivative(self, *args, **kwargs):
        return self.xgrid.derivative(*args, **kwargs)
    
    # vertical grid
    def get_z(self, *args, **kwargs):
        # should add mechanism to automatically load h from grid
        return gop.get_z(self, *args, **kwargs)
    
    # relative vorticity
    def get_relative_vorticity(self, u, v):
        xi = (-self.xgrid.derivative(u, 'eta') 
               + self.xgrid.derivative(v, 'xi')
              ).rename('vorticity')
        return xi
    
    def get_ertel_pv(self, ds, z=None, time=None, typ='ijk'):
        """
        #
        #   epv    - The ertel potential vorticity with respect to property 'lambda'
        #
        #                                       [ curl(u) + f ]
        #   -  epv is given by:           EPV = --------------- . del(lambda)
        #                                            rho
        #
        #   -  pvi,pvj,pvk - the x, y, and z components of the potential vorticity.
        #
        #   -  Ertel PV is calculated on horizontal rho-points, vertical w-points.
        #
        #
        #   tindex   - The time index at which to calculate the potential vorticity.
        #   depth    - depth
        #
        # Adapted from rob hetland.
        #
        """

        # Grid parameters
        ds = ds.isel(time=time) if time is not None else ds
        grid = self.xgrid

        f = self['grid'].f
        rho0 = self.params_input['rho0']

        # 3D variables
        z=z if z is not None else self.get_z(zeta=ds.ssh)
        dz = grid.diff(z,'s').persist()
        u = ds['u'].persist()
        v = ds['v'].persist()
        w = ds['w'].persist()

        try:
            rho = ds['rho'].persist()
        except Exception:
            print('rho not in history file')
            return

        if 'k' in typ:

            # Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*drho/dz
            # Compute d(v)/d(xi) at PSI-points.
            dvdxi = grid.derivative(v,'xi')

            # Compute d(u)/d(eta) at PSI-points.
            dudeta = grid.derivative(u,'eta')

            # Compute d(rho)/d(z) at horizontal RHO-points and vertical W-points
            drhodz = grid.diff(rho,'s') / dz

            #  Compute Ertel potential vorticity <k hat> at horizontal RHO-points and
            #  vertical W-points. 
            omega = dvdxi - dudeta
            omega = f + gop.x2rho(omega,grid)
            pvk = grid.interp(omega,'s') * drhodz
            del dvdxi, dudeta, drhodz, omega
        else:
            pvk = 0.

        if 'i' in typ:

            #  Ertel potential vorticity, term 2: (dw/dy - dv/dz)*(drho/dx)
            #  Compute d(w)/d(y) at horizontal V-points and vertical RHO-points
            dwdy = grid.derivative(w,'eta')

            #  Compute d(v)/d(z) at horizontal V-points and vertical W-points
            dz_v = grid.interp(dz,'eta')
            dvdz = grid.diff(v,'s') / dz_v

            #  Compute d(rho)/d(xi) at horizontal U-points and vertical RHO-points
            drhodx = grid.derivative(rho,'xi')

            #  Add in term 2 contribution to Ertel potential vorticity at horizontal RHO-points and
            #  vertical W-points.
            pvi = (gop.x2w(dwdy, grid)-gop.x2w(dvdz, grid)) * gop.x2w(drhodx, grid)
            del dwdy, dz_v, dvdz, drhodx
        else:
            pvi = 0.

        if 'j' in typ:

            #  Ertel potential vorticity, term 3: (du/dz - dw/dx)*(drho/dy)
            #  Compute d(u)/d(z) at horizontal U-points and vertical W-points
            dz_u = grid.interp(dz, 'xi')
            dudz = grid.diff(u,'s') / dz_u

            #  Compute d(w)/d(x) at horizontal U-points and vertical RHO-points
            dwdx = grid.derivative(w,'xi')

            #  Compute d(rho)/d(eta) at horizontal V-points and vertical RHO-points
            drhodeta = grid.derivative(rho,'eta')

            #  Add in term 3 contribution to Ertel potential vorticity at horizontal RHO-points and
            #  vertical W-points..
            pvj =  (gop.x2w(dudz,grid)-gop.x2w(dwdx,grid)) * gop.x2w(drhodeta,grid)
            del dz_u, dudz, dwdx, drhodeta

        else:
            pvj = 0.

        #
        #
        # Sum potential vorticity components, and divide by rho0
        #
        pvi = pvi / rho0
        pvj = pvj / rho0
        pvk = pvk / rho0
        pv = pvi + pvj + pvk

        z_w = self.get_z(vgrid='w').fillna(0.)
        pv = pv.assign_coords(coords={"z":z_w})

        return pv.squeeze()
    
    def get_dtdz(self, ds, z=None, time=None):
        """
        Compute dT/dz at horizontal rho point/vertical w point
        ds : dataset, containing T field
        z : xarray.DataArray, z in meters at rho points
        time : int, time index 
        """
        
        # keep time index if not None
        ds=ds.isel(time=time) if time is not None else ds

        # Grid parameters
        grid = self.xgrid

        # compute z coordinates
        z=z if z is not None else self.get_z(zeta=ds.ssh)
        z_w = self.get_z(zeta=ds.ssh, vgrid='w').fillna(0.)

        dtdz = (grid.diff(ds['temp'],'s') / grid.diff(z,'s')).squeeze()
        dtdz = dtdz.assign_coords(coords={"z":z_w})
        return dtdz
    
    def get_richardson(self, ds, z=None, time=None):
        """
        Ri is given by:      N²/((du/dz)² - (dv/dz)²)
             with N = sqrt(-g/rho0 * drho/dz)
        Ri is calculated at RHO-points and w level

        ds : dataset, which contains 3D u,v and rho fields
        z : xarray datarray, z depths at rho levels
        time : int, time index
        """

        # Grid parameters
        ds=ds.isel(time=time) if time is not None else ds
        grid = self.xgrid
        try:
            g = self.params_output['g']
        except :
            g = 9.81
        rho0 = self.params_input['rho0']

        # compute Z
        z=z if z is not None else self.get_z(zeta=ds.ssh)
        z_w = self.get_z(zeta=ds.ssh, vgrid='w').fillna(0.)

        N2 = self.get_N2(ds.rho, z, g=g)
        dudz = grid.diff(ds.u,'s')/grid.diff(gop.x2u(z,grid),'s')
        dvdz = grid.diff(ds.v,'s')/grid.diff(gop.x2v(z,grid),'s')

        Ri = xr.ufuncs.log10(N2 / (grid.interp(dudz,'xi')**2 +  grid.interp(dvdz,'eta')**2))
        Ri = Ri.assign_coords(coords={"z":z_w})
        return(Ri)
    
    def get_streamfunction(self,pm,pn,pv,verbo=False):
        """
        Compute the stream function from the relative vorticity
        Invert the laplacian to solve the poisson equation Ax=b
        A is the horizontal laplacian, b is the vorticity
        Input:
            - pm : (DataArray) 1/dx metric
            - pn : (DataArray) 1/dy metric
            - pv : (DataArray) relative vorticity
            - verbo : (Boolean) verbose mode
        Output:
            (DataArray) the computed streamfunction 
        """
        
        if np.any(np.isnan(pv)): 
            print("Can't inverse the laplacian, non compact domain, pv contains nan values")
            return None

        #######################################################
        #Create matrix A
        #######################################################
        if verbo: print('creating matrix A')
        A = gop.poisson_matrix(pm.values,pn.values)

        #######################################################
        #Solve matrix A
        A = A.tocsr()
        #######################################################

        if verbo: print('creating matrix b')
        b = -1. * pv.values.flatten() # right hand side
        ml = ruge_stuben_solver(A)                # construct the multigrid hierarchy
        if verbo: print(ml)                             # print hierarchy information
        x = ml.solve(b, tol=1e-8)                       # solve Ax=b to a tolerance of 1e-8     
        #x = solve(A,b,verb=False,tol=1e-8)

        if verbo: print("residual: ", np.linalg.norm(b-A*x))          # compute norm of residual vector

        chi = xr.DataArray(
            data=x.reshape(pm.shape),
            dims=["y_rho", "x_rho"],
            coords={'nav_lon_rho':pv.nav_lon_rho, 'nav_lat_rho':pv.nav_lat_rho}
            )
        return chi
    
    # buoyancy frequency
    def get_N2(self, *args, **kwargs):
        return gop.get_N2(self, *args, **kwargs)
    
    # pressure
    def get_p(self, *args, **kwargs):
        return gop.get_p(self.xgrid, *args, **kwargs)
    
    # vertical modes
    #   store:
    def store_vmodes(self, 
                     name, 
                     vmodes, 
                     projections=None, 
                     directory='diagnostics/',
                     **kwargs):
        """ store Vmode object and projections in a diagnostic directory
        
        Parameters
        ----------
        name: str
            Name of the vertical modes
        vmodes: Vmodes object
            Vertical mode object
        projections: xarray.Dataset
            Projections to load
        directory: str, optional
            Directory where diagnostics will be stored (absolute or relative to output directory).
            Default is 'diagnostics/'
        **kwargs: passed to to_zarr
        """
        _dir = _check_diagnostic_directory(directory, self.dirname, create=True)
        file_path = path.join(_dir, name+'.zarr')
        return vmodes.store(file_path, projections=projections, **kwargs)
        
    #   load:
    def load_vmodes(self, name, persist=None, directory='diagnostics/'):
        """ load Vmode object and projections from a diagnostic directory
        
        Parameters
        ----------
        name: str
            Name of the vertical modes
        persist: boolean
            Turns data loading on/off, default is False
        directory: str, optional
            Directory where diagnostics will be stored (absolute or relative to output directory).
            Default is 'diagnostics/'        
        """
        from .vmodes import load_vmodes as load_vm
        _dir = _check_diagnostic_directory(directory, self.dirname, create=False)
        file_path = path.join(_dir, name+'.zarr')
        return load_vm(file_path, self.xgrid, persist=persist)
  

    def _compute_metrics_curvilinear(self):
        from xgcm import Grid
    
        # curvilinear grid
        # Create xgcm grid without metrics
        ds = self.ds[self.grid_name]    
        coords={'xi': {'center':'x_rho', 'inner':'x_u'}, 
                'eta': {'center':'y_rho', 'inner':'y_v'}, 
                's': {'center':'s_rho', 'outer':'s_w'}}
        self.xgrid = Grid(ds, 
                  periodic=self.grid_periodicity,
                  coords=coords,
                  boundary='extend')
        
        # drop f from coordinates
        for s in self.outputs:
            try:
                self.ds[s] = self.ds[s].reset_coords(['f'])
            except:
                pass
            
        for s in self.outputs:
            # compute horizontal coordinates
            self.ds[s]['nav_lon_u'] = self.xgrid.interp(self.ds[s].nav_lon_rho,'xi')
            self.ds[s]['nav_lat_u'] = self.xgrid.interp(self.ds[s].nav_lat_rho,'xi')
            self.ds[s]['nav_lon_v'] = self.xgrid.interp(self.ds[s].nav_lon_rho,'eta')
            self.ds[s]['nav_lat_v'] = self.xgrid.interp(self.ds[s].nav_lat_rho,'eta')
            self.ds[s]['nav_lon_psi'] = self.xgrid.interp(self.ds[s].nav_lon_v,'xi')
            self.ds[s]['nav_lat_psi'] = self.xgrid.interp(self.ds[s].nav_lat_u,'eta')
            # set as coordinates in the dataset
            _coords = ['nav_lon_u','nav_lat_u',
                       'nav_lon_v','nav_lat_v',
                       'nav_lon_psi','nav_lat_psi',
                      ]
            self.ds[s] = self.ds[s].set_coords(_coords)
        
        for s in [s for s in self.outputs if s!=self.grid_name]:
            # compute z coordinate at rho/w points
            if 'ssh' in [v for v in self.ds[s].data_vars] and \
               's_rho' in [d for d in self.ds[s].dims.keys()] and \
                self.ds[s]['s_rho'].size>1:
                z_r = self.get_z(zeta=self.ds[s].ssh).fillna(0.)
                z_w = self.get_z(zeta=self.ds[s].ssh, vgrid='w').fillna(0.)
                self.ds[s]['z_r'] = z_r
                self.ds[s]['z_w'] = z_w
                self.ds[s]['z_u'] = self.xgrid.interp(z_r,'xi')
                self.ds[s]['z_v'] = self.xgrid.interp(z_r,'eta')
                self.ds[s]['z_psi'] = self.xgrid.interp(self.ds[s].z_u,'eta')
                self.ds[self.grid_name]['z_r'] = z_r
                self.ds[self.grid_name]['z_w'] = z_w
                self.ds[self.grid_name]['z_u'] = self.xgrid.interp(z_r,'xi')
                self.ds[self.grid_name]['z_v'] = self.xgrid.interp(z_r,'eta')
                self.ds[self.grid_name]['z_psi'] = \
                          self.xgrid.interp(self.ds[self.grid_name].z_u,'eta')
                # set as coordinates in the dataset
                _coords = ['z_r','z_w','z_u','z_v','z_psi']
                self.ds[s] = self.ds[s].set_coords(_coords)
                self.ds[self.grid_name] = self.ds[self.grid_name].set_coords(_coords)

        ds = self.ds[self.grid_name]    

        # add horizontal metrics for u, v and psi point
        if 'pm' in ds and 'pn' in ds:
            ds['dx_r'] = 1/ds['pm']
            ds['dy_r'] = 1/ds['pn']
        else: # backward compatibility, hack
            dlon = self.xgrid.interp(self.xgrid.diff(ds.nav_lon_rho,'xi'),'xi')
            dlat =  self.xgrid.interp(self.xgrid.diff(ds.nav_lat_rho,'eta'),'eta')
            ds['dx_r'], ds['dy_r'] = dll_dist(dlon, dlat, ds.nav_lon_rho, ds.nav_lat_rho)
        dlon = self.xgrid.interp(self.xgrid.diff(ds.nav_lon_u,'xi'),'xi')
        dlat = self.xgrid.interp(self.xgrid.diff(ds.nav_lat_u,'eta'),'eta')
        ds['dx_u'], ds['dy_u'] = dll_dist(dlon, dlat, ds.nav_lon_u, ds.nav_lat_u)
        dlon = self.xgrid.interp(self.xgrid.diff(ds.nav_lon_v,'xi'),'xi')
        dlat = self.xgrid.interp(self.xgrid.diff(ds.nav_lat_v,'eta'),'eta')
        ds['dx_v'], ds['dy_v'] = dll_dist(dlon, dlat, ds.nav_lon_v, ds.nav_lat_v)
        dlon = self.xgrid.interp(self.xgrid.diff(ds.nav_lon_psi,'xi'),'xi')
        dlat = self.xgrid.interp(self.xgrid.diff(ds.nav_lat_psi,'eta'),'eta')
        ds['dx_psi'], ds['dy_psi'] = dll_dist(dlon, dlat, ds.nav_lon_psi, ds.nav_lat_psi)

        # add vertical metrics for u, v, rho and psi points
        if 'z_r' in [v for v in ds.coords]:
            ds['dz_r'] = self.xgrid.diff(ds.z_r,'s')
            ds['dz_w'] = self.xgrid.diff(ds.z_w,'s')
            ds['dz_u'] = self.xgrid.diff(ds.z_u,'s')
            ds['dz_v'] = self.xgrid.diff(ds.z_v,'s')
            ds['dz_psi'] = self.xgrid.diff(ds.z_psi,'s')

        # add areas metrics for rho,u,v and psi points
        ds['rAr'] = ds.dx_psi * ds.dy_psi
        ds['rAu'] = ds.dx_v * ds.dy_v
        ds['rAv'] = ds.dx_u * ds.dy_u
        ds['rAf'] = ds.dx_r * ds.dy_r

        self.ds[self.grid_name] = ds
        
        # create new xgcmgrid with vertical metrics
        coords={'xi': {'center':'x_rho', 'inner':'x_u'}, 
                'eta': {'center':'y_rho', 'inner':'y_v'}}
        if 'z_r' in ds:
            coords.update({'s': {'center':'s_rho', 'outer':'s_w'}})
        metrics = {
               ('xi',): ['dx_r', 'dx_u', 'dx_v', 'dx_psi'], # X distances
               ('eta',): ['dy_r', 'dy_u', 'dy_v', 'dy_psi'], # Y distances
               ('xi', 'eta'): ['rAr', 'rAu', 'rAv', 'rAf'] # Areas
              }
        if 'z_r' in ds:
            metrics.update({('s',): ['dz_r', 'dz_u', 'dz_v', 'dz_psi', 'dz_w']}), # Z distances
        
        return ds, coords, metrics

def _compute_metrics(ds):
    """ Compute metrics from data available in grid.nc
    This code should be update for realistic curvilinear grids
    """
    if 'pm' in ds:
        ds['dx_rho'] = 1/ds['pm']
    else: # backward compatibility, hack
        ds['dx_rho'] = ds['x_rho'].shift(x_rho=-1)-ds['x_rho']
        ds['dx_rho'] = ds['dx_rho'].fillna(ds['dx_rho'].shift(x_rho=1))
    if 'pn' in ds:
        ds['dy_rho'] = 1/ds['pn']
    else: # backward compatibility, hack
        ds['dy_rho'] = ds['y_rho'].shift(y_rho=-1)-ds['y_rho']
        ds['dy_rho'] = ds['dy_rho'].fillna(ds['dy_rho'].shift(y_rho=1))
    ds['dx_rho2d'] = 0.*ds['dy_rho'] + ds['dx_rho']
    ds['dy_rho2d'] = ds['dy_rho'] + 0.*ds['dx_rho']     
    # code below should be updated for curvilinear grids
    ds['dx_u'] = ((ds['x_rho'].shift(x_rho=-1)-ds['x_rho'])[:-1]
                        .drop('x_rho')
                        .rename({'x_rho': 'x_u'})
                  )
    ds['dy_u'] = (ds['y_rho'].shift(y_rho=-1)-ds['y_rho'])
    ds['dx_v'] = (ds['x_rho'].shift(x_rho=-1)-ds['x_rho'])
    ds['dy_v'] = ((ds['y_rho'].shift(y_rho=-1)-ds['y_rho'])[:-1]
                     .drop('y_rho')
                     .rename({'y_rho': 'y_v'})
                  )
    ds['dx_u2d'] = 0.*ds['y_rho'] + ds['dx_u']
    ds['dy_u2d'] = ds['dy_u'] + 0.*ds['x_u']
    ds['dx_v2d'] = 0.*ds['y_v'] + ds['dx_v']
    ds['dy_v2d'] = ds['dy_v'] + 0.*ds['x_rho']
    ds['rA'] = ds['dx_rho']*ds['dy_rho']
    ds['rAu'] = ds['dx_u']*ds['dy_u']
    ds['rAv'] = ds['dx_v']*ds['dy_v'] 
    
    metrics = {
               ('xi',): ['dx_rho', 'dx_u', 'dx_rho2d', 'dx_u2d'], # X distances
               ('eta',): ['dy_rho', 'dy_v', 'dy_rho2d', 'dy_v2d'], # Y distances
               ('xi', 'eta'): ['rA', 'rAu', 'rAv'] # Areas
              }
    return ds, metrics
    
def dll_dist(dlon, dlat, lon, lat):
    """
    Converts lat/lon differentials into distances in meters
    PARAMETERS
    ----------
    dlon : xarray.DataArray longitude differentials 
    dlat : xarray.DataArray latitude differentials 
    lon : xarray.DataArray longitude values
    lat : xarray.DataArray latitude values
    RETURNS
    -------
    dx : xarray.DataArray distance inferred from dlon 
    dy : xarray.DataArray distance inferred from dlat 
    """
    distance_1deg_equator = 111000.0
    dx = dlon * xr.ufuncs.cos(xr.ufuncs.deg2rad(lat)) * distance_1deg_equator 
    dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
    return dx, dy
   
def _check_file_overwrite(file, overwrite):
    """ Check whether one can overwrite a file, return False otherwise
    """
    _isfile = path.isfile(file)
    if not _isfile or (_isfile and overwrite):
        return True
    else:
        print('Cannot overwrite {}'.format(file))
        return
    

def _check_diagnostic_directory(directory, dirname, 
                                create=False):
    """ Check existence of a directory and create it if necessary
    """
    # create diagnostics dir if not present
    if path.isdir(directory):
        # directory is an absolute path
        _dir = directory
    elif path.isdir(path.join(dirname, directory)):
        # directory is relative
        _dir = path.join(dirname, directory)
    else:
        if create:
            # need to create the directory
            _dir = path.join(dirname, directory)
            os.mkdir(_dir)
            print('Create new diagnostic directory {}'.format(_dir))
        else:
            raise OSError('Directory does not exist')
    return _dir
    
def _move_singletons_as_attrs(ds, ignore=[]):
    """ change singleton variables and coords to attrs
    This seems to be required for zarr archiving
    """
    for c,co in ds.coords.items():
        if co.size==1 and ( len(co.dims)==1 and co.dims[0] not in ignore or len(co.dims)==0 ):
            ds = ds.drop_vars(c).assign_attrs({c: ds[c].values})
    for v in ds.data_vars:
        if ds[v].size==1 and ( len(v.dims)==1 and v.dims[0] not in ignore or len(v.dims)==0 ):
            ds = ds.drop_vars(v).assign_attrs({v: ds[v].values})
    return ds

def _check_chunks_sizes(da):
    """ checks that chunk sizes are above the _chunk_size_threshold
    """
    averaged_chunk_size, total_size = _get_averaged_chunk_size(da)
    assert averaged_chunk_size==total_size or averaged_chunk_size>_chunk_size_threshold, \
        '{} chunks are two small, rechunk such that chunk sizes'.format(da.name) \
        + ' exceed {} elements on average,'.format(_chunk_size_threshold) \
        + ' there are currently ' \
        + '{} points per chunks on average'.format(averaged_chunk_size)

def _get_averaged_chunk_size(da):
    """ returns the averaged number of elements in the dataset
    """
    # total number of elements
    total_size = int(np.array(list(da.sizes.values())).prod())
    # averaged number of chunks along each dimension:
    if da.chunks:
        chunk_size_dims = np.array([np.max(d) for d in da.chunks])
        chunk_size = int(chunk_size_dims.prod())
    else:
        chunk_size = total_size
    return chunk_size, total_size

def _auto_rechunk_da(da):
    """ Automatically rechunk a DataArray such as chunks number of elements
    exceeds _chunk_size_threshold
    """
    dims = ['x_rho', 'x_u', 'y_rho', 'y_v', 
            's_rho', 's_w', 'mode', 
            'time']
    for d in dims:
        # gather da number of elements and chunk sizes
        averaged_chunk_size, total_size = _get_averaged_chunk_size(da)
        # exit if there is one chunk
        if averaged_chunk_size==total_size:
            break
        # rechunk along dimenion d
        if (d in da.dims) and averaged_chunk_size<_chunk_size_threshold:
            dim_chunk_size = np.max(da.chunks[da.get_axis_num(d)])
            # simple rule of 3
            factor = max(1, np.ceil(_chunk_size_threshold/averaged_chunk_size))
            new_chunk_size = int( dim_chunk_size * factor )
            # bounded by dimension size
            new_chunk_size = min(da[d].size, new_chunk_size)
            da = da.chunk({d: new_chunk_size})
    return da

def _auto_rechunk(ds):
    """ Wrapper around _auto_rechunk_da for datasets
    Accepts DataArrays as well however.
    """
    if isinstance(ds, xr.DataArray):
        return _auto_rechunk_da(ds)
    for k, da in ds.items(): # data_vars
        ds = ds.assign(**{k: _auto_rechunk_da(da)})
    for k, da in ds.coords.items():
        ds = ds.assign_coords(**{k: _auto_rechunk_da(da)})
    return ds

def _get_dir_size(dir_path):
    ''' Returns the size of a directory in bytes
    '''
    process = os.popen('du -s '+dir_path)
    size = int(process.read().split()[0]) # du returns kb
    process.close()
    return size*1e3

def print_zarr_archive_info(zarr_archive):
    """ Print basic information about a zarr archive
    """
    print('  location: {} '.format(zarr_archive))
    # get archive size
    arch_size = _get_dir_size(zarr_archive)
    print('  size:     {:.1f} GB'.format(arch_size/1e9))
    #
    ds = xr.open_zarr(zarr_archive)
    # get largest item typical chunks
    n_dim_max = 0
    for v in ds:
        if ds[v].chunks and ds[v].ndim>n_dim_max:
            _size = list(ds[v].sizes.values())
            _chunks = [np.max(d) for d in ds[v].chunks]
            n_dim_max = ds[v].ndim
    if n_dim_max>0:
        print('  typical chunks: ('
              +','.join('{}'.format(c) for c in _chunks)
              +') for size ('
              +','.join('{}'.format(c) for c in _size)
              +')'
             )
    else:
        print('  data is not chunked')
        
def _fix_nochunk_encoding(da):
    ''' Fix in place the encoding for nonchunked arrays such that zarr 
    writing does not automatically chunked arrays.
    
    Parameters
    ----------
    da: xr.DataArray, xr.Dataset
        variable or dataset to fix
    '''
    if isinstance(da, xr.Dataset):
        for v in da:
            _fix_nochunk_encoding(da[v])
    if isinstance(da, xr.DataArray):
        if not da.chunks:
            da.encoding['chunks'] = -1
            
def _reset_chunk_encoding(ds):
    ''' Delete chunks from variables encoding. 
    This may be required when loading zarr data and rewriting it with different chunks
    
    Parameters
    ----------
    ds: xr.DataArray, xr.Dataset
        Input data
    '''
    if isinstance(ds, xr.DataArray):
        return _reset_chunk_encoding(ds.to_dataset()).to_array()
    #
    for v in ds.coords:
        if 'chunks' in ds[v].encoding:
            del ds[v].encoding['chunks']
    for v in ds:
        if 'chunks' in ds[v].encoding:
            del ds[v].encoding['chunks']
    return ds
