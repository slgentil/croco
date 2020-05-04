"""
Module postp
 - Run class holds metadata about a CROCO run
"""

import os, fnmatch
from os.path import join, isfile, isdir
import numpy as np
from numpy import ma
import pandas as pd
import xarray as xr
from glob import glob

second2day = 1./86400.
grav = 9.81

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
                 dirname, prefix='',
                 outputs=[],
                 tdir_max=0,
                 grid_params={},
                 grid_periodicity=False, 
                 chunk_time=1,
                 verbose=False):
        """ Run object that gathers output and grid information

        Parameters
        ----------
        dirname: str
            Path to base directory where model output lies
        prefix: str, optional
            Prefix to all netcdf output files (e.g. 'file_')
        outputs: list, optional
            List of outputs to load, file names should look like: 
            [prefix+nc+'*.nc' for nc in outputs]
            Default is to load no outputs, i.e. outputs is empty
        tdir_max: int, optional
            Maximum run iteration loaded, default is 0
        grid_params: dict, optional
            Relevant grid parameters: y0, beta, yr_beta
        grid_periodicity: boolean, list, optional
            Passed to xgcm:
            Whether the grid is periodic (i.e. "wrap-around"). If a list is
            specified (e.g. ``['xi', 'eta']``), the axis names in the list will be
            be periodic and any other axes founds will be assumed non-periodic.
        chunk_time: int, optional
            Time chunk size, default is 1.
        verbose: boolean, optional
            Prints information if true. False by default
        """
        self.dirname = os.path.expanduser(dirname)
        self.verbose = verbose
        self.prefix = prefix
        if isfile(self.dirname+'t1/{}grid.nc'.format(prefix)):
            _nc_default = ['grid']
        else:
            # for backward compatibility
            _nc_default = ['his']
        self.outputs = list(set(_nc_default + outputs))
        self.zarr_dir = join(self.dirname, 'zarr')
        self.tdir_max = tdir_max # limits the number of t directories
        #
        self._grid_params = grid_params
        self.grid_periodicity = grid_periodicity
        #
        if isinstance(chunk_time,dict):
            self._chunk_time = chunk_time
        else:
            self._chunk_time = {nc:chunk_time for nc in self.outputs}
        #
        self._explore_tree()   # Find files that we know how to handle
        self._read_input_params()  # Scan croco.in for parameters
        self._read_output_stats()  # Scan output.mpi for parameters and stats
        self._open_data_files()   # Open the NetCDF files as scDatasets
        self._read_grid()    # Read the horizontal/vertical grid

    def __repr__(self):
        return ("Run: "+self.dirname+"\n"
                "  output keys: "+" / ".join(s for s in self.outputs)+"\n"
        )
        
    #def _repr_html_(self):
    #    return ("<b>"+self.dirname+"</b>"
    #            "<p>"+"".join(self.ds[s].__repr__() for s in self.outputs)+"</p>"
    #    )
        
    def __del__(self):
        """ Close any files linked to the datasets
        """
        for s in self.outputs:
            self.ds[s].close()

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

    def _explore_tree(self):
        if self.verbose:
            print("Analysing directory " + self.dirname)

        # Find the list of segments (t1, t2, ...), a.k.a. chains
        _segs = []
        _segs_avail = []
        for item in os.listdir(self.dirname):
            if isdir(join(self.dirname, item)):
                if item != "t0" and item[0] == "t" and item[1:].isdigit():
                    i = int(item[1:])
                    _segs_avail.append(i)
                    # Drop the 't' prefix so we can sort by increasing integer
                    if (self.tdir_max>0 and i<=self.tdir_max) or (self.tdir_max==0):
                        _segs.append(i)
        # Sort and restore the 't'
        self.segs_avail = ['t' + str(x) for x in sorted(_segs_avail)]
        self.segs = ['t' + str(x) for x in sorted(_segs)]
        if self.verbose:
            print("Found " + str(len(self.segs)) + " segments")

        # Now loop over segments in sequential order
        self.log_files = []
        self.nc_files = {s: [] for s in self.outputs}
        for segname in self.segs:
            if isfile(join(self.dirname, segname, "output.mpi")):
                self.log_files.append((segname, "output.mpi"))
            for key in self.outputs:
                # We use intermediate lists for each segment so we can sort them
                filename = []
                for cfile in os.listdir(join(self.dirname, segname)):
                    if fnmatch.fnmatchcase(cfile, self.prefix+"*.nc"):
                        if fnmatch.fnmatchcase(cfile, self.prefix+key+"*.nc"):
                            filename.append((segname, cfile))
                # Append sorted intermediate lists, such that the entire list is
                # in the correct order for accumulation with scDataset
                self.nc_files[key].extend(sorted(filename))
                
        if self.verbose:
            for key in self.outputs:
                print("Found " + str(len(self.nc_files[key])) + " " + key + " files")

    def _open_data_files(self):
        """ 
        Open data files (preferentially zarr archives, otherwise ncfiles)
        """
        # Open datasets found in list self.outputs
        if self.verbose:
            print("Opening datasets: ", self.outputs)
        self.ds = {}
        for key in self.outputs:
            # check exsistence of zarr archives
            zarr_archive = join(self.zarr_dir, key+'.zarr')
            if isdir(zarr_archive):
                _ds = xr.open_zarr(zarr_archive)
            else:
                _ds = self._create_xrDataset(key)
            #
            self.ds[key] = _ds

    def _create_xrDataset(self, key):
        # Helper function to synthesise inputs and call xarray
        assert len(self.nc_files[key]) > 0, \
                'No output with key {} to be found'.format(key)
        ncset = self.nc_files[key]
        tdir = [x[0] for x in ncset]
        offsets = self.t0.copy()
        files = [join(self.dirname, x[0], x[1]) for x in ncset]
        datasets = []
        for f, td in zip(files, tdir):
            if 'grid' in f:
                ds = xr.open_dataset(f, drop_variables=["x_rho","y_rho", \
                                                         "x_psi","y_psi"])
                ds = self._adjust_grid(ds)
                return ds
            else:
                try:
                    # chunks should be an option
                    ds = xr.open_dataset(f, chunks={'time_counter': self._chunk_time[key],
                                                    's_rho': 1})
                except ValueError:
                    ds = xr.open_dataset(f, chunks={'time_counter': self._chunk_time[key]})

                t0 = pd.Timestamp(ds.time_counter.time_origin).to_datetime64()
                _timevars = (t for t in ['time_counter', 'time_center', 'time_instant'] if t in ds)
                for tvar in _timevars:
                    t = ds[tvar]
                    t = ((t-t0)/np.timedelta64(1, 's') + offsets[td])*second2day
                    ds[tvar] = t
                # drop coordinates for easier concatenation
                #ds = ds.drop([k for k in ds.coords \
                #                if k not in ['time_counter','time_centered']])
                datasets.append(ds)
        ds = xr.concat(datasets, dim='time_counter',
                        coords='minimal', compat='override')
        ds = ds.rename_dims({'time_counter':'time'})
        ds = ds.set_coords([c for c in ds.data_vars if 'time' in c])
        ds = self._adjust_grid(ds)
        return ds

    def _read_input_params(self):
        """
        Short function to find parameters from a croco run directory.
        Currently we only examine croco.in.
        """
        # Read croco.in to extract parameters
        romsfile=join(self.dirname, self.segs[0], "croco.in")
        if isfile(romsfile):
            f = open(romsfile)
            pline=[] #previous line
            params = {}
            for line in iter(f):
                if 'time_stepping:' in pline:
                    params['dt']=tofloat(line.split()[1])
                    if self.verbose:
                        print("Detected time step of " + str(params['dt']) + " s")
                elif 'S-coord:' in pline:
                    tmp = [tofloat(x) for x in line.split()]
                    params['theta_s'], params['theta_b'] = tmp[0], tmp[1]
                    params['Hc'] = tmp[2]
                    if self.verbose:
                        print("Detected theta_s = " + str(params['theta_s']))
                        print("Detected theta_b = " + str(params['theta_b']))
                        print("Detected Hc = " + str(params['Hc']) + " m")
                elif 'rho0:' in pline:
                    params['rho0']=tofloat(line.split()[0])
                    if self.verbose:
                        print("Detected rho0 = " + str(params['rho0']) + " kg/m^3")
                pline=line
            f.close()
            self.params_input = params
        else:
            print("File not found: "+romsfile)
            self.params_input = None

    def _read_output_stats(self):
        # Now read each output.mpi and get the energy diagnostics
        self.t0 = dict()
        self.params_output = dict()
        n=0
        nbstats=None
        for ii, cfile in enumerate(self.log_files):
            f = open(join(self.dirname, cfile[0], cfile[1]))
            search = False
            firstline = True
            skipfirstline = ii>0
            pline=[] # previous line
            for line in iter(f):
                if ii==0 and 'hmax' in pline and 'grdmin' in pline:
                    self.H=tofloat(line.split()[1])
                    if self.verbose:
                        print("Detected H = " + str(self.H) + " m")
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
                    if self.verbose:
                        print("Found " + str(len(statnames)) + " columns in output.mpi:")
                        print(statnames)
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
        _dims = (d for d in ['x_v', 'y_u', 'x_w', 'y_w'] if d in ds.dims)
        for d in _dims:
            ds = ds.rename({d: d[0]+'_rho'})
        #
        _coords = [d for d in [d for d in ds.data_vars.keys()] if "nav_" in d]
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
        if 'grid' in self.outputs:
            ds = self.ds['grid']
        elif 'his' in self.outputs: # backward compatibility
            ds = self.ds['his']
        self.L = ds.sizes['x_rho']
        self.M = ds.sizes['y_rho']
        self.Lm = self.L - 1
        self.Mm = self.M - 1
        self.N  = ds.sizes['s_rho']
        self.Np = self.N + 1
        
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

        if self.verbose:
            print("Grid size: (L ,M, N) = (" + str(self.L) + ", " + str(self.M) + ", " + str(self.N) + ")")

        # Create xgcm grid
        # axis
        coords={'xi': {'center':'x_rho', 'inner':'x_u'}, 
                'eta': {'center':'y_rho', 'inner':'y_v'}, 
                's': {'center':'s_rho', 'outer':'s_w'}}
        # add metrics terms
        ds = _compute_metrics(ds)
        metrics = {
                   ('xi',): ['dx_rho', 'dx_u', 'dx_rho2d', 'dx_u2d'], # X distances
                   ('eta',): ['dy_rho', 'dy_v', 'dy_rho2d', 'dy_v2d'], # Y distances
                   ('xi', 'eta'): ['rA', 'rAu', 'rAv'] # Areas
                  }
        # generate xgcm grid
        self.xgrid = Grid(ds, 
                          periodic=self.grid_periodicity,
                          coords=coords, 
                          metrics=metrics)
        
        if 'grid' not in self.outputs: # backward compatibility
            self.grid = ds
    
    ### store data to zarr archives
    def _is_zarr_archive(self, key):
        """ Utils, test existence of a zarr archive
        """
        zarr_archive = join(self.zarr_dir, key+'.zarr')
        return isdir(zarr_archive)

    def store_zarr(self, chunks={}, **kwargs):
        """ writes data in zarr archives
        
        Parameters
        ----------
        chunks: dict, optional
            Dictionary with output keys as keys and dimension chunk sizes 
            as values
        **kwargs:
            Passed to xr.to_zarr method
        """
        assert len(self.segs)==len(self._segs_avail), \
            'Cannot store zarr archives if the full dataset is not loaded.' \
            +'Reload without limiting tdir_max.'
        D = {}
        for key in self.outputs: # filter out grid?
            ds = self[key]
            ds = _move_singletons_as_attrs(ds)
            #
            if key in chunks:
                ds = ds.chunk(chunks[key])
            # check each avariable average chunk sizes and stops if necessary
            n_threshold = 4000*3000
            # loop around vars, coords
            for k in ds: # .data_vars
                averaged_chunk_size = 1
                assert averaged_chunk_size > size_threshold, \
                    '{} chunks are two small, rechunk such that chunk sizes' \
                    +' exceed 4000x3000 elements on average'
            for k in ds.coords:
                averaged_chunk_size = None
                assert averaged_chunk_size > size_threshold, \
                    '{} chunks are two small, rechunk such that chunk sizes' \
                    +' exceed 4000x3000 elements on average'
            # datasets are first collected in a list
            D[key] = ds
        for key, ds in D.items():
            zarr_archive = join(self.zarr_dir, key+'.zarr')            
            ds.to_zarr(zarr_archive, **kwargs)
            print('{} stored in {}'.format(key,zarr_archive))
            
    def delete_nc(self, outputs=None):
        """ Delete netcdf files if zarr corresponding zarr archives have
        been produced
        
        Parameters
        ----------
        outputs: list of str, optional
            Output keys to delete, all keys available otherwise
        """
        if outputs is None:
            outputs = self.outputs
        for key in outputs:
            files = [join(self.dirname, x[0], x[1]) for x in self.nc_files[key]]
            for f in files:
                if _is_zarr_archive(key) and isfile(f):
                    os.remove(f)
                    print('{} deleted'.format(f))
    
    ### store/load diagnostics
    def store_diagnostic(self, name, data, 
                         overwrite=False,
                         file_format=None, 
                         directory='diagnostics/',
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
        **kwargs:
            Any keyword arguments that will be passed to the file writer
        """
        # create diagnostics dir if not present
        _dir = _check_diagnostic_directory(directory, self.dirname, 
                                           create=True)
        #
        if isinstance(data, xr.DataArray):
            self.store_diagnostic(name, data.to_dataset(),
                                  overwrite=overwrite,
                                  file_format=file_format,
                                  directory=directory,
                                  **kwargs
                                 )
        elif isinstance(data, xr.Dataset):
            success=False
            if file_format is None or file_format.lower() in ['zarr', '.zarr']:
                _file = join(_dir, name+'.zarr')
                write_kwargs = dict(kwargs)
                if overwrite:
                    write_kwargs.update({'mode': 'w'})
                _move_singletons_as_attrs(data).to_zarr(_file, **write_kwargs)
                success=True
            elif file_format.lower() in ['nc', 'netcdf']:
                _file = join(_dir, name+'.nc')
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
        _file = glob(join(_dir,name+'.*'))
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
        file_path = join(_dir, name+'.zarr')
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
        file_path = join(_dir, name+'.zarr')
        return load_vm(file_path, self.xgrid, persist=persist)

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
    return ds

def _check_file_overwrite(file, overwrite):
    """ Check whether one can overwrite a file, return False otherwise
    """
    _isfile = isfile(file)
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
    if isdir(directory):
        # directory is an absolute path
        _dir = directory
    elif isdir(join(dirname, directory)):
        # directory is relative
        _dir = join(dirname, directory)
    else:
        if create:
            # need to create the directory
            _dir = join(dirname, directory)
            os.mkdir(_dir)
            print('Create new diagnostic directory {}'.format(_dir))
        else:
            raise OSError('Directory does not exist')
    return _dir
    
def _move_singletons_as_attrs(ds):
    """ change singleton variables and coords to attrs
    This seems to be required for zarr archiving
    """
    for c in ds.coords:
        if ds[c].size==1:
            ds = ds.drop_vars(c).assign_attrs({c: ds[c].values})
    for v in ds.data_vars:
        if ds[v].size==1:
            ds = ds.drop_vars(v).assign_attrs({v: ds[v].values})
    return ds

