"""
Module postp
 - Run class holds metadata about a CROCO run
"""
#rom lpolib.grid import CGrid, s_coordinate
#rom lpolib.utils import tofloat, get_gradp, rho2u, rho2v
#rom lpolib.scDataset import scDataset
import os, fnmatch
import numpy as np
from numpy import ma
import pandas as pd
import xarray as xr
from glob import glob
#rom IPython import embed

second2day = 1./86400.

from crocosi.gridop import *

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
    def __init__(self, dirname, verbose=False, prefix='', open_nc=[],
                 tdir_max=0, grid_params={}, chunk_time=1):
        """
        Constructor; we inspect a run directory and assemble scDataset
                     classes for ave, his, etc., construct a CGrid class, and
                     read the energy diagnostics from output.mpi
        Parameters: dirname - directory to search for model output
                    verbose - if True, prints diagnostics as it works
                    open_nc - a list of nc files to open in addition to 'grid',
                              which may be an empty list.
        """
        self.dirname = os.path.expanduser(dirname)
        self.verbose = verbose
        self.prefix = prefix
        self.open_nc = ['grid'] + open_nc
        self.tdir_max = tdir_max # limits the number of t directories
        self._grid_params = grid_params
        if isinstance(chunk_time,dict):
            self._chunk_time = chunk_time
        else:
            self._chunk_time = {nc:chunk_time for nc in self.open_nc}
        #
        self._findfiles()   # Find files that we know how to handle
        self._readparams()  # Scan croco.in for parameters
        self._readstats()   # Scan output.mpi for parameters and stats
        self._openfiles()   # Open the NetCDF files as scDatasets
        self._readgrid()    # Read the horizontal/vertical grid

    def __repr__(self):
        return ("Run: "+self.dirname+"\n"
                "  datasets suffixes: "+" / ".join(s for s in self.open_nc)+"\n"
        )
        
    #def _repr_html_(self):
    #    return ("<b>"+self.dirname+"</b>"
    #            "<p>"+"".join(self.ds[s].__repr__() for s in self.open_nc)+"</p>"
    #    )
        
    def __del__(self):
        """ Close any files linked to the datasets
        """
        for s in self.open_nc:
            self.ds[s].close()

    def __getitem__(self, key):
        """ Load data set by providing suffix
        """
        # assert key in self.open_nc
        if key in self.open_nc:
            return self.ds[key]
        elif key in self.params_output.keys():
            return self.params_output[key]
        elif key in self.params_input.keys():
            return self.params_input[key]

    def __setitem__(self, key, item):
        """
        Load data set by providing suffix
        """
        if key in self.open_nc:
            self.ds[key] = item
        elif key in self.params_input.keys():
            self.params_input[key] = item
        elif key in self.params_output.keys():
            self.params_output[key] = item
    
    def __iter__(self):
        """ this allows looping over datasets
        """
        for item in self.open_nc:
            yield item, self.ds[item]

    def _findfiles(self):
        if self.verbose:
            print("Analysing directory " + self.dirname)

        # Find the list of segments (t1, t2, ...), a.k.a. chains
        self.segs = []
        for item in os.listdir(self.dirname):
            if os.path.isdir(os.path.join(self.dirname, item)):
                if item != "t0" and item[0] == "t" and item[1:].isdigit():
                    # Drop the 't' prefix so we can sort by increasing integer
                    if (self.tdir_max>0 and int(item[1:])<=self.tdir_max) or (self.tdir_max==0):
                        self.segs.append(int(item[1:]))
        # Sort and restore the 't'
        self.segs = ['t' + str(x) for x in sorted(self.segs)]
        if self.verbose:
            print("Found " + str(len(self.segs)) + " segments")

        # Now loop over segments in sequential order
        self.filename={}
        for suffix in self.open_nc:
            self.filename[suffix] = []
            self.flog = []
            # self.fave, self.fhis, self.finst, self.fsta1, self.fsta2, self.flog = [], [], [], [], [], []
            for segname in self.segs:
                if os.path.isfile(os.path.join(self.dirname, segname, "output.mpi")):
                    self.flog.extend([(segname, "output.mpi")])
                # We use intermediate lists for each segment so we can sort them
                filename = []
                for cfile in os.listdir(os.path.join(self.dirname, segname)):
                    if fnmatch.fnmatchcase(cfile, self.prefix+"*.nc"):
                        if fnmatch.fnmatchcase(cfile, self.prefix+suffix+"*.nc"):
                            filename.extend([(segname, cfile)])
                # Append sorted intermediate lists, such that the entire list is
                # in the correct order for accumulation with scDataset
                self.filename[suffix].extend(sorted(filename))
            if self.verbose:
                print("Found " + str(len(self.filename[suffix])) + " " + suffix + " files")

    def _openfiles(self):
        """
        Constructs xarray datasets for each list in self.open_nc
        """
        # Open datasets found in list self.open_nc
        if self.verbose:
            print("Opening NC datasets: ", self.open_nc)
        self.ds = {}
        for suffix in self.open_nc:
            if len(self.filename[suffix]) > 0 :
                self.ds[suffix] = self._create_xrDataset(self.filename[suffix], suffix)

    def _create_xrDataset(self, ncset, suffix):
        # Helper function to synthesise inputs and call xarray
        tdir = [x[0] for x in ncset]
        offsets = self.t0.copy()
        files = [os.path.join(self.dirname, x[0], x[1]) for x in ncset]
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
                    ds = xr.open_dataset(f, chunks={'time_counter': self._chunk_time[suffix],
                                                    's_rho': 1})
                except ValueError:
                    ds = xr.open_dataset(f, chunks={'time_counter': self._chunk_time[suffix]})

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

    def _readparams(self):
        """
        Short function to find parameters from a croco run directory.
        Currently we only examine croco.in.
        """
        # Read croco.in to extract parameters
        romsfile=os.path.join(self.dirname, self.segs[0], "croco.in")
        if os.path.isfile(romsfile):
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

    def _readstats(self):
        # Now read each output.mpi and get the energy diagnostics
        self.t0 = dict()
        self.params_output = dict()
        n=0
        nbstats=None
        for ii, cfile in enumerate(self.flog):
            f = open(os.path.join(self.dirname, cfile[0], cfile[1]))
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
                ds = ds.assign_coords(f=ds.f_rho) # shortcut?
        return ds

    def _readgrid(self, check=False):
        """ !!! to document !!!
        """
        from xgcm import Grid

        def s_coordinate(sc, theta_s, theta_b):
            '''
            Allows use of theta_b > 0 (July 2009)
            '''
            one64 = np.float64(1)

            if theta_s > 0.:
                csrf = ((one64 - np.cosh(theta_s * sc)) /
                        (np.cosh(theta_s) - one64))
            else:
                csrf = -sc ** 2
            sc1 = csrf + one64
            if theta_b > 0.:
                Cs = ((np.exp(theta_b * sc1) - one64) /
                      (np.exp(theta_b) - one64) - one64)
            else:
                Cs = csrf
            return Cs

        # Store grid sizes
        if 'grid' in self.open_nc:
            ds = self.ds['grid']
        elif 'his' in self.open_nc: 
            ds = self.ds['his']
        self.L = ds.sizes['x_rho']
        self.M = ds.sizes['y_rho']
        self.Lm = self.L - 1
        self.Mm = self.M - 1
        self.N  = ds.sizes['s_rho']
        self.Np = self.N + 1

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

        # Add topography in dataset if not in
        if 'h' not in list(ds.data_vars):
            ds['h']=(['y_rho','x_rho'],  self.H*np.ones((self.M,self.L)))

        if self.verbose:
            print("Grid size: (L ,M, N) = (" + str(self.L) + ", " + str(self.M) + ", " + str(self.N) + ")")

        # Create xgcm grid
        coords={'xi':{'center':'x_rho', 'inner':'x_u'}, 
                'eta':{'center':'y_rho', 'inner':'y_v'}, 
                's':{'center':'s_rho', 'outer':'s_w'}}
        ds.attrs['xgcm-Grid'] = Grid(ds, coords=coords)
        self.xgrid = ds.attrs['xgcm-Grid']
