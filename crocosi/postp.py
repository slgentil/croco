"""
Module postp
 - CROCOrun class holds metadata about a CROCO run
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

class CROCOrun(object):
    """
    Class CROCOrun contains several xarray classes for each netCDF file type
    (ave, his, etc.), a grid object, and online output diagnostics (e.g. energy, ...).
    """
    def __init__(self, dirname, verbose=False, prefix=None, open_nc=[],
                 tdir_max=0, grid_params={}):
        """
        Constructor; we inspect a run directory and assemble scDataset
                     classes for ave, his, etc., construct a CGrid class, and
                     read the energy diagnostics from output.mpi
        Parameters: dirname - directory to search for model output
                    verbose - if True, prints diagnostics as it works
                    open_nc - a list of nc files to open in addition to 'his',
                              which may be an empty list.
        """
        self.dirname = os.path.expanduser(dirname)
        self.verbose = verbose
        self.prefix = prefix
        self.open_nc = ['his'] + open_nc   # Ensure that we always open 'his'
        self.tdir_max = tdir_max # limits the number of t directories
        self._grid_params = grid_params
        #
        self._findfiles()   # Find files that we know how to handle
        self._readparams()  # Scan croco.in for parameters
        self._readstats()   # Scan output.mpi for parameters and stats
        self._openfiles()   # Open the NetCDF files as scDatasets
        self._readgrid()    # Read the horizontal/vertical grid

    def __del__(self):
        """
        Destructor: close the netcdf files upon destruction
        """
        # self.ds[:].close()

    def __getitem__(self,item):
        """
        Load data set by providing suffix
        """
        assert item in self.open_nc
        return self.ds[item]

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
            self.ds[suffix] = self._create_xrDataset(self.filename[suffix])

    def _create_xrDataset(self, ncset):
        # Helper function to synthesise inputs and call xarray
        tdir = [x[0] for x in ncset]
        offsets = self.t0.copy()
        files = [os.path.join(self.dirname, x[0], x[1]) for x in ncset]
        datasets = []
        deltat=0
        for f, td in zip(files, tdir):
            try:
                ds = xr.open_dataset(f, chunks={'time_counter': 1, 's_rho': 1})
            except ValueError:
                ds = xr.open_dataset(f, chunks={'time_counter': 1})
            if deltat==0:
                dt = offsets[td]
                deltat = (ds.time_counter[-1] - ds.time_counter[0]) \
                         /pd.Timedelta('1s')/(ds.time_counter.size-1) \
                         *ds.time_counter.size
            if 'time_counter' in ds:
                t = ds['time_counter']
                t0 = t.isel(time_counter=0)
                t = ((t-t0)/np.timedelta64(1, 's') + dt)*second2day
                ds['time_counter'] = t
            if 'time_center' in ds:
                t = ds['time_center']
                t = ((t-t0)/np.timedelta64(1, 's') + dt)*second2day
                ds['time_center'] = t
            dt = dt + deltat
            # drop coordinates for easier concatenation
            #ds = ds.drop([k for k in ds.coords \
            #                if k not in ['time_counter','time_centered']])
            datasets.append(ds)
        _ds = xr.concat(datasets, dim='time_counter',
                        coords='minimal', compat='override')
        _ds = _ds.rename({'time_counter':'time'})
        _ds = self._adjust_grid(_ds)
        return _ds

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
                elif 'tidal_diag:' in pline:
                    params['omega']=tofloat(line.split()[0])
                    if self.verbose:
                        print("Detected omega = " + str(params['omega']) + " 1/s")
                elif 'jet_ywidth' in pline:
                    params['jet_ywidth'] = tofloat(line.split()[2])
                    params['jet_weight'] = tofloat(line.split()[3])
                    if self.verbose:
                        print("Detected jet_ywidth = " + str(params['jet_ywidth']) + " m")
                        print("Detected jet_weight = " + str(params['jet_weight']))
                elif 'y_itide' in pline:
                    params['y_itide'] = tofloat(line.split()[4])
                    if self.verbose:
                        print("Detected y_itide = " + str(params['y_itide']) + " m")
                pline=line
            f.close()
            self.params = params
        else:
            print("File not found: "+romsfile)
            self.params = None

    def _readstats(self):
        # Now read each output.mpi and get the energy diagnostics
        self.t0 = dict()
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
                if search and ii==0 and 'STEP' in line:
                    # Found header; save titles and create empty storage array
                    statnames = line.split()
                    nbstats = len(statnames)
                    statdata  = np.empty([5000,nbstats])
                    if self.verbose:
                        print("Found " + str(len(statnames)) + " columns in output.mpi:")
                        print(statnames)
                elif search and len(line.split())==nbstats and not 'STEP' in line:
                    # The 114 condition may be fragile. ToDo: better condition
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
        eta_suff={}
        for c in ds.coords:
            new_c = c.replace('nav_lat','eta').replace('nav_lon','xi')
            ds = ds.rename({c:new_c})
            if 'eta_rho' in new_c:
                eta_suff[new_c] = new_c.lstrip('eta_rho')
        # fills in grid parameters, f, f0, beta
        if 'f0' in self._grid_params:
            ds = ds.assign_attrs(f0=self._grid_params['f0'])
        if 'beta' in self._grid_params:
            ds = ds.assign_attrs(beta=self._grid_params['beta'])
            for c, suff in eta_suff.items():
                ds = ds.assign_coords(**{'f'+suff: ds.beta*ds[c]+ds.f0})
        return ds

    def _readgrid(self, check=False):
        """ !!! old code, update or delete !!!
        """
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
        self.L = self.ds['his'].sizes['x_rho']
        self.M = self.ds['his'].sizes['y_rho']
        self.Lm = self.L - 1
        self.Mm = self.M - 1
        self.N  = self.ds['his'].sizes['s_rho']
        self.Np = self.N + 1

        # add S-coordinate stretching curves at RHO-points in dataset if not in
        if 'sc_r' not in list(self.ds['his'].data_vars):
            sc = ((np.arange(1, self.N + 1, dtype=np.float64)) - self.N - 0.5) / self.N
            self.ds['his']["sc_r"]=(['s_rho'],  sc)
        if 'Cs_r' not in list(self.ds['his'].data_vars):
            cs = s_coordinate(sc, self.params['theta_s'], self.params['theta_b'])
            self.ds['his']["Cs_r"]=(['s_rho'],  cs)

        # Add topography in dataset if not in
        if 'h' not in list(self.ds['his'].data_vars):
            self.ds['his']['h']=(['y_rho','x_rho'],  self.H*np.ones((self.M,self.L)))

        if self.verbose:
            print("Grid size: (L ,M, N) = (" + str(self.L) + ", " + str(self.M) + ", " + str(self.N) + ")")

    def getZ(self,pt,zeta=None):
        """ !!! need update !!!
        getZ: Simple function to return vertical grid at specified level
              Valid arguments for pt are 'w','rho','u','v','t','s'
              If zeta is specified, it must be on the rho points
              and will be interpolated to u or v points as needed
        """
        if zeta is None:
            # We'll just return the nominal grid
            h=self.H*np.ones([1,1])
        else:
            h=self.H*np.ones(zeta.shape)
        s=s_coordinate(h, self.theta_b, self.theta_s, self.Hc, self.N, zeta=zeta)

        pt=pt.upper()
        if (pt=='RHO') or (pt=='T') or (pt=='S') or  ((pt=='U') and zeta is None) or ( (pt=='V') and zeta is None) :
            return s.z_r[:]
        elif (pt=='U'):
            # assumes zonal periodicity
            return rho2u(s.z_r[:])
        elif (pt=='V'):
            return rho2v(s.z_r[:])
        elif (pt=='W'):
            return s.z_w[:]

    def get_type(self, v):
        """  !!! old code, update or delete !!!
        Returns the grid point type. Queries the hgrid class for horizontal
        point type, and then checks the vertical dimension to identify 'W' arrays.
        """
        pt = self.hgrid.get_htype(v)
        if (pt == 'RHO') and (v.shape[0] == self.Np):
            return 'W'
        else:
            return pt
