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
        #self._readgrid(**grid_params)    # Read the horizontal grid

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
        offsets = [self.t0[x[0]] for x in ncset]
        files = [os.path.join(self.dirname, x[0], x[1]) for x in ncset]
        datasets = []
        for f, dt in zip(files, offsets):
            try:
                ds = xr.open_dataset(f, chunks={'time_counter': 1, 's_rho': 1})
            except ValueError:
                ds = xr.open_dataset(f, chunks={'time_counter': 1})
            if 'time_counter' in ds:
                t = ds['time_counter']
                t0 = t.isel(time_counter=0)
                t = ((t-t0)/np.timedelta64(1, 's') + dt)*second2day
                ds['time_counter'] = t
            if 'time_center' in ds:
                t = ds['time_center']
                t = ((t-t0)/np.timedelta64(1, 's') + dt)*second2day
                ds['time_center'] = t
            # drop coordinates for easier concatenation
            #ds = ds.drop([k for k in ds.coords \
            #                if k not in ['time_counter','time_centered']])
            datasets.append(ds)
        _ds = xr.concat(datasets, dim='time_counter',
                        coords='minimal', compat='override')
        _ds = self._adjustgrid(_ds)
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

    def _adjustgrid(self, ds):
        for c in ds.coords:
            new_c = c.replace('nav_lat','eta').replace('nav_lon','xi')
            ds = ds.rename({c:new_c})
        # fills in grid parameters, f, f0, beta
        if 'f0' in self._grid_params:
            ds['f0'] = self._grid_params['f0']
        if 'beta' in self._grid_params:
            ds['beta'] = self._grid_params['beta']
            ds = ds.assign_coords(f=ds.beta*ds.eta_rho+ds.f0)
        return ds

    def _readgrid(self, check=False):
        # Synthesize the x_vert and y_vert that the pyroms CGrid class requests
        # Get 1D vector of x and y points on cell box sides
        xx=self.his.variables['nav_lon_dom_U'][0,:]
        yy=self.his.variables['nav_lat_dom_V'][:,0]
        # Extrapolate by one at each end (assumes equal grid spacing at the ends)
        xxx=np.append(np.append(2*xx[0]-xx[1],xx),2*xx[-1]-xx[-2])
        yyy=np.append(np.append(2*yy[0]-yy[1],yy),2*yy[-1]-yy[-2])

        x_vert, y_vert = np.meshgrid(xxx,yyy)
        self.hgrid = CGrid(x_vert, y_vert)

        if check:
            # Verify that CGrid class regenerated the right grid
            a=np.max(np.abs(self.hgrid.x_u-self.his.variables['nav_lon_dom_U'][:]))
            b=np.max(np.abs(self.hgrid.y_u-self.his.variables['nav_lat_dom_U'][:]))
            c=np.max(np.abs(self.hgrid.x_v-self.his.variables['nav_lon_dom_V'][:]))
            d=np.max(np.abs(self.hgrid.y_v-self.his.variables['nav_lat_dom_V'][:]))
            e=np.max(np.abs(self.hgrid.x_rho-self.his.variables['nav_lon_dom_T'][:]))
            f=np.max(np.abs(self.hgrid.y_rho-self.his.variables['nav_lat_dom_T'][:]))
            if np.any([a,b,c,d,e,f]):
                print("CGrid synthesis failure!")
            else:
                print("CGrid synthesis successful!")

        # Store grid sizes
        self.L = len(xx)
        self.M = len(yy)
        self.Lm = self.L - 1
        self.Mm = self.M - 1
        sr=self.his.variables['s_r'][:]
        self.N  = len(sr)
        self.Np = self.N + 1

        # Copy the horizontal dimenions to the hgrid object
        self.hgrid.L  = self.L
        self.hgrid.Lm = self.Lm
        self.hgrid.M  = self.M
        self.hgrid.Mm = self.Mm

        # fills in grid parameters, f, f0, beta
        #for _p, _v in params():
        #    setattr(self,_p,_v)
        #if 'beta' in params():
        #    self.hgrid.f = beta*(self.hgrid.y_rho-np.mean(self.hgrid.y_rho[:,0]))
        #if 'f0' in params():
        #    self.hgrid.f += self.f0

        if self.verbose:
            print("Grid size: (L ,M, N) = (" + str(self.L) + ", " + str(self.M) + ", " + str(self.N) + ")")

    def getXY(self, pt):
        """
        getXY: Simple function to return x and y grid points at specified cell
               position. Valid arguments are 'T','rho','U','V','psi'
        """
        pt=pt.upper()
        if (pt=='T') or (pt=='S') or (pt=='RHO'):
            return self.hgrid.x_rho, self.hgrid.y_rho
        elif (pt=='U'):
            return self.hgrid.x_u, self.hgrid.y_u
        elif (pt=='V'):
            return self.hgrid.x_v, self.hgrid.y_v
        elif (pt=='PSI'):
            return self.hgrid.x_psi, self.hgrid.y_psi
        else:
            return [], []

    def getZ(self,pt,zeta=None):
        """
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
        """
        Returns the grid point type. Queries the hgrid class for horizontal
        point type, and then checks the vertical dimension to identify 'W' arrays.
        """
        pt = self.hgrid.get_htype(v)
        if (pt == 'RHO') and (v.shape[0] == self.Np):
            return 'W'
        else:
            return pt

    def getstats(self, *args):
        """
        getstats: Returns a tuple of diagnostics from the roms log (output.mpi)
                  Arguments are strings indicating the column titles
        """
        out=[]
        for carg in args:
            idx = [i for i, x in enumerate(self.statnames) if x==carg]
            out.append(self.statdata[:,idx])
        if len(out)==1:
            return out[0]
        else:
            return tuple(out)

    def getini(self, config=1):
        """
        getini: return jet initial conditions
        """

        # init grid, grid will be updated later on
        u=np.zeros((self.N,self.M+1,self.L))
        v=np.zeros((self.N,self.M,self.L+1))
        zeta=np.zeros((self.M+1,self.L+1))
        z_r=self.getZ('RHO', zeta=zeta)
        z_w=self.getZ('W', zeta=zeta)

        # a building function
        def asymth(zz,zasym,dzasym):
            return (zz-zasym)*(1.0+0.5*((zz-zasym)+np.abs(zz-zasym))/dzasym) \
                    +zasym

        # from ana_initial.F
        # southern profile
        rhomaxs=27.75;
        bgdrhos=9.8e-6;
        zs1=-1000;
        dzs=700;
        # drhos=1.4;
        # northern profile
        rhomaxn=27.7573;
        bgdrhon=9.8e-6;
        zn1=-400;
        dzn=300;
        # surface Charney mode
        # drhosfs=1.5;
        # drhosfn=0.0;
        # z0=-300;
        drhosf=0.00;
        # z0p=-110;
        # aponte jet
        # z0p=z0;
        # alpha1=0.0075;

        if config==1:
            # cfg1 ~ setup 1 (albeit for alpha1 and z0p)
            # strong charney mode
            print('Configuration 1: strong charney mode (~ setup 1)')
            drhos=1.4;
            drhosfs=1.5;
            drhosfn=0.0;
            alpha1=0.0075;
            alpha2=1.;
            z0=-300;
            z0p=z0;
            #
            drhosfs1=drhosfs*alpha1;
            drhosfs2=drhosfs*alpha2;
            drhosfn1=drhosfn*alpha1;
            drhosfn2=drhosfn*alpha2;
        elif config==2:
            # cfg2 ~ setup 2 (albeit for alpha1
            # Phillips, weak charney mode
            print('Configuration 2: Phillips, weak charney mode (~ setup 2)')
            drhos=1.4;
            drhosfs=1.5;
            drhosfn=0.0;
            alpha1=0.0075;
            alpha2=0;
            z0=-300;
            z0p=z0;
            #
            drhosfs1=drhosfs*alpha1;
            drhosfs2=drhosfs*alpha2;
            drhosfn1=drhosfn*alpha1;
            drhosfn2=drhosfn*alpha2;
        elif config==3:
            # cfg 3
            # setup 1: Phillips>Charney
            print('Configuration 3: Phillips>Charney (setup 1)')
            drhos=1.4;
            drhosfs1=0; # drhosf
            drhosfn1=0; # drhosf
            drhosfs2=1.5; # drhosfs
            drhosfn2=0.0; # drhosfn
            #drhosfs=0.0; # test
            #drhosfn=1.5; # test
            z0=-300;
            z0p=-110;
        elif config==4:
            # cfg 4
            # setup 2: Phillips
            print('Configuration 4: Phillips (setup 2)')
            drhos=1.4;
            drhosfs1=0.; # drhosf
            drhosfn1=0.; # drhosf
            drhosfs2=0.; # drhosfs
            drhosfn2=0.; # drhosfn
            z0=-300;
            z0p=-110;
        elif config==5:
            # cfg 5
            # setup 4: Phillips + surface temperature anomaly
            print('Configuration 5: Phillips + surface temperature anomaly (setup 4)')
            drhos=1.4;
            drhosfs1=0.71; # drhosf
            drhosfn1=0.71; # drhosf
            drhosfs2=1.5; # drhosfs
            drhosfn2=0.; # drhosfn
            z0=-300;
            z0p=-110;
        elif config==6:
            # cfg 6
            # setup 5: strong charney
            print('Configuration 6: strong charney (setup 5)')
            drhos=0.;
            drhosfs1=0.; # drhosf
            drhosfn1=0.; # drhosf
            drhosfs2=2.5; # drhosfs
            drhosfn2=1.0; # drhosfn
            z0=-600;
            z0p=-110;

        # zonal perturbation
        cff_perturb=0.02;

        ### ! --- Build northern and southern density profiles ---

        # ! First compute background density profiles associated with
        # ! a gentle stratification that does not depend
        # ! on jet side. It is there to ensure static stability.

        i=0;j=0; # zk = z_r(istr,jstrR,k)
        h = self.H

        rhoprof = np.zeros((self.N,2))
        rhoprof[:,0]=rhomaxs-bgdrhos*(z_r[:,j,i]+h);
        rhoprof[:,1]=rhomaxn-bgdrhon*(z_r[:,j,i]+h);
        rho0=rhoprof[:]

        # ! Second, get main north/south contrast with a distorded
        # ! tanh shape of variable amplitude.Distorsion of the tanh
        # ! is done with the asym function that increases
        # ! stratification in the upper ocean

        dzs_a=1.3*dzs;
        zk=asymth(z_r[:,j,i],zs1,dzs_a);
        rhoprof[:,0]= rhoprof[:,0]-drhos*(0.5+0.5*np.tanh((zk-zs1)/dzs));

        dzn_a=1.3*dzn;
        drhon=-(rhoprof[self.N-1,0]-rhoprof[self.N-1,1])/(0.5+0.5*np.tanh((zk[-1]-zn1)/dzn));
        zk=asymth(z_r[:,j,i],zn1,dzn_a);
        rhoprof[:,1]= rhoprof[:,1]-drhon*(0.5+0.5*np.tanh((zk-zn1)/dzn));

        rho1 = rhoprof[:];

        zk=z_r[:,j,i];
        rhoprof[:,0] += -drhosfs1*(np.exp((zk-z0p)/np.abs(z0p)))/(np.exp(1.));
        #      -drhosfs*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.));
        rhoprof[:,0] += -drhosfs2*0.5*(1+np.tanh((zk-z0)/np.abs(z0)))/np.tanh(1.);
        #      -drhosfs*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.);
        rhoprof[:,1] += -drhosfn1*(np.exp((zk-z0p)/np.abs(z0p)))/(np.exp(1.));
        #      -drhosfn*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.));
        rhoprof[:,1] += -drhosfn2*0.5*(1+np.tanh((zk-z0)/np.abs(z0)))/np.tanh(1.);
        #      -drhosfn*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.);

        # averages the profiles in order to have smaller difference
        # and ultimately a weaker jet
        print("Jet ywidth is ", self.jet_ywidth/1e3, "km")
        print("Jet weight is ", self.jet_weight)
        rhoprof[:,1] = self.jet_weight * rhoprof[:,1] \
                        + ( 1 - self.jet_weight) * rhoprof[:,0];

        rho2 = rhoprof[:];

        # horizontal indices
        x = (np.arange(-1,self.Lm+1)+0.5)/(self.Lm-2.0)
        y = (np.arange(-1,self.Mm+1)+0.5)/self.Mm -0.5
        y = np.tile(y.reshape((1,y.size,1)),(self.N,1,self.L+1))
        x = np.tile(x.reshape((1,1,x.size)),(self.N,self.M+1,1))


        # if flag_jet_perturb:
        y = y + cff_perturb*np.exp(z_r/1000.) \
                          *np.exp(-(x-0.5)**2/0.05) \
                          *( 0.5*np.sin(2.*np.pi*x) \
                            +0.5*np.sin(6.*np.pi*x) );

        y = y*np.pi*self.hgrid.el/self.jet_ywidth + np.pi/2;

        Fyz=0.5-(y-np.sin(y)*np.cos(y))/np.pi;
        Fyz[np.where(y<0)]=0.5;
        Fyz[np.where(y>np.pi)]=-0.5;
        Fyz[:]=Fyz[:]+0.5;

        rhos = np.tile(rhoprof[:,0].reshape((self.N,1,1)),(1,self.M+1,self.L+1))
        rhon = np.tile(rhoprof[:,1].reshape((self.N,1,1)),(1,self.M+1,self.L+1))
        rho = Fyz*rhos + (1.-Fyz)*rhon;

        # compute pressure and pressure gradients
        flag_flux=False;
        dpdx,dpdy,P = get_gradp(rho,z_r,z_w,self,flag_flux);

        # adjust sea level in order to have zero pressure signal at the bottom
        zeta = -P[0,:,:]/9.81; # P is already divided by rho0
        zeta += -zeta[:].mean();

        # recompute depth and pressure gradients
        z_r=self.getZ('RHO', zeta=zeta)
        z_w=self.getZ('W', zeta=zeta)
        dpdx,dpdy,P = get_gradp(rho,z_r,z_w,self,flag_flux);

        # get Coriolis parameter
        f=self.hgrid.f[:]

        # compute geostrophic velocities
        #print dpdx.shape, dpdy.shape
        u[:,1:-1,:] = (dpdy[:,:-1,:-1]+dpdy[:,:-1,1:]+ \
                       dpdy[:,1:,:-1]+dpdy[:,1:,1:])*.25;
        u[:,0,:] = u[:,1,:]
        u[:,-1,:] = u[:,-2,:]
        u = u/f[:,:-1]
        #u = np.concatenate((u[:,[0],:],u,u[:,[-1],:]),axis=1)/f[:,:-1];
        v[:,:,1:-1] = (dpdx[:,:-1,:-1]+dpdx[:,:-1,1:]+ \
             dpdx[:,1:,:-1]+dpdx[:,1:,1:])*.25;
        # assumes zonal periodicity
        v[:,:,0]=v[:,:,-2]
        v[:,:,-1]=v[:,:,1]
        #print v.shape, f.shape
        #v = np.concatenate((v[:,[0],:],v),axis=1)/f[:,:-1];

        return zeta, rho, u, v

def getrnames(argv,rpath=None):
    """
    Helper function to select run names
    Usage: Put this in file.py:
        rpaths, rnames = getrnames(sys.argv)
    And run it as
        python file.py 2km
        python file.py seq
    """
    if rpath is not None:
        rnames = argv
        rpaths  = [os.path.join(rpath,x) for x in rnames]   # complete paths
        return rpaths, rnames

    # Find the lists of run names and paths (rnames/rpaths)
    rpath="/home2/pharos/othr/aponte/roms_ird/"
    rnames=[]
    if "allruns" in argv: # All runs
        rnames  = [os.path.join("caparmor",x) for x in os.listdir(os.path.join(rpath, "caparmor"))]
        rnames += [os.path.join("fermi",x) for x in os.listdir(os.path.join(rpath, "fermi"))]
        rnames += ["fermi/jet_cfg1_wp5_4km/followup"]
        rnames  = [x for x in rnames if "jet_cfg1" in x and not "1km" in x] # include cfg1, exclude 1km cases

    elif "seq" in argv:  # The official sequence for the paper
        rnames.append("fermi/jet_cfg1_wp5_4km")              #    0-1500 d
        rnames.append("fermi/jet_cfg1_wp5_2km_k1.0e8")       # 1500-2000 d
        rnames.append("fermi/jet_cfg1_wp5_2km_decay_itide")  # 2000-4000 d

    elif "4km" in argv:  # 4km for rapid testing
        rnames.append("caparmor/jet_cfg1_wp5_4km_decay_itide")

    elif "2km" in argv: # shorthand for the usual case
        rnames.append("fermi/jet_cfg1_wp5_2km_decay_itide")

    elif '/' in argv[1]: # specify the complete rname
        if os.path.isdir(os.path.join(rpath,argv[1])):
            rnames.append(argv[1])

    elif '_' in argv[1]: # specify the suffix following jet_cfg1_
        tmp = "fermi/jet_cfg1_" + argv[1]
        if os.path.isdir(os.path.join(rpath,tmp)):
            rnames.append(tmp)

    else:
        print("Unknown case specified")

    rpaths  = [os.path.join(rpath,x) for x in rnames]   # complete paths
    return rpaths, rnames
