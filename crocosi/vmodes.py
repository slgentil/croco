""" vmodes.py python module
WARNING this is a WIP, see ipython notebook for the current version under development
this version is not working with the class Vmodes yet. """
import scipy.sparse as sp
import scipy.sparse.linalg as la
import numpy as np
import xarray as xr
import crocosi.gridop as gop

from .postp import g_default

# default values
_sig = .1
_nmodes = 10

### N.B.: norm is H

# This should be a class
class Vmodes(object):
    """ Description of the class
    """
    def __init__(self, xgrid, nmodes,
                 zc, zf, N2,
                 free_surf=True,
                 persist=False,
                 g=g_default, sigma=_sig):
        """ Create a Vmode object
        
        Parameters
        ----------
            xgrid: xgcm.Grid
                xgcm grid object required for grid manipulations
            ...
            
        """
        self.xgrid = xgrid
        self.nmodes = nmodes
        self._znames = {"zc": zc.name, "zf": zf.name}
        self.g = g
        self.free_surf = free_surf
        self.sigma=_sig
        
        # create dataset
        N2name = N2.name if N2.name else "N2"
        # time_instant is conflicting if specifying coords in Dataset... weird
        self.ds = xr.Dataset({N2name:N2}).assign_coords({zc.name:zc})
        if zf.name not in self.ds.coords:
            self.ds = self.ds.assign_coords({zf.name:zf})
        self._compute_vmodes()
        
        self.ds['H'] = np.abs(zf.sel(s_w=-1,method='nearest'))
        self.ds['dz'] = xgrid.diff(zf, "s").rename("dz")
        
        if persist:
            self.ds.persist()
     
    def __getitem__(self, item):
        """ Enables calls such as vm['N2']
        """
        if item in ['zc', 'zf']:
            return self.ds[self._znames[item]]
        return self.ds[item]

    def _compute_vmodes(self):
        """ compute vertical modes and store the results into the dataset """
        dm = get_vmodes(self['zc'],
                        self['zf'],
                        self.ds.N2,
                        nmodes=self.nmodes,
                        free_surf=self.free_surf,
                        sigma=self.sigma)
        self.ds = xr.merge([self.ds, dm], compat="override")
        pass
    
#############################################################################
####### projection and reconstruction #######################################
    
    def project(self, data, vartype="p", **kwargs):
        """ Project a variable on vertical modes (p-modes or w-modes)
        Internally call project_puv or project_w
        
        Parameters:
        ___________
        data: xarray.DataArray
            array containing the data to be projected. 
        vartype: str, optional (default: "p", i.e. pressure modes)
            string specifying whether projection should be done w-modes ("w"), buoyancy modes ("b") or pressure modes (any other value)
        zz: xarray.DataArray or str or bool, optional (default: False)
            array containing the z-grid of the data, or string containing the name of the z coord, or boolean saying wether we should interpolate (finding the z-coord by its own). The data will be interpolated onto the vmodes grid points prior to projection.
        isel: Dict, optional (default: None)
            indices applied to the vmodes dataset prior to projection
        align: bool, optional (default: True)
            wether alignment between the data and vmodes DataArray should be performed before projecting
        
        Returns
        _______
        xarray.DataArray
            Projected array
            
        See also
        ________
        project_puv, project_w, project_b, reconstruct
        
        """
        if vartype == "w":
            return self.project_w(data, **kwargs)
        elif vartype == "b":
            return self.project_b(data, **kwargs)
        else:
            return self.project_puv(data, **kwargs)

    def project_puv(self, data, zz=False, isel=None, align=True):
        """ projection on p-mode 
        compute int_z(phi_n * field) using sum
        
        Parameters:
        ___________
        data: xarray.DataArray
        
        Returns:
        ________
        xarray.Datarray
        
        See also:
        _________
        project, reconstruct_p
        
        """

        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)

        if align:
            data, dm = xr.align(data, dm, join="inner")
        if not( zz is None or zz is False ):
            if zz is True:
                zz = _get_z_coord(data)
            if isinstance(zz, str):
                zz = data.coords[zz]
            elif align:
                data, zz = xr.align(data, zz, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], zz, data)
        res = (dm.dz*data*dm.phi).sum("s_rho")/dm.norm

        return res
    
    def project_w(self, data, zz=False, isel=None, align=True): 
        """ projection on w-mode (varphi = -c**2/N2 * dphidz)
        for reconstructing, use w = wn*varphi (see reconstruct_w)
        
        interpolation uses linear interpolation, but midpoints should be OK 
            (since it gives something equivalent to trapezoidal integration upon integration)
            
        Parameters:
        ___________
        data: xarray.DataArray
        
        Returns:
        ________
        xarray.Datarray
        
        See also:
        _________
        project, reconstruct_w
        
        """
        
        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)
        if align:
            data, dm = xr.align(data, dm, join="inner")    
        if not( zz is None or zz is False ):
            if zz is True:
                zz = _get_z_coord(data)
            if isinstance(zz, str):
                zz = data.coords[zz]
            elif align:
                data, zz = xr.align(data, zz, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], zz, data)
        zf, zc = self._znames['zf'], self._znames["zc"]
        prov = (data * self._w2rho(-dm.dphidz, zc=dm[zc], zf=dm[zf]) * dm.dz).sum(dim="s_rho")
        if self.free_surf:
            prov += self.g * ( -dm.dphidz / dm.N2 
                              * self.xgrid.interp(data, "s", boundary="extrapolate") 
                             ).isel(s_w=-1).drop(self._znames["zf"])
       
        return prov/dm.norm
    
    def project_b(self, data, zz=False, isel=None, align=True): 
        """ projection on b-mode (dphidz)
        for reconstructing, use -c**2*bn*dphidz (see reconstruct b)
        
        interpolation uses linear interpolation, but midpoints should be OK 
            (since it gives something equivalent to trapezoidal integration upon integration)
            
        Parameters:
        ___________
        data: xarray.DataArray
        
        Returns:
        ________
        xarray.Datarray
        
        See also:
        _________
        project, project_w, reconstruct_b
        
        """
        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)
        if align:
            data, dm = xr.align(data, dm, join="inner")    
        if not( zz is None or zz is False ):
            if zz is True:
                zz = _get_z_coord(data)
            if isinstance(zz, str):
                zz = data.coords[zz]
            elif align:
                data, zz = xr.align(data, zz, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], zz, data)
        zf, zc = self._znames['zf'], self._znames["zc"]
        prov = (data * self._w2rho(dm.dphidz/dm.N2, zc=dm[zc], zf=dm[zf])* dm.dz).sum("s_rho")
        if self.free_surf:
            prov += self.g * (self.xgrid.interp(data, "s", boundary="extrapolate")
                          * dm.dphidz / dm.N2**2
                         ).isel(s_w=-1).drop(self._znames["zf"])
        return -prov/dm.norm 

    def reconstruct(self, modamp, vartype=None, **kwargs):
        """ Reconstruct a variable from modal amplitudes
        Internally call reconstruct_puv or reconstruct w or reconstruct_b
        
        Parameters:
        ___________
        modamp: xarray.DataArray
            array containing the modal amplitudes 
        vartype: str, optional, {"p", "u", "v", "w", "b"} (default: "p", i.e. pressure modes)
            string specifying whether reconstruction should be done using w-modes ("w"), buoyancy modes ("b") or pressure modes ("p", "u" or "v")
        isel: Dict, optional (default: None)
            indices applied to the vmodes dataset prior to reconstruction
        
        Returns
        _______
        xarray.DataArray
            Reconstructed field array
            
        See also
        ________
        reconstruct_puv, reconstruct_w, reconstruct_b, project
        
        """
        
        if vartype is None:
            vartyps = "puvbw"
            if sum(s in modamp.name.lower() for s in vartyps)==1:
                vartype = next((s for s in vartyps if s in modamp.name.lower()))
            else: 
                raise ValueError("unable to find what kind of basis to use for reconstruction")
                
        if vartype in "puv":
            return self.reconstruct_puv(modamp, **kwargs)
        elif vartype == "w":
            return self.reconstruct_w(modamp, **kwargs)
        elif vartype == "b":
            return self.reconstruct_b(modamp, **kwargs)
        
    def reconstruct_puv(self, modamp, isel=None, align=True):
        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)
        if align:
            modamp, dm = xr.align(modamp, dm, join="inner")    
        return (modamp * dm.phi).sum("mode")

    def reconstruct_w(self, modamp, isel=None, align=True):
        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)
        if align:
            modamp, dm = xr.align(modamp, dm, join="inner")    
        return (-dm.c**2 / dm.N2 * dm.dphidz * modamp).sum("mode")
    
    def reconstruct_b(self, modamp, isel=None, align=True):
        if isel is None:
            dm = self.ds
        else:
            dm = self.ds.isel(isel)
        if align:
            modamp, dm = xr.align(modamp, dm, join="inner")    
        return (-dm.c**2 * dm.dphidz * modamp).sum("mode")
    ### utilitaries 
    
    def _w2rho(self, data, zf=None, zc=None, align=True):
        """ routine to interpolate fields from mean rho point to mean w point
        adapt to different ensemble of points if specified
        wrapper to gop.w2rho
        could be replaced by xgrid.interp one correct metrics are provided
        """
        if zf is None:
            if align:
                zf, data = xr.align(self["zf"], data, join="inner")
            else:
                zf= self["zf"]
        if zc is None:
            if align:
                zc, data = xr.align(self["zc"], data, join="inner")
            else:
                zc = self["zc"]
        return gop.w2rho(data, self.xgrid, zc, zf)

def _get_z_coord(ds):
    """ return name of z coordinate """
    if "s_rho" in ds.dims:
        zname = next((z for z in ['z_r', 'z_rho'] if z in ds.coords), None)
    elif "s_w" in ds.dims:
        zname = "z_w" if "z_w" in ds.coords else None
    elif "z" in ds.dims:
        zname = "z" if "z" in ds.coords else None
    else:
        raise ValueError("could not find z dimension")
    if zname is None:
        raise ValueError("could not find z coordinate")
    return zname

def get_vmodes(zc, zf, N2, nmodes=_nmodes, **kwargs):
    """ wrapper for calling compute_vmodes with apply_ufunc. Includes unstacking of result
    this is what you should need in your scripts, unless you are using numpy arrays only 
    input:
        - zc: z-levels at center of cells
        - zf: z-levels at cell faces
        - N2: BVF at cell faces 
        - nmodes (default:10): number of vertical modes (+ barotropic)
    z levels must be in ascending order (first element is at bottom, last element is at surface)
    output: xarray dataset with phase speed, pressure-modes (at center) and b-modes (at faces)
    kwargs:
        - free_surf (bool): use free surface boundary condition (default:True)
        - sigma (scalar or None): for shift-invert in eig (default: 0.1)
    """
    kworg = {"nmodes":nmodes, "stacked":True}
    if kwargs is not None:
        kworg.update(kwargs)
    N = zc.s_rho.size
    res = xr.apply_ufunc(compute_vmodes, 
                         zc.chunk({"s_rho":-1}), 
                         zf.chunk({"s_w":-1}),
                         N2.chunk({"s_w":-1}), 
                         kwargs=kworg, 
                         input_core_dims=[["s_rho"],["s_w"],["s_w"]],
                         dask='parallelized', 
                         output_dtypes=[np.float64],
                         output_core_dims=[["s_stack","mode"]],
                         output_sizes={"mode":nmodes+1,"s_stack":2*(N+1)}
                        )
    res['mode'] = np.arange(nmodes+1)
    # unstack variables
    c = res.isel(s_stack=0).rename('c')
    phi = (res.isel(s_stack=slice(1,N+1))
           .rename('phi')
           .rename({'s_stack': 's_rho'})
           #.assign_coords(z_rho=zc)
          )
    dphidz = (res.isel(s_stack=slice(N+1,2*N+2))
              .rename('dphidz')
              .rename({'s_stack': 's_w'})
            #  .assign_coords(z_w=zf)
             )
    # merge data into a single dataset
    other_dims = tuple([dim for dim in zc.dims if dim!="s_rho"]) # extra dims    
    dm = (xr.merge([c, phi, dphidz, -zf.isel(s_w=0, drop=True).rename('norm')])
          .transpose(*('mode','s_rho','s_w')+other_dims)
          #.assign_coords(norm=-zf.isel(s_w=0, drop=True)) #, N2=N2
         )
    return dm  ### hard-coded norm = H

def compute_vmodes(zc_nd, zf_nd, N2f_nd, 
                   nmodes=_nmodes, free_surf=True,
                   g=g_default,
                   sigma=_sig, stacked=True,
                   **kwargs):
    """
    wrapper for vectorizing compute_vmodes_1D over elements of axes other than vertical dim
    that's not elegant, nor efficient
    here z is last axis (because it is core dim)
    you can use this if you are using numpy (but make sure z is last dim)
    """
    assert zc_nd.ndim==zf_nd.ndim==N2f_nd.ndim
    assert zf_nd.shape==N2f_nd.shape
    if "nmodes" in kwargs:
        nmodes = kwargs["nmodes"]
    if "free_surf" in kwargs:
        free_surf = kwargs["free_surf"]
        if "stacked" in kwargs:
            stacked = kwargs["stacked"]

    if zc_nd.ndim>1:
        nxy = zc_nd.shape[:-1]
        nn = np.prod(nxy)
        zc_nd = zc_nd.reshape(nn,zc_nd.shape[-1])
        zf_nd = zf_nd.reshape(nn,zf_nd.shape[-1])
        N2f_nd = N2f_nd.reshape(nn,N2f_nd.shape[-1])
        ii = 0
        cn, phin, dphi = compute_vmodes_1D(zc_nd[ii,:], zf_nd[ii,:], N2f_nd[ii,:], nmodes, \
                                free_surf=free_surf, g=g, sigma=sigma)
        cn, phin, dphi = cn[None,:], phin[None,...], dphi[None,...]
        for ii in range(1,nn):
            res = compute_vmodes_1D(zc_nd[ii,:], zf_nd[ii,:], N2f_nd[ii,:], nmodes, \
                                free_surf=free_surf, g=g, sigma=sigma)
            cn = np.vstack([cn, res[0][None,...]])
            phin = np.vstack([phin,res[1][None,...]])
            dphi = np.vstack([dphi,res[2][None,...]])
        if stacked:
            return np.hstack([cn[:,None,:],phin,dphi]).reshape(nxy+(-1,nmodes+1))
        else:
            return cn.reshape(nxy+(nmodes,)), phin.reshape(nxy+(-1,nmodes)), \
                        dphi.reshape(nxy+(-1,nmodes))
    else:
        return compute_vmodes_1D(zc_nd, zf_nd, N2f_nd, nmodes, \
                                free_surf=free_surf, g=g, sigma=sigma)

            
def compute_vmodes_1D(zc, zf, N2f, 
                      nmodes=_nmodes, free_surf=True, 
                      g=g_default, sigma=_sig):
    """ compute vertical modes: solution of SL problem (phi'/N^2)'+k*phi=0'
    returns phi at rho points, dphi at w points and c=1/sqrt(k) 
    normalization such that int(phi^2)=H, w-modes=d(phi)/dz
    copy-pasted from M. Dunphy's vmodes_MD.py
    TODO: correct rigid lid & barotropic mode """
    # Precompute a few quantities
    assert zc.ndim==zf.ndim==N2f.ndim==1
    assert len(zc)+1==len(N2f)==len(zf)
    dzc=np.diff(zc)
    dzf=np.diff(zf)
    Np=len(zf) #self.zf)
    N20=N2f[-1] #self.N2f[-1]  # N^2(z=0)
    H = abs(zf[0])
    
    # Build Dz, C-to-F grid 
    v1=-1.0/np.concatenate([dzc,np.ones(1)])
    v2= 1.0/np.concatenate([np.ones(1),dzc])
    v12=np.stack([v1, v2])
    Dz=sp.spdiags(v12,[-1, 0],Np,Np-1,format="lil")
    # Adjust matrix for BCs
    Dz[0,:]=0
    Dz[-1,:]=0
    if free_surf:
        Dz[-1,-1]=np.divide(-N20, g + N20*(zf[-1] - zc[-1]))

    # Build Dz2, F-to-C grid
    v1=-1.0/np.concatenate([dzf,np.ones(1)])
    v2= 1.0/np.concatenate([np.ones(1),dzf])
    v12=np.stack([v1,v2])
    Dz2=sp.spdiags(v12,[0, 1],Np-1,Np,format="lil")
    
    # Construct A, solve eigenvalue problem
    iN2=sp.spdiags(1.0/N2f,0,Np,Np)
    A=-Dz2*iN2*Dz
    ev,ef = la.eigs(A.tocsc(),nmodes+1,sigma=sigma)
    ev,ef = np.real(ev), np.real(ef)
    
    # Convert eigvenvalues to c_e, sort appropriately
    c=1.0/np.sqrt(np.real(ev))
    ii=(-c).argsort()
    c=c[ii]               # c_e value
    phic=ef[:,ii]         # phi at cell centres

    # Normalize and set surface value positive
    for mi in range(nmodes+1):
        fn=phic[:,mi]        # current phi
        s=np.sign(fn[-1])         # sign at surface
        if s==0:
            s=1;
        tmp = np.sum((fn**2)*dzf)/H
        phic[:,mi] = s*fn/np.sqrt(tmp) # (1/H)*\int_{-H}^{0} \phi_m(z)^2 dz = 1
        
    # dphi at cell faces: phi'=dphidz (buoyancy modes)
    dphif = Dz*phic
    # this would give w-modes: np.r_[np.zeros((1,nmodes+1)),(dzf[:,None]*phic).cumsum(axis=0)]
    return c, phic, dphif

