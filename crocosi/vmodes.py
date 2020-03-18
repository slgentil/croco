import scipy.sparse as sp
import scipy.sparse.linalg as la
import numpy as np
import xarray as xr

### default values
_g = 9.81
_sig = .1
_nmodes = 10
### N.B.: norm is H

# This should be a class

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
           .assign_coords(z_rho=zc)
          )
    dphidz = (res.isel(s_stack=slice(N+1,2*N+2))
              .rename('dphidz')
              .rename({'s_stack': 's_w'})
              .assign_coords(z_w=zf)
             )
    # merge data into a single dataset
    other_dims = tuple([dim for dim in zc.dims if dim!="s_rho"]) # extra dims    
    dm = (xr.merge([c, phi, dphidz])
          .transpose(*('mode','s_rho','s_w')+other_dims)
          .assign_coords(N2=N2, norm=-zf.isel(s_w=0, drop=True))
         )
    return dm  ### hard-coded norm = H

def compute_vmodes(zc_nd, zf_nd, N2f_nd, nmodes=_nmodes, free_surf=True, g=_g, 
                   sigma=_sig, stacked=True, **kwargs):
    """
    wrapper for vectorizing compute_vmodes_1D over elements of axes other than vertical dim
    that's not elegant, nor efficient
    here z is last axis (because it is core dim)
    you can use this if you are using numpy (but make sure z is last dim)
    """
    assert zc_nd.ndim==zf_nd.ndim==N2f_nd.ndim
    assert zf_nd.shape==N2f_nd.shape
    if "nmodes" in kwargs:
        nmodes = kwarg["nmodes"]
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

            
def compute_vmodes_1D(zc, zf, N2f, nmodes=_nmodes, free_surf=True, g=_g, sigma=_sig):
    """ compute vertical modes: solution of SL problem (phi'/N^2)'+k*phi=0'
    returns phi at rho points, dphi at w points and c=1/sqrt(k) 
    normalization such that int(phi^2)=H, w-modes=d(phi)/dz
    copy-pasted from M. Dunphy's vmodes_MD.py
    TODO: correct rigid lid & barotropic mode """
    print("computing {0} modes", nmodes)
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

#def project(...)
