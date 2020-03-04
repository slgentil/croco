import sys
import xarray as xr
import numpy as np

# ------------------------------------------------------------------------------------------------
_g = 9.81

def add_coords(ds, var, coords):
    for co in coords:
        var.coords[co] = ds.coords[co]

def rho2u(v, ds):
    """
    interpolate horizontally variable from rho to u point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'xi')
    add_coords(ds, var, ['xi_u','eta_u'])
    var.attrs = v.attrs
    return var.rename(v.name)

def u2rho(v, ds):
    """
    interpolate horizontally variable from u to rho point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'xi')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    return var.rename(v.name)

def v2rho(v, ds):
    """
    interpolate horizontally variable from rho to v point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'eta')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    return var.rename(v.name)

def rho2v(v, ds):
    """
    interpolate horizontally variable from rho to v point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'eta')
    add_coords(ds, var, ['xi_v','eta_v'])
    var.attrs = v.attrs
    return var.rename(v.name)

def rho2psi(v, ds):
    """
    interpolate horizontally variable from rho to psi point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'xi')
    var = grid.interp(var,'eta')
    var.attrs = v.attrs
    return var.rename(v.name)

def psi2rho(v, ds):
    """
    interpolate horizontally variable from rho to psi point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'xi')
    var = grid.interp(var,'eta')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    return var.rename(v.name)

def get_z(run, zeta=None, h=None, vgrid='r', hgrid='r', vtrans=None):
    ''' compute vertical coordinates
        zeta should have the size of the final output
        vertical coordinate is first in output
    '''

    try:
        ds = run.ds['grid']
    except Exception:
        ds = run.ds['his']
    N = run.N
    hc = run.params_input['Hc']

    _h = ds.h if h is None else h
    _zeta = 0*ds.h if zeta is None else zeta

    # swith horizontal grid if needed (should we use grid.interp instead?)
    if hgrid in ['u','v']:
        funtr = eval("rho2"+hgrid)
        if zeta is None:
            _zeta = funtr(_zeta, ds)
        _h = funtr(_h, ds)
    
    # determine what kind of vertical corrdinate we are using (NEW_S_COORD)
    if vtrans is None:
        vtrans = ds.Vtransform.values
    else:
        if isinstance(vtrans, str):
            if vtrans.lower()=="old":
                vtrans = 1
            elif vtrans.lower()=="new":
                vtrans = 2
            else:
                raise ValueError("unable to understand what is vtransform")
                
    sc=ds['sc_'+vgrid]
    cs=ds['Cs_'+vgrid]

    if vtrans == 2:
        z0 = (hc * sc + _h * cs) / (hc + _h)
        z = _zeta + (_zeta + _h) * z0
    else:
        z0 = hc*sc + (_h-hc)*cs
        z = z0 + _zeta*(1+z0/_h)
        
    zdim = "s_"+vgrid.replace('r','rho')
    if z.dims[0] != zdim:
        z = z.transpose(*(zdim,)+_zeta.dims)
    return z.rename('z_'+vgrid)

def get_p(grid,rho,zw,zr=None,g=_g):
    """ compute (not reduced) pressure by integration from the surface, 
    taking rho at rho points and giving results on w points (z grid)
    with p=0 at the surface. If zr is not None, compute result on rho points """
    if zr is None:
        dz = grid.diff(zw, "s")
        p = grid.cumsum((rho*dz).sortby("s_rho",ascending=False), "s",                             to="outer", boundary="fill").sortby("s_w",ascending=False)
    else:
        """ it is assumed that all fields are from bottom to surface"""
        rna = {"s_w":"s_rho"}
        dpup = (zr - zw.isel(s_w=slice(0,-1)).drop("s_w").rename(rna))*rho
        dpdn = (zw.isel(s_w=slice(1,None)).drop("s_w").rename(rna) - zr)*rho
        p = (dpup.shift(s_rho=-1, fill_value=0) + dpdn).sortby(rho.s_rho, ascending=False)                .cumsum("s_rho").sortby(rho.s_rho, ascending=True).assign_coords(z_r=zr)
    return _g *p.rename("p")


def get_uv_from_psi(psi, ds):
    # note that u, v are computed at rho points
    x, y = ds.xi_rho, ds.eta_rho
    #
    u = - 0.5*(psi.shift(eta_rho=1)-psi)/(y.shift(eta_rho=1)-y) \
        - 0.5*(psi-psi.shift(eta_rho=-1))/(y-y.shift(eta_rho=-1))
    #
    v =   0.5*(psi.shift(xi_rho=1)-psi)/(x.shift(xi_rho=1)-x) \
        + 0.5*(psi-psi.shift(xi_rho=-1))/(x-x.shift(xi_rho=-1))
    return u, v

def get_pv(u, v, rho, rho_a, f, f0, zr, zw, ds):

    # relative vorticity
    xi_v =  ( (v.diff('xi_rho')/ds.xi_rho.diff('xi_rho'))
             .rename({'xi_rho':'xi_psi', 'eta_v':'eta_psi'})
             .assign_coords(xi_psi=ds.xi_u.rename({'xi_u':'xi_psi'}),
                            eta_psi=ds.eta_v.rename({'eta_v':'eta_psi'})))

    xi_u = -( (u.diff('eta_rho')/ds.eta_rho.diff('eta_rho'))
             .rename({'eta_rho':'eta_psi', 'xi_u':'xi_psi'})
             .assign_coords(xi_psi=ds.xi_u.rename({'xi_u':'xi_psi'}),
                            eta_psi=ds.eta_v.rename({'eta_v':'eta_psi'})))
    xi = psi2rho(xi_u+xi_v, ds)

    # stretching term
    drho = (rho.shift(z_r=1)+rho)*0.5 * zr.diff('z_r')/rho_a.diff('z_r')
    drho = drho.rename({'z_r':'z_w'}).assign_coords(z_w=zw[:-1])
    Sint = drho.diff('z_w')/zw[:-1].diff('z_w')
    Sint = Sint.rename({'z_w':'z_r'}).assign_coords(z_r=zr[1:-1])

    # note that top and bottom values are not used in the solver
    S = f0 * xr.concat([0.*rho.isel(z_r=0), Sint, 0.*rho.isel(z_r=-1)], dim='z_r') #.transpose('z_r','eta_rho','xi_rho')

    # assemble pb
    q = (xi + S + f - f0 ).rename('q') # order of variable conditions dimension order

    return q

def interp2z_3d(z0, z, v, extrap):
    import crocosi.fast_interp3D as fi  # OpenMP accelerated C based interpolator
    # check v and z have identical shape
    assert v.ndim==z.ndim
    # test if temporal dimension is present
    if v.ndim == 1:
        lv = v.squeeze()[:,None,None]
        lz = z.squeeze()[:,None,None]
    elif v.ndim == 2:
        lv = v[...,None]
        lz = z[...,None]
    else:
        lz = z[...]
        lv = v[...]
    #
    if extrap:
        zmin = np.min(z0)-1.
        lv = np.concatenate((lv[[0],...], lv), axis=0)
        lz = np.concatenate((zmin+0.*lz[[0],...], lz), axis=0)
    #
    return fi.interp(z0.astype('float64'), lz.astype('float64'), lv.astype('float64'))

def interp2z(z0, z, v, extrap):
    ''' interpolate vertically
    '''
    # check v and z have identical shape
    assert v.ndim==z.ndim
    # test if temporal dimension is present
    if v.ndim == 4:
        vi = [interp2z_3d(z0, z[...,t], v[...,t], extrap)[...,None] for t in range(v.shape[-1])]
        return np.concatenate(vi, axis=0) # (50, 722, 258, 1)
        #return v*0 + v.shape[3]
    else:
        return interp2z_3d(z0, z, v, extrap)


def N2Profile(run, strat, z, g=9.81):
    """
    Method to compute the N2 profile : 
    """
    
    grid = run.ds['his'].attrs['xgcm-Grid']
    N2 = -g/run.params_input['rho0'] * grid.diff(strat,'s') / grid.diff(z,'s')
    N2.isel(s_w=0).values = N2.isel(s_w=1).values
    N2.isel(s_w=-1).values = N2.isel(s_w=-2).values
    # if np.any(N2<0):
    #     print("Unstable N2 profile detected")
    return (N2)

def hinterp(ds,var,coords=None):
    import pyinterp
    #create Tree object
    mesh = pyinterp.RTree()

    L = ds.dims['x_r']
    M = ds.dims['y_r']
    N = ds.dims['s_r']
    z_r = get_z(ds)
    #coords = np.array([coords])

    # where I know the values
    z_r = get_z(ds)
    vslice = []
    #lon_r = np.tile(ds.lon_r.values,ds.dims['s_r']).reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
    lon_r = np.tile(ds.lon_r.values,(ds.dims['s_r'],1,1)).reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
    lat_r = np.tile(ds.lat_r.values,(ds.dims['s_r'],1,1)).reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
    z_r = get_z(ds)
    mesh.packing(np.vstack((lon_r.flatten(), lat_r.flatten(), z_r.values.flatten())).T,
                            var.values.flatten())

    # where I want the values
    zcoord = np.zeros_like(ds.lon_r.values)
    zcoord[:] = coords
    vslice, neighbors = mesh.inverse_distance_weighting(
        np.vstack((ds.lon_r.values.flatten(), ds.lat_r.values.flatten(), zcoord.flatten())).T,
        within=True)    # Extrapolation is forbidden)


    # The undefined values must be set to nan.
    print(ds.mask_rho.values)
    mask=np.where(ds.mask_rho.values==1.,)
    vslice[int(ds.mask_rho.values)] = float("nan")

    vslice = xr.DataArray(np.asarray(vslice.reshape(ds.dims['y_r'], ds.dims['x_r'])),dims=('x_r','y_r'))
    yslice = ds.lat_r
    xslice = ds.lon_r

    return[xslice,yslice,vslice]
