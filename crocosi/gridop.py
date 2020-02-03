import sys
import xarray as xr
import numpy as np

# ------------------------------------------------------------------------------------------------

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

def get_z(run, zeta=None, h=None, vgrid='r', hgrid='r'):
    ''' compute vertical coordinates
        zeta should have the size of the final output
        vertical coordinate should be first
    '''

    try:
        ds = run.ds['grid']
    except Exception:
        ds = run.ds['his']
    N = run.N
    hc = run.params_input['Hc']

    if h is not None:
        _h = h
    else:
        _h = ds.h

    if zeta is not None:
        _zeta = zeta.values
    else:
        _zeta = np.zeros(h.shape)

    #
    if hgrid is 'u':
        _zeta = rho2u(_zeta, ds)
        _h = rho2u(_h, ds)
    elif hgrid is 'v':
        _zeta = rho2v(_zeta, ds)
        _h = rho2v(_h, ds)
    #
    sc=ds['sc_'+vgrid]
    cs=ds['Cs_'+vgrid]

    #
    z0 = (hc * sc + _h * cs) / (hc + _h)
    z = np.squeeze(np.repeat(_zeta[:, :, np.newaxis], sc.shape[0], axis=2) + (_zeta + _h) * z0)
    # manually swap dims, could also perform align with T,S
    if z.ndim ==4:
        z = z.transpose(z.dims[0], z.dims[3], z.dims[1], z.dims[2])
    elif z.ndim == 3:
        z = z.transpose(z.dims[2], z.dims[0], z.dims[1])
    return z.rename('z_'+vgrid)

def get_p(rho, zeta, rho0, rho_a=None):
    #
    if rho_a is None:
        _rho_a = 0.
    else:
        _rho_a = rho_a
    #
    _rho = g*(rho.shift(z_r=1)+rho)*.5*(rho.z_r.shift(z_r=1)-rho.z_r)
    _rho = _rho.sortby(_rho.z_r, ascending=False).shift(z_r=1).fillna(0.)
    p0 = g*(rho0+_rho_a.sel(z_r=0,method='nearest')+rho.sel(z_r=0,method='nearest'))*zeta
    p = _rho.cumsum('z_r').sortby(_rho.z_r, ascending=True) + p0
    p = p.rename('p')
    return p

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
