import sys
import xarray as xr
import numpy as np

# ------------------------------------------------------------------------------------------------

def rho2u(v, ds):
    """
    interpolate horizontally variable from rho to u point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'X')
    add_coords(ds, var, ['xi_u','eta_u'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def u2rho(v, ds):
    """
    interpolate horizontally variable from u to rho point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'X')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def v2rho(v, ds):
    """
    interpolate horizontally variable from rho to v point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'Y')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def rho2v(v, ds):
    """
    interpolate horizontally variable from rho to v point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'Y')
    add_coords(ds, var, ['xi_v','eta_v'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def rho2psi(v, ds):
    """
    interpolate horizontally variable from rho to psi point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'X')
    var = grid.interp(var,'Y')
    add_coords(ds, var, ['xi_u','eta_v'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def psi2rho(v, ds):
    """
    interpolate horizontally variable from rho to psi point
    """
    grid = ds.attrs['xgcm-Grid']
    var = grid.interp(v,'X')
    var = grid.interp(var,'Y')
    add_coords(ds, var, ['xi_rho','eta_rho'])
    var.attrs = v.attrs
    var.name = v.name
    return var

def get_z(run, zeta=None, h=None, vgrid='r', hgrid='r'):
    ''' compute vertical coordinates
        zeta should have the size of the final output
        vertical coordinate should be first
    '''

    ds = run.ds['his']
    N = run.N
    hc = run.params['Hc']

    if zeta is not None:
        _zeta = zeta
    else:
        _zeta = run.ds['his'].ssh
    if h is not None:
        _h = h
    else:
        _h = run.ds['his'].h
    #
    if hgrid is 'u':
        _zeta = rho2u(_zeta, ds)
        _h = rho2u(_h, ds)
    elif hgrid is 'v':
        _zeta = rho2v(_zeta, ds)
        _h = rho2v(_h, ds)
    #
    sc=run.ds['his']['sc_'+vgrid]
    cs=run.ds['his']['Cs_'+vgrid]

    #
    z0 = (hc * sc + _h * cs) / (hc + _h)
    z = np.squeeze(_zeta + (_zeta + _h) * z0)
    # manually swap dims, could also perform align with T,S
    if z.ndim == 3:
        z = z.transpose(sc.dims[0], _zeta.dims[0], _zeta.dims[1])
    elif z.ndim == 2:
        z = z.transpose(sc.dims[0], _zeta.dims[0])
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

def interp2z(z0, z, v, extrap):
    ''' Interpolate on a horizontally uniform grid
    '''
    import fast_interp3D as fi  # OpenMP accelerated C based interpolator
    #
    if v.ndim == 1 or z.ndim == 1 :
        lz = z.squeeze()[:,None,None]
        lv = v.squeeze()[:,None,None]
    elif v.ndim == 2 :
        lz = z[...,None]
        lv = v[...,None]
    else:
        lz = z[...]
        lv = v[...]
    #
    if extrap:
        zmin = np.min(z0)-1.
        lv = np.vstack((lv[[0],...], lv))
        lz = np.vstack((zmin+0.*lz[[0],...], lz))
    #
    vi = fi.interp(z0.astype('float64'), lz.astype('float64'), lv.astype('float64'))
    return vi

def interp2z_xr(z0, z, v,  hgrid='rho', extrap=True):
    _xmap = {'rho': 'rho', 'u': 'u', 'v': 'rho'}
    _ymap = {'rho': 'rho', 'u': 'rho', 'v': 'v'}
    _xgrid, _ygrid = _xmap[hgrid], _ymap[hgrid]
    return (xr.DataArray(interp2z(z0.values, z.values, v.values, extrap),
                         dims=['z_r','eta_'+_ygrid,'xi_'+_xgrid],
                         coords={'z_r': z0.values, 'xi_'+_xgrid: v['xi_'+_xgrid],
                                 'eta_'+_ygrid: v['eta_'+_ygrid]}) )
