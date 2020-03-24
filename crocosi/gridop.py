import sys
import xarray as xr
import numpy as np
from collections import OrderedDict

from .postp import g_default

# ------------------------------------------------------------------------------

def _get_spatial_dims(v):
    """ Return an ordered dict of spatial dimensions in the s/z, y, x order
    """
    dims = OrderedDict( (d, next((x for x in v.dims if x[0]==d), None))
                        for d in ['s','y','x'] )
    return dims

def x2rho(v, grid):
    """ Interpolate from any grid to rho grid
    """
    dims = _get_spatial_dims(v)
    vout = v.copy()
    if dims['x'] == 'x_u':
        vout = grid.interp(vout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
    return vout

def x2u(v, grid):
    """ Interpolate from any grid to u grid
    """
    dims = _get_spatial_dims(v)
    vout = v.copy()
    if dims['x'] == 'x_rho':
        vout = grid.interp(vout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
    return vout

def x2v(v, grid):
    """ Interpolate from any grid to u grid
    """
    dims = _get_spatial_dims(v)
    vout = v.copy()
    if dims['x'] == 'x_rho':
        vout = grid.interp(vout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
    return vout

def x2x(v, grid, target):
    if target is 'rho':
        return x2rho(v, grid)
    elif target is 'u':
        return x2u(v, grid)
    elif target is 'v':
        return x2v(v, grid)
    
# ------------------------------------------------------------------------------

def get_z(run, zeta=None, h=None, vgrid='r', 
          hgrid=None, vtransform=None):
    ''' Compute vertical coordinates
        Spatial dimensions are placed last, in the order: s_rho/s_w, y, x
        
        Parameters
        ----------
        run: crocosi.gridop.run
            Simulation output object
        zeta: xarray.DataArray, optional
            Sea level data, default to 0 if not provided
            If you use slices, make sure singleton dimensions are kept, i.e do:
                zeta.isel(x_rho=[i])
            and not :
                zeta.isel(x_rho=i)
        h: xarray.DataArray, optional
            Water depth, searche depth in grid if not provided
        vgrid: str, optional
            Vertical grid, 'r'/'rho' or 'w'. Default is 'rho'
        hgrid: str, optional
            Any horizontal grid: 'r'/'rho', 'u', 'v'. Default is 'rho'
        vtransform: int, str, optional
            croco vertical transform employed in the simulation.
            1="old": z = z0 + (1+z0/_h) * _zeta  with  z0 = hc*sc + (_h-hc)*cs
            2="new": z = z0 * (_zeta + _h) + _zeta  with  z0 = (hc * sc + _h * cs) / (hc + _h)
    '''

    grid, xgrid = run['grid'], run['xgrid']

    _h = grid.h if h is None else h
    _zeta = 0*grid.h if zeta is None else zeta

    # switch horizontal grid if needed
    if hgrid in ['u','v']:
        _h = x2x(_h, xgrid, hgrid)
        _zeta = x2x(_zeta, xgrid, hgrid)
    
    # align datasets (zeta may contain a slice along one dimension for example)
    _h, _zeta  = xr.align(_h, _zeta, join='inner')
    
    # determine what kind of vertical corrdinate we are using (NEW_S_COORD)
    if vtransform is None:
        vtransform = grid.Vtransform.values
    else:
        if isinstance(vtransform, str):
            if vtransform.lower()=="old":
                vtransform = 1
            elif vtransform.lower()=="new":
                vtransform = 2
            else:
                raise ValueError("unable to understand what is vtransform")

    if vgrid in ['r', 'rho']:
        vgrid = 'rho'
        sc = grid['sc_r']
        cs = grid['Cs_r']
    else:
        sc = grid['sc_'+vgrid]
        cs = grid['Cs_'+vgrid]

    hc = run['Hc']
    if vtransform == 1:
        z0 = hc*sc + (_h-hc)*cs
        z = z0 + (1+z0/_h) * _zeta
    else:
        z0 = (hc * sc + _h * cs) / (hc + _h)
        z = z0 * (_zeta + _h) + _zeta
    
    # reorder spatial dimensions and place them last
    sdims = list(_get_spatial_dims(z).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in z.dims if d not in sdims) + sdims
    z = z.transpose(*reordered_dims)
    
    return z.rename('z_'+vgrid)

# ------------------------------------------------------------------------------

def w2rho(data, grid, z_w=None, z_r=None):
    """ Linearly interpolates from w vertical grid to rho one.

    Parameters
    ----------
    data: xarray.DataArray
       Array to be interpolated. May contain `z_w` and `z_r` has coordinates.
    grid: xgcm.Grid
        Grid object containing grid information required for the interpolation
    z_w: xarray.DataArray, optional
        Vertical grid at cell faces (w). If not provided will search for variable `z_w` in `data`.
    z_r: xarray.DataArray, optional
        Vertical grid at cell centers (rho). If not provided will search for variable `z_r` in `data`.

    Returns
    -------
    xarray.DataArray
        Interpolated array

    Notes
    -----
    The linear interpolation leverages `xgcm`_ library.

    .. _xgcm:
       https://xgcm.readthedocs.io/en/latest/
    This routine should give the same result as xgcm.interp with updated metrics
    """
    
    if z_w is None and "z_w" in data:
        z_w = data.z_w
        
    if z_r is not None:
        dzr = grid.diff(z_w, "s") #.diff("s_w")
        idn, iup = slice(0,-1), slice(1,None)
        rna = {"s_w":"s_rho"}
        z_w, data = z_w.drop("s_w"), data.drop("s_w")
        w1 = (z_w.isel(s_w=iup).rename(rna) - z_r)/dzr
        w2 = (z_r - z_w.isel(s_w=idn).rename(rna))/dzr
        w_i = (w1*data.isel(s_w=idn).rename(rna) + w2*data.isel(s_w=iup).rename(rna))
        w_i = w_i.assign_coords(z_rho=z_r)
    else:
        w_i = grid.interp(data, "s")
        if z_w is not None:
            w_i = w_i.assign_coords(z_rho=grid.interp(z_w, "s"))
        
    return w_i

def interp2z_3d(z0, z, v, b_extrap=2, t_extrap=2):
    """
    b_extrap, t_extrap:
        0 set to NaN
        1 set to nearest neighboor
        2 linear extrapolation
    """
    import crocosi.fast_interp3D as fi  # OpenMP accelerated C based interpolator
    # check v and z have identical shape
    assert v.ndim==z.ndim
    # add dimensions if necessary
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
    return fi.interp(z0.astype('float64'), lz.astype('float64'), lv.astype('float64'), 
                     b_extrap, t_extrap).squeeze()

def interp2z(z0, z, v, b_extrap, t_extrap):
    ''' interpolate vertically
    '''
    # check v and z have identical shape
    assert v.ndim==z.ndim
    # test if temporal dimension is present
    if v.ndim == 4:
        vi = [interp2z_3d(z0, z[...,t], v[...,t], b_extrap, t_extrap)[...,None] 
                  for t in range(v.shape[-1])]
        return np.concatenate(vi, axis=0) # (50, 722, 258, 1)
        #return v*0 + v.shape[3]
    else:
        return interp2z_3d(z0, z, v, b_extrap, t_extrap)

# ------------------------------------------------------------------------------
    
def get_N2(run, rho, z, g=g_default):
    """ Compute square buoyancy frequency N2 
    ... doc to be improved
    """
    grid = run['xgrid']
    N2 = -g/run['rho0'] * grid.diff(rho, 's', boundary='fill', fill_value=np.NaN) \
            / grid.diff(z, 's', boundary='fill', fill_value=np.NaN)
    # cannot find a solution with xgcm, weird
    N2 = N2.fillna(N2.shift(s_w=-1))
    N2 = N2.fillna(N2.shift(s_w=1))
    return N2

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

def get_p(grid, rho, zw, zr=None, g=g_default):
    """ Compute (not reduced) pressure by integration from the surface, 
    taking rho at rho points and giving results on w points (z grid)
    with p=0 at the surface. If zr is not None, compute result on rho points
    
    Parameters
    ----------
    ...
    
    
    """
    if zr is None:
        dz = grid.diff(zw, "s")
        p = (grid.cumsum((rho*dz).sortby("s_rho",ascending=False), 
                         "s", to="outer", boundary="fill")
             .sortby("s_w", ascending=False)
            )
    else:
        # it is assumed that all fields are from bottom to surface
        rna = {"s_w":"s_rho"}
        dpup = (zr - zw.isel(s_w=slice(0,-1)).drop("s_w").rename(rna))*rho
        dpdn = (zw.isel(s_w=slice(1,None)).drop("s_w").rename(rna) - zr)*rho
        p = ((dpup.shift(s_rho=-1, fill_value=0) + dpdn)
             .sortby(rho.s_rho, ascending=False)                
             .cumsum("s_rho")
             .sortby(rho.s_rho, ascending=True)
             .assign_coords(z_r=zr)
            )
    return _g*p.rename("p")

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
