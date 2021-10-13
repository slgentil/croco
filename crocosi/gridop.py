#import sys # NL: imported but not used
import xarray as xr
import numpy as np
from collections import OrderedDict
from scipy.sparse import lil_matrix

from .postp import grav

### Default values and parameters
_z_dim_database = ['z', 's_rho', 's_w']
_z_coord_database = {"s_rho": ["z_r", "z_rho"], 
                     "s_w": ["z_w"], 
                     "z": ["z"]
                    }

# ---------------------- horizontal grid manipulations -------------------------

def _get_spatial_dims(v):
    """ Return an ordered dict of spatial dimensions in the s/z, y, x order
    """
    dims = OrderedDict( (d, next((x for x in v.dims if x[0]==d), None))
                        for d in ['s','y','x'] )
    return dims

def _get_spatial_coords(v):
    """ Return an ordered dict of spatial dimensions in the s/z, y, x order
    """
    coords = OrderedDict( (d, next((x for x in v.coords if d in x), None))
                        for d in ['z','lat','y','lon','x'] )
    #coords = {k: coords[k] for k in coords.keys() if coords[k]}
    return coords
        
def x2rho(v, grid, run=None):
    """ Interpolate from any grid to rho grid
    """
    dims = _get_spatial_dims(v)
    coords = _get_spatial_coords(v)
    vout = v.copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u':
        vout = grid.interp(vout, 'xi')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'eta')
    if dims['s'] == 's_w':
        vout = grid.interp(vout, 's')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 's')
    # assign coordinates
    if run:
        if not run.grid_regular:
            vout = vout.assign_coords(coords={'nav_lon_rho':run['grid'].nav_lon_rho})
            vout = vout.assign_coords(coords={'nav_lat_rho':run['grid'].nav_lat_rho})
            vout = vout.assign_coords(coords={'z_r':zout})
                            
    return vout

def x2u(v, grid, run=None):
    """ Interpolate from any grid to u grid
    """
    dims = _get_spatial_dims(v)
    coords = _get_spatial_coords(v)
    vout = v.copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_rho':
        vout = grid.interp(vout, 'xi')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'eta')
    if dims['s'] == 's_w':
        vout = grid.interp(vout, 's')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 's')
    if run:
        if not run.grid_regular:
            vout = vout.assign_coords(coords={'nav_lon_u':run['grid'].nav_lon_u})
            vout = vout.assign_coords(coords={'nav_lat_u':run['grid'].nav_lat_u})
            vout = vout.assign_coords(coords={'z_u':zout})
    return vout

def x2v(v, grid, run=None):
    """ Interpolate from any grid to v grid
    """
    dims = _get_spatial_dims(v)
    coords = _get_spatial_coords(v)
    vout = v.copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u':
        vout = grid.interp(vout, 'xi')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'xi')
    if dims['y'] == 'y_rho':
        vout = grid.interp(vout, 'eta')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'eta')
    if dims['s'] == 's_w':
        vout = grid.interp(vout, 's')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 's')
    if run:
        if not run.grid_regular:
            vout = vout.assign_coords(coords={'nav_lon_v':run['grid'].nav_lon_v})
            vout = vout.assign_coords(coords={'nav_lat_v':run['grid'].nav_lat_v})
            vout = vout.assign_coords(coords={'z_v':zout})
    return vout

def x2w(v, grid, run=None):
    """ Interpolate from any grid to w grid
    """
    dims = _get_spatial_dims(v)
    coords = _get_spatial_coords(v)
    vout = v.copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u':
        vout = grid.interp(vout, 'xi')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'xi')
    if dims['y'] == 'y_v':
        vout = grid.interp(vout, 'eta')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'eta')
    if dims['s'] == 's_rho':
        vout = grid.interp(vout, 's')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 's')
    if run:
        if not run.grid_regular:
            vout = vout.assign_coords(coords={'nav_lon_rho':run['grid'].nav_lon_rho})
            vout = vout.assign_coords(coords={'nav_lat_rho':run['grid'].nav_lat_rho})
            vout = vout.assign_coords(coords={'z_w':zout})
    return vout

def x2psi(v, grid, run=None):
    """ Interpolate from any grid to psi grid
    """
    dims = _get_spatial_dims(v)
    coords = _get_spatial_coords(v)
    vout = v.copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_rho':
        vout = grid.interp(vout, 'xi')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'xi')
    if dims['y'] == 'y_rho':
        vout = grid.interp(vout, 'eta')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 'eta')
    if dims['s'] == 's_w':
        vout = grid.interp(vout, 's')
        if run:
            if not run.grid_regular:
                zout = grid.interp(zout, 's')
    if run:
        if not run.grid_regular:
            vout = vout.assign_coords(coords={'nav_lon_f':run['grid'].nav_lon_f})
            vout = vout.assign_coords(coords={'nav_lat_f':run['grid'].nav_lat_f})
            vout = vout.assign_coords(coords={'z_f':zout})
    return vout

def x2x(v, grid, target, run=None):
    if target == 'rho':
        return x2rho(v, grid, run=run)
    elif target == 'u':
        return x2u(v, grid, run=run)
    elif target == 'v':
        return x2v(v, grid, run=run)
    elif target == 'w':
        return x2w(v, grid, run=run)
    elif target == 'psi':
        return x2psi(v, grid, run=run)

def _get_grid_point(var):
    dims = var.dims
    if "x_u" in dims:
        if "y_rho" in dims:
            return 'u'
        else:
            return 'psi'
    elif "y_v" in dims:
        return 'v'
    else:
        if 's_rho' in dims:
            return 'rho'
        else:
            return 'w'
        
        
def poisson_matrix(pm,pn):    
    """
    Initialize the elliptic equation matrix (d_xx + d_yy)
    Input :
        - pm : (ndarray) 1/dx coefficents
        - pn : (ndarray) 1/dy coefficents
    Output:
        return a sparse matrix of the laplacian operator
    """
    # elliptic equation matrix:  d_xx + d_yy
    [nx,ny] = pm.shape
    
    ndim = ny*nx;
    i_s  = ny;
    js  = 1;
    A=lil_matrix((ndim,ndim))
    ############################
    
    for i in range(nx): 
        for j in range(ny): 
                idx = i*ny + j;
                diag = 0.;
                if j>0:
                    dy2i = 0.5*(pn[i,j]+pn[i,j-1])*pn[i,j]
                    A[idx,idx-js] = dy2i;
                    diag -= dy2i;
                if i>0:
                    dx2i = 0.5*(pm[i,j]+pm[i-1,j])*pm[i,j]
                    A[idx,idx-i_s] = dx2i;
                    diag -= dx2i;
                if i<nx-1:
                    dx2i = 0.5*(pm[i,j]+pm[i+1,j])*pm[i,j]
                    A[idx,idx+i_s] = dx2i;
                    diag -= dx2i;
                if j<ny-1:
                    dy2i = 0.5*(pn[i,j]+pn[i,j+1])*pn[i,j]
                    A[idx,idx+js] = dy2i;
                    diag -= dy2i;
                A[idx,idx] = diag
    
    return A
# -----------------------  grid interpolation ------------------------

def find_nearest_above(my_array, target, axis=0):
    diff = target - my_array
    diff = diff.where(diff>0,np.inf)
    return xr.DataArray(diff.argmin(axis=axis))

def get_slice(run, var, z, longitude=None, latitude=None, depth=None):
        """
        #
        #
        # This function interpolate a 3D variable on a slice at a constant depth or
        # constant longitude or constant latitude. 
        # This function uses interp2z over z axe and calculate linear interpolation over x and y axes
        #
        # On Input:
        #
        #    ds      dataset to find the grid
        #    var     (dataArray) Variable to process (3D matrix).
        #    z       (dataArray) Depths at the same point than var (3D matrix).
        #    longitude   (scalar) longitude of the slice (scalar meters, negative).
        #    latitude    (scalar) latitude of the slice (scalar meters, negative).
        #    depth       (scalar) depth of the slice (scalar meters, negative).
        #
        # On Output:
        #
        #    vnew    (dataArray) Horizontal slice
        #
        #
        """
        if longitude is None and latitude is None and depth is None:
            "Longitude or latitude or depth must be defined"
            return None
        
        if z.shape != var.shape:
            print('slice: var and z shapes are different')
            return

        # Find dimensions of the variable
        dims = _get_spatial_dims(var)
        grid_point = _get_grid_point(var)
        N = len(var[dims['s']])

        # Find horizontal coordinates of the variable
        x = [var.coords[s] for s in var.coords if "lon_" in s or "x_" in s][0]
        y = [var.coords[s] for s in var.coords if "lat_" in s or "y_" in s][0] 


        # Find the indices of the grid points just below the longitude/latitude/depth
        if longitude is not None:
            axe = x.get_axis_num(dims['x'])
            indices = find_nearest_above(x, longitude, axis=axe)
        elif latitude is not None:
            axe = y.get_axis_num(dims['y'])
            indices = find_nearest_above(y, latitude, axis=axe)
            
        # Initializes the 2 slices around the longitude/latitude/depth
        if longitude is not None:
            x1 = x.isel({dims['x']:indices})
            x2 = x.isel({dims['x']:indices+1})
            if not run.grid_regular:
                y1 = y.isel({dims['x']:indices})
                y2 = y.isel({dims['x']:indices+1})
            z1 = z.isel({dims['x']:indices})
            z2 = z.isel({dims['x']:indices+1})
            v1 = var.isel({dims['x']:indices})
            v2 = var.isel({dims['x']:indices+1})
        elif latitude is not None:
            if not run.grid_regular:
                x1 = x.isel({dims['y']:indices})
                x2 = x.isel({dims['y']:indices+1})
            y1 = y.isel({dims['y']:indices})
            y2 = y.isel({dims['y']:indices+1})
            z1 = z.isel({dims['y']:indices})
            z2 = z.isel({dims['y']:indices+1})
            v1 = var.isel({dims['y']:indices})
            v2 = var.isel({dims['y']:indices+1})

        # Do the linear interpolation
        if longitude is not None:
            xdiff = x1 - x2
            if not run.grid_regular:
                ynew =  (((y1 - y2) * longitude + y2 * x1 - y1 * x2) / xdiff)
            znew =  (((z1 - z2) * longitude + z2 * x1 - z1 * x2) / xdiff).fillna(0)
            vnew =  (((v1 - v2) * longitude + v2 * x1 - v1 * x2) / xdiff)
        elif latitude is not None:
            ydiff = y1 - y2
            if not run.grid_regular:
                xnew =  (((x1 - x2) * latitude + x2 * y1 - x1 * y2) / ydiff)
            znew =  (((z1 - z2) * latitude + z2 * y1 - z1 * y2) / ydiff).fillna(0)
            vnew =  (((v1 - v2) * latitude + v2 * y1 - v1 * y2) / ydiff)
        elif depth is not None:
            try:
                zmask = x2x(run['grid'].mask_rho.where(run['grid'].h>abs(depth),np.nan) ,
                        run.xgrid,grid_point)    
            except:
                zmask = var.isel({dims['s']:0}) * 0. +1
            swapdim ='s_rho' if dims['s']=='s_w' else 's_w'
            zt = (var.isel({dims['s']:[0]}).fillna(0.)*0. + depth).rename({dims['s']:swapdim})
            vnew = zmask*interp2z(zt, z, var)

        # Add the coordinates to dataArray
        if longitude is not None:
            if not run.grid_regular:
                ynew = ynew.expand_dims({dims['s']: N})
            else:
                ynew = y
            vnew = vnew.assign_coords(coords={"z":znew})
            vnew = vnew.assign_coords(coords={y.name:ynew})

        elif latitude is not None:
            if not run.grid_regular:
                xnew = xnew.expand_dims({dims['s']: N})
            else:
                xnew = x
            vnew = vnew.assign_coords(coords={"z":znew})
            vnew = vnew.assign_coords(coords={x.name:xnew})

        elif depth is not None:
            vnew = vnew.assign_coords(coords={y.name:y})
            vnew = vnew.assign_coords(coords={x.name:x})

        return vnew.squeeze()

def get_slices(run, var, z, longitude=None, latitude=None, depth=None):
    """
    #
    #
    # This function interpolate a 3D variable on slices at constant depths/longitude/latitude
    # This function use xcgm transform method and needs xgcm.Grid to be defined over the 3 axes. 
    # !!! For now, it works only with curvilinear coordinates !!!
    #
    # On Input:
    #
    #    run      dataset to find the grid
    #    var     (dataArray) Variable to process (3D matrix).
    #    z       (dataArray) Depths at the same point than var (3D matrix).
    #    longitude   (scalar,list or ndarray) longitude of the slice (scalar meters, negative).
    #    latitude    (scalar,list or ndarray) latitude of the slice (scalar meters, negative).
    #    depth       (scalar,list or ndarray) depth of the slice (scalar meters, negative).
    #
    # On Output:
    #
    #    vnew    (dataArray) Horizontal slice
    #
    #
    """
    from matplotlib.cbook import flatten

    if longitude is None and latitude is None and depth is None:
        "Longitude or latitude or depth must be defined"
        return None

    # check typ of longitude/latitude/depth
    longitude = longitude.tolist() if isinstance(longitude,np.ndarray) else longitude
    longitude = [longitude] if (isinstance(longitude,int) or isinstance(longitude,float)) else longitude

    latitude = latitude.tolist() if isinstance(latitude,np.ndarray) else latitude
    latitude = [latitude] if (isinstance(latitude,int) or isinstance(latitude,float)) else latitude

    depth = depth.tolist() if isinstance(depth,np.ndarray) else depth
    depth = [depth] if (isinstance(depth,int) or isinstance(depth,float)) else depth

     # Find dimensions and coordinates of the variable
    dims = _get_spatial_dims(var)
    coords = _get_spatial_coords(var)

    if longitude is not None:
        axe = 'xi'
        coord_ref = coords['x'] if coords['x'] else coords['lon']
        coord_x = coords['y'] if coords['y'] else coords['lat']
        coord_y = coords['z']
        slices_values = longitude
    elif latitude is not None:
        axe = 'eta'
        coord_ref = coords['y'] if coords['y'] else coords['lat']
        coord_x = coords['x'] if coords['x'] else coords['lon']
        coord_y = coords['z']
        slices_values = latitude
    else:
        axe = 's'
        coord_ref = coords['z']
        coord_x = coords['x'] if coords['x'] else coords['lon']
        coord_y = coords['y'] if coords['y'] else coords['lat']
        slices_values = depth

    # Recursively loop over time if needed
    if len(var.dims) == 4:
        vnew = [get_slices(run, var.isel(time=t), z.isel(time=t), 
                      longitude=longitude, latitude=latitude, depth=depth) 
                      for t in range(len(var.time))]
        vnew = xr.concat(vnew, dim='time')
    else:
        vnew = run.xgrid.transform(var, axe, slices_values, 
                               target_data=var[coord_ref])

    # Do the linear interpolation
    if not depth:
        x = run.xgrid.transform(var[coord_x], axe, slices_values, 
                                   target_data=var[coord_ref])\
                     .expand_dims({dims['s']: len(var[dims['s']])})
        y = run.xgrid.transform(var[coord_y], axe, slices_values, 
                                   target_data=var[coord_ref])

        # Add the coordinates to dataArray
        vnew = vnew.assign_coords(coords={coord_x:x})
        vnew = vnew.assign_coords(coords={coord_y:y})
    else:
        # Add the coordinates to dataArray
        vnew = vnew.assign_coords(coords={coord_x:var[coord_x]})
        vnew = vnew.assign_coords(coords={coord_y:var[coord_y]})

    return vnew.squeeze().unify_chunks()
    
def hinterp(ds,var,coords=None):
    """ This is needs proper documentation:
    https://numpydoc.readthedocs.io/en/latest/format.html
    """
    import pyinterp
    #create Tree object
    mesh = pyinterp.RTree()

    #L = ds.dims['x_r'] # NL: assigned but not used
    #M = ds.dims['y_r']
    #N = ds.dims['s_r']
    z_r = get_z(ds)
    #coords = np.array([coords])

    # where I know the values
    z_r = get_z(ds)
    vslice = []
    #lon_r = np.tile(ds.lon_r.values,ds.dims['s_r']).reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
    lon_r = (np.tile(ds.lon_r.values,(ds.dims['s_r'],1,1))
            .reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
            )
    lat_r = (np.tile(ds.lat_r.values,(ds.dims['s_r'],1,1))
             .reshape(ds.dims['s_r'], ds.dims['y_r'], ds.dims['x_r'])
            )
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
    #mask=np.where(ds.mask_rho.values==1.,) # NL: assigned but not used
    vslice[int(ds.mask_rho.values)] = float("nan")

    vslice = xr.DataArray(np.asarray(vslice.reshape(ds.dims['y_r'], ds.dims['x_r'])),
                          dims=('x_r','y_r'))
    yslice = ds.lat_r
    xslice = ds.lon_r

    return[xslice,yslice,vslice]

# --------------------------- depth coordinate ---------------------------------

def get_cs_sigma(vtransform, theta_s, theta_b, N, vgrid, xarray):
    """ get croco vertical grid auxilliary variables:
    vertical stretching and sigma
    
    https://www.myroms.org/wiki/Vertical_S-coordinate
    """
    if xarray:
        if vgrid=='r':
            sigma = xr.DataArray(np.arange(-1.+1./(N+1),0., 1/(N+1)), dims='s_rho')
        elif vgrid=='w':
            sigma = xr.DataArray(np.arange(-1.,0., 1/(N+1)), dims='s_w')
    else:
        if vgrid=='r':
            sigma = np.arange(-1.+1./(N+1),0., 1/(N+1))
        elif vgrid=='w':
            sigma = np.arange(-1.,0., 1/(N+1))
    #
    if vtransform == 1:
        cs = (1-theta_b)*np.sinh(theta_s*sigma)/np.sinh(theta_s) \
             + theta_b*0.5 \
               *(np.tanh((sigma+0.5)*theta_s)/np.tanh(0.5*theta_s)-1.)
    elif vtransform == 2:
        #_sc = (sigma-N)/N
        _sc = sigma
        if theta_s>0:
            csf = (1-np.cosh(theta_s*sigma)) / (np.cosh(theta_s)-1.)
        else:
            csf = -sigma**2
        if theta_b>0:
            cs = (np.exp(theta_s*csf)-1.) / (1-np.exp(-theta_b))
        else:
            cs = csf
    return cs, sigma
    
def get_z_core(zeta, h, vtransform, hc, 
               sigma=None, cs=None,
               theta_s=None, theta_b=None, N=None, vgrid=None,
               xarray=True):
    """ Computes z grid
    
    Parameters
    ----------
    zeta: np.ndarray, xr.DataArray
        sea level
    h: np.ndarray, xr.DataArray
        water depth
    vtransform: int
        vertical transform
    hc: float
        thickness parameter
    sc: float, np.ndarray, xr.DataArray
        sigma array (-1,0)
    cs: float, np.ndarray, xr.DataArray
        vertical stretching
        
    https://www.myroms.org/wiki/Vertical_S-coordinate
    """

    if cs is None and sigma is None:
        cs, sigma = get_cs_sigma(vtransform, theta_s, theta_b, N, vgrid, xarray)

    if vtransform == 1:
        z0 = hc*sigma + (h-hc)*cs
        z = z0 + (1+z0/h) * zeta
    elif vtransform == 2:
        z0 = (hc * sigma + h * cs) / (hc + h)
        z = z0 * (zeta + h) + zeta

    return z

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
        vtransform = int(grid.Vtransform)
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
    
    z = get_z_core(_zeta, _h, vtransform, hc, sigma=sc, cs=cs)
    
    # reorder spatial dimensions and place them last
    sdims = list(_get_spatial_dims(z).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in z.dims if d not in sdims) + sdims
    z = z.transpose(*reordered_dims, transpose_coords=True)
    
    return z.rename('z_'+vgrid)

def get_z_dim(ds, z_dim_database=_z_dim_database):
    """ return list of recognized vertical dimensions according to a reference list 
    """
    dims = [d for d in ds.dims if d in z_dim_database]
    dims = dims if dims else None # replace empty list by None
    return dims

def get_z_coord(ds, coords=_z_coord_database, dims=_z_dim_database):
    """ return list of names of z coordinate in an xarray DataSet or DataArray"""
    #if "s_rho" in ds.dims:
        #zname = next((z for z in ['z_r', 'z_rho'] if z in ds.coords), None)
    #elif "s_w" in ds.dims:
        #zname = "z_w" if "z_w" in ds.coords else None
    #elif "z" in ds.dims:
        #zname = "z" if "z" in ds.coords else None
    zname = []
    for dim in dims:
        if dim in ds.dims:
            zname += [z for z in coords[dim] if z in ds.coords]
    if not zname:
        raise ValueError("could not find z coordinate")
    return zname

def get_xgrid_ax_name(xgrid, sdim):
    """ return name of ax in xgrid object that match the list of dimensions 'sdim' """
    return next(( lax.name for lax in xgrid.axes.values() \
                    if all([dim in lax.coords.values() for dim in sdim]) ))

# ---------------------------- vertical interpolation --------------------------
# 
def w2rho(data, grid, z_r=None, z_w=None, s_dims=["s_rho", "s_w"]):
    """ Linearly interpolates from w vertical grid to rho one.

    Parameters
    ----------
    data: xarray.DataArray
       Array to be interpolated. May contain `z_w` and `z_r` as coordinates.
    grid: xgcm.Grid
        Grid object containing grid information required for the interpolation
    z_r: xarray.DataArray, optional 
        Vertical grid at cell centers (rho). If not provided will search for variable `z_r` in `data`.
    z_w: xarray.DataArray, optional
        Vertical grid at cell faces (w). If not provided will search for variable `z_w` in `data`.
    s_dims: list of two str, optional
        names of vertical dimensions in, resp. , z_r and z_w (default: s_rho, s_w)

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
    
    if z_w is None:
        z_w = data.z_w
        
    if z_r is None:
        z_r = data.z_r
    s_rho, s_w = s_dims
    xgrid_s = get_xgrid_ax_name(grid,s_dims)
    dzr = grid.diff(z_w, xgrid_s) #.diff("s_w")
    idn, iup = slice(0,-1), slice(1,None)
    rna = {s_w:s_rho}
    z_w, data = z_w.drop(s_w), data.drop(s_w)
    # TODO use shift instead of isel (if the routine is maintained) 
    w1 = (z_w.isel({s_w:iup}).rename(rna) - z_r)/dzr
    w2 = (z_r - z_w.isel({s_w:idn}).rename(rna))/dzr
    w_i = (w1*data.isel({s_w:idn}).rename(rna) + w2*data.isel({s_w:iup}).rename(rna))

    return w_i.assign_coords(z_rho=z_r)


def _guess_dim_mapping(x, y):
    # guess mapping
    x_shape = {d+1: dsize for d, dsize in zip(range(x.ndim-1),x.shape[1:])}
    y_shape = y.shape
    ndim = y.ndim
    dim_map = np.arange(ndim)
    find_dim = lambda dsize, shape: \
            next((d for d, size in shape.items() if size==dsize), -1)
    for d, dsize in zip(np.arange(1,ndim), y_shape[1:]):
        xd = find_dim(dsize, x_shape)
        xd1 = find_dim(1, x_shape)
        if xd>-1:
            # use first matching element
            dim_map[d] = xd
            del x_shape[xd]
        elif xd1>-1:
            # or default to first singleton dimension
            dim_map[d] = xd1
            del x_shape[xd1]
        else:
            _msg = 'align_broadcast_but_dim0 could not align and broadcast' \
                  +'arrays with shape: {} and {}'.format(x.shape, y_shape)
            raise ValueError(_msg)
    return dim_map

def _align_dims_but_one(x,y, exclude=None):
    """ Align dimensions of x against that of y but one
    
    Parameters
    ----------
    x, y: ndarray
    dims: tuple, optional
        Dimensions that need to be excluded from the alignment
        By default exclude the first dimension
    """
    if exclude:
        # swaps dimensions to first spot
        _x = _align_dims_but_one(x.swapaxes(0,exclude[0]),
                                 y.swapaxes(0,exclude[1]))
        return _x.swapaxes(0, exclude[1])
    # add dimensions of necessary to x
    if x.ndim<y.ndim:
        _x = x.reshape(x.shape + tuple(1 for d in range(y.ndim-x.ndim)))
    else:
        _x = x
    # get mapping or crash
    dim_map = _guess_dim_mapping(_x, y)
    # transpose and broadcast
    _x = np.transpose(_x, dim_map)
    # could also broadcast
    #_x = np.broadcast_to(_x, (_x.shape[0],)+y.shape[1:])
    return _x

def interp2z_np_3d(zt, z, v, b_extrap=2, t_extrap=2):
    """ Interpolates 3D arrays from one vertical grid to another
    
    Parameters
    ----------
    zt: ndarray
        Target vertical grid, may be 1D/2D/3D, ascending order
    z: ndarray
        Initial vertical grid, may be 1D/2D/3D, ascending order
    v: ndarray
        Initial data, may be 1D/2D/3D
    b_extrap: int, optional
        Bottom boundary condition. 1: closest value, 2: linear extrapolation (default), otherwise: NaN
    t_extrap: int, optional
        Top boundary condition. 1: closest value, 2: linear extrapolation (default), otherwise: NaN
        
    Returns
    -------
    vt: ndarray
    
    Notes
    -----
    The vertical dimension must be the first one.
    zt and z must be in ascending order along their vertical dimension.
    zt, z, v must have the same number of dimensions.
    zt may have singleton dimensions along horizontal dimensions.
    z and v must have the same shape.
    If you bring any modification to this code, make sure it does not break
    tutorials/dev/interp2z.ipynb.
    """
    import crocosi.fast_interp3D as fi  # OpenMP accelerated C based interpolator
    # check shapes and number of dimensions
    assert v.shape==z.shape, \
        'mismatch: v_shape={} but z_shape={}'.format(v.shape, z.shape)
    assert v.ndim==zt.ndim, \
        'mismatch: v.ndim={} but zt.ndim={}'.format(v.ndim, zt.ndim)    
    # add dimensions to v and z to be 3D
    if v.ndim == 1:
        _v = v.squeeze()[:,None,None]
        _z = z.squeeze()[:,None,None]
        _zt = zt.squeeze()[:,None,None]
    elif v.ndim == 2:
        _v = v[...,None]
        _z = z[...,None]
        _zt = zt[...,None]
    else:
        _z = z[...]
        _v = v[...]
        _zt = zt[...]
    # broadcast zt along all dimensions but the first one
    _zt = np.broadcast_to(_zt, (_zt.shape[0],)+_v.shape[1:])
    #
    out = fi.interp(_zt.astype('float64'), 
                    _z.astype('float64'),
                    _v.astype('float64'), 
                    b_extrap, t_extrap)
    if v.ndim==1:
        return out[:,0,0]
    elif v.ndim==2:
        return out[...,0]
    else:
        return out

def interp2z_np(zt, z, v, zdim=None, zdimt=None, **kwargs):
    ''' Interpolates arrays from one vertical grid to another
    
    Parameters
    ----------
    zt:  ndarray
        Target vertical grid, may be 1D/2D/3D but NOT 4D, ascending order
    z:   ndarray
        Initial vertical grid, may be 1D/2D/3D/4D, ascending order
    v: ndarray
        Initial data, may be 1D/2D/3D/4D
    zdim: tuple
        Position and size of the vertical dimension for v and z 
        as a tuple: (pos, size)
    zdimt: tuple
        Position and size of the vertical dimension for zt 
        as a tuple: (pos, size)
    **kwargs: passed to interp2z_np_3D

    Returns
    -------
    vt: ndarray
    Preserves v and z shapes, except for the vertical dimension whose size
    will match that of zt

    Notes
    -----
    v and z must have the same shape.
    zt and z must be in ascending order along their vertical dimension.
    If you bring any modification to this code, make sure it does not break
    tutorials/dev/interp2z.ipynb.
    '''
    # check v and z have identical shape
    assert v.shape==z.shape, \
            'v and z have different shapes: {}, {}'.format(v.shape, z.shape)
    assert v.shape[zdim[0]]==zdim[1], \
            'v shape and zdim are not consistent: {}, {}'.format(v.shape, zdim)
    assert zt.shape[zdimt[0]]==zdimt[1], \
            'zt shape and zdimt are not consistent: {}, {}'.format(zt.shape, zdimt)
    # zdim and zdimt are not optional parameters in fact:
    assert zdim, 'zdim was not provided'
    assert zdimt, 'zdimt was not provided'
    # align zt dimensions against that of v
    _zt = _align_dims_but_one(zt, v, exclude=(zdimt[0], zdim[0]))
    if v.ndim == 4:
        # move vertical dimension in second position
        _z = z.swapaxes(1,zdim[0])
        _v = v.swapaxes(1,zdim[0])
        _zt = _zt.swapaxes(1,zdim[0])
        #
        vi = [interp2z_np_3d(_zt[t,...], _z[t,...], _v[t,...], **kwargs)[None,...]
                      for t in range(_v.shape[0])]
        vi = np.concatenate(vi, axis=0)
        #
        return vi.swapaxes(1, zdim[0])
    else:
        _z = z.swapaxes(0, zdim[0])
        _v = v.swapaxes(0, zdim[0])
        _zt = _zt.swapaxes(0, zdim[0])
        #
        vi = interp2z_np_3d(_zt, _z, _v, **kwargs)
        #
        return vi.swapaxes(0, zdim[0])

def interp2z(zt, z, v, zt_dim=None, z_dim=None, 
             override_dims=False, output_dims=None,
             **kwargs):
    ''' Interpolates xarrays from one vertical grid to another
    
    Parameters
    ----------
    zt:  xarray.DataArray
        Target vertical grid, ascending order
    z:   xarray.DataArray
        Initial vertical grid, ascending order
    v:   xarray.DataArray
        Initial data
    zt_dim: str, optional
        Name of the target vertical dimension 
    z_dim: str, optional
        Name of the initial vertical dimension
    override_dims: boolean, optional
        If True, allow zt and z to have same vertical dimension name but different
        dimension values
    output_dims: list, optional
        This list will adjust the order of output dimensions
    **kwargs: passed to interp2z_np_3D

    Returns
    -------
    vt: xarray.DataArray
    Preserves v dimension order

    Notes
    -----
    When zt_dim or z_dim are not provided, we search through a database
    of known vertical dimension names.
    zt and z must be in ascending order along their vertical dimension.
    If you bring any modification to this code, make sure it does not break
    tutorials/dev/interp2z.ipynb.
    '''
    # search for vertical dimensions names
    _zdims_database = ['z', 's_rho','s_w'] # erase 'z' eventually
    if zt_dim is None:
        zt_dim = next((d for d in zt.dims if d in _zdims_database), None)
        assert zt_dim, 'Could not find a vertical dimension for zt'
    if z_dim is None:
        z_dim = next((d for d in v.dims if d in _zdims_database), None)
        assert z_dim, 'Could not find a vertical dimension for z'
    # test if z_dim and zt_dim are equal but refer to dimensions with different
    # values
    if z_dim==zt_dim and not (v[z_dim] == zt[zt_dim]).all():
        if override_dims:
            _zt_dim = 'zt_swap'
            _zt = zt.rename({zt_dim: _zt_dim})
        else:
            _msg = 'z and zt dimensions have same name but different values.' \
                  +' This is not allowed unless override_dims=True'
            raise ValueError(_msg)
    else:
        _zt_dim = zt_dim
        _zt = zt
    # broadcast variables against each other
    _v, _z, _zt = xr.broadcast(v, z, _zt, exclude=[_zt_dim, z_dim])
    #
    pos_size = lambda v, dim: (v._get_axis_num(dim), v[dim].size)
    z_pos, z_size = pos_size(_v, z_dim)
    zt_pos, zt_size = pos_size(_zt, _zt_dim)
    #
    # adjust chunks in the vertical dimension
    _zt = _zt.chunk({_zt_dim: -1})
    _z = _z.chunk({z_dim: -1})
    _v = _v.chunk({z_dim: -1})
    # provide core dimensions if necessary
    ufunc_kwargs = {}
    if _zt_dim!=z_dim or _zt[_zt_dim].size!=_z[z_dim].size:
        # core dimensions are moved to the last position
        z_pos, zt_pos = -1, -1
        ufunc_kwargs.update(**{'input_core_dims': [[_zt_dim], [z_dim], [z_dim]],
                               'output_core_dims': [[_zt_dim]]})
    #
    interp_kwargs = {'zdim': (z_pos, z_size), 'zdimt': (zt_pos, zt_size)}
    interp_kwargs.update(kwargs)
    vout = xr.apply_ufunc(interp2z_np, _zt, _z, _v,
                          kwargs=interp_kwargs,
                          dask='parallelized',
                          output_dtypes=[np.float64],
                          **ufunc_kwargs)
    # reorder dimensions to match that of v (modulo change along vertical)
    if not output_dims:
        output_dims = [_zt_dim if d==z_dim else d for d in list(_zt.dims)]
    vout = vout.transpose(*output_dims)
    # swap vertical dim name back to original value
    if _zt_dim is not zt_dim:
        vout = vout.rename({_zt_dim: zt_dim})
    return vout

# ----------------------------- physics & dynamics -----------------------------
    
def get_N2(run, rho, z, g=grav):
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

def get_oldp(grid, rho, zw, zr=None, g=grav):
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
    return g*p.rename("p")

# !!! code below same as croco
def get_p(grid, rho, z_w, z_r=None, g=None, rho0=None):
    """ 
    Compute (not reduced) pressure by integration from the surface, 
    taking rho at rho points and giving results on rho points (z grid).
    
    Parameters
    ----------
    grid : xgcm grid
    rho  : Density (DataArray)
    z_w  : depth at vertical w points (DataArray)
    z_r  : depth at vertical rho points (DataArray)
    g    : acceleration of gravity (float)
    rho0 : mean density (float)
    
    """
    # useful parameters
    eps = 1.0e-10
    if g is None: g=9.81
    if rho0 is None: rho0=1000.
    GRho=g/rho0
    HalfGRho=0.5*GRho
    N = rho.sizes['s_rho']
    OneFifth = 1.0/5.0
    OneTwelfth = 1.0/12.0
    
    # compute z_r if None
    if z_r is None: z_r = grid.interp(z_w, 's')
        
    # dz and drho on w levels
    dR = grid.diff(rho,'s', boundary='extrapolate').rename('dR')
    dZ = grid.diff(z_r,'s', boundary='extrapolate').rename('dZ')

    # modified dz and dr on w levels
    dZi = 2. * ( dZ * dZ.shift(s_w=1,fill_value=0) 
            / (dZ + dZ.shift(s_w=1,fill_value=0)) )
    dZi = xr.concat([dZ.isel(s_w=0), dZi.isel(s_w=slice(1,None))], dim='s_w')
    dRi = 2. * ( dR * dR.shift(s_w=1,fill_value=0) 
            / (dR + dR.shift(s_w=1,fill_value=0)) )
    dRi = xr.concat([dR.isel(s_w=0), dRi.isel(s_w=slice(1,None))], dim='s_w')

    # Pressure at the surface on rho level
    Pn = (g*z_w.isel(s_w=-1) + GRho*( rho.isel(s_rho=-1) 
                                   + 0.5*(rho.isel(s_rho=-1)-rho.isel(s_rho=-2)) 
                                   * (z_w.isel(s_w=-1)-z_r.isel(s_rho=-1)) 
                                   / (z_r.isel(s_rho=-1)-z_r.isel(s_rho=-2))
                                  ) * (z_w.isel(s_w=-1)-z_r.isel(s_rho=-1))
         )
    
    # Pressure of each slice
    rhoz = (rho + rho.shift(s_rho=-1))*(z_r.shift(s_rho=-1) - z_r)
    dRz = ((grid.diff(dRi,'s'))*(z_r.shift(s_rho=-1) - z_r
                                 -OneTwelfth*(grid.interp(dZi, 's'))) )
    dZr = ((grid.diff(dZi,'s'))*(rho.shift(s_rho=-1) - rho
                                 -OneTwelfth*(grid.interp(dRi,'s'))) )
    Pk = (HalfGRho*( rhoz - OneFifth*(dRz - dZr)))

    # replace pressure at the surface
    Pk = xr.concat([Pk.isel(s_rho=slice(0,-1)), Pn], dim='s_rho')

    # integrate from the surface to bottom
    P = (Pk
            .sortby(Pk.s_rho, ascending=False)                
            .cumsum("s_rho")
            .sortby(Pk.s_rho, ascending=True)
            .assign_coords(z_r=z_r)
        )

    return P.rename('P')

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
    xi = x2rho(xi_u+xi_v, ds) # NL: psi2rho replaced by x2rho, not checked

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

# ----------------------------- grid rechunk -----------------------------

def ds_hor_chunk(ds, keep_complete_axe=None, wanted_chunk=100):
    """
    Rechunk Dataset or DataArray such as each partition size is about 100Mb
    Input:
        - ds : (Dataset or DataArray) object to rechunk
        - keep_complete_axe : (character) Horizontal axe to keep with no chunk (x or y)
        - wanted_chunk : (integer) size of each partition in Mb
    Output:
        - object rechunked
    """
    
    #check input parameters
    if not isinstance(ds, (xr.Dataset,xr.DataArray)):
        print('argument must be a xarray.DataArray or xarray.Dataset')
        return
    if keep_complete_axe and keep_complete_axe!= 'x' and keep_complete_axe!='y':
        print('keep_complete_axe must equal x or y')
        return
    
    # get horizontal dimensions of the Dataset/DataArray
    chunks = {}
    hor_dim_names = {}
    if isinstance(ds,xr.Dataset):
        s_dim = max([ds.dims[s] for s in ds.dims.keys() if 's_' in s])
        for key,value in ds.dims.items():
            if 'x' in key or 'y' in key:
                hor_dim_names[key] = value
    else:
        s_dim = max([len(ds[s]) for s in ds.dims if 's_' in s])
        for key in ds.dims:
            if 'x' in key or 'y' in key:
                hor_dim_names[key] = len(ds[key])
                           
    if keep_complete_axe:
        # get the maximum length of the dimensions along the argument axe
        # set the chunks of those dimensions to -1 (no chunk)
        h_dim=[]
        for s in hor_dim_names.keys():
            if keep_complete_axe in s:
                h_dim.append(hor_dim_names[s])
                hor_dim_names[s] = -1
        h_dim = max(h_dim) 
        # compute the chunk on the dimensions along the other horizontal axe
        size_chunk = int(np.ceil(wanted_chunk*1.e6 / 4. / s_dim / h_dim))
    else : 
        # compute the chunk on all horizontal dimensions
        size_chunk = int(np.ceil(np.sqrt(wanted_chunk*1.e6 / 4. / s_dim )))
    
    # Initialize the dctionnary of chunks
    for s in hor_dim_names.keys():
        chunks[s]=min([size_chunk,hor_dim_names[s]])
        
    return ds.chunk(chunks)