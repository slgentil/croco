#!/usr/bin/env python
# coding: utf-8
""" cswbib.py
contains routines and wrapper usefull for DiagDynCSW analysis. Also some different version of routines that are present in gridop.py and should be merged together in the future
Copy pasted from DiagDynCSW_proj_v3.ipynb (2020/02/05), probably modified later on.
Requirements:
 - grid (xgcm-grid)
"""

import numpy as np
import xarray as xr
import crocosi.postp as pp
import crocosi.gridop as gp
from crocosi.vmodes import get_vmodes

_g = 9.81

################################################################################
##################### Miscellaneous routines ###################################
################################################################################
def N2Profile(run, strat, z):
    """
    Method to compute the N2 profile : 
    N.B.: s_w must be correct, i.e. with values -1 at bottom and 0 at top
    """   
    grid = run.ds['grid'].attrs['xgcm-Grid']
    rho0 = run.params_input["rho0"]
    N2 = -_g/rho0 * grid.diff(strat,'s') / grid.diff(z,'s')
    
    # put 0 at bottom and top
    N2 = N2.where(N2.s_w > -1, other=0)#.values = 0*N2.isel(s_w=1).values
    N2 = N2.where(N2.s_w < 0, other=0)
    # copy nearest values at bottom and top
    N2bot = N2.shift(s_w=-1).where(N2.s_w == -1, other=0)
    N2top = N2.shift(s_w=1).where(N2.s_w == 0, other=0)
    N2 += N2bot + N2top
    return N2.rename("bvf")


###########################################################################
########################## Interpolation routines #########################
###########################################################################

### this should be in gridop 
def zi_w2rho(grid, data, z_w=None, z_r=None):
    """ interpolate linearly from z_w grid to z_r grid
    warning: this version uses grid (xgcm)
    N.B.: z_r, z_w can the grid at any time or at rest (zeta=0) 
    if z_w is None, will try to get the values from data
    if z_r is None, will interpolate at midpoints using grid.interp
    if z_w is None but z_r is not, will fail """
    
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
    ''' interpolate vertically WARNING z must be first dim in z0
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
    
def interp_coredim(z0, z, v, b_extrap, t_extrap):
    ''' interpolate vertically
    when using core_dims=..., the core dims are moved to last axis
    here I transpose the fields to overcome this
    .not robust. can't plot field for instance, because 
    error "replacement data must match the Variable's shape"
    '''
    # check v and z have identical shape
    assert v.ndim==z.ndim
    # test if temporal dimension is present
    if v.ndim == 4: # ! WARNING not tested !!!
        vi = [interp2z_3d(z0.T, z[...,t,:].T, v[...,t,:].T, b_extrap, t_extrap).T[...,None] 
                  for t in range(v.shape[-1])]
        return np.concatenate(vi, axis=0) # (50, 722, 258, 1)
        #return v*0 + v.shape[3]
    else:
        return interp2z_3d(z0.T, z.T, v.T, b_extrap, t_extrap).T # transpose

def z2zmoy(data, zz, zmoy):
    # if not re-chunking, nanny restarts worker
    if "s_rho" in data.dims: data = data.chunk({"s_rho":-1})
    elif "s_w" in data.dims: data = data.chunk({"s_w":-1})
    prov = xr.apply_ufunc(interp2z, zmoy, zz, data, 2, 2,
                    dask='parallelized', output_dtypes=[np.float64])
    return prov.assign_coords(z_rho=zmoy)

def z2zwoy(data, zz, zmoy):
    # if not re-chunking, nanny restarts worker
    # to avoid broadcasting when interpolating from rho to w, need to specify core dims
    if "s_rho" in data.dims: 
        prov = xr.apply_ufunc(interp_coredim, zmoy.chunk({"s_w":-1}), zz.chunk({"s_rho":-1}),                           data.chunk({"s_rho":-1}), 2, 2, dask='parallelized', 
                        output_dtypes=[np.float64], 
                    input_core_dims=[["s_w"],["s_rho"],["s_rho"],[],[]], 
                    output_core_dims = [["s_w"]])
    else:
        prov = xr.apply_ufunc(interp2z, zmoy.chunk({"s_w":-1}), zz.chunk({"s_w":-1}),                           data.chunk({"s_w":-1}), 2, 2,                     dask='parallelized', output_dtypes=[np.float64])
    return prov.assign_coords(z_w=zmoy)


################################################################################
#####################  Vertical integration  ###################################
################################################################################
def trapzw(grid, data, zz=None):
    """ trapz for a function on w points (zz is also on w points)
    (b-a)*.5*(f(a)+f(b)) = sum(dz*grid.interp(f, "s"))"""
    if zz is None:
        zz = data.zz
    return (grid.diff(zz, "s")*grid.interp(data, "s")).sum("s_rho")

# def trapz(grid, data, zz):
#    """ wrapper of np.trapz. Does not work, everything crashes """
#    if "s_w" in data.dims:
#        indim = "s_w"
#    elif "s_rho" in data.dims:
#        indim = "s_rho"
#    else:
#        raise ValueError("core dim not found")
#    #print("found core dim:", indim, data.dims, zz.dims)
#    res = xr.apply_ufunc(np.trapz, data, zz,dask='parallelized', input_core_dims = [[indim],[indim]], 
#                output_dtypes=[np.float64], exclude_dims=set((indim,)))
    

################################################################################
#####################  Projection on vertical modes  ###########################
################################################################################
def proj_p(data, dm, zz=None):
    """ projection on p-mode """
    grid = dm.attrs["xgcm-Grid"]
    if zz is not None:
        data = z2zmoy(data, zz, dm.z_rho)    
    phin = dm.phi
    dz = grid.diff(dm.z_w, "s")
    return (dz*data*phin).sum("s_rho")/dm.norm

def proj_w(data, dm, zz=None): 
    """ using what I call "w-modes".
    for reconstructing, use w = wn*modw 
    interpolation uses linear interpolation, but midpoints should be OK 
        (since it gives something equivalent to trapezoidal integration upon integration)
    """
    grid = dm.attrs["xgcm-Grid"]
    dz = grid.diff(dm.z_w, "s")
    # this should be similar to -N2/c**2*dphidz
    phiw = grid.cumsum(dz*dm.phi, "s", to="outer", boundary="fill")\
            .chunk({"s_w":-1})
    if zz is not None:
        data = z2zmoy(data, zz, dm.z_rho)        
    prov = (data*zi_w2rho(grid, phiw*dm.N2, dm.z_w, dm.z_rho)*dz).sum(dim="s_rho") \
            + _g* (phiw*grid.interp(data, "s", boundary="extrapolate")).isel(s_w=-1)
    return prov/dm.norm/dm.c**2

def proj_b(data, dm, zz=None): #zz, phib=dm.dphidz, dz=dz):
    """ N.B.: for reconstructing b, use -c**2*bn*dphidz """
    grid = dm.attrs["xgcm-Grid"]
    phib, N2 = dm.dphidz, dm.N2
    dz = grid.diff(dm.z_w, "s")
    if zz is not None:
        data = z2zmoy(data, zz, dm.z_rho)    
    prov = (data*grid.interp(phib/N2, "s")*dz).sum("s_rho") \
            + _g* (grid.interp(data, "s", boundary="extrapolate")*phib/N2**2).isel(s_w=-1)
    return -prov/dm.norm #data.isel(s_rho=-1)

if False: # routines not updated
    def projww(data, zz, phiw=modw, dz=dz, cn=dm.c): # not very cool
        """ using what I call "w-modes" and w data on w-points.
        for reconstructing, use w = wn*modw """
        prov = z2zwoy(data, zz)
        prov = (grid.interp(prov*phiw*N2, "s")*dz).sum(dim="s_rho")             + grav* (phiw*prov).isel(s_w=-1)
        return prov/ds.h/cn**2

    def projwb(data, zz, phiw=modw, dz=dz, cn=dm.c): # not very cool
        """ using what I call "w-modes".
        for reconstructing, use w = wn*modw """
        prov = z2zwoy(data, zz)
        prov = (prov*phiw*N2).integrate(dim="s_w")             + grav* (phiw*prov).isel(s_w=-1)
        return prov/ds.h/cn**2

    def projw_wb(data, zz, phib=dm.dphidz, dz=dz):
        """ this one is using varphi modes. You can verify that it does not change a lot 
        for reconstructing w, use -c^2/N^2*wn*dphidz """
        prov = z2zmoy(data, zz)
        prov = (prov*grid.interp(phib, "s")*dz).sum(dim="s_rho")             + grav* (phib/N2*grid.interp(prov, "s")).isel(s_w=-1)
        return -prov/ds.h

################################################################################
#####################  Horizontal derivatives  #################################
################################################################################

def diff_xy(grid, data, axis, pm, ongrid=True, docor=False, zr=None, zw=None, doverb=False):
    """ diff_xy
    general routine (wrapper) for horizontal derivative. xarray implementation, using xgcm-Grid.
    First order staggered grid or second order on same grid
    can have correction for sigma-to-z if z_rho, z_w are present in data-array
    WARNINGS: 
       - not tested on non-cartesian grid, although it should be OK
       - passages from one grid to another (and diff) use "extrapolate" argument in grid
     """
    if ongrid:
        return diff_xy_ongrid(grid, data, axis, pm=pm, docor=docor, zr=zr, zw=zw, doverb=doverb)
    else:
        return diff_xy_offgrid(grid, data, axis, pm=pm, docor=docor, zr=zr, zw=zw, doverb=doverb)
    
    
def diff_xy_ongrid(grid, data, axis, pm=None, docor=False, zr=None, zw=None, doverb=False):
    if axis in ["x","xi"]:
        axis, coord = "x", "xi"
    elif axis in ["y","eta"]:
        axis, coord = "y", "eta"
    else:
        raise ValueError("could not recognize axis", axis)
    bnd = "extrapolate" # set how grid treats boundary. Overrided by periodic of activated 
    
    res = grid.interp(diff_xy_offgrid(grid, data, axis, pm=pm, docor=False, doverb=doverb), coord, boundary=bnd)

    if docor: 
        res += diffxy_corsz(grid, data, axis, ongrid=True, pm=pm, zr=zr, zw=zw, doverb=doverb)
    
    return res

def diff_xy_offgrid(grid, data, axis, pm=None, docor=False, zr=None, zw=None, doverb=False):
    if axis in ["x","xi"]:
        axis, coord = "x", "xi"
    elif axis in ["y","eta"]:
        axis, coord = "y", "eta"
    else:
        raise ValueError("could not recognize axis", axis)
    bnd = "extrapolate" # set how grid treats boundary. Overrided by periodic of activated 
    
    res = grid.diff(data, coord, boundary=bnd)
    
    if pm is not None:
        if doverb: print("using pm")
        if not isinstance(pm, float):
            dim = [dim for dim in pm.dims if dim.startswith(axis+"_")][0]
            if doverb: print("recognized pm dim", dim)
            if dim not in res.dims and not isinstance(pm, float):
                pm = grid.interp(pm, coord, boundary=bnd)
        res *= pm
    else:
        dim = [dim for dim in data.dims if dim.startswith(axis+"_")][0]
        if doverb: print("recognized data dim", dim)
        if dim in data.coords:
            if doverb: print("using coord for computing dl")
            dl = grid.diff(data.coords[dim], coord, boundary=bnd)
            res /= dl
        else:
            print('no dl or pm found, returning a non-dimensional diff')
    
    if docor:
        if zr is None:
            zr = data.z_rho
        if zw is None:
            zw = data.z_w
        res += diffxy_corsz(grid, data, axis, ongrid=False, pm=pm, zr=zr, zw=zw, doverb=doverb)
    return res

def diffxy_corsz(grid, data, axis, pm=None, zr=None, zw=None, ongrid=False, doverb=False):
    """ Compute correction associated with horizontal derivative of z levels
    for ongrid or offgrid 
    Configurations available: 
       - if ongrid, data and z are on the same grid
       - if offgrid, data and z are on staggered grids, with z on target grid
    """
    if axis in ["x","xi"]:
        axis, coord = "x", "xi"
    elif axis in ["y","eta"]:
        axis, coord = "y", "eta"
    else:
        raise ValueError("could not recognize axis", axis)
    bnd = "extrapolate" # set how grid treats boundary. Overrided by periodic of activated 
    
    dsig = grid.interp(grid.diff(data, "s", boundary=bnd), "s", boundary=bnd) # on data pnts
    hz = grid.diff(zw, "s") # on z points = data-scat
    dlz = grid.interp(grid.diff(zr, coord, boundary=bnd), coord, boundary=bnd) # on z points
    
    ### results on scattered points
    if ongrid:
        res = -dsig / hz*dlz
    else:
        res = -grid.interp(dsig, coord, boundary=bnd)/hz*dlz # 
    
    if pm is not None:
        if doverb: print("using pm")
        if not isinstance(pm, float):
            dim = [dim for dim in pm.dims if dim.startswith(axis+"_")][0]
            if doverb: print("recognized pm dim", dim)
            if dim not in res.dims and not isinstance(pm, float):
                pm = grid.interp(pm, coord, boundary=bnd)
        res *= pm
    else:
        dim = [dim for dim in data.dims if dim.startswith(axis+"_")][0]
        if doverb: print("recognized data dim", dim)
        if dim in data.coords:
            if doverb: print("using coord for computing dl")
            dl = grid.diff(data.coords[dim], coord, boundary=bnd)
            res /= dl
        elif doverb:
            print('no dl or pm found, returning a non-dimensional diff')
    
    return res

# wrappers
def diff_x(grid, data, dx=1, ongrid=True, docor=False, zr=None, zw=None): 
    return diff_xy(grid, data, "x", dx, ongrid=ongrid, order=order, docor=docor, zr=zr, zw=zw)

def diff_y(grid, data, dy=1, ongrid=True, docor=False, zr=None, zw=None): 
    return diff_xy(grid, data, "y", dy, ongrid=ongrid, order=order, docor=docor, zr=zr, zw=zw)



################################################################################
############################  Plots  ###########################################
################################################################################
def plt_imodx(data, ax=None, x=None, y=None, **kwargs):
    """ plotting x (or y) vs modes 2D quantity, using pcolormesh through xarray
    data must be a xarray DataArray """
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    args = {}
    if ax is not None:
        args['ax'] = ax
    if x is not None:
        args['x'] = x
    if y is not None:
        args['y'] = y
    if "mode" in data.dims:
        nmod = data['mode'].size-1
    elif "modm" in data.dims:
        nmod = data['modm'].size-1
    else:
        nmod = Nmodes 
        print("WARNING Nmodes is hard coded")
    for key, val in kwargs.items():
        args[key] = val
        
    hpc = data.plot(add_colorbar=False, **args)
    if ax is None:
        ax = plt.gca()
    ax.set_yticks(np.arange(nmod+1))
    ax.set_yticks(np.arange(nmod)+.5, minor=True)
    ax.set_ylim([-.5, nmod+.5])
    ax.grid(True, axis="x")
    ax.grid(False, which="major", axis="y")
    ax.grid(True, which="minor", axis="y")
    cbax = inset_axes(ax, width="3%", height="80%", loc="lower right") 
    cbl = plt.colorbar(hpc, cax=cbax)
    return hpc
