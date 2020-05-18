""" vmodes.py python module
Defines class Vmodes and routines for computing vertical modes given a stratification profile
the corresponding Sturm-Liouville problem is (phi'/N2)' + lam^2 phi = 0, with proper BCs. 
"""
import os
from glob import glob

import scipy.sparse as sp
import scipy.sparse.linalg as la
import numpy as np
import xarray as xr

from matplotlib import pyplot as plt

from . import gridop as gop
from .postp import grav, _move_singletons_as_attrs
from . import plot as cplt

# default values
_sig = .1
_nmodes = 10
_free_surf = True

# core attribute and variables
_core_variables = ['phi', 'dphidz', 
                   'c', 'norm']
_core_attrs = ['nmodes', '_xgrid_z', 'g', 
               'free_surf', 'sigma']

class Vmodes(object):
    """ Class for vertical modes computation and manipulation (projection and reconstruction)
    The vertical modes 'phi' are defined following the standard Sturm-Liouville problem 
    (e.g. Gill, 1982 or any standard GFD textbook), see below.
    Those are referred as 'pressure-like' modes. 
    Their vertical derivatives, named 'dphidz', are referred as 'w-like' modes.

    Parameters:
    ___________
    xgrid: xgcm.Grid 
        xgcm grid object associated with the dataset, required for grid manipulations
    zc: xarray.Datarray
        vertical grid at cell centers
    zf: xarray.Datarray
        vertical grid at cell interfaces
    N2: xarray.Datarray
        Brunt-vaisala frequency, at cell interfaces
    nmodes: int, optional
        number of baroclinic modes to compute (barotropic mode will be added)
    free_surf: bool, optional
        whether using a free-surface condition or rigid lid
    persist: bool, optional
        persist the dask dataset containing the modes or not
    grav: float, optional
        gravity constant
    sigma: float, optional
        parameter for shift-invert method in scipy.linalg.eig (default: _sig)
        
    Attributes:
    ___________
    xgrid: xgcm.Grid
    nmodes: float
    ds: xarray.Dataset
        dataset containing the vertical modes and associated data
    g: float
    free_surf: float
    sigma: float

    Methods:
    ________
    project(data, vartype="p", **kwargs)
        project 'data' on vertical modes
    project_puv(data, z=False, sel=None, align=True)
        project 'data' on pressure-like vertical modes
    project_w(data, z=False, sel=None, align=True) 
        project 'data' on w-like vertical modes
    project_b(data, z=False, sel=None, align=True)
        project 'data' on b-like vertical modes
    reconstruct(projections, vartype=None, **kwargs)
        reconstruct field from projections coefficient (modal amplitudes)
    reconstruct_puv(projections, sel=None, align=True)
        reconstruct field from projections coefficient using pressure-like modes
    reconstruct_w(projections, sel=None, align=True)
        reconstruct field from projections coefficient using w-like modes
    reconstruct_b(projections, sel=None, align=True)
        reconstruct field from projections coefficient using w-like modes

    Notes:
    ______
    zc, zf, N2 should be ordered from bottom to top, and zf and N2 have the same z dimension
    The vertical modes are definied following the equation:
    .. math:: (\phi'/N^2)' + k^2\phi=0 
    with boundary conditions
        :math:`\phi'=0` at the bottom
        :math:`g\phi' + N^2\phi=0` at the surface (or :math:`\phi'=0` for a rigid lid condition). 
    Computation of the vertical modes is performed using second order finite difference
    The eigenvalue used in this class is c=1/k.
    """
    def __init__(self, xgrid, zc, zf, N2, 
                 nmodes=_nmodes,
                 free_surf=_free_surf,
                 persist=False,
                 g=grav, sigma=_sig):
        """ Create a Vmode object
        
        Parameters:
        ___________
        xgrid: xgcm.Grid 
            xgcm grid object associated with the dataset, required for grid manipulations
        zc: xarray.Datarray
            vertical grid at cell centers
        zf: xarray.Datarray
            vertical grid at cell interfaces
        N2: xarray.Datarray
            Brunt-vaisala frequency, at cell interfaces
        nmodes: int, optional
            number of baroclinic modes to compute (barotropic mode will be added)
        free_surf: bool, optional
            whether using a free-surface condition or rigid lid
        persist: bool, optional
            persist the dask dataset containing the modes or not
        grav: float, optional
            gravity constant
        sigma: float, optional
            parameter for shift-invert method in scipy.linalg.eig (default: _sig)
            
        Notes:
        ______
        zc, zf, N2 should be ordered from bottom to top, and zf and N2 have 
        the same z dimension
        """
        self.xgrid = xgrid
        self.nmodes = nmodes
        self._znames = {"zc": zc.name, "zf": zf.name}
        sdim = {"zc": gop.get_z_dim(zc)[0], "zf": gop.get_z_dim(zf)[0]}
        self._zdims = sdim
        xgrid_z = gop.get_xgrid_ax_name(xgrid,sdim.values())
        self._xgrid_z = xgrid_z
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
        
        self.ds['H'] = np.abs(zf.isel({sdim['zf']:0}))
        self.ds['dz'] = xgrid.diff(zf, xgrid_z).rename("dz")
        
        # in order to use a consistent color code
        self._colors = cplt.get_cmap_colors(self.nmodes)
        
        if persist:
            self.ds = self.ds.persist()
     
    def __getitem__(self, item):
        """ Enables calls such as vm['N2']
        """
        if item in ['zc', 'zf']:
            return self.ds[self._znames[item]]
        return self.ds[item]

    def __repr__(self):
        strout = 'Vmode object with dimensions {}\n'.format(tuple(self.ds.dims.keys()))
        strout += '  Number of modes = {}\n'.format(self.nmodes)
        strout += '  N2, min/max = {0:.1e}, {1:.1e}\n'.format(\
                                self['N2'].min().values, self['N2'].max().values)
        strout += '  Options / parameters: g={0:.2f}, free_surf={1}, eig_sigma={2:.1e}\n'\
                            .format(self.g, self.free_surf, self.sigma)
        return strout
    
    def plot(self, mode_max=3):
        """ Show modal structure along with celerities
        
        Parameters
        ----------
        mode_max: int, optional
            Maximum mode number to show
        """    
        ds = self.ds.sel(mode=slice(mode_max))
        
        if mode_max<self.nmodes:
            colors = cplt.get_cmap_colors(mode_max)
        else:
            colors = self._colors

        fig, axes = plt.subplots(1,2)
        #
        ax = axes[0]
        (ds['N2']*1e5).plot.line(ax=ax, y=self._znames['zf'], 
                                 color='0.5', lw=2)
        for n, col in zip(ds['mode'], colors):
            c = float(ds.c.sel(mode=n))
            ds['phi'].sel(mode=n).plot.line(ax=ax, 
                                            y=self._znames['zc'],
                                            color=col,
                                            label='{:.1f} m/s'.format(c)
                                           )
        ax.set_title('phi')
        ax.grid()
        ax.legend()
        #
        ax = axes[1]
        (ds['N2']*1e3).plot.line(ax=ax, y=self._znames['zf'], 
                                 color='0.5', lw=2)
        for n, c in zip(ds['mode'], colors):
            ds['dphidz'].sel(mode=n).plot.line(ax=ax, 
                                               y=self._znames['zf'], 
                                               color=c
                                              )
        ax.set_title('dphidz')
        ax.grid()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        #return fig, axes

    def _compute_vmodes(self):
        """ compute vertical modes and store the results into the dataset 
        wrapper of external function get_vmodes """
        dm = get_vmodes(self['zc'],
                        self['zf'],
                        self.ds.N2,
                        nmodes=self.nmodes,
                        free_surf=self.free_surf,
                        sigma=self.sigma, g=self.g, 
                        z_dims=[self._zdims['zc'], self._zdims['zf']])
        self.ds = xr.merge([self.ds, dm], compat="override")
    
    # projections
    def project(self, data, vartype="p", **kwargs):
        """ Project a variable on vertical modes (p-modes or w-modes)
        Internally calls project_puv or project_w
        
        Parameters:
        ___________
        data: xarray.DataArray
            array containing the data to be projected. 
        vartype: str, optional (default: "p", i.e. pressure modes)
            string specifying whether projection should be done w-modes ("w"), 
            buoyancy modes ("b") or pressure modes (any other value)
        z: xarray.DataArray or str or bool, optional (default: False)
            array containing the z-grid of the data, or string containing the name of the z coord, 
            or boolean saying whether we should interpolate (finding the z-coord by its own). 
            The data will be interpolated onto the vmodes grid points prior to projection.
        sel: Dict, optional (default: None)
            indices applied to the vmodes dataset prior to projection
        align: bool, optional (default: True)
            whether alignment between the data and vmodes DataArray should be performed before projecting.
            Note that any mismatch between data and vertical modes horizontal grids will raise an error
            (horizontal interpolation needs to be performed prior to projection)
        
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

    def project_puv(self, data, z=False, sel=None, align=True):
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

        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
            
        if align:
            data, dm = xr.align(data, dm, join="inner")
            
        _check_hdim_mismatch(data, dm)
            
        if not( z is None or z is False ):
            if z is True:
                z, = gop.get_z_coord(data)
            if isinstance(z, str):
                z = data.coords[z]
            elif align:
                data, z = xr.align(data, z, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], z, data)
            
        res = (dm.dz*data*dm.phi).sum(self._zdims['zc'])/dm.norm

        return res
    
    def project_w(self, data, z=False, sel=None, align=True): 
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
        
        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
            
        if align:
            data, dm = xr.align(data, dm, join="inner")
            
        _check_hdim_mismatch(data, dm)

        if not( z is None or z is False ):
            if z is True:
                z, = gop.get_z_coord(data)
            if isinstance(z, str):
                z = data.coords[z]
            elif align:
                data, z = xr.align(data, z, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], z, data)
        
        zf, zc = self._znames['zf'], self._znames["zc"]
        prov = (data * self._w2rho(-dm.dphidz, zc=dm[zc], zf=dm[zf]) * dm.dz).sum(self._zdims['zc'])
        
        if self.free_surf:
            prov += self.g * ( -dm.dphidz / dm.N2 
                              * self.xgrid.interp(data, self._xgrid_z, boundary="extrapolate") 
                             ).isel({self._zdims["zf"]:-1}).drop(self._znames["zf"])
       
        return prov/dm.norm
    
    def project_b(self, data, z=False, sel=None, align=True): 
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
        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
            
        if align:
            data, dm = xr.align(data, dm, join="inner")
            
        _check_hdim_mismatch(data, dm)
        
        if not( z is None or z is False ):
            if z is True:
                z, = gop.get_z_coord(data)
            if isinstance(z, str):
                z = data.coords[z]
            elif align:
                data, z = xr.align(data, z, join="inner")
            data = gop.interp2z(dm[self._znames['zc']], z, data)

        zf, zc = self._znames['zf'], self._znames["zc"]
        prov = (data * self._w2rho(dm.dphidz/dm.N2, zc=dm[zc], zf=dm[zf])* dm.dz).sum("s_rho")
        
        if self.free_surf:
            prov += self.g * (self.xgrid.interp(data, self._xgrid_z, boundary="extrapolate")
                          * dm.dphidz / dm.N2**2
                         ).isel({self._zdims['zf']:-1}).drop(self._znames["zf"])
        
        return -prov/dm.norm 

    # reconstructions
    def reconstruct(self, projections, vartype=None, **kwargs):
        """ Reconstruct a variable from modal amplitudes
        Internally call reconstruct_puv or reconstruct w or reconstruct_b
        
        Parameters:
        ___________
        projections: xarray.DataArray
            array containing the modal projection coefficients (modal amplitudes)
        vartype: {"p", "u", "v", "w", "b"}, optional 
            string specifying whether reconstruction should be done using w-modes ("w"), buoyancy modes ("b") or pressure modes ("p", "u" or "v"). Default is "p".
        sel: Dict, optional (default: None)
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
            if sum(s in projections.name.lower() for s in vartyps)==1:
                vartype = next((s for s in vartyps if s in projections.name.lower()))
            else: 
                raise ValueError("unable to find what kind of basis to use for reconstruction")
                
        if vartype in "puv":
            return self.reconstruct_puv(projections, **kwargs)
        elif vartype == "w":
            return self.reconstruct_w(projections, **kwargs)
        elif vartype == "b":
            return self.reconstruct_b(projections, **kwargs)
        
    def reconstruct_puv(self, projections, sel=None, align=True):
        """ Reconstruct a variable from modal amplitudes, using pressure-like modes
            
        See also
        ________
        reconstruct
        
        """

        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
        if align:
            projections, dm = xr.align(projections, dm, join="inner")    
        return (projections * dm.phi).sum("mode")

    def reconstruct_w(self, projections, sel=None, align=True):
        """ Reconstruct a variable from modal amplitudes, using w-like modes and w-like normalization
            
        See also
        ________
        reconstruct
        
        """
        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
        if align:
            projections, dm = xr.align(projections, dm, join="inner")    
        return (-dm.c**2 / dm.N2 * dm.dphidz * projections).sum("mode")
    
    def reconstruct_b(self, projections, sel=None, align=True):
        """ Reconstruct a variable from modal amplitudes, using w-like modes and b-like normalization
            
        See also
        ________
        reconstruct
        
        """
        if sel is None:
            dm = self.ds
        else:
            dm = self.ds.sel(sel)
        if align:
            projections, dm = xr.align(projections, dm, join="inner")    
        return (-dm.c**2 * dm.dphidz * projections).sum("mode")
    
    # utilitaries 
    def _w2rho(self, data, zf=None, zc=None, align=True):
        """ routine to interpolate fields from mean rho point to mean w point
        adapt to different ensemble of points if specified
        wrapper to gop.w2rho
        could be replaced by xgrid.interp once correct metrics are provided
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
        return gop.w2rho(data, self.xgrid, zc, zf, s_dims=[self._zdims["zc"],self._zdims["zf"]])

    # store/load 
    def store(self, file_path, projections=None, **kwargs):
        """ Store Vmodes object along with projections into zarr archives
        
        Parameters
        ----------
        file_path: str
            File path where data will be stored, e.g. '/home/roger/vmodes.zarr'
        projections: xarray.Dataset
            Dataset containing projections
        **kwargs:
            kwargs passed to `to_zarr`
        """
        ds = self._wrap_in_dataset()
        _file = file_path
        # undesired singleton time coordinates
        ds = _move_singletons_as_attrs(ds)
        #
        ds.to_zarr(file_path, **kwargs)
        print('Store vertical modes in {}'.format(file_path))
        #
        if projections:
            _file = file_path.strip('.zarr')+'_projections.zarr'
            ds = _move_singletons_as_attrs(projections)
            ds.to_zarr(_file, **kwargs)
            print('Store projections in {}'.format(_file))
    
    def _wrap_in_dataset(self):
        """ wrap critical data in a xarray dataset
        """
        # we won't store xgrid as the intend is to store this in a zarr archive
        #
        # One should try to figure out whether minimal information can be extracted
        # from the xgrid object and store as attribute to the xarray dataset
        ds = self.ds.assign_attrs(**{a: getattr(self, a) for a in _core_attrs})
        ds = ds.assign_attrs(zc_name=self._znames['zc'],
                             zf_name=self._znames['zf'])    
        return ds


def load_vmodes(file_path, xgrid, persist=False):
    """ load vertical modes from a file
    
    Parameters
    ----------
    file_path: str
        Path to vmode datafile (zarr archive)
    xgrid: xgcm.Grid object
        Required for grid manipulations
    """
    ds = xr.open_zarr(file_path)
    vm = Vmodes(xgrid,
                ds[ds.zc_name], ds[ds.zf_name], 
                ds['N2'], 
                nmodes=ds.nmodes,
                free_surf=ds.free_surf,
                persist=persist,
                g=ds.g, sigma=ds.sigma)
    # transfer mode from ds to vm
    for v in _core_variables:
        if v in ds:
            vm.ds[v] = ds[v]
        elif v in ds.attrs:
            # singleton cases (norm)
            vm.ds[v] = ds.attrs[v]
    if persist:
        vm.ds = vm.ds.persist()
    # search for projections:
    _pfile = file_path.strip('.zarr')+'_projections.zarr'
    if os.path.isdir(_pfile):
        projections = xr.open_zarr(_pfile)
        if persist:
            projections = projections.persist()
        return vm, projections
    else:
        return vm

def _check_hdim_mismatch(data, dm):
    # test if mismatching horizontal dimensions will trigger broadcasting
    dims_data = gop._get_spatial_dims(data)
    dims_vm = gop._get_spatial_dims(dm)
    for d in ['x', 'y']:
        assert (dims_vm[d]==None) or (dims_data[d]==None) or (dims_vm[d]==dims_data[d]), \
            'Data and vertical mode have mismatching horizontal dimensions ' \
            +' ({},{}), '.format(dims_vm[d], dims_data[d]) \
            +'consider horizontally interpolating data first'
    
def get_vmodes(zc, zf, N2, nmodes=_nmodes, **kwargs):
    """ compute vertical modes
    Wrapper for calling `compute_vmodes` with DataArrays through apply_ufunc. 
    z levels must be in ascending order (first element is at bottom, last element is at surface)
    
    Parameters:
    ___________
    zc: DataArray
        z-levels at center of cells
    zf: DataArray
        z-levels at cell faces (vertical dim is one element longer than zc)
    N2: DataArray
        Brunt-Vaisala Frequency at cell faces
    nmodes: int, optional
        number of vertical baroclinic modes (barotropic is added)
    free_surf: bool, optional
        whether to use free surface boundary condition or not
    sigma: scalar or None, optional
        parameter for shift-invert method in scipy.linalg.eig (default: _sig)
    g: scalar, optional
        gravity constant
    z_dims: list of str, optional
        vertical dimension names in zc, zf (default: "s_rho", "s_w")

    Returns:
    ________
    xarray.DataSet: vertical modes (p and w) and eigenvalues
    
    See Also:
    _________
    compute_vmodes: routine for computing vertical modes from numpy arrays
    Vmodes: class for computing and manipulating vertical modes
   
    """
    kworg = {"nmodes":nmodes, "stacked":True}
    z_dims = None
    if kwargs is not None:
        z_dims = kwargs.pop("z_dims", None)
        kworg.update(kwargs)
    if z_dims:
        s_rho, s_w = z_dims
    else:
        s_rho, s_w = "s_rho", "s_w"

    N = zc[s_rho].size
    res = xr.apply_ufunc(compute_vmodes, 
                         zc.chunk({s_rho:-1}), 
                         zf.chunk({s_w:-1}),
                         N2.chunk({s_w:-1}), 
                         kwargs=kworg, 
                         input_core_dims=[[s_rho],[s_w],[s_w]],
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
           .rename({'s_stack': s_rho})
           #.assign_coords(z_rho=zc)
          )
    dphidz = (res.isel(s_stack=slice(N+1,2*N+2))
              .rename('dphidz')
              .rename({'s_stack': s_w})
            #  .assign_coords(z_w=zf)
             )
    # merge data into a single dataset
    other_dims = tuple([dim for dim in zc.dims if dim!=s_rho]) # extra dims
    norm = -zf.isel({s_w:0}, drop=True).rename('norm')
    dm = (xr.merge([c, phi, dphidz, norm])
          .transpose(*('mode',s_rho,s_w)+other_dims)
         )
    return dm  ### hard-coded norm = H

def compute_vmodes(zc_nd, zf_nd, N2f_nd, 
                   nmodes=_nmodes, free_surf=_free_surf,
                   g=grav,
                   sigma=_sig, stacked=True,
                   **kwargs):
    """
    Compute vertical modes from stratification.
    wrapper for vectorizing compute_vmodes_1D over dimensions of numpy arrays other than the vertical one

    Parameters:
    ___________
    zc_nd: (...,N) ndarray
        vertical grid at cell centers
    zf_nd: (...,N+1) ndarray
        vertical grid at cell interfaces
    N2f_nd: (...,N+1) ndarray
        Brunt-Vaisala frequency at cell interfaces
    nmodes: int, optional
        number of baroclinic modes to compute (barotropic mode will be added)
    free_surf: bool, optional
        whether using a free-surface condition or rigid lid
    g: float, optional
        gravity constant
    sigma: float, optional
        parameter for shift-invert method in scipy.linalg.eig (default: _sig)
    stacked: bool
        whether the result should be returned as a stacked array of size (..., 2N+2). Default is True

    Returns:
    res: ndarray (...,2N+2,nmodes) 
        stacked array of eigenvalues, p-like modes and w-like modes
    OR
    c: (...,nmodes) ndarray
        eigenvalues (pseudo phase speed, c=1/sqrt(k))
    phi: (...,N,nmodes) ndarray
        p-like modes at cell centers
    dphidz: (...,N+1,nmodes) ndarray
        w-like modes at cell interfaces

    Notes:
    ______
    z dimension must be the last axis 
    you can use this if you are using numpy (but make sure z is last dim)
    The vertical modes are definied following the equation:
    .. math:: (\phi'/N^2)' + k^2\phi=0 
    with boundary conditions:
        :math:`\phi'=0` at the bottom
        :math:`g\phi' + N^2\phi=0` at the surface (or :math:`\phi'=0` for a rigid lid condition). 
    Computation of the vertical modes is performed using second order finite difference
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
        cn, phin, dphi = compute_vmodes_1D(zc_nd, zf_nd, N2f_nd, nmodes, \
                                free_surf=free_surf, g=g, sigma=sigma)
        if stacked:
            return np.vstack([cn[None,:],phin,dphi]) #.reshape(nxy+(-1,nmodes+1))   
        else:
            return cn, phin, dphi
            
def compute_vmodes_1D(zc, zf, N2f, 
                      nmodes=_nmodes, free_surf=True, 
                      g=grav, sigma=_sig):
    """
    Compute vertical modes from stratification.

    Parameters:
    ___________
    zc: (N) ndarray
        vertical grid at cell centers
    zf: (N+1) ndarray
        vertical grid at cell interfaces
    N2f: (N+1) ndarray
        Brunt-Vaisala frequency at cell interfaces
    nmodes: int, optional
        number of baroclinic modes to compute (barotropic mode will be added)
    free_surf: bool, optional
        whether using a free-surface condition or rigid lid
    g: float, optional
        gravity constant
    sigma: float, optional
        parameter for shift-invert method in scipy.linalg.eig (default: _sig)

    Returns:
    ________
    c: (nmodes) ndarray
        eigenvalues (pseudo phase speed, c=1/sqrt(k))
    phi: (N,nmodes) ndarray
        p-like modes at cell centers
    dphidz: (N+1,nmodes) ndarray
        w-like modes at cell interfaces

    Notes:
    ______
    copy-pasted from M. Dunphy's vmodes_MD.py
    WARNING maybe the barotropic mode for rigid lid is wrong 
    The vertical modes are definied following the equation:
    .. math:: (\phi'/N^2)' + k^2\phi=0 
    with boundary condition :math:`\phi'=0` at the bottom and :math:`g\phi' + N^2\phi=0` at the surface (or :math:`\phi'=0` for a rigid lid condition). 
    Computation of the vertical modes is performed using second order finite difference

    """ 
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

# --------------------------------- model modes -----------------------------------

def open_vmode_file(f, r):
    ds = xr.open_dataset(f)
    for d in ds.dims:
        if ds[d].size==r.N:
            ds = ds.rename({d: 's_rho'})
            ds['s_rho'] = r['grid']['s_rho']
        elif ds[d].size==r.N+1:
            ds = ds.rename({d: 's_w'})
            ds['s_w'] = r['grid']['s_w']
        else:
            ds = ds.rename({d: 'mode'})
            ds['mode'] = np.arange(0,ds['mode'].size)
    ds = ds.rename(array=f.split('/')[-1].strip('vmodes_').replace('.nc',''))
    return ds

def load_vmodes_from_files(r, tnum=1):
    """ Load vertical modes from online diagnostic files and recompute modes
    
    Parameters
    ----------
    r: postp.run
    tnum: int, optional
        t directory that is loaded, default is 1.
    """
    # load online outputs
    _files = glob(os.path.join(r.dirname,'t{}/vmodes*.nc'.format(tnum)))
    vmodes_croco = (xr.merge([open_vmode_file(f,r)
                             for f in _files])
                    .rename({'zr': 'z_rho', 'zw': 'z_w'})
                   )
    # recompute vertical modes from output stratification
    vmodes = Vmodes(r['xgrid'],
                    vmodes_croco.z_rho, vmodes_croco.z_w, 
                    vmodes_croco.N2, 
                    vmodes_croco.mode.size, 
                    persist=True)
    return vmodes_croco, vmodes
