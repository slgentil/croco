import os

import numpy as np
import xarray as xr

from dask import delayed
import threading

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import scipy.io

# -------------------------------- various utils -------------------------------

def get_cmap_colors(Nc, cmap='plasma'):
    """ load colors from a colormap to plot lines
    
    Parameters
    ----------
    Nc: int
        Number of colors to select
    cmap: str, optional
        Colormap to pick color from (default: 'plasma')
    """
    scalarMap = cmx.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=Nc),
                                   cmap=cmap)
    return [scalarMap.to_rgba(i) for i in range(Nc)]

def DefCmap():
    """ construct a python colormap from matlab mat file, stored in the same directory as this module"""
    matfile = scipy.io.loadmat(os.path.join(os.path.dirname(__file__),'map_64_wc.mat'))
    return array2cmap(np.array(matfile['cm']))


def array2cmap(X):
    N = X.shape[0]

    r = np.linspace(0., 1., N + 1)
    r = np.sort(np.concatenate((r, r)))[1:-1]

    rd = np.concatenate([[X[i, 0], X[i, 0]] for i in range(N)])
    gr = np.concatenate([[X[i, 1], X[i, 1]] for i in range(N)])
    bl = np.concatenate([[X[i, 2], X[i, 2]] for i in range(N)])

    rd = tuple([(r[i], rd[i], rd[i]) for i in range(2 * N)])
    gr = tuple([(r[i], gr[i], gr[i]) for i in range(2 * N)])
    bl = tuple([(r[i], bl[i], bl[i]) for i in range(2 * N)])

    cdict = {'red': rd, 'green': gr, 'blue': bl}
    return colors.LinearSegmentedColormap('my_colormap', cdict, N)#

# -------------------------------- movies --------------------------------------

def movie_figure(atom_plot, *da, i=0, test=False,
                 overwrite=True, fig_suffix=None, fig_dir=None, 
                 figsize=(4,5), ax_kwargs={}, **atom_kwargs):
    """ Generate a figure, make a plot and store figure. 
    May be distributed across workers
    
    Parameters
    ----------
    atom_plot: method
        Atomic plot method used, should look like:
        def atom_plot(ax, *da, test=..., ax_kwargs=..., **plt_kwargs):
            ...
    *da: xarray DataArray, or other
        Variables containing data that will be plotted
    i: integer, optional
        Index that will be used in filename. Default is 0
    test: boolean, optional
        Flag to test figure generation. Default is False
    overwrite: boolean, optional
        Turn overwrite existing figure files. Default is True
    fig_suffix: str, optional
        Suffix of figure files. Default to da name
    fig_dir: str, optional
        Directory where figures will be stored. 
        Default to figs in user scratch directory
    figsize: tuple, optional
        Figure size. Default to (4,5)
    ax_kwargs: dict, optional
        Dictionnary of keyword arguments passed on the axis
    **atom_kwargs:
        Keyword arguments passed on xarray atomic plot call.
    """
    
    if fig_dir is None:
        fig_dir = os.environ['SCRATCH']+'/figs/'
        
    if test:
        _da = []
        for d in da:
            if isinstance(d, xr.DataArray):
                _da.append(d.isel(time=i))
            else:
                _da.append(d)
        da = tuple(_da)
        
    if fig_suffix is None:
        fig_suffix = '_'.join(d.name for d in da if hasattr(d,'name'))
        
    figname = fig_dir+fig_suffix+'_t%05d' %(i)+'.png'
    
    if not os.path.isfile(figname) or overwrite:
        #
        MPL_LOCK = threading.Lock()
        with MPL_LOCK:
            if not test:
                plt.switch_backend('agg')    
            #
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1,1,1)
            #
            atom_plot(ax, *da, test=test, ax_kwargs=ax_kwargs, **atom_kwargs)
            #
            if not test:
                fig.savefig(figname, bbox_inches = 'tight')
                plt.close()
            #
            m = 1.
    else:
        m = -1.
    if not test:
        return m

def movie_wrapper(atom_plot, client, *da, Nt=None, Nb=None, **kwargs):
    """ wrap figure generation in delayed objects that can be distributed to
    workers. Send by batches of Nb figures.
    
    Parameters
    ----------
    atom_plot: method
        Atomic plot method used, should look like:
        def atom_plot(ax, *da, test=..., ax_kwargs=..., **plt_kwargs):
            ...
    *da: xarray.DataArray's
        Contain data to be plotted
    Nt: int, optional
        Number of figures generated. Default to entire dataset
    Nb: int, optional
        Number of figures to generate in each batch. Default to number of 
        available threads
    **kwargs: optional
        Keyword arguments passed on to movie_figure
    """
    if Nt is None:
        Nt = max([d.time.size for d in da])
    if Nb is None:
        Nb = len(client.nthreads())
    #
    rg = range(0,Nt)
    II = np.array_split(rg,len(rg)/Nb)
    print('%d batches to be done'%len(II))
    #
    delayed_fig = delayed(movie_figure)
    #
    for I in II:
        print(' batch %d-%d'%(I[0],I[-1]))
        values = [delayed_fig(atom_plot, 
                              *(d.isel(time=i) for d in da),
                              i=i, **kwargs) 
                  for i in I]
        futures = client.compute(values)
        results = client.gather(futures)