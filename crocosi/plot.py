import os

import numpy as np
from dask import delayed
import threading

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

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

# -------------------------------- movies --------------------------------------

def movie_figure(da, atom_plot, i=0, test=False,
                 overwrite=True, fig_suffix=None, fig_dir=None, 
                 figsize=(4,5), ax_kwargs={}, **plt_kwargs):
    """ Generate a figure, make a plot and store figure. 
    May be distributed across workers
    
    Parameters
    ----------
    da: xarray DataArray
        Variable containing data that will be plotted
    atom_plot: method
        Atomic plot method used, should look like:
        def atom_plot(da, ax, ax_kwargs, **plt_kwargs):
            ...
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
    **plt_kwargs:
        Keyword arguments passed on xarray plot call.
    """
    
    if fig_dir is None:
        fig_dir = os.environ['SCRATCH']+'/figs/'
        
    if test:
        da = da.isel(time=i)
        
    if fig_suffix is None:
        fig_suffix = da.name
        
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
            atom_plot(da, ax, ax_kwargs, **plt_kwargs)
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

def movie_wrapper(da, atom_plot, client, Nt=None, Nb=None, **get_kwargs):
    """ wrap figure generation in delayed objects that can be distributed to
    workers. Send by batches of Nb figures.
    
    Parameters
    ----------
    
    da: Dask array
        contain data to be plotted.
        
    Nt: int, optional
        Number of figures generated. Default to entire dataset.

    Nb: int, optional
        Number of figures to generate in each batch. Default to number of 
        available threads.
        
    **kwargs: optional
        Keyword arguments passed on to get_fig.
    
    """
    if Nt is None:
        Nt = da.time.size
    if Nb is None:
        Nb = len(client.nthreads())
    #
    rg = range(0,Nt)
    II = np.array_split(rg,len(rg)/Nb)
    #
    print('%d batches to be done'%len(II))
    for I in II:
        print(' batch %d-%d'%(I[0],I[-1]))
        values = [delayed(movie_figure)(da.isel(time=i), atom_plot, 
                                        i, **get_kwargs) for i in I]
        futures = client.compute(values)
        results = client.gather(futures)