
    #!/usr/bin/python
# -*- coding:Utf-8 -*-

import sys
import glob
import os
import numpy as np
import pandas as pd
#import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.collections as col
import numpy.ma as ma


def plot_colorline(x,y,c):
    """
    Plots XY Plot of given trajectory, with color as a function of c
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(c.min(), c.max())
    lc = LineCollection(segments, cmap='RdYlBu', norm=norm)
    # Set the values used for colormapping
    lc.set_array(c)
    lc.set_linewidth(2)

    ax = plt.gca()
    line = ax.add_collection(lc)
    fig.colorbar(line, ax=ax)
    # ax.set_xlim(x.min(), x.max())
    # ax.set_ylim(y.min(), y.max())
    ax.set_xlim(1, 256)
    ax.set_ylim(1, 720)
    # ax.set_xlim(1, 128)
    # ax.set_ylim(1, 360)

    return

def get_floats(dpath):
#
# load float data, written by wrt_floats_mpi.F
#
    # float files in directory and sub-directories t* of dpath
    files = glob.glob(dpath+'/float.????')
    # Files in subdirectories t1,t2,...
    if len(files) == 0:
        files=[]
        tpath = glob.glob(dpath+'/t[1-9]*/')
        for idir in range(len(tpath)):
            dirfiles = glob.glob(tpath[idir]+'/float.????')
            files.append(dirfiles)
        # flatten list
        files = [item for sublist in files for item in sublist]

    # read every float files in a pandas dataframe
    df = pd.DataFrame()
    if len(files)>0: 
        for file in files:
                tmp = pd.read_csv(file, delim_whitespace=True, names=['float','time','xgrid','ygrid','zgrid','depth','temp'])
                df = df.append(tmp)

    # sort the dataframe on float number and time
    df = df.sort_values(by = ['float', 'time'])

    # find how many floats in the dataframe
    nfloats = df.float.nunique()
    print('Number of floats:',nfloats)

    return nfloats,df


# dpath = "/home1/dunree/slgentil/models/croco/croco/CONFIGS/Run_JETN/t1"
# dpath = "/home/datawork-lops-osi/slgentil/croco/jetn/jet_cfg1_wp75_4km_1500a2500j_float/t1"
# dpath = "/home/datawork-lops-osi/slgentil/croco/jetn/jet_cfg1_wp75_4km_1500a2500j_float_z5/t1"
dpath = "/home1/scratch/slgentil/jet_cfg1_wp75_4km_1500a2000j_floats_lev50/t1/"

nfloats,df = get_floats(dpath)

# for iflt in range(nfloats):
# for iflt in range(10,900,100):
for iflt in range(210,211):

    # Get trajectory for float number iflt
    xgrid = df[df['float'] == iflt+1]['xgrid'].values
    ygrid = df[df['float'] == iflt+1]['ygrid'].values
    depth = df[df['float'] == iflt+1]['depth'].values
    temp = df[df['float'] == iflt+1]['temp'].values
    time = df[df['float'] == iflt+1]['time'].values

    # insert None point if float jump from East/West
    for index in range(xgrid.size-1):
        if abs(xgrid[index+1] - xgrid[index]) > 250 :
            xgrid = ma.masked_invalid(np.insert(xgrid,index+1,np.NaN))
            ygrid = ma.masked_invalid(np.insert(ygrid,index+1,np.NaN))
            depth = ma.masked_invalid(np.insert(depth,index+1,np.NaN))
            temp = ma.masked_invalid(np.insert(temp,index+1,np.NaN))
            time = ma.masked_invalid(np.insert(time,index+1,np.NaN))

    fig=plt.figure()

    # plot 2D line(x,y)
    # plt.plot(xgrid[::24],ygrid[::24])
    # plt.axis([1, 128, 1, 360]) 


    # plot 2D line (x,y) colored with temp
    inctime=1
    plot_colorline(xgrid[::inctime],ygrid[::inctime],temp[::inctime])

    # plot 3D line (x,y,z)
    # ax = fig.gca(projection='3d')
    # ax.plot(xgrid,ygrid,depth)

    title = 'float {:4d}, initial position ({:6.1f},{:6.1f}),time {:6.1f} to {:6.1f} days'.format(int(iflt),\
            xgrid[0],ygrid[0],time[0],time[-1])
    plt.title(title)

plt.show()
