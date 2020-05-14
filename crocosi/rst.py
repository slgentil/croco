#!/usr/bin/python
# -*- coding:Utf-8 -*-

import sys
import os
import time

import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset


# ------------------------------ join ------------------------------------------

def join(run, tdir):
    """ join rst files
    
    Parameters
    ----------
    run: str
        Path to simulation
    tdir: str
        directory where restart files are, e.g. 't5/'
    """

    start = time.time()
    #
    
    # input dir
    inDir = os.path.join(run, tdir)
    inFname = os.popen('ls '+inDir+'jetn_rst.*.nc').readlines()

    # output file
    outDir = os.path.join(run,'data_rst')
    if not os.path.exists(outDir) :
        os.mkdir(outDir)
    outFname = os.path.join(outDir,'jetn_rst.nc')

    # Ouverture d'un fichier Netcdf
    nc = Dataset(inFname[0][:-1], 'r') 

    #---------------------------------------------------
    # get global variables
    # get time
    scrum_time = nc.variables['scrum_time'][:]
    scrum_time_units = nc.variables['scrum_time'].units
    # get time steps
    time_step = nc.variables['time_step'][:]
    # get partition
    partition = nc.getncattr('partition')
    np_xi = partition[2]
    np_eta = partition[3]

    #---------------------------------------------------
    # get grid dimensions-lpython2.6 -lpthread -lm -lutil -ldl
    # local dimensions, valid for an interior grid
    xi_rho = len(nc.dimensions['xi_rho'])
    eta_rho = len(nc.dimensions['eta_rho'])
    s_rho = len(nc.dimensions['s_rho'])
    # compute output global sizes
    out_Lm = (xi_rho - 1) * np_xi
    out_Mm = (eta_rho - 1) * np_eta
    out_N = s_rho

    nc.close()

    # create output file-lpython2.6 -lpthread -lm -lutil -ldl
    ncid = Dataset(outFname, 'w', format='NETCDF4')

    # define dimensions
    ncid.createDimension('xi_rho', out_Lm + 2)
    ncid.createDimension('xi_u', out_Lm + 1)
    ncid.createDimension('eta_rho', out_Mm + 2)
    ncid.createDimension('eta_v', out_Mm + 1)
    ncid.createDimension('s_rho', out_N)
    ncid.createDimension('s_w', out_N + 1)
    out_time = ncid.createDimension('time', None)
    ncid.createDimension('auxil', 4)

    # define variables
    out_scrum_time = ncid.createVariable('scrum_time', 'f8', ('time',))
    out_scrum_time.units = scrum_time_units
    out_timestep = ncid.createVariable('time_step', 'i4', ('auxil',))

    out_zeta = ncid.createVariable('zeta', 'f8', ('time', 'eta_rho', 'xi_rho',))
    out_ubar = ncid.createVariable('ubar', 'f8', ('time', 'eta_rho', 'xi_u',))
    out_vbar = ncid.createVariable('vbar', 'f8', ('time', 'eta_v', 'xi_rho',))

    out_temp = ncid.createVariable('temp', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho',))
    out_u = ncid.createVariable('u', 'f8', ('time', 's_rho', 'eta_rho', 'xi_u',))
    out_v = ncid.createVariable('v', 'f8', ('time', 's_rho', 'eta_v', 'xi_rho',))

    # global variable for ncjoin
    ncid.partition = partition

    # store time related values
    out_scrum_time[:] = scrum_time
    out_timestep[:] = time_step

    # loop around out grid subblocks
    j0r = 0
    j0v = 0
    for jj in range(0, np_eta):

        i0r = 0
        i0u = 0 
        for ii in range(0, np_xi):

           mynode = (jj * np_xi) + ii

           # open mynode input file
           print('Tile ' + str(mynode))
           nc = Dataset(inFname[mynode][:-1], 'r') 

           # get grid dimensions
           xi_rho = len(nc.dimensions['xi_rho'])
           xi_u = len(nc.dimensions['xi_u'])
           eta_rho = len(nc.dimensions['eta_rho'])
           eta_v = len(nc.dimensions['eta_v'])
           s_rho = len(nc.dimensions['s_rho'])

           # load and store in output file

           # zeta
           zeta = nc.variables['zeta'][:]
           out_zeta[:, j0r:j0r + eta_rho, i0r:i0r + xi_rho] = zeta
           # ubar
           ubar = nc.variables['ubar'][:]
           out_ubar[:, j0r:j0r + eta_rho, i0u:i0u + xi_u] = ubar
           # vbar
           vbar = nc.variables['vbar'][:]
           out_vbar[:, j0v:j0v + eta_v, i0r:i0r + xi_rho] = vbar
           # u
           u = nc.variables['u'][:]
           out_u[:, :, j0r:j0r + eta_rho, i0u:i0u + xi_u] = u
           # v
           v = nc.variables['v'][:]
           out_v[:, :, j0v:j0v + eta_v, i0r:i0r + xi_rho] = v
           # temp
           temp = nc.variables['temp'][:]
           out_temp[:, :, j0r:j0r + eta_rho, i0r:i0r + xi_rho] = temp

           nc.close()

           # update indices
           i0r = i0r + xi_rho
           i0u = i0u + xi_u

        # update indices
        j0r = j0r + eta_rho
        j0v = j0v + eta_v

    # close new tile file
    ncid.close()
    print(outFname + ' created')
    print('Elapsed time is ', str(time.time() - start))

# ------------------------------ interp  ---------------------------------------

#
# Interpolate roms rst data to double resolution
#
# assumes x periodicity
# assumes one time record in rst
# 
# the grid arrangement is not obvious, choices have to be made
#
# variables horizontal global sizes:
#    zeta, temp    Lm+2 x Mm+2
#    ubar, u       Lm+1 x Mm+2
#    vbar, v       Lm+2 x Mm+1
#
# Inner block:
#    zeta, temp    Lm/np_xi   x  Mm/np_eta
#    ubar, u       Lm/np_xi   x  Mm/np_eta
#    vbar, v       Lm/np_xi   x  Mm/np_eta
#
# West bdy
#    zeta, temp    Lm/np_xi+1
#    ubar, u       Lm/np_xi
#    vbar, v       Lm/np_xi+1
# East bdy
#    zeta, temp    Lm/np_xi+1
#    ubar, u       Lm/np_xi+1
#    vbar, v       Lm/np_xi+1
#
# South bdy
#    zeta, temp    Mm/np_eta+1
#    ubar, u       Mm/np_eta+1
#    vbar, v       Mm/np_eta
# North bdy
#    zeta, temp    Mm/np_eta+1
#    ubar, u       Mm/np_eta+1
#    vbar, v       Mm/np_eta+1
#
# Periodicity:
#    eta(1)=eta(Lm+1)   and    eta(2)=eta(Lm+2)
#      u(1)=u(Lm+1)
#      u(1) is in between eta(1) and eta(2) !!!! check ????
#

def interp(run, out_np_xi, out_np_eta):

    print('NP_XI = ',out_np_xi,' NP_ETA = ',out_np_eta)
        
    start = time.time()

    # input dir
    inDir = os.path.join(run,'data_rst')
    inFname = os.path.join(inDir,'jetn_rst.nc')

    # output file
    outDir = inDir
    outFname = os.path.join(outDir,'jetn_rst.')
    #outFname='./data_rst/jetn_rst.'

    # interpolation method
    interp_method='linear';

    # get global grid dimensions
    nc = Dataset(inFname, 'r') 
    in_xi_rho = len(nc.dimensions['xi_rho'])
    in_eta_rho = len(nc.dimensions['eta_rho'])
    in_s_rho = len(nc.dimensions['s_rho'])

    # store
    in_LLm=in_xi_rho-2
    in_MMm=in_eta_rho-2
    in_N=in_s_rho

    # compute output global Lm dimensions
    out_LLm=in_LLm*2
    out_MMm=in_MMm*2
    out_Lm=int(out_LLm/out_np_xi)
    out_Mm=int(out_MMm/out_np_eta)
    out_N=in_N*2

    # compute global grid
    # zonal and meridional direction are different:
    # zonal direction: some rho points are collocated
    # meridional direction: some v points are collocated
    in_xr=np.linspace(-0.5,in_LLm+0.5,in_LLm+2)
    in_xu=np.linspace(0,in_LLm,in_LLm+1)
    in_yr=np.linspace(-0.5,in_MMm+0.5,in_MMm+2)
    in_yv=np.linspace(0,in_MMm,in_MMm+1)
    in_sw=np.linspace(-1,0,in_N+1)
    in_sr=np.linspace(-0.5,in_N+0.5,in_N+2)


    out_xr=(np.linspace(-0.25,in_LLm+0.25,out_LLm+2))
    out_xu=np.linspace(0,in_LLm,out_LLm+1)
    out_yr=(np.linspace(-0.25,in_MMm+0.25,out_MMm+2))
    out_yv=np.linspace(0,in_MMm,out_MMm+1)
    out_sw=np.linspace(-1,0,out_N+1)
    out_sr=np.linspace(0.25,in_N-0.25,out_N)

    # get time
    scrum_time = nc.variables['scrum_time'][:]
    scrum_time_units = nc.variables['scrum_time'].units
    # get time steps
    time_step = nc.variables['time_step'][:]

    # loop around out grid subblocks

    for jj in range(0,out_np_eta):

        for ii in range(0,out_np_xi):
            
            mynode = (jj * out_np_xi) + ii
            print('Tile = {} / {}'.format(mynode, out_np_eta*out_np_xi))
            
            # calculate offset of out subblock in global domain
            iminmpiu=(ii*out_Lm)
            imaxmpiu=iminmpiu+out_Lm-1  
            iminmpir=1+(ii*out_Lm)
            imaxmpir=iminmpir+out_Lm-1          
            if (ii==0): 
                iminmpir = iminmpir-1
            if (ii==out_np_xi-1):
                imaxmpiu=imaxmpiu+1
                imaxmpir=imaxmpir+1
                
            jminmpiv=(jj*out_Mm)
            jmaxmpiv=jminmpiv+out_Mm-1
            jminmpir=1+(jj*out_Mm)
            jmaxmpir=jminmpir+out_Mm-1   
            if (jj==0): 
                jminmpir = jminmpir-1
            if (jj==out_np_eta-1):
                jmaxmpiv=jmaxmpiv+1
                jmaxmpir=jmaxmpir+1        
                   
            # calculate lenghts of subblock 
            out_xi_rho=imaxmpir-iminmpir+1
            out_xi_u=imaxmpiu-iminmpiu+1
            out_eta_rho=jmaxmpir-jminmpir+1
            out_eta_v=jmaxmpiv-jminmpiv+1
     
            # create netcdf file
            fmt='{:0'+str(np.int(np.floor(np.log10(out_np_xi*out_np_eta))+1))+'d}'
            ncid = Dataset(outFname+fmt.format(mynode)+'.nc', 'w', format='NETCDF4')
            
            # define dimensions
            ncid.createDimension('xi_rho', out_xi_rho)
            ncid.createDimension('xi_u', out_xi_u)
            ncid.createDimension('eta_rho', out_eta_rho)
            ncid.createDimension('eta_v', out_eta_v)
            ncid.createDimension('s_rho', out_N)
            ncid.createDimension('s_w', out_N+1)
            ncid.createDimension('time', None)
            ncid.createDimension('auxil', 4)  
             
            #  define variables
            out_scrum_time = ncid.createVariable('scrum_time','f8',('time',))
            out_scrum_time.units = scrum_time_units
            out_timestep = ncid.createVariable('time_step','i4',('auxil',))
            
            out_zeta = ncid.createVariable('zeta','f8',('time','eta_rho','xi_rho',))
            out_ubar = ncid.createVariable('ubar','f8',('time','eta_rho','xi_u',))
            out_vbar = ncid.createVariable('vbar','f8',('time','eta_v','xi_rho',))
            
            out_temp = ncid.createVariable('temp','f8',('time','s_rho','eta_rho','xi_rho',))
            out_u = ncid.createVariable('u','f8',('time','s_rho','eta_rho','xi_u',))
            out_v = ncid.createVariable('v','f8',('time','s_rho','eta_v','xi_rho',))        

            # global variable for ncjoin
            ncid.partition = [ii,jj,out_np_xi,out_np_eta]

            #  store time related values
            
            # store time related values
            out_scrum_time[:] = scrum_time
            out_timestep[:] = time_step        

            # Load/interp/store 2D/3D fields
            
            #zeta
            zeta = nc.variables['zeta'][:]
            i0 = np.where(in_xr<=out_xr[iminmpir])[-1][-1]
            i1 = np.where(in_xr>=out_xr[imaxmpir])[0][0]
            j0 = np.where(in_yr<=out_yr[jminmpir])[-1][-1]
            j1 = np.where(in_yr>=out_yr[jmaxmpir])[0][0] 
            out_zeta[0,:,:] = interp2d( in_yr[j0:j1+1], in_xr[i0:i1+1], zeta[0,j0:j1+1,i0:i1+1],
                              out_yr[jminmpir:jmaxmpir+1], out_xr[iminmpir:imaxmpir+1])
            
            #temp
            temp = nc.variables['temp'][:]
            #deal with lower and upper bdy: dTdz=0
            temp=np.concatenate((temp[0,0:1,j0:j1+1,i0:i1+1], temp[0,:,j0:j1+1,i0:i1+1],temp[0,-1:,j0:j1+1,i0:i1+1]),axis=0)
            out_temp[0,:,:,:] = interp3d( in_sr, in_yr[j0:j1+1], in_xr[i0:i1+1], temp,
                              out_sr, out_yr[jminmpir:jmaxmpir+1], out_xr[iminmpir:imaxmpir+1])  
            
            #ubar
            ubar = nc.variables['ubar'][:]
            i0 = np.where(in_xu<=out_xu[iminmpiu])[-1][-1]
            i1 = np.where(in_xu>=out_xu[imaxmpiu])[0][0]
            j0 = np.where(in_yr<=out_yr[jminmpir])[-1][-1]
            j1 = np.where(in_yr>=out_yr[jmaxmpir])[0][0] 
            out_ubar[0,:,:] = interp2d( in_yr[j0:j1+1], in_xu[i0:i1+1], ubar[0,j0:j1+1,i0:i1+1],
                              out_yr[jminmpir:jmaxmpir+1], out_xu[iminmpiu:imaxmpiu+1]) 
            
            #u
            u = nc.variables['u'][:]
            #deal with lower and upper bdy: dudz=0
            u=np.concatenate((u[0,0:1,j0:j1+1,i0:i1+1], u[0,:,j0:j1+1,i0:i1+1],u[0,-1:,j0:j1+1,i0:i1+1]),axis=0);
            out_u[0,:,:,:] = interp3d( in_sr, in_yr[j0:j1+1], in_xu[i0:i1+1], u,
                              out_sr, out_yr[jminmpir:jmaxmpir+1], out_xu[iminmpiu:imaxmpiu+1])      
            
            #vbar
            vbar = nc.variables['vbar'][:]
            i0 = np.where(in_xr<=out_xr[iminmpir])[-1][-1]
            i1 = np.where(in_xr>=out_xr[imaxmpir])[0][0]
            j0 = np.where(in_yv<=out_yv[jminmpiv])[-1][-1]
            j1 = np.where(in_yv>=out_yv[jmaxmpiv])[0][0] 
            out_vbar[0,:,:] = interp2d( in_yv[j0:j1+1], in_xr[i0:i1+1], vbar[0,j0:j1+1,i0:i1+1],
                              out_yv[jminmpiv:jmaxmpiv+1], out_xr[iminmpir:imaxmpir+1])   
            
            #v
            v = nc.variables['v'][:]
            #deal with lower and upper bdy: dvdz=0
            v=np.concatenate((v[0,0:1,j0:j1+1,i0:i1+1], v[0,:,j0:j1+1,i0:i1+1],v[0,-1:,j0:j1+1,i0:i1+1]),axis=0);
            out_v[0,:,:,:] = interp3d( in_sr, in_yv[j0:j1+1], in_xr[i0:i1+1], v,
                              out_sr, out_yv[jminmpiv:jmaxmpiv+1], out_xr[iminmpir:imaxmpir+1]) 
            
            ncid.close()

    # close input file
    nc.close()
    print('./data_rst/jetn_rst.[0-' + str(out_np_xi*out_np_eta-1) + '].nc created')
    print('Elapsed time is ',str(time.time() - start))

# ------------------------------ interp mpi  -----------------------------------

#
# Interpolate roms rst data to double resolution
#
# assumes x periodicity
# assumes one time record in rst
# 
# the grid arrangement is not obvious, choices have to be made
#
# variables horizontal global sizes:
#    zeta, temp    Lm+2 x Mm+2
#    ubar, u       Lm+1 x Mm+2
#    vbar, v       Lm+2 x Mm+1
#
# Inner block:
#    zeta, temp    Lm/np_xi   x  Mm/np_eta
#    ubar, u       Lm/np_xi   x  Mm/np_eta
#    vbar, v       Lm/np_xi   x  Mm/np_eta
#
# West bdy
#    zeta, temp    Lm/np_xi+1
#    ubar, u       Lm/np_xi
#    vbar, v       Lm/np_xi+1
# East bdy
#    zeta, temp    Lm/np_xi+1
#    ubar, u       Lm/np_xi+1
#    vbar, v       Lm/np_xi+1
#
# South bdy
#    zeta, temp    Mm/np_eta+1
#    ubar, u       Mm/np_eta+1
#    vbar, v       Mm/np_eta
# North bdy
#    zeta, temp    Mm/np_eta+1
#    ubar, u       Mm/np_eta+1
#    vbar, v       Mm/np_eta+1
#
# Periodicity:
#    eta(1)=eta(Lm+1)   and    eta(2)=eta(Lm+2)
#      u(1)=u(Lm+1)
#      u(1) is in between eta(1) and eta(2) !!!! check ????
#

def interp_mpi(out_np_xi, out_np_eta):
    
    from mpi4py import MPI

    print('NP_XI = ',out_np_xi,' NP_ETA = ',out_np_eta)

    start = time.time()

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    myrank = comm.Get_rank()

    # input variable
    # restart file joined with rst_join.py
    inFname='./data_rst/jetn_rst.nc'

    # output file
    outFname='./data_rst/jetn_rst.'

    # interpolation method
    interp_method='linear';

    # get global grid dimensions
    nc = Dataset(inFname, 'r') 
    in_xi_rho = len(nc.dimensions['xi_rho'])
    in_eta_rho = len(nc.dimensions['eta_rho'])
    in_s_rho = len(nc.dimensions['s_rho'])

    # store
    in_LLm=in_xi_rho-2
    in_MMm=in_eta_rho-2
    in_N=in_s_rho

    # compute output global Lm dimensions
    out_LLm=in_LLm*2
    out_MMm=in_MMm*2
    out_Lm=out_LLm/out_np_xi
    out_Mm=out_MMm/out_np_eta
    out_N=in_N*2

    # compute global grid
    # zonal and meridional direction are different:
    # zonal direction: some rho points are collocated
    # meridional direction: some v points are collocated
    in_xr=np.linspace(-0.5,in_LLm+0.5,in_LLm+2)
    in_xu=np.linspace(0,in_LLm,in_LLm+1)
    in_yr=np.linspace(-0.5,in_MMm+0.5,in_MMm+2)
    in_yv=np.linspace(0,in_MMm,in_MMm+1)
    in_sw=np.linspace(-1,0,in_N+1)
    in_sr=np.linspace(-0.5,in_N+0.5,in_N+2)


    out_xr=(np.linspace(-0.25,in_LLm+0.25,out_LLm+2))
    out_xu=np.linspace(0,in_LLm,out_LLm+1)
    out_yr=(np.linspace(-0.25,in_MMm+0.25,out_MMm+2))
    out_yv=np.linspace(0,in_MMm,out_MMm+1)
    out_sw=np.linspace(-1,0,out_N+1)
    out_sr=np.linspace(0.25,in_N-0.25,out_N)

    # get time
    scrum_time = nc.variables['scrum_time'][:]
    scrum_time_units = nc.variables['scrum_time'].units
    # get time steps
    time_step = nc.variables['time_step'][:]

    # loop around out grid subblocks

    for jj in range(0,out_np_eta):

        for ii in range(0+myrank,out_np_xi,size):
            
            mynode = (jj * out_np_xi) + ii
    #        print 'Tile =',mynode
    #        if (np.mod(mynode+1,size)==myrank):
            print('Tile =',mynode,myrank)

            # calculate offset of out subblock in global domain
            iminmpiu=(ii*out_Lm)
            imaxmpiu=iminmpiu+out_Lm-1  
            iminmpir=1+(ii*out_Lm)
            imaxmpir=iminmpir+out_Lm-1          
            if (ii==0): 
                iminmpir = iminmpir-1
            if (ii==out_np_xi-1):
                imaxmpiu=imaxmpiu+1
                imaxmpir=imaxmpir+1
                
            jminmpiv=(jj*out_Mm)
            jmaxmpiv=jminmpiv+out_Mm-1
            jminmpir=1+(jj*out_Mm)
            jmaxmpir=jminmpir+out_Mm-1   
            if (jj==0): 
                jminmpir = jminmpir-1
            if (jj==out_np_eta-1):
                jmaxmpiv=jmaxmpiv+1
                jmaxmpir=jmaxmpir+1        
                   
            # calculate lenghts of subblock 
            out_xi_rho=imaxmpir-iminmpir+1
            out_xi_u=imaxmpiu-iminmpiu+1
            out_eta_rho=jmaxmpir-jminmpir+1
            out_eta_v=jmaxmpiv-jminmpiv+1
     
            # create netcdf file
            fmt='{:0'+str(np.int(np.floor(np.log10(out_np_xi*out_np_eta))+1))+'d}'
            ncid = Dataset(outFname+fmt.format(mynode)+'.nc', 'w', format='NETCDF4')
            
            # define dimensions
            ncid.createDimension('xi_rho', out_xi_rho)
            ncid.createDimension('xi_u', out_xi_u)
            ncid.createDimension('eta_rho', out_eta_rho)
            ncid.createDimension('eta_v', out_eta_v)
            ncid.createDimension('s_rho', out_N)
            ncid.createDimension('s_w', out_N+1)
            ncid.createDimension('time', None)
            ncid.createDimension('auxil', 4)  
             
            #  define variables
            out_scrum_time = ncid.createVariable('scrum_time','f8',('time',))
            out_scrum_time.units = scrum_time_units
            out_timestep = ncid.createVariable('time_step','i4',('auxil',))
            
            out_zeta = ncid.createVariable('zeta','f8',('time','eta_rho','xi_rho',))
            out_ubar = ncid.createVariable('ubar','f8',('time','eta_rho','xi_u',))
            out_vbar = ncid.createVariable('vbar','f8',('time','eta_v','xi_rho',))
            
            out_temp = ncid.createVariable('temp','f8',('time','s_rho','eta_rho','xi_rho',))
            out_u = ncid.createVariable('u','f8',('time','s_rho','eta_rho','xi_u',))
            out_v = ncid.createVariable('v','f8',('time','s_rho','eta_v','xi_rho',))        

            # global variable for ncjoin
            ncid.partition = [ii,jj,out_np_xi,out_np_eta]

            #  store time related values
            
            # store time related values
            out_scrum_time[:] = scrum_time
            out_timestep[:] = time_step        

            # Load/interp/store 2D/3D fields
            
            #zeta
            zeta = nc.variables['zeta'][:]
            i0 = np.where(in_xr<=out_xr[iminmpir])[-1][-1]
            i1 = np.where(in_xr>=out_xr[imaxmpir])[0][0]
            j0 = np.where(in_yr<=out_yr[jminmpir])[-1][-1]
            j1 = np.where(in_yr>=out_yr[jmaxmpir])[0][0] 
            out_zeta[0,:,:] = interp2d( in_yr[j0:j1+1], in_xr[i0:i1+1], zeta[0,j0:j1+1,i0:i1+1],
                              out_yr[jminmpir:jmaxmpir+1], out_xr[iminmpir:imaxmpir+1])
            
            #temp
            temp = nc.variables['temp'][:]
            #deal with lower and upper bdy: dTdz=0
            temp=np.concatenate((temp[0,0:1,j0:j1+1,i0:i1+1], temp[0,:,j0:j1+1,i0:i1+1],temp[0,-1:,j0:j1+1,i0:i1+1]),axis=0)
            out_temp[0,:,:,:] = interp3d( in_sr, in_yr[j0:j1+1], in_xr[i0:i1+1], temp,
                              out_sr, out_yr[jminmpir:jmaxmpir+1], out_xr[iminmpir:imaxmpir+1])  
            
            #ubar
            ubar = nc.variables['ubar'][:]
            i0 = np.where(in_xu<=out_xu[iminmpiu])[-1][-1]
            i1 = np.where(in_xu>=out_xu[imaxmpiu])[0][0]
            j0 = np.where(in_yr<=out_yr[jminmpir])[-1][-1]
            j1 = np.where(in_yr>=out_yr[jmaxmpir])[0][0] 
            out_ubar[0,:,:] = interp2d( in_yr[j0:j1+1], in_xu[i0:i1+1], ubar[0,j0:j1+1,i0:i1+1],
                              out_yr[jminmpir:jmaxmpir+1], out_xu[iminmpiu:imaxmpiu+1]) 
            
            #u
            u = nc.variables['u'][:]
            #deal with lower and upper bdy: dudz=0
            u=np.concatenate((u[0,0:1,j0:j1+1,i0:i1+1], u[0,:,j0:j1+1,i0:i1+1],u[0,-1:,j0:j1+1,i0:i1+1]),axis=0);
            out_u[0,:,:,:] = interp3d( in_sr, in_yr[j0:j1+1], in_xu[i0:i1+1], u,
                              out_sr, out_yr[jminmpir:jmaxmpir+1], out_xu[iminmpiu:imaxmpiu+1])      
            
            #vbar
            vbar = nc.variables['vbar'][:]
            i0 = np.where(in_xr<=out_xr[iminmpir])[-1][-1]
            i1 = np.where(in_xr>=out_xr[imaxmpir])[0][0]
            j0 = np.where(in_yv<=out_yv[jminmpiv])[-1][-1]
            j1 = np.where(in_yv>=out_yv[jmaxmpiv])[0][0] 
            out_vbar[0,:,:] = interp2d( in_yv[j0:j1+1], in_xr[i0:i1+1], vbar[0,j0:j1+1,i0:i1+1],
                              out_yv[jminmpiv:jmaxmpiv+1], out_xr[iminmpir:imaxmpir+1])   
            
            #v
            v = nc.variables['v'][:]
            #deal with lower and upper bdy: dvdz=0
            v=np.concatenate((v[0,0:1,j0:j1+1,i0:i1+1], v[0,:,j0:j1+1,i0:i1+1],v[0,-1:,j0:j1+1,i0:i1+1]),axis=0);
            out_v[0,:,:,:] = interp3d( in_sr, in_yv[j0:j1+1], in_xr[i0:i1+1], v,
                              out_sr, out_yv[jminmpiv:jmaxmpiv+1], out_xr[iminmpir:imaxmpir+1]) 
            
            ncid.close()

    # close input file
    nc.close()
    print('./data_rst/jetn_rst.[0-' + str(out_np_xi*out_np_eta-1) + '].nc created')
    print('Elapsed time is ',str(time.time() - start))



# ------------------------------ utils -----------------------------------------

def interp3d(x, y, z, v, xi, yi, zi, method='linear'):
    """Interpolation on 3-D. x, y, z, xi, yi, zi should be 1-D
    and v.shape == (len(x), len(y), len(z))"""
    q = (x, y, z)
    qi = (xi, yi, zi)
    for j in range(3):
        v = interp1d(q[j], v, axis=j, kind=method)(qi[j])
    return v

def interp2d(x, y, v, xi, yi, method='linear'):
    """Interpolation on 2-D. x, y, xi, yi should be 1-D
    and z.shape == (len(x), len(y))"""
    q = (x, y)
    qi = (xi, yi)
    for j in range(2):
        v = interp1d(q[j], v, axis=j, kind=method)(qi[j])
    return v

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
                      )
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim )
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars