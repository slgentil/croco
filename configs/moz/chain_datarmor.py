#!/usr/bin/python
# -*- coding:Utf-8 -*-

"""
  script python pour lancer une simulation en plusieurs chaînages sur caparmor
 
  [syntaxe] : romsc_chain workdir nbchain elaptim resolution jobname
      arguments:
          workdir : répertoire qui sera créé sous $WORK
          nb_chain : nombre de chainages
          elap_time : temps elapsed pour le job HH:MM:SS
          resolution : 128, 256, 512 ou 1024
          jobname : nom générique des batchs
"""

import sys
import os
import shutil


# check number of arguments

if  len(sys.argv) < 7 :
    print '--------------------------------------------------------------------------'
    print '[syntaxe] : romsc_chain workdir nbchain elaptim resolution jobname restart' 
    print '--------------------------------------------------------------------------'
    print 'workdir : répertoire qui sera créé sous $WORK'
    print 'nbchain : nombre de chainages'
    print 'elaptim : temps elapsed pour chaque chainage HH:MM:SS  '
    print 'resolution : resolution 256, 512, 1024, 1152 '
    print 'jobname : nom générique des batchs'   
    print 'restart : 0 (initial) or 1 (restart)'
    quit()

workdir=sys.argv[1]
nbchain=int(sys.argv[2])
elaptim=sys.argv[3]
resolution=int(sys.argv[4])
jobname=sys.argv[5]
restart=int(sys.argv[6])
prefix='moz'

if resolution==128:
    nbproc_roms=2*8
    # nbproc_xios=0
    nbproc_xios=1
elif resolution==256:
    # nbproc_roms=2*16
    nbproc_roms=2*8
    # nbproc_xios=0
    nbproc_xios=1
elif resolution==512:
    nbproc_roms=16*32
    nbproc_xios=2
elif resolution==768:
    nbproc_roms=12*64
    # nbproc_xios=0
    nbproc_xios=4
elif resolution==1024:
    # nbproc_roms=16*32
    # nbproc_xios=16
    nbproc_roms=8*32
    nbproc_xios=16
elif resolution==1152:
    nbproc_roms=12*64
    nbproc_xios=4
    # nbproc_xios=10
else:
    print 'resolution not implemented yet'
    quit()
    
nb_cores = nbproc_roms+nbproc_xios    
nb_nodes = int((nb_cores-1)/28)+1
    
#==============================================================================
    
# Création du répertoire workdir dans /work

startdir=os.getcwd()
HOME = os.getenv('HOME')
USER = os.getenv('USER')
#WORK = os.getenv('DATAWORK')
WORK = os.getenv('SCRATCH')
RPATH = WORK+'/'+workdir
if os.path.exists(RPATH) :    
    os.system('rm -Rf '+RPATH)
os.mkdir(RPATH)
os.chdir(RPATH)

#==============================================================================

# Création des sous_répertoire t0,t1,..tn dans workdir

for t in range(0,nbchain+1):
    tdir='t'+str(t)
    os.mkdir(tdir)
    # os.mkdir(tdir+'/CROCO_FILES')

    if t != 0 :
     
        # recopie des exécutables et input files dans chaque répertoire t
        shutil.copy(startdir+'/croco',tdir)
        shutil.copy(startdir+'/croco.in',tdir)
        # shutil.copy(startdir+'/CROCO_FILES/'+prefix+'_grd.nc',tdir+'/CROCO_FILES')
        # shutil.copy(startdir+'/CROCO_FILES/'+prefix+'_bry.nc',tdir+'/CROCO_FILES')
        # shutil.copy(startdir+'/CROCO_FILES/'+prefix+'_ini.nc',tdir+'/CROCO_FILES')
        shutil.copy(startdir+'/field_def.xml',tdir) 
        shutil.copy(startdir+'/domain_def.xml',tdir) 
        shutil.copy(startdir+'/axis_def.xml',tdir)
        shutil.copy(startdir+'/grid_def.xml',tdir)
        shutil.copy(startdir+'/iodef.xml',tdir)   
#       shutil.copy(startdir+'/floats.in',tdir)
#       shutil.copy('/appli/NETCDF/4.3.3.1-mpt-intel/XIOS-rsync/bin/xios_server.exe',tdir)
        shutil.copy(HOME+'/lib/XIOS/bin/xios_server.exe',tdir)

#==============================================================================

# Copie un backup des fichiers sources dans WORK

if (os.path.exists('backup') == False):
    os.mkdir('backup')
    cmd='ls '+startdir+"/*.[Fh]"
    listfiles = os.popen(cmd).readlines() 
    for ifile in (listfiles): 
      shutil.copy(ifile[:-1],'backup')
        
#==============================================================================

# change NRREC in file croco.in depending on restart  

for t in range(1,nbchain+1):
    tdir='t'+str(t)
    os.chdir(RPATH+'/'+tdir)
    
    if t == 1: 
       with open('croco.in','r') as f:
           lines = f.readlines()
           for index, line in enumerate(lines):
               if 'NRREC' in line:
                  if restart==0 :
                     lines[index+1]='         0\n'
                  else:
                     lines[index+1]='         1\n'
                  break
    else:
        with open('croco.in','r') as f:
           lines = f.readlines()
           for index, line in enumerate(lines):
               if 'NRREC' in line:
                  lines[index+1]='         1\n' 
                  break  
    f.close()   
    f = open('croco.in','w')
    f.writelines(lines)

#==============================================================================


# Création du script jobname* dans chaque workdir/t*

for t in range(1,nbchain+1):
    tdir='t'+str(t)   
    os.chdir(RPATH+'/'+tdir)
    
    fo = open('job_datarmor','w')
    fo.write('#!/bin/csh\n')
    fo.write('#PBS -N '+jobname+str(t)+'\n')
    fo.write('#PBS -q mpi\n')
    fo.write('#PBS -l select='+str(nb_nodes)+':ncpus=28:mpiprocs=28:mem=120G\n')
    fo.write('#PBS -l walltime='+elaptim+'\n')
    fo.write('\n') 
    fo.write('# get the path for command module,mpirun\n')
    fo.write('source /usr/share/Modules/3.2.10/init/csh\n') 
    fo.write('module purge\n')  
    fo.write('module load   NETCDF/4.3.3.1-mpt-intel2016  intel-cmkl-16/16.0.4.258\n')  
    fo.write('\n')    
    fo.write('# cd to the directory you submitted your job\n')
    fo.write('cd $PBS_O_WORKDIR\n')
    fo.write('\n')     
    fo.write('date\n') 
    if t!=1 or restart==1:
       fo.write('cp ../t'+str(t-1)+'/'+prefix+'_grd.nc .\n')  
       fo.write('cp ../t'+str(t-1)+'/'+prefix+'_ini.nc .\n')    
       fo.write('cp ../t'+str(t-1)+'/'+prefix+'_bry.nc .\n')    
    # fo.write('time $MPI_LAUNCH -n '+str(nbproc_roms)+' roms >& output.mpi\n')   
    fo.write('time $MPI_LAUNCH -n '+str(nbproc_roms)+' croco : -n '+str(nbproc_xios)+ ' xios_server.exe >& output.mpi\n')  

    fo.close()

#==============================================================================

# Commandes pour lancer le batch

for t in range(1,nbchain+1):    
    tdir='t'+str(t)   
    os.chdir(RPATH+'/'+tdir)
    
    if (t == 1):
        commande='qsub -h job_datarmor'  
    else:
        f=os.popen('qselect -N '+jobname+str(t-1)+' -u '+USER)
        numjob=f.read()
        commande='qsub -W depend=afterany:'+str(numjob[:-1])+' job_datarmor'
    os.system(commande)
os.chdir(RPATH+'/t1')
f=os.popen('qselect -N '+jobname+str(1)+' -u '+USER)
numjob=f.read()
commande='qrls '+numjob[:-1]     

if restart==1:
   print 'Put the files (grd,ini,bry) in '+RPATH+'/t0'
print 'Change directory: cd '+RPATH+'/t1'
print 'Run commande : '+commande
