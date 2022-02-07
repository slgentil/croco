#!/usr/bin/python
# -*- coding:Utf-8 -*-

import sys
import os
import shutil

###########################################################################
# workdir : répertoire qui sera créé sous $SCRATCHDIR
# nb_chain : nombre de chainages
# elap_time : Elapsed time limit in seconds
# resolution :  512 ou 1024
# projid : numero de projet gen6130(Claire) ou gen7515(Aurelien)
# restart : 0 (initial) or 1 (restart)
###########################################################################


# test number of arguments

if  len(sys.argv) < 7 :
    print '--------------------------------------------------------------------------'
    print '[syntaxe] : chain_irene.py workdir nbchain elaptim resolution projid restart' 
    print '--------------------------------------------------------------------------'
    print 'workdir : répertoire qui sera créé sous $WORK'
    print 'nbchain : nombre de chainages'
    print 'elaptim : temps elapsed pour chaque chainage (secondes)'
    print 'resolution : resolution 128, 256, 512, 1024, 1620, 3000'
    print 'projid :numero de projet gen12938(Claire)'   
    print 'restart : 0 (initial) or 1 (restart)'
    quit()
    
workdir=sys.argv[1]
nbchain=int(sys.argv[2])
elaptim=sys.argv[3]
resolution=int(sys.argv[4])
projid= sys.argv[5]
restart=int(sys.argv[6])

###########################################################################

# Configuration according to resolution
if resolution==128:
    nbproc_roms=2*8
    nbproc_xios=1
elif resolution==256:
    nbproc_roms=2*8
    nbproc_xios=1
elif resolution==512:
    nbproc_roms=4*32
    nbproc_xios=2
elif resolution==1024:
    nbproc_roms=16*32
    nbproc_xios=64
elif resolution==1620:
    nbproc_roms=10*47
    nbproc_xios=188
elif resolution==3000:
    nbproc_roms=50*100
    nbproc_xios=0
else:
    print 'resolution not implemented yet'
    quit()

nb_cores = nbproc_roms+nbproc_xios    

###########################################################################

# Création du répertoire workdir dans WORK

startdir=os.getcwd()
HOME = os.getenv('HOME')
WORK = os.getenv('CCCWORKDIR')
WORKG = os.getenv('ALL_CCCWORKDIR')
#WORK = os.getenv('CCCSCRATCHDIR')
#WORK = os.getenv('ALL_CCCSCRATCHDIR')
SCRATCHG = os.getenv('ALL_CCCSCRATCHDIR')
RPATH = SCRATCHG+'/'+workdir
if os.path.exists(RPATH) :
    os.system('rm -Rf '+RPATH)
os.mkdir(RPATH)
os.chdir(RPATH)

#==============================================================================

# Création des sous_répertoire t0,t1,..tn dans workdir

for t in range(0,nbchain+1):
    tdir='t'+str(t)
    os.mkdir(tdir)
    if t != 0 :
     
        # recopie des exécutables et input files dans chaque répertoire t
        if os.path.exists(startdir+"/croco"):
            shutil.copy(startdir+'/croco',tdir)
        else:
            print "croco not found in current directory"
        if os.path.exists(startdir+"/croco.in"):
            shutil.copy(startdir+'/croco.in',tdir)
        else:
            print "croco.in not found in current directory"      	
        if os.path.exists(startdir+"/iodef.xml"):
         	shutil.copy(startdir+'/iodef.xml',tdir)
        else:
        	print "iodef.xml not found in current directory"      	
        if os.path.exists(startdir+"/domain_def.xml"):
     	    shutil.copy(startdir+'/domain_def.xml',tdir)
        else:
     	    print "domain_def.xml not found in current directory"      	
        if os.path.exists(startdir+"/domain_def_low.xml"):
     	    shutil.copy(startdir+'/domain_def_low.xml',tdir)
        else:
     	    print "domain_def_low.xml not found in current directory"      	
        if os.path.exists(startdir+"/field_def.xml"):
     	    shutil.copy(startdir+'/field_def.xml',tdir)
        else:
     	    print "field_def.xml not found in current directory"      	
        if os.path.exists(startdir+"/field_def_low.xml"):
     	    shutil.copy(startdir+'/field_def_low.xml',tdir)
        else:
     	    print "field_def_low.xml not found in current directory"      	
        if os.path.exists(startdir+"/axis_def.xml"):
     	    shutil.copy(startdir+'/axis_def.xml',tdir)
        else:
     	    print "axis_def.xml not found in current directory"      	
        if os.path.exists(startdir+"/grid_def.xml"):
     	    shutil.copy(startdir+'/grid_def.xml',tdir)
        else:
     	    print "grid_def.xml not found in current directory"      	
        if os.path.exists(startdir+"/grid_def_low.xml"):
     	    shutil.copy(startdir+'/grid_def_low.xml',tdir)
        else:
     	    print "grid_def_low.xml not found in current directory"      	
        #if os.path.exists(startdir+"/xios_server.exe"):
     	#    shutil.copy(startdir+'/xios_server.exe',tdir)
        #else:
     	#    print "xios_server.exe not found in current directory"

#==============================================================================

# Copie un backup des fichiers sources dans WORK

if (os.path.exists('backup') == False):
    os.mkdir('backup')
    cmd='ls '+startdir+"/*.[Fh]"
    listfiles = os.popen(cmd).readlines()  
    for ifile in (listfiles): 
        shutil.copy(ifile[:-1],'backup')
        
#==============================================================================

# change NRREC in file croco.in.MEDDY depending on restart  

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
# Création du script job_curie dans chaque tdir


for t in range(1,nbchain+1):
    tdir=RPATH+'/t'+str(t)
    os.chdir(tdir)
    fo = open('job_irene','w')
    fo.write('#!/bin/bash\n')
    fo.write('#MSUB -r moz'+str(t)+'\n')
    fo.write('#MSUB -o output.mpi\n')
    fo.write('#MSUB -e log.err\n')
    fo.write('#MSUB -q skylake\n')
    fo.write('#MSUB -A '+projid+'\n')
    fo.write('#MSUB -T '+elaptim+'\n')
    fo.write('#MSUB -n '+str(nb_cores)+'\n')
    fo.write('#MSUB -x\n')
    fo.write('#MSUB -m work,scratch\n')

# Création du corps du script

    fo.write('cd ${BRIDGE_MSUB_PWD}\n')
    #fo.write('cat << END > app.conf\n')
    #fo.write(str(nbproc_roms)+' ./croco\n')
    #fo.write(str(nbproc_xios)+' ./xios_server.exe\n')
    #fo.write('END\n')
    if t!=1 or restart==1:
       #fo.write('cp ../t'+str(t-1)+'/moz_grd.nc .\n')  
       fo.write('ln -s '+WORKG+'/init_3000/moz_grd.*.nc .\n')    
       fo.write('cp ../t'+str(t-1)+'/moz_ini.*..nc .\n')    
       #fo.write('cp ../t'+str(t-1)+'/moz_bry.nc .\n')    
       fo.write('ln -s '+WORKG+'/init_3000/moz_bry.*.nc .\n')    
#    fo.write('ccc_mprun -f app.conf\n')
    fo.write('ccc_mprun ./croco\n')
    if t != nbchain:
        fo.write('cd ../t'+str(t+1)+'\n')
        fo.write('ccc_msub job_irene\n')
      
    fo.close()

###########################################################################

# Instruction pour lancer le batch
if restart==1:
    print 'Put the restart files in '+RPATH+'/t0'
print 'Change directory: cd '+RPATH+'/t1'
print 'Run commande : qsub job_irene'
