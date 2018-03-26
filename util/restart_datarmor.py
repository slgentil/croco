#!/usr/bin/python
# -*- coding:Utf-8 -*-

"""
  script python pour relancer une simulation de plusieurs chaînages sur datarmor
 
  [syntaxe] : restart_datarmor tmin tmax jobname
      arguments:
          tmin : chainage min à relancer
          tmax : chainage max à relancer
          jobname : nom générique des batchs
          ex : restart_datarmor 2 10 monjob
"""

import sys
import os

if  len(sys.argv) < 3 :
    print '----------------------------------------------------------------------'
    print '[syntaxe] : restart_qg.py tmin tmax jobname' 
    print '----------------------------------------------------------------------'
    print 'tmin : chainage min à relancer'
    print 'tmax : chainage max à relancer'
    print 'jobname : nom générique des batchs' 
    print '   se placer dans le répertoire qui contient t1,t2,...'
    print "   ex: restart_datamor.py 2 10 monjob" 
    quit()

tmin=int(sys.argv[1])
tmax=int(sys.argv[2])
jobname=sys.argv[3]

startdir=os.getcwd()
USER = os.getenv('USER')

for t in range(tmin,tmax+1):    
    tdir='t'+str(t)   
    os.chdir(startdir+'/'+tdir)
    
    if (t == tmin):
        commande='qsub -h job_datarmor'  
    else:
        f=os.popen('qselect -N '+jobname+str(t-1)+' -u '+USER)
        numjob=f.read()
        commande='qsub -W depend=afterany:'+str(numjob[:-1])+' job_datarmor'
    os.system(commande)
os.chdir(startdir+'/t'+str(tmin))
f=os.popen('qselect -N '+jobname+str(tmin)+' -u '+USER)
numjob=f.read()
commande='qrls '+numjob[:-1]     
print 'Run commande : '+commande