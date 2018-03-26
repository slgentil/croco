#!/usr/bin/python
# -*- coding:Utf-8 -*-

"""
  script python pour tuer une simulation de plusieurs chaînages sur datarmor
 
  [syntaxe] : kill_datarmor.py jobname
      arguments:
          jobname : nom générique des batchs
"""

import sys
import os

if  len(sys.argv) < 1 :
    print '----------------------------------------------------------------------'
    print '[syntaxe] : kill_datarmor.py jobname' 
    print '----------------------------------------------------------------------'
    print 'jobname : nom générique des batchs' 
    print "   ex: kill_datarmor.py monjob" 
    quit()

jobname=sys.argv[1]

startdir=os.getcwd()
USER = os.getenv('USER')
substring='.datarmor'

f=os.popen('qstat | grep '+jobname+'| grep '+USER)
lines = f.readlines()
for index, line in enumerate(lines):
    if substring in line:
		ind=line.find(substring)
		commande = 'qdel -W force '+line[:ind]
		os.system(commande)
