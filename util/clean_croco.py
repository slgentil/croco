#!/usr/bin/python
# -*- coding:Utf-8 -*-

"""
  script python pour supprimer les xios_client* et xios_server* (sauf 0)
  dans le répertoire courant dans dans les deux niveaux inférieurs

  example for 2048 procs for croco and 64 procs for xios
  rm -f xios_client_{0001..2047}.out
  rm -f xios_server_{0001..0063}.out
  rm -f xios_*.err

 
  [syntaxe] : clean_croco.py

"""

import sys
import os
import shutil
import numpy as np
import glob

def clean_files():
  """
  fonction pour définir et supprimer les fichiers
  """

  # find number of procs for croco (number of files xios_client_*.out) 
  # and xios (number of files xios_server_*.out)
  nbcroco = len(glob.glob("xios_client_*.out"))
  nbxios = len(glob.glob("xios_server_*.out"))
  nbcroco=max(int(nbcroco)-1,1)
  nbxios=max(int(nbxios)-1,1)

  # get exponent for cmd format (xios_client_01.out or xios_client_001.out)
  logcroco = int(np.floor(np.log10(nbcroco)))+1
  logxios = int(np.floor(np.log10(nbxios)))+1

  cmd1 = "rm -f xios_client_{"+str(1).zfill(logcroco)+".."+str(nbcroco).zfill(logcroco)+"}.out"
  cmd2 = "rm -f xios_server_{"+str(1).zfill(logcroco)+".."+str(nbxios).zfill(logcroco)+"}.out"
  cmd3 = "rm -f xios_*.err"

  if nbcroco>1 :
    print os.getcwd()
    print cmd1
    print cmd2
    print cmd3

  os.system(cmd1)
  os.system(cmd2)
  os.system(cmd3)

# Main

PWD=os.getcwd()

# Clean files in local directory
clean_files()   

# Clean files in sub-directory
for dir1 in os.listdir(PWD): 
  if os.path.isdir(PWD+'/'+dir1):
    os.chdir(PWD+'/'+dir1)
    clean_files()

    # Clean files in sub-sub-directory
    for dir2 in os.listdir(PWD+'/'+dir1):    
      if os.path.isdir(PWD+'/'+dir1+'/'+dir2):
        os.chdir(PWD+'/'+dir1+'/'+dir2)
        clean_files()
