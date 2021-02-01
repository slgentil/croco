#!/usr/bin/python
# -*- coding:Utf-8 -*-

# """
#   script python pour supprimer les xios_client* et xios_server* (sauf le
#   premier de chaque)
#   dans le répertoire courant dans dans les deux niveaux inférieurs
# 
#   [syntaxe] : clean_croco.py
#
# """

import sys
import os
import shutil
import glob

def clean_files():
  """
  fonction pour définir et supprimer les fichiers
  """

  # find all croco clients (xios_client_*.out) 
  # and xios servers (xios_server_*.out)
  croco = glob.glob("xios_client_*.out")
  xios = glob.glob("xios_server_*.out")
  # sort lists
  croco.sort()
  xios.sort()
  # remove first element of list
  del croco[0]
  del xios[0]
  # get number of files in list
  nbcroco = len(croco)
  nbxios = len(xios)
  # convert list to string
  croco = ' '.join(croco)
  xios = ' '.join(xios)

  # # get exponent for cmd format (xios_client_01.out or xios_client_001.out)
  # logcroco = int(np.floor(np.log10(nbcroco)))+1
  # logxios = int(np.floor(np.log10(nbxios)))+1
  # cmd1 = "rm -f xios_client_{" + str(1).zfill(logcroco) + ".." + str(nbcroco).zfill(logcroco) + "}.out"
  # cmd2 = "rm -f xios_server_{" + str(nbcroco+2).zfill(logcroco) + ".." + str(nbcroco+nbxios+1).zfill(logcroco) + "}.out"
  # cmd3 = "rm -f xios_*.err"

  if nbcroco >=1:
    cmd1 = "rm " + croco
    # print(cmd1)
    os.system(cmd1)
  if nbxios >=1:
    cmd2 = "rm " + xios 
    # print(cmd2) 
    os.system(cmd2)
  cmd3 = "rm -f xios_*.err"
  # print(cmd3)
  os.system(cmd3)

  return

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
