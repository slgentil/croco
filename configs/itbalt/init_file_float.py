# -*- coding:Utf-8 -*-

"""
  python script to create the floats.in file which contains the initial positions of the floats
 
  [syntaxe] : python init_file_floats.py Fcoor xmin xmax ymin ymax nx ny
      arguments:
          Fcoor : type of coordinates (0=grid point, 1=meter or degree)
          xmin : min coordinate on x axis where to release floats (in kilometers)
          xmax : max coordinate on x axis
          ymin : min coordinate on y axis
          ymax : max coordinate on y axis
          nx : number of floats on x axis
          ny : number of floats on y axis
"""

import sys
import os
import shutil

# check number of arguments

if  len(sys.argv) < 9 :
    print('--------------------------------------------------------------------------------')
    print('[syntaxe] : python init_file_floats.py Fcoor xmin xmax ymin ymax nx ny z')
    print('--------------------------------------------------------------------------------')
    print('Fcoor : type of coordinates (0=grid point, 1=meter or degree)')
    print('xmin : min coordinate on x axis where to release floats (in meters)')
    print('xmax : max coordinate on x axis')
    print('ymin : min coordinate on y axis')
    print('ymax : max coordinate on y axis')
    print('nx : number of floats on x axis')
    print('ny : number of floats on y axis')
    print('z  : z depth (mètres,<0) or k level(>=1)')
    quit()

# Input parameters
# ----------------

# Type of coordinates : 0 = grid points, 1 = meters or degrees (depending on SPHERICAL flag)
Fcoor=int(sys.argv[1])
# Min and max coordinate of the distribution of floats 
xmin=float(sys.argv[2])
xmax=float(sys.argv[3])
ymin=float(sys.argv[4])
ymax=float(sys.argv[5])
nx=int(sys.argv[6])
ny=int(sys.argv[7])
Fz0=float(sys.argv[8])

# # Type of coordinates : 0 = grid points, 1 = meters or degrees (depending on SPHERICAL flag)
# Fcoor = 1

# # Min and max coordinate of the distribution of floats in meters
# xmin = 50000
# xmax = 975000
# ymin = 1040000
# ymax = 1840000
# nx = 30
# ny = 30

# ----------------

Fgrd = 0
dx = 0
if nx>1:
    dx = (xmax - xmin) / (nx - 1)
dy = 0
if ny>1:
    dy = (ymax - ymin) / (ny - 1)
Ft0 = 0.0
Ftype = 0
Fcount = 1
Fdt = 0.0
Fdx = 0.0
Fdy = 0.0
Fdz = 0.0

# Create file floats.in and write first 3 lines
fo = open('floats.in','w')
fo.write('1  Ftitle (a80)\n')
fo.write('CROCO - Initial Drifters Locations\n')
fo.write('2      Ft0,      Fx0,      Fy0,   Fz0, Fgrd,Fcoor,Ftype,Fcount,     Fdt,      Fdx,      Fdy,      Fdz\n')

for ifloat in range(nx):
    for jfloat in range(ny):

        Fx0 = xmin + (dx * (ifloat))
        Fy0 = ymin + (dy * (jfloat))

        fo.write('{0:10.1f}{1:10.1f}{2:10.1f}{3:7.1f}{4:6d}{5:6d}'\
            '{6:6d}{7:6d}{8:10.1f}{9:10.1f}{10:10.1f}{11:10.1f}\n'.\
            format(Ft0, Fx0, Fy0, Fz0, Fgrd, Fcoor, Ftype, Fcount, Fdt, Fdx, Fdy, Fdz))

fo.write('99 END of float input data')      
fo.close()