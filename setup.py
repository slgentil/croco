#!/usr/bin/env python

from distutils.core import setup

INSTALL_REQUIRES = ['xarray >= 0.10.6', ]

setup(name='croco_osi',
      description='Tools for croco pre/post processing, OSI-LOPS',
      url='https://github.com/slgentil/croco',
      packages=['croco_osi'])
