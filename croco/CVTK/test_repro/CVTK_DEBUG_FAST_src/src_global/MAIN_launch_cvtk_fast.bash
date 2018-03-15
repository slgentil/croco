#!/bin/bash

#=========================================================
# This scripts launch the CVTK FAST procedure
# on DATARMOR for the 3 type of tests :
# - VORT  (as VORTEX)
# - KTEST (as TEST_CASES)
# - REG   (as REGIONAL)
#
# Usage : ./Lauch_CVTK_FAST.bash > CVTK_FAST.out 2>&1 
# It will launch
# 1) setup_env_cvtk.bash
# 2) exec_cvtk.bash  (if setup_env_cvtk.bash finish OK)
# 3) gather_all.bash (if exec_cvtk.bash finish OK)
# The scripts are using qsub on serial queue
#
#  Gildas Cambon, LOPS (gildas.cambon@ird.fr) February 2018
#=========================================================

export TESTROOTDIR /home7/datawork/crocodev/CVTK_DEBUG_REG_DEV/TESTROOT

#-1 Set-up the CVTK_FAST environement for reg, ana and vort type
qsub -N setup_env_cvtk setup_env_cvtk.bash

#-2 Lauch the various CVTK_FAST for reg, ana and vort
qsub -N exec_cvtk -W depend=afterok:setup_env_cvtk exec_cvtk.bash

#-3 Get the recap for each CVTK_FAST type 
qsub -N gather_all -W depend=afterok:exec_cvtk gather_all.bash

#-End 
qrls `qselect -N setup_env_cvtk`
