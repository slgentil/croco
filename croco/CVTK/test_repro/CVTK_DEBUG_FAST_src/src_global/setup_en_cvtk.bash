#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -M gildas.cambon@ird.fr -m abe
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================

#-1 Set-up the CVTK_FAST environement for reg, ana and vort type
./Scripts_reg/create_link_master_reg.sh
./Scripts_ana/create_link_master_ana.sh
./Scripts_vort/create_link_master_vort.sh
