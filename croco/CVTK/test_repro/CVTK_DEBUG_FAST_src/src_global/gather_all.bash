#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -M gildas.cambon@ird.fr -m abe
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================

cd ${TESTROOTDIR}/REG
./gather_recap.bash REG
cd -

cd ${TESTROOTDIR}/VORT
./gather_recap.bash VORT
cd -

cd ${TESTROOTDIR}/KTEST
./gather_recap.bash KTEST
cd -
