#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -l mem=20gb
#PBS -j oe 

cd $PBS_O_WORKDIR
#===================================
source configure_file
par1='SERIAL'

./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log 
mv croco croco_${par1}.exe

#Run (and produce check file)
time ./croco_${par1}.exe $CROCOIN > serial_${TEST_NAME}.log

