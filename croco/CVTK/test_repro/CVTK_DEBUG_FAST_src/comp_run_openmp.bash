#!/bin/bash
#===================================
#PBS -q omp
#PBS -l ncpus=4
#PBS -l walltime=02:00:00
#PBS -l mem=20gb
#PBS -j oe 
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================
source configure_file

par1='OPENMP'
# Compile
time ./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log
/bin/mv croco croco_${par1}.exe

# Run
export OMP_NUM_THREADS=$NBPROCS
./croco_${par1}.exe $CROCOIN > openmp_${NBPROCS}_${TEST_NAME}.log
