#!/bin/bash
#===================================
#PBS -q omp
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mpiprocs=4
#PBS -l mem=20gb
#PBS -j oe 
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================
source configure_file

par1='MPI'
# Compile
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log
/bin/mv croco croco_${par1}.exe

# Run
$MPI_LAUNCH -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log

