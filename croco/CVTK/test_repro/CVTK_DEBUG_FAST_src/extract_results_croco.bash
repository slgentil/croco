#!/bin/bash
#===================================
#PBS -q sequentiel
#PBS -l walltime=02:00:00
#PBS -j oe 
#PBS -M gildas.cambon@ird.fr -m abe
cd $PBS_O_WORKDIR
echo $PBS_O_LOGNAME
#===================================
#set -x

#==
source configure_file
#==

numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
today=`date +%Y%m%d`
mv Recap_${TEST_NAME}.git${numrev} Recap_${TEST_NAME}_${today}.git${numrev}

#============================================
filein_openmp=openmp_${NBPROCS}_${TEST_NAME}.log
filein_mpi=mpi_${NBPROCS}_${TEST_NAME}.log

fileout_openmp=openmp_checkbugbin_${TEST_NAME}.txt
fileout_mpi=mpi_checkbugbin_${TEST_NAME}.txt

rm -Rf $fileout_openmp $fileout_mpi

touch $fileout_openmp
touch $fileout_mpi

GREP_CMD='grep -m 1'

#=============================================
# OPENMP
echo "===================" > $fileout_openmp
echo 'OPENMP CHECK (BUGBIN detection)' >> $fileout_openmp 
${GREP_CMD} BUGBIN $filein_openmp >> $fileout_openmp
res_omp=`${GREP_CMD} BUGBIN $filein_openmp`
echo 'res_omp='$res_omp >> $fileout_openmp
if [ -z $res_omp ] ; then 
echo 'check [passed]'  >> $fileout_openmp
else
echo 'check [failed]'  >> $fileout_openmp
fi

#MPI
echo "===================" > $fileout_mpi
echo 'MPI CHECK (BUGBIN detection)' >> $fileout_mpi
${GREP_CMD} BUGBIN $filein_mpi >> $fileout_mpi
res_mpi=`${GREP_CMD} BUGBIN $filein_mpi`
echo 'res_mpi='$res_mpi >> $fileout_mpi
if [ -z $res_mpi ] ; then 
echo 'check mpi [passed]'  >> $fileout_mpi
else
echo 'check mpi [failed]'  >> $fileout_mpi
fi

# MERGE
cat $fileout_mpi >>  $fileout_openmp
mv $fileout_openmp checkbugbin_${TEST_NAME}.txt

#============================================
file=Recap_${TEST_NAME}_${today}.git${numrev}
#echo '  ' >> $file
echo ' '
echo '>> oooooooooooooo <<' >> $file
echo 'REVISION GIT :' $numrev >> $file
echo 'DATE         :' $today >> $file
cat checkbugbin_${TEST_NAME}.txt >> $file 
echo '---------------------------' >> $file
echo ' ' >> $file
cp $file ../Log_Details
#echo '@@@@@@@@@@@@@@@@@@@@ ' >> $file

# Cleaning
rm -f $fileout_openmp $fileout_mpi
