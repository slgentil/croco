#!/bin/bash
#
#======================================================================
# ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al) 
# and Rutgers University (Arango et al) are under MIT/X style license.
# ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
# 
# ROMS_AGRIF website : http://www.croco-ocean.org
#======================================================================
#
#---------------------------------------------------------------------
# Script to Run CVTK DEBUG procedure managing parallelization type 
# AND AGRIF nesting type (No nesting, Nesting 1-way, Nesting 2-ways) : 
# VORTEX and REGIONAL
#--------------------------------------------------------------------
#set -x

echo "=============================================="
echo "=> CONFIG "$mytest
echo " "
echo "=> MPIRUN COMMAND: "$MPIRUN
echo "Remove *.exe* *.log* "
[ ! -z "$(ls *.exe* 2>/dev/null)" ] && /bin/rm *.exe*
[ ! -z "$(ls *.log* 2>/dev/null)" ] &&/bin/rm *.log*
echo "Remove the CHECKFILE"
[ -f check_file ] && /bin/rm check_file
echo "Remove AGRIF_FixedGrids.in"
/bin/rm -f AGRIF_FixedGrids.in 
echo " " 

#=============================================================================================
# Sources
#
#sed -n -e '/SOURCE=/p' jobcomp_rvtk.bash > tmp1
#sed -n '$p' tmp1 > tmp2
#eval "SOURCE=`sed -n -e '/SOURCE=/ s/.*\= *//p' tmp2`"
#rm -f tmp1 tmp2
#echo 'Sources code: '$SOURCE
#===
source CONFIGURE_GLOBAL
#===
SOURCE_CVTK=${SOURCE_CROCO}/../CVTK/test_repro/CVTK_DEBUG_FAST_src
echo 'Sources CVTK tests: '$SOURCE_CVTK
#
# Get updated files
#
#/bin/cp ${SOURCE_CVTK}/Config_files/cppdefs_dev_cvtk.h cppdefs_dev_cvtk.h
#/bin/cp ${SOURCE_CVTK}/Config_files/cppdefs_cvtk.h cppdefs_bak1.h.SERIAL
#/bin/cp ${SOURCE_CVTK}/Config_files/param_cvtk.h param_bak0.h.SERIAL
#/bin/cp ${SOURCE_CVTK}/Config_files/cppdefs_cvtk.h cppdefs_bak1.h.OPENMP
#/bin/cp ${SOURCE_CVTK}/Config_files/param_cvtk.h param_bak0.h.OPENMP
#/bin/cp ${SOURCE_CVTK}/Config_files/cppdefs_cvtk.h cppdefs_bak1.h.MPI
#/bin/cp ${SOURCE_CVTK}/Config_files/param_cvtk.h param_bak0.h.MPI
#
/bin/cp ${SOURCE_CROCO}/cppdefs_dev.h cppdefs_dev_cvtk.h
sed 's/'undef\ \ \*RVTK_DEBUG'/'define\ RVTK_DEBUG'/' < cppdefs_dev_cvtk.h > cppdefs_dev_cvtk.h.tmp
mv cppdefs_dev_cvtk.h.tmp cppdefs_dev_cvtk.h

/bin/cp ${SOURCE_CROCO}/cppdefs.h cppdefs_bak1.h.SERIAL
/bin/cp ${SOURCE_CROCO}/param.h param_bak0.h.SERIAL

/bin/cp ${SOURCE_CROCO}/cppdefs.h cppdefs_bak1.h.OPENMP
/bin/cp ${SOURCE_CROCO}/param.h param_bak0.h.OPENMP

/bin/cp ${SOURCE_CROCO}/cppdefs.h cppdefs_bak1.h.MPI
/bin/cp ${SOURCE_CROCO}/param.h param_bak0.h.MPI
#==

# Source config file
#LIST_KEY0='PSOURCE PSOURCE_NCFILE FRC_BRY CLIMATOLOGY TIDES AGRIF AGRIF_2WAY BULK_FLUX MPI OPENMP'
#LIST_KEY_PHYS='REGIONAL FRC_BRY CLIMATOLOGY TIDES BULK_FLUX'
#LIST_KEY_PAR='OPENMP MPI' ; FLAG_MPI=1 ; FLAG_OPENMP=1
#LIST_KEY_NEST=''
#KEY_DEBUG='RVTK_DEBUG'
#LIST_WORDS='ETALON difference: ABNORMAL ERROR BUGBIN GRID#'
#CONFIG_NAME='BENGUELA_VHR'
#/bin/ln -sf AGRIF_FixedGrids.in.REGIONAL AGRIF_FixedGrids.in

#List of test cases with only one points in one direction
LIST_2DV_X='GRAV_ADJ IGW INNERSHELF INTERNAL SHOREFACE SWASH THACKER'
LIST_2DV_Y='OVERFLOW SHELFRONT'

source configure_file

Is2DV_X=0
Is2DV_Y=0
[ -n "$(echo $LIST_2DV_X |grep "${CONFIG_NAME}")" ] && Is2DV_X=1
[ -n "$(echo $LIST_2DV_Y |grep "${CONFIG_NAME}")" ] && Is2DV_Y=1

##############################################################################
#
#   FILL THE CPPDEFS.H 
# Title
echo TESTS OF $LIST_CONFIG
#
# 1- UNDEF ALL THE KEYS
#
echo 'undef '$LIST_KEY0
for par in SERIAL OPENMP MPI ; do 
    echo 'PARA=' $par
    for EXAMPLE in $LIST_KEY0
    do
	sed 's/'define\ \ \*$EXAMPLE'/'undef\ $EXAMPLE'/' < cppdefs_bak1.h.$par > cppdefs_bak2.h.$par
	/bin/mv cppdefs_bak2.h.$par cppdefs_bak1.h.$par
    done
    
    # 2- DEFINE THE TYPE OF DEBUG TEST 
    sed 's/'undef\ \ \*$KEY_DEBUG'/'define\ $KEY_DEBUG'/' < cppdefs_bak1.h.$par > cppdefs_bak2.h.$par
    /bin/mv cppdefs_bak2.h.$par cppdefs_bak1.h.$par
    
    
    # 3- DEFINE THE NAME OF THE CONFIG
    sed 's/'undef\ \*BENGUELA_LR'/'define\ $CONFIG_NAME'/' < cppdefs_bak1.h.$par > cppdefs_bak2.h.$par
    /bin/mv cppdefs_bak2.h.$par cppdefs_bak1.h.$par
    
    # 4- DEFINE THE VARIOUS CPPKEYS
    #=4.1
    for EXAMPLE in $LIST_KEY_PHYS
    do
	sed 's/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs_bak1.h.$par > cppdefs_bak2.h.$par
	/bin/mv cppdefs_bak2.h.$par cppdefs_bak1.h.$par
    done
    #==4.2
    for EXAMPLE in $LIST_KEY_NEST
    do
	echo $EXAMPLE
	sed 's/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs_bak1.h.$par > cppdefs_bak2.h.$par
	/bin/mv cppdefs_bak2.h.$par cppdefs_bak1.h.$par
    done
done

#=====================================================================================================
#=====================================================================================================
#####
echo 'FLAG OPENMP' ${FLAG_OPENMP}
echo 'FLAG MPI' ${FLAG_MPI}
####
#echo ' '
echo '==============================='
echo 'START TESTING ...             '
#echo '==============================='

# 3- 'MPI'
##############################################################################
# Serial runs
##############################################################################

#echo ''
par1='SERIAL'
echo "SERIAL NPP=1 TEST $mytest"
[ -f check_file ] && /bin/rm check_file
sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=1'/' < param_bak0.h.$par1 > param_bak1.h.$par1
sed 's/'NPP=1'/'NPP=1'/' < param_bak1.h.$par1 > param_bak2.h.$par1
/bin/mv param_bak2.h.$par1 param_bak1.h.$par1 ; rm param_bak0.h.$par1

rm -Rf Compile_$par1 ; mkdir Compile_$par1
cp param_bak1.h.$par1 Compile_$par1/param.h.OK
cp cppdefs_bak1.h.$par1 Compile_$par1/cppdefs.h.OK

#echo "qsub -h -N ${TEST_NAME}_SE comp_run_serial.bash"
qsub -h -N ${TEST_NAME}_SE comp_run_serial.bash
myjobid_serial="`qselect -N ${TEST_NAME}_SE -u $USER`"

# 4- 
##############################################################################
# Openmp runs
##############################################################################
if [ ${FLAG_OPENMP} = 1 ]; then 
    
    par1='OPENMP'
    if [ $Is2DV_Y == 1 ]; then
	echo "OPEN-MP 1X2 NPP=2 TEST $mytest"
	#export OMP_NUM_THREADS=4
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=2'/' < param_bak0.h.$par1 > param_bak1.h.$par1
	sed 's/'NPP=4'/'NPP=2'/' < param_bak1.h.$par1 > param_bak2.h.$par1
	
    elif [ $Is2DV_X == 1 ]; then
	echo "OPEN-MP 2x1 NPP=2 TEST $mytest"
	#export OMP_NUM_THREADS=4
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=2,\ NSUB_E=1'/' < param_bak0.h.$par1 > param_bak1.h.$par1
	sed 's/'NPP=4'/'NPP=2'/' < param_bak1.h.$par1 > param_bak2.h.$par1
	
    else
	echo "OPEN-MP 2X2 NPP=4 TEST $mytest"
	#export OMP_NUM_THREADS=4
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=2,\ NSUB_E=2'/' < param_bak0.h.$par1 > param_bak1.h.$par1
	sed 's/'NPP=4'/'NPP=4'/' < param_bak1.h.$par1 > param_bak2.h.$par1
	
    fi	 
    
    /bin/mv param_bak2.h.$par1 param_bak1.h.$par1 ; rm param_bak0.h.$par1
    sed 's/'undef\ \ \*${par1}'/'define\ ${par1}'/' < cppdefs_bak1.h.$par1 > cppdefs_bak2.h.$par1
    /bin/mv cppdefs_bak2.h.$par1 cppdefs_bak1.h.$par1
    
    rm -Rf Compile_$par1 ; mkdir Compile_$par1
    cp param_bak1.h.$par1 Compile_${par1}/param.h.OK
    cp cppdefs_bak1.h.$par1 Compile_${par1}/cppdefs.h.OK
    
    #   echo "qsub -N ${TEST_NAME}_OM -W depend=afterok:$myjobid_serial comp_run_openmp.bash" 
    qsub -N ${TEST_NAME}_OM -W depend=afterok:$myjobid_serial comp_run_openmp.bash
    
    #   echo ""myjobid_openmp=`qselect -N ${TEST_NAME}_OM -u $USER`""
    myjobid_openmp="`qselect -N ${TEST_NAME}_OM -u $USER`"
fi 

###############################################################################

# 4- RVTK_DEBUG_REG_DEV
##############################################################################
# Mpi runs
##############################################################################
if [ ${FLAG_MPI} = 1 ]; then 
    
    par1='MPI'
    if [ $Is2DV_Y == 1 ]; then
	echo "MPI 1X2 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=1,\ NP_ETA=2'/' < param_bak0.h.$par1 > param_bak1.h.$par
	
    elif [ $Is2DV_X == 1 ]; then
	echo "MPI 2X1 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=2,\ NP_ETA=1'/' < param_bak0.h.$par1 > param_bak1.h.$par
    
    else	
	echo "MPI 2X2 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=2,\ NP_ETA=2'/' < param_bak0.h.$par1 > param_bak1.h.$par1
    fi
    
    rm param_bak0.h.$par1
    sed 's/'undef\ \ \*$par1'/'define\ $par1'/' < cppdefs_bak1.h.$par1 > cppdefs_bak2.h.$par1
    /bin/mv cppdefs_bak2.h.$par1 cppdefs_bak1.h.$par1
    rm -Rf Compile_$par1 ; mkdir Compile_$par1
    cp param_bak1.h.$par1 Compile_${par1}/param.h.OK
    cp cppdefs_bak1.h.$par1 Compile_${par1}/cppdefs.h.OK
    
    #  echo "qsub -N mpi_${TEST_NAME}_MP -W depend=afterok:$myjobid_serial comp_run_mpi.bash"
    qsub -N ${TEST_NAME}_MP -W depend=afterok:$myjobid_serial comp_run_mpi.bash
    
    #  echo "myjobid_mpi="`qselect -N ${TEST_NAME}_MP -u $USER`""
    myjobid_mpi="`qselect -N ${TEST_NAME}_MP -u $USER`"
fi

#########################################################################################################
# 5- Extract results
##############################################################################
#  runs
#echo ' '
echo "EXTRACTION $mytest"
if [[ ${FLAG_MPI} = 1 &&  ${FLAG_OPENMP} = 1 ]]; then 
    #echo "qsub -N ${TEST_NAME}_EX -W depend=afterok:${myjobid_mpi}:${myjobid_openmp} extract_results_croco.bash"
    qsub -N ${TEST_NAME}_EX -W depend=afterany:${myjobid_mpi}:${myjobid_openmp} extract_results_croco.bash
    
elif [ ${FLAG_OPENMP} = 1 ]; then 
    #echo "qsub -N ${TEST_NAME}_EX -W depend=afterok:${myjobid_openmp} extract_results_croco.bash"
    qsub -N ${TEST_NAME}_EX -W depend=afterany:${myjobid_openmp} extract_results_croco.bash
    
elif [ ${FLAG_MPI} = 1 ]; then 
    #echo "qsub -N ${TEST_NAME}_EX -W depend=afterok:${myjobid_mpi} extract_results_croco.bash"
    qsub -N ${TEST_NAME}_EX -W depend=afterany:`qselect -N mpi_${TEST_NAME} -u $USER` extract_results_croco.bash
fi
#########################################################################################################


qrls `qselect -N ${TEST_NAME}_SE`
