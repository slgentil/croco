#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD

source ../CONFIGURE_REG

#echo 'common scripts
ln -sf $dir_home/../CONFIGURE_GLOBAL .
ln -sf $dir_home/../CONFIGURE_REG .

[ -d TEST_CASES ] && rm TEST_CASES
ln -sf $dir_home/../TEST_CASES_CVTK TEST_CASES
ln -sf $dir_home/../TEST_CASES_CVTK/VHR/croco.in.VHR croco.in
ln -sf $dir_home/../TEST_CASES_CVTK/VHR/croco.in.1.VHR croco.in.1
ln -sf $dir_home/../TEST_CASES_CVTK/VHR/AGRIF_FixedGrids.in.REGIONAL.VHR AGRIF_FixedGrids.in
ln -sf $dir_datafile/VHR/*.nc* .

ln -sf $dir_home/../jobcomp_rvtk.bash .
ln -sf $dir_home/../extract_results_croco.bash .
ln -sf $dir_home/../comp_run_*.bash .

# specific scripts
ln -sf $dir_home/src_test/test_croco_reg.sh .
ln -sf $dir_home/src_test/rvtk_fast_qsub_REGIONAL.bash .
#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
