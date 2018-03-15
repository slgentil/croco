#!/bin/bash
#echo '============================================================='
echo 'Create the link between TESTROOT . and '$PWD

source ../CONFIGURE_VORT

#echo 'common scripts
ln -sf $dir_home/../CONFIGURE_GLOBAL .
ln -sf $dir_home/../CONFIGURE_VORT .

[ -d TEST_CASES ] && rm TEST_CASES
ln -sf $dir_home/../TEST_CASES_CVTK TEST_CASES
ln -sf $dir_home/../TEST_CASES_CVTK/croco.in.Vortex croco.in
ln -sf $dir_home/../TEST_CASES_CVTK/croco.in.Vortex.1 croco.in.1
ln -sf $dir_home/../TEST_CASES_CVTK/AGRIF_FixedGrids.in.VORTEX AGRIF_FixedGrids.in
ln -sf $dir_home/../TEST_CASES_CVTK/vortex_grd.nc* .
ln -sf $dir_home/../TEST_CASES_CVTK/vortex_ini.nc* .
ln -sf $dir_home/../TEST_CASES_CVTK/vortex_clm.nc .

ln -sf $dir_home/../jobcomp_rvtk.bash .
ln -sf $dir_home/../extract_results_croco.bash .
ln -sf $dir_home/../comp_run_*.bash .

# specific scripts
ln -sf $dir_home/src_test/test_croco_vort.sh .
ln -sf $dir_home/src_test/rvtk_fast_qsub_VORT.bash .
#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
