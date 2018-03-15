#!/bin/bash
#echo '============================================================='
echo 'Create the link between CVTK_DEBUG_FAST_src and CVTK_DEBUG/test dir'

ROOTDIR="/home7/datahome/crocodev/GIT/croco/CVTK/test_repro/CVTK_DEBUG_FAST_src"
source "$ROOTDIR/CONFIGURE_REG"

rm -Rf $dir_test
mkdir $dir_test
mkdir $dir_test/Junk
ln -sf $dir_home/Configure_Test_reg $dir_test/

# configure files
ln -sf $dir_home/../CONFIGURE_GLOBAL $dir_test/
ln -sf $dir_home/../CONFIGURE_REG $dir_test/

# common scripts and programms
ln -sf $dir_home/../Log_Details $dir_test/
ln -sf $dir_home/../gather_recap.bash $dir_test/
ln -sf $dir_home/../gitinfo.sh $dir_test/
ln -sf $dir_home/../git_process.bash $dir_test/
ln -sf $dir_home/../mk_CLEANALL.bash $dir_test/
ln -sf $dir_home/..//mk_CHECKALL.bash $dir_test/
ln -sf $dir_web/Log_Summary $dir_test/

# ana specific and programms
ln -sf $dir_home/src_test/mk_TestDIR_reg.bash $dir_test/
ln -sf $dir_home/src_test/mk_TESTALL_reg.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_reg Configure_Test
cd -
#
echo 'Well done: Finish linking'
