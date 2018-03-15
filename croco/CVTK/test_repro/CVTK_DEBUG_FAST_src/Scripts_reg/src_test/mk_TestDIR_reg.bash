#!/bin/bash
#set -x
source CONFIGURE_REG

[ ! -d gitinfos ] && ./git_process.bash

# Directory creation
echo "======================================================"
[ ! -d $1 ] && mkdir $1
echo 'Create and setup dir. : ' $1

rm -Rf Configure_Test ; ln -sf  Configure_Test_reg Configure_Test

cp Configure_Test_reg/$1 $1/configure_file 
cp gitinfos $1
ln -sf $dir_home/src_test/create_link_for_Test_reg.sh $1
#===========================================================
cd $1 ;  
export mytest=$1
./create_link_for_Test_reg.sh
./test_croco_reg.sh
cd -

echo "======================================================"
echo "  "

