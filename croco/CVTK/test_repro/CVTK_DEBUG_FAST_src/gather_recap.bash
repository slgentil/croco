#!/bin/bash
#################################
# Gather all the log file to put them in a global one
# to have the echo on set -x
#set -x
# Use input argument as
# gather_recap KTEST,VORT, REG
#=================================

# $1: type of cas test
# $2: date
# usage : gather_recap.bash $1 $2
#---------------------------------
type_test=$1
today=$2
[ ! -n "$(echo "$2")" ] && today=`date +%Y%m%d`
#==
ligne=`grep -n revision gitinfos | cut -d: -f1`
ligne2=$((ligne + 1))
numrev=`head -$ligne2 gitinfos | tail -1 | tr -d '\n' | tr -d ' '`
#==
for testREGO in `ls Configure_Test` ; do  
    testREG=`echo $testREGO | cut -d/ -f2-`
    echo $testREG
    #echo "$testREG/Recap_${testREG}_${today}_*.git* Log_Summary"
    cp $testREG/Recap_${testREG}_${today}.git* Junk
done

cd Junk
for i in `ls -1 Recap_*${today}.git*` ; do 
    echo $i 
    cat $i >>  gather_recap_tmp
done
cd -

mv Junk/gather_recap_tmp "Log_Summary/$1_gather_recap_${today}_git${numrev}"
cp "Log_Summary/$1_gather_recap_${today}_git${numrev}" .
##rm -Rf Junk/*
