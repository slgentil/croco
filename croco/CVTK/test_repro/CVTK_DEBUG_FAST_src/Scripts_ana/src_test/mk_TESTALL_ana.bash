#!/bin/bash

./git_process.bash

for testconf in `ls -1 ./Configure_Test/ `;do
echo $testconf
rm -rf $testconf
./mk_TestDIR_ana.bash $testconf
done
#

exit


