#!/bin/bash
tutorials=`ls | grep "^[0-9]"`
execute="../../bin/hopr"
echo 'EXECUTE ALL TUTORIALS:'
for tutorial in $tutorials
do
  echo
  echo
  cd $tutorial
  inifiles=`ls *.ini`
  for inifile in $inifiles
  do
    echo '===> EXECUTING TUTORIAL: ' $tutorial ', WITH INIFILE:' $inifile
    outfile='calculation_'$tutorial'_'$inifile'.log'
    $execute $inifile > $outfile
    tail -n 4  $outfile
  done
  cd ..
done
