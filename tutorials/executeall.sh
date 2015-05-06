#!/bin/bash
tutorials=`ls |grep [0-9]`
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
    outfile='calculation_'$inifile'.log'
    ../../bin/hopr $inifile > $outfile
    tail -n 4  $outfile
  done
  cd ..
done
