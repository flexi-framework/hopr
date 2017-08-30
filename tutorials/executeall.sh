#!/bin/bash

if [ -z "$1" ]; then
  HOPREXEC=$(readlink -f ../build/bin/hopr)
else
  HOPREXEC=$(readlink -f "$1")
fi
if [ ! -f "$HOPREXEC" ]; then
  echo "Executable $HOPREXEC does not exist!"
  exit 1
fi

success=0
tutorials=`ls | grep "^[0-9]"`
echo 'EXECUTE ALL TUTORIALS:'
for tutorial in $tutorials
do
  echo
  cd $tutorial
  inifiles=`ls *.ini`
  for inifile in $inifiles
  do
    echo '===> EXECUTING TUTORIAL: ' $tutorial ', WITH INIFILE:' $inifile
    outfile='calculation_'$tutorial'_'$inifile'.log'
    $HOPREXEC $inifile > $outfile
    tmp=$?
    if [ $tmp != 0 ]; then
      success=$tmp
    fi
    tail -n 10  $outfile
  done
  cd ..
done
if [ $success != 0 ]; then
  echo "Errors have occurred for one or more example."
fi
exit $success
