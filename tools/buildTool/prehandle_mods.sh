#!/bin/sh
for mf in "$@"
do
  if [ -e $mf ]; then
    cp -p $mf .tmp_$mf
  fi
done
