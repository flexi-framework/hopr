#!/bin/sh
for mf in "$@"
do
  if [ -e .tmp_$mf ]; then
    if perl -w compare_module_file.pl -compiler PLACEHOLDER_COMPILER_ID $mf .tmp_$mf; then
      touch -r .tmp_$mf $mf
    fi
    rm -f .tmp_$mf
  fi
done
