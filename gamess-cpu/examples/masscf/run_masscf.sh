#!/bin/bash

silent=0
# Get command-line options
while getopts "sp" flag $@
do
  case $flag in
      s) silent=1;;
      p) continue;;
      # ? means an invalid option
      ?) echo "$0: Problem with option: -$flag" | tee $logfile; exit 1;;
      # : means an option without an argument - aparently not...
      #:) echo "$0: Option $OPTARG requires an argument!" | tee $logfile; exit 1;;
  esac
done

export RUN_SILENT=${RUN_SILENT:-0}
if [ $silent -eq 1 ]
then
  export RUN_SILENT=1
  export GAMESS_VALIDATE_QUIETLY=1
fi
if [ $RUN_SILENT -eq 1 ]
then
  echo "Running masscf"
fi

binary=../../bin/gamess

date
for job in `cat ./jobs.list | grep -v "#"`
do
    if [ $RUN_SILENT -eq 1 ]
    then
      echo " *** Running job: $job  ***"
    fi
    $binary < ${job}.in > LOG/${job}.log
done
date

# Now validate & check for any errors
export GAMESS_VTAB=./REF/masscf.vtab

error=0
for job in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $job LOG/${job}.log || error=1
done

# Cleanup
rm -f ed*

exit $error
