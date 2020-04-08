#!/bin/bash

echo "Running drf"

date
# Loop over jobs and run them
for job in `cat ./jobs.list | grep -v "#"`
do
  ./RUN ./$job
  # Need to keep files for cd013 jobs
  if [ ! ${job:0:6} = "cd013_" ]
      then
      /bin/rm -f da* neq* mfg* gamess*
    fi
done
date
/bin/rm -f adapt diag sort tran mfg*
/bin/rm -f rpa* tda* tm* spectrum residues roots ft* table.tex fort.*
/bin/rm -f con4

# Now validate & check for any errors
export GAMESS_VTAB=./drf.vtab

error=0
for job in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $job LOG/${job}.log || error=1
done

exit $error