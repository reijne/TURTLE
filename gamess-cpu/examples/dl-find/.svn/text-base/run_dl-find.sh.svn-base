#!/bin/bash

echo "Running dl-find"

date
for job in `cat ./jobs.list | grep -v "#"`
do
    echo " *** Running job: $job  ***"
    ../../bin/gamess < ${job}.in > LOG/${job}.log
done
date

# Now validate & check for any errors
export GAMESS_VTAB=./REF/dl-find.vtab

error=0
for job in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $job LOG/${job}.log || error=1
done

# Cleanup
rm -f neb*  path*

exit $error
