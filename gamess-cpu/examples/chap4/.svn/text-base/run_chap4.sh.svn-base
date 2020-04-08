#!/bin/bash

echo "Running chap4"

date
# just run the existing script
./run_chap4
date

# Now validate & check for any errors
export GAMESS_VTAB=./chap4.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
