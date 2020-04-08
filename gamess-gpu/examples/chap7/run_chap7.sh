#!/bin/bash

echo "Running chap7"


date
# just run the existing script
./run_chap7
date

# Now validate & check for any errors
export GAMESS_VTAB=./chap7.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
