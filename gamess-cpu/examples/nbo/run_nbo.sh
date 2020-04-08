#!/bin/bash

echo "Running nbo"

date
# just run the existing script
./run_nbo
date

# Now validate & check for any errors
export GAMESS_VTAB=./nbo.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
