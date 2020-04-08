#!/bin/bash

echo "Running mopac"

export PATH=.:$PATH

date
# just run the existing script
./run_mopac7
date

# Now validate & check for any errors
export GAMESS_VTAB=./mopac7.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
