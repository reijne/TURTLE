#!/bin/bash

echo "Running benchmark"

date
# just run the existing script
./run_bench
date

# Now validate & check for any errors
export GAMESS_VTAB=./bench.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
