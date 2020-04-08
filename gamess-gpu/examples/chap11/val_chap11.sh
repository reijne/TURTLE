#!/bin/bash

export PATH=.:$PATH

# Now validate & check for any errors
export GAMESS_VTAB=./chap11.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
