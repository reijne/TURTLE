#!/bin/bash


# Now validate & check for any errors
export GAMESS_VTAB=./drf.vtab

error=0
for job in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $job LOG/${job}.log || error=1
done

exit $error
