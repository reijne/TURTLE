#!/bin/bash

echo "Running benchmark_99"

date
# Loop over jobs and run them
for job in `cat ./jobs.list | grep -v "#"`
do
  ./RUN ./$job
done
date
/bin/rm -rf ft30 ft57 ft60 ft61 ft62 ft63 ft64 ft65 ft66 ft69 ft78 ft90 options.dft

# Now validate & check for any errors
export GAMESS_VTAB=./bench.vtab

error=0
for job in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $job LOG/${job}.log || error=1
done

exit $error
