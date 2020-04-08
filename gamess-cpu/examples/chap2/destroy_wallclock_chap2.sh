#!/bin/sh
#
# update for chap2
#
# This script updates all tolerance files where the test result
# did not match the reference result. As this script works 
# indiscrimately it should be used with great care!
#

# Set the paths to directories and programs
chap2dir=`pwd`

# Work out the toldiff directory
toldiffdir=$chap2dir/../../utilities/toldiff

destroy=./destroy_wallclock.py

# Set defaults for variables
parallel=0
error=0 # Flag if we encounter any errors (sets exit code of script)

for option
do
  case "$option" in
    -p) parallel=1;;
    *)  echo "unknown option $option"; exit 1;;
  esac
done

if ( test $parallel -eq 0 )
    # Serial jobs
    then
    # Set list of jobs
    job_list=`cat $chap2dir/jobs.serial | grep -v "#"`
    echo 
    echo "destroying the wall clock timings for the serial jobs"
    echo "====================================================="
    echo 
else
    # Set list of jobs
    job_list=`cat $chap2dir/jobs.parallel | grep -v "#"`
    echo 
    echo 
    echo "destroying the wall clock timings for the parallel jobs"
    echo "======================================================="
    echo 
fi

if ( ! test -e LOG-wallclock)
then
  mkdir LOG-wallclock
fi
for job in $job_list
do
  if ( test -e $chap2dir/LOG/${job}.log )
  then
      echo "destroying wall clock times in LOG/${job}.log"
      $destroy ${chap2dir}/LOG/${job}.log ${chap2dir}/LOG-wallclock/${job}.log
  else
      error=1
      echo "no output file LOG/${job}.log"
  fi
done
echo

# Exit with error code
exit $error
