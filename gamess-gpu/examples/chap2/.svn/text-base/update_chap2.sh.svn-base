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

toldiff=${toldiffdir}/toldiff
#toldiff="python $toldiffdir/toldiff.py"
#toldiff="/usr/local/packages/python/Python2.4/bin/python2.4 $toldiffdir/toldiff.py"

# Environment variables required for using tkdiff with toldiff
export TOLDIFF_EXE=`which diff`

# General toldiff environment variables
export TOLDIFF_OUTPUT=summary
export TOLDIFF_SUMMARY="OK:OK:Updating"
export TOLDIFF_EXIT="0:0:1"


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

echo "You are about to UPDATE ALL tolerance files in this directory."
echo "Are you sure about this? (y/n)"
read input
case "$input" in
  y) echo "Starting update";;
  n) echo "Aborting update"
     exit 0;;
  *) echo "Please answer y or n"
     exit 1;;
esac

echo

if ( test $parallel -eq 0 )
    # Serial jobs
    then
    # Set list of jobs
    job_list=`cat $chap2dir/jobs.serial | grep -v "#"`
    echo 
    echo "updating the tolerances for the serial jobs"
    echo "==========================================="
    echo 
else
    # Set list of jobs
    job_list=`cat $chap2dir/jobs.parallel | grep -v "#"`
    echo 
    echo 
    echo "updating the tolerances for the parallel jobs"
    echo "============================================="
    echo 
fi

for job in $job_list
  do
  if ( test -e $chap2dir/LOG/${job}.log )
      then
      result=`$toldiff $chap2dir/REF/${job}.log $chap2dir/LOG/${job}.log`
      status=$?
      printf "%30s %8s\n" ${job} ${result}
      if ( test \( ${status} -eq 1 \) )
	  then
          $toldiff --update ${chap2dir}/REF/${job}.log ${chap2dir}/LOG/${job}.log
      fi
  else
      error=1
      echo "no output file LOG/${job}.log"
  fi
done
echo

# Exit with error code
exit $error
