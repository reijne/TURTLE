#!/bin/sh
#
# Diff or compare the tolerances for all files.
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
export TOLDIFF_SUMMARY="N/A:N/A:N/A"
export TOLDIFF_EXIT="0:0:0"

unset TKDIFFRC


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
    echo "diffing the tolerances for the serial jobs"
    echo "=========================================="
    echo 
else
    # Set list of jobs
    job_list=`cat $chap2dir/jobs.parallel | grep -v "#"`
    echo 
    echo 
    echo "diffing the tolerances for the parallel jobs"
    echo "============================================"
    echo 
fi

for job in $job_list
  do
  if ( test -e $chap2dir/REF/${job}.tol )
      then
      ./diff_tolerances.sh "${job}"
  else
      error=1
      echo "no tolerance file REF/${job}.tol"
  fi
done
echo

# Exit with error code
exit $error
