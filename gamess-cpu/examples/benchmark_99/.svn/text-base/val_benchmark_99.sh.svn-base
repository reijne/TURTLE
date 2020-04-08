#!/bin/sh
#
# validation for benchmark_99
#
# By default only a brief report is generated.
# If the -g flag is presented tktoldiff is started if the test
# result did not match the reference data.
#

# Set the paths to directories and programs
benchmark_99dir=`pwd`

# Work out the toldiff directory
toldiffdir=$benchmark_99dir/../../utilities/toldiff

toldiff=${toldiffdir}/toldiff
#toldiff="python $toldiffdir/toldiff.py"
#toldiff="/usr/local/packages/python/Python2.4/bin/python2.4 $toldiffdir/toldiff.py"

# Need to set the path to the toldiff file for when we use tktoldiff
# as this executes a shell and expects to find toldiff in it's path
export PATH=$toldiffdir:$PATH


# Environment variables required for using tkdiff with toldiff
export TOLDIFF_EXE=`which diff`
export TKDIFFRC=${toldiffdir}/tktoldiffrc

# General toldiff environment variables
export TOLDIFF_OUTPUT=summary
export TOLDIFF_SUMMARY="OK:OK:Failed"
export TOLDIFF_EXIT="0:0:1"

# For when running the old vtab-based validation
export GAMESS_VTAB=./bench.vtab

# Set list of jobs
job_list=`cat $benchmark_99dir/jobs.list | grep -v "#"`

# Set defaults for variables
graphical=0
text=0
oldval=0 # Use old vtab-based validation
error=0 # Flag if we encounter any errors (sets exit code of script)

# Parse the command-line
for option
do
  case "$option" in
    -g) graphical=1;;
    -t) text=1;;
    -o) oldval=1;;
    *)  echo "unknown option $option"; exit 1;;
  esac
done


# Print header
echo 
echo "The results of the benchmark_99 examples"
echo "========================================"
echo 

for job in $job_list
  do
  if ( test -e $benchmark_99dir/LOG/${job}.log )
      then
      if ( test $oldval -eq 1 )
	  then
          #
	  # Using old vtab stuff
          # 
	  ../../utilities/validate $job $benchmark_99dir/LOG/${job}.log || error=1
      else
          #
	  # New toldiff approach
          #
	  result=`$toldiff $benchmark_99dir/REF/${job}.log $benchmark_99dir/LOG/${job}.log`
	  status=$?
	  printf "%30s %8s\n" ${job} ${result}
	  if ( test \( ${status} -eq 1 \) )
	      then
	      error=1
	      if ( test ${graphical} -eq 1 )
		  then
		  tkdiff ${benchmark_99dir}/REF/${job}.log ${benchmark_99dir}/LOG/${job}.log
	      elif ( test  ${text} -eq 1 )
		  then
		  $toldiff --output full ${benchmark_99dir}/REF/${job}.log ${benchmark_99dir}/LOG/${job}.log
	      fi
	  fi
      fi
  else
      error=1
      echo "no output file LOG/${job}.log"
  fi
done
echo

# Exit with error code
exit $error
