#!/bin/bash
#

#
# Run every example we've got...
#

# Need to have . in the path in order to use the old scripts
export PATH=.:$PATH

norerun=0
silent=0
parallel=0
# Get command-line options
while getopts "nps" flag $@
do
  case $flag in
      n) norerun=1;;
      s) silent=1;;
      p) parallel=1;;
      # ? means an invalid option
      ?) echo "$0: Problem with option: -$OPTARG"; exit 1;;
      # : means an option without an argument - aparently not...
      #:) echo "$0: Option $OPTARG requires an argument!" | tee $logfile; exit 1;;
  esac
done

# Set if we rerun all jobs or just failing ones
if [ $norerun -eq 1 ]
then
   export NO_RERUN=1
fi

# Set if we want to minimize the output (particularly under the control of 
# buildbot). The environment variable RUN_SILENT is used in RUN.
# Also tell validate to shhh it unless an error is discovered.
if [ $silent -eq 1 ]
then
   export RUN_SILENT=1
   export GAMESS_VALIDATE_QUIETLY=1
fi

# Loop over all cases and run them
error=0
echo "Starting to validate all jobs at:"
date
for directory in `cat ./testdirs | grep -v "#"`
do
    pushd ./$directory > /dev/null
    if [ $parallel -eq 1 ]
	then
	./run_${directory}.sh -p
    else
	./run_${directory}.sh
    fi
    [ $? -ne 0 ] && error=1
    popd > /dev/null
done
echo "Finished testing all jobs at:"
date

# Let the user know if any jobs failed

# Create file called ERROR if any failed to let the automatic testing
# framework know what's happening.
[ $error -ne 0 ] && touch ./ERROR

exit $error
