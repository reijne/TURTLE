#!/bin/bash
#

#
# Test jobs from Chapter 2 of GAMESS Manual
#

# Set default values
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
      ?) echo "$0: Problem with option: -$OPTARG" | tee $logfile; exit 1;;
      # : means an option without an argument - aparently not...
      #:) echo "$0: Option $OPTARG requires an argument!" | tee $logfile; exit 1;;
  esac
done

# Set if we rerun all jobs or just failing ones
if [ $norerun -eq 1 ]
then
   export NO_RERUN=1
fi

# Set if we want to silence the running of the test jobs (RUN_SILENT used
# in RUN)
# Also tell validate to shhh it unless an error is discovered.
if [ $silent -eq 1 ]
then
   export RUN_SILENT=1
   export GAMESS_VALIDATE_QUIETLY=1
fi

# Set which jobs we're running
joblist=./jobs.list

# Create ecplib file
if [ ! -f ../../libs/ecplib ]
then
   pushd ../../libs;./make_lib > /dev/null ;popd > /dev/null
fi

# Create TABLE file for mrdci jobs (not if running in parallel under buildbot)
if [ $parallel -ne 1 ]
then
   if [ ! -f ../../libs/TABLE ]
   then
      pushd ../../libs;./make_table > /dev/null ;popd > /dev/null
   fi
fi

echo "Running chap2 examples."
date
# Loop over all cases and run them
for job in `cat $joblist | grep -v "#"`
do
    ./RUN $job
done
date

/bin/rm -rf adapt diag sort tran mfg*
/bin/rm -rf rpa* tda* tm* spectrum residues roots ft* table.tex fort.*
/bin/rm -rf table options.dft


#./val_chap2.sh -t
./val_chap2.sh -o
