#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l place=scatter:excl
#PBS -o drf.pbs.log
#PBS -e drf.pbs.err
#PBS -N DRF
#PBS -l walltime=1:00:00
#PBS -q workq
#PBS -P PR39

# latest intel compilers, mkl and intel-mpi

module purge
module load compiler/intel-15.0

ulimit -s unlimited
ulimit -c 0

export OMP_NUM_THREADS=1

cat $PBS_NODEFILE

code=${HOME}/MFG/GAMESS-UK/bin/gamess
MYPATH=$HOME/MFG/GAMESS-UK/examples/drf
MYDATA=$HOME/MFG/GAMESS-UK/examples/drf
WDPATH=/scratch/$USER/GAMESS-UK-drf

#
# Test jobs from DFT Chapter of GAMESS Manual
#
date

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

# run jobs in examples directory
#
cd $MYPATH
# Set which jobs we're running
joblist=$MYDATA/jobs.list

env

echo "Running drf examples."
date
# Loop over all cases and run them
for job in `cat $joblist | grep -v "#"`
do
    $MYPATH/RUN $job
if [ ! ${job:0:6} = "cd013_" ]
      then
      /bin/rm -f da* neq* mfg* gamess*
    fi
done
date
/bin/rm -f adapt diag sort tran mfg*
/bin/rm -f rpa* tda* tm* spectrum residues roots ft* table.tex fort.*
/bin/rm -f con4

$MYPATH/val_drf.sh -o
#
