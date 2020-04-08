#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l place=scatter:excl
#PBS -o dft.pbs.log
#PBS -e dft.pbs.err
#PBS -N DFT
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
MYPATH=$HOME/MFG/GAMESS-UK/examples/dft
MYDATA=$HOME/MFG/GAMESS-UK/examples/dft
WDPATH=/scratch/$USER/GAMESS-UK-dft

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

# Create ecplib file
if [ ! -f ../../libs/ecplib ]
then
   pushd ../../libs;./make_lib > /dev/null ;popd > /dev/null
fi

env

echo "Running dft examples."
date
# Loop over all cases and run them
for job in `cat $joblist | grep -v "#"`
do
    $MYPATH/RUN $job
done
date

$MYPATH/val_dft.sh -o
#
/bin/rm -rf mfged2 mfged3 ftn058 options.dft sort
