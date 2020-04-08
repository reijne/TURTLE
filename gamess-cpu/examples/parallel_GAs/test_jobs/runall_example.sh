#!/bin/bash
# This script was created under an LSF environment - 
# you will need to edit it accordingly for your queue system.
#BSUB -n 32
#BSUB -J valgamess
#BSUB -o valgamess.out.%J
#BSUB -W 04:00

# The base location of the GAMESS-UK directory
gukroot=/home/ngs0425/GAMESS-UK_dev

# The location of the binary to be validated
binary=$gukroot/bin/gamess-uk

# Where the input files can be found
inputdir=$gukroot/examples/parallel_GAs/input_files

# The name of the file holding the list of jobs
jobs_list=$inputdir/jobs.list

# The directory that the job runs in.
workdir=$LS_SUBCWD

# Create the list of nodes
nproc=0
hostlist=hostlist
> $hostlist
for h in $LSB_HOSTS
do
  echo $h >> $hostlist
  ((nproc++))
done

# Loop over all the jobs listed in the jobs_list file and submit them:
#for job in "/HF.crno4"
for job in `cat $jobs_list | grep -v "#"`
do
   mpirun -np $nproc -machinefile $hostlist $binary < $inputdir/${job}.in  > $workdir/${job}.out
done
