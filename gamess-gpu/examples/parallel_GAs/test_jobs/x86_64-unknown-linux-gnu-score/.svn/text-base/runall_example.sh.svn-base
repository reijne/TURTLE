#!/bin/bash

# Grid Engine stuff
#REM: nnodes + 1
#$ -pe score 4
#$ -masterq master.q@dl1.nw-grid.ac.uk
#$ -cwd
#$ -V
#$ -j y
#$ -o runtest.log

# How many processors per node
PPN=4

# Main GUK directory - change for your system
root=$SGE_O_WORKDIR/../../../../

# The location of the binary to be validated
binary=$root/bin/gamess-uk

# Where the input files can be found
inputdir=$root/examples/parallel_GAs/input_files

# The name of the file holding the list of jobs
jobs_list=$inputdir/jobs.list

# Loop over all the jobs listed in the jobs_list file and submit them:
for job in `cat $jobs_list | grep -v "#"`
do
  time scout -wait -F $HOME/.score/ndfile.$JOB_ID \
      -e /tmp/scrun.$JOB_ID -nodes=$((NSLOTS-1))x$PPN \
      $binary  < $inputdir/${job}.in > ${job}.out
done

# Now validate the jobs & create file depending on the result
../validate.sh && touch DONE || touch ERROR
