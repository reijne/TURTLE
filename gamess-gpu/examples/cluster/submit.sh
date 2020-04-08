#!/bin/bash
# Parallel job submission script:
# Use qsub -pe score (nodes + 1) ./submit.sh
# where nodes is the number of nodes to be used
# to submit the job. 
#$ -S /bin/bash
# define the parallel queue for ccp1 (where job is run)
#$ -masterq ccp1.q
#$ -j y
# define the parallel runtime environment and number of nodes
# Note: number of nodes is one more than needed as one copy resides
# on the master node
#$ -pe score 5
# Use location that job was submitted as working directory
#$ -cwd
# export all environment variables to slave jobs
#$ -V
#$ -l h_rt=20:00:00
# If you stick to under 6 hours you will be able to access all queues
#$ -o par_6b.log
# tell sge where to put output

# run the job:
# scout is the SCORE remote shell runtime environment
# -wait: wait until all remote hosts are locked via MessageBoard (msgbserv)
# -F specify file in which the host names are listed line by line
# Given by SGEEEE
# -e create scout enviroment (final command of scout)
#  /tmp/scrun.$JOB_ID is a symbolic link to 
# /opt/score/bin/bin.i386-redhat7-linux2_4/scrun.exe
# All options which follow are those of scrun:
# -nodes NSLOTS will be one more than number of nodes as one
# resides on the master node; There are also 2 CPUs per node.
# Final command is the executable to be run

cp par_6b.in datain

scout -wait -F $HOME/.score/ndfile.$JOB_ID -e /tmp/scrun.$JOB_ID \
 -nodes=$((NSLOTS-1))x2 ../../bin/gamess-uk > par_6b.out.8.4X2
