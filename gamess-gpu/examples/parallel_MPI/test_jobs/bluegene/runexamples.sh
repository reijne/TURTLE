#!/bin/bash
# LoadLeveler is configured in such a way that we have no control
# over the mpirun command. So we have arrange everything beforehand.

binary=../../../../bin/gamess-uk
inputdir=../../input_files
jobs_list=`cat $inputdir/jobs.list | grep -v "#"`
me=`whoami`
mycwd=`pwd`
nproc=$1

for job in $jobs_list
do
  if (! test -e $mycwd/${job}.output.${nproc})
  then
    #
    # First make sure that this is going to be our only
    # job in the queue.
    #
    myqueue=`llq | grep $me`
    while ((${#myqueue}!=0)) 
    do
       sleep 180
       myqueue=`llq | grep $me`
    done
    #
    # OK, no jobs of mine in the queue so I can safely submit one
    # use: mpirun -h for more options (apparently queue limit is 12 hours)
    #
    cp $inputdir/${job}.in datain
    cat <<EOF > job.ll
#@ job_name = ${job}
#@ arguments = -np ${nproc} -exe ${binary} -cwd ${mycwd} -timeout 21600
#@ input = /dev/null
#@ output = $mycwd/${job}.output.${nproc}
#@ error = $mycwd/${job}.error.${nproc}
#@ notification = complete
#@ queue
EOF
   llsubmit job.ll
   sleep 180
  fi
done

/bin/rm -rf ed3-cyclo ed3-tot12 ed3-morphine ed3-valino ftn058
/bin/rm -rf options.dft datain

