#
# This file contains utility functions for buildbot and is executed by the buildbot_parallel and
# buildbot_serial scripts.
#
# Currently, it just has functions for submitting to different queueing systems.
#

################################################################################################
#
# Sun Grid Engine
#
################################################################################################
run_sge(){

# Submit job to an SGE queueing system, get the job id and poll till it 
# has finished.
#
# Check for a file called ERROR to see is we return with an 
# an error condition or not

local maxcount=240 # how many times to loop
local mysleep=1  # Polling interval in minutes

# Submit
qsub $1 > job.id || return 1

# Get job id
local jobid=`cat job.id | awk '{print $3}'`

local count=0
echo -n "Polling..."
while [ 1 ]
do

  sleep ${mysleep}m
  echo -n "$count ..."
  [ $count -ge $maxcount ] && return 1
  
  # See if we've finished
  qstat -j $jobid > /dev/null || return 0

  count=$((count+=1))
done
}

################################################################################################
#
# Loadleveller
#
################################################################################################

run_ll(){

# Submit job to a LoadLeveller queueing system, get the job id and poll till it 
# has finished.
#
# Check for a file called ERROR to see is we return with an 
# an error condition or not

local maxcount=360 # how many times to loop
local mysleep=60  # Polling interval in seconds

# Submit
llsubmit $1 > job.id || return 1


# Get job id
local jobid=`cat job.id  | grep "llsubmit: The job" | cut -d '"' -f2`
local count=0
echo -n "Polling..."
while [ 1 ]
do

  sleep $mysleep
  echo -n "$count ..."
  if [ $count -ge $maxcount ]
  then
	echo "Error - Job ran out of time!"
        llcancel $jobid
        touch ./ERROR
	return 1
  fi
  
  # See if we've finished
  llq -j $jobid | grep "llq: There is currently no job status to report." > /dev/null && return 0

  count=$((count+=1))
done
}


################################################################################################
#
# Sun Grid Engine
#
################################################################################################
run_lsf(){

# Submit job to an LSF queueing system, get the job id and poll till it 
# has finished.
#
# Check for a file called ERROR to see is we return with an 
# an error condition or not

local maxcount=240 # how many times to loop
local mysleep=1  # Polling interval in minutes

# Submit
bsub < $1 > job.id || return 1

# Get job id
local jobid=`cat job.id | cut -d "<"  -f2 | cut -d ">" -f1`

local count=0
echo -n "Polling..."
while [ 1 ]
do

  sleep ${mysleep}m
  echo -n "$count ..."
  if [ $count -ge $maxcount ]
  then
	echo "Error - Job ran out of time!"
        bkill $jobid
        touch ./ERROR
	return 1
  fi

  
  # See if we've finished
  # Output looks like:
  # JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
  # 105691  jmht93  DONE  n4_1       hapu33      4*lsfhost.l gamess-uk  Jul 28 13:53
  bjobs $jobid | awk '{getline; if ( $3 == "DONE" || $3 == "EXIT" ) exit 1}'  || return 0

  count=$((count+=1))
done
}
