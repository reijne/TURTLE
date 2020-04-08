#!/bin/sh 
# helper script for rg_exe.quadrics
# saves writing a shell to do the cd (which fails on
# our system with stale nfs handles.
scratch=$1
jobname=$2
exec=$3
home=$4
mkdir -p $scratch/$jobname
cd $scratch/$jobname
if test $RMS_RANK -eq 0
then
echo Start at `date`  directory `pwd`
fi
$exec  < $home/$jobname.in 
