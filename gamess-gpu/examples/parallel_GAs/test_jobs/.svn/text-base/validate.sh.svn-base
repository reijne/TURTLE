#!/bin/sh

# Simple script to validate the base outputs
# Default locations of the validate program, vtab file and
# list of jobs is set below - change these as necessary

# Work out where the base test directory is
if [ "${0:0:1}" = "/" ]
then
    # we've been called with an absolute path so just get the directory
    testdir=`dirname $0`
else
    # We've been called with a relative path so work out where the script lives
    # relvative to us, go there and print the directory name to get its full path
    path_to_script=`dirname $0`
    testdir=`(cd $path_to_script; pwd)`
fi

# Set the locate of the vtab file
export GAMESS_VTAB=$testdir/validate.vtab

# Set the validate command
val_cmd=`(cd ${testdir}/../../../utilities; pwd)`/validate
 
# Set the location of the list of jobs to validate
jobs_list=`(cd ${testdir}/../input_files/ ; pwd)`/jobs.list

# Function to state script usage
usage(){
echo "Usage is $0 [<directory>]"
echo "Where <directory> is an optional argument specifying the directory"
echo "containing the outputs to be validated."
echo "If no directory keyword is given the outputs are expected to reside"
echo " in the current directory."
exit 1
}

# See how we are being run
if [ $# -eq 1 ]
    then
    directory=`(cd $1; pwd)`
elif [ $# -eq 0 ]
    then
    directory=`pwd`
else
    usage
    exit 1
fi

# Check we can find the directory
if [ ! -d $directory ]; then
    echo "Cannot find directory: $directory!"
    usage
else
    echo "Validating the outputs in the directory: $directory"
fi

# Now cd to the directory and validate all the outputs
cd $directory
error=0
for job in `cat $jobs_list | grep -v "#"`
do
output=${job}.out
$val_cmd $job $output || error=1
done

# Exit with error code
exit $error
