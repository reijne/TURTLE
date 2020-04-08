#!/bin/sh
# Jens' script to scour the cvs archives to find
# when a particular change was made to a file 14/6/05

usage(){
echo
echo " This script will search through all the cvs revisions of a file"
echo " to find when a particular change occured."
echo
echo " Usage is as follows:"
echo
echo " $0 \"<string_to_search_for>\" <file>"
echo
exit 1
}

# Check we've got the correct number of command-line arguments
[ $# -ne 2 ] && usage

# Set the variables
change=$1
infile=$2
outfile=$infile.cvsdiff

# Check we can find the input file
if [ ! -f $infile ]; then
echo
echo " ### Error! Cannot find find file $infile ###"
usage
fi

# Get the revision numnbers
revisions=`cvs log $infile | grep "revision\>" | cut -d" " -f2`

# Dirty hack to simulate an array
i=0
for r in $revisions
do
  i=`expr $i + 1`
  eval rev$i="$r"
done


# find out how many elements we have
narray=$i

# The following loop starts with the most recent revision, and diffs it with 
# the preceeding revision. The output of the diff is sent into $outfile
# and also piped to grep, if grep find the string, the cvs log for that change
# is appended to $outfile and the loop stopped. Otherwise the loop continues

# ( Is a bit dirty, but \$rev$i in the below translates to (e.g. $rev35)
# ie the 35th element of the array holding the revision numbers )

i=1
found=0
while [ $found -ne 1 ]
  do 
  eval r1=\$rev$i
  i=`expr $i + 1`
  # if we are at the end of the revisions crash out...
  [ $i -gt $narray ] && echo "Could not find change" && exit
  eval r2=\$rev$i

  echo "cvs diffing $r1 and $r2 ..."
  cvs diff -r $r1 -r $r2 $infile | tee $outfile | grep "$change" > /dev/null
  if [ $? -eq 0 ]; then
      # Found the change - say so and then append the log info
      # for that change to the output
      echo
      echo "Found change between revisions $r1 and $r2"
      cvs log -r$r1 $infile >> $outfile
      echo
      echo "Please see the file $outfile for the results"
      found=1
  fi
done
