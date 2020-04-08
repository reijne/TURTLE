#!/bin/bash

# Simple script for working out what the compiler is breaking
# Compiler the code with no optimisation and then create a file
# with a list of objects to recompile in it called obj.list
# The script will then recompile the object with the flags specified
# in the variable compile and run the tests listed in "tests" and stop
# if one of the tests fails

filelist=`cat obj.list`
tests="c2028_e c2028_f"
compile="xlf -O3 -c -qnosave -qarch=pwr4    -qEXTNAME -qxlf77=leadzero"

#
# Function to run the tests in chap2
#
root=/hpcx/home/c01/c01/jmht/bbot/GAMESS-UK
chap2=$root/examples/chap2
export GAMESS_VTAB=$chap2/chap2.vtab
validate=$root/utilities/validate
runtest(){
    pushd $chap2
    for t in $tests
      do
      echo RUNNING TEST $t
      $chap2/$t > $chap2/LOG/$t.log
      $validate $t $chap2/LOG/$t.log
      if [ $? -ne 0 ]
	  then
	  echo "FAILED TEST: $t"
	  popd
	  return -1
      fi
    done
    popd
}

#runtest
#exit

#
# Loop over the objects
#
for file in $filelist
do
  echo "CHECKING FILE: $file"
  name=`basename $file .o`
  echo "$compile $name.f"
  $compile $name.f || { echo "ERROR COMPILLING FILE: $name"; exit 1; }
  make > make.log 2>&1 || { echo "ERROR IN MAKE"; exit 1; }
  runtest || { echo "FOUND ERROR WITH FILE: $name"; exit 1; }
done
