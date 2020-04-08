#!/bin/bash
#
# bash script to check when a particular change was introduced into the code
#
# Given a starting version, the script checks out the code as of a particular
# revision (specified by the variable thisrev)
# subsequent to a given starting date (which is set by the variable nextdate)
#
# The script then configures the code with the -s option (so a recent version
# of the configure script must be used and the file configure.dat must contain
# the build and compiler option) and builds the code.
#
# The script then runs the test specified in the "testcode" variable
# (running the function run_chap2_eg in this case, which runs a
# chap2 test as specified by the local variable "c2test" in that function).
# If the test passes, the source
# code is update to the day of the next change, the code rebuilt and the
# test re-run
#
# The root variable must be set to the root of the GAMESS-UK directory
#
# If a file called date.list exists in the root directory, it is expected
# to contain a list of dates on which the code should be tested. If the
# file doesn't exist one is created that contains an entry for each date
# the the repository was altered.

start=5262 # A known working revision
end=5293   # 0 implies current repository head -  otherwise put a failing version here

testcode=run_chap2_eg # The test to run
makefailok=0          # Are makes allowed to fail?
root=`pwd`            # The root of the GAMESS-UK directory

# Create a file to log everything
log=$root/check.log; > $log


run_chap2_eg(){
# See if a particular chap2 example is validating or not

# Specify the test to run
local c2test=c2001_a

# Need to put . in the path
export PATH=.:$PATH

c2dir=$root/examples/chap2
export GAMESS_VTAB=$c2dir/chap2.vtab
cd $c2dir
./RUN $c2test
[ $? -ne 0 ] && { echo "RUN $c2test FAILED!" ; cd - >/dev/null; return 1 ; }
mv LOG/$c2test.log LOG/${c2test}.log.$nextdate
echo $root/utilities/validate $c2test  LOG/$c2test.log.$nextdate
$root/utilities/validate $c2test  LOG/$c2test.log.$nextdate
ret=$?
echo ret is $ret
cd - >/dev/null
return $ret
}

build_code(){
  # Clean, configure and build the code
  cd $root/m4
  # Should be able to get away without making clean so that only new files are remade
  #make realclean
  ./configure -s > make.log.$thisrev 2>&1 || { echo "*** CONFIGURE FAILED!!!! ***" >>$log; cd ->/dev/null; return 1; }
  make >> make.log.$thisrev 2>&1 || { echo "***MAKE FAILED!!! ***" >> $log; cd - >/dev/null; return 1; }
  cd - >/dev/null
}


# If we don't have a specified end point, get the head of the repository
if [ $end -eq 0 ]
then
     end=`svn info -r HEAD http://ccpforge.cse.rl.ac.uk/svn/gamess-uk/trunk | grep "Revision:" | awk '{print $2}'`
fi
thisrev=$end

while [ 1 ]
  do

  # Get date of next change to repository
  let thisrev=($end-$start)/2+$start

  # See if we've got the version or something's gone wrong
  if [ $thisrev -lt $start ]
      then
      echo "ERROR GETTING NEXT REVISION: $nexrev"
      exit 1
  elif [ $thisrev -eq $start ]
      then
      echo " FOUND FAILING CHECKIN AT REVISION: $thisrev"| tee -a $log
      echo " See $log for more details" | tee -a $log
      exit 0
  fi

  echo -n "CHECKING REVISION $thisrev: " | tee -a $log

  # Need to remove the newscf directory
  rm -rf newscf

  # Update code to this date
  svn update -r $thisrev >> $log 2>&1


  # Hack need to use new configure and ignore any effects of the update
  cp $root/../configure .

  # Try to build the code
  build_code

  # See if it worked and if we need to bail
  if [ $? -ne 0 ]
      then
      echo "*** MAKE FAILED ***" | tee -a $log
      [ $makefailok -eq 0 ] && exit 1 || continue
  fi

   # Run the test
   echo "Running test..." >> $log
   $testcode >> $log 2>&1

   if [ $? -ne 0 ]
       then
       echo " *** TEST FAILED *** " | tee -a $log
       end=$thisrev
   else
       echo " ***   TEST OK   *** " | tee -a $log
       start=$thisrev
   fi

done
