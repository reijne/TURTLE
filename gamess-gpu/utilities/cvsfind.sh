#!/bin/bash
#
# bash script to check when a particular change was introduced into the code
#
# The script checks out the code on each day that the code was altered 
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


nextdate=2007-09-01  # Set date to look for change from (format: YYYY-MM-DD)
testcode=run_chap2_eg # The test to run
makefailok=1          # Are makes allowed to fail?
root=`pwd`            # The root of the GAMESS-UK directory

# The name of the file containing the list of dates - if this doesn't
# exist it is created from the cvs history
datelist=$root/date.list
log=$root/check.log; > $log


get_next_date(){
    # Given a date, find the first change to the repository after this date
    local startdate=$1

    # If the history file doesn't exist, create it
    [ ! -f $datelist ] && cvs history -cal -D $startdate | grep GAMESS-UK | awk '{print $2}' | sort -u > $datelist

    # Get the next date from the file
    nextdate=`head -1 $datelist`

    # Need to check we got a date
    check_date $nextdate
    [ $? -ne 0 ] && { echo "get_next_date failed at check_date!!!!"; nextdate=""; return 1; }

    #[ $nextdate = $startdate ] && { echo "NEXTDATE IS SAME AS STARTDATE!!!"; nextdate=""; return 1; }

    # Date is o.k. so remove it from the file
    mv $datelist ${datelist}.old
    sed '1 d' < ${datelist}.old > $datelist
    rm  ${datelist}.old
}

check_date(){
# Check a date is valid
[ $# -ne 1 ] && echo "check_date needs a date!" && return 1
[ ! $1 ] && return 1
echo $1 | grep "[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}" > /dev/null 2>&1
}

run_chap2_eg(){
# See if a particular chap2 example is validating or not

# Specify the test to run
local c2test=c2020_i

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
  ./configure -s > make.log.$nextdate 2>&1 || { echo "*** CONFIGURE FAILED!!!! ***" >>$log; cd ->/dev/null; return 1; }
  make >> make.log.$nextdate 2>&1 || { echo "***MAKE FAILED!!! ***" >> $log; cd - >/dev/null; return 1; }
  cd - >/dev/null
}

while [ 1 ]
  do

  # Get date of next change to repository
  get_next_date $nextdate

  [ ! $nextdate ] && echo "COULD NOT FIND NEXT DATE!!!" | tee -a $log && exit 1

  echo -n "CHECKING DATE $nextdate: " | tee -a $log

  # Update code to this date
  cvs update -d -R -D $nextdate >> $log 2>&1 

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
  [ $? -ne 0 ] && { echo " *** TEST FAILED *** - see file $log for more details" | tee -a $log ; exit 1; } || { echo " ***   test ok   *** " | tee -a $log ; }

done
