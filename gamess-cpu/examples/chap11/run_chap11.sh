#!/bin/bash

echo "Running chap11"

export PATH=.:$PATH
date
for job in `cat ./jobs.list | grep -v "#"`
do

  # Run the test case
  ./RUN ./$job

  # Delete any files
  if [ $job = "c11000a" ]
      then
      /bin/rm -f hclmain hcldump
  elif [ $job = "c11001a" ]
      then
      /bin/rm -f mfged2 mfged3
  elif [ $job = "c11006" ]
      then
      /bin/rm -f namgmain namgdump
  elif [ $job = "c11007" ]
      then
      /bin/rm -f nicomain nicodump
  elif [ $job = "c11008" ]
      then
      /bin/rm -f gvbmain gvbdump
  elif [ $job = "c11011_a" ]
      then
      /bin/rm -f hpsimain hpsidump
  elif [ $job = "c11011_b" ]
      then
      /bin/rm -f  mfged2 mfged3
  elif [ $job = "c11012_a" ]
      then
      /bin/rm -f  c2h4main c2h4dump
  elif [ $job = "c11012_b" ]
      then
      /bin/rm -f  c2h4main c2h4dump
  elif [ $job = "c11012_c" ]
      then
      /bin/rm -f  c2h4main c2h4dump
  elif [ $job = "c11012_d" ]
      then
      /bin/rm -f  pyred
  elif [ $job = "c11016" ]
      then
      /bin/rm -f h2oed1 h2omain h2odump h2oed4 h2oed6 h2oed9 h2oed10 h2oed11
  elif [ $job = "c11017" ]
      then
      /bin/rm -f beomain beodump beoed1 beoed4 beoed5 beoed6 beoed9
  elif [ $job = "c11018" ]
      then
      /bin/rm -f beoed10 beomain beodump beoed4 beoed6 beoed13
  elif [ $job = "c11018a" ]
      then
      /bin/rm -f beoed10 beomain beodump beoed4 beoed6 beoed13
  elif [ $job = "c11020" ]
      then
      /bin/rm -f nh3main nh3dump nh3tran nh3sel01 nh3sel02 nh3ham nh3diag
  elif [ $job = "c11023" ]
      then
      /bin/rm -f  pyred2 pyred3 pyrtran
  elif [ $job = "c11024" ]
      then
      /bin/rm -f  pyred2 pyred3 pyrtran
  elif [ $job = "c11025" ]
      then
      /bin/rm -f  pyred2 pyred3 pyrtran
  elif [ $job = "c11025a" ]
      then
      /bin/rm -f  pyred2 pyred3 table
  elif [ $job = "c11026" ]
      then
      /bin/rm -f  h2omain h2odump file8
  fi
done

/bin/rm -f mfged* table ftn058
/bin/rm -f mged2 mged3 pyred2 pyred3
date

# Now validate & check for any errors
export GAMESS_VTAB=./chap11.vtab

error=0
for test in `cat ./jobs.list | grep -v "#"`
do
    ../../utilities/validate $test LOG/${test}.log || error=1
done

exit $error
