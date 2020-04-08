#!/bin/sh
echo "AUTO_UPDATE"
echo "==========="
echo
echo "This script automatically runs all test cases and updates the tolerance"
echo "files. This eliminates the critical inspection of the validation results"
echo "in favour of quickly discovering trivial differences."
echo
echo "This procedure is only safe if the exact executable used has been"
echo "validated correctly once before. Only then is it safe to assume that any"
echo "differences discovered in this procedure are more than likely to be"
echo "trivial."
echo
echo "To make absolutely sure that no rogue tolerances are introduced the"
echo "script ./diff_all_tolerances.sh graphically diffs the new tolerances"
echo "against the old ones. This allows you to check that new tolerances only"
echo "appear in places where you expect them."
echo
echo "You are strongly encouraged to perform this check before committing"
echo "any updated tolerance files."
echo
echo "Finally, the script writes the output the run script and the validation"
echo "script to a log file every iteration. If there are any problems these"
echo "log files are there to help you."
echo
echo "Are you sure you want to proceed? (y/n/m)"

read input
case "$input" in
  y) echo "Yes. Starting auto update in 30 seconds..." 
     sleep 30 ;;
  n) echo "No. Aborting auto update..." 
     exit 0 ;;
  m) echo "Maybe? Maybe is not nearly good enough!"
     echo "Make sure you understand exactly what this procedure means!"
     exit 1 ;;
  *) echo "Please answer y, n, or m"
     exit 1 ;;
esac

iteration=1

./run_chap2.sh 2>&1 > "./result_auto_update_${iteration}.log"
until (./val_chap2.sh -t 2>&1 >> "./result_auto_update_${iteration}.log")
do
  yes | ./update_chap2.sh  > /dev/null
  iteration=$((iteration+1))
  ./run_chap2.sh 2>&1 > "./result_auto_update_${iteration}.log" 
done
