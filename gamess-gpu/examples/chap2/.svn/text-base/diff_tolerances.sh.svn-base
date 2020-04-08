#!/bin/sh
# This script compares the current tolerances against the previous tolerances
# as stored in SVN. To this end it uses the toldiff --show option to obtain
# the tolerances in a human readable form. Tkdiff is used to display the
# differences in context.
#
pid=$$
job=$1
reference_file=REF/${job}.log
current_tolerance_file=REF/${job}.tol
old_tolerance_file=/tmp/old_${job}_${pid}.tol
tolerance_diff=/tmp/tolerance_diff_${job}_${pid}.dff
current_tolerance_show=/tmp/current_${job}_${pid}.shw
old_tolerance_show=/tmp/old_${job}_${pid}.shw
toldiff=../../utilities/toldiff/toldiff.py
#
# Get the diff of the current tolerance file against the old one
#
svn diff ${current_tolerance_file} > ${tolerance_diff}
#
# Reconstruct the old tolerance file from the current one and the diff
# by reverse patching the current tolerance file.
#
cp -f ${current_tolerance_file} ${old_tolerance_file}
patch --reverse ${old_tolerance_file} ${tolerance_diff}
#
# Generate the tolerances in the context of the reference file to provide
# a human readable form of them.
#
$toldiff --show --tolerance ${old_tolerance_file} ${reference_file} > ${old_tolerance_show}
$toldiff --show --tolerance ${current_tolerance_file} ${reference_file} > ${current_tolerance_show}
#
# Display the differences between the tolerances.
#
diff ${old_tolerance_show} ${current_tolerance_show} > /dev/null
status=$?
if ( test \( ${status} -eq 1 \) )
then
    tkdiff ${old_tolerance_show} ${current_tolerance_show}
else
    echo "${job}: no change in tolerances"
fi
#
# Clean up
#
rm -f ${tolerance_diff} ${old_tolerance_file}
rm -f ${old_tolerance_show} ${current_tolerance_show}
exit 0
