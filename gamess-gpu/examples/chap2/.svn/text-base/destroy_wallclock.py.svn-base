#!/usr/bin/env python
#
# destroy_wallclock
# =================
#
# This is a script to find the wall clock time prints in a GAMESS-UK output
# and replace them with garbage so that the tolerance file can be update to
# ignore these bloody timings.
# 
import os
import sys

narg = len(sys.argv)

input_fnm  = "-"
output_fnm = "-"

if narg >= 2:
  input_fnm = sys.argv[1]
if narg == 3:
  output_fnm = sys.argv[2]

if input_fnm == "-":
  fpin = sys.stdin
else: 
  fpin = open(input_fnm,"r")
data = fpin.readlines()
fpin.close()

nlines = len(data)
iline = 0
while (iline < nlines):
  if data[iline].find("diagonalizer timing analysis") > -1:
    #
    # Found integral/interaction timing analysis, now deal with the
    # wall clock times
    #
    tline = iline + 5
    while (data[tline].find("other") == -1):
      line = data[tline][0:26]+"###.###"+data[tline][33:43]+"#.###\n"
      data[tline] = line
      tline = tline + 1
    line = data[tline][0:26]+"###.###\n"
    data[tline] = line
  if data[iline].find("gamess timing analysis") > -1:
    #
    # Found a timing output block, now deal with the wall clock times
    #
    tline = iline + 3
    while (data[tline].find("****") == -1):
      line = data[tline][0:39]+"#.##"+data[tline][43:47]+"###.##"+data[tline][53:65]+"#.##"+data[tline][69:73]+"###.##\n"
      data[tline] = line
      tline = tline + 1
    tline = tline + 1
    line = data[tline][0:39]+"#.##"+data[tline][43:65]+"#.##\n"
    data[tline] = line
  iline = iline + 1

if output_fnm == "-":
  fpout = sys.stdout
else:
  fpout = open(output_fnm,"w")
iline = 0
while (iline < nlines):
  fpout.write(data[iline])
  iline = iline + 1
fpout.close()
