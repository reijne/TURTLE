#!/usr/bin/python
#
# This is a little program to rewrap code that Maple has generated for
# Density Functional Repository.
#
# Maple will wrap lines to 72 characters irrespective of the contents. 
# In particular Maple may wrap lines in the middle of numeric constants or
# identifiers. For the Fortran compilers this is not a problem but for 
# preprocessors that need to expand macro's and such this is a complete 
# disaster.
#
# This script unwraps the Maple lines (which use "#" as the continuation
# character) and rewraps them at sensible positions. The script works like
# a filter, i.e. it reads the input file from standard input and writes its
# results to standard output.
#
import sys
from string import rfind

def rewrap_line(longline):
  # note: longline has a newline character at the end which is counted by 
  #       len as well but is not included in the 72 characters that are allowed
  #       in Fortran
  while len(longline)-1 > 72:
    i = -1
    # wrap before * / ( ) + or -
    i = max(i,rfind(longline,"*",0,71))
    i = max(i,rfind(longline,"/",0,71))
    i = max(i,rfind(longline,"(",0,71))
    i = max(i,rfind(longline,")",0,71))
    i = max(i,rfind(longline,"+",0,71))
    # wrap before - but not in the middle of a numerical constant...
    j = rfind(longline,"-",0,71)
    k = rfind(longline,"D-",0,71)
    if j-1 == k:
      j = rfind(longline,"-",0,k)
    i = max(i,j)
    if i == -1:
      sys.stderr.write("No sensible break point found in:\n")
      sys.stderr.write(longline)
      exit(1)
    sys.stdout.write(longline[:i]+"\n")
    longline = "     #" + longline[i:]
  sys.stdout.write(longline)

longline = ""
line = sys.stdin.readline()
while line:
  if line[:6] == "     #":
    length = len(longline)
    if len > 0:
      longline = longline[:length-1]
    longline += line[6:]
  else:
    rewrap_line(longline)
    longline = line
  line = sys.stdin.readline()
rewrap_line(longline)

