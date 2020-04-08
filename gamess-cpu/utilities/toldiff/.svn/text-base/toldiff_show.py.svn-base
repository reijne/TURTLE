# Copyright (C) 2006 Huub van Dam, Science and Technology Facilities Council,
# Daresbury Laboratory.
# All rights reserved.
#
# Developed by:        Huub van Dam
#                      Science and Technology Facilities Council
#                      Daresbury Laboratory
#                      Computational Science and Engineering Department
#                      Computational Chemistry Group
#                      http://www.cse.clrc.ac.uk/ccg
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"),
# to deal with the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions: 
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimers. 
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution. 
# Neither the names of the Science and Technology Facilities Council,
# Daresbury Laboratory, the Computational Science and Engineering Department,
# the Computational Chemistry Group, nor the names of its contributors may be
# used to endorse or promote products derived from this Software without
# specific prior written permission. 
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS WITH THE SOFTWARE. 

import toldiff_lcs
import sys
import copy
import math
import string

def show_tolerance(out_fp,ref_txt,Nb,Ne,chg_txt,add_txt,del_txt,ferr):
  """Prints out the reference with the tolerance marked on it. This function
     is essentially an informational aid to help the user understand what
     toldiff is doing to his data.

     For every line we print
     1. the number of lines that may be added after it if that exceeds zero
     2. a 'X' if the line may be deleted and a ' ' otherwise
     3. the line of the reference with the characters that may change replaced
        by a '#'"""
  try:
    #
    # First find out how many characters are needed to print the numbers of
    # potentially added lines. Also find out whether the first addition might
    # appear after line zero.
    #
    line_zero = add_txt.has_key(Nb-1)
    maxno = 1
    keys = add_txt.keys()
    while (len(keys) > 0):
      key = keys.pop()
      if add_txt[key] > maxno:
        maxno = add_txt[key]
    numdigit = int(1.0+math.log10(maxno))
    #
    # Use numdigit to create the right format string
    #
    fmtn = "%"+str(numdigit)+"d%c%s"
    fmte = ""
    i = 0
    while (i < numdigit):
      fmte = fmte + " "
      i = i + 1
    fmte = fmte + "%c%s"
    if line_zero:
      nadd = add_txt[Nb-1]
      cdel = " "
      line = "\n"
      out_fp.write(fmtn % (nadd,cdel,line))
    #
    # Loop over all the lines in the reference and sort it out
    #
    i = Nb
    while (i <= Ne):
      if chg_txt.has_key(i):
        line = copy.deepcopy(ref_txt[i])
        line = toldiff_lcs.tol_decode(chg_txt[i],line)
        if string.find(line,"\n") == -1:
          line = line + "\n"
      else:
        line = ref_txt[i]
      if del_txt.has_key(i):
        cdel = "X"
      else:
        cdel = " "
      if add_txt.has_key(i):
        nadd = add_txt[i]
        out_fp.write(fmtn % (nadd,cdel,line))
      else:
        out_fp.write(fmte % (cdel,line))
      i = i + 1
  except IOError, e:
    (errno,errmsg) = e
    try:
      ferr.write("toldiff: error writing the marked reference file\n")
      ferr.write("toldiff: error message: ")
      ferr.write(errmsg)
      ferr.write("\n")
    except IOError, e:
      pass
    sys.exit(5)
