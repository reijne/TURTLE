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

import sys
import string

# now follows a random string that can be changed to force subversion to
# update the version number of this file (sigh...): aaa

def version_toldiff(fp,errfp):
  """Print out the version information to the specified file object."""
  revision = "$Revision: 65 $"
  date     = "$Date: 2007-10-01 22:27:10 +0100 (Mon, 01 Oct 2007) $"
  rlist    = string.split(revision)
  dlist    = string.split(date)
  version  = "toldiff version "+rlist[1]+" of "+dlist[1]+"\n"
  try:
    fp.write(version)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing version information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def load_plain_text(fp,data,nline):
  """Reads a plain text file refered to by the file descriptor 'fp'. The 
     contents of the file are stored in the dictionary 'data'. The number of
     lines read are returned in 'nline'. 
     The line number is used as the key in the dictionary, the lines themselves
     are stored as the entries. The lines include the newline character. 
     Finally the first real line is stored as line number 1."""
  #
  # IOError exceptions are handled by the calling routine.
  #
  line = fp.readline()
  while line:
    nline = nline + 1
    data[nline] = line
    line = fp.readline()
  return (data,nline)

def load_tolerances(fp):
  """Reads a tolerance file from the file refered to by the file descriptor
     'fp'. The contents of the file are stored in three sections:

     Section 1: The change section
     -----------------------------

     Contains a specification of characters that are allowed to differ from the
     reference file. This information is stored in a dictionary of with line 
     numbers as keys. The format of the section is: 

     "<line number> <splice 1>,<splice 2>, ... ,<splice n>"

     Every splice is defined as "<char begin>:<char end>"

     Section 2: The add section
     --------------------------

     Contains a specification of the number of lines that might be added at
     a particular point in the reference file. The information will be stored
     in a dictionary with the reference file line numbers as key. Each entry
     records the maximum number of lines that may be added at that point. 
     The format of this section is:

     "<line number> <maximum number of added lines>"

     Section 3: The delete section
     -----------------------------

     Contains as specification of the lines in reference file that may be
     missing from the data file. The information will once again be stored in 
     a dictionary, but the entries in the dictionary will be empty. Instead
     the fact that a reference file line is mentioned in the dictionary alone
     means that it might be deleted in the data file. The format of this 
     section simply is:

     "<line number>"

     
     Finally, to separate the three sections the following constructs are used:
    
     <section type>
     ...
     end

     Where section type can be any of the keywords "change", "add" or "delete".

     Apart from the tolerance sections every tolerance file will also contain
     a version line. This line specifies the version of the toldiff script
     that last wrote the tolerance file. The version line is produced by
     the function version_toldiff and follows the associated format.
     """
  #
  # IOError exceptions are handled by the calling routine.
  #
  change = { }
  add    = { }
  delete = { }
  line = fp.readline()
  while line:
    list = string.split(line)
    if   list[0] == "change":
      line = fp.readline()
      while (line != "end\n"):
        list           = string.split(line)
        number         = list[0]
        tolerances     = list[1]
        number         = int(number)
        change[number] = tolerances
        line = fp.readline()
      line = fp.readline()
    elif list[0] == "add":
      line = fp.readline()
      while (line != "end\n"):
        list        = string.split(line)
        number      = list[0]
        numlines    = list[1]
        number      = int(number)
        numlines    = int(numlines)
        add[number] = numlines
        line = fp.readline()
      line = fp.readline()
    elif list[0] == "delete":
      line = fp.readline()
      while (line != "end\n"):
        number = int(line)
        delete[number] = None
        line = fp.readline()
      line = fp.readline()
    elif list[0] == "toldiff":
      # this is the version line which is currently not needed
      # so skip this.
      line = fp.readline()
    elif list[0] == "#":
      line = fp.readline()
    else:
      print line
      print "Unknown tolerance file section encountered\n"
      exit(10)
  return (change,add,delete)

def save_tolerances(fp,change,add,delete,errfp):
  """Saves the tolerances specified in the change, add and delete dictionaries.
     To file pointer fp is used to write the data to. See load_tolerances
     for a specification of the file format."""
  version_toldiff(fp,errfp)
  #
  # Beyond this point IOError exceptions are handled by the calling routine
  # above this one.
  #
  fp.write("change\n")
  list = change.keys()
  list.sort()
  for line in list:
    tol = change[line]
    if tol != "":
      fp.write("%d %s\n" % (line,tol))
  fp.write("end\n")
  fp.write("add\n")
  list = add.keys()
  list.sort()
  for line in list:
    numlines = add[line]
    fp.write("%d %d\n" % (line,numlines))
  fp.write("end\n")
  fp.write("delete\n")
  list = delete.keys()
  list.sort()
  for line in list:
    fp.write("%d\n" % line)
  fp.write("end\n")
