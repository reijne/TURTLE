#!/usr/bin/python
#
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

import os
import sys
import string
import toldiff_files
import toldiff_lcs
import toldiff_diff
import toldiff_update
import toldiff_transfer
import toldiff_show

def license_toldiff(fp,errfp):
  """Print out the license information to the specified file object."""
  try:
    fp.write("""
    Copyright (C) 2006 Huub van Dam, Science and Technology Facilities Council,
    Daresbury Laboratory.
    All rights reserved.

    Developed by:        Huub van Dam
                         Science and Technology Facilities Council
                         Daresbury Laboratory
                         Computational Science and Engineering Department
                         Computational Chemistry Group
                         http://www.cse.clrc.ac.uk/ccg

    Permission is hereby granted, free of charge, to any person obtaining 
    a copy of this software and associated documentation files (the "Software"),
    to deal with the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    and/or sell copies of the Software, and to permit persons to whom the 
    Software is furnished to do so, subject to the following conditions: 

    Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimers. 
    Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimers in the documentation 
    and/or other materials provided with the distribution. 
    Neither the names of the Science and Technology Facilities Council,
    Daresbury Laboratory, the Computational Science and Engineering Department,
    the Computational Chemistry Group, nor the names of its contributors may be
    used to endorse or promote products derived from this Software without
    specific prior written permission. 

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS WITH THE SOFTWARE. 
    \n""")
    sys.exit(1)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing license information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def usage_toldiff(fp,errfp):
  """Print out the usage information to the specified file object."""
  try:
    fp.write("""
  Usage:

     toldiff [[--diff] <reference file> <data file>]
             [--update <reference file> <data file>]
             [--transfer <reference file> <new reference file>]
             [--show <reference file>]
             [--tolerance <tolerance file name>]
             [--new-tolerance <new tolerance file name>]
             [--diff-exe <diff executable>]
             [--output full|summary|none]
             [--summary <identical>:<equivalent>:<different>]
             [--exit <identical>:<equivalent>:<different>]
             [--[no]exact] [--[no]tolerant] [--[no]best]
             [--help] [--license] [--version]

  Toldiff is a script that compares two files allowing for tolerable
  differences. Tolerable differences often arise in meta data like who ran
  the test and on which date, timing data, and which machines and how many
  processors were used. In scientific/technical codes additional variations
  may result from numerical accuracy limitations.

  Toldiff is designed to assist in software testing by suppressing tolerable
  or trivial differences and highlighting only the significant ones. This
  facilitates checking whether an output a program has just produced matches
  the reference result obtained in the past.

  The toldiff script knows of the following files:

  A. The reference file:
     - THE correct file
  B. The data file: 
     - A file the correctness of which is to be tested against the reference
       file. Once its correctness has been established it may be used to update
       the tolerances.
  C. The tolerance file:
     - This file records where all allowed differences can occur, if any.
  D. The new reference file:
     - The file that is to replace the reference file after a change has
       taken place that outdates the reference file
  E. The new tolerance file:
     - This file records where allowed differences can occur relative to the
       new reference file instead of the current reference file.
    
  The script offers three processes:

  1. The diff process:
     - This process reports all differences between the reference file and
       the data file that are not explicitly tolerated.
  2. The update process:
     - This process updates the tolerances file adding all differences between 
       the reference file and the data file that were not tolerated before.
  3. The transfer process:
     - If the current reference file needs to be replaced by a new one this 
       process will carry as many as possible known tolerances relative to the
       current reference file over to the new reference file.

  There are various command line options to control toldiff. In cases where
  environment variables can be used as an alternative to command line options
  the precedence is handled as:
  - environment variables take precedence over default settings
  - command line options take precedence over environment variables.

  There are three categories of options this script will recognise:

  1. Process options:

    1.1 --diff <reference file name> <data file name>

        This triggers the script to perform the default diff process of
        comparing the data file against the reference file.

    1.2 --update <reference file name> <data file name>

        This requests the update process to be performed updating the
        tolerance file to allow for any differences between the reference and
        data files.

    1.3 --tolerance <tolerance file name>

        This option allows explicit specification of the tolerance file name.
        If omitted the script will construct a name for the tolerance file 
        from the name of the reference file.

    1.4 --transfer <reference file name> <new reference file name>

        This option invokes the transfer process to migrate as many tolerances
        as possible from the current reference file over to the new one.

    1.5 --new-tolerance <new tolerance file name>

        This option allows for the explicit specification of the name of the
        new tolerance file. If this is omitted the script will construct a name
        for the new tolerance file from the new reference file name.

    1.6 --diff-exe <diff executable>

        This option enbles replacing some of the Python diff implementation
        by invoking a binary diff program. This greatly improves the
        performance without changing the functionality. As an alternative 
        mechanism the environment variable TOLDIFF_EXE may be set to specify
        the diff program. In case both the command line option and the
        environment variable are provided the command line option has 
        precedence.

    1.7 --output full|summary|none

        This option controls the amount of output toldiff produces. The default
        setting "full" results in printing a full diff output. The setting
        "summary" suppresses the diff output and replaces it with a short 
        string for files being identical, equivalent or different. The values
        of these strings can be specified with the --summary option. Finally,
        setting "none" suppresses all output. Other than the --output option
        setting the TOLDIFF_OUTPUT environment variable does the same.

    1.8 --summary <identical>:<equivalent>:<different>

        This option allows the specification of short results for toldiff. The
        first string is reported if the reference file and data file are 
        identical. The second string is reported if the reference and data files
        are not identical but all differences are tolerated. The last string
        is reported if there are differences that are not tolerated. The
        default strings are "identical", "equivalent", and "different". Finally,
        these settings can be specified by setting the TOLDIFF_SUMMARY
        environment variable. In both ways the values are colomn separated.

    1.9 --exit <identical>:<equivalent>:<different>

        This option specifies the exit codes for toldiff. The first value is
        reported if the reference file and data file are identical. The second
        value is reported if the reference and data files are not identical but
        all differences are tolerated. The last value is reported if there are
        differences that are not tolerated. The default values are 0, 0, and 1.
        Finally, these settings can be specified by setting the TOLDIFF_EXIT
        environment variable. In both ways the values are colomn separated.

  2. Information options:

    2.1 --help

        Print this information on how to use this scripts.

    2.2 --show <reference file name>

        Prints the reference file marking all the known tolerances on it.
        This allows checking how the program has resolved differences through
        the tolerances chosen.
        The tolerances are marked on each line in the following order:
        1. The number of lines that may be inserted after this line.
        2. Whether this line may be deleted in which case it will be marked by
           a 'X', otherwise white space indicates that the line has to be 
           present.
        3. The contents of the line are shown with those characters that may
           change replaced by '#'.

    2.3 --version

        Print the version number of the toldiff script you are using.

    2.4 --license

        Print the license conditions under which this script is distributed.

  3. Debug options:

    These options are normally set automatically based on the requirements of
    the selected process. The default settings aim to complete the selected
    process with the highest efficiency. However, for debugging purposes it
    is possible to override these settings. You are free to try them to your
    own peril.

    3.1 --[no]exact

        Enable or disable the file differencing procedure that is based on
        exact line matches.

    3.2 --[no]tolerant

        Enable or disable the file differencing procedure that uses a line
        comparison which allows for tolerable differences between lines.

    3.3 --[no]best

        Enable or disable the file differencing procedure that matches lines
        based on maximum similarity.

  Copyright 2006, Huub van Dam, Science and Technology Facilities Council,
  Daresbury Laboratory\n""")
    sys.exit(1)
  except IOError, e:
    (errno,errmsg) = e
    try:
      errfp.write("toldiff: error writing usage information\n")
      errfp.write("toldiff: error message: ")
      errfp.write(errmsg)
      errfp.write("\n")
    except IOError, e:
      pass
    sys.exit(5)

def load_file(filename,err_fp):
  """Open and load a file. Returns the file text and the number of lines.
     The routine also handles I/O errors. I.e. it reports the error to the
     user and terminates the program."""
  text  = { }
  lines = 0
  try:
    file_fp = open(filename,"r")
    (text,lines) = toldiff_files.load_plain_text(file_fp,text,lines)
    file_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      err_fp.write("toldiff: I/O error on file: ")
      err_fp.write(filename)
      err_fp.write("\n")
      err_fp.write("toldiff: I/O error message: ")
      err_fp.write(errmsg)
      err_fp.write("\n")
    except IOError, e:
      pass
    sys.exit(10)
  return (text,lines)

def store_tolerance(tol_fnm,chg_txt,add_txt,del_txt,err_fp):
  """Open and write the tolerance file. The routine handles any I/O errors.
     I.e. it reports the error to the user and terminates the program."""
  try:
    tol_fp = open(tol_fnm,"w")
    toldiff_files.save_tolerances(tol_fp,chg_txt,add_txt,del_txt,err_fp)
    tol_fp.close()
  except IOError, e:
    (errno,errmsg) = e
    try:
      err_fp.write("toldiff: I/O error encountered attempting to write: ")
      err_fp.write(tol_fnm)
      err_fp.write("\n")
      err_fp.write("toldiff: I/O error message: ")
      err_fp.write(errmsg)
      err_fp.write("\n")
    except IOError, e:
      pass
    sys.exit(30)

def run_diff(diff_exe,ref_fnm,dat_fnm,fp):
  """This routine starts off an external diff program. It is assumed that
     the diff program ran successfully if the error file remains empty. 
     In this case the file descriptor of the standard output of the diff
     process is returned.
     If the error file is not empty then the contents of the error file is
     echoed to the user and the program terminates."""
  cmd = diff_exe+" "+ref_fnm+" "+dat_fnm
  try:
    (in_fp,out_fp,err_fp) = os.popen3(cmd)
  except IOError, e:
    (errno,errmsg) = e
    try:
      fp.write("toldiff: I/O error on external diff standard error file\n")
      fp.write("toldiff: I/O error message: ")
      fp.write(errmsg)
      fp.write("\n")
    except IOError, e:
      pass
    sys.exit(25)
  in_fp.close()
  return (out_fp,err_fp)
      
def find_overall_lcs(lexact,ltol,lbest,tol,ref_fnm,dat_fnm,diff_exe,err_fp):
  """Find the overall LCS including the tolerances. The general procedure is
     simply to establish the exact LCS, then try to resolve as much of the
     mismatches by considering the tolerances, then try to match the remaining
     differences to minimize the mismatches.

     This routine will read in the reference file and the data file as well.
     The reason for this is that this is more efficient in case an external
     diff program is used for the first phase.

     The routine returns the overall LCS, the reference file text, the data
     file text and beginning and ending line numbers of both files.

     This routine allows each phase to be disabled explicitly through a
     flag passed in as an argument:
     - lexact: if false skip the exact matching
     - ltol  : if false skip the tolerant matching
     - lbest : if false skip the minimal difference matching.
     """
  lcs = [ ]
  if lexact:
    if diff_exe == "":
      Nb = 1
      Mb = 1
      (ref,Ne) = load_file(ref_fnm,err_fp)
      (dat,Me) = load_file(dat_fnm,err_fp)
      lcs = toldiff_lcs.find_lcs1(ref,Nb,Ne,dat,Mb,Me)
    else:
      error = false
      (diff_out_fp,diff_err_fp) = run_diff(diff_exe,ref_fnm,dat_fnm,err_fp)
      (ref,Nb,Ne,dat,Mb,Me,lcs) = toldiff_diff.diff_to_lcs(diff_out_fp,err_fp)
      diff_out_fp.close()
      try:
        line = diff_err_fp.readline()
        while line:
          error = true
          err_fp.write("toldiff:"+line)
          line = diff_err_fp.readline()
        diff_err_fp.close()
      except IOError, e:
        (errno,errmsg) = e
        try:
          err_fp.write("toldiff: I/O error on external diff standard error file\n")
          err_fp.write("toldiff: I/O error message: ")
          err_fp.write(errmsg)
          err_fp.write("\n")
        except IOError, e:
          pass
        sys.exit(25)
      if error:
        sys.exit(20)
  else:
    Nb = 1
    Mb = 1
    (ref,Ne) = load_file(ref_fnm,err_fp)
    (dat,Me) = load_file(dat_fnm,err_fp)
  if (len(tol) <= 0) or (not ltol):
    #
    # No tolerances were specified or this phase is explicitly suppressed
    #
    pass
    #
  else:
    #
    # Consider all the differences and try to resolve as many as possible.
    #
    if (len(lcs) <= 0):
      #
      # Then the new LCS is simply the result of the tolerant diff
      #
      lcs = toldiff_lcs.find_lcs2(tol,ref,Nb,Ne,dat,Mb,Me)
      #
    else:
      #
      # First consider whether there is anything to compare before the first
      # snake
      #
      lcs1 = lcs
      (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
      if (xbot1 > Mb) and (ybot1 > Nb):
        lcs = toldiff_lcs.find_lcs2(tol,ref,Nb,ybot1-1,dat,Mb,xbot1-1)
      else:
        lcs = [ ]
      xtop0 = xtop1
      ytop0 = ytop1
      lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      while (len(lcs1) > 0 ):
        (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
        if (xbot1 > xtop0+1) and (ybot1 > ytop0+1):
          lcs2 = toldiff_lcs.find_lcs2(tol,ref,ytop0+1,ybot1-1,dat,xtop0+1,xbot1-1)
          lcs = lcs + lcs2
        xtop0 = xtop1
        ytop0 = ytop1
        lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      if (Ne >= ytop0+1) and (Me >= xtop0+1):
        #
        # The some more stuff at the end left to do
        #
        lcs2 = toldiff_lcs.find_lcs2(tol,ref,ytop0+1,Ne,dat,xtop0+1,Me)
        lcs = lcs + lcs2

  if (not lbest):
    #
    # This phase is explicitly suppressed
    #
    pass
    #
  else:
    #
    # Consider all the differences and try to match different lines as best as
    # possible minimizing the number of differences.
    #
    if (len(lcs) <= 0):
      #
      # Then the new LCS is simply the result of the best match diff,
      # which will probably hurt as this will get very expensive.
      #
      lcs = toldiff_lcs.find_lcs3(tol,ref,Nb,Ne,dat,Mb,Me)
      #
    else:
      #
      # First consider whether there is anything to compare before the first
      # snake
      #
      lcs1 = lcs
      (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
      if (xbot1 > Mb) and (ybot1 > Nb):
        lcs = toldiff_lcs.find_lcs3(tol,ref,Nb,ybot1-1,dat,Mb,xbot1-1)
      else:
        lcs = [ ]
      xtop0 = xtop1
      ytop0 = ytop1
      lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      while (len(lcs1) > 0 ):
        (xbot1,ybot1,xtop1,ytop1,type1) = lcs1.pop(0)
        if (xbot1 > xtop0+1) and (ybot1 > ytop0+1):
          lcs2 = toldiff_lcs.find_lcs3(tol,ref,ytop0+1,ybot1-1,dat,xtop0+1,xbot1-1)
          lcs = lcs + lcs2
        xtop0 = xtop1
        ytop0 = ytop1
        lcs.append((xbot1,ybot1,xtop1,ytop1,type1))
      if (Ne >= ytop0+1) and (Me >= xtop0+1):
        #
        # There is some more stuff at the end left to do
        #
        lcs2 = toldiff_lcs.find_lcs3(tol,ref,ytop0+1,Ne,dat,xtop0+1,Me)
        lcs = lcs + lcs2
  return (lcs,ref,Nb,Ne,dat,Mb,Me)

def construct_tolerance_filename(ref_fnm,dat_fnm,tol_fnm):
  if tol_fnm == "":
    #
    # No tolerance file name given so try and construct one
    #
    i = string.rfind(ref_fnm,".")
    if i == -1:
      #
      # The reference file name has no extension.
      #
      ref_ext = ""
      #
    else:
      #
      # The reference file name has an extension extract it.
      #
      ref_ext = ref_fnm[i:]
      #
    j = string.rfind(dat_fnm,".")
    if j == -1:
      #
      # The data file name has no extension.
      #
      dat_ext = ""
      #
    else:
      #
      # The data file name has an extension extract it.
      #
      dat_ext = dat_fnm[j:]
      #
    tol_ext = ".tol"
    if (tol_ext == ref_ext) or (tol_ext == dat_ext):
      tol_ext = ".tlr"
      if (tol_ext == ref_ext) or (tol_ext == dat_ext):
        tol_ext = ".tlc"
    tol_fnm = ref_fnm[:i]+tol_ext
  return tol_fnm

true = (0 == 0)
false = not true
#
# Set up the default comparison options
#
diff     = 1
update   = 2
transfer = 3
show     = 4
process  = diff
lexact   = true
ltol     = true
lbest    = false
#
# Set up default comparison results
#
identical  = 1
equivalent = 2
different  = 3
#
# Set up default comparison exit codes
#
exit_identical  = 0
exit_equivalent = 0
exit_different  = 1
#
# Set up default comparison summary texts
#
text_identical  = "identical"
text_equivalent = "equivalent"
text_different  = "different"
#
# Set up output options and default output option
#
output_full    = 3
output_summary = 2
output_none    = 1
output         = output_full

lcs = [ ]
diff_exe    = ""
tol_fnm     = ""
tol_new_fnm = ""
ref_fnm     = ""
dat_fnm     = ""
narg = len(sys.argv)
iarg = 1

if os.environ.has_key("TOLDIFF_EXE"):
  diff_exe = os.environ["TOLDIFF_EXE"]

if os.environ.has_key("TOLDIFF_OUTPUT"):
  output = os.environ["TOLDIFF_OUTPUT"]
  if output == "FULL" or output == "full":
    output = output_full
  elif output == "SUMMARY" or output == "summary":
    output = output_summary
  elif output == "NONE" or output == "none":
    output = output_none

if os.environ.has_key("TOLDIFF_EXIT"):
  exit_codes = os.environ["TOLDIFF_EXIT"]
  exit_codes = string.split(exit_codes,":")
  if len(exit_codes) == 3:
    exit_identical  = int(exit_codes[0])
    exit_equivalent = int(exit_codes[1])
    exit_different  = int(exit_codes[2])

if os.environ.has_key("TOLDIFF_SUMMARY"):
  text_summaries = os.environ["TOLDIFF_SUMMARY"]
  text_summaries = string.split(text_summaries,":")
  if len(text_summaries) == 3:
    text_identical  = text_summaries[0]
    text_equivalent = text_summaries[1]
    text_different  = text_summaries[2]

if narg == 1:
  usage_toldiff(sys.stdout,sys.stderr)
while iarg < narg:
  if sys.argv[iarg]   == "--exact":
    lexact = true
  elif sys.argv[iarg] == "--noexact":
    lexact = false
  elif sys.argv[iarg] == "--tolerant":
    ltol   = true
  elif sys.argv[iarg] == "--notolerant":
    ltol   = false
  elif sys.argv[iarg] == "--best":
    lbest  = true
  elif sys.argv[iarg] == "--nobest":
    lbest  = false
  elif sys.argv[iarg] == "--tolerance":
    iarg = iarg + 1
    if iarg < narg:
      tol_fnm = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing tolerance file name\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--new-tolerance":
    iarg = iarg + 1
    if iarg < narg:
      tol_new_fnm = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing new tolerance file name\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--diff-exe":
    iarg = iarg + 1
    if iarg < narg:
      diff_exe = sys.argv[iarg]
    else:
      try:
        sys.stderr.write("toldiff: missing diff executable specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--diff":
    process = diff
  elif sys.argv[iarg] == "--update":
    process = update
  elif sys.argv[iarg] == "--transfer":
    process = transfer
  elif sys.argv[iarg] == "--show":
    process = show
  elif sys.argv[iarg] == "--version":
    toldiff_files.version_toldiff(sys.stdout,sys.stderr)
    sys.exit(0)
  elif sys.argv[iarg] == "--help":
    usage_toldiff(sys.stdout,sys.stderr)
  elif sys.argv[iarg] == "--license":
    license_toldiff(sys.stdout,sys.stderr)
  elif sys.argv[iarg] == "--exit":
    iarg = iarg + 1
    if iarg < narg:
      exit_codes = sys.argv[iarg]
      exit_codes = string.split(exit_codes,":")
      if len(exit_codes) == 3:
        exit_identical  = int(exit_codes[0])
        exit_equivalent = int(exit_codes[1])
        exit_different  = int(exit_codes[2])
    else:
      try:
        sys.stderr.write("toldiff: missing exit codes specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--summary":
    iarg = iarg + 1
    if iarg < narg:
      text_summaries = sys.argv[iarg]
      text_summaries = string.split(text_summaries,":")
      if len(text_summaries) == 3:
        text_identical  = text_summaries[0]
        text_equivalent = text_summaries[1]
        text_different  = text_summaries[2]
    else:
      try:
        sys.stderr.write("toldiff: missing summaries specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  elif sys.argv[iarg] == "--output":
    iarg = iarg + 1
    if iarg < narg:
      output = sys.argv[iarg]
      if output == "FULL" or output == "full":
        output = output_full
      elif output == "SUMMARY" or output == "summary":
        output = output_summary
      elif output == "NONE" or output == "none":
        output = output_none
      else:
        sys.stderr.write("toldiff: unknown output specification: %s\n" % output)
        sys.exit(5)
    else:
      try:
        sys.stderr.write("toldiff: missing output specification\n")
      except IOError, e:
        pass
      sys.exit(5)
  else:
    argstr = sys.argv[iarg]
    if (process < show) and (iarg == narg-2):
      ref_fnm = sys.argv[iarg]
      iarg = iarg + 1
      dat_fnm = sys.argv[iarg]
    elif (process == show) and (iarg == narg-1):
      ref_fnm = sys.argv[iarg]
    elif argstr[0:1] == "-":
      try:
        sys.stderr.write("toldiff: unknow option encountered: ")
        sys.stderr.write(argstr)
        sys.stderr.write("\n")
      except IOError, e:
        pass
      sys.exit(8)
    else:
      sys.stderr.write("toldiff: missing reference or data files?\n")
      sys.exit(9)
  iarg = iarg + 1

if ref_fnm == "":
  sys.stderr.write("toldiff: error: no reference filename given\n")
  sys.exit(5)

if (process < show) and (dat_fnm == ""):
  sys.stderr.write("toldiff: error: no data filename given\n")
  sys.exit(6)

tol_fnm = construct_tolerance_filename(ref_fnm,dat_fnm,tol_fnm)
if process == transfer:
  tol_new_fnm = construct_tolerance_filename(dat_fnm,ref_fnm,tol_new_fnm)

ref_txt   = { }
dat_txt   = { }
chg_txt   = { }
add_txt   = { }
del_txt   = { }
ref_lines = 0
dat_lines = 0

try:
  tol_fp  = open(tol_fnm,"r")
  (chg_txt,add_txt,del_txt) = toldiff_files.load_tolerances(tol_fp)
  tol_fp.close()
except IOError, e:
  #
  # If an exception was thrown it is assumed that there is no valid 
  # tolerance file present. Hence proceed as if there is no tolerance 
  # information.
  #
  pass

if process == diff:
  (lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines) = find_overall_lcs(lexact,ltol,lbest,chg_txt,ref_fnm,dat_fnm,diff_exe,sys.stderr)
  lcs = toldiff_lcs.filter_lcs(lcs,Nb,ref_lines,Mb,dat_lines,add_txt,del_txt)
  analysis = toldiff_diff.lcs_analysis(Nb,ref_lines,Mb,dat_lines,lcs,identical,equivalent,different)
  if output == output_full:
    toldiff_diff.lcs_to_diff(ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,lcs,sys.stdout,sys.stderr)
  elif output == output_summary:
    if   analysis == identical:
      sys.stdout.write("%s" % text_identical)
    elif analysis == equivalent:
      sys.stdout.write("%s" % text_equivalent)
    elif analysis == different:
      sys.stdout.write("%s" % text_different)
    else:
      sys.stderr.write("illegal value of analysis")
  elif output == output_none:
    pass
  else:
    sys.stderr.write("illegal value of output")
  if   analysis == identical:
    sys.exit(exit_identical)
  elif analysis == equivalent:
    sys.exit(exit_equivalent)
  elif analysis == different:
    sys.exit(exit_different)
  else:
    sys.stderr.write("illegal value of analysis")
elif process == update:
  (lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines) = find_overall_lcs(true,true,true,chg_txt,ref_fnm,dat_fnm,diff_exe,sys.stderr)
  chg_txt = toldiff_update.lcs_to_change(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,chg_txt)
  add_txt = toldiff_update.lcs_to_add(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,add_txt)
  del_txt = toldiff_update.lcs_to_delete(lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines,del_txt)
  store_tolerance(tol_fnm,chg_txt,add_txt,del_txt,sys.stderr)
elif process == transfer:
  (lcs,ref_txt,Nb,ref_lines,dat_txt,Mb,dat_lines) = find_overall_lcs(true,true,false,chg_txt,ref_fnm,dat_fnm,diff_exe,sys.stderr)
  (chg_new,add_new,del_new) = toldiff_transfer.transfer_tol(lcs,Nb,ref_lines,Mb,dat_lines,chg_txt,add_txt,del_txt)
  store_tolerance(tol_new_fnm,chg_new,add_new,del_new,sys.stderr)
elif process == show:
  Nb = 1
  (ref_txt,Ne) = load_file(ref_fnm,sys.stderr)
  toldiff_show.show_tolerance(sys.stdout,ref_txt,Nb,Ne,chg_txt,add_txt,del_txt,sys.stderr)
else:
  try:
    sys.stderr.write("toldiff: internal error: invalid process")
  except IOError, e:
    pass
  sys.exit(999)
