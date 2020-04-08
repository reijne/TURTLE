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
import copy
import string

true = (0 == 0)
false = not true

class search_path:
  """The search_path class is used to record the progress of a search
     path through the two files. There are three components that should be 
     tracked
     1. A list of snakes each snake is recorded by its starting and end point.
     2. As the last point is not recorded if there is no match the last point
        explored in the search path needs to be recorded separately.
     3. For the purpose of eliminating failing paths it is useful to keep 
        track of the length longest common subsequence sofar.
     """
  def __init__(self):
    """Initialize the search path.
       1. Create an empty snake list.
       2. Create an initial last point (although this will require resetting
          for backward searches).
       3. Set the length sofar to zero.
       """
    self.snakelist = []
    self.lastpoint = (0,0)
    self.lcslength = 0

  def set_lastpoint(self,i,j):
    """Set the last point searched to (i,j)"""
    self.lastpoint = (i,j)

  def increase_length(self,n=1):
    """Increase the LCS length sofar by n"""
    self.lcslength = self.lcslength + n

  def add_snake(self,i1,j1,i2,j2,itype):
    """Add a new snake defined by points (i1,j1) and (i2,j2)
       also store the type of the snake."""
    self.snakelist.append((i1,j1,i2,j2,itype))

class search_path_linked:
  """The search_path class is used to record the progress of a search
     path through the two files. However the search path class as implemented
     has as a severe disadvantage that if a search splits into two different
     directions a new search has to be created that is a copy of the one it 
     split off from. This copy operation is extremely expensive. 

     So in the search_path_linked implementation the snake list is stored as
     a linked list instead. Every snake links back to the previous snake in the
     search path. Thus avoiding the whole snakelist to be copied if a search
     splits. 

     To generate the final snake list that is returned from the search the
     linked snake list will need to be transformed.

     The linked search path has three components that should be tracked
     1. A linked list of snakes each snake is recorded by its starting and end
        point.
     2. As the last point is not recorded if there is no match the last point
        explored in the search path needs to be recorded separately.
     3. For the purpose of eliminating failing paths it is useful to keep 
        track of the length longest common subsequence sofar.
     """
  def __init__(self):
    """Initialize the search path.
       1. Create an empty linked snake list.
       2. Create an initial last point (although this will require resetting
          for backward searches).
       3. Set the length sofar to zero.
       """
    self.snakelist = None
    self.lastpoint = (0,0)
    self.lcslength = 0

  def set_lastpoint(self,i,j):
    """Set the last point searched to (i,j)"""
    self.lastpoint = (i,j)

  def increase_length(self,n=1):
    """Increase the LCS length sofar by n"""
    self.lcslength = self.lcslength + n

  def add_snake(self,i1,j1,i2,j2,itype,previous_snake=None):
    """Add a new snake defined by points (i1,j1) and (i2,j2)
       also store the type of the snake."""
    if (previous_snake == None):
      self.snakelist = (i1,j1,i2,j2,itype,self.snakelist)
    else:
      self.snakelist = (i1,j1,i2,j2,itype,previous_snake)

def min(a,b):
  """Return the minimum of two values"""
  if (a<b):
    return a
  else:
    return b

def max(a,b):
  """Return the maximum of two values"""
  if (a>b):
    return a
  else:
    return b

def transform_snake_list(linked_snake):
  """Transforms a linked snake list to a canonical snake list, i.e.
     from something of the form:

       (x1,y1,x2,y2,itype,previous_snake)

     to something of the form:

       [ ..., (x1,y1,x2,y2,ityp), ... ]

     This is needed to go from the linked list snake search to a snake list
     that is more suitable to other parts of the program."""

  lcs = []
  while (linked_snake):
    (x1,y1,x2,y2,itype,linked_snake) = linked_snake
    lcs.insert(0,(x1,y1,x2,y2,itype))
  return lcs
    
def find_lcs1(ref,Nb,Ne,dat,Mb,Me):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns the longest common subsequence (LCS) in 'lcs'. The LCS
     is stored as a list of snakes. A snake is a sequence of line pairs
     (Xi,Yj) to (Xi+p,Yj+p) where the lines X and Y in every pair match.
     Whatever happens between two snakes in a path is irrelevant.
     As this routine looks for exact matches it produces type 1 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however I cannot guarantee that understood it well enough to reproduce
     the actual published algorithm.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs = { }
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # NF - counter for generating forward search path keys
  #NF = 1
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  # the length of the longest LCS sofar
  max_len_lcs = 0

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      opt_len = s.lcslength + min(Me-xi+1,Ne-yj+1)
      if (opt_len > max_len_lcs):
        #
        # There (still) is hope that this search can beat the best search sofar
        #
        #
        # First try whether we are onto a snake
        #
        xi1 = xi+1
        yj1 = yj+1
        if yj1 <= Ne and xi1 <= Me:
          if ref[yj1] == dat[xi1]:
            # yes, we are onto a snake
            xi2 = xi1 + 1
            yj2 = yj1 + 1
            while (yj2 <=Ne and xi2 <= Me):
              if ref[yj2] == dat[xi2]:
                xi2 = xi2 + 1
                yj2 = yj2 + 1
              else:
                break
            xi2 = xi2 - 1
            yj2 = yj2 - 1
            s.add_snake(xi1,yj1,xi2,yj2,1)
            s.increase_length(xi2-xi1+1)
            s.set_lastpoint(xi2,yj2)
            xi = xi2
            yj = yj2
            finished = (yj2 == Ne and xi2 == Me)
            if finished:
              lcs = transform_snake_list(s.snakelist)
              break
        #
        # update the maximum LCS length
        #
        max_len_lcs = max(max_len_lcs,s.lcslength)
        #
        # now explore the way forward, horizontal first
        #
        keep_horizontal = false
        xih = xi+1
        yjh = yj
        if xih <= Me:
          if FP.has_key((xih,yjh)):
            if FP[(xih,yjh)] < s.lcslength:
              keep_horizontal = true
          else:
            keep_horizontal = true
          if xih+yjh < flij:
            flij = xih+yjh
          finished = (yjh == Ne and xih == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        #
        # now explore the vertical direction
        #
        keep_vertical = false
        xiv = xi
        yjv = yj+1
        if yjv <= Ne:
          if FP.has_key((xiv,yjv)):
            if FP[(xiv,yjv)] < s.lcslength:
              keep_vertical = true
          else:
            keep_vertical = true
          if xiv+yjv < flij:
            flij = xiv+yjv
          finished = (yjv == Ne and xiv == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        if keep_vertical:
          if keep_horizontal:
            # Keeping both horizontal and vertical search direction
            # So generate a new search path
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
            FP[(xiv,yjv)] = sa.lcslength
            #
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping only the vertical search direction
            # So simply update the current search path
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
            FP[(xiv,yjv)] = s.lcslength
        else:
          if keep_horizontal:
            # Keeping only the horizontal search direction
            # So simply update the current search path
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping neither the horizontal or vertical search direction
            # So remove the current path from the search list
            pass
      else:
        pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
  return lcs

def tol_decode(tol_line,txt_line):
  """Decode a single line of the tolerance change section and mark the
     relevant characters in the text line. Then return the modified 
     text line.
     """
  clist = string.split(tol_line,',')
  c1 = 0
  c2 = len(clist)
  while (c1 < c2):
    p = clist[c1]
    (s1,s2) = string.split(p,':')
    s1 = int(s1)
    s2 = int(s2)
    while (s1 <= s2):
      txt_line = txt_line[:s1] + "#" + txt_line[s1+1:]
      s1 = s1 + 1
    c1 = c1 + 1
  return txt_line

def tol_encode(txt_line):
  """We assume that the text line contains hash character in all places where
     differences should be tolerated. This function will find the hash 
     characters and encode their locations in the format of a tolerance 
     change section line.
     """
  #
  # s1,s2 - beginning and end of a tolerance
  # intol - whether we are in a tolerance part
  # ntol  - the number of tolerance parts found
  #
  tol_line = ""
  true  = (0==0)
  false = not true
  ntol  = 0
  i     = 0
  n     = len(txt_line)
  intol = false
  while (i < n):
    if txt_line[i] == "#":
      if intol:
        s2 = i
      else:
        s1 = i
        s2 = i
        intol = true
        ntol  = ntol + 1
    else:
      if intol:
        if ntol > 1:
          tol_line = tol_line + ","
        tol_line = tol_line + str(s1)+":"+str(s2)
        intol = false
      else:
        pass
    i = i + 1
  if intol:
    if ntol > 1:
      tol_line = tol_line + ","
    tol_line = tol_line + str(s1)+":"+str(s2)
  return tol_line

def tol_compare(tol,ref,yj,dat,xi):
  """Compares two lines of text taking into account tolerable differences
     stored in 'tol' if any.
     """
  if tol.has_key(yj):
    tmp_ref = copy.deepcopy(ref[yj])
    tmp_dat = copy.deepcopy(dat[xi])
    tmp_ref = tol_decode(tol[yj],tmp_ref)
    tmp_dat = tol_decode(tol[yj],tmp_dat)
    result = (tmp_ref == tmp_dat)
  else:
    result = (ref[yj] == dat[xi])
  return result

def find_lcs2(tol,ref,Nb,Ne,dat,Mb,Me):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns the longest common subsequence (LCS) in 'lcs'. The LCS
     is stored as a list of snakes. A snake is a sequence of line pairs
     (Xi,Yj) to (Xi+p,Yj+p) where the lines X and Y in every pair match.
     Whatever happens between two snakes in a path is irrelevant.
     In this particular routine the string comparison is modified based
     on the information held 'tol'. For every relevant line the tol
     dictionary holds a string of splices of characters where differences
     should be tolerated. As this routine uses a tolerant comparison it
     generates type 2 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however I cannot guarantee that understood it well enough to reproduce
     the actual published algorithm.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs = { }
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # NF - counter for generating forward search path keys
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  # the length of the longest LCS sofar
  max_len_lcs = 0

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      opt_len = s.lcslength + min(Me-xi+1,Ne-yj+1)
      if (opt_len > max_len_lcs):
        #
        # There (still) is hope that this search can beat the best search sofar
        #
        # First try whether we are onto a snake
        #
        xi1 = xi+1
        yj1 = yj+1
        if yj1 <= Ne and xi1 <= Me:
          # line comparison 1
          if tol_compare(tol,ref,yj1,dat,xi1):
            # yes, we are onto a snake
            xi2 = xi1 + 1
            yj2 = yj1 + 1
            while (yj2 <=Ne and xi2 <= Me):
              # line comparison 2
              if tol_compare(tol,ref,yj2,dat,xi2):
                xi2 = xi2 + 1
                yj2 = yj2 + 1
              else:
                break
            xi2 = xi2 - 1
            yj2 = yj2 - 1
            s.add_snake(xi1,yj1,xi2,yj2,2)
            s.increase_length(xi2-xi1+1)
            s.set_lastpoint(xi2,yj2)
            xi = xi2
            yj = yj2
            finished = (yj2 == Ne and xi2 == Me)
            if finished:
              lcs = transform_snake_list(s.snakelist)
              break
        #
        # update the maximum LCS length
        #
        max_len_lcs = max(max_len_lcs,s.lcslength)
        #
        # now explore the way forward, horizontal first
        #
        keep_horizontal = false
        xih = xi+1
        yjh = yj
        if xi1 <= Me:
          if FP.has_key((xih,yjh)):
            if FP[(xih,yjh)] < s.lcslength:
              keep_horizontal = true
          else:
            keep_horizontal = true
          if xih+yjh < flij:
            flij = xih+yjh
          finished = (yjh == Ne and xih == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        #
        # now explore the vertical direction
        #
        keep_vertical = false
        xiv = xi
        yjv = yj+1
        if yjv <= Ne:
          if FP.has_key((xiv,yjv)):
            if FP[(xiv,yjv)] < s.lcslength:
              keep_vertical = true
          else:
            keep_vertical = true
          if xiv+yjv < flij:
            flij = xiv+yjv
          finished = (yjv == Ne and xiv == Me)
          if finished:
            lcs = transform_snake_list(s.snakelist)
            break
        if keep_vertical:
          if keep_horizontal:
            # Keeping both horizontal and vertical search directions
            # So generate a new search path
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
            FP[(xiv,yjv)] = sa.lcslength
            #
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping only the vertical search direction
            # So simply update the current search path
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
            FP[(xiv,yjv)] = s.lcslength
        else:
          if keep_horizontal:
            # Keeping only the horizontal search direction
            # So simply update the current search path
            s.set_lastpoint(xih,yjh)
            FV[(xih,yjh)] = s
            FP[(xih,yjh)] = s.lcslength
          else:
            # Keeping neither the horizontal or vertical search direction
            # So simply let the current search path die
            pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
  return lcs


def num_same(tol,ref,yj,dat,xi):
  """Compares two lines of text taking into account tolerable differences
     stored in 'tol' if any and computing the number of characters that 
     are the same in the two lines.
     """
  if tol.has_key(yj):
    tmp_ref = copy.deepcopy(ref[yj])
    tmp_dat = copy.deepcopy(dat[xi])
    tmp_ref = tol_decode(tol[yj],tmp_ref)
    tmp_dat = tol_decode(tol[yj],tmp_dat)
  else:
    tmp_ref = ref[yj]
    tmp_dat = dat[xi]
  len_ref = len(tmp_ref)
  len_dat = len(tmp_dat)
  if len_ref < len_dat:
    length = len_ref
  else:
    length = len_dat
  result = 0
  i = 0
  while (i < length):
    if tmp_ref[i] == tmp_dat[i]:
      result = result + 1
    i = i + 1
  return result

def find_lcs3(tol,ref,Nb,Ne,dat,Mb,Me):
  """Compares the data stored in 'dat' against the data in 'ref', 
     and returns a common subsequence where lines are matched up so as
     to maximise the number of equal characters. The common subsequence
     is stored as a list of snakes. As the matching is based on a very
     much relaxed line comparison it generates type 3 snakes.

     The algorithm used here is inspired by:

       E. Myers, 'An O(ND) Difference Algorithm and Its Variations'
       Algorithmica 1, 2 (1986), 251-266
       http://www.cs.arizona.edu/people/gene/PAPERS/diff.ps

     however as the algorithm has to find a maximum number of matching 
     characters the algorithm has to be substantially modified. In particular
     the greedy approach does not work, neither does the first across the line
     criterion.

     The algorithm we use is to progress in three directions from every point.
     1. Horizontally,
     2. Vertically
     3. Diagonally
     As diagonal moves are the only ones that can up the similarity score it
     is expected that they are most likely to survive. Therefore we use the
     deepcopy for the horizontal and vertical moves. Finally, as there are no
     truely matching lines all type 3 snakes are only of length 1.
     
     Huub van Dam, SciTech Daresbury Laboratory, June 2006.
     """

  lcs_t = None
  # FP - Forward Pij
  # Records the maximum number of diagonal lines of all candidates paths that
  # passed through node (i,j). P is a dictionary with tuples (i,j) as keys and
  # the maximum number as data.
  FP = { }

  # FV - Forward search path vector
  # Stores the forwards search paths. 
  FV = { }
  # 
  s = search_path_linked()
  s.set_lastpoint(Mb-1,Nb-1)
  FV[(Mb-1,Nb-1)] = s

  # flij - forward last i+j
  # foij - forward old i+j
  # lij is the smallest i+j of the end point of any search path in the current
  # pass.
  # oij is the value of lij of the previous pass. These values will be used
  # to eliminate any entries in P that are no longer nessecary.
  flij = Me+Ne
  foij = Mb+Nb

  finished = (0 != 0)

  # D - the number of non-diagonal steps
  D = 0
  while (D < (Me+Ne) and not(finished)):
    #
    # Work on the forward searches
    #
    D = D + 1
    flij = Me+Ne
    Fkeys = FV.keys()
    Fkeys.sort()
    while (len(Fkeys) > 0):
      key = Fkeys.pop()
      s = FV[key]
      del FV[key]
      (xi,yj) = s.lastpoint
      #
      # now explore the way forward, horizontal first
      #
      keep_horizontal = false
      xih = xi+1
      yjh = yj
      if xih <= Me:
        if FP.has_key((xih,yjh)):
          if FP[(xih,yjh)] < s.lcslength:
            FP[(xih,yjh)] = s.lcslength
            if (yjh == Ne) and (xih == Me):
              #
              # This search path has reached the end point so remove it from
              # the active search paths
              #
              lcs_t = s.snakelist
            else:
              #
              # This search path has not hit the end yet and is still viable
              # so keep it
              #
              keep_horizontal = true
        else:
          FP[(xih,yjh)] = s.lcslength
          if (yjh == Ne) and (xih == Me):
            #
            # This search path has reached the end point so remove it from
            # the active search paths
            #
            lcs_t = s.snakelist
          else:
            #
            # This search path has not hit the end yet and is still viable
            # so keep it
            #
            keep_horizontal = true
        if xih+yjh < flij:
          flij = xih+yjh
      #
      # now explore the vertical direction
      #
      keep_vertical = false
      xiv = xi
      yjv = yj+1
      if yjv <= Ne:
        if FP.has_key((xiv,yjv)):
          if FP[(xiv,yjv)] < s.lcslength:
            FP[(xiv,yjv)] = s.lcslength
            if (yjv == Ne) and (xiv == Me):
              #
              # This search path has reached the end point so remove it from
              # the active search paths
              #
              lcs_t = s.snakelist
            else:
              #
              # This search path has not hit the end yet and is still viable
              # so keep it
              #
              keep_vertical = true
        else:
          FP[(xiv,yjv)] = s.lcslength
          if (yjv == Ne) and (xiv == Me):
            #
            # This search path has reached the end point so remove it from
            # the active search paths
            #
            lcs_t = s.snakelist
          else:
            #
            # This search path has not hit the end yet and is still viable
            # so keep it
            #
            keep_vertical = true
        if xiv+yjv < flij:
          flij = xiv+yjv
      #
      # finally explore the diagonal direction
      #
      keep_diagonal = false
      xid = xi+1
      yjd = yj+1
      if yjd <= Ne and xid <= Me:
        length_diag = s.lcslength + num_same(tol,ref,yjd,dat,xid)
        if FP.has_key((xid,yjd)):
          if FP[(xid,yjd)] < length_diag:
            keep_diagonal = true
            FP[(xid,yjd)] = length_diag
        else:
          keep_diagonal = true
          FP[(xid,yjd)] = length_diag
        if xid+yjd < flij:
          flij = xid+yjd
      if keep_diagonal:
        if keep_horizontal:
          # Keeping at least both the horizontal and diagonal search directions
          # So create a new search path for the horizontal direction
          sa = copy.copy(s)
          sa.set_lastpoint(xih,yjh)
          FV[(xih,yjh)] = sa
        if keep_vertical:
          # Keeping at least both the vertical and diagonal search directions
          # So create a new search path for the vertical direction
          sa = copy.copy(s)
          sa.set_lastpoint(xiv,yjv)
          FV[(xiv,yjv)] = sa
        # Keeping the diagonal search direction as well
        # So update the search path
        s.set_lastpoint(xid,yjd)
        s.add_snake(xid,yjd,xid,yjd,3)
        s.lcslength = length_diag
        if (yjd == Ne) and (xid == Me):
          #
          # This search path has reached the end point so remove it from
          # the active search paths. The diagonal case is special because
          # the snake needed to be added first.
          #
          lcs_t = s.snakelist
        else:
          FV[(xid,yjd)] = s
      else:
        if keep_horizontal:
          if keep_vertical:
            # Keeping both the horizontal and vertical search directions
            # So create a new search path for the vertical direction
            sa = copy.copy(s)
            sa.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = sa
          # Keeping the horizontal direction in any case
          # So update the search path for the horizontal direction
          s.set_lastpoint(xih,yjh)
          FV[(xih,yjh)] = s
        else:
          if keep_vertical:
            # Keeping only the vertical search direction
            # So update the search path for the vertical direction
            s.set_lastpoint(xiv,yjv)
            FV[(xiv,yjv)] = s
          else:
            # Not keeping anything so let the search path die
            pass
    #
    # now tidy up FP
    #
    ij = foij - (Mb+Nb)
    while (ij <= flij - (Mb+Nb)):
      xi = Mb + ij
      yj = Nb
      while (xi >= Mb):
        if FP.has_key((xi,yj)):
          del FP[(xi,yj)]
        xi = xi - 1
        yj = yj + 1
      ij = ij + 1
    foij = flij
    finished = (len(FV) == 0)
  lcs = transform_snake_list(lcs_t)
  return lcs

def filter_lcs(lcs,Nb,Ne,Mb,Me,add,delete):
  """Take the list of snakes and the add and delete dictionaries and remove
     all tolerable adds and/or deletes from the snake list."""
  lcsout = [ ]
  nsnake = len(lcs)
  if nsnake == 0:
    #
    # Check added lines
    #
    ok_add = false
    if Me-Mb+1 <= 0:
      ok_add = true
      xi1 = Mb
      xi2 = Me
    elif add.has_key(Nb-1):
      if add[Nb-1] >= Me-Mb+1:
        ok_add = true
        xi1 = Mb
        xi2 = Me
      else:
        xi1 = Mb
        xi2 = Mb
    else:
      xi1 = Mb
      xi2 = Mb
    #
    # Check deleted lines
    #
    ok_del = true
    iline = Nb
    while (iline <= Ne):
      ok_del = ok_del and delete.has_key(iline)
      iline = iline + 1
    if ok_del:
      yj1 = Nb
      yj2 = Ne
    else:
      yj1 = Nb
      yj2 = Nb
    if ok_add and ok_del:
      if xi2 < xi1 and yj2 < yj1:
        pass
      else:
        lcsout.append((xi1,yj1,xi2,yj2,2))
  else:
    (xi1,yj1,xi2,yj2,ityp) = lcs.pop(0)
    #
    # Check added lines
    #
    ok_add = false
    if xi1-Mb <= 0:
      xi0 = Mb
      ok_add = true
    elif add.has_key(Nb-1):
      if add[Nb-1] >= xi1-Mb:
        xi0 = Mb
        ok_add = true
      else:
        xi0 = xi1-1
    else:
      xi0 = xi1-1
    #
    # Check deleted lines
    #
    iline = Nb
    ok_del = true
    while (iline < yj1):
      ok_del = ok_del and delete.has_key(iline)
      iline = iline + 1
    if ok_del:
      yj0 = Nb
    else:
      yj0 = yj1-1
    if ok_add and ok_del:
      if xi1 == Mb and yj1 == Nb:
        pass
      else:
        lcsout.append((xi0,yj0,xi1-1,yj1-1,2))
    lcsout.append((xi1,yj1,xi2,yj2,ityp))
    xi3 = xi2
    yj3 = yj2
    while (len(lcs) > 0):
      (xi1,yj1,xi2,yj2,ityp) = lcs.pop(0)
      #
      # Check added lines
      #
      ok_add = false
      if xi1 - xi3 - 1 <= 0:
        xi0 = xi3+1
        ok_add = true
      elif add.has_key(yj3):
        if add[yj3] >= xi1 - xi3 - 1:
          xi0 = xi3+1
          ok_add = true
        else:
          xi0 = xi1-1
      else:
        xi0 = xi1-1
      #
      # Check deleted lines
      #
      iline = yj3+1
      ok_del = true
      while (iline < yj1):
        ok_del = ok_del and delete.has_key(iline)
        iline = iline + 1
      if ok_del:
        yj0 = yj3+1
      else:
        yj0 = yj1-1
      if ok_add and ok_del:
        lcsout.append((xi0,yj0,xi1-1,yj1-1,2))
      lcsout.append((xi1,yj1,xi2,yj2,ityp))
      xi3 = xi2
      yj3 = yj2
    #
    # final bit
    #
    # Check added lines
    #
    ok_add = false
    if Me - xi3 <= 0:
      ok_add = true
      xi4 = Me
    elif add.has_key(yj3):
      if add[yj3] >= Me - xi3:
        ok_add = true
        xi4 = Me
      else:
        xi4 = xi3+1
    else:
      xi4 = xi3+1
    #
    # Check added lines
    #
    iline = yj3+1
    ok_del = true
    while (iline <= Ne):
      ok_del = ok_del and delete.has_key(iline)
      iline = iline + 1
    if ok_del:
      yj4 = Ne
    else:
      yj4 = yj3+1
    if ok_add and ok_del:
      if xi3 == Me and yj3 == Ne:
        pass
      else:
        lcsout.append((xi3+1,yj3+1,xi4,yj4,2))
  return lcsout

def print_lcs(fp,lcs):
  """Print the snake list.
     """
  fp.write("=== printing the snake list ===\n")
  i = 0
  nsnake = len(lcs)
  while (i < nsnake):
    (xi1,yj1,xi2,yj2,ityp) = lcs[i]
    fp.write("(%d,%d)->(%d,%d) : %d\n" % (xi1,yj1,xi2,yj2,ityp))
    i = i + 1
  fp.write("=== end of the snake list =====\n")
