c
c <HTML>
c <HEAD><TITLE>Vector Service Program Manual</TITLE></HEAD>
c
c <BODY>
c <PRE>
      subroutine servec(cr,caller)
c
c <H1>The Vector Service Program</H1>
c
c
c <H2>Tabel of Contents</H2></PRE><UL>
c <LI><a href="#notation"><B>Notation Conventions</B></a>
c <LI><a href="#directives"><B>The Directives</B></a>
c <UL>
c <LI><a href="#title">The TITLE Directive</a>
c <LI><a href="#read">The READ Directive</a>
c <LI><a href="#write">The WRITE Directive</a>
c <LI><a href="#init">The INIT Directive</a>
c <LI><a href="#overlap">The S Directive</a>
c <LI><a href="#overlap">The overlap Directive</a>
c <LI><a href="#screen">The SCREEN Directive</a>
c <LI><a href="#print">The PRINT Directive</a>
c <LI><a href="#schmidt">The SCHMIDT Directive</a>
c <LI><a href="#lowdin">The LOWDIN Directive</a>
c <LI><a href="#normalise">The NORMALISE Directive</a>
c <LI><a href="#check">The CHECK Directive</a>
c <LI><a href="#sdiag">The SDIAG Directive</a>
c <LI><a href="#locmp">The LOCMP</a>
c <LI><a href="#moperm">The MOPERM Directive</a>
c <LI><a href="#aoperm">The AOPERM Directive</a>
c <LI><a href="#combine">The COMBINE Directive</a>
c <LI><a href="#extra">The EXTRA Directive</a>
c <LI><a href="#summary">The SUMMARY Directive</a>
c <LI><a href="#transform">The TRANSFORM Directive</a>
c <LI><a href="#occup">The OCCUP Directive</a>
c <LI><a href="#eigen">The EIGEN Directive</a>
c <LI><a href="#crit">The CRIT Directive</a>
c <LI><a href="#ndim">The NDIM Directive</a>
c <LI><a href="#nmdim">The NMDIM Directive</a>
c <LI><a href="#log">The LOG Directive</a>
c <LI><a href="#addvec">The ADDVEC Directive</a>
c <LI><a href="#symmetry">The SYMMETRY Directive</a>
c <LI><a href="#maxdim">The MAXDIM Directive</a>
c <LI><a href="#ATOMMO">The ATOMMO Directive</a>
c <LI><a href="#DUMPFILE">The DUMPFILE Directive</a>
c <LI><a href="#stop">The STOP Directive</a>
c </UL>
c <LI><a href="#error"><B>Error Monitoring</B></a>
c </UL><PRE>
c
c
c
c
c..             atmol vector - service program
c..      -- j.h. v. lenthe  december 1982 / may 1987 / may 1988 ---
c..              moved to GAMESS november 2002
c
c***********************************************************************
c                 i n p u t   -  description
c***********************************************************************
c
c-----------------------------------------------------------------------
c
c<a NAME="notation"><H2>   Notation Conventions</H2></a>
c
c     .. istring   =  string of integer numbers finished by 'end'
c                     'to' convention is allowed i.e. 2 to 5 = 2 3 4 5
c     .. tapnam    =  atmol file 
c                     next parameters are iblock (start-block
c                     for dumpfile) and isect (section-number)
c     .. 'a'/'b'   =  choice between textstrings 'a' and 'b'
c     .. 'flop'    =  text-string  (directive is always text)
c     .. ('aa')    =  optional text-string
c     ..  **       =  after ** comment is given
c     ..  ex.      =  ex. denotes an example (starting in column 1)
c-----------------------------------------------------------------------
c
c ** note  atmol input is free format / maximal 72 columns           **
c **       tapnam etc. must be on same card                          **
c **       istring may go over many cards                            **
c
c ** note  servec works in non symmetry-adapted mode                 **
c **       no basis-set information is (can be) written              **
c **       a getq is possible if section is written to the standard  **
c **       dumpfile, in which case the result may be adapted         **
c **       GAMESS requires a gamess keyword on the 'write' directive **
c **       **HOWEVER** no real checking of vectors is performed      **
c
c-----------------------------------------------------------------------
c
c<a NAME="directives"><H2>   The Directives</H2></a>
c
c *************** each directive is immediately obeyed  ****************
c
c
c
c<a NAME="title"><B>  The TITLE Directive</B></a>
c
c title                 ** give on next line the new title
c       ex.  title
c       ex.  this is the title that will appear on vector-files
c
c<a NAME="read"><B>  The READ Directive</B></a>
c
c read   tapnam  ('notran')   ** read vectors
c                             ** if 'notran' is specified the vectors
c                             ** are not transformed to the original
c                             ** basis / all ctrans (symmetry adaption) information
c                             ** is lost however (**not retained**)
c       ex.   read ed3  1 1
c
c read  'input' ndim nmdim    ** read vectors from input freeform
c       ex    read input 3 2
c read  'free' ndim nmdim     ** read vectors from input using read *
c
c read  'format'    ** read vectors and occupations etc. formatted
c                   ** for exact form of input *see* the print vectors
c                   ** 'format' directive. this option is used to
c                   ** transport vectors as text-file
c
c       ex.   read format
c       ex.   .....................
c       ex.   ............................
c
c
c<a NAME="write"><B>  The WRITE Directive</B></a>
c
c write tapnam  ('gamess') ** write vectors (see read)
c                          ** if 'gamess' is given vectors are gamess-compatible
c
c
c<a NAME="init"><B>  The INIT Directive</B></a>
c
c init   tapnam         ** produce new (empty) dumpfile
c       ex.   init  ed3 1
c
c
c<a NAME="S"><B>  The S Directive</B></a>
c
c s  'print'             ** read overlap integrals 
c   or
c s tapnam (without section) ndims 'print' ** read s from foreign dumpfile
c                        **  if print specified the s-matrix is printed
c       ex.   s  
c       ex.   s ed5 1 122
c
c<a NAME="overlap"><B>  The Overlap Directive</B></a>
c
c  overlap 'format'      ** print the overlap matrix (multiplied by 100) between the orbitals
c
c       ex.  overlap
c       ex.  overlap  (i5,(t8,20f7.2),1x)
c
c
c<a NAME="screen"><B>  The SCREEN Directive</B></a>
c
c  screen 'crit'   ** make all s-matrix elements < crit exactly zero
c                  ** if crit is omitted the current criterium is used
c
c<a NAME="print"><B>  The PRINT Directive</B></a>
c
c  print  'on'/'off'    ** switches print flag
c  print  'vectors'  nv  ('occ'/'eig')  ** print first nv vectors
c                               **  default nv = nmdim
c                               **  if 'occ' => give occupations
c                               **  if 'eig' => give orbital energies
c       ex.   print off
c      2ex.   print vectors 7
c
c  print  'vectors'  nv  'format' format  ** print vectors formatted
c           ** according to format(format). nv=0 means nv = nmdim
c           **  default format =  5f16.9
c           ** this option is used to transport vectors as text-file
c           ** see also *read format* / the output looks as follows
c           *** 'vectors' ndim  nv     format
c           *** ((vc(i,j),i=1,ndim),j=1,nv)  (format(format))
c           *** 'occ'  nn         ** nn = 0 for no occupations
c           *** (occ(i),i=1,nn)   (format(format))
c           *** 'eig'  nn         ** nn = 0 for no eigen-values
c           *** (eig(i),i=1,nn)   (format(format))
c           *** 'end of format-print'
c           ** the input for read format should look the same
c
c      ex.   print vectors 0 format   5e16.10
c
c<a NAME="schmidt"><B>  The SCHMIDT Directive</B></a>
c
c  schmidt  istring     ** schmidt orthogonalize the orbitals specified
c       ex.   schmidt   1 5 7 3 end
c
c  schmidt  'set' ('norm'/'nonorm') istring istring 
c                     ** schmidt first istring onto second
c                     ** (no internal orthogonalisations)
c                     ** if 'norm' is specified the set is normalised
c                     == The default is nonorm (the string may be omitted)
c       ex.   schmidt  set  1 4 end 2 3 end
c
c<a NAME="lowdin"><B>  The LOWDIN Directive</B></a>
c
c  lowdin   istring     ** lowdin orthogonalise the orbitals specified
c       ex.   lowdin    2 3 4 5  end
c
c  **  accuracy of schmidt : crit   / accuracy of lowdin : crit*10
c
c<a NAME="normalise"><B>  The NORMALISE Directive</B></a>
c
c  normalise istring   ** normalise the orbitals specified
c
c<a NAME="check"><B>  The CHECK Directive</B></a>
c
c  check  ('print')    **  check vectors on orthonormality  within crit.
c                      **  if in the overlap-matrix of the vectors an
c                      **  element > crit is found the matrix is printed
c                      **  if 'print' is selected it is printed anyway
c
c<a NAME="sdiag"><B>  The SDIAG Directive</B></a>
c
c  sdiag ('set') ('istring') **  diagonalise s-matrix and leave its eigen-
c                            **  vectors as current vector-set, ordered by
c                            **  decreasing eigen-value / needs s-matrix
c                            **  see  's' and 'crit' directives
c                            **  ** if set is specified the s-matrix
c                            **  for that set of unit vectors is diagp-
c                            **  nalised and the eigenvectors are of
c                            **  the original dimension
c
c<a NAME="locmp"><B>  The LOCMP Directive</B></a>
c
c  locmp iset nset     ** localise within the set of nset orbitals
c                      ** starting at orbital iset.
c                      ** V. Magnasco, A. Perico, JCP 47,971 (1967)
c                      ** DOI:<a href="http://dx.doi.org/10.1063/1.1712065">10.1063/1.1712065</a>
c  istring             ** a string of numbers finished by an end indicates
c                      ** the ao's to localise an mo on. when all mo's
c                      ** to be localised are specified, give 'end'
c  'end'               ** # mo's to be localised  .le. nset
c                      ** an moperm may be useful before applying locmp
c                      ** also the 'to' keyword may be used to specify
c                      ** the orbital numbers, i.e. '1 2 3 4 5',
c                      ** '1 to 5', and 'to 5' are all equivalent
c                      **   *** locmp needs s-integrals and vectors so
c                      **   *** the 's' and a 'read' directive must
c                      **   *** be present
c       ex.   locmp 3  7
c       ex.   1 4 5  end
c       ex.   1 2 3  end
c       ex.   end
c
c<a NAME="moperm"><B>  The MOPERM Directive</B></a>
c
c  moperm   istring     ** permute mo's
c       ex.   moperm   1 to 10 13 12  end
c
c<a NAME="aoperm"><B>  The AOPERM Directive</B></a>
c
c  aoperm   istring     ** permute ao's
c       ex.   aoperm  17 18 end
c
c<a NAME="combine"><B>  The COMBINE Directive</B></a>
c
c  combine  tapnam  tapnam ** identical to atmolscf combine option
c                          ** except that no orthogonalisation is done
c                          ** combine first and second vector - set
c        ** next cards define the origin of the ao's **
c           ifr ito 'a'/'b' ii ** ao's ifr to ito (inclusive) in new set
c           ... ...   ...   .. ** are ao's ii onwards of vectorset a/b
c           ... ...   ...   .. ** go on until all ao's are assigned
c           'end'
c        ** next cards define the origin of the mo's **
c           ifr ito 'a'/'b' ii ** same as for ao's / finish with 'end'
c           'end'
c       ex.   combine ed6 1 1 ed4 200 3
c       ex.   1 5 a 1
c       ex.   6 10 b 1
c       ex.   11 20 a 6
c       ex.   21 26 b 6
c       ex.   end
c       ex.   1  4  a 1
c       ex.   5  7  b 1
c       ex.   8 18  a 5
c       ex.   19 26  b 4
c       ex.   end
c        **  Alternatively a simplified form may be used just specifying quantities
c        **  and on the same line
c       ex.  combine  ed6 1 1 ed4 200 3
c       ex.  a 5 b 5  a 10  b 6  end
c       ex.  a 4 b 3  a 11  b 8  end
c
c<a NAME="extra"><B>  The EXTRA Directive</B></a>
c
c  extra    istring     ** add ao's and corresponding unit-vectors
c                       ** (in existing vectors zeros are inserted)
c                       ** warning  : ndim **must** equal nmdim
c           ** e.g. one performed a dz-calculation on h2 using ao's
c           ** 1,2 and 3,4 . now one adds polarisation functions (p)
c           ** on each h-atom (p's will be 3,4,5 and 8,9,10). the start-
c           ** orbitals for this new calculation are produced (starting
c           ** from the old 4*4 set) by putting :
c       ex.   extra 3 4 5 8 9 10 end
c
c<a NAME="summary"><B>  The SUMMARY Directive</B></a>
c
c  summary              ** print summary of files and parameters
c
c<a NAME="transform"><B>  The TRANSFORM Directive</B></a>
c
c  transform tapnam ('notran')   ** transform nmdim*ndim vectors in core
c                                ** with a nmdim*nmdim set on tapnam
c                                ** notran as in read
c
c<a NAME="occup"><B>  The OCCUP or OCC Directive</B></a>
c
c  occup   nval (occ(i),i=1,nval) ** read occupation numbers (free form)
c      ex.    occup 3  2.0 2.0 1.5
c      ex2.   occ 2  2 1 
c           ** the occupation-numbers will appear on any tape4 or dump-
c           ** file written , following this directive
c
c<a NAME="eigen"><B>  The EIGEN Directive</B></a>
c
c  eigen   nval              ** read nval orbital-energies on following
c          (eig(i),i=1,nval) ** according to format(5e15.5)
c      ex.    eigen 2
c      ex.     -20.5673456     -11.43456
c
c<a NAME="crit"><B>  The CRIT Directive</B></a>
c
c  crit    acr  npower ** set criterium to acr.10**(-npower)
c                      ** (default  1.0e-14)
c      ex.    crit    1.0   20
c
c<a NAME="ndim"><B>  The NDIM Directive</B></a>
c
c  ndim    nd   ** reduce ao-dimension to nd (loose last ao's)
c      ex.    ndim    22
c
c<a NAME="nmdim"><B>  The NMDIM Directive</B></a>
c
c  nmdim   nmd  ** reduce mo-dimension to nmd (loose last mo's)
c
c<a NAME="log"><B>  The LOG Directive</B></a>
c
c  log  'on'/'of' ** switch printflag for printing of directive-cards
c                 ** in interactive sessions log on is more forgiving
c
c<a NAME="addvec"><B>  The ADDVEC Directive</B></a>
c
c  addvec  tapnam ** add vectors (same dimension as current vectors)
c                 ** from specified file ,as last orbitals
c      ex.    addvec ed3 1 2
c  addvec  'unit' n1 n2 .. .. 'end' ** add specified unit vectors
c                                   ** (i.e. with 1.0 in positions n1,n2, ..
c                                   **  respectively ) // 'to' may be used
c      ex.    addvec unit 11 7 1 end
c      ex.    ndim 10 \ addvec unit 1 to 10 end
c      ex.    extra 1 to 10 end
c
c<a NAME="symmetry"><B>  The SYMMETRY Directive</B></a>
c
c  symmetry  nao nirr n1 n2 n3 .. ** 'symmetry' selection of mo's
c           ** nao = number of ao's in each mo to look at (max 8)
c           ** nirr= number of representations to select  (max 12)
c           ** n1,n2,n3 etc : dimensions of subsets within mo's
c            ** next cards define selection (nirr cards) **
c               (iao(i),ich(i),i=1,nao) ** iao = number of ao
c               .. .. .. .. .. .. .. .. ** ich = 1,0 or -1 = sign of ao
c               .. .. .. .. .. .. .. .. **       in mo  (parity is ok)
c       ex.   symmetry 3 2  3 5
c       ex.   1 1    3 0    12  -1
c       ex.   1 0    4 0    10   1
c
c  symmetry 'h'            ** the symmetry is determined from the
c                          ** (1-electron) h-matrix 
c          ** following directives may be issued / finish with 'end'
c      'sets'  is1 is2 ... ** dimensions of subsets in orbital-space
c     'rename' ir1 ir2 ... ** new numbers of the representations
c      'end'               **
c       ex.   symmetry h tape3
c       ex.   order 1 4 2 3
c       ex.   sets  1 4 8 1
c       ex.   end
c      ** this sorts the orbitals for water in a dz basis , keeping
c      ** the 1s and the 1s* , the occupied and the virtual spaces
c      ** apart / the a2 symmetry does not occur in the occupied
c      ** space ,is therefore originally assigned # 4 and is reordered
c      ** to position 2
c
c *******************************************************************
c **** a directive : crit 5 5 : is recommended before symmetry   ****
c **** selection, since orbitals are usually only accurate to    ****
c **** 5 figures                                                 ****
c *******************************************************************
c
c<a NAME="maxdim"><B>  The MAXDIM Directive</B></a>
c
c  maxdim 'low'/'medium'/'high' / set maximum dimension  (default low)
c                              ** low : all actions possible
c                              ** medium : no orthogonalisation possible
c                              ** high  : and no combine possible
c       ex.   maxdim   high
c
c<a NAME="atommo"><B>  The ATOMMO Directive</B></a>
c
c  atommo 'atom' istring 'mo' istring / sorts the mo's listed in the order of the atoms
c
c       ex. atommo atom 1 2 3 end mo 1 2 3 4 5  end
c
c<a NAME="dumpfile"><B>  The DUMPFILE Directive</B></a>
c
c  dumpfile tapnam (without section)  / set the standard dumpfile to tapnam 
c               ** the dumpfile must exist
c               ** default the GAMESS dumpfile is used
c               ** if not available the first dumpfile vectors were read from
c
c       ex. dumpfile ed4 1
c
c<a NAME="stop"><B>  The STOP / FINISH Directive</B></a>
c
c  stop/finish     **  ends calculation
c       ex.   stop
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
      logical lprint,otran_save
      dimension cr(*)
      character*(*) caller
      data ntri/6/
c
      if (lchbas(dum,dum,idum,'check').ne.0) 
     1   call caserr('vb - basis confusion')
c
c...  save ctrans
c
      call save_ctrans
c
      call servec_gamess(cr,caller)
c
      ncore = igmem_max_memory()
      ncore = min(ncore-1000,ntri*maxorb*(maxorb+1)/2)
      kcr =  igmem_alloc(ncore)
c
      mm = ncore/ntri
      maxdim = (dsqrt(1.0d0+8.0d0*mm)-1.0d0)/2.0d0
c     maxdim = min(maxdim,255)
c
      call driver_servec (cr(kcr),cr(kcr+maxdim**2),
     1             cr(kcr+maxdim**2+maxdim*(maxdim+1)/2),
     2             cr(kcr+2*maxdim**2+maxdim*(maxdim+1)/2))
c
      call gmem_free(kcr)
c
c...  restore ctrans block
c
      call rest_ctrans
c
      end
      subroutine servec_gamess(q,caller)
c
c...   communication between servec and gamess
c...   mostly like vb
c...   a bit cumbersome ......
c
      implicit real*8 (a-h,o-z), integer (i-n)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
      common /junkc/ zcom_serv(19),ztitle_serv(10)
      character*8 zcom_serv,ztitle_serv
      dimension q(*)
      character*(*) caller
c
      do i=1,19
         zcom_serv(i) = zcom(i)
      end do
      zcom_serv(4) = 'servec-g'
      do i=1,10
         ztitle_serv(i) = ztitle(i)
      end do
c
c...  calculate 1-electron ints if they are not present 
c...  may be needed  .............
c
      call cgamtovb(numed3,iblk3)
      call ini_set_dump
      call stractlt(ionsec,num)
      nbasisg = num
c
      if (ionsec.le.0) then
         ifl = 0
      else
         call sectst(ionsec,ifl)
      end if
c
      if (ifl.eq.0) then
         if (caller.eq.'start'.and.lpseud.ne.1) then
c...      ispchk set6s opbas,odbas etc .... needed for symmetrizing
            ibasis = ispchk()
            call standv(0,q)
         else if (caller.eq.'server'.or.lpseud.eq.1.or.
     1            caller.eq.'driver') then
            write(iwr,*) 
     1 '**warning possibly no 1-electron integrals in servec**'
            if (lpseud.eq.1) write(iwr,*) 'because of pseudo-potentials'
            if (caller.eq.'server'.or.caller.eq.'driver') 
     1           write(iwr,*) 'because no basis info'
         else
            call caserr('unknown caller to servec')
         end if
      end if
c
      return
      end
      subroutine driver_servec(vc,s,sq1,tri1)
c
c...  vector-service program  driver
c
      implicit real*8 (a-h,o-z)
      dimension vc(*),sq1(*),s(*),tri1(*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
      character*8 ztitle,zcomm
      common /junkc/ zcomm(19),ztitle(10)
c
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
c      for communicating occ
c...  standard blkorbs (from getqvb :-) 
c...  including this everywhere too tiresome
c
      real*8 deig, dpop, etott
      integer nba, new, ncol, jeig, jpop, ipadblko
c
      common/blkorbs/deig(maxorb),dpop(maxorb),
     *               etott,nba,new,ncol,jeig,jpop,ipadblko
c...    was local/blkin ...
c...    /junk/ may be 3*maxorb reals + 16*maxorb integers (cf adapt)
      common /junk/ieig,iocc,eig(maxorb),occ(maxorb),hulp(maxorb),
     1             scrat(maxorb),ip(maxorb),ipb(maxorb),ijun(maxorb,9)
c                  ijun10-ijun12 available
      logical lprint
      logical lprinl
      character*8 mon,moff,mormat,mocc,meig,mset,mectors,minput,
     1     mblan,mmoper,maoper,munit,mh,mrewin,mretur,notran,
     1     itest,ites,mfree,zcomm4,zcomm5,zformat
      character*30 nformat
c
      parameter (nww=35)
      character*4  iww(nww),ytest
c     data iww/'title','read','write','init','s','overlap','print',
c    1         'schmidt','lowdin','moperm','aoperm','combine','extra',
c    2         'summary','occup','occ','eigen','crit','stop','finish',
c    3         'ndim','log','nmdim','addvec','symmetry',
c    4                'maxdim','check', 'locmp','',
c    5         'sdiag','tranform','screen','normalis'/
      data iww/'titl','read','writ','init','s   ','over','prin',
     1         'schm','lowd','mope','aope','comb','extr',
     2         'summ','occu','occ ','eige','crit','stop','fini',
     3         'ndim','log ','nmdi','addv','symm',
     4                'maxd','chec', 'locm',' ',
     5         'sdia','tran','scre','norm','atom','dump'/
c
      data lprinl/.false./,ndims,ndim,nmdim/3*0/
      data crit/1.0d-14/
c 
      data mon,moff,mormat,mocc,meig,mset,mectors,minput,mfree,
     1     mblan,mmoper,maoper,munit,mh,mrewin,mretur,notran
     1 /'on','off','format','occ','eig','set','vectors','input','free',
     1  ' ','moperm','aoperm','unit','h','rewind','return','notran'/
c
      maxdio = maxdim
      lprint = .true.
      ieig = 0
      iocc = 0
c
      write(iwr,2) maxdim
2     format(/,10x,'*******************************************',
     1       /,10x,'gamess vector service module',' maxdim ',i5,
     1       /,10x,'*******************************************'/)
c
3     call reset_dump
1     call input
c
      if (nbasisg.le.0) nbasisg = ndim
      i492 = 492
      call stractlt(i492,nbasisg)
c
      call inpa4(ytest)
      if (lprinl.or.lprint) call outrec
      ipp = locatc(iww,nww,ytest) 
      if (ipp.le.0) then
         if (lprinl) then
            write(iwr,*) ' *** input not recognised / try again ***'
            go to 1
         end if
         call caserr('directive not recognised')
      end if
      go to (10,20,30,40,50,55,60,70,80,90,100,110,120,130,140,140,150,
     1       160,1000,1000,170,180,190,200,210,240,250,270,280,
     1       290,300,310,320,330,340)
     2       ipp
c...  title
10    read(ird,11) ztitle
11    format(10a8)
      if (lprint) write(iwr,12) ztitle
12    format('  current title ',10a8,/)
      go to 1
c...  read
20    call inpa(ites)
      if (ites.ne.minput.and.ites.ne.mfree) go to 21
c...  read input
      call rdvcin(vc,ndim,nmdim,ites)
      go to 1
c...  read format
21    if (ites.ne.mormat) go to 25
      call form_servec(vc,ndim,nmdim,occ,eig,iocc,ieig,-1)
      go to 1
c...  read tapnm
25    call cinput(jrec,jump,-1) 
      call set_dump(isect,.false.)
c...   check for notran keyword
      call inpa(ites)
      ntra = 0
      if (ites.eq.notran) call caserr('notran not active')
      call getqvb(vc,ndim,nmdim,isect,'print')
      if (jpop.eq.1) call dcopy(nmdim,dpop,1,occ,1)
      if (jeig.eq.1) call dcopy(nmdim,deig,1,eig,1)
      iocc = jpop
      ieig = jeig
      go to 3
c...  write
30    call set_dump(isect,.false.)
      call inpa4(ytest)
      if (ytest.eq.'game'.or.nmdim.eq.ndim) then
        zcomm(5) = 'servgg'
        write(iwr,'(a)') ' Vectors are written (unadapted) for GAMESS'
      else
        zcomm(5) = 'servec'
        write(iwr,'(a)') ' ** Vectors are not for GAMESS **'
      endif
      call ini_ctrans
      call putq(zcomm,ztitle,eig,occ,ndim,ndim,nmdim,
     1          ieig,iocc,vc,isect,iddd)
      go to 3
c...  init
40    call set_dump(isect,.true.)
      go to 3
c...  s
50    call set_dump(nb,.false.)
c...  read s integrals from standard or foreign dumpfile 
      call get1e(s,dummy,'s',s)
      ndims = nbasisg
      if (nb.ne.0) ndims = nb
c...   print ?
      call inpa(itest)
      if (itest.eq.'print') then
         write(iwr,51) ndims
51       format(/' S-matrix of dimension',i4,' as read ..')
         call prtrs(s,ndims,iwr)
      end if
      go to 3
c...  overlap
55    nformat = char1(istrt(jrec+1):istrt(jrec+1)+29)
      if (nformat.eq.' ') nformat = '(i5,(t8,20f7.2),1x)'
      if (nmdim.lt.1.or.ndims.le.0) then
         write(iwr,*) 'no vectors or S-matrix'
         go to 3
      end if
      call fmos(tri1,s,vc,hulp,nmdim,ndim,crit)
      do i=1,nmdim*(nmdim+1)/2
         tri1(i) = tri1(i)*100.0d0
      end do
      write(iwr,*) ' ------------ MO Overlap (* 100) -------------'
      write(iwr,*) 'format',nformat
      k = 1
      do i=1,nmdim
         je = k + i - 1
         write(iwr,nformat) i,(tri1(j),j=k,je)
c56       format(i5,(t8,20f7.2),1x)
         k = k + i
      end do
      go to 3
c...  print
60    call inpa(itest)
      if (itest.eq.mon) lprint=.true.
      if (itest.eq.moff) lprint=.false.
      if (itest.ne.mectors) go to 63
      call inpi(nn)
      if (nn.le.0) nn = nmdim
      call inpa(itest)
      if (itest.eq.mormat) go to 61
      if (itest.eq.mocc) then
         call prvc(vc,nn,ndim,occ,'v','l')
      else if (itest.eq.meig) then
         call prvc(vc,nn,ndim,eig,'v','l')
      else 
         call prvc(vc,nn,ndim,eig,'o','l')
      end if
      go to 63
c...  print vectors  nv 'format' format
61    call form_servec(vc,ndim,nn,occ,eig,iocc,ieig,1)
c
63    go to 1
c...  schmidt
70    iort=-1
      call inpa(itest)
      if (itest.eq.mset) go to 75
      if (itest.ne.mblan)  call cinput(jrec,jump,-1) 
      go to 85
c.    set
75    iort = 0
c
      call inpa(itest)
      if (itest.eq.'norm') then
         ites = 'norm'
         write(iwr,'(a)') ' Schmidt set orbitals are Normalised '
      else if (itest.eq.'nonorm') then
         ites = 'nonorm'
         write(iwr,'(a)') ' Schmidt set orbitals are **not** Normalised'
      else
         write(iwr,'(a)') ' Schmidt set orbitals are **not** Normalised'
         ites = 'nonorm'
         call cinput(jrec,jump,-1)
      end if
c
      call rdistr(ip,n1,maxorb)
      call rdistr(ip(n1+1),n2,maxorb-n1)
      if (n2.eq.255) n2 = nmdim - n1
      if ((n1+n2).gt.nmdim) call caserr('# mos in schmidt set > nmdim')
      go to 87
c...  lowdin
80    iort= 1
c
85    call rdistr(ip,nn,maxorb)
      nn = min(nn,nmdim)
87    if (ndims.lt.ndim) call caserr('s-matrix dim wrong')
      call coperm(ip,hulp,nmdim,iwr)
      call pinver(ip,ipb,nmdim)
      call change(ip,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      if (iort.eq.-1) call normvc(vc,s,sq1,ndim,nn,crit)
      if (iort.eq.0)
     1 call schmids(vc,vc(n1*ndim+1),s,hulp,n1,n2,ndim,crit,ites)
      if (iort.eq.1) call lowdins(vc,s,tri1,sq1,hulp,scrat,nn,ndim,crit)
      call change(ipb,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      go to 1
c...  moperm
90    call rdistr(ip,nn,maxorb)
      call coperm(ip,hulp,nmdim,iwr)
      if (lprint) write(iwr,95) mmoper,(ip(i),i=1,nmdim)
      if (lprint) write(iwr,'(1x)')
      call change(ip,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      go to 1
95    format(' ',a6,'utation :',(t15,25i4))
c...  aoperm
100   call rdistr(ip,nn,maxorb)
      call coperm(ip,hulp,ndim,iwr)
      if (lprint) write(iwr,95) maoper,(ip(i),i=1,ndim)
      if (lprint) write(iwr,'(1x)')
      call change(ip,scrat,scrat,vc,ndim,1,hulp,nmdim,ndim)
      go to 1
c...  combine
110   maxsq1 = maxdio**2+maxdio*(maxdio+5)/2 - 4*maxdim
      call combin(vc,sq1,ndim,nmdim,
     1            ijun(1,1),ijun(1,2),ijun(1,3),ijun(1,4),ijun(1,5),
     2            ijun(1,6),ijun(1,7),ijun(1,8),ijun(1,9),maxsq1)
      go to 1
c...  extra
120   call extra_servec(vc,occ,eig,hulp,ndim,nmdim)
      go to 1
c...  summary
130   call secsum
      call sumvecs
      call whtps
      write(iwr,131) lprint,maxdim,ndim,nmdim
131   format(/' lprint ',l5,'  maxdim ',i4,'  ndim ',i4,'  nmdim ',i4,/)
      go to 1
c...  occup  ... occ
140   call vclr(occ,1,maxorb)
      call inpi(n)
      do 141 i=1,n
141   call rinpf(occ(i))
      iocc = 1
      if (lprint) write(iwr,146) (occ(i),i=1,nmdim)
      if (lprint) write(iwr,'(1x)')
      go to 1
155   format(5e15.5)
146   format(/' occupations :',(t18,8f13.6))
147   format(/' eigenvalues :',(t18,8f13.6))
c...  eigen
150   call vclr(eig,1,maxorb)
      call inpi(n)
      read(ird,155) (eig(i),i=1,n)
      ieig = 1
      if (lprint) write(iwr,147) (eig(i),i=1,nmdim)
      if (lprint) write(iwr,'(1x)')
      go to 1
c...  crit
160   call inpf(acrit)
      call inpi(icrit)
      crit = acrit*10.0d0**(-icrit)
      if (lprint) write(iwr,161) crit
161   format('  criterion set to ',e12.5,/)
      go to 1
c...  ndim
170   call inpi(ndimi)
      ko = ndim + 1
      kn = ndimi + 1
      do 175 i=2,nmdim
         call fmove(vc(ko),vc(kn),ndimi)
         ko = ko + ndim
         kn = kn + ndimi
175   continue
      ndim = ndimi
      go to 1
c...  log
180   call inpa(itest)
      if (itest.eq.mon) lprinl = .true.
      if (itest.eq.moff) lprinl = .false.
      go to 1
c...  nmdim
190   call inpi(nmdim)
      if (nmdim.le.0) ndim = 0
      go to 1
c...  addvec
200   call inpa(itest)
      if (itest.eq.munit) go to 205
c...  addvec from dumpfile or tape4
      call cinput(jrec,jump,-1)
      call set_dump(isect,.false.)
      mtype = 0
      call secget(isect,mtype,k)
      if (mtype.eq.33) then
         call getqvb(vc(nmdim*ndim+1),nd2,nmd2,isect,tag)
      else
         call getq(vc(nmdim*ndim+1),eig(nmdim+1),
     1             occ(nmdim+1),nd2,nmd2,3,ii,ii,isect)
      end if
      if (nd2.ne.ndim) call caserr('dim. of extra set <> resident set')
      nmdim = nmdim + nmd2
      if (nmdim.gt.maxorb) call caserr('too many vectors')
      go to 3
c...  addvec  / unit vectors
205   call rdistr(ip,nn,maxorb)
      do i=1,nn
         if (ip(i).eq.0) exit
         call vclr(hulp,1,ndim)
         if (ip(i).gt.ndim) 
     1     call caserr('ndim not set big enough - addvec')
         hulp((ip(i))) = 1.0d0
         call fmove(hulp,vc(nmdim*ndim+1),ndim)
         nmdim = nmdim + 1
      end do
c
206   if (lprint) write(iwr,209) ndim,nmdim
209   format(' addvec results in ndim=',i4,' / nmdim=',i4,/)
      go to 1
c...  symmetry
210   call inpa(ites)
      if (ites.ne.mh) go to 215
c...  symmetry check on h-matrix (get from standard dumpfile)
      call get1e(tri1,dummy,'h',tri1)
      call symh_servec(tri1,sq1,vc,occ,eig,hulp,nmdim,ndim,crit,lprint)
      go to 1
c...  symmetry check on vector-coefficients
215   call cinput(jrec,jump,-1)
      call symmetr(vc,eig,occ,ip,hulp,ipb,ndim,nmdim,crit)
      go to 1
c...  maxdim
240   call smaxd(maxdim,ncore,maxdio)
      go to 1
c...  check
250   if (ndims.lt.ndim) call caserr('s-matrix dimension error')
      call checks(vc,s,tri1,hulp,ndim,nmdim,crit)
      go to 1
c...  locmp
270   if (ndims.lt.ndim) call caserr('s-matrix dimension error')
      call inpi(ill)
      call inpi(nll)
      write(iwr,271) ill,nll
271   format(/' // localise in set starting at mo ',i3, ',of dim. ',i3)
      if (ill+nll-1.gt.nmdim) call caserr('# vectors too small /local')
      call locmp(vc(ndim*(ill-1)+1),s,hulp,sq1,nll,ndim,50,crit)
      go to 1
c...  
280   go to 1
c...  sdiag 'set' 'istring'
290   ndim = ndims
      if (ndim.le.0) call caserr('no or too small s-matrix')
      call inpa(ites)
      if (ites.eq.mset) then
         call rdistr(ip,nn,maxorb)
         write(iwr,293) nn
293      format(/' SDIAG SET dimension ',i5)
         write(iwr,95) mset,(ip(i),i=1,nn)
         kk = 0
         do 194 i=1,nn
            do 194 j=1,i
               kk = kk + 1
               imi = max(ip(i),ip(j))
               jmj = min(ip(i),ip(j))
               ks = imi*(imi-1)/2 + jmj
               tri1(kk) = s(ks)
194      continue
      else
         nn = ndim
         call fmove(s,tri1,ndim*(ndim+1)/2)
      end if
      call ibasgn(nn+1,0,nn,ipb)
      call jacobi(tri1,iky,nn,vc,ipb,nn,eig,2,3,crit)
      ieig = 1
      iocc = 0
      nmdim = ndim
      write(iwr,291) (eig(i),i=1,nn)
291   format(//' eigenvalues of s-matrix ',/,(7e11.3))
      write(iwr,292)
292   format(/'  eigenvectors of s-matrix are the current vectors ')
c
c...  make nn * ndim vectors if set was specified
c
      if (ites.eq.mset) then
         nmdim = nn
         call fmove(vc,sq1,nn*nn)
         call vclr(vc,1,ndim*ndim)
         do 195 i=1,nmdim
            do 195 j=1,nmdim
            vc(ip(j)+(i-1)*ndim)=sq1(j+(i-1)*nmdim)
195      continue
      end if
      go to 1
c...  transform
300   call set_dump(isect,.false.)
c...   check for notran keyword
      call inpa(ites)
      ntra = 0
      if (ites.eq.notran) call caserr('notran not implemented') 
      call getq(sq1,eig,occ,nd2,nmd2,ntra,ieig,iocc,isect)
      if (nd2.ne.nmd2) call caserr('transf. matrix must be square')
      if (nd2.ne.nmdim) call caserr('dim. transform = nmdim')
      call matmpl(vc,sq1,tri1,nmdim,ndim)
      ndims = 0
      go to 3
c...  screen crit
310   if (ndims.eq.0) call caserr('s-matrix must be present for screen')
      call inpf(acrit)
      call inpi(icrit)
      if (acrit.ne.0.0d0) then
         critsc = acrit*10.0d0**(-icrit)
      else
         critsc = crit
      end if
      write(iwr,311) critsc
311   format(/'  S-matrix screened for elements < ',1pe12.5)
      smax = 0.0d0 
      nmax = 0
      do 312 i=1,ndims*(ndims+1)/2
         if (dabs(s(i)).lt.critsc) then
            if (dabs(s(i)).gt.dabs(smax)) smax = s(i)
            if (s(i).ne.0.0d0) nmax = nmax + 1
            s(i) = 0.0d0
         end if
312   continue
      write(iwr,313) nmax,smax
313   format('  # elements zeroed',i6,' max elemt zeroed ',1pe12.5)
      go to 1
c...  normalise
320   if (ndims.lt.ndim) call caserr('s-matrix needed for normalise')
      call rdistr(ip,nn,maxorb)
      nn = min(nn,nmdim)
      call coperm(ip,hulp,nmdim,iwr)
      call pinver(ip,ipb,nmdim)
      call change(ip,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      do 321 i=1,nn
321   call normvc(vc(ndim*(i-1)+1),s,sq1,ndim,1,crit)
      call change(ipb,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      go to 1
c. ..  atommo
330   if (ndim.le.0.or.nmdim.le.0) call caserr('atommo needs orbitals')
      call inpa4(ytest)
      if (ytest.ne.'atom') call caserr('atommo atom first')
      call rdistr(ijun(1,1),natom,maxorb)
      call inpa4(ytest)
      if (ytest.ne.'mo') call caserr('atommo mos next')
      call rdistr(ip,nmo,maxorb)
         call coperm(ip,hulp,nmdim,iwr)
         call pinver(ip,ipb,nmdim)
         call change(ip,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      call atommo(vc,ndim,nmo,ijun(1,1),natom,ijun(1,3),hulp,
     1            ijun(1,4),ijun(1,2))
      call change(ijun(1,2),occ,eig,vc,1,ndim,hulp,ndim,nmo)
         call change(ipb,occ,eig,vc,1,ndim,hulp,ndim,nmdim)
      go to 1
c...  dumpfile (set master dumpfile)
340   call set_dump(isec,.false.)
      call ini_set_dump
      go to 1
c...  stop
1000  write(iwr,999)
999   format(/,10x,'*******************************************',
     1       /,10x,'   end of gamess vector service module',
     1       /,10x,'*******************************************'/)
c
      return
      end
      subroutine set_dump(isect,ocreate)
c
c...   set the dumpfile as specified in the input (return the section)
c...   if ocreate=true initialse the dumpfile
c
      implicit none
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      logical ocreate, oint
      character*4 ytest
      character*8 ztest
      integer iblkdu_real,numdu_real
      integer iblk,iunit,irep,isect,jrec,jump
      save iblkdu_real,numdu_real
      data iblkdu_real,numdu_real/2*0/
c
      call inpa4(ytest)
      iblk  = 1
      call filec(ytest,iblk,iunit,irep)
      if (irep.ne.0) then
c
c       ytest is not a file name, so it is either:
c       - another keyword
c       - a section number
c
        call intval(ytest,isect,oint)
        if (.not.oint) then
c
c         ytest is not a section number so it must be a keyword
c
          irep  = 0
          isect = 0
          iblk  = iblkdu_real
          iunit = numdu_real
          numdu = numdu_real
          call cinput(jrec,jump,-1) 
        endif
        return
      else
c
c       ytest is a file name so get the number of the starting block
c
        call inpi(iblk)
        call filec(ytest,iblk,iunit,irep)
      endif
      if (numdu.gt.0) call revind
      if (ocreate) then
         call vclr(apos,1,508)
         iblkdu = iblk
         numdu = iunit
         orevis(1) = .false.
         if (numdu.gt.0) call revind
      end if
      if (numdu_real.le.0) then
         numdu_real = iunit
         iblkdu_real = iblk
         write(iwr,1) yed(numdu_real),iblkdu_real
      end if
      call secini(iblk,iunit)
      call inpa(ztest)
      call intval(ztest,isect,oint)
      if (.not.oint) then
        isect = 0
        call cinput(jrec,jump,-1) 
      endif
      return
c
      entry reset_dump
      if (numdu.gt.0) call revind
      if (numdu_real.gt.0) call secini(iblkdu_real,numdu_real)
      return
c
      entry ini_set_dump
      iblkdu_real = iblkdu
      numdu_real = numdu
      if (numdu_real.gt.0) then
         write(iwr,1) yed(numdu_real),iblkdu_real
1        format(1x,/,1x,' **** master dumpfile set to ',a4,
     1          ' block',i6,' ****')
      else
         write(iwr,2)
2        format(1x,/,1x,' **** master dumpfile not yet set ****')
      end if
c
      return
      end
      subroutine smaxd(maxdim,ncore,mxdio)
c
c...  set new maxdim
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*4 name(4),itest
      data name/' ','low','medi','high'/
c
      call inpa4(itest)
c
      do 10 i=1,4
      if (itest.eq.name(i)) go to 20
10    continue
      call caserr('maxdim specification not recognised')
c
20    go to (1,1,3,4),i
1     maxdim = mxdio
      write(iwr,61) maxdim
61    format(' maxdim restored to old value ',i3,/)
      return
3     maxdim = mxdio + dsqrt(dble(mxdio))
      maxdim = min(maxdim,255)
      write(iwr,63) maxdim
63    format(' maxdim set to',i4,' no orthog possible',/)
      return
4     maxdim = dsqrt(dble(ncore))
      maxdim = min(maxdim,255)
      write(iwr,64) maxdim
64    format(' maxdim set to',i5,' no orthog or combine or trans',/)
      return
      end
      subroutine extra_servec(vc,occ,eig,hulp,ndim,nmdim)
      implicit real*8 (a-h,o-z)
      dimension vc(*),occ(*),eig(*),hulp(*)
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
      logical lprint
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/junk/ijunk(2),rjunk(4*maxorb),jug(maxorb),ind1(maxorb)
c
      if (ndim.ne.nmdim) call caserr('nmdim ne ndim at start of extra')
      call rdistr(jug,mextra,maxorb)
2     if(mextra)9,8,9
9     if (lprint) write(iwr,10)
10    format(/' list of extra basis functions')
      if (lprint) write(iwr,11)(jug(i),i=1,mextra)
11    format(30i4)
      if (lprint) write(iwr,11)
      goto 12
8     write(iwr,13)
13    format(/' no extra basis functions specified',/)
      return
c...  add unit matrix and permute ao's
12    call extr(vc,vc,nmdim,nmdim+mextra)
      call izero(ind1,nmdim+mextra,1)
      call vclr(occ(nmdim+1),1,mextra)
      call vclr(eig(nmdim+1),1,mextra)
      call icopy(mextra,jug,1,ind1(nmdim+1),1)
      nmdim = nmdim + mextra
      ndim = nmdim
      call coperm(ind1,jug,ndim,iwr)
      call pinver(ind1,jug,ndim)
      call change(jug,hulp,hulp,vc,ndim,1,ind1,nmdim,ndim)
c
      return
      end
      subroutine combin(q,scr,nbasis,newbas,
     1                  lga,lgb,lgab,ilifa,ilifb,iga,igb,igab,itab,
     2                  maxsq1)
      implicit real*8 (a-h,o-z)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/blkorbs/deig(maxorb),dpop(maxorb),
     1             dumc,idumc(3),jeig,jpop
      common /junk/ieig,iocc,eig(maxorb),occ(maxorb)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical iftran,nopr,mugg,mug,iprin,jchek,ipunch
      logical lga,lgb,lgab
      dimension iab(2),q(*),scr(*)
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
      logical lprint
      dimension lga(*),lgb(*),lgab(*),ilifa(*),
     *          ilifb(*),iga(*),igb(*),igab(*),itab(*)
      character*4 nam,ytest,itest,iab,iend,itab,iblank
      logical old
      data iab/'a','b'/
      data iend/'end'/,iblank/' '/
c
      nwb = 0
      newa = 0
      iocc = 1
      ieig = 1
c... read in a-set vectors
555     call set_dump(isec,.false.)
        call getqvb(scr(nwb+1),nbb,newb,isec,'print')
        iocc = min(iocc,jpop)
        ieig = min(ieig,jeig)
        if (iocc.eq.1) call dcopy(newb,dpop,1,occ(newa+1),1)
        if (ieig.eq.1) call dcopy(newb,deig,1,eig(newa+1),1)
        call reset_dump
        if (nbb.gt.maxdim) call caserr('dimensions wrong in combine') 
        nw=nbb*newb
        if (nwb.eq.0) then
          newa=newb
          nba=nbb
          nwb=nbb*newb
          go to 555
        end if
c... read in b-set vectors
      newbas = newa+newb
      nbasis = nba + nbb
      if ((nwb+nw).gt.maxsq1) call caserr('dimension > maxdim /combine')
      if (nbasis.gt.maxorb.or.newbas.gt.maxorb) 
     1    call caserr('maxorb exceeded in combin')
c
      if (iocc.eq.1) call dcopy(newbas,occ,1,dpop,1)
      if (ieig.eq.1) call dcopy(newbas,eig,1,deig,1)
c
      old = .true.
      n = 0
      ka = 0
      kb = 0
      la=0
      do 4 ia=1,newbas
      ilifa(ia)=(ia-1)*nba
      ilifb(ia)=(ia-1)*nbb+nwb
      lga(ia)=.false.
      lgb(ia)=.false.
4     lgab(ia)=.false.
c... process lcbf params
888   if (old) callinput
      callinpa4(itest)
      if (itest.eq.iblank) then
         call input
         call inpa4(itest)
      end if
      if(itest.eq.iend)goto999
      if (itest.eq.iab(1).or.itest.eq.iab(2)) then
         old = .false.
         call inpi(nn)
         m = n + 1
         n = n + nn
         if (itest.eq.iab(1)) then
            k = ka + 1
            ka = ka + nn
         else
            k = kb + 1
            kb = kb + nn
         end if
      else
         if (.not.old) call caserr('do not mix old and new')
         call cinput(jrec,jump,-1) 
         callinpi(m)
         callinpi(n)
         callinpa4(itest)
         callinpi(k)
      end if
      if(m.lt.1.or.n.gt.newbas.or.k.lt.1.or.m.gt.n)
     *call caserr('ao-specification wrong')
      if(itest.eq.iab(2))goto8
      if(itest.ne.iab(1)) call caserr('ao-specification wrong')
      do7ia=m,n
      if(k.gt.newa.or.lga(k).or.lgab(ia)) call caserr('ao-spec wrong')
      lga(k)=.true.
      lgab(ia)=.true.
      igab(ia)=k
      itab(ia)=iab(1)
      iga(k)=ia
      k=k+1
7     la=la+1
      goto888
8     do9 ia=m,n
      if(k.gt.newb.or.lgb(k).or.lgab(ia)) call caserr('ao-spec wrong')
      lgb(k)=.true.
      lgab(ia)=.true.
      igab(ia)=k
      itab(ia)=iab(2)
      igb(k)=ia
      k=k+1
9     la=la+1
      goto888
999   if(la.ne.newbas) call caserr('ao-spec wrong in combine')
      if (.not.lprint) go to 43
      write(iwr,300)
300   format(/' origins of lcbf of ab system'/,1x)
c
      nbatch = (newbas-1)/24 + 1
      do 10 kk=1,nbatch
         ke = min(kk*24,newbas)
         kb = (kk-1)*24 + 1
         write(iwr,301) (k,k=kb,ke)
         write(iwr,302) (itab(k),k=kb,ke)
         write(iwr,303) (igab(k),k=kb,ke)
10    continue
301   format(/'  ab  lcbf ',24i5)
302   format('  subsystem',24(4x,a1))
303   format('    lcbf   ',24i5)
c
c... construct trial mos
43    nbsq = nbasis**2
      la=0
      do 44 ia=1,newbas
      lga(ia)=.false.
      lgb(ia)=.false.
44    lgab(ia)=.false.
c... process mo specification lines
      call vclr(q,1,nbsq)
      call vclr(occ,1,newbas)
      call vclr(eig,1,newbas)
      n = 0
      ka = 0
      kb = 0
8888  if (old) callinput
      callinpa4(itest)
      if (itest.eq.iblank) then
         call input
         call inpa4(itest)
      end if
      if(itest.eq.iend)goto9999
      if (itest.eq.iab(1).or.itest.eq.iab(2)) then
         if (old) call caserr(' do not mix old and new ') 
         call inpi(nn)
         m = n + 1
         n = n + nn
         if (itest.eq.iab(1)) then
            k = ka + 1
            ka = ka + nn
         else
            k = kb + 1
            kb = kb + nn
         end if
      else
         if (.not.old) call caserr('do not mix old and new')
         call cinput(jrec,jump,-1)  
         callinpi(m)
         callinpi(n)
         callinpa4(itest)
         callinpi(k)
      end if
      if(m.lt.1.or.n.gt.newbas.or.k.lt.1.or.m.gt.n)
     * call caserr('mo-spec wrong in combine')
      if(itest.eq.iab(2))goto88
      if(itest.ne.iab(1))call caserr('mo-spec wrong in combine')
      do77 ia=m,n
      if(k.gt.newa.or.lga(k).or.lgab(ia))call caserr('mo-spec wrong')
      lga(k)=.true.
      lgab(ia)=.true.
      igab(ia)=k
      itab(ia)=iab(1)
      la=la+1
      ja=ilifa(k)
      nn= (ia-1)*nbasis
      if (iocc.gt.0) occ(ia) = dpop(k)
      if (ieig.gt.0) eig(ia) = deig(k)
      k=k+1
      do77mm=1,newa
77    q(nn+iga(mm))=scr(ja+mm)
      goto8888
88    do99ia=m,n
      if(k.gt.newb.or.lgb(k).or.lgab(ia))call caserr('mo-spec wrong')
      lgb(k)=.true.
      lgab(ia)=.true.
      igab(ia)=k
      itab(ia)=iab(2)
      la=la+1
      ja=ilifb(k)
      nn= (ia-1)*nbasis
      if (iocc.gt.0) occ(ia) = dpop(k+newa)
      if (ieig.gt.0) eig(ia) = deig(k+newa)
      k=k+1
      do99mm=1,newb
99    q(nn+igb(mm))=scr(ja+mm)
      goto8888
9999  newbas = la
cjvl  if(la.ne.newbas)call caserr('mo-spec wrong')
      if (.not.lprint) return
      write(iwr,400)
400   format(/' origins of trial mos for ab system'/,1x)
c
      nbatch = (newbas-1)/24 + 1
      do 110 kk=1,nbatch
         ke = min(kk*24,newbas)
         kb = (kk-1)*24 + 1
         write(iwr,304) (k,k=kb,ke)
         write(iwr,302) (itab(k),k=kb,ke)
         write(iwr,305) (igab(k),k=kb,ke)
110   continue
304   format(/'  ab  mo   ',24i5)
305   format('     mo    ',24i5)
      write(iwr,306)
306   format(1x)
c
      return
      end
      subroutine rdvcin(vc,ndim,nmdim,ites)
c
c...  read vectors from input
c
      implicit real*8 (a-h,o-z)
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
      dimension vc(*)
      character*8 ites,mfree
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      data mfree/'free'/
c
      call inpi(ndim)
      call inpi(nmdim)
      if (ndim.gt.maxdim.or.nmdim.gt.maxdim)call caserr('core overflow')
      nn = ndim*nmdim
c
      if (ites.eq.mfree) then
         read (ird,*) (vc(i),i=1,nn)
      else
         do 10 i=1,nn
10       call rinpf(vc(i))
      end if
c
      return
      end
      subroutine change(idm,dm,em,vc,nspc,nspr,vcin,ndim,nmdim)
      implicit real*8 (a-h,o-z)
      dimension dm(*),em(*),idm(*),vc(*),vcin(*)
c
c     perform permutations of vectors and occupation-numbers/ eigenval
c     according to idm   / nspc,nspr  : column,row-spacing
c
      do 100 k=1,nmdim
      if (idm(k).eq.k) go to 100
c
      l = k
      ll = (l-1)*nspr  +  1
c
c     save vector k
c
      do 50 i=1,ndim
   50 vcin(i) = vc(ll+(i-1)*nspc)
      din = dm(l)
      ein = em(l)
c
      n = idm(l)
   60 nn = (n-1)*nspr  +  1
      do 70 i=1,ndim
   70 vc(ll+(i-1)*nspc) = vc(nn+(i-1)*nspc)
      dm(l) = dm(n)
      em(l) = em(n)
c
      idm(l) = l
      l = n
      ll = nn
      n = idm(l)
      if (idm(n).ne.n) go to 60
c
      do 80 i=1,ndim
   80 vc(ll+(i-1)*nspc) = vcin(i)
      dm(l) = din
      em(l) = ein
      idm(l) = l
c
  100 continue
      return
      end
      subroutine symmetr(vc,eig,occ,ip,hulp,ndsub,ndim,nmdim,crit)
c
c...  select mo's by 'symmetry'-relations
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension vc(ndim,nmdim),eig(nmdim),occ(nmdim),ip(ndim),hulp(ndim)
      dimension ipos(8,12),ich(8,12),nsbs(12),ndsub(nmdim)
c
      common/cntrl_servec/ lprint,maxdim,ncore,nbasisg
      logical lprint
c
      call inpi(nao)
      call inpi(nirr)
      if (nao.gt.8.or.nao.le.0.or.nirr.gt.12.or.nirr.le.0)
     1   call caserr('# aos >8 or <1 or # symm > 12 or < 1')
c
      nsub = 0
      nn = 0
10    nsub = nsub + 1
      call inpi(ndsub(nsub))
      if (ndsub(nsub).le.0) ndsub(nsub) = nmdim - nn
      nn = nn + ndsub(nsub)
      if (nn.lt.nmdim) go to 10
      if (nn.gt.nmdim) call caserr('subset wrong in symmetry')
c
      do 30 irr=1,nirr
         call input
         do 20 iao=1,nao
            call inpi(ipos(iao,irr))
            if (ipos(iao,irr).lt.0.or.ipos(iao,irr).gt.nmdim)
     1         call caserr('ao spec wrong in symmetry')
            call inpi(ich(iao,irr))
20       continue
30    continue
c
      if (.not.lprint) go to 50
      write(iwr,602) nao,nirr,nsub
      do 40 irr=1,nirr
40    write(iwr,603) irr,(ipos(iao,irr),ich(iao,irr),iao=1,nao)
602   format(/'  symmetry selection on',i3,' ao s for',i3,
     1       ' representations in',i4,' subsets',//,
     2       '   representation     ao;ich  etc ....')
603   format(1x,i7,4x,'/',8(i5,i5,3x,'/'))
50    continue
c
c...  loop over subsets
c
      kkv = 0
      do 100 is=1,nsub
         nd = ndsub(is)
         call izero(ip,nd,1)
         kk = 0
         do 80 irr=1,nirr
            do 60 kv=1,nd
               if (icomps(vc(1,kkv+kv),ipos(1,irr),ich(1,irr),nao,crit)
     1             .eq.0) go to 60
               kk = kk + 1
               ip(kk) = kv
60          continue
80       nsbs(irr) = kk
         call coperm(ip,hulp,nd,iwr)
c         if (lprint) write(iwr,604) is,(ip(i)+kkv,i=1,nd)
c604      format('  permutation for subset ',i4,/,(5x,16i4))
         call change(ip,occ(kkv+1),eig(kkv+1),vc(1,kkv+1),1,ndim,
     1               hulp,ndim,nd)
         do 70 i=2,nirr
70       nsbs(nirr-i+2) = nsbs(nirr-i+2) - nsbs(nirr-i+1)
         if (lprint) write(iwr,600) is,nd,(nsbs(i),i=1,nirr)
600      format('  subset ',i3,' (dim ',i3,' ) repr.dim : ',12i4)
100   kkv = kkv + nd
c
      write(iwr,601)
601   format(1x)
c
      return
      end
      integer function icomps(v,ipos,ich,nao,crit)
c
c...  see if v is in representation described by ipso/ich
c
      implicit real*8 (a-h,o-z)
      dimension v(*),ipos(*),ich(*)
c
      icomps = 0
      ipar = 0
c
      do 100 i=1,nao
         if (ich(i)) 10,20,30
c...  -1
10       if (ipar.eq.0) go to 15
         if (ipar*v(ipos(i)).lt.-crit) go to 100
         return
15       ipar = 1
         if (v(ipos(i)).gt.crit) ipar = -1
         go to 10
c...   0
20       if (dabs(v(ipos(i))).lt.crit) go to 100
         return
c...  +1
30       if (ipar.eq.0) go to 35
         if (ipar*v(ipos(i)).gt.crit) go to 100
         return
35       ipar = 1
         if (v(ipos(i)).lt.-crit) ipar = -1
         go to 30
100   continue
c
      icomps = 1
c
      return
      end
      subroutine checks(vc,s,smo,hulp,ndim,nmdim,crit)
c
c
c...  check orthogonality of vectors vc and print s-matrix (mo-basis)
c...  if not a unit-matrix within threshold crit
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension vc(*),s(*),smo(*),hulp(*)
c
      character*4 itest,mprin
      data mprin/'prin'/
c
      call fmos(smo,s,vc,hulp,nmdim,ndim,crit)
c
c     prepare diagonal  (substract  1.0 ) and check
c
      k = 0
      kkd = 0
      ddev = 0.0d0
      do 10 i=1,nmdim
         k = k + i
         smo(k) = smo(k) - 1.0d0
         if (dabs(smo(k)).lt.crit.or.dabs(smo(k)).le.ddev) go to 10
            ddev = dabs(smo(k))
            kkd = i
10    continue
c
c     check off-diagonal
c
      kkoffx = 0
      kkoffy = 0
      offdev = 0.0d0
      if (nmdim.le.1) go to 100
      k = 0
      do 20 i=2,nmdim
         k = k + 1
         je = i - 1
         do 20 j=1,je
            k = k + 1
            if (dabs(smo(k)).lt.crit.or.dabs(smo(k)).le.offdev) go to 20
            offdev = dabs(smo(k))
            kkoffx = i
            kkoffy = j
20    continue
c
100   continue
c...  print non-orth-warnings
      if (kkd.ne.0) write(iwr,1) kkd,ddev
1     format(' **********    check ********** ',/,
     1       ' max. deviation from normality at',i4,' of ',e12.5)
      if (kkoffx.ne.0) write(iwr,2) kkoffx,kkoffy,offdev
2     format(' **********    check ********** ',/,
     1       ' max. deviation from orthogonality at',i4,' /',i4,
     2       ' of ',e12.5)
c
c     print if there's a 'print' on the directive
c
      call inpa4(itest)
      if (itest.ne.mprin) return
c
      write(iwr,3) nmdim,ndim
3     format(/' the (S-I) matrix for ',i4,' vectors (dim',i4,') is ')
c
      call prtrs(smo,nmdim,iwr)
c
      return
      end
      subroutine prtrs(s,n,iwr)
c
c     triangle print for s-matrices
c
      implicit real*8 (a-h,o-z)
      dimension s(*)
      k = 1
      do 10 i=1,n
         je = k + i - 1
         write(iwr,1) i,(s(j),j=k,je)
10    k = k + i
c
1     format(i5,(t6,10e12.5))
c
      return
      end
      subroutine symh_servec(
     1           ihao,hmo,vc,occ,eig,isy,nmdim,ndim,crit,lprint)
c
c     symmetry division using h-matrix blocking
c     ihao   :  h-matrix (ao-basis) + scratch space (ndim*(ndim+1)/2)
c     hmo    :  scratch       (nmdim*(nmdim+1)/2)
c     vc     :  vector-set
c     occ/eig   occupations + eigenvalues (must be permuted as well)
c     isy    :  scratch vector to store symmetries  (ndim)
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension vc(ndim,nmdim),ihao(*),hmo(*),isy(ndim),occ(*),eig(*)
      dimension nd(20),ior(20)
      logical lprint
c 
      character*8 mend,mrenam,msets,ites
      data mend,mrenam,msets /'end','rename','sets'/
c...  transform h-matrix to mo-basis
c
      call fmos(hmo,ihao,vc,isy,nmdim,ndim,crit)
c
c...  sort symmetries out
c
      call symdiv(hmo,nmdim,nsym,isy,crit)
c
      if (lprint) write(iwr,1) (isy(i),i=1,nmdim)
1     format(/' symmetry selection on h-matrix blocking  / ',
     1 ' global assignment :',/,(20i3))
      if (lprint) write(iwr,7)
      nset = 1
      nd(1) = nmdim
c
c...  try to read more instructions
c
5     call input
      call inpa(ites)
c
      if (ites.eq.mend) go to 100
      if (ites.eq.mrenam) go to 10
      if (ites.eq.msets) go to 20
      call cinput(jrec,jump,-1) 
      go to 100
c
c...  rename
c
10    nn = 0
11    nn = nn + 1
      call inpi(ior(nn))
      if (ior(nn).ne.0) go to 11
      if (nn.ne.nsym+1) call caserr('wrong # symmetries')
c...  change names of representations
      do 12 i=1,nmdim
12    isy(i) = ior(isy(i))
c
      if (lprint) write(iwr,2) (isy(i),i=1,nmdim)
2     format(' new global assignment :',/,(20i3))
      if (lprint) write(iwr,7)
c
      go to 5
c
c...  sets
c
20    nset = 0
      nn = 0
21    nset = nset + 1
      if (nset.gt.20) call caserr('too many sets')
      call inpi(nd(nset))
      nn = nn + nd(nset)
      if (nd(nset).gt.0) go to 21
      nd(nset) = nmdim - nn
      if (nd(nset).eq.0) nset = nset - 1
      go to 5
c
c...  end or not recognized word
c
100   if (lprint) write(iwr,3) nset,(nd(i),i=1,nset)
3     format(' division of orbital space into',i4,' sets / dimensions ',
     1       10i4)
      kk = 0
      ks = 0
      do 39 i=1,nset
         nnd = nd(i)
         do 35 j=1,nsym
            ior(j) = 0
            do 33 k=1,nnd
               if (isy(k+kk).ne.j) go to 33
               ior(j) = ior(j) + 1
               ks = ks + 1
               ihao(ks) = kk + k
33          continue
35       continue
         if (lprint) write(iwr,6) i,(ior(j),j=1,nsym)
6     format(' subset',i3,' dim. irrep. repr. :',20i4)
39    kk = kk + nnd
c
      if (lprint) write(iwr,4) (ihao(i),i=1,nmdim)
4     format(1x,/,' total symmetry permutation :',/,(5x,30i4))
      if (lprint) write(iwr,7)
7     format(1x)
c
      call change(ihao,occ,eig,vc,1,ndim,hmo,ndim,nmdim)
c
      return
      end
      subroutine symdiv(hmat,ndim,nrep,isy,critsym)
c
c...  find symmetry classification from 1-electron matrix
c...  results (symmetry-numbers) are returned in isy
c
      implicit real*8 (a-h,o-z)
      dimension hmat(*),isy(ndim)
c
      logical ch
      ch(i,j) = dabs(hmat(j*(j-1)/2+i)).le.critsym
c
      call izero(isy,ndim,1)
c
      nrep = 0
c
10    continue
      do 100 i=1,ndim
         if (isy(i).ne.0) go to 30
         nrep = nrep + 1
         isy(i) = nrep
30       if (i.eq.ndim) return
         jb = i+1
         do 50 j=jb,ndim
            if (ch(i,j)) go to 50
            if (isy(j).le.0.or.isy(j).eq.isy(i)) go to 40
            call comsym_servec(isy(i),isy(j),isy,ndim,nrep)
            go to 50
40          isy(j) = isy(i)
50       continue
100   continue
c
c
      return
      end
      subroutine comsym_servec(ir1,ir2,iro,norb,nrep)
c
c...  symmetry ir1 and ir2 are identical really
c...  change the highest to the lowest one
c...  a simplified version of the direct-ci routine
c
      implicit real*8 (a-h,o-z)
      dimension iro(norb)
c
      irf = max(ir1,ir2)
      irt = min(ir1,ir2)
      nrep = nrep - 1
c
10    do 20 i=1,norb
20    if (iro(i).eq.irf) iro(i) = irt
c
c...   shift everything to get a consecutive list of symm.numbers again
c
      if (irf.gt.nrep) return
      irt = irf
      irf = irf + 1
      go to 10
c
      end
      subroutine form_servec(vc,ndim,nmdim,occ,eig,iocc,ieig,ipr)
c
c...  print (ipr=1) or read (ipr=-1) a complete vector file 'formatted'
c...   this option is intended for transport of vectors as text
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension vc(*),occ(*),eig(*)
      character*8 dumbo(3),mblan,mvecto,mocc,meig,dumb,ites
      data dumbo/'(',' ',')'/,dumb/'5f16.9'/,mblan/' '/
      data mvecto,mocc,meig/'vectors','occ','eig'/
c
      if (ipr.lt.0) go to 1000
c...  print / read format first
      call inpa(dumbo(2))
      if (dumbo(2).eq.mblan) dumbo(2) = dumb
c
c... begin printing
c
      write(iwr,1)
1     format(//' read format')
      write(iwr,2) ndim,nmdim,dumbo(2)
2     format(/' vectors ',2i10,5x,a8)
      nnnd = ndim*nmdim
      write(iwr,dumbo)  (vc(i),i=1,nnnd)
      nn = nmdim
      if (iocc.le.0) nn = 0
      write(iwr,3) nn
3     format(/' occ ',i10)
      if (nn.gt.0) write(iwr,dumbo) (occ(i),i=1,nmdim)
      nn = nmdim
      if (ieig.le.0) nn = 0
      write(iwr,4) nn
4     format(/' eig ',i10)
      if (nn.gt.0) write(iwr,dumbo) (eig(i),i=1,nmdim)
      write(iwr,5)
5     format(/' end of format-print')
c
c...  end of print
c
      return
c
c...  read
c
1000  call rchar(ites)
      if (ites.ne.mvecto) go to 1000
      call inpi(ndim)
      call inpi(nmdim)
      call inpa(dumbo(2))
      if (dumbo(2).eq.mblan) dumbo(2) = dumb
      nnnd = ndim*nmdim
      read(ird,dumbo) (vc(i),i=1,nnnd)
      call rchar(ites)
      if (ites.ne.mocc) call caserr('error on read format, no occ')
      call inpi(nn)
      if (nn.gt.0) read(ird,dumbo) (occ(i),i=1,nmdim)
      iocc = 0
      if (nn.gt.0) iocc = 1
      call rchar(ites)
      if (ites.ne.meig) call caserr('error on read format, no eig')
      call inpi(nn)
      if (nn.gt.0) read(ird,dumbo) (eig(i),i=1,nmdim)
      ieig = 0
      if (nn.gt.0) ieig = 1
      write(iwr,6)
6     format(/' end of format-read')
c
      return
      end
      subroutine locmp(u,sover,iloc,g,ntot,norb,nrmx,var)
c
      implicit real*8 (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer g(*),iloc(*)
      real*8 u(*),sover(*)
      character*8 ites,iend
c
c      localization
c     v.magnasco,a.perico jcp 47,971 (1967)
c
c   ntot  =  size of subset
c   norb  =  dimension of orbital basis
c   nrmx=maximum number of cycles to localize
c   var =accuracy of localization
c
c     u orbitals (ntot*norb)
c     sover  overlap matrix (norb*(norb+1)/2)
c     iloc  scratch array (norb)
c     g     scratch array (ntot*norb)
c
      data simin/1.0d-7/,iend/'end'/
c
      nmx=ntot*norb
      do 10 i=1,nmx
 10   g(i) = 0
      nrot = 0
      nmoc = 0
c
      write(iwr,601)
601   format(/' magnasco/perico localisation ')
c
 20   call rchar(ites)
      if (ites.ne.iend) then
         call cinput(jrec,jump,-1) 
         call rdistr(iloc,max,norb)
         nmoc = nmoc + 1
         write(iwr,602) nmoc,(iloc(j),j=1,max)
 602     format(/,1x,i3,(2x,40i3))
         idx=norb*(nmoc-1)
         do j=1,max
            jdx=idx+iabs(iloc(j))
            g(jdx)=isign(1,iloc(j))
         enddo
         go to 20
      end if
c...   print if it fits on the screen easily
      if (norb.le.10) then
         il = 0
         do 50 j=1,nmoc
            is = il + 1
            il = il + norb
50       write(iwr,603) (g(i), i = is,il)
603      format(1x,(/1x,40i3))
      end if
c
c...   start rotating the orbitals
c
60    sd = 0.0d0
         nrot=nrot+1
         do 90 i=1,nmoc
            idx=norb*(i-1)
            iplus=i+1
            if (iplus.gt.ntot) go to 100
            do 90 j=iplus,ntot
               jdx=norb*(j-1)
               sa=0.0d0
               sb=0.0d0
               ndx=0
               do 70 m=1,norb
                  im=idx+m
                  jm=jdx+m
                  do 70 n=1,m
                     ndx=ndx+1
                     in=idx+n
                     jn=jdx+n
                     if (j.le.nmoc) then
                        fac=g(im)*g(in)-g(jm)*g(jn)
                     else
                        fac=g(im)*g(in)
                     end if
                     if (fac.ne.0.0d0) then
                         fac=fac*sover(ndx)
                         if (m.ne.n) fac=fac+fac
                         sa=sa+fac*(u(im)*u(in)-u(jm)*u(jn))
                         sb=sb+fac*(u(im)*u(jn)+u(jm)*u(in))
                     end if
70             continue
               sc= dsqrt(sa*sa+sb*sb)
               if (sc.ne.0.0d0) then
                  cosine=(sa+sc)/(sc+sc)
                  if(cosine.lt.0.0d0) cosine=0.0d0
                  sine=1.0d0-cosine
                  if (sine.ge.simin) then
                     sine = dsign(dsqrt(sine),sb)
                  else
                     sine = sb/(sa+sa)
                     sine = sine*(1-1.5d0*sine**2)
                  end if
                  cosine = dsqrt(cosine)
                  do 80 m=1,norb
                     im=idx+m
                     jm=jdx+m
                     sa=u(im)
                     sb=u(jm)
                     uim = sa*cosine + sb*sine
                     ujm = -sa*sine + sb*cosine
                     sd = sd + dabs(u(im)-uim) + dabs(u(jm)-ujm)
                     u(im) = uim
                     u(jm) = ujm
80                continue
               end if
90       continue
100      continue
      if (var.lt.sd.and.nrot.lt.nrmx) goto 60
c
c...   end of calc
c
      write(iwr,604) nrot,sd
604   format(5h0****,' number of rotations :',i3,'  convergence : ',
     1 e12.5,5h ****)
      if (nrot.gt.nrmx) write(iwr,605)
605   format(/' #####   warning  no convergence ##### '/)
      write(iwr,606)
606   format(/' *******  orbitals localised *****'/)
c
      return
      end
      subroutine ini_ctrans
c
      implicit none
c
c...  create a unit ctrans
c...  and save/restore original trans
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      real*8 ctran2
      integer ilifc2, ntran2, itran2
      logical otran2,otri2
      common/junk2/ilifc2(maxorb),ntran2(maxorb),itran2(mxorb3),
     +             ctran2(mxorb3),otran2,otri2
      integer i
c
      do i=1,maxorb
         ntran(i) = 1
         itran(i) = i
         ctran(i) = 1.0d0
         ilifc(i) = i-1
      end do
      otran = .false.
c
      return
c
      entry save_ctrans
c...    save ctrans-block so that servec may freely mess around
      call dcopy(mxorb3,ctran,1,ctran2,1)
      call icopy(2*maxorb+mxorb3,ilifc,1,ilifc2,1)
      otran2 = otran
      otri2 = otri
c
      return
c
      entry rest_ctrans
c...    restore ctrans-block at the end of servec 
      call dcopy(mxorb3,ctran2,1,ctran,1)
      call icopy(2*maxorb+mxorb3,ilifc2,1,ilifc,1)
      otran = otran2
      otri = otri2
c
      return
      end
      subroutine atommo(vc,ndim,nmo,iat,natom,maskv,vv,iato,iperm)
c
c...  sort mo's in de order of atoms
c     call atommo(vc,ndim,nmo,ijun(1,1),natom,ijun(1,3),hulp,
c    1            ijun(1,4),ijun(1,2))
c
      implicit real*8 (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension vc(ndim,nmo),maskv(ndim),vv(natom)
      dimension iat(natom),iato(natom),iperm(nmo)
c
      iii = 0
      do i=1,ndim
         maskv(i) = 0
         if (nat.lt.100) then 
            read(zbflab(i)(1:2),'(i2)') iatom
         else if (nat.ge.100.and.nat.lt.1000) then 
            read(zbflab(i)(1:3),'(i3)') iatom
         else if (nat.ge.1000.and.nat.lt.10000) then 
            read(zbflab(i)(1:4),'(i4)') iatom
         else if (nat.ge.10000.and.nat.lt.100000) then 
            read(zbflab(i)(1:5),'(i5)') iatom
         else 
            call caserr(' nat > 100000 not allowed in servec')
         end if
         do j=1,natom
            if (iat(j).eq.iatom) maskv(i) = j
         end do
         iii = iii + maskv(i)
      end do
      if (iii.eq.0) call caserr('geometry/basis is needed in atommo')
c
c...  asign mo's to atoms 
c
      do j=1,nmo
         iato(j) = 0
         do i=1,natom
            vv(i) = 0.0d0
         end do
         do i=1,ndim
           if (maskv(i).gt.0) vv(maskv(i))=vv(maskv(i))+vc(i,j)*vc(i,j)
         end do
         dm = 0.0d0
         do i=1,natom
            if (vv(i).gt.dm) then
               dm = vv(i)
               iato(j) = i
            end if
         end do
      end do
c
c...  translate into permutation
c
      kk = 0
      do i=1,natom
         do j=1,nmo
            if (iato(j).eq.i) then
               kk = kk + 1
               iperm(kk) = j
            end if
         end do
         iat(i) = kk
      end do
      do i=natom,2,-1
         iat(i) = iat(i) - iat(i-1)
      end do
      write(iwr,10) nmo,natom,(iato(i),i=1,nmo)
10    format(' ******** assigning',i4,' mos to',i4,' atoms ********',
     1       (/3x,20i5),1x)
      write(iwr,11) (iat(i),i=1,natom)
11    format('   # orbitals / atom ',(t21,20i5),1x)
      if (kk.ne.nmo) call caserr('atommo error')
c
c...  make order sort of unique
c
      ne = 0
      do i=1,natom
         nb = ne + 1
         ne = ne + iat(i)
         do j=nb,ne
            do k=1,ndim
               if (maskv(k).eq.i.and.
     1             dabs(vc(k,iperm(j))).gt.1.0d-2)  go to 20
            end do
            call caserr('atommo fails to find significant component')
20          do jj=j+1,ne
               if (dabs(vc(k,iperm(jj))).gt.dabs(vc(k,iperm(j)))) then
                  itemp = iperm(j)
                  iperm(j) = iperm(jj)
                  iperm(jj) = itemp
               end if
            end do
         end do
         nb = ne + 1
      end do
c
      return 
      end
c</BODY>
c</HTML>
