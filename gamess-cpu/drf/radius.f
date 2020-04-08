      function radius(name,iradex,charge,afct)
c=======================================================================
c     this function assigns a radius to an atom
c     value is in bohr according to :
c
c     bondi:j.phys.chem.(1964) 68 441
c     table i:  rw
c
cmas  changed function feb 7, 1996
c=======================================================================
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c * * * this is for communication with hondo8
c
      character *8 errmsg(3)
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      parameter (ndif=106)
      character*2 names(ndif),name*4
      dimension rad(ndif)
      dimension afct(2)
c
      parameter ( nfrec = 18 )
      dimension rfrec0(nfrec),rfrec1(nfrec),rfrec2(nfrec),rfrec3(nfrec)
c
      data errmsg/'program','stop in','-radius-'/
      data names/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     *           'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     *           'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     *           'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     *           'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     *           'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     *           'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     *           'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     *           'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     *           'pa','u ','np','pu','am','cm','bk','cf','es','fm',
     *           'md','no','lr','e','qq','xx'/
c
c    -----  note: e refers to a center without nuclear charge (bsse)
c
      data rad  /
c         h        he       li       be       b        c
     1 0.2267d1,0.2646d1,0.0000d1,0.0000d1,0.0000d1,0.3213d1,
c         n        o        f        ne       na       mg
     1 0.2929d1,0.2872d1,0.2778d1,0.2910d1,0.0000d1,0.0000d1,
c         al       si       p        s        cl       ar       k
     1 0.0000d1,0.3968d1,0.3401d1,0.3402d1,0.3307d1,0.3553d1,0.0000d1,
c         ca       sc       ti       v        cr       mn
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         fe       co       ni       cu       zn       ga
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         ge       as       se       br       kr       rb
     1 0.0000d1,0.3500d1,0.3590d1,0.3500d1,0.3817d1,0.0000d1,
c         sr       y        zr       nb       mo       tc
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         ru       rh       pd       ag       cd       in
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         sn       sb       te       i        xe       cs
     1 0.0000d1,0.0000d1,0.3892d1,0.3741d1,0.3817d1,0.0000d1,
c         ba       la       ce       pr       nd       pm
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         sm       eu       gd       tb       dy       ho
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         er       tm       yb       lu       hf       ta
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         w        re       os       ir       pt       au
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         hg       tl       pb       bi       po       at
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         rn       fr       ra       ac       th       pa
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         u        np       pu       am       cm       bk
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         cf       es       fm       md       no       lr
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         e        qq       xx
     1 0.0000d1,0.0000d1,0.1000d1/
c
c=======================================================================
c     data obtained by fitting of polynoma (3rd degree) to rvdw-values
c     for ions with charge ranging from q=-2 to q=+2
c     by monte carlo sampling of ab initio charge densities
c
c                                        2         3
c     rvdw(q) = rvdw(0) + r1 * q + r2 * q  + r3 * q
c
c     values are in a (angstroms) and e (electron charge)
c
c     fitting done by vladimir frecer
c=======================================================================
c-----------------------------------------------------------------------
c     h
c-----------------------------------------------------------------------
      data rfrec0(1),rfrec1(1),rfrec2(1),rfrec3(1)/
     .       1.520d0, -1.054d0, -0.293d0,  0.000d0/
c-----------------------------------------------------------------------
c     he
c-----------------------------------------------------------------------
      data rfrec0(2),rfrec1(2),rfrec2(2),rfrec3(2)/
     .       1.330d0, -0.300d0, -0.129d0, -0.037d0/
c-----------------------------------------------------------------------
c     li
c-----------------------------------------------------------------------
      data rfrec0(3),rfrec1(3),rfrec2(3),rfrec3(3)/
     .       2.210d0, -0.996d0,  0.003d0,  0.096d0/
c-----------------------------------------------------------------------
c     be
c-----------------------------------------------------------------------
      data rfrec0(4),rfrec1(4),rfrec2(4),rfrec3(4)/
     .       2.180d0, -0.335d0, -0.103d0, -0.051d0/
c-----------------------------------------------------------------------
c     b
c-----------------------------------------------------------------------
      data rfrec0(5),rfrec1(5),rfrec2(5),rfrec3(5)/
     .       2.060d0, -0.298d0, -0.007d0,  0.009d0/
c-----------------------------------------------------------------------
c     c
c-----------------------------------------------------------------------
      data rfrec0(6),rfrec1(6),rfrec2(6),rfrec3(6)/
     .       1.930d0, -0.234d0, -0.013d0,  0.008d0/
c-----------------------------------------------------------------------
c     n
c-----------------------------------------------------------------------
      data rfrec0(7),rfrec1(7),rfrec2(7),rfrec3(7)/
     .       1.800d0, -0.170d0, -0.015d0,  0.003d0/
c-----------------------------------------------------------------------
c     o
c-----------------------------------------------------------------------
      data rfrec0(8),rfrec1(8),rfrec2(8),rfrec3(8)/
     .       1.700d0, -0.140d0, -0.014d0,  0.002d0/
c-----------------------------------------------------------------------
c     f
c-----------------------------------------------------------------------
      data rfrec0(9),rfrec1(9),rfrec2(9),rfrec3(9)/
     .       1.610d0, -0.096d0, -0.012d0,  0.004d0/
c-----------------------------------------------------------------------
c     ne
c-----------------------------------------------------------------------
      data rfrec0(10),rfrec1(10),rfrec2(10),rfrec3(10)/
     .       1.510d0, -0.098d0, -0.005d0,  0.003d0/
c-----------------------------------------------------------------------
c     na
c-----------------------------------------------------------------------
      data rfrec0(11),rfrec1(11),rfrec2(11),rfrec3(11)/
     .       2.240d0, -0.835d0, -0.009d0,  0.094d0/
c-----------------------------------------------------------------------
c     mg
c-----------------------------------------------------------------------
      data rfrec0(12),rfrec1(12),rfrec2(12),rfrec3(12)/
     .       2.420d0, -0.375d0, -0.062d0, -0.039d0/
c-----------------------------------------------------------------------
c     al
c-----------------------------------------------------------------------
      data rfrec0(13),rfrec1(13),rfrec2(13),rfrec3(13)/
     .       2.410d0, -0.365d0,  0.021d0,  0.005d0/
c-----------------------------------------------------------------------
c     si
c-----------------------------------------------------------------------
      data rfrec0(14),rfrec1(14),rfrec2(14),rfrec3(14)/
     .       2.330d0, -0.291d0,  0.007d0,  0.006d0/
c-----------------------------------------------------------------------
c     p
c-----------------------------------------------------------------------
      data rfrec0(15),rfrec1(15),rfrec2(15),rfrec3(15)/
     .       2.260d0, -0.275d0,  0.015d0,  0.006d0/
c-----------------------------------------------------------------------
c     s
c-----------------------------------------------------------------------
      data rfrec0(16),rfrec1(16),rfrec2(16),rfrec3(16)/
     .       2.150d0, -0.240d0,  0.020d0,  0.004d0/
c-----------------------------------------------------------------------
c     cl
c-----------------------------------------------------------------------
      data rfrec0(17),rfrec1(17),rfrec2(17),rfrec3(17)/
     .       2.080d0, -0.208d0,  0.017d0,  0.001d0/
c-----------------------------------------------------------------------
c     ar
c-----------------------------------------------------------------------
      data rfrec0(18),rfrec1(18),rfrec2(18),rfrec3(18)/
     .       1.970d0, -0.251d0,  0.010d0,  0.017d0/
      data zero,one,two,three/0.0d0,1.0d0,2.0d0,3.0d0/
      data bohr/0.529177249d0/
c-----------------------------------------------------------------------
c     begin
c-----------------------------------------------------------------------
cafc
cafc  gamess-names for dummy-atoms
      if (name(1:2) .eq. 'bq') then
        radius = 0.0d0
        return
      endif
      do i=1,ndif
        if (name(:2).eq.names(i))  goto 20
      enddo
c-----------------------------------------------------------------------
c     error : symbol not found
c-----------------------------------------------------------------------
      write(iwr,*) ' no radius found for ',name
      call hnderr(3,errmsg)
c
   20 if (iradex.eq.4.and.i.le.nfrec) then
c-----------------------------------------------------------------------
c     use data obtained by vladimir frecer
c     remember that radius is then in angstroms
c-----------------------------------------------------------------------
         radius = rfrec0(i) + rfrec1(i)*charge + rfrec2(i)*charge**two
     .                      + rfrec3(i)*charge**three
         radius = radius / bohr
      elseif (iradex.eq.4) then
c-----------------------------------------------------------------------
c     try to use frecer's model for an atom which was not fitted
c-----------------------------------------------------------------------
         write (iwr,*) 'no frecer values for ' ,name
         call hnderr(3,errmsg)
      elseif (iradex.eq.0) then
c-----------------------------------------------------------------------
c     use value from table
c-----------------------------------------------------------------------
         radius = rad(i)
      else
c-----------------------------------------------------------------------
c     compute radius from polarizability if iradex = 1 or 2
c-----------------------------------------------------------------------
         radius = afct(iradex)*alfa(name,iradex,ier)**(one/three)
      endif
      return
      end
      subroutine drfzfa(xexp)
c------
c      calculates and updates the (expanded) field due to
c      external charges at the expansion centra of the
c      quantum mechanical system
c
c      the interaction energy between the external charges
c      and the nuclei (extnuc) is also calculated
c
c      penetration (overlap) effects
c      resulting from close contacts
c      can be accounted for (modxza flag)
c
c      --------  p.th. van duijnen, ibm-kingston 1985, and
c               groningen, dec. 1991.
c
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  dummy variables
c
      dimension xexp(3,nexp)
c
c-----  common blocks
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
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
      character*8 bflab
      common/mollab/mwanam(128),mwbnam(128),bflab(1024)
c
      character*8 zruntp, zguess, zconf, zstate
      character*8 zcom,title,ztagg,zsymm
      character*8 zorb,zpseud
      character*8 anam
      character*10 bnam
      common/runlab/zcom(19),title(10),anam(maxat),bnam(maxorb),
     +           ztagg(maxat),zsymm(7),zorb(maxorb),zpseud(maxat)
c
c
cmw      include 'drf/dimpar'
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 zfa
      common/drfzfa1/zfa(4*(mxat +1),mxgran)
      real*8 zpa
      common/drfzfa2/zpa(mxat +1,mxgran)
      real*8 dzfa
      common/derzfa1/dzfa(12*(mxat+1),mxgran)
      real*8 dpseu
      common/derzfa2/dpseu(3,(mxat+1),mxgran)
c
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
c-----  local variables
c
      dimension p(3), q(3), pq(3), w(3)
      dimension b(3,3)
c
c
      character*80 card
      character*16 namj, nami
      logical group
      logical ipol
c
      data third/.33333333333333333333333333333d00/
      data sixth/.16666666666666666666666666667d00/
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,one,pt5,pt75/0.0d00,1.0d00,0.5d00,0.75d00/
      data small /1.0d-03/
c
c-----  begin
c
      ier = 0
c
c-----  read expansion centers.
c       = # internal  atoms
c       (+ center of internal nuclear charge)
c
      call daread(idafdrf,iodadrf,xexp,3*nexp,1)
c
      if (mcupdt) then
c 1-----
c  -----  if update of zfa is required, read zfa from disk
c
        call daread(idafdrf,iodadrf,zfa,4*nzfa,3)
        call daread(idafdrf,iodadrf,zpa,nzfa,4)
        call daread(idafdrf,iodadrf,dzfa,12*nzfa,5)
        call daread(idafdrf,iodadrf,dpseu,3*nzfa,6)
c
c  -----  update field at expansion centra according to
c         monte carlo move
c
        call expupd(xexp,extnuc(iactst),repmod(iactst),
     1              extnucg,repmodg,iactst)
c
        call dawrit(idafdrf,iodadrf,zfa,4*nzfa,76,navdrf)
        call dawrit(idafdrf,iodadrf,zpa,nzfa,77,navdrf)
        call dawrit(idafdrf,iodadrf,dzfa,12*nzfa,78,navdrf)
        call dawrit(idafdrf,iodadrf,dpseu,3*nzfa,79,navdrf)
c 1-----
      else
c 1-----
c  -----  calculate zfa, dzfa (derivative) and model repulsion
c
        pseufac = pt75
c
c  -----  initialize pseudo repulsion and nuclear-external charge
c         interaction
c
        epseu = zero
        extnuc(iactst) = zero
        repmod(iactst) = zero
        do 400, igr = 1, ngran
          extnucg(igr,iactst) = zero
          repmodg(igr,iactst) = zero
  400   continue
c
c  -----  nspec counts the external points close to internal centra,
c         for which a model repulsion is calculated, and for which
c         the "reverse" dispersion is calculated, if required
c
        nspec = 0
c
        call clear(zfa,4*nzfa)
        call clear(zpa,nzfa)
        call clear(dzfa,12*nzfa)
        call clear(dpseu,3*nzfa)
c
c  -----  loop over the expansion centra
c
        do 500, j = 1, nexp
c   2-----
c    -----  exclude ambiguous atom: its interaction with the
c           discrete classical system is treated classically
c           in subr. drfzfp
c
          if (nambpt(j) .ne. 0) goto 500
c
          if (j .le. nat) then
c     3-----
c      -----  current expansion centre is at a nucleus
c
            za = czan(j)
            namj = anam(j)
            alfj = alfat(j)
c           alfj = alfa(namj,0,ier)
c
c      -----  # of  valence electrons
c
            call drfnval(namj,nvalj,znuc)
c
c      -----  scaling factor for slater-kirkwood estimate
c             of dispersion and repulsion energy
c
            if (alfj .ne. zero) then
              aovernj= sqrt(alfj/nvalj)
            endif
c     3-----
          else
c     3-----
c      -----  non-nuclei are given a polarizability of 1.
c
            alfj = one
c     3-----
          endif
c
c    -----  pointers to arrays w.r.t. expansion center -j-
c
          jf = (j-1)*4
          nzf = jf*3
c
c    -----  position vector -j- into -q-
c
          q(1) = xexp(1,j)
          q(2) = xexp(2,j)
          q(3) = xexp(3,j)
c
c    -----  pointer to external points in -xtpts-
c           with special properties
c           imp: polarizable external points
c
          imp = 1
c
c    -----  loop over external points
c
          do 300, ii = 1, nxtpts
c     3-----
c      -----  check whether external point is polarizable
c
            ipol = ii .eq. mpol(imp)
c
c      -----  check whether external point is an ambiguous
c             atom. the interaction with the qm system is taken care of
c             by inclusion into the qm system (whereby it is treated exa
c
            if (ncutpt(ii) .ne. 0) goto 290
c
c      -----  position vector of external point -ii- in -p-
c
            p(1) = xpts(1,ii)
            p(2) = xpts(2,ii)
            p(3) = xpts(3,ii)
c
            call distab(q,p,pq,dist)
c
            if (dist .lt. small) then
              write(iwr,95) anam(j), nxcent(ii), dist
  95          format(/,1x,'WARNING: interaction ',
     1   'between ', a16, ' and ', a16, ' at', e15.8,
     2   ' bohr distance: skipped')
              goto 300
            endif
c
            dmind1 = one/dist
            dmin3 = dmind1*(dmind1**2)
            factp = one
            factf = one
            v = one
c
            alfi = alfext(ii)
            if (alfi .eq. zero) alfi = one
c
c           if (chrg(ii) .eq. zero) goto 250
c
c      -----  get index of analysis group of external point -ii-
c
            igr = igranl(ii)
c
c      -----  account for penetration effects (optional)
c
            if (modxza .ne. 0) then
c       4-----
c        -----  set distance criterium for penetration effects
c               depending on polarizabilities of -j- and -ii-
c
              s = (alfj*alfi)**sixth
              v = dist/s
c
c        -----  check distance criterium
c               modify interactions (thole's cone model)
c
              if (ithole .eq. 1) then
c         5-----
                if (v .le. afact) then
                  av = v/afact
                  factp = av**4 - two*av**3 + two*av
                  factf = four*av**3 - three*av**4
                endif
              else if (ithole .eq. 2) then
                au = afact*v
                factp = (one - (pt5*au + one)*exp(-au))
                factf = (one - (pt5*au**2 + au + one)*exp(-au))
c         5-----
              endif
c       4-----
            endif
            factp = factp*dmind1
            factf = factf*dmin3
c
c      -----  update coulomb interaction with internal nuclei
c
            if ((field(:4) .ne. ' ') .and. (j .le. nat)) then
              extnuca = za*chrg(ii)*factp
              extnucg(igr,iactst) = extnucg(igr,iactst) + extnuca
              extnuc(iactst) = extnuc(iactst) + extnuca
            endif
c
c      -----  update -zfa-
c
c             -zfa- : sum(p) {f(p;q) + v(p;q) - f(p;q).q} * q(p)
c
c             the field operator and potential of all external charges
c             in p (-p-) at expansion centre q (-j- = -q-)
c
            sum = zero
            do 210, k = 1, 3
              zfa(jf+k,igr) = zfa(jf+k,igr) + chrg(ii)*(pq(k)*factf)
              sum = sum - q(k)*(pq(k)*factf)
  210       continue
c
c      -----  the expanded term f(p;q).q is not to be added if
c             the field of the external charges is treated
c             as a perturbation and the qm charge distribution
c             is modelled by point charges (and dipoles)
c
            zfa(jf+4,igr) = zfa(jf+4,igr) + chrg(ii)*(sum+factp)
c
c      -----  update -zpa- : the potential of the external charges
c             at the expansion centra
c
            zpa(j,igr) = zpa(j,igr) + chrg(ii)*factp
c
c      -----  classical field gradient tensor into -b-
c
            call clear(b,9)
c
c      -----  b(3,3) = t(p;q) : field gradient of charge in -p- at -q-
c             the gradient is scaled according to v
c
            call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c      -----  multiply position vector of -j- with -b- to give -w-
c
            call matvec(b,q,w,3,.false.)
c
c      -----  update derivative of -zfa-
c
            do 214, k = 1, 3
              do 213, l = 1, 3
                dzfa(nzf+(k-1)*4+l,igr) = dzfa(nzf+(k-1)*4+l,igr)
     1                               - b(k,l)*chrg(ii)
  213         continue
              dzfa(nzf+(k-1)*4+4,igr) = dzfa(nzf+(k-1)*4+4,igr)
     1                              + w(k)*chrg(ii)
  214       continue
c
  250       continue
c
c      -----  non bonding repulsion taken from charmm:
c             j.comp.chem. 4 (1983) 213
c
c      -----  this is called the pseudo- or model repulsion throughout
c             the program
c
            if ((field(:4) .ne. ' ') .and.
     1          (j .le. nat) .and.
     1          (iqmclr .eq. 1) .and. 
     1          (namj(:2) .ne. 'bq') .and.
     1          (nxcent(ii)(:2) .ne. 'qq') .and.
     2          (nxcent(ii)(:2) .ne. 'e ') .and.
     3          (nxcent(ii)(:2) .ne. '  ')) then
c       4-----
c        -----  this might be a close atom: apply distance
c               criterium for pauli repulsion
c
              if (dist .lt. dstmin) then
c         5-----
c          -----  store external point counter in array relating to
c                 atom, for calculation of exact integrals
c
c          -----  calculate pauli repulsion and
c                 derivative, using some model potential
c
c               if (irepopt .eq. 0) then
c
c                 the necessary ingredients:
c
c            -----  # valence electrons
c
                  nami=nxcent(ii)
c
c            -----  check if the point represents a group
c
                  if (nami(:5) .eq. 'group') then
c           6-----
c              -----  non-bonded interactions evaluated only between
c                     internal atoms and individual members
c                     of a group:
c                     skip the group-representing point
c
                    goto 290
                  else
c
c              -----  get the number of valence electrons of the
c                     external point
c
                    call drfnval(nxcent(ii),nvali,znuc)
c
c              -----  check if it was already found
c                     close to another internal atom
c
c                   do 220, nsp = 1, nspec
c                     if (nspecl(nsp) .eq. ii) goto 230
c 220               continue
c
c              -----  this external point has not been found
c                     close to an internal atom:
c                     add it to the list
c
c                   nspec = nspec + 1
c                   nvalel(nspec) = nvali
c                   polars(nspec) = alfa(nxcent(ii),0,ier)
c                   polars(nspec) = alfext(ii)
c                   nspecl(nspec) = ii
c 230               continue
c
c              -----  1/r12 term for selected atoms
c
c              -----  get polarizability of external atom
c
c                   if (ipol) then
c                     alfi = polar(imp)
c                   else
c                     alfi = alfa(nami,0,ier)
c                   endif
                    alfi = alfext(ii)
                                                
                    if (alfi .ne. zero) then
                      aoverni = sqrt(alfi/nvali)
                    endif
c
c              -----  calculate equilibrium lj distance from
c                     "size" of internal and external atom
c
c              -----  size of external atom: ri
c
                    ri = radext(ii)*rfact
c
c              -----  size of external atom: rjx
c
                    rjx = radat(j)*rfact
c
c              -----  in case of a h-bond, the sizes need to be
c                     modified
c
c              -----  check h-bond criterium
c
                    if ((ihbond .eq. 1) 
     1                  .and. (dist .le. hbondl)) then
                      if(namj(:2).eq.'h'.and.
     1                  (nami(:2).eq.'n'.or.
     2                   nami(:2).eq.'o'.or.
     3                   nami(:2).eq.'f')) rjx=hbondr
                      if(nami(:2).eq.'h'.and.
     1                  (namj(:2).eq.'n'.or.
     2                   namj(:2).eq.'o'.or.
     3                   namj(:2).eq.'f')) ri=hbondr
                    endif
cmas  changed feb 7, 1996
c
c              -----  calculate the 1/r12 lj term for the close contact
c
                    if (alfj .ne. zero .and. alfi .ne. zero) then
                      pseu = pseufac*alfj*alfi/(aovernj+aoverni)
                      pseu1 = (ri+rjx)**6
                      repmodg(igr,iactst) = repmodg(igr,iactst) +
     1                                      pseu*pseu1*(dmind1**12)
                      epseu = epseu + pseu*pseu1*(dmind1**12)
c
c                -----  and derivative
c
                      do 270, k = 1, 3
                        dpseu(k,j,igr) = dpseu(k,j,igr) +
     1                         twelve*pq(k)*pseu*pseu1*(dmind1**14)
  270                 continue
                    endif
c          6-----
                  endif
c               else if (irepopt .eq. 1) then
c
c            -----  set some necessary blabla for integrals
c
c                  call drfpseu
c
c               endif
c         5-----
              endif
c       4-----
            endif
  290       if (ipol) imp = imp + 1
c     3-----
  300     continue
c   2-----
  500   continue
c
c  -----  write zfa, zpa, dzfa and dpseu to disk.
c
        call dawrit(idafdrf,iodadrf,zfa,4*nzfa,3,navdrf)
        call dawrit(idafdrf,iodadrf,zpa,nzfa,4,navdrf)
        call dawrit(idafdrf,iodadrf,dzfa,12*nzfa,5,navdrf)
        call dawrit(idafdrf,iodadrf,dpseu,3*nzfa,6,navdrf)
c
        repmod(iactst) = epseu
c 1-----
      endif
c
      if (idrfout .eq. 5) then
        do 600, igr = 1, ngran
          call hatout(dzfa(1,igr),4,3*nexp,5,'dzfa')
  600   continue
      endif
c
      if (idrfout .eq. 3) then
        do 700, igr = 1, ngran
          call hatout(zfa(1,igr),4,nexp,5,'zfa')
          call hatout(zpa(1,igr),1,nexp,5,'zpa')
  700   continue
      endif
c
      return
      end
