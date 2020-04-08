c
c     Routines interfacing between HONDO and GAMESS-UK
c     ================================================
c
c     These routines are needed by:
c     - drf
c     - nmr
c
c-----------------------------------------------------------------------
c
      subroutine dsxyz
      implicit none
c
c     dsxyz -> dsxyzh [analg.m] : gauss-hermite quadrature using 
c                                 minimum point formula
c
      call dsxyzh
      end
c
c-----------------------------------------------------------------------
c
      subroutine dvxyz
      implicit none
c
c     dvxyz -> dvxyzh [analg.m] : gauss-hermite quadrature using 
c                                 minimum point formula
c
      call dvxyzh
      end
c
c-----------------------------------------------------------------------
c
      subroutine fiflsh(iunit)
      implicit none
c
c     filflsh -> flushn [util3] : simply flush a fortran unit number
c
      integer iunit
      call flushn(iunit)
      end
c
c-----------------------------------------------------------------------
c
      subroutine hnderr(n,a)
      implicit none
c
c     hnderr -> caserr [util1.m] : print error message and exit
c
      integer n
      character*8 a(n)
      integer i
      character*132 errmsg
      write(errmsg,*)(a(i),i=1,n)
      call caserr(errmsg)
      end
c
c-----------------------------------------------------------------------
c
      subroutine initpk(a)
      implicit none
c
c     initpk -> caserr : initpk is a stub in HONDO and should not be 
c                        called.
c
      character*(*) a
      call caserr2("initpk: should not get here")
      end
c
c-----------------------------------------------------------------------
c
      subroutine initpn(a)
      implicit none
c
c     initpn -> caserr : initpn is a stub in HONDO and should not be 
c                        called.
c
      character*(*) a
      call caserr2("initpn: should not get here")
      end
c
c-----------------------------------------------------------------------
c
      subroutine isoout
      implicit none
      call caserr("isoout: have a look at this")
      end
c
c-----------------------------------------------------------------------
c
      subroutine lcpadd(a,b,n)
      implicit none
c
c     lcpadd -> pg_dgop [parallel.m] : global sum (ignore scratch space
c                                      in argument b)
c
      integer n
      real*8 a(n),b(n)
      call pg_dgop(6001,a,n,'+')
      end
c
c-----------------------------------------------------------------------
c
      subroutine lcpsnc
      implicit none
c
c     lcpsnc -> pg_synch [parallel.m] : synchronize processes
c
      call pg_synch(6002)
      end
c
c-----------------------------------------------------------------------
c
      real*8 function prpdot(a,b,n)
      implicit none
c
c     prpdot -> tracep [util1.m] : trace of product of 2 symmetric
c                                  matrices
c
      integer n
      real*8 a(n), b(n)
      real*8 tracep
      prpdot = tracep(a,b,n)
      end
c
c-----------------------------------------------------------------------
c
      subroutine prtr(d,n)
      implicit none
c
c     prtr -> prtri : print a triangular matrix
c
      integer n
      real*8 d(n)
      call prtri(d,n)
      end
c
c-----------------------------------------------------------------------
c
      subroutine rewfil(iunit)
      implicit none
c
c     rewfil -> rewftn [machscf.m] : rewind fortran stream
c
      integer iunit
      call rewftn(iunit)
chvd  rewind(iunit)
      end
c
c-----------------------------------------------------------------------
c
      subroutine root4
      implicit none
c
c     root4 -> roots4 [intega.m] : root for Gauss-Rys quadrature
c
      call roots4
      end
c
c-----------------------------------------------------------------------
c
      subroutine root5
      implicit none
c
c     root5 -> roots5 [intega.m] : root for Gauss-Rys quadrature
c
      call roots5
      end
c
c-----------------------------------------------------------------------
c
      subroutine sxyz
      implicit none
c
c     sxyz -> sxyzh [analg.m] : gauss-hermite quadrature using 
c                               minimum point formula
c
      call sxyzh
      end
c
c-----------------------------------------------------------------------
c
      subroutine daopen(idafh,ioda,navh,ndar)
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
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
       character *132 filatt
       logical ok
c
c ---  lrecl = logical record length
c
c ----- open a direct access file -----

c        idafh ..... fortran file number
c        ioda ..... index array of length -ndar-
c        navh  ..... associated variable
c        ndar ..... number of -da- index words
c                   it must be equal to the number of logical records to
c                   be written on -idafh-.
c
      parameter (lrec1=1024,lrec2=4096)
c
      integer isingl, nbits
      common /machin1/ isingl,nbits
c
      integer ixddaf
      integer ndar10, ndar20, ndar31, ndar41, ndar42, ndar43, ndar44
      integer ndar45, ndar46, ndar47, ndar48, ndar49
      common /dafindx/ ixddaf(8192),ndar10,ndar20,
     +                 ndar31,ndar41,ndar42,ndar43,ndar44,ndar45,
     +                 ndar46,ndar47,ndar48,ndar49
c
      integer lrec10, lrec20, lrec31, lrec41, lrec42, lrec43
      integer lrec44, lrec45, lrec46, lrec47, lrec48, lrec49
      common /dafrec/ lrec10,lrec20,lrec31,lrec41,lrec42,lrec43,
     +                lrec44,lrec45,lrec46,lrec47,lrec48,lrec49
c
      common/junk2/navhh,navhhh,ibuff(1)
      dimension ioda(2,ndar)
      dimension iunits(12)
      data iunits /
     +    110, 120,  131,   141,   142,  143,    144,   145,
c        da10 da20  da31 neqinf neqsta neqcst neqind neqino
c       =mt0  =mt1  =mt2   =mt3   =mt4   =mt5   =mt6   =mt7
c
     +    146,   147,    148,   149 /
c      neqpol neqdis  neqrep  neqrqm
c       =mt8    =mt9   =mt10   =mt12
c
      data m1/1/
c
c     read the directory and the first available (file) record
c     of an existing da-file.
c     if it does not exist(new file), or it is too short, the directory
c     will be initialized with zero's
c
c     save dimension of index array
c
      lrec10 = lrec1
      lrec20 = lrec2
      lrec31 = lrec1
      lrec41 = lrec1
      lrec42 = lrec1
      lrec43 = lrec1
      lrec44 = lrec1
      lrec45 = lrec1
      lrec46 = lrec1
      lrec47 = lrec1
      lrec48 = lrec1
      lrec49 = lrec1
      if (idafh.eq.110) then
        ndar10=ndar
      else if (idafh.eq.120) then
        ndar20=ndar
      else if (idafh.eq.131) then
        ndar31=ndar
      else if (idafh.eq.141) then
        ndar41=ndar
      else if (idafh.eq.142) then
        ndar42=ndar
      else if (idafh.eq.143) then
        ndar43=ndar
      else if (idafh.eq.144) then
        ndar44=ndar
      else if (idafh.eq.145) then
        ndar45=ndar
      else if (idafh.eq.146) then
        ndar46=ndar
      else if (idafh.eq.147) then
        ndar47=ndar
      else if (idafh.eq.148) then
        ndar48=ndar
      else if (idafh.eq.149) then
        ndar49=ndar
      else
      endif
c
  100 continue
c
      nav = lenwrd()
c
c     map idafh onto internal fortran stream
c
      itmp = locat1(iunits,12,idafh)
      if (itmp.le.0) call caserr(
     +   'invalid fortran stream number in daopen')
      idafi = 20 + itmp
c
      nword = ndar + ndar + 2
      nword2 = nword / nav
c
      if (idafi.le.0.or.idafi.gt.maxlfn) then
         call caserr("daopen: file number idafi out of range")
      endif
      filatt = zedfil(idafi)
      inquire (file=filatt, exist=ok)
      if (ooprnt)  then
        call strtrm(filatt,lencha)
        write(6,*)'daopen: exist,iposun = ', idafi, 
     +   ok,iposun(idafi)
        write(6,*)'daopen: name = ', filatt(1:lencha)
      endif
      if (ok.or.iposun(idafi).ne.0) then
c
c     idafi exists or is already opened
c     restore header blocks
c
       call search (m1,idafi)
       if(ooprnt) write(6,*) 'header words = ', nword
       call readi(navhh,nword,m1,idafi)
       call icopy(2*ndar,ibuff,1,ioda,1)
       navh = navhh
c
      else
c
c     open data set through search
        call search (m1,idafi)
        do i=1,ndar
         ioda(1,i)=0
         ioda(2,i)=0
        enddo
c
c   write directory and the first available record address in the first
c   (file) block(s)
c
        navh = lensec(nword2) + 1
        navhh = navh
        call icopy(2*ndar,ioda,1,ibuff,1)
        if(ooprnt) write(6,*) 'header words = ', nword
        call wrt3i(navhh,nword,m1,idafi)
c
      endif
      if(ooprnt)  then
       write(6,*)'daopen: next block = ', navh
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine daclos(idafh,navh,iw)
      implicit real*8  (a-h,o-z),integer  (i-n)
      dimension iunits(12)
      data iunits /
     +    110, 120, 131, 141, 142, 143, 144, 145, 146, 147, 
     +    148, 149 /
c
c     map idafh onto internal fortran stream
c
      itmp = locat1(iunits,12,idafh)
      if (itmp.le.0) call caserr(
     +   'invalid fortran stream number in daopen')
      idafi = 20 + itmp
c
      close(unit=idafi,status='keep')
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine daread(idafh,ioda,x,nx,idar)
c
c  this routine reads from da-file
c
      implicit real*8  (a-h,o-z),integer  (i-n)
c
c
      integer isingl, nbits
      common /machin1/ isingl,nbits
c
      integer ixddaf
      integer ndar10, ndar20, ndar31, ndar41, ndar42, ndar43, ndar44
      integer ndar45, ndar46, ndar47, ndar48, ndar49
      common /dafindx/ ixddaf(8192),ndar10,ndar20,
     +                 ndar31,ndar41,ndar42,ndar43,ndar44,ndar45,
     +                 ndar46,ndar47,ndar48,ndar49
c
      integer lrec10, lrec20, lrec31, lrec41, lrec42, lrec43
      integer lrec44, lrec45, lrec46, lrec47, lrec48, lrec49
      common /dafrec/ lrec10,lrec20,lrec31,lrec41,lrec42,lrec43,
     +                lrec44,lrec45,lrec46,lrec47,lrec48,lrec49
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      integer mxioda
      parameter(mxioda=1000)
      dimension ioda(2,mxioda), buff(512)
      dimension x(*)
      dimension iunits(12)
      data iunits /
     +    110, 120, 131, 141, 142, 143, 144, 145, 146, 147, 
     +    148, 149 /
c
c     map idafh onto internal fortran stream
c
      itmp = locat1(iunits,12,idafh)
      if (itmp.le.0) call caserr(
     +   'invalid fortran stream number in daread')
      idafi = 20 + itmp
c
      nav = lenwrd()
c
      if (idar.le.0.or.idar.gt.mxioda) then
         write(*,*)"*** daread: idar   = ",idar
         write(*,*)"*** daread: mxioda = ",mxioda
         call caserr("daread: idar out of range")
      endif
c
      ldar =   nx
      iblock=ioda(1,idar)
      ierr = 0 
c
      nw = ioda(2,idar)
      if(nw.le.0) then
       write(6,*)'daread: invalid number of words in section', nw
       ierr = ierr + 1
      endif
      if (nx.le.0) then
       write(6,*)'daread: attempting to read invalid no. of words',
     + nx
       ierr = ierr + 1
      endif
      if (nw .lt.nx) then
       write(6,*)'daread: warning - more words read than present',
     +     idafi, iblock, nx
       ierr = ierr + 1
      endif
c
      if (ooprnt)
     + write(6,*) 'daread: idafh,iblock,nword = ', idafi, iblock, nx
c
      if (ierr.gt.0) call caserr(
     +'error detected in daread')
c
c
c     current rdedx cannot be used for this, for that routine
c     assumes that there is sufficient space in the receiving
c     array to hold the material as written. daread will 
c     sometimes recall fewer words than written, allocating space
c     for just the former.
c
      call search(iblock,idafi)
      j=1
 20   if(j.gt.ldar)go to 30
      call find(idafi)
      call get(buff,l)
      m = l
      if(j+l-1.gt.ldar) m = ldar-j +1
      if(m.gt.0)call dcopy(m,buff,1,x(j),1)
      j=j+l
      go to 20
30    return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dawrit(idafh,ioda,x,nx,idar,navh)
c
c  this routine writes to da-file
c
      implicit real*8  (a-h,o-z),integer  (i-n)
c
c
      integer isingl, nbits
      common /machin1/ isingl,nbits
c
      integer ixddaf
      integer ndar10, ndar20, ndar31, ndar41, ndar42, ndar43, ndar44
      integer ndar45, ndar46, ndar47, ndar48, ndar49
      common /dafindx/ ixddaf(8192),ndar10,ndar20,
     +                 ndar31,ndar41,ndar42,ndar43,ndar44,ndar45,
     +                 ndar46,ndar47,ndar48,ndar49
c
      integer lrec10, lrec20, lrec31, lrec41, lrec42, lrec43
      integer lrec44, lrec45, lrec46, lrec47, lrec48, lrec49
      common /dafrec/ lrec10,lrec20,lrec31,lrec41,lrec42,lrec43,
     +                lrec44,lrec45,lrec46,lrec47,lrec48,lrec49
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      common/junk2/navhh,navhhh,ibuff(1)
      integer mxioda
      parameter (mxioda=1000)
      dimension ioda(2,mxioda),x(*)
      dimension iunits(12)
      data iunits /
     +    110, 120, 131, 141, 142, 143, 144, 145, 146, 147, 
     +    148, 149 /
      data m1/1/
c
c     map idafh onto internal fortran stream for gamess i/o
c
      itmp = locat1(iunits,12,idafh)
      if (itmp.le.0) call caserr(
     +   'invalid fortran stream number in dawrit')
      idafi = 20 + itmp
c
      nav = lenwrd()
      ldar =   nx
      nrecs = lensec(ldar)
c
c ---    check whether record exists already (ioda(1,idar).ne.0)
c        if so, check whether it is large enough to contain new data.
c        if not,delete entry in -ioda- and write new record at -navh-.
c
      if (idar.le.0.or.idar.gt.mxioda) then
         write(*,*)"*** dawrit: idar   = ",idar
         write(*,*)"*** dawrit: mxioda = ",mxioda
         call caserr("dawrit: idar out of range")
      endif
      if((ioda(1,idar).ne.0).and.(ldar.le.ioda(2,idar))) then
        irec1=ioda(1,idar)
      else
        irec1=navh
      endif
c
      iblock=irec1
      if (ooprnt) 
     +   write(6,*) 'dawrit: idafi,iblock,nword = ', idafi, iblock, nx
      call wrt3(x,nx,iblock,idafi)
c
c ----- if new record is written, update directory on file
c
      if(iblock.eq.navh) then
        if(idafh.eq.110) ndar=ndar10
        if(idafh.eq.120) ndar=ndar20
        if(idafh.eq.131) ndar=ndar31
        if(idafh.eq.141) ndar=ndar41
        if(idafh.eq.142) ndar=ndar42
        if(idafh.eq.143) ndar=ndar43
        if(idafh.eq.144) ndar=ndar44
        if(idafh.eq.145) ndar=ndar45
        if(idafh.eq.146) ndar=ndar46
        if(idafh.eq.147) ndar=ndar47
        if(idafh.eq.148) ndar=ndar48
        if(idafh.eq.149) ndar=ndar49
c
        ioda(1,idar)=iblock
        ioda(2,idar)=nx
        navh=navh+nrecs
c
        nword = 2 * (ndar + 1)
        navhh = navh
        call icopy(2*ndar,ioda,1,ibuff,1)
        call wrt3i(navhh,nword,m1,idafi)
        if (ooprnt)
     +  write(6,*) 'dawrit: idafi, dirctory pointer = ', idafi, navhh
      endif
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function igms_int_cutoff()
      implicit none
c
c     Return the -log10 of the integral cutoff for storing integrals 
c     in GAMESS-UK. 
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
      igms_int_cutoff = icut
      return
      end

