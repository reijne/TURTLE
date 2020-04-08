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
      REAL a(n),b(n)
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
      REAL function prpdot(a,b,n)
      implicit none
c
c     prpdot -> tracep [util1.m] : trace of product of 2 symmetric
c                                  matrices
c
      integer n
      REAL a(n), b(n)
      REAL tracep
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
      REAL d(n)
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
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/utilc)
INCLUDE(../m4/common/disc)
INCLUDE(../m4/common/discc)
INCLUDE(../drf/comdrf/sizesrf)
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
INCLUDE(../drf/comdrf/darw)
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
      implicit REAL  (a-h,o-z),integer  (i-n)
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
      implicit REAL  (a-h,o-z),integer  (i-n)
c
INCLUDE(../drf/comdrf/darw)
INCLUDE(../m4/common/utilc)
_IF(ga)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/disc)
_ENDIF
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
_IF(ga)
      if( iwhere( idafi ) .eq. 6 .or. iwhere(idafi) .eq. 7 ) then 
         call rdedx_ga( x, ldar, iblock, idafi)
      else
_ENDIF
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
_IF(ga)
      end if
_ENDIF
      end
c
c-----------------------------------------------------------------------
c
      subroutine dawrit(idafh,ioda,x,nx,idar,navh)
c
c  this routine writes to da-file
c
      implicit REAL  (a-h,o-z),integer  (i-n)
c
INCLUDE(../drf/comdrf/darw)
INCLUDE(../m4/common/utilc)
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
INCLUDE(common/restar)
      igms_int_cutoff = icut
      return
      end

