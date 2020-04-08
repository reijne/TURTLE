c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/intege.m,v $
c  $State: Exp $
c  
      subroutine cangt1(buf,ibuf,lbuf,inbrel,map,nzero,q,ifort)
      implicit real*8  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
c
c     reads in canonical list, version 1 - produces
c     map and condensed triangle
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real read from buf/ibuf
c     inbint = position of last integer read from buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      if (inbrel.eq.-1) then
         call reads(buf,lbuf,ifort)
         inbrel = 0
      end if
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     get value of nzero in buffer
c
 20   nzero = ibuf(lenrel(inbrel)+1)
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be read
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            map(nfirst+i) = ibuf(inbint+i)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            map(nfirst+i) = ibuf(inbint+i)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            q(nfirst+i) = buf(inbrel+i)
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            q(nfirst+i) = buf(inbrel+i)
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call reads(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60) , iret
c
 100  return
      end
      subroutine cangt2(buf,ibuf,lbuf,inbrel,map,nzero,q,ifort)
      implicit real*8  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
c
c     reads in canonical list
c     version 2 - produces map and expanded list of triangle
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real read from buf/ibuf
c     inbint = position of last integer read from buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      if (inbrel.eq.-1) then
         call reads(buf,lbuf,ifort)
         inbrel = 0
      end if
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     get value of nzero in buffer
c
 20   nzero = ibuf(lenrel(inbrel)+1)
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be read
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            map(nfirst+i) = ibuf(inbint+i)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            map(nfirst+i) = ibuf(inbint+i)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            q(map(nfirst+i)) = buf(inbrel+i)
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            q(map(nfirst+i)) = buf(inbrel+i)
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call reads(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60) , iret
c
 100  return
      end
      subroutine canonc(q,maxq,num,iblki,ifili,ifort,ltin)
c
c-------------------------------------------------------------
c     put list of integrals into canonical order
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      logical ltin
      dimension q(maxq)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/craypk/labs(1360)
c
c     ibl5 = no of integrals in block of sortfile
c     only machine dependent feature should be structure of /bufa/
c
      iwrt = 6
      ibl5 = nsz340
      nav = lenwrd()
      iilen = nsz340*nav/2
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = num*(num+1)/2
      n2 = nij
c
      maxt = (maxq-n2-lenint(n2))/nij
      nword = (maxq/(1+lenrel(1)))*lenrel(1)
c
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
c
c     nbuck is the number of buckets required
c
      nbuck = (nij/maxt) + 1
      nadd = min(maxt,nij)
      maxa = nbuck*(ibl5+ibl5i)
c
c
      if (nbuck.gt.maxb) then
         write (iwrt,6010) maxq , maxa
         call caserr('stop')
      end if
c
c     read through original file producing sorted file
c
      call vclr(q,1,maxa)
      call setsto(numlab,0,labs)
      call setbfa
      call canrd1(q,q,maxa,iblki,ifili,ltin)
c
c     read through the sort file to give final result
c
      maxqq = nij*nadd
      lbuf = nsz*512
      i1 = 1
      i2 = i1 + lbuf
      i3 = i2 + lenint(n2)
      call canwt1(q(i1),lbuf,q(i2),q(i3),maxqq,ifort)
c
      call closbf
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine canput(buf,ibuf,lbuf,inbrel,map,nzero,q,last,ifort)
      implicit real*8  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
      logical last
c
c     writes out canonical list
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real written into buf/ibuf
c     inbint = position of last integer written into buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     put value of nzero in buffer
c
 20   ibuf(lenrel(inbrel)+1) = nzero
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be written
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            ibuf(inbint+i) = map(i+nfirst)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            ibuf(inbint+i) = map(i+nfirst)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            buf(inbrel+i) = q(map(i+nfirst))
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            buf(inbrel+i) = q(map(i+nfirst))
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call wrt3s(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60,110) , iret
c
 100  if (last) then
         iret = 4
         go to 90
      end if
 110  return
      end
      subroutine canrd1(a,ia,maxa,iblki,ifili,ltin)
c
c     does the sorting part to get coulomb matrices
c     lower triangles only are produced
c     adaption for use in sorting to canonical order ,
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      logical ltin
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
c     note arrays a and ia actually overlap
c
      dimension a(maxa),ia(maxa)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/blkin/gin(510),nint
      common/stak/btri,mlow,nstack,iblock
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/bufb/nkk1,mkk1,g(1)
      common/craypk/labs(1360)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension iaad(340),ibbu(340)
      data lastb/999999/
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i)
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      call vclr(g,1,nsz340+nsz170)
c
      iblock = 0
c     ninb = no of elements in bucket(coreload)
      ninb = nadd*nij
c
      call search(iblki,ifili)
 30   call find(ifili)
      call get(gin,nw)
      if (nw.eq.0) then
c
c     empty anything remaining in buckets
c
         do 40 ibuck = 1 , nbuck
            nwb = nwbuck(ibuck)
            if (nwb.ne.0) then
               call stopbk
               mkk1 = mark(ibuck)
               nkk1 = nwb
               call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
               call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 40      continue
c
c
         call stopbk
         return
      else
         if (ltin) then
            call unpack(gin(num2e+1),lab816,labs,numlab)
            do 50 int = 1 , nint
               n4 = int + int + int + int
               ij = iky(labs(n4-2)) + labs(n4-3)
               kl = iky(labs(n4  )) + labs(n4-1)
c
               iaddr = (ij-1)*nij + kl
c
c     iaddr is address of integral in final sequence
c
               ibuck = (iaddr-1)/ninb
               iaad(int) = iaddr - ninb*ibuck
               ibbu(int) = ibuck + 1
 50         continue
         else
            call unpack(gin(num2e+1),16,labs,680)
            do 60 int = 1 , nint
               n4 = int + int
c
               iaddr = (labs(n4-1)-1)*nij + labs(n4)
c
c     iaddr is address of integral in final sequence
c
               ibuck = (iaddr-1)/ninb
               iaad(int) = iaddr - ninb*ibuck
               ibbu(int) = ibuck + 1
 60         continue
         end if
c
c     element goes in bucket ibuck with modified address
c
         do 70 int = 1 , nint
            ibuck = ibbu(int)
            nwb = nwbuck(ibuck) + 1
            a(ibase(ibuck)+nwb) = gin(int)
            ia(ibasen(ibuck)+nwb) = iaad(int)
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - empty
c
               call stopbk
               mkk1 = mark(ibuck)
               nkk1 = nwb
               call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
               call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
c
 70      continue
         go to 30
      end if
      end
      subroutine canwt1(buf,lbuf,ibuff,q,maxqq,ifort)
c----------------------------------------------------------
c     this reads back down the back-chained sort file
c     produced by canrd1 to give a final file
c     containing the coulomb matrices arranged
c     sequentially
c----------------------------------------------------------
      implicit real*8  (a-h,o-z)
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
      logical last
      dimension q(maxqq),buf(lbuf),ibuff(n2)
      common/sortpk/labin(1)
      common/stak/btri,mlow,nstack,iblock
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/bufb/nkk,mkk,g(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call rewedz(ifort)
c
      min = 1
      max = nadd
      inbrel = 0
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxqq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     triangles min thru max are in core - clear them out
c
            itri = 1
            do 30 n = min , max
c
c       get map of non-zero elements
c
               call dlstne(n,q(itri),1,0.0d0,nzero,ibuff)
               last = n.eq.nij
               call canput(buf,buf,lbuf,inbrel,ibuff,nzero,q(itri),last
     +  ,ifort)
               itri = itri + nij
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nij) max = nij
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
            call unpack(g(nsz341),32,labin,ibl5)
            call dsctr(nkk,g,labin,q)
c
            go to 20
         end if
 40   continue
      return
      end
      subroutine ver_intege(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/intege.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
