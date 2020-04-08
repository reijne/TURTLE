c 
c  $Author: jmht $
c  $Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
c  $Locker:  $
c  $Revision: 6090 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util5.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   util5(ci) =
c ******************************************************
c ******************************************************
      subroutine pinver(ip,ipb,n)
c
c     invert permutation ip ==> ipb
c
      implicit real*8  (a-h,o-z),integer  (i-n)
      dimension ip(n),ipb(n)
c
      do 10 i=1,n
10    ipb(ip(i)) = i
c
      return
      end
      subroutine anti1(r,a,n)
      implicit real*8  (a-h,o-z),integer  (i-n)
      dimension r(*),a(n,*)
      m=1
      do 1 i=2,n
      im1=i-1
      do 1 j=1,im1
      r(m)=a(j,i)-a(i,j)
    1 m=m+1
      return
      end
      subroutine zbasgn(n,ibegin,iadd,iarr)
      implicit real*8  (a-h,o-z),integer  (i-n)
      integer iarr
      dimension iarr(*)
      incr=ibegin
      do 1 loop=1,n
      iarr(loop)=incr
    1 incr=incr+iadd
      return
      end
      subroutine smask
      implicit real*8  (a-h,m-z),integer    (i-l)
      integer mask,z1
      common/maskc/mask(64)
      data z1/z'8000000000000000'/
      mask(1)=z1
      do 1 i=2,64
    1 mask(i)=ior(mask(i-1),ishft(z1,1-i))
      return
      end
      subroutine isquar(ir,ia,mrowr,n)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ir(*),ia(*)
c... convert triangle ia to square ir
      iii=0
      jjj=1
      k=1
      do 40 i=1,n
      jj=jjj
      do 30 j=1,i
      ir(iii+j)=ia(k)
      ir(jj)=ia(k)
      k=k+1
30    jj=jj+mrowr
      iii=iii+mrowr
40    jjj=jjj+1
      return
      end
      subroutine sqsitr(r,s,t,nk)
      implicit real*8  (a-h,o-z), integer  (i-n)
      dimension s(*),t(*)
      dimension r(nk,*)
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
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
      do 1 i=1,nk
    1 r(i,i)=s(iky(i+1))
      m=1
      do 2 j=2,nk
      l=iky(j)
      j1=j-1
      do 2 i=1,j1
      r(i,j)=s(l+i)+t(m)
      r(j,i)=s(l+i)-t(m)
   2  m=m+1
      return
      end
      subroutine symm1a(r,a,nk)
      implicit real*8  (a-h,o-z),integer  (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=1,nk
      do 1 j=1,i
      r(m)=r(m)+a(j,i)+a(i,j)
    1 m=m+1
      return
      end
      subroutine anti1s(r,a,nk)
      implicit real*8  (a-h,o-z),integer (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=2,nk
      i1=i-1
      do 1 j=1,i1
      r(m)=r(m)+a(i,j)-a(j,i)
    1 m=m+1
      return
      end
      subroutine anti1a(r,a,nk)
      implicit real*8  (a-h,o-z),integer  (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=2,nk
      i1=i-1
      do 1 j=1,i1
      r(m)=r(m)-a(i,j)+a(j,i)
    1 m=m+1
      return
      end
      subroutine transe(qq,imap,local,maxtrs)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qq(*),imap(*),local(*)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      common/scra/ijklin(2,3412),index(3412)
c
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=((imap(ijklin(1,loop))-mloww)*maxtrs+
     +  local(ijklin(1,loop)))* lenbas+ijklin(2,loop)
      enddo
      call dsctr(nkk,gin,index,qq)
      return
      end
      subroutine transb(qq)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qq(*)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      common/scra/ijklin(2,3412),index(3412)
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=(ijklin(1,loop)-mloww)*lenbas + 
     +              ijklin(2,loop)
      enddo
      call dsctr(nkk,gin,index,qq)
      return
      end
      subroutine jsortp(map)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
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
      common/scra/ibuk(3400)
     * ,             index(3400)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      dimension map(*)
      mstack=1
      call igthr(nstack,map,index,itx)
      do 1 i=1,nstack
      ibuk(i)=(index(i)-mlow)*btri+1
 1    continue
      return
      end
      subroutine mcsrto (g,nijkl)
c
c...  subroutine to handle sortfile i/o. hopefully most machine dependnt
c     features are isolated here.
      implicit real*8  (a-h,o-z)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
      common /bufb/ nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      dimension g(nwb),nijkl(*)
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g,1,gout,1)
      call pack(gout(nsz341),lab1632,nijkl,nsz680)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
      nwb=0
      nwbuck(ibuck)=0
      return
      end
      subroutine mcsrti(mk,nk)
c
c...  subroutine to handle sortfile input. hopefully most machine
c     dependent features are isolated here.
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
      common /bufb  / nwbnwb,lnklnk,gout(5118)
      common /junk  /klin(2,3412)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      iblock = mk
      call rdbak(iblock)
      call stopbk
      mk=lnklnk
      nk=nwbnwb
      call unpack(gout(nsz341),lab1632,klin,nsz680)
      if(odebug(34)) then
       write(6,*)'mcsrti: gout, kin, lin'
       do int = 1,nwbnwb
        write(6,100)int,gout(int),klin(1,int),klin(2,int)
100     format(1x,i5,f20.9,2(1x,i7))
       enddo
      endif
c
      return
      end
      subroutine mult3 (q,h,r,p,nactiv,nbasis)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension q(*),h(*),r(*),p(*)
c
      integer nbnb, nanb, nanb2, nbasd2, nana, nananb
      common /mlngth/ nbnb,nanb,nanb2,nbasd2,nana,nananb
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
      do 1 loop=1,nbasis
      m=iky(loop+1)
    1 h(m)=h(m)*0.5d0
      nacnba=nactiv*nbasis
      if(nbasd2.ge.nactiv)go to 2
      call vclr(p,1,nacnba+nactiv*nactiv)
      call mxmm(q(nanb+1),nsa4,h,p,nactiv,nactiv,nbasis)
      call mxmb(p,1,nactiv,q,1,nbasis,p(nacnba+1),1,nactiv,
     * nactiv,nbasis,nactiv)
      call symm1(r,p(nacnba+1),nactiv)
      return
    2 continue
      call vclr(p,1,nacnba)
      n=1
      l=1
      do   4 loop=1,nactiv
      m=1
      do   5 moop=1,nbasis
      call daxpy(moop,q(n),h(m),1,p(l),1)
      n=n+1
   5  m=m+moop
   4  l=l+nbasis
      m=1
      n=1
      do 6 loop=1,nactiv
      l=1
      do   7 moop=1,loop
      temp1=ddot(nbasis,q(n),1,p(l),1)
      r(m)=temp1+ddot(nbasis,q(l),1,p(n),1)
      m=m+1
   7  l=l+nbasis
   6  n=n+nbasis
      return
      end
      subroutine mrgle(n,scalar,r)
      implicit real*8  (a-h,o-z), integer (i-n)
      dimension r(*)
      do 1 loop=1,n
      if(dabs(r(loop)).ge.scalar) go to 1
      r(loop)=dsign(scalar,r(loop))
    1 continue
      return
      end
      subroutine mrgge(n,scalar,r)
      implicit real*8  (a-h,o-z), integer (i-n)
      dimension r(*)
      do 1 loop=1,n
      if(dabs(r(loop)).lt.scalar) go to 1
      r(loop)=dsign(scalar,r(loop))
    1 continue
      return
      end
      subroutine stackr(g,nijkl)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*)
      dimension nijkl(*)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
      common/junk/nwbuck(1500),
     * itx(3400),ktx(3400),gtx(3400)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/scra/ibu(3400),itxktx(3400)
      call jsort1
    1 ibuck=ld340t(g,nijkl)
      if(ibuck)3,2,3
    3 call stopbk
      lnklnk=mark(ibuck)
      ib=ibase(ibuck)+1
      ibn=ibasen(ibuck)+1
      isz=nsz340
      call dcopy(isz,g(ib),1,gout,1)
      call pack(gout(nsz341),lab1632,nijkl(ibn),nsz680)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
      mstack=mstack+1
      if(mstack-nstack)1,1,2
   2  nstack=0
      return
      end
      function isort1(itri,ktri,gtri)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/integ(1360)
      common/blkin/gin(510),mword
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ipad
      dimension itri(*),ktri(*),gtri(*)
c
      call unpack(gin(num2e+1),lab816,integ,numlab)
      n=0
      iword=1
      do 1 loop=1,mword
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
      itx=iky(i)+j
      ktx=iky(k)+l
      gtx=gin(loop)
      if(mloww.ge.itx.or.mhi.lt.itx)go to 2
      n=n+1
      itri(n)=itx
      ktri(n)=ktx
      gtri(n)=gtx
    2 if(mloww.ge.ktx.or.mhi.lt.ktx.or.itx.eq.ktx)go to 1
      n=n+1
      itri(n)=ktx
      ktri(n)=itx
      gtri(n)=gtx
    1 iword=iword+4
*     write(6,121)i,j,k,l,itri(n),ktri(n),gtri(n)
*121  format('isort1: i,j,k,l,itri,ktri,G = ', 4i3,1x,2i6,2x,f10.5)
      isort1=n
      return
      end
      function isort2(itri,ktri,gtri)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/blkin/gin(510),nword
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/junke/maxt,ires,ipas,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ntri
      common/craypk/ij205(2,340)
      dimension itri(*),ktri(*),gtri(*)
      call unpack(gin(num2ep+1),lab1632,ij205,numlabp)
      n=0
      do 1 loop=1,nword
      ij=ij205(1,loop)
      if(mloww.ge.ij.or.mhi.lt.ij)go to 1
      n=n+1
      gtri(n)=gin(loop)
      itri(n)=ij
      ktri(n)=ij205(2,loop)
    1 continue
      isort2=n
      return
      end
      function isort3(imap,itri,ktri,gtri)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension itri(*),ktri(*),gtri(*),imap(*)
      common/blkin/gin(510),nword
      common/junke/maxt,ires,ipas,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ntri
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/ij205(2,340)
c
      call unpack(gin(num2ep+1),lab1632,ij205,numlabp)
      n=0
      do 1 loop=1,nword
      ij=ij205(1,loop)
      itrish=imap(ij)
      if(mloww.ge.itrish.or.mhi.lt.itrish)go to 1
      n=n+1
      gtri(n)=gin(loop)
      itri(n)=ij
      ktri(n)=ij205(2,loop)
    1 continue
      isort3=n
      return
      end
      subroutine jsort1
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
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
      common/scra/ibuk(3400)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      mstack=1
      do 1 i=1,nstack
    1 ibuk(i)=(itx(i)-mlow)*btri+1
      return
      end
      function ld340t(buf,ijklbuf)
      implicit real*8  (a-h,o-z),integer (i-n)
      common/bufb/nkk,mkk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
c     common/three/mark(1500),lnumb(1500)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/scra/ib(3400)
      common/junk/nwbuf(1500),itx(3400),ktx(3400),v(3400)
      dimension buf(*),ijklbuf(*)
      mstak=mstack
    2 ld340t=ib(mstak)
      vv=v(mstak)
      nw=nwbuf(ld340t)+1
      ln=ibase(ld340t)+nw
      buf(ln)=vv
      ijklbuf(ln+ln-1)=itx(mstak)
      ijklbuf(ln+ln  )=ktx(mstak)
*     write(6,121)ln,buf(ln),ijklbuf(ln+ln-1),ijklbuf(ln+ln  )
*121  format('ld340t: ',1x,i5,f15.5,2x,2i10)
      if(nw.ge.nsz340)go to 4
      nwbuf(ld340t)=nw
      mstak=mstak+1
      if(mstak.le.nstack)go to 2
      ld340t=0
      return
    4 nwbuf(ld340t)=0
      mstack=mstak
      return
      end
      subroutine transc(qq)
      implicit real*8  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      common/junke/maxt(8),mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
      common/scra/ijklin(2,3412),index(3412)
      dimension qq(*)
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=(ijklin(1,loop)-mloww)*lenb4 + 
     +              ijklin(2,loop)
      enddo
      call dsctr(nkk,gin,index,qq)
      return
      end
      function locate(ref,nref,nw,pat)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer ref,pat
      dimension ref(nw,*)
      do 1 loop=1,nref
      if(ref(1,loop).eq.pat)goto2
    1 continue
      locate=0
      return
    2 locate=loop
      return
      end
      subroutine upak8w(gijkl,iiii,mapper)
      implicit real*8  (a-h,o-z),integer  (i-n)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/blkin/gout(510),nwb
      common/craypk/i205(4,340)
      dimension icray(1360)
      dimension mapper(*),gijkl(*),iiii(*)
c
      call unpack(gijkl,lab816,i205,numlab)
      call igthr(nwb*4,mapper,icray,i205)
      int4=1
      do 1 iw=1,nwb
      j=icray(int4 )
      i=icray(int4+1)
      l=icray(int4+2)
      k=icray(int4+3)
      int4=int4+4
      if(i.ge.j)goto 2
      m=i
      i=j
      j=m
2     if(k.ge.l)goto 3
      m=k
      k=l
      l=m
3      if((i*256+j).ge.(k*256+l))goto 4
      m=i
      i=k
      k=m
      m=j
      j=l
      l=m
  4   i205(1,iw)=j
      i205(2,iw)=i
      i205(3,iw)=l
      i205(4,iw)=k
 1    continue
      end
      subroutine upak8s(gijkl,iiii)
      implicit real*8  (a-h,o-z),integer  (i-n)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
      common/junk2/index(4,3412)
      common/bufb/nint,nlink,mij(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
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
      dimension gijkl(*),iiii(*)
c     write(6,*) 'upak8s'
c     j = nsz340/2 + nsz341
c     write(6,99996) (gijkl(i),i=nsz341,j)
c99996 format(1x,4(1x,z16))
      call setsto(13648,0,index)
      call unpack(gijkl(nsz341),lab816,index,nsz340*4)

      do 1 i=1,(nint+1)/2
      i205(2*i-1)=index(2,2*i-1)
      j205(2*i-1)=index(1,2*i-1)
      k205(2*i-1)=index(4,2*i-1)
      l205(2*i-1)=index(3,2*i-1)
      i205(2*i)=index(2,2*i)
      j205(2*i)=index(1,2*i)
      k205(2*i)=index(4,2*i)
   1  l205(2*i)=index(3,2*i)
      return
      end
      subroutine vvdv(n,r,a,b)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      dimension a(*),b(*),r(*)
      do 1 loop=1,n
    1 r(loop)=a(loop)/b(loop)
      return
      end
      subroutine gtrian(n,scalar,r,a,b)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),b(*)
      if(n.le.0)go to 2
      do 1 loop=1,n
    1 r(loop)=scalar*b(loop)-a(loop)
    2 return
      end
      subroutine wto (buf,lbuf)
      character*(*) buf
      character*80 obuf
      obuf=buf(1:lbuf)
      call sptchk(obuf)
      return
      end
      subroutine fixorb(v,ndim,nrow,ncol)
      implicit none
c
c     This routine canonicalises the sign of the vectors.
c     This is useful in comparing vectors from different runs for
c     debugging purposes.
c     The canonicalisation is done such that the sum of all coefficients
c     in every vector is positive.
c
      integer ndim      ! the leading dimension of the vectors
      integer nrow      ! the length of the vectors
      integer ncol      ! the number of vectors
      real*8 v(ndim,ncol) ! the vectors
c
      integer i
      real*8 sum, dsum
c
      do i = 1, ncol
         sum = dsum(nrow,v(1,i),1)
         if (sum.lt.0.0d0) then
            call dscal(nrow,-1.0d0,v(1,i),1)
         endif
      enddo
      end
      subroutine dchksm(n,x,ix,string)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension x(*)
      character*(*) string
c
c real check sum routine for debug
c
      sum = 0.0d0
      asum =0.0d0
      if (n.gt.0) then
         sum = dsum(n,x,ix)
         asum = dasum(n,x,ix)
      endif
      write(iwr,1) string,n,ix,sum,asum
 1    format(1x,a,2i11,2(2x,f23.15))
      return
      end
      subroutine pchksm(n,k,ik,string)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension k(2,*),ix(2),io(2),ia(2)
      dimension iix(2),ioo(2),iaa(2)
      character*(*) string

c
c this is the Sun format - replace code as required
c new machines
c
      data iix / z'00000000',z'00000000'/
      data ioo / z'00000000',z'00000000'/
      data iaa / z'ffffffff',z'ffffffff'/
c
c packed integer check sum routine for debug
c

      do 20 loop=1,2
      ix(loop)=iix(loop)
      io(loop)=ioo(loop)
 20   ia(loop)=iaa(loop)
      if (n.gt.0) then
         do 10 i = 1,1+(n-1)*ik,ik
            ix(1) = ieor(ix(1),k(1,i))
            ix(2) = ieor(ix(2),k(2,i))
            io(1) = or(io(1),k(1,i))
            io(2) = or(io(2),k(2,i))
            ia(1) = and(ia(1),k(1,i))
            ia(2) = and(ia(2),k(2,i))
 10      continue
c         write(6,2) (k(i),i=1,n)
c2        format(4(1x,z16))
      endif
c
      write(iwr,1) string,n,ik,ix,io,ia
 1    format(1x,a,2i6,3(2x,2z8))
      return
      end
      subroutine ichksm(n,k,ik,string)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension k(*)
      character*(*) string
c
c integer check sum routine for debug
c
      isum = 0
      iasum = 0
      if (n.gt.0) then
         do 10 i = 1,1+(n-1)*ik,ik
            isum = isum + k(i)
            iasum = iasum + iabs(k(i))
 10      continue
c         write(6,2) (k(i),i=1,n)
c2        format(7(1x,i10))
      endif
c
      write(iwr,1) string,n,ik,isum,iasum
 1    format(1x,a,2i11,2(2x,i11))
      return
      end
      subroutine wrtsor( text,i205,j205,k205,l205,ggg,nwb)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character *(*) text
      dimension i205(*),j205(*),k205(*),l205(*),ggg(*)
      dimension f(4),ia(4),ja(4),ka(4),la(4)
      write(iwr,1) text,nwb
1     format(/1x,a6,' sorted integrals = ',i5)
      iseq = 0
      n = 0
      do loop=1,nwb
      n = n + 1
      f(n)  = ggg(loop)
      ia(n) = i205(loop)
      ja(n) = j205(loop)
      ka(n) = k205(loop)
      la(n) = l205(loop)
      if(n.eq.4) then
       iseq = iseq + 1
       write(iwr,200)iseq,(ia(k),ja(k),ka(k),la(k),
     +               f(k),k=1,n)
       n = 0
      endif
      enddo
      if(n.ne.0) then
       iseq=iseq+1
       write(iwr,200)iseq,(ia(k),ja(k),ka(k),la(k),
     +               f(k),k=1,n)
      endif
      return
200   format(i5,2x,4(4i4,f11.6))
      end
      subroutine dbgvec(text,array,n)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character *(*) text
      dimension array(*)
      write(iwr,200)text,n,(array(i),i=1,n)
200   format(1x,a8,i8/(8f12.5))
      return
      end
      subroutine dbgvecv(text,array,n)
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character *(*) text
      dimension array(*)
      dimension f(6),ia(6)
      data thresh /1.0d-5/
c
      write(iwr,200)text,n
      iseq = 0
      m = 0
      do loop=1,n
       if(abs(array(loop)).ge.thresh) then
        m = m + 1
        f(m)  = array(loop)
        ia(m) = loop
        if(m.eq.6) then
        iseq = iseq + 1
        write(iwr,201) iseq,(ia(k),f(k),k=1,m)
        m = 0
        endif
       endif
      enddo
      if(m.ne.0) then
       iseq=iseq+1
       write(iwr,201)iseq,(ia(k),f(k),k=1,m)
      endif
      return
200   format(/1x,'**** ',a10,2x,i8/)
201   format(i5,2x,6(i6,f12.5))
      end
      subroutine wtfort(q,nword,irec,ifort)
c
c     routines for input/output on fortran streams
c     if called with nword.eq.0 no i/o takes place, but the file is
c     positioned ready to read or write at record irec.
c
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension q(*)
      rewind ifort
      assign 10 to label
      if(irec-1)10,10,30
10    if(nword.gt.0)call wtfors(q,nword,ifort)
      return
c
      entry rdfort(q,nword,irec,ifort)
      rewind ifort
      assign 20 to label
      if(irec-1)20,20,30
20    if(nword.gt.0)call rdfors(q,nword,ifort)
      return
c
30    irecm=irec-1
      do 40 irec=1,irecm
40    read(ifort,err=60,end=50)
      goto label,(10,20)
60    write(iwr,70)ifort,irec
      call caserr('error on attempting to read fortran dataset.')
      return
50    write(iwr,70)ifort,irec
      call caserr('end of file reached on reading fortran dataset.')
      return
70    format(' fault on fortran stream',i3,' at record',i7)
      end
      subroutine wtfors(q,nword,ifort)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(nword)
      if(nword.gt.0)write(ifort)q
      return
      end
      subroutine rdfors(q,nword,ifort)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(nword)
      if(nword.gt.0)read(ifort,end=60,err=50)q
      return
60    call caserr('end of file reached on reading fortran dataset.')
50    call caserr(' error on reading fortran dataset.')
      return
      end
      subroutine ver_util5(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util5.m,v $
     +     "/
      data revision /"$Revision: 6090 $"/
      data date /"$Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
