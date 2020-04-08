c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/machci.m,v $
c  $State: Exp $
c  
      subroutine closbf2(idel)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *1024 zsort,blank
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz2/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      if(nsmax.ge.0) then
       call closecc(nsort)
       if(ooprnt)write(iwr,1)
       nsmax = -1
       call gtnv('psort',zsort)
       blank = ' '
       if(zsort.eq.blank) zsort = 'psort'
       call strtrm(zsort,length)
       call delcc(zsort,length,ierrio)
       if(ierrio.ne.0)call ioerr('delete',0,zsort)
       if(ooprnt)write(iwr,2)
      endif
1     format(/1x,'*** p-sortfile closed')
2      format(/1x,'*** p-sortfile deleted')
      return
      end
      subroutine closbf3
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/blksi3/nsz(12),nsmax,nsort
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      if(nsmax.gt.0)
     *call closecc(nsort)
      if(ooprnt)write(iwr,1)
1     format(/
     * ' tablefile closed')
      return
      end
      subroutine closbf(idel)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *1024 zsort,blank
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      if(nsmax.ge.0) then
       call closecc(nsort)
       nsmax = -1
       call gtnv('sort',zsort)
       blank = ' '
       if(zsort.eq.blank) zsort = 'sort'
       call strtrm(zsort,length)
       call delcc(zsort,length,ierrio)
       if(ierrio.ne.0)call ioerr('delete',0,zsort)
       if(ooprnt) then
        write(iwr,1)
1       format(/1x,'*** sortfile closed')
        write(iwr,2)
2       format(/1x,'*** sortfile deleted')
       endif
      endif
      return
      end
      subroutine rdbak(iblock)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
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
      common/bufb/buffer(5120)
      call srchcc(nsort,(iblock+1),ierrio)
      if(ierrio.ne.0)call ioerr('search',0,'sort')
      call getccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('read',0,'sort')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine rdbak2(iblock)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz2/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/bufc/buffer(512)
      call srchcc(nsort,(iblock+1),ierrio)
      if(ierrio.ne.0)call ioerr('search',0,'psort')
      call getccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('read',0,'psort')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine rdbak3(iblock)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/blksi3/nsz,nsz512,nszz(8),nsstat
     *,nslen,nsmax,nsort
      common/bufd/buffer(512)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      if (ooprnt) then
       ibl = iblock + 1
       write(iwr,10) nsort, ibl
 10    format(1x,'**** get from table:', i3,
     +           ' at block',i5)
      endif
      call srchcc(nsort,(iblock+1),ierrio)
      if(ierrio.ne.0)call ioerr('search',0,'table')
      call getccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('read',0,'table')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine setbfa
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *1024 zsort,blank
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
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
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
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/disc/ispz(3),irep
c
c     sort file settings
c
      call setsrtp(nsz,o255i,
     +  nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     +  nsz510,nsz341,nsz342,nsz680)
      nsstat=0
c
      call gtnv('sort',zsort)
      blank = ' '
      if(zsort.eq.blank) zsort = 'sort'
      call strtrm(zsort,length)
      irep = 0
      if (nsmax.lt.0) then
       call opencc(zsort,length,nsort,ierrio)
       if(ierrio.ne.0)call ioerr('open',0,zsort)
       if(irep.ne.0)
     * call caserr(
     * 'error creating/opening sort file')
       if(ooprnt)write(iwr,1)zsort(1:length),nsz
       nsmax = 0
      endif
c
      return
1     format(/
     * ' sortfile ',a,' allocated '/
     * '          : blockfactor =',i3, ' block(s) ')
      end
      subroutine setbfb
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *1024 zsort,blank
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
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
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
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz2/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/disc/ispz(3),irep
c
c     sort file settings
c
      call setsrtp(nsz,o255i,
     +  nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     +  nsz510,nsz341,nsz342,nsz680)
      nsstat=0
c
      call gtnv('psort',zsort)
      blank = ' '
      if(zsort.eq.blank) zsort = 'psort'
      call strtrm(zsort,length)
      irep=0
      if (nsmax.lt.0) then
       call opencc(zsort,length,nsort,ierrio)
       if(ierrio.ne.0)call ioerr('open',0,zsort)
       if(irep.ne.0)
     * call caserr(
     * 'error creating/opening p-sort file')
       if(ooprnt)write(iwr,1)zsort(1:length),nsz
       nsmax = 0
      endif
c
      return
1     format(/
     * ' p-sortfile ',a,' allocated '/
     * '          : blockfactor =',i3, ' block(s) ')
      end
      subroutine setbfc
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *1024 zsort,blank
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
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      common/blksi3/nsz,nsz512,nsz340,nsz170,nsz85,
     * nszij,nszkl,nsz510,nsz341,nsz342,nsstat,nslen,nsmax,nsort
c
c     table file settings
c
      nsz512=nsz*512
      nsz170=(nsz512-2)/3
      nsz85 =nsz170/2
      nszkl=nsz170+1
      nsz340=nsz170+nsz170
      nsz341=nsz340+1
      nsz510=nsz340+nsz170
      nszij=nsz340+1
      nsz342=nsz341+nsz85
      nsstat=0
c
      call gtnv('table',zsort)
      blank = ' '
      if(zsort.eq.blank) zsort = 'table'
      call strtrm(zsort,length)
      irep=0
      call opencc(zsort,length,nsort,ierrio)
      if(ierrio.ne.0)call ioerr('open',0,zsort)
      if(irep.ne.0)
     *call caserr(
     *'error creating/opening table-ci file')
      if(ooprnt)write(iwr,1)zsort(1:length),nsz
1     format(/
     * ' table-ci file ',a,' allocated '/
     * '          : blockfactor =',i3, ' block(s) ')
c
      return
      end
      subroutine stopbk2
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
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
      common/disc/ispz(3),irep
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz2/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak2/btri,mlow(2),iblock
      if(nsstat.eq.0)go to 2
      irep = 0
      nsstat=0
 2    if(nsmax.lt.iblock)nsmax=iblock
      return
      end
      subroutine stopbk3
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
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
      common/disc/ispz(3),irep
      common/blksi3/nsz(10),nsstat,
     *nslen,nsmax,nsort
      common/stak3/btri,mlow(2),iblock
      if(nsstat.eq.0)go to 2
      irep = 0
      nsstat=0
 2    if(nsmax.lt.nslen)nsmax=nslen
      return
      end
      subroutine stopbk
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
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
      common/disc/ispz(3),irep
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
      if(nsstat.eq.0)go to 2
      irep = 0
      nsstat=0
 2    if(nsmax.lt.iblock)nsmax=iblock
      return
      end
      subroutine sttout2
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4(y)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz2/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/bufc/buffer(1)
      common/stak2/btri,mlow,nstack,iblock
      if (iblock.ne.nslen) then
         call srchcc(nsort,(iblock+1),ierrio)
         if(ierrio.ne.0)call ioerr('search',0,'psort')
      endif
      call putccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('write',0,'psort')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine sttout3(jrec)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4(y)
      common/blksi3/nsz,nsz512,nsz340(8),nsstat
     * ,nslen,nsmax,nsort
      common/bufd/buffer(1)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      common/stak3/btri,mlow,nstack,iblock
      iblock=jrec
      if (ooprnt) then
       ibl = iblock + 1
       write(iwr,10) nsort, ibl
 10    format(1x,'**** write to table:', i3,
     +           ' at block',i5)
      endif
      if (iblock.ne.nslen) then
         call srchcc(nsort,(iblock+1),ierrio)
         if(ierrio.ne.0)call ioerr('search',0,'table')
      endif
      call putccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('write',0,'table')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine sttout
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4(y)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/bufb/buffer(1)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
      if (iblock.ne.nslen) then
         call srchcc(nsort,(iblock+1),ierrio)
         if(ierrio.ne.0)call ioerr('search',0,'sort')
      endif
      call putccn(nsort,buffer,nsz,ierrio)
      if(ierrio.ne.0)call ioerr('write',0,'sort')
      nsstat=nsstat+1
      nslen=iblock+nsz
      return
      end
      subroutine setsrtp(nsz,o255i,
     +  nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     +  nsz510,nsz341,nsz342,nsz680)
c
c     sort file settings 
c
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      logical o255i
c
      nsz512=nsz*512
      if (o255i) then
       nsz170=(nsz512-2)/3
       nsz85 =nsz170/2
       nszkl=nsz170+1
       nsz340=nsz170+nsz170
       nsz680=nsz340+nsz340
       nsz341=nsz340+1
       nsz510=nsz340+nsz170
       nszij=nsz340+1
       nsz342=nsz341+nsz85
      else    
       nsz340=(nsz512-2)/2
       nsz680=nsz340+nsz340
       nsz341=nsz340+1
      endif   
c
      return  
      end 
      subroutine ver_machci(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/machci.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
