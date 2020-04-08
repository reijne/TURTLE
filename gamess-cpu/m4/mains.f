c 
c  $Author: jmht $
c  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
c  $Locker:  $
c  $Revision: 6317 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mains.m,v $
c  $State: Exp $
c  
c
c use unix memory allocation from C
c
      program gamess
      implicit none
      integer ierr
*
* for MPI code in parallel.m (IJB)
*
*  Maximum Number of processors we can run on
*
      Integer    max_processors
      Parameter( max_processors = 1024 )
*
* The workers communicator - why cant this be called something
* like shop_steward ?
*
* For newscf_f90 need a BLACS context, possibly two
*
      Integer MPI_COMM_WORKERS, blacs_context, blacs_uhf_context
      Integer MPI_MAXABS ! handle for MPI_ALLREDUCE function
      Logical is_farm
      Common /mpi_hack_job/ MPI_COMM_WORKERS, blacs_context, 
     +	blacs_uhf_context, is_farm, MPI_MAXABS

*
* The global GAMESS communicator - will be MPI_COMM_WORLD unless we
* are running under task farming
*
      Integer MPI_COMM_GAMESS
      Common /mpi_comms/ MPI_COMM_GAMESS

      Save   /mpi_hack_job/
      Save   /mpi_comms/
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)


      integer ioff
      integer idum
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer iret
      logical opg_root
      integer ipg_nodeid



      call pg_begin
      call reset_time_periods




       iret=0
       call initj(nmaxly)

c
c print out Node information
c placed here so that GA data server processes don't get listed
c
      call pg_pproc
c
c
c If we are not using MA tools for memory, allocate
c
      call mallocc(qq(1),qq(2),nmaxly,ioff,memhandle)
      if(ioff .eq. 0)call caserr('memory problem')

      maxq=nmaxly


      call gmem_set(qq(1 + ioff))
      call init_segm
      call driver(qq(1 + ioff))

      call gmem_highwater(IGMEM_QUIET)
      call gmem_summarize('mains.m','gamess','after_driver',IGMEM_DEBUG)


c
c flag unconverged scf
c     write(6,*)'irest', irest
c
      if (irest .ne. 0) iret = 1

      call clenup2
      call pg_end(iret)
      write(iwr,*)'pg_end should not return'
      end
c
      integer function icoresize()
c
c default for workstation version
      icoresize=4000000
      end
c
      subroutine clenup2
      implicit none
      integer code
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer idum
c
c flush standard out
      call flushn(6)
c release core memory
c note problem on C200
      call freec(qq(1), memhandle, idum)
      end
c




**==io.f
      block data iodata
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),llll(maxlfn),iop(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),istatf(maxfrt),
     *  lrec(maxfrt),nrec(maxfrt),oform(maxfrt)
c
c /bufa/
c
c
      common/bufa/aaaaa(8190)
      common/blksiz/nsz,nsz512(8),nsznav,nsstat,
     *  nslen,nsmax
      common/blksiz2/ns2,ns2512(10),npznav,npstat,
     * n2len,n2max
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     * n3len,n3max
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      data maxb,iblkla/9999999,1/
      data jrec/-1/
      data jwidth/132/,oswit/.false./,oterm/.false./,oflush/.true./
      data nerr/999/
      data nline,noline/0,0/
      data iwidth/132/
      end
**==all.f
      block data all_data
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      common/nottwi/obeen,obeen2,obeen3,obeen4
      parameter (ncmax=65)
      parameter (mcen33 = 9 * maxat)
      parameter (ntr  = 1900)
      common/junk2/pcap(13,ntr)
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
c
      integer natom, nres, nbonh, nbona, ntheth, ntheta, nphih
      integer nphia, nnb, ntypes, nrcp, nconp, mema, memb, maxnb
      integer mbona, mtheta, mphia, natc, npro, numbnd, numang
      integer numphi, natyp, nphb, imodf, jmodf
      common /modfov/ natom,nres,nbonh,nbona,ntheth,ntheta,nphih,
     +                nphia,nnb,ntypes,nrcp,nconp,mema,memb,maxnb,
     +                mbona,mtheta,mphia,natc,npro,numbnd,numang,
     +                numphi,natyp,nphb,imodf(25),jmodf(50)
c
c
      common/bufb/rrbufb(1500+mcen33),iibufb(200+mxgrps)
      common/maskc/rmaskc(64)
      character *8 zdate,ztime,zaccno,zanam
      common/jinfo/zdate,ztime,zaccno,zanam
c     common/cslosc/closc(7),icslo(1084)
      common/craypk/icall(14+mxcsf*mxnshl+mxcsf+mxcsf),isymm(mxroot),
     *cofff(mxcsf*mxroot),t0,t1
c = 256 bfn / 100 centres
      common/blkin/ctxxx(maxat*12+3),itcxxx(4)
      common/junk/ccxyz(36776),iicxyz(101)
c     common/linkmc/
c = 256  bfn / 1000 centres
c
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
      common /rcombd/ win(36),pxin(36),pyin(36),pzin(36)
c     dimension of ydd2
      common/direc/zrunt(50),yrunt(50),ydd1(100),ydd2(120),
     *             zdd3(70) ,zdd4(50)
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
      common/hermit/h(45)
      common/wermit/w(45)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
c
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      real*8 slder, alimin
      common/trntim/slder,alimin
c
c (common/timanb)
c
      real*8 timbeg,  tmlast,  timsec2
      real*8 wtimbeg, wtmlast, wtimsec
      real*8 cpwa
      integer nodel
      common/timanb/timbeg,tmlast,timsec2(14),
     *wtimbeg,wtmlast,wtimsec(14),cpwa,nodel
c
c
      common/junkc/ylabel(26),ztype(35),zlist(300),znumb(100),zhead(90),
     * zheadr(20),zunit(33),ztittt(100),ztitl(120),zspin(30),
     * zsta,zhz,zgauss,zcont(10),zconta(maxat),zmul(2,maxorb)
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
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/symtry/invt(206)
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
      common/foropt/vibsiz,nvib
      common/vibrtn/nxyz,iatom,icoord,ivib
      common/transf/psmal,qsmal,rsmal,pnew,qnew,rnew,pp,qp,rp
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,p0,q0,r0
      common/czmat/ianz(maxnz,5),bla(maxnz,3),lbla(maxnz,3),nz,nnvar
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
      character *8 zvar
      common/csubch/zvar(maxvar)
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      common/shlnos/qq4,litt(21)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      real*8 dij, dkl
      common /denss/ dij(225),dkl(225)
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      common/setint/ikn(24),bp01(16),ij12(4)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinf/ ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      common/pkfil/opkfl(3)
c
      real*8 p, q, r
      common /bshell/ p(mxshel),q(mxshel),r(mxshel)
c
      common/root/suw(85),nroots(2)
      common/flips/ib123(13+4*mxprms)
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
      common/geom/axyz(34)
      common/pqgeom/abcd(24)
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
      common/type/jtype(2)
      common/pgeom/gep(8*mxprms*mxprms)
      common/qgeom/acx(6)
      common/maxc/cmax(mxprim),cmaxx(4*mxprms+2),ismlp(mxprms*mxprms+2)
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
      common/astore/qq,theta,n
      common/miscg/mab,mcd,ngangb
      common/lt/labcdt(4)
c
      real*8 conp
c
c     conp = prefactors from pairs of primitives
c
      common /const/ conp(mxprms*mxprms)
c
      common/shllfo/shllf(12*mxprms+8)
c
      real*8 en, accd
      integer ig, ierec, igrec, iterat, jpoint
      common /cntl2/ en(200),ig(200),ierec(200),igrec(200),
     +               accd(6,200),iterat,jpoint
c
      common/limy/mpppp(23+mxroot)
      common/reso/iresco(400)
      common/prpspc/centx(3),mopic(3)
c
      real*8 stos1, stos2, stos3, stos4, stos5, stos6, stos7
      real*8 rleesf
      real*8 c11al, c11be, c31al, c31be, c33al, c33be
      integer nbfs, minf, maxf, nangm
c
      common /datbas/ stos1(54),stos2(54),stos3(54),stos4(54),
     + stos5(54),stos6(54),stos7(54),
     + rleesf(36,2),
     + nbfs(24),minf(24),maxf(24),nangm(24),
     + c11al(12),c11be(12),c31al(12),c31be(12),c33al(12),c33be(12)
c
c
      real*8 alphas,chargat,spinat
      integer nswapa, nswapb, nswap, next,isecat
      integer natconf,nnatconf,iatcon,iatstates,natdiff
      parameter (nnatconf=20)
      logical odiff,oexc,oground,oatdo,oalway,uhfatom,nonrel,forcat
      character*8 zatconf,atmode,zatdiff
      character*60 string_ch
      character*2 zatstate
c
      common /datgue/ alphas(maxorb),nswapa,nswapb,nswap(80),
     +                next(maxorb),oground,oatdo,isecat,oalway,
     +                odiff,oexc
      common /datgue2/ chargat(nnatconf),spinat(2,4,nnatconf),
     +                 natconf,uhfatom,iatcon,nonrel(nnatconf),
     +                 forcat(nnatconf),iatstates(nnatconf),
     +                 natdiff
      common /datgue3/ string_ch(2,nnatconf),zatconf(nnatconf),atmode,
     +                 zatstate(nnatconf),zatdiff(nnatconf)
c
      real*8 bias, ajm, ajm1, ajm2
      integer nvar, mope, ncs, nos, npopul
c
      common/datplt/bias(8),nvar(14),mope(2),ajm(110),ajm1(110),
     + ajm2(110),ncs(90),nos(90),npopul(90)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      common /miscop/ cmin(16,maxvar),
     2 spppp(31),fspa(3,maxat),fspaa(maxvar),ia(maxat+1)
      common /copt / rcvopt
      common /sver  / tver  (18)
c
      integer mapsie
      common /linkan/ mapsie(maxorb)
c
      common /dshlno/ dshlno(22)
c
      real*8 de
      common/grad2/de(3,maxat)
c
      common /dshlnf/ dshlnf(24*mxprms + mxprms*mxprms + 15),idshln(4)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      common /setd  / setd  (ncmax,10),setdd(6),isetd(34)
      common /ffq/ fq0,fq1,fq2,fq3,fq4,fq5
c...   for UHF natorbs
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
c
      integer ix, iy, iz
      common /inxblk/ ix(35),iy(35),iz(35)
c
      common/rtdata/whi(45),wlow(45),rhi(45),rlow(45),rfac(45),
     1amps(9),ipoint(9),madd(20)
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
      common/restrz/zgam(50)
      common/restrr/realg(30)
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      common/restrl/logg(50)
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
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      real*8 fieldx, fieldy, fieldz
      logical ofield
      common /gms_field/  fieldx, fieldy, fieldz, ofield
c
c ...
c     simons & jorgensen optimisation.
c ...
      common /zjorgs/ystat(3)
      common /jorgs/maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     *outjor,outdeb,maxsta
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
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
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
c The oblock array determines what is written out to the punchfile:
c
c Below sourced from server.m
c oblock(1)   - coords - single point or last point
c keyword: coor
      integer PUN_COORD
      parameter (PUN_COORD=1)
c oblock(2)   - coordinates at each optimisation step
c keyword: opti
      integer PUN_OPT_COORD
      parameter (PUN_OPT_COORD=2)
c oblock(3)   - starting coordinates
c keyword: init
      integer PUN_START_COORD
      parameter (PUN_START_COORD=3)
c oblock(4)   - atom connectivity
c keyword: conn
      integer PUN_ATOM_CONN
      parameter (PUN_ATOM_CONN=4)
c oblock(5)   - molecular orbital occupations
c keyword: occu
      integer PUN_MO_OCC
      parameter (PUN_MO_OCC=5)
c oblock(6)   - eigenvectors
c keyword: vect
      integer PUN_EIGV
      parameter (PUN_EIGV=6)
c oblock(7)   - run type
c keyword: type
      integer PUN_RUNT
      parameter (PUN_RUNT=7)
c oblock(8)   - gvb pair data
c keyword: gvb
      integer PUN_GVBPAIR
      parameter (PUN_GVBPAIR=8)
c oblock(9)   - force constant matrix
c keyword: forc
      integer PUN_FCMAT
      parameter (PUN_FCMAT=9)
c oblock(10)  - vibrational frequencies
c keyword: vibr
      integer PUN_VIBFREQ
      parameter (PUN_VIBFREQ=10)
c oblock(11)  - basis set
c keyword: basi
      integer PUN_BASIS
      parameter (PUN_BASIS=11)
c oblock(12)  - orbital energies
c keyword: eige
      integer PUN_ORB_ENER
      parameter (PUN_ORB_ENER=12)
c oblock(13)  - total scf energies
c keyword: scfe
      integer PUN_SCF_ENER
      parameter (PUN_SCF_ENER=13)
c oblock(14)  - job title
c keyword: titl
      integer PUN_JOB_TITLE
      parameter (PUN_JOB_TITLE=14)
c oblock(15)  - overlap matrix
c keyword: over
      integer PUN_OVLP_MAT
      parameter (PUN_OVLP_MAT=15)
c oblock(16)  - grid data from dumpfile section(s)...
c keyword: grid
      integer PUN_GRID_DATA
      parameter (PUN_GRID_DATA=16)
c oblock(17)  - mulliken analysis
c keyword: mull
      integer PUN_MULLIKEN
      parameter (PUN_MULLIKEN=17)
c oblock(18)  - lowdin analysis
c keyword: lowd
      integer PUN_LOWDIN
      parameter (PUN_LOWDIN=18)
c oblock(19)  - spin densities
c keyword: spin
      integer PUN_SPIN_DENS
      parameter (PUN_SPIN_DENS=19)
c oblock(20)  - scf type
c keyword: leve
      integer PUN_SCF_TYPE
      parameter (PUN_SCF_TYPE=20)
c oblock(21)  - normal coordinates for vibrational modes
c keyword: norm
      integer PUN_VIB_MODES
      parameter (PUN_VIB_MODES=21)
c oblock(22)  - two electron integrals
c keyword: twoe
      integer PUN_2E_INT
      parameter (PUN_2E_INT=22)
c oblock(23)  - gradients
c keyword: grad
      integer PUN_GRADIENTS
      parameter (PUN_GRADIENTS=23)
c oblock(24)  - transformation matrix
c keyword: tran
      integer PUN_TRANS_MAT
      parameter (PUN_TRANS_MAT=24)
c oblock(25)  - potential derived charges
c keyword: pdc
      integer PUN_POT_DERV_CHG
      parameter (PUN_POT_DERV_CHG=25)
c oblock(26)  - grid symmetry array
c keyword: gsym
      integer PUN_GRID_SYMM
      parameter (PUN_GRID_SYMM=26)
c oblock(27)  - cartesian hessian matrix
c keyword: secd
      integer PUN_HESSIAN
      parameter (PUN_HESSIAN=27)
c oblock(28)  - distributed multipole analysis
c keyword: dma
      integer PUN_MPOLE_ANAL
      parameter (PUN_MPOLE_ANAL=28)
c oblock(29)  - density matrix
c keyword: dens
      integer PUN_DENS_MAT
      parameter (PUN_DENS_MAT=29)
c oblock(30)  - total energy
c keyword: ener
      integer PUN_TOT_ENER
      parameter (PUN_TOT_ENER=30)
c oblock(31)  - internal coordinate hessian matrix
c keyword: hess
      integer PUN_ZMAT_HESS
      parameter (PUN_ZMAT_HESS=31)
c oblock(32)  - dipole moment
c keyword: dipo
      integer PUN_DIPOLE
      parameter (PUN_DIPOLE=32)
c oblock(33)  - infrared intensities 
c keyword: infr
      integer PUN_INFRARED
      parameter (PUN_INFRARED=33)
c oblock(34)  - raman intensities
c keyword: rama
      integer PUN_RAMAN
      parameter (PUN_RAMAN=34)
c CURRENTLY ONLY USED FOR XML
c oblock(35)  - metadata (user, compilation date etc)
c keyword: meta
      integer PUN_METADATA
      parameter (PUN_METADATA=35)
c CURRENTLY ONLY USED FOR XML
c oblock(36)  - input parameters
c keyword: inpa
      integer PUN_INPUT_PARAM
      parameter (PUN_INPUT_PARAM=36)


      integer LENPUN
      parameter(LENPUN=60)

      logical oblock, opunold
      integer nfblck, nblsec, iblsec, nblgri, iblgri
      integer itwoe, ntwoe, iblfmt
      common/blocks/nfblck,nblsec,iblsec(10),nblgri,iblgri(10),
     +              itwoe,ntwoe,iblfmt(LENPUN),oblock(LENPUN),opunold
c

c
c symmetry assignment
c
      common/atmol3o/osim1,osim2
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
c
c
      real*8 func0
      integer ncoord, npts, nserch, iupdat, icode, iupcod
      common /seerch/ func0,ncoord,npts,nserch,iupdat,icode,iupcod
c
      common/l/lockk
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     * no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/scfopt/maxit,mconv,nconv,npunch,
     * accdi1,accdi2,odiis,icoupl(3),dmpcut,acccc(6),iaccv(2),
     * rshift(4),iextn,junks,dampoo(7)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     * ctran(mxorb3),otran,otri
c
      integer ilifd, ntrad, itrad, loc, ncomp
      real*8 ctrad
      common/tranao/ilifd(maxorb),ntrad(maxorb),itrad(mxorb3),
     * ctrad(mxorb3),loc(maxorb),ncomp(maxorb)
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
      real*8 tr, trx, try, trz, prmoms, aprev
      integer indmx, jaxis, igroup, igrp80
      logical oprmom
      common /molsym/ tr(3,3),trx,try,trz,indmx,jaxis,igroup,igrp80,
     +                prmoms(3),aprev(maxat3,3),oprmom
c
      common/dnfnw/oden(3),iiden(3),vdnfnw(8,200),rdnfnw(1275),
     *  tdnfnw(50,3),ndnfnw(50),eltotxx,nomoxx
c
c     ------ scrf common blocks
c
      real*8 aradix,aradiy,aradiz, dielec
      logical oscrf
      common/scrf/aradix,aradiy,aradiz,dielec,oscrf
      common/dplrfc/tele(3),tnuc(3),tmol(3),dtot
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
c
c ----- Pisa solvation common blocks
c
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
      common/surfgn/numbq(maxat),irejec(maxat),cavsze(maxat)
     * ,nocent,itmax,damper
c
c
      real*8 radvdw
      common/vdwrad/radvdw(105)
c
c approximate covalent radius table in au
      real*8 cov
      common/coval/cov(100)
      common/indxsv/indx
      common/orb/nint(9)
c
c ***** hondo common blocks ****
c
      common/blk1/realx(510),irealx(4)
      common/sortpk/spk(3412)
      common/tapes/iko(750)
cjmht does not appear to be needed. It's now in the common directory
c     as an include file if it should be.
cjmht      common/funcon/bohr(9)
      common/cigrad/eds(130)
      common/specal/ipqrs(14)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/jspec/nspj(maxorb+1)
      common/maxlen/maxq,maxqq,omem(2*maxlfn+2)
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
      logical oqskip,omixed
      common/qercom/frqq(31),nfrqq,npoleb,iq6(36),
     + iqsec(50),oqskip(50),iqang(50),npb,ijunq(4),omixed,ifileb,iblkb
      character *8 qnames
      common/qercmx/qnames(50)
      common/out/aj(512)
      common/scfblk/al(34)
c
      real*8 clat, zlat
      integer nlat, maxlat
      parameter (maxlat=216)
      common /lattic/ clat(3,maxlat),zlat(maxlat),nlat
c
      character *8 groupy
      common/indsyx/groupy
      common/symmos/imos(43,maxorb),immos(162)
      common/defunk/cobas(maxorb,3),idatbf(2,maxat),
     + ijdd1(8,numspl),ijdd2(11)
c
      real*8 adiis
      common/worksp/adiis(286)
c
      common/zlabs/izlabs(408)
      common/small/dipd(3,maxat*3),pold(3,3,maxat*3) ,vcd(maxat*3,6)
c     common/mcgrad/av(8)
      common/incrs/bc(60)
      common/crnams/ iangc(40),iangs(24)
      character *8 pnamc,ac6,pnams
      character *4 bfnam,atnam
      common/crnamx/ac6(6),pnamc(40),pnams(24),bfnam(20),atnam(37)
c flag to avoid recursive calls to caserr
      common/iotrp/ocaser

c *************************************
c ***** end of hondo common blocks ***
c *************************************
cdrf start
      common/c_of_m/pcm,qcm,rcm
cdrf end
c
      character*2 zelem(0:103)
      integer nelem
      parameter(nelem=103)
      common/periodic/zelem
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c
      real*8 aoccc
      integer isetc, iync, ns_canon, isec_canon
      logical opr_canon,o_canon,oset_canon
      common/canon_nat/aoccc(10),isetc(10),iync(10),
     1                 ns_canon,isec_canon,
     2                 opr_canon,o_canon,oset_canon
c
c...  common to controL fiddling of (inverse) hessians e.g. resetting
      INTEGER icall_h,ifreq_h,khes7_h
      common/fiddle_hesc/icall_h,ifreq_h,khes7_h
c
      data field/' '/
      data version/' '/
      data pcm,qcm,rcm/3*0.0d0/
      data obeen,obeen2,obeen3 /3*.false./
c
      data ac6/'c6','c7','c8','c9','c10','cn'/
      data iangs/1,1,1,  2,2,2,2,2,   3,3,3,3,3,3,3,
     1   4,4,4,4,4,4,4,4,4/
      data iangc/1,1,1,  2,2,2,2,2,2,  3,3,3,3,3,3,3,3,3,3,
     1   4,4,4,4, 4,4,4,4,4, 4,4,4,4,4,4   ,0,0,0,0,0,0/
      data pnams/'q10','q11c','q11s','q20','q21c','q21s','q22c','q22s',
     1    'q30','q31c','q31s','q32c','q32s','q33c','q33s',
     2 'q40','q41c','q41s','q42c','q42s','q43c','q43s','q44c','q44s'/
      data pnamc/'x','y','z','xx','yy','zz','xy','xz','yz','xxx',
     1         'yyy','zzz','xyy','yzz','xxz','xzz','xxy','yyz',
     2         'xyz','xxxx','yyyy','zzzz','xxyy','yyzz','xxzz',
     3         'xxxy','xxxz','yyyz','xyyy','xzzz','yzzz','xyzz',
     4         'xyyz', 'xxyz','lx','ly','lz','px','py','pz'/
c
c     runlab
c
      data zbflab/maxorb*'        '/
c
c
c     set data in /rtdata/ used by integration routines
c
      data madd/17,17,18,18,19,19,20,20,21,21,22,22,
     &   23,23,24,24,25,25,26,26/
c
c     roots and weights generated using nag routine d01bcf
c     rhi contains squares of gauss-hermite roots of order
c     2,4,6.......18
c     whi are corresponding weights
c     rlow contains squares of gauss-legendre roots of order
c     2,4,6,......,18
c     wlow are corresponding weights
c     rfac contains ratios rlow/rhi
c
      data rhi/
     &0.50000000000000d+00,0.27247448713916d+01,0.27525512860841d+00,
     &0.55253437422632d+01,0.17844927485433d+01,0.19016350919349d+00,
     &0.85886356890120d+01,0.39269635013583d+01,0.13390972881264d+01,
     &0.14530352150332d+00,0.11807189489972d+02,0.64147297336620d+01,
     &0.30859374437175d+01,0.10745620124369d+01,0.11758132021178d+00,
     &0.15129959781108d+02,0.91242480375311d+01,0.51961525300544d+01,
     &0.25525898026682d+01,0.89830283456962d+00,0.98747014068481d-01,
     &0.18528277495852d+02,0.11989993039824d+02,0.75540913261018d+01,
     &0.43897928867310d+01,0.21805918884504d+01,0.77213792004278d+00,
     &0.85115442997594d-01,0.21984272840963d+02,0.14972627088426d+02,
     &0.10093323675221d+02,0.64831454286271d+01,0.38094763614849d+01,
     &0.19051136350314d+01,0.67724908764929d+00,0.74791882596818d-01,
     &0.25485979166099d+02,0.18046505467729d+02,0.12771825354869d+02,
     &0.87697567302685d+01,0.56944233429577d+01,0.33691762702432d+01,
     &0.16923950797932d+01,0.60323635708175d+00,0.66702230958194d-01/
      data whi/
     &0.88622692545276d+00,0.81312835447246d-01,0.80491409000551d+00,
     &0.45300099055090d-02,0.15706732032286d+00,0.72462959522439d+00,
     &0.19960407221138d-03,0.17077983007414d-01,0.20780232581489d+00,
     &0.66114701255824d+00,0.76404328552331d-05,0.13436457467813d-02,
     &0.33874394455481d-01,0.24013861108232d+00,0.61086263373533d+00,
     &0.26585516843565d-06,0.85736870435883d-04,0.39053905846291d-02,
     &0.51607985615884d-01,0.26049231026416d+00,0.57013523626248d+00,
     &0.86285911681258d-08,0.47164843550191d-05,0.35509261355194d-03,
     &0.78500547264582d-02,0.68505534223466d-01,0.27310560906425d+00,
     &0.53640590971209d+00,0.26548074740116d-09,0.23209808448653d-06,
     &0.27118600925380d-04,0.93228400862421d-03,0.12880311535510d-01,
     &0.83810041398987d-01,0.28064745852853d+00,0.50792947901661d+00,
     &0.78281997721164d-11,0.10467205795793d-07,0.18106544810936d-05,
     &0.91811268679300d-04,0.18885226302685d-02,0.18640042387545d-01,
     &0.97301747641316d-01,0.28480728566998d+00,0.48349569472546d+00/
      data rlow/
     &0.33333333333333d+00,0.74155574714581d+00,0.11558710999705d+00,
     &0.86949939491826d+00,0.43719785275109d+00,0.56939115967007d-01,
     &0.92215660849206d+00,0.63467747623464d+00,0.27618431387246d+00,
     &0.33648268067507d-01,0.94849392628837d+00,0.74833462838728d+00,
     &0.46159736149627d+00,0.18783156765245d+00,0.22163568807218d-01,
     &0.96346127870282d+00,0.81742801326687d+00,0.59275012773154d+00,
     &0.34494237942742d+00,0.13530001165525d+00,0.15683406607401d-01,
     &0.97275575129749d+00,0.86199133320339d+00,0.68426201565315d+00,
     &0.47237153700448d+00,0.26548115726894d+00,0.10183270400277d+00,
     &0.11675871940146d-01,0.97891421016235d+00,0.89222197421380d+00,
     &0.74931737854740d+00,0.57063582016217d+00,0.38177105339712d+00,
     &0.20977936861551d+00,0.79300559811486d-01,0.90273770256471d-02,
     &0.98320148322563d+00,0.91359942257427d+00,0.79673916319752d+00,
     &0.64594166107702d+00,0.47843096553757d+00,0.31334338332122d+00,
     &0.16953901896600d+00,0.63446670693112d-01,0.71868028362264d-02/
      data wlow/
     &0.10000000000000d+01,0.34785484513745d+00,0.65214515486255d+00,
     &0.17132449237917d+00,0.36076157304814d+00,0.46791393457269d+00,
     &0.10122853629038d+00,0.22238103445337d+00,0.31370664587789d+00,
     &0.36268378337836d+00,0.66671344308689d-01,0.14945134915058d+00,
     &0.21908636251598d+00,0.26926671931000d+00,0.29552422471475d+00,
     &0.47175336386513d-01,0.10693932599532d+00,0.16007832854335d+00,
     &0.20316742672307d+00,0.23349253653835d+00,0.24914704581340d+00,
     &0.35119460331752d-01,0.80158087159761d-01,0.12151857068790d+00,
     &0.15720316715819d+00,0.18553839747794d+00,0.20519846372130d+00,
     &0.21526385346316d+00,0.27152459411755d-01,0.62253523938648d-01,
     &0.95158511682493d-01,0.12462897125553d+00,0.14959598881658d+00,
     &0.16915651939500d+00,0.18260341504492d+00,0.18945061045507d+00,
     &0.21616013526485d-01,0.49714548894970d-01,0.76425730254889d-01,
     &0.10094204410629d+00,0.12255520671148d+00,0.14064291467065d+00,
     &0.15468467512627d+00,0.16427648374583d+00,0.16914238296314d+00/
      data rfac/
     &0.66666666666667d+00,0.27215603006790d+00,0.41992718021791d+00,
     &0.15736566546394d+00,0.24499839134006d+00,0.29942188282334d+00,
     &0.10736939391571d+00,0.16162041638918d+00,0.20624663818033d+00,
     &0.23157228207121d+00,0.80331896688366d-01,0.11665879303696d+00,
     &0.14958091987121d+00,0.17479825778177d+00,0.18849566212812d+00,
     &0.63679037660486d-01,0.89588534847421d-01,0.11407481291265d+00,
     &0.13513427777031d+00,0.15061737139021d+00,0.15882410982598d+00,
     &0.52501143266839d-01,0.71892563268415d-01,0.90581644583619d-01,
     &0.10760679357614d+00,0.12174729195090d+00,0.13188408619684d+00,
     &0.13717689210026d+00,0.44527932183337d-01,0.59590208781963d-01,
     &0.74238913033865d-01,0.88018358749513d-01,0.10021614971993d+00,
     &0.11011383507948d+00,0.11709216189089d+00,0.12069995716395d+00,
     &0.38578132580971d-01,0.50624727552267d-01,0.62382560132156d-01,
     &0.73655596266152d-01,0.84017456504924d-01,0.93002965172432d-01,
     &0.10017697462623d+00,0.10517713322195d+00,0.10774456465678d+00/
c
c
      data amps/24.0d0,30.0d0,35.0d0,39.0d0,43.0d0,47.0d0,51.0d0,55.0d0,
     c60.0d0/
      data ipoint/0,1,3,6,10,15,21,28,36/
      data ix/0,
     +        1,0,0,
     +        2,0,0,1,1,0,
     +        3,0,0,2,2,1,0,1,0,1,
     +        4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data iy/0,
     +        0,1,0,
     +        0,2,0,1,0,1,
     +        0,3,0,1,0,2,2,0,1,1,
     +        0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data iz/0,
     +        0,0,1,
     +        0,0,2,0,1,1,
     +        0,0,3,0,1,0,1,2,2,1,
     +        0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
c ...
      data  jsym/1,2,1,3,4,1,4,3,2,1,5,6,7,8,1,6,5,8,7,2,1,7,8
     *          ,5,6,3,4,1,8,7,6,5,4,3,2,1/
      data nix/0,
     +         1,0,0,
     +         2,1,1,0,0,0,
     +         3,2,2,1,1,0,0,1,0,0,
     +         4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data niy/0,
     +         0,1,0,
     +         0,1,0,2,1,0,
     +         0,1,0,1,2,3,2,0,1,0,
     +         0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data niz/0,
     +         0,0,1,
     +         0,0,1,0,1,2,
     +         0,0,1,1,0,0,1,2,2,3,
     +         0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
      data isx /1,0,0,1,1,0,1/
      data isy /0,1,0,1,0,1,1/
      data isz /0,0,1,0,1,1,1/
      data irrep/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/
      data idh /-1,1,-1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,
     11,1,1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,
     2-1,1,1,-1/
      data jdh /-1,1,-1,1,-1,-1,-1,-1,1/
c
c      list of all possible runtype's
c
      data zrunt/
     *'scf','trans','ci','saddle','optimize','optxyz',
     *'gradient','force','adapt','analyse','gf','tda',
     *'jorgense','branch','integral','nmr','********',
     *'irc','rpagrad','rpaoptim',
     *'index4', 'polariza', 'fockder','hessian','dipder',
     *'polder','magnetiz', 'intensit','test', 'lagrange',
     *'hyperpol','lazzeret','raman'  ,'infrared',6*'--------',
     *'mopac','response','montec',7*'********'/
      data yrunt/
     *'scf ','tran','ci  ','sadd','opti','optx','grad',
     *'forc','adap','anal','gf'  ,'tda' ,'jorg','bran',
     *'inte','nmr' ,'****','irc ','rpag','rpao',
     *'inde','pola','fock','hess','dipd',
     *'pold','magn','inte','test','lagr',
     *'hype','lazz','rama','infr',6*'----',
     *'mopa','resp','mont',7*'****'/
c
c-weights
c ams data deleted, all atomic masses should now be handled by
c amass_get function
c

c
cjvl
cjvl  add 7,8,9 points hermite quadrature
cjvl   (hfff is dimensioned 7+8+9=24)
cjvl  allows up to g-functions in dma analysis
cjvl  M.Abramowitz, I.A.Stegun, Handbook of Mathematical Fucntions
cjvl  9th printing, (Dover publications, New York) (1972) page 924
c
      data h /0.0d0,
     +-.707106781186548d0,.707106781186548d0,
     +-1.22474487139159d0,0.0d0,
     + 1.22474487139159d0,-1.65068012388578d0,-.524647623275290d0,
     +.524647623275290d0,1.65068012388578d0,
     +-2.02018287045609d0,-.958572464613819d0,0.0d0,
     + .958572464613819d0,2.02018287045609d0,-2.350604973674d0,
     + -1.335849074014d0,-0.436077411928d0,0.436077411928d0,
     +  1.335849074014d0, 2.350604973674d0,
     + -2.651961356835d0, -1.673551628767d0,
     + -0.81628788285896d0,   0.0d0,  0.81628788285896d0,
     +  1.673551628767d0 , 2.651961356835d0,
     + -2.930637420257d0, -1.981656756696d0,
     + -1.157193712447d0, -0.3811869902073d0,
     +  0.3811869902073d0, 1.157193712447d0,
     +  1.981656756696d0,  2.930637420257d0,
     + -3.190993201781528d0, -2.266580584531843d0,
     + -1.468553289216668d0, -0.723551018752838d0, 0.0d0,
     +  0.723551018752838d0,  1.468553289216668d0,
     +  2.266580584531843d0,  3.190993201781528d0  /
c
      data w /1.77245385090552d0,
     +.8862269254528d0,.8862269254528d0,
     +.2954089751509d0,1.181635900604d0,
     + .2954089751509d0, 8.131283544725d-02,8.049140900055d-01,
     + 8.049140900055d-01,8.131283544725d-02,
     + 1.995324205905d-02,3.936193231522d-01,
     + 9.453087204829d-01,3.936193231522d-01,1.995324205905d-02,
     + 4.530009905509d-03,
     + 1.570673203229d-01, 7.246295952244d-01,7.246295952244d-01,
     + 1.570673203229d-01, 4.530009905509d-03,
     + 0.9717812450995d-03, 0.5451558281913d-01,
     + 0.4256072526101d0,   0.8102646175568d0, 0.4256072526101d0,
     + 0.5451558281913d-01, 0.9717812450995d-03, 0.1996040722114d-03, 
     + 0.1707798300741d-01, 0.2078023258149d0,   0.6611470125582d0,
     + 0.6611470125582d0,   0.2078023258149d0,   0.1707798300741d-01, 
     + 0.1996040722114d-03, 0.3960697726326d-04, 0.4943624275537d-02,
     + 0.8847452739438d-01, 0.4326515590026d0, 0.7202352156061d0,
     + 0.4326515590026d0,   0.8847452739438d-01,
     + 0.4943624275537d-02, 0.3960697726326d-04  /
c
      data ydd1 /
     *'safe','minb','bypa','nosy','zmat','ipri',
     *'nopr','angs','ang ','adap','****','main',
     *'mfil','seco','sfil','ffil','tfil','dmfi',
     *'twop','cifi','afil','maxb','loop','rloo',
     *'dump','inte','titl','prin','accu','supe',
     *'atom','geom','basi','gtos','stos','gaus',
     *'symm','char','norm','mult','rest','pseu',
     *'ecp' ,'form','fiel','comb','mapp','free',
     *'hfil','punc','denf','scrf','pssc','orie',
     *'symt','elec','dege','grou','vers','xfie',
     *'moro','reac','harm','zora','bqbq','blur',
     *'news','chm ','bqfo','cres','vb  ','depe',
     *'serv','ghos','disp','keys','xmli','xmlo',
     *'copq','cpwa','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****'
     * /
c     do add keywords to ydd2, three other places
c     have to be changed: file m4/common/direc, 
c     this file, common/direc/zrun ...ydd2(120),
c     this file, data ... limit2 ...
      data ydd2 /
     *'maxc','thre','npun','damp','lock',
     *'swap','conv','leve','valu','xtol',
     *'step','lsea','tolm','tols','tans','minm',
     *'diis','orbs','noci','stat','simu','sort',
     *'augm','onep','hess','cive','cano','newt',
     *'line','lagr','jkop','ciac','supa','trac',
     *'casp','supe','pass','cipa','stop','open',
     *'popl','upda','nucl','prop','cent','sele',
     *'maxj','reca','powe','bfgs','bfgx','rfo ',
     *'cuto','optp','optd','mine','maxe','pref',
     *'inne','less',
     *'ente','runt','conf','loca','boys','lmo',
     *'grap','core','irc' ,'acti','dire','i.p.',
     *'ip'  ,'mull','mult','mcsc','dma', 'mrdc',
     *'full','mode','rpa' ,'sacc','mopa','vect',
     *'scft','nato','potf','omit','ccsd','qcis',
     *'dpa', 'save','hmat','nbo', 'cdft','dft' ,
     *'smea','aver','dies','vdwa','dlfi','mass',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****',
     *'****','****','****','****','****','****' /
      data zdd3 /
     1'output','cards','vcd','gtens','field','skip','imptprin','open',
     2'ncoorb','keepfock','jacobi','oscf','grhf','contract','weights',
     7'uhf','canonica','open-sin','sections','density',
     8'mos','order','epsnes','ignore','salvage','shells','orbcore',
     b 'gauge','copyfile','reorder','restrict','irest2',
     c'irest3','irest4','irest5','mp2','seczero','closed',
     d'spdoc', 'spuoc','helfeygr','mprest','mp3','pwftol'
     e,'origin','fcm','pump2','ci','icflag'
     f,'frozen','minfc','adaptint','cpf',
     g'optorb','mpflag','nonuclea','frequenc',
     h'onelec','cutoff','chfconv','perturba','mixed',
     i'noprint','irestp','irest1','debug','mpasses','irest6','ordermos',
     j'********'/
cjvl  disabled zdd4 'check', 'intermol', 'ghost1', 'ghost2', 'supermol',
cjvl  disabled zdd4 'dispersi', 'multipol',
      data zdd4 / 'hfscf','gradient','gradone','optimize','force',
     1    '*****', '*****', '*****', '*****', '*****',
     2    '*****', '*****', 'index4', 'polariza',
     3    '       ','      ','magnetiz', 'saddle','integral',
     4    'fockder','hessian','dipder','polder','adapt',
     5    'scf','potentia','cichf','property','intensit',
     6    'optimise','energy','moguess','test','drt',
     7    'polarisa','magnetis','cisort',
     8    'lagrange','lazzeret','hyperpol',
     9     'allinfo','raman','infrared',7*'        '/
c
      data cpwa/0.0d0/
      data nbfs  /1,1,3,4,1,3,4, 6,1,3,4, 6,1,3, 6,1,3, 6,10,15,
     +   1,4, 6,4/
      data minf  /1,1,2,1,1,2,1, 5,1,2,1, 5,1,2, 5,1,2, 5,11,21,
     +   1,1, 5,1/
      data maxf  /1,1,4,4,1,4,4,10,1,4,4,10,1,4,10,1,4,10,20,35,
     +   1,4,10,4/
      data nangm /1,1,2,2,1,2,2, 3,1,2,2, 3,1,2, 3,1,2, 3, 4, 5,
     +   1,2, 3,2/
c     1s
      data stos1 /
     +1.24d0,1.69d0,
     +2.69d0,3.68d0,4.68d0,5.67d0,6.67d0,7.66d0,8.65d0,9.64d0,
     +10.61d0,11.59d0,12.56d0,13.53d0,14.50d0,15.47d0,16.43d0,17.40d0,
     +18.61d0,19.58d0,
     +20.56d0,21.54d0,22.53d0,23.52d0,24.50d0,25.49d0,26.47d0,27.46d0,
     +28.44d0,29.43d0,
     +30.42d0,31.40d0,32.39d0,33.37d0,34.36d0,35.34d0,
     +36.32d0,37.31d0,38.29d0,39.27d0,40.26d0,41.24d0,
     +42.22d0,43.21d0,44.19d0,45.17d0,46.15d0,47.14d0,
     +48.12d0,49.10d0,50.08d0,51.07d0,52.05d0,53.03d0/
c     2sp
      data stos2 /
     +0.00d0,0.00d0,
     +0.80d0,1.15d0,1.50d0,1.72d0,1.95d0,2.25d0,2.55d0,2.88d0,
     +3.48d0,3.90d0,4.36d0,4.83d0,5.31d0,5.79d0,6.26d0,6.74d0,
     +7.26d0,7.74d0,
     +8.22d0,8.70d0,9.18d0,9.66d0,10.13d0,10.61d0,11.09d0,11.56d0,
     +12.04d0,12.52d0,
     +12.99d0,13.47d0,13.94d0,14.40d0,14.87d0,15.34d0,
     +15.81d0,16.28d0,16.72d0,17.19d0,17.66d0,18.12d0,
     +18.59d0,19.05d0,19.51d0,19.97d0,20.43d0,20.88d0,
     +21.33d0,21.79d0,22.25d0,22.71d0,23.17d0,23.63d0/
c     3sp
      data stos3 /
     +10*0.0d0,
     +1.75d0,1.70d0,1.70d0,1.75d0,1.90d0,2.05d0,2.10d0,2.33d0,
     +2.75d0,3.01d0,
     +3.21d0,3.44d0,3.67d0,3.89d0,4.11d0,4.33d0,4.55d0,4.76d0,
     +4.98d0, 5.19d0,
     +5.26d0,5.58d0,5.90d0,6.22d0,6.54d0,6.86d0,
     +7.18d0, 7.49d0, 7.97d0, 8.21d0, 8.51d0, 8.82d0,
     +9.14d0, 9.45d0, 9.77d0,10.09d0,10.41d0,10.74d0,
     +11.08d0,11.39d0,11.71d0,12.03d0,12.35d0,12.66d0/
c     3d
      data stos4 /
     +20*0.0d0,
     +1.10d0,1.90d0,2.55d0,3.05d0,3.45d0,3.75d0,4.10d0,4.35d0,
     +4.60d0, 4.90d0,
     +5.26d0,5.58d0,5.90d0,6.22d0,6.54d0,6.86d0,
     +7.18d0, 7.49d0, 7.97d0, 8.21d0, 8.51d0, 8.82d0,
     +9.14d0, 9.45d0, 9.77d0,10.09d0,10.41d0,10.74d0,
     +11.08d0,11.39d0,11.71d0,12.03d0,12.35d0,12.66d0/
c      4sp
      data stos5 /
     +18*0.0d0,
     +1.43d0,1.36d0,
     +1.60d0,1.70d0,1.70d0,1.75d0,1.65d0,1.55d0,1.55d0,1.60d0,
     +1.60d0,1.90d0,
     +1.80d0,2.00d0,2.12d0,2.22d0,2.38d0,2.54d0,
     +3.02d0, 3.16d0, 3.29d0, 3.48d0, 3.67d0, 3.87d0,
     +4.05d0, 4.24d0, 4.41d0, 4.59d0, 4.76d0, 4.93d0,
     +4.65d0, 4.89d0, 5.12d0, 5.36d0, 5.59d0, 5.82d0/
c     4d
      data stos6 /
     +36*0.0d0,
     +0.00d0, 0.00d0, 1.40d0, 1.95d0, 2.40d0, 2.70d0,
     +3.00d0, 3.20d0, 3.45d0, 3.60d0, 3.75d0, 3.95d0,
     +4.65d0, 4.89d0, 5.12d0, 5.36d0, 5.59d0, 5.82d0/
c     5sp
      data stos7 /
     +36*0.0d0,
     +1.90d0, 1.80d0, 1.80d0, 1.90d0, 1.90d0, 1.95d0,
     +1.85d0, 1.75d0, 1.75d0, 1.80d0, 1.80d0, 2.10d0,
     +2.05d0, 2.15d0, 2.20d0, 2.28d0, 2.42d0, 2.57d0/
c
      data rleesf/1.20d+00,1.00d+00,1.03d+00,1.03d+00,1.03d+00,1.00d+00,
     1            0.99d+00,0.99d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     2            1.00d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     3         12*1.00d+00,6*0.00d+00,
     4            1.15d+00,1.00d+00,1.12d+00,1.12d+00,1.12d+00,1.04d+00,
     5            0.98d+00,0.98d+00,1.00d+00,1.00d+00,1.00d+00,1.00d+00,
     6            1.00d+00,1.00d+00,1.02d+00,1.01d+00,1.01d+00,1.00d+00,
     7         12*1.00d+00,6*0.00d+00/
      data c11al/2*0.375d0,2*-0.125d0,2*0.125d0,2*-0.125d0,2*0.375d0,
     c2*0.125d0/
      data c11be/-0.25d0,2*0.25d0,4*-0.25d0,0.25d0,-0.75d0,0.75d0,0.25d0
     c,0.75d0/
      data c31al/2*0.125d0,2*0.625d0,2*0.375d0,2*0.625d0,2*0.125d0,
     c2*0.375d0/
      data c31be/-0.25d0,0.75d0,5*-0.25d0,-1.25d0,3*-0.25d0,0.75d0/
      data c33al/2*1.375d0,2*0.875d0,2*1.125d0,2*0.875d0,2*1.375d0,
     c2*1.125d0/
      data c33be/-0.75d0,2*-0.25d0,2*-0.75d0,-0.25d0,-0.75d0,-0.25d0,
     *-1.25d0,0.25d0,-0.25d0,-0.75d0/
      data ncom/1,6,3,6,3,6,6,7,10,10,3,10,10,9,1,3,3,1,1,6,3,3/
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,
     +         0,2,2,0,2,1,1,5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1,
     +         6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/
      data kin/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,
     +         1,2,0,2,1,2,1,0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2,
     +         0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2/
      data nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,3,
     +         3,0,2,2,1,1,2,0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2,
     +         0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2/
      data ngmx,nfmx,nomx,ncmx,ntmx,ntcol,icon,npmx/
     *1000,1024,1024,100,20,20,24*0,100/
      data cconv/
     *9.07618d0,17.7497d0,.171524d0,.324123d0,2.54154d0,
     *1.344911d0,.79198d0,.2800116d0,.711688d0,.14818d0,.14818d0/
c
c-weights
c atwt data deleted, all atomic masses should now be handled by
c amass_get function
c
      data itps,maxnt/5,22/
      data cspin/
     *2.792780d0,-2.1275d0,1.085467d0,-0.392533d0,0.896167d0,
     *0.7024d0,0.2018d0,-0.37874d0,2.6288d0,-0.2206d0,
     *.739167d0,-.17106d0,.72828d0,-.5553d0,1.1317d0,
     *.2144333d0,.273943d0,-.185714d0,.1304667d0,-.188143d0,
     *.679486d0,-.15766d0,.7355714d0,-.1581333d0,.6888d0,
     *.0902d0,.6641429d0,-.2495667d0,.742d0,.1751d0/
c
      data bias,nvar/
     *5.0d0,.0001d0,.001d0,.1d0,.000001d0,.001d0,.0d0,.0d0,
     *0,10,10,50,0,0,0,0,1,1,1,0,1,1/
      data ncs/
     *0,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,
     *0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/
      data nos/
     *1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1/
      data npopul/
     *1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,1,2,3,4,5,0,0,0,1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,0 /
      data ajm/
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,-1.6666666667d0,.1333333333d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,-.6666666667d0,-0.0666666667d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.33333333333d0,-.1333333333d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.1666666667d0,
     *-0.0166666667d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-0.0666666667d0,0.0053333333d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ajm1/
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     * 0.0d0,0.0d0,-1.6666667d0, 0.1333333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.6666667d0,-0.0666667d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.3333333d0,-0.1333333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.1666667d0,-0.0166667d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,0.0d0,-0.0666667d0, 0.0053333d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *0.0d0,0.0d0,0.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-1.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ajm2/
     *.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-1.8d0,.05714286d0,.05714286d0,
     *0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,-.8d0,-.10612245d0,
     *.0367347d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *-.4666666667d0,-.07891156d0,-.0154195d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,.0d0,-.3d0,-0.05d0,-0.05d0,.0d0,.0d0,.0d0,.0d0,.0d0,
     *.0d0,.0d0,0.0d0,-.2d0,-0.05714286d0,-.05714286d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.13333333d0,-0.022222222d0,-0.022222222d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.08571429d0,-0.01449396d0,-0.002832153d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.05d0,-0.00663265d0,0.002295918d0,0.0d0,
     $0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *-0.022222222d0,0.000705467d0,0.000705467d0,0.0d0,
     * 0.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0,.0d0/
      data ird,iwr,ipu/ 5,6,7/
      data ntapes/1,2,3,4,8/
      data ooprnt,opun7/.false.,.false./
      data omax,maxbrd/.true.,0/
      data ologf,logf/.false.,0/
      data ti,tx/ 0.0d0,0.0d0/
c     limit2= dimension of ydd2
      data limit1,limit2,limit3,mxtask/
     *   80,120, 120, 50/
      data prmoms/3*-1.0d0/
      data indx/0/
      data cov/
     * 0.47d0,3.80d0,    
     * 2.76d0,1.99d0, 1.62d0,1.33d0,1.23d0,1.14d0,0.95d0,3.80d0,
     * 3.42d0,2.85d0, 2.38d0,2.09d0,1.90d0,1.90d0,1.90d0,3.80d0,
     * 4.18d0,3.42d0,
     * 3.04d0,2.66d0,2.57d0,2.66d0,2.66d0,
     * 2.66d0,2.57d0,2.57d0,2.57d0,2.57d0,
     * 2.47d0,2.38d0,2.18d0,2.18d0,2.18d0,3.80d0,
     * 4.46d0,3.80d0,               
     * 3.42d0,2.94d0,2.76d0,2.76d0,2.57d0,
     * 2.47d0,2.57d0,2.66d0,3.04d0,2.94d0,
     * 2.94d0,2.76d0,2.76d0,2.66d0,2.66d0,3.80d0,
     * 4.94d0,4.09d0,
     * 3.71d0,3.52d0,3.52d0,3.52d0,3.52d0,
     * 3.52d0,3.52d0,3.42d0,3.33d0,3.33d0,
     * 3.33d0,3.33d0,3.33d0,3.33d0,
     * 3.33d0,2.94d0,2.76d0,2.57d0,2.57d0,
     * 2.47d0,2.57d0,2.57d0,2.57d0,2.85d0,
     * 3.61d0,3.42d0,3.04d0,3.61d0,3.61d0,3.80d0,
     * 4.94d0,4.09d0,
     * 3.71d0,3.42d0,3.42d0,3.33d0,3.33d0,
     * 3.33d0,3.33d0,3.23d0,3.13d0,3.13d0,
     * 3.13d0,3.13d0 /
     
c      data cov /0.57d0, 2.50d0,
c     * 2.53d0, 2.00d0, 1.53d0, 1.46d0, 1.42d0, 1.38d0, 1.34d0,3.00d0,
c     * 2.91d0, 2.69d0, 2.46d0, 2.23d0, 2.08d0, 1.93d0, 1.87d0,3.50d0,
c     * 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,2.00d0,
c     * 2.00d0, 3.00d0, 2.00d0, 3.00d0,
c     * 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,
c     * 64*2.0d0/
c
c ----- van der waals radii data in au, atomic number ordering
c
      data radvdw/
     *     2.26767d0,1.87083d0, 2*0.0d0, 4.10071d0, 3.495995d0, 
     *     2.83459d0, 2.645617d0, 2.551131d0, 3*0.0d0, 4.78101d0, 
     *     4.23299d0, 3.590481d0, 3.495995d0, 3.40151d0, 14*0.0d0,
     *     3.81725d0, 3.77945d0,3.77945d0,3.684967d0, 15*0.0d0,
     *     4.15740d0, 4.15740d0,  4.062913d0, 26*0.0d0,
     *     2.83459d0, 25*0.0d0 /
      data ocaser/.false./
      data nseto,ncores,npair,nset,nopen,nope/6*0/
      data nfblck/-1/
      data timsec,walsec,tstart,estart/102*0.0d0/
c
c common/periodic
c
      data zelem/
     $         'bq',
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','es','fm','md','no','lw'   /
c
c  common/harmon
c
      data oharm,opharm/.false.,.false./,odepen/.false./
c
c  common/overrule
c
      data i_depen_check,i_global_diag/-1,-999/
      data n_depen_check,o_depen_print/0,.false./
c
c  common/zorac
c
      data ozora,onlys,onlyp,onlyd,onlyf,onlyg,oscalz,opzora,
     1     op2zora,ocontr,oioraz
     2    /5*.false.,         .true., 5*.false./
      data oint_zora,o1e_zora,osmall_zora/3*.false./
      data icoul_z,isexch_z/2,0/
      data niter_z/100/,opre_zora/.false./,int_zora/1/
      data onoext_z,oso/2*.false./
      data cspeed/137.0359895d0/
      data is_z/0/,igauge_z/0/,nat_z/2/
      data critat_z/1.0d-6/
      data oatint_z,oatscf_z,oscalatz/2*.false.,.false./
      data irest_z,numdu_z/2*0/
c
c  common canon_nat
c
      data o_canon,opr_canon/2*.false./
c
c  common datgue
c
      data natconf/0/
c
c  common fiddle_hesc
c
      data icall_h,ifreq_h/0,987654/
c
      end
      block data savscf
c
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...  block data to initialise variables that could be better
c...  handled by save statements in subroutines
c...  unfortunately some compilers have trouble with those
c...  (variables have been checked not to cause interfeence
c...   in any of the other routines mentioned below)
c
c...  subroutines affected :
c     randd   (anala)   : u
c     adapt   (guess)   : idadap
c     chekpp  (input)   : ifirst
c     sinfo   (integc)  : ijato,klato,ijshlo
c     fj      (integb)  : common/fjsav/dp,dm,tab,hm,ep,em
c     fm      (integb)  : common/fmsav/t,est,ect
c     fsi     (integb)  : common/fsisav/errfs,dawfs
c     extrpd/extrpm(scf): common/extsav/odampr,oextpr

c
c...  machine specific uses are left alone in:
c     mxmb (FPS-MAX) (machscf) : ienter,ndots,lenofp,lenofl,ifun
c     ioerr (Apollo) (machscf) : igain
c     timef (atpghs) (machscf) : ienter,tst
c     get (io for ?) (machscf) : irbuf,iselb
c     oipsci (ipsc)  (machscf) : icount
c     proc21 (FPS)   (util3)   : arg
c     jksup1 (FPS)   (util3)   : arg
c
c..   in model are saved but not pre-initialised :
c     model (model) : nr,nr3,len173
c
c     for fixnag (util3/scf) seperate commons nagsav/nagsat are used
c
      common /saveco/ u,idadap,ifirst,ijato,klato,ijshlo
c
      parameter (maxf=10)
      character*20 fail(maxf)
      integer nfail(maxf),ifails(maxf)
      common/nagsat/ fail
      common/nagsav/ nfail,nf,ifails
c
      common/fjsav/dp,dm,tab,hm,ep,em
      common/fmsav/t,est,ect
      common/fsisav/errfs,dawfs
c
      common/extsav/odampr,oextpr
c
      data nfail/maxf*0/,nf/0/,ifails/maxf*0/
c
      end
      block data  fock_critical_data
      implicit none
c
c
c...   storage for fock contributions
c
      real*8 fockc
      integer*4 iposc,maxfock,nwfc,nwfcc
      logical ofipc,centurion,oshmem_inuse,ofipc_debug
      logical ofipc_shmp,ofipc_nofock,ofipc_root
      parameter (maxfock=2)
      common/fock_cratic/ fockc(maxfock),iposc(maxfock),nwfc,nwfcc,
     +                    ofipc,oshmem_inuse,ofipc_debug,ofipc_shmp,
     +                    ofipc_nofock,ofipc_root
      data ofipc,oshmem_inuse,nwfc/.false.,.false.,0/
      data ofipc_debug,ofipc_shmp,ofipc_nofock,ofipc_root/4*.false./
c
      end

      subroutine ver_mains(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mains.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
