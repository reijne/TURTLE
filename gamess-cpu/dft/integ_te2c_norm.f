c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/dft/integ_te2c_norm.m,v $
c  $State: Exp $
c
      subroutine te2c_rep_norm(icd_tag,gout,rmatrix)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
***************************************************************
* common blocks used
**************************************************************
**** from main-line code
*********************************************************************
* cslosc * files * infoa *  parallel * parcntl 
* prints * restar * sizes * statis * timez
***************************************************************
* local common blocks
*
* tabinx, junkx, mapperx, flipsx, flip7x, ijlabx
* pqgeomx, qgeomx, pgeomx, geomx, shlg70x, piconx, miscgx, astorex
* ginfx, bshellx, constx, shllfox, maxcx, savecx, typex
* shlinfx
* shltx, miscx, setintx, rootx, denssx, shlnosx, indezx, rtdatx
* inxblkx
********************************************************************
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
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
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

c
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
ccINCLUDE(common/dft_auxvar)
      integer iwr
      common/dft_iofile/iwr
c
      dimension gout(*)
      dimension rmatrix(*)
c
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      call cpuwal(begin,ebegin)
c ***
c *** establish arrays normally set in main line code
c ***
c
c ... /mapperx/
c
      do 30 i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
 30   continue
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue
c
c ... /auxvarx/
c
c      ofast = .false.
c
c /cslosc/
c
      call setsto(10,0,intcut)
      call setsto(1060,0,intmag)
      nopk = 1
c
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_dft
      nopkr = nopk
      iofrst = iofsym
      nindmx = 1
      call jkint_gamess_norm(icd_tag,gout,rmatrix)
      call final_dft
      return
      end
      subroutine jkint_gamess_norm(icd_tag,gout,rmatrix)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
c
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common /shlg70x/ ishell,jshell,kshell,lshell,
     +                 inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      integer iwr
      common/dft_iofile/iwr
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c                      1296   if d or m shells
c                     10000   if f shells
c                     50625   if g shells
c
c     ----- this version can handle g shells   -----
c
      common/junkx/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      dimension ib(4,4)
      dimension gout(*)
      dimension rmatrix(*)
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
c
c     ----- two-electron integrals -----
c
      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
        igt(loop) = ib(1,ii)
        jgt(loop) = ib(2,ii)
        kgt(loop) = ib(3,ii)
        lgt(loop) = ib(4,ii)
 30   continue
c
c     ----- set some parameters -----
c
c
      imc=0
c
c     l2 = iky(numorb(1))+numorb(1)
c
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     are the rotated axis integrals to be used only?
c
      ogauss = .true.
c
cc      if (intg76.ne.0) call filmax_2c
c
c     ----- ishell -----
c
      jnn=0
      do 150 ii = 1 , nshell(icd_tag)
c
c     ----- print intermediate restart data -----
c
         ishell = ii
         dt0 = time - tim0
         dt1 = time - tim1
         tim1 = time
         ikyii = iky(ii)
c
c     ----- jshell -----
c
         j0 = jst0
         jj=1
         jst0 = 1
         minj = 1
         maxj = 1
         itrij = ikyii + jj
c
c     ----- get information about i-shell and j-shell -----
c
         jshell = 0
         locj = 1
         idum = 0

         call shells_dft(gout,1,ishell,jshell,ishell,jshell,
     &        icd_tag,-1,icd_tag,-1,2,1,idum)

         call ijprim_dft(2)
         if (nij.eq.0) go to 170
c
c     ----- kshell -----
c
               k0 = kst0
               kmc=0
               kk=ii
               kst0 = 1
               ikykk = iky(kk)
               kshell = kk
               mink = mini
               maxk = maxi
               lock = loci
               itrjk = iky(max(jj,kk)) + min(jj,kk)
c
c     ----- lshell ----
c
               maxll = kk
               maxll = jj
               ll=1
               q4 = 1.0d0
c
c     ----- (ii,jj//kk,ll) -----
c
               lshell = 0
               qq4 = q4
C *
C *For shell quartet with no sp shell(s)
C *
c
c     ----- get information about ksh and lsh -----
c
                    idum = 0
                    call shells_dft(gout,2,
     &                   ishell,jshell,ishell,jshell,
     &                   icd_tag,-1,icd_tag,-1,2,1,idum)

                    call genral_dft(gout,2)
                    inn=0
                    do i=mini,maxi
                      jnn=jnn+1
                      inn=inn+1
                      knn=0
                      nn=ijgt(inn)+klgt(inn)
                      rmatrix(jnn)=gout(nn)
                    enddo
C *
C *For shells quartet containing sp shell(s)
C *
c
c
 170       continue
 150     continue
c
c
      return
      end
c
      subroutine te2c_rep_schwarz(icd_tag,gout,rmatrix)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...  computes the 2 centre 2 electron integrals needed for 
c...  applying the Schwarz inequality to 3 centre 2 electron
c...  integrals.
c
***************************************************************
* common blocks used
**************************************************************
**** from main-line code
*********************************************************************
* cslosc * files * infoa *  parallel * parcntl 
* prints * restar * sizes * statis * timez
***************************************************************
* local common blocks
*
* tabinx, junkx, mapperx, flipsx, flip7x, ijlabx
* pqgeomx, qgeomx, pgeomx, geomx, shlg70x, piconx, miscgx, astorex
* ginfx, bshellx, constx, shllfox, maxcx, savecx, typex
* shlinfx
* shltx, miscx, setintx, rootx, denssx, shlnosx, indezx, rtdatx
* inxblkx
********************************************************************
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
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
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)

c
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
ccINCLUDE(common/dft_auxvar)
      integer iwr
      common/dft_iofile/iwr
c
      dimension gout(*)
      dimension rmatrix(*)
c
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      call cpuwal(begin,ebegin)
c ***
c *** establish arrays normally set in main line code
c ***
c
c ... /mapperx/
c
      do 30 i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
 30   continue
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue
c
c ... /auxvarx/
c
c      ofast = .false.
c
c /cslosc/
c
      call setsto(10,0,intcut)
      call setsto(1060,0,intmag)
      nopk = 1
c
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_dft
      nopkr = nopk
      iofrst = iofsym
      nindmx = 1
      call jkint_gamess_schwarz(icd_tag,gout,rmatrix)
      call final_dft
      return
      end
      subroutine jkint_gamess_schwarz(icd_tag,gout,rmatrix)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms * mxprms)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
c
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common /shlg70x/ ishell,jshell,kshell,lshell,
     +                 inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      integer iwr
      common/dft_iofile/iwr
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c                      1296   if d or m shells
c                     10000   if f shells
c                     50625   if g shells
c
c     ----- this version can handle g shells   -----
c
      common/junkx/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      dimension ib(4,4)
      dimension gout(*)
      dimension rmatrix(*)
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
c
c     ----- two-electron integrals -----
c
      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
        igt(loop) = ib(1,ii)
        jgt(loop) = ib(2,ii)
        kgt(loop) = ib(3,ii)
        lgt(loop) = ib(4,ii)
 30   continue
c
c     ----- set some parameters -----
c
c
      imc=0
c
c     l2 = iky(numorb(1))+numorb(1)
c
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     are the rotated axis integrals to be used only?
c
      ogauss = .true.
c
cc      if (intg76.ne.0) call filmax_2c
c
c     ----- ishell -----
c
      jnn=0
      do 150 ii = 1 , nshell(icd_tag)
c
c     ----- print intermediate restart data -----
c
         ishell = ii
         dt0 = time - tim0
         dt1 = time - tim1
         tim1 = time
         ikyii = iky(ii)
c
c     ----- jshell -----
c
         j0 = jst0
         jj=1
         jst0 = 1
         minj = 1
         maxj = 1
         itrij = ikyii + jj
c
c     ----- get information about i-shell and j-shell -----
c
         jshell = 0
         locj = 1

         call shells_dft(gout,1,ishell,jshell,ishell,jshell,
     &        icd_tag,-1,icd_tag,-1,2,1,idum)

         call ijprim_dft(2)
         if (nij.eq.0) go to 170
c
c     ----- kshell -----
c
               k0 = kst0
               kmc=0
               kk=ii
               kst0 = 1
               ikykk = iky(kk)
               kshell = kk
               mink = mini
               maxk = maxi
               lock = loci
               itrjk = iky(max(jj,kk)) + min(jj,kk)
c
c     ----- lshell ----
c
               maxll = kk
               maxll = jj
               ll=1
               q4 = 1.0d0
c
c     ----- (ii,jj//kk,ll) -----
c
               lshell = 0
               qq4 = q4
C *
C *For shell quartet with no sp shell(s)
C *
c
c     ----- get information about ksh and lsh -----
c
                    idum = 0
                    call shells_dft(gout,2,
     &                   ishell,jshell,ishell,jshell,
     &                   icd_tag,-1,icd_tag,-1,2,1,idum)
                    call genral_dft(gout,2)
                    inn=0
                    rmax = 0.0d0
                    do i=mini,maxi
                      jnn=jnn+1
                      inn=inn+1
                      knn=0
                      nn=ijgt(inn)+klgt(inn)
                      rmax=max(rmax,gout(nn))
                    enddo
                    rmatrix(ii)=sqrt(rmax)
C *
C *For shells quartet containing sp shell(s)
C *
c
c
 170       continue
 150     continue
c
c
      return
      end
      subroutine ver_dft_integ_te2c_norm(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/integ_te2c_norm.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
