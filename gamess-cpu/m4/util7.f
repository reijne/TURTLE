c 
c  $Author: mrdj $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util7.m,v $
c  $State: Exp $
c  
      subroutine mpsrtj(q,iq,ncorb,iblki,ifili,ifort)
c
c-------------------------------------------------------------
c     sort out coulomb integrals
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      common/craypk/labs(1360)
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
      dimension q(*),iq(*)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      logical swop1,swop2
c     logical uhf
c     character *8 open
c     data open/'open'/
c
c
c
c     ibl5 = no of integrals in block of sortfile
c
c     uhf = scftyp.eq.open
c
      character*7 fnm
      character*6 snm
      data fnm,snm/'util7.m','mpsrtj'/
c
      ibl5 = nsz340
      iilen = nsz340*mach12/2
      ibl52 = iilen
      call setsto(1360,0,labs)
      ibl5i = lenint(ibl5)
      nij = ncorb*(ncorb+1)/2
      n2 = ncorb*ncorb
c
c       are going to sort the integrals from file ifili and
c       starting block iblki so that for ij (i.ge.j) all kl
c       integrals are available in a square on stream ifort
c
c       maxt is the number of triangles (squares for exchange ints)
c       which can be held in core (allowing n2 wkspace for reading back)
c       which is the number in each bucket
c
      i1  = igmem_alloc_all_inf(maxq,fnm,snm,'i1',IGMEM_DEBUG)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq-n2)/nij
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
      maxa = nbuck*(ibl5+ibl5i) + n2
c
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxq , maxa
         call caserr('stop')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      swop1 = .true.
      swop2 = .true.
      ipss = 1
      call rdsrtj(q(i1),iq(ii1),iblki,ifili,swop1,swop2,ipss)
c
c       read through the sort file to give final result
c
      maxqq = nij*nadd
      call wtsrtj(q(i1),q(i1+n2),maxqq,ncorb,ifort,ipss)
c
      call closbf(0)
      call gmem_free_inf(i1,fnm,snm,'i1')
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine rdsrtj(a,ia,iblki,ifili,swop1,swop2,ipss)
c
c     does the sorting part to get coulomb matrices
c     lower triangles only are produced
c-----------------------------------------------------------------
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
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk1,mkk1,g(5118)
c
      integer ixoxxa, ixoxxb
      integer ispin, iblkzz, iblkz, npassm, intblk, mupblk
      common /uhfspn/ ispin,iblkzz,iblkz,npassm,intblk(4),mupblk(8),
     +      ixoxxa,ixoxxb
c
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/maxlen/maxq
      common/blkin/gin(510),nint
      common/craypk/labs(1360)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      dimension ia(*),a(*)
      logical swop1,swop2
      data lastb/999999/
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c
      if (ipss.le.1) then
         ibl5i = lenint(ibl5)
         do 20 ibuck = 1 , nbuck
            nwbuck(ibuck) = 0
            mark(ibuck) = lastb
            i = (ibuck-1)*(ibl5+ibl5i)
            ibase(ibuck) = i
            ibasen(ibuck) = lenrel(i+ibl5)
 20      continue
c
         call vclr(g,1,nsz340+nsz170)
c
         iblock = 0
      end if
c      ninb = no of elements in bucket(coreload)
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
         call unpack(gin(num2e+1),lab816,labs,numlab)
         do 50 int = 1 , nint
            n4 = int + int + int + int
            i1 = labs(n4-2)
            j1 = labs(n4-3)
            k1 = labs(n4  )
            l1 = labs(n4-1)
            val = gin(int)
            ij = iky(i1) + j1
            kl = iky(k1) + l1
c
c
            if (swop1) then
               iaddr = (ij-1)*nij + kl
c
c     iaddr is address of integral in final sequence
c
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
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
            end if
c
c
c     now check if klij integral is different
c
            if (swop2) then
               iaddr = (kl-1)*nij + ij
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
               nwbuck(ibuck) = nwb
               if (nwb.eq.ibl5) then
c
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
c
            end if
c
c
 50      continue
         go to 30
      end if
      end
      subroutine rdsrtk(a,ia,iblki,ifili,ncoorb)
c-------------------------------------------------------------
c     reads the integral file to produce the back-chained
c     sortfile of exchange integrals
c---------------------------------------------------------------
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
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk1,mkk1,g(5118)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/maxlen/maxq
      common/blkin/gin(510),nint
      common/craypk/labs(1360)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      dimension ia(*),a(*)
      data lastb/999999/
c
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
      ninb = nadd*n2
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
         call unpack(gin(num2e+1),lab816,labs,numlab)
         do 50 int = 1 , nint
            n4 = int + int + int + int
            i1 = labs(n4-2)
            j1 = labs(n4-3)
            k1 = labs(n4  )
            l1 = labs(n4-1)
            val = gin(int)
            i11 = (i1-1)*ncoorb
            j11 = (j1-1)*ncoorb
            k11 = (k1-1)*ncoorb
            l11 = (l1-1)*ncoorb
c----------------------------------------------------------------
c     do the contribution to k([ik],j,l)
c
            ik = iky(i1) + k1
            jl = l11 + j1
c
c
            iaddr = (ik-1)*n2 + jl
c
c     iaddr is address of integral in final sequence
c
            ibuck = (iaddr-1)/ninb
            iaddr = iaddr - ninb*ibuck
            ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
            nwb = nwbuck(ibuck) + 1
            a(ibase(ibuck)+nwb) = val
            ia(ibasen(ibuck)+nwb) = iaddr
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - so empty it
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
c----------------------------------------------------------------
c     if i1=k1 then do k([ik],l,j) - if i1 .ne. k1
c     then this contribution arises when (il/kj) is
c     processed
c
            if (i1.eq.k1 .and. j1.ne.l1) then
               ik = iky(i1) + k1
               jl = j11 + l1
               iaddr = (ik-1)*n2 + jl
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
               nwbuck(ibuck) = nwb
               if (nwb.eq.ibl5) then
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
            end if
c
c
c-----------------------------------------------------------------
c
c     contribution to k([il],j,k)
c
            if (k1.ne.l1) then
c
c       if k1.eq.l1 contribution already covered
c       in k([ik],j,l) above
c
               il = iky(i1) + l1
               jk = k11 + j1
c
               iaddr = (il-1)*n2 + jk
c
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
c
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
               nwbuck(ibuck) = nwb
               if (nwb.eq.ibl5) then
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
            end if
c
c-------------------------------------------------------------------
c     contribution to k([jl],i,k) or k([jl],k,i)
c
            jl = iky(j1) + l1
            ik = k11 + i1
c
            if (j1.lt.l1) then
               jl = iky(l1) + j1
               ik = i11 + k1
            end if
c
            iaddr = (jl-1)*n2 + ik
c
c
            ibuck = (iaddr-1)/ninb
            iaddr = iaddr - ninb*ibuck
            ibuck = ibuck + 1
c
            nwb = nwbuck(ibuck) + 1
            a(ibase(ibuck)+nwb) = val
            ia(ibasen(ibuck)+nwb) = iaddr
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
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
c--------------------------------------------------------------
c     if j1.eq.l1 then include k([jl],k,i)
c     if j1.ne.l1 then this is included when the integral
c     (il/kj) is processed
c
            if (j1.eq.l1 .and. i1.ne.k1) then
               jl = iky(j1) + l1
               ik = i11 + k1
               iaddr = (jl-1)*n2 + ik
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
               nwbuck(ibuck) = nwb
               if (nwb.eq.ibl5) then
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
            end if
c
c-------------------------------------------------------------
c     contribution to k([jk],i,l) or k([jk],l,i)
c
            if (i1.ne.j1) then
c
c       if i1.eq.j1 these contributions covered by
c       earlier terms
c
               jk = iky(j1) + k1
               il = l11 + i1
c
               if (j1.lt.k1) then
                  jk = iky(k1) + j1
                  il = i11 + l1
               end if
c
               iaddr = (jk-1)*n2 + il
c
               ibuck = (iaddr-1)/ninb
               iaddr = iaddr - ninb*ibuck
               ibuck = ibuck + 1
c
c
               nwb = nwbuck(ibuck) + 1
               a(ibase(ibuck)+nwb) = val
               ia(ibasen(ibuck)+nwb) = iaddr
               nwbuck(ibuck) = nwb
               if (nwb.eq.ibl5) then
c
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
c-------------------------------------------------------------
c     if j1.eq.k1 then also include k([jk],l,i)
c     if j1.ne.k1 this comes from integral (ik/jl)
c
               if (j1.eq.k1 .and. i1.ne.l1) then
                  jk = iky(j1) + k1
                  il = i11 + l1
c
                  iaddr = (jk-1)*n2 + il
c
                  ibuck = (iaddr-1)/ninb
                  iaddr = iaddr - ninb*ibuck
                  ibuck = ibuck + 1
c
                  nwb = nwbuck(ibuck) + 1
                  a(ibase(ibuck)+nwb) = val
                  ia(ibasen(ibuck)+nwb) = iaddr
                  nwbuck(ibuck) = nwb
                  if (nwb.eq.ibl5) then
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
               end if
            end if
c
c
 50      continue
         go to 30
      end if
      end
      subroutine mpdag0(ifort1,ifort2,nov,n2,b,c)
      implicit real*8  (a-h,o-z)
c
c     matrix transposition ( middle part of transformation )
c     in-core version
      dimension b(nov*n2),c(nov)
c
      m1 = 1
      call search(m1,ifort1)
c
c     read the matrix into core
c
      ib = 1
      do 20 ixx = 1 , nov
         call rdedz(b(ib),n2,ifort1)
         ib = ib + n2
 20   continue
c
c
c     rewind file
c     ( can be overwriting original )
c
      m1 = 1
      call search(m1,ifort2)
c
c     write out transposed matrix
c     (could cause a few paging problems here ! on a virtual
c      machine )
c
      do 40 iel = 1 , n2
         ipos = iel
         do 30 iij = 1 , nov
            c(iij) = b(ipos)
            ipos = ipos + n2
 30      continue
         call wtedz(c,nov,ifort2)
 40   continue
      return
      end
      subroutine wtsrtj(buf,q,maxqq,ncoorb,ifort,ipss)
c----------------------------------------------------------
c     this reads back down the back-chained sort file
c     produced by edsrd1 to give a final file
c     containing the coulomb matrices arranged
c     sequentially
c----------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension q(maxqq),buf(n2)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk,mkk,g(5118)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      if (ipss.le.1) call rewedz(ifort)
c
      min = 1
      max = nadd
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
            j = 1
            do 30 n = min , max
               call squr(q(j),buf,ncoorb)
               call wtedz(buf,n2,ifort)
               j = j + nij
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
            call unpack(g(nsz341),32,labs,ibl5)
            call dsctr(nkk,g,labs,q)
            go to 20
         end if
 40   continue
      return
      end
      subroutine wtsrtk(q,maxqq,ifort,ipss)
c
c-------------------------------------------------------------
c     reads the back chained sortfile produced by edsrd2
c     to give a sequential file of exchange matrices
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension q(maxqq)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk,mkk,g(5118)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      if (ipss.le.1) call rewedz(ifort)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxqq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     squares min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),n2,ifort)
               j = j + n2
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
            call unpack(g(nsz341),32,labs,ibl5)
            call dsctr(nkk,g,labs,q)
            go to 20
         end if
 40   continue
      return
      end
      subroutine mpsrtk(q,iq,ncoorb,iblki,ifili,ifort,ipss)
c
c -------------------------------------------------------------
c     this produces matrices of exchange integrals
c--------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/craypk/labs(1360)
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
      dimension q(*),iq(*)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
c
c
c     ibl5 = no of integrals in block of sortfile
c     only machine dependent feature should be structure of /bufa/
c
      character*7 fnm
      character*6 snm
      data fnm,snm/'util7.m','mpsrtk'/
c
      ibl5 = nsz340
      iilen = nsz340*lenwrd()/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = ncoorb*(ncoorb+1)/2
      n2 = ncoorb*ncoorb
c
c       are going to sort the integrals from file ifili and
c       starting block iblki so that for ij (i.ge.j) all kl
c       integrals are available in a square on stream ifort
c
c       maxt is the number of triangles (squares for exchange ints)
c       which can be held in core (allowing n2 wkspace for reading back)
c       which is the number in each bucket
c
      i1  = igmem_alloc_all_inf(maxq,fnm,snm,'i1',IGMEM_DEBUG)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq-n2)/n2
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
      maxa = nbuck*(ibl5+ibl5i) + n2
c
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxq , maxa
         call caserr('stop')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      call rdsrtk(q(i1),iq(ii1),iblki,ifili,ncoorb)
c
c       read through the sort file to give final result
c
      maxqq = nadd*n2
      call wtsrtk(q(i1+n2),maxqq,ifort,ipss)
      call closbf(0)
      call gmem_free_inf(i1,fnm,snm,'i1')
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine mpsrt1(q,iq,iblki,ifili,ifort,nrows,ncols,nocc)
c
c -------------------------------------------------------------
c     this produces k(ij) = (ai|bj) - (aj|bi)  for all i,j occ
c     and for a strict lower triangle of ab ( antisymmetric matrix )
c--------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/craypk/labs(1360)
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
      dimension q(*),iq(*)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
c
c
c     ibl5 = no of integrals in block of sortfile
c     only machine dependent feature should be structure of /bufa/
c
      ibl5 = nsz340
      iilen = nsz340*lenwrd()/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = nrows
      n2 = ncols
c
c       are going to sort the integrals from file ifili and
c       starting block iblki so that for ij (i.ge.j) all kl
c       integrals are available in a square on stream ifort
c
c       maxt is the number of triangles (squares for exchange ints)
c       which can be held in core (allowing n2 wkspace for reading back)
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq-n2)/n2
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
      maxa = nbuck*(ibl5+ibl5i) + n2
c
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxq , maxa
         call caserr('stop')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      call rdsrt1(q(i1),iq(ii1),iblki,ifili,nocc)
c
c       read through the sort file to give final result
c
      maxqq = maxq
      ipss = 1
      call wtsrt1(q(i1),maxqq,ifort,ipss)
      call closbf(0)
      call gmem_free(i1)
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine rdsrt1(a,ia,iblki,ifili,nocc)
c-------------------------------------------------------------
c     reads the integral file to produce the back-chained
c     sortfile of elements (ai|bj) and -(aj|bi)
c---------------------------------------------------------------
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
      dimension a(*),ia(*)
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
      common/bufb/nkk1,mkk1,g(5118)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/maxlen/maxq  
      common/blkin/gin(510),nint
      common/craypk/labs(1360)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      data lastb/999999/
c
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
      ninb = nadd*n2
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
         call unpack(gin(num2e+1),lab816,labs,numlab)
         do 50 int = 1 , nint
            n4 = int + int + int + int
            i1 = labs(n4-2)
            j1 = labs(n4-3)
            k1 = labs(n4  )
            l1 = labs(n4-1)
            if (i1.ge.nocc) then
               if (j1.le.nocc) then
                  if (k1.ge.nocc .and. k1.ne.i1) then
                     if (l1.le.nocc) then
                        val = gin(int)
                        l11 = (l1-1)*nocc
                        j11 = (j1-1)*nocc
                        i11 = i1 - nocc
                        k11 = k1 - nocc
c----------------------------------------------------------------
c    do the contribution to k([jl],ab)
c
                        ik = iky(i11) + k11
                        jl = l11 + j1
c
c
                        iaddr = (jl-1)*n2 + ik
c
c     iaddr is address of integral in final sequence
c
                        ibuck = (iaddr-1)/ninb
                        iaddr = iaddr - ninb*ibuck
                        ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                        nwb = nwbuck(ibuck) + 1
                        a(ibase(ibuck)+nwb) = val
                        ia(ibasen(ibuck)+nwb) = iaddr
                        nwbuck(ibuck) = nwb
                        if (nwb.eq.ibl5) then
c
c     this block full - so empty it
c
                           call stopbk
                           mkk1 = mark(ibuck)
                           nkk1 = nwb
                           call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
                           call pack(g(nsz341),32,ia(ibasen(ibuck)+1)
     +                           ,ibl5 )
                           call sttout
                           nwbuck(ibuck) = 0
                           mark(ibuck) = iblock
                           iblock = iblock + nsz
                        end if
c
                        jl = j11 + l1
                        iaddr = (jl-1)*n2 + ik
                        ibuck = (iaddr-1)/ninb
                        iaddr = iaddr - ninb*ibuck
                        ibuck = ibuck + 1
                        nwb = nwbuck(ibuck) + 1
                        a(ibase(ibuck)+nwb) = -val
                        ia(ibasen(ibuck)+nwb) = iaddr
                        nwbuck(ibuck) = nwb
                        if (nwb.eq.ibl5) then
c
c     this block full - so empty it
c
                           call stopbk
                           mkk1 = mark(ibuck)
                           nkk1 = nwb
                           call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
                           call pack(g(nsz341),32,ia(ibasen(ibuck)+1)
     +                               , ibl5)
                           call sttout
                           nwbuck(ibuck) = 0
                           mark(ibuck) = iblock
                           iblock = iblock + nsz
                        end if
                     end if
                  end if
               end if
            end if
c
c
 50      continue
         go to 30
      end if
      end
      subroutine rdsrt2(a,ia,iblki,ifili,nocc1,nocc2,ipss,
     1    noc,nvr)
c-------------------------------------------------------------
c     reads the integral file to produce the back-chained
c     sortfile of elements (ai|bj) all ab for a given ij (mixed)
c---------------------------------------------------------------
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
      common/bufb/nkk1,mkk1,g(5118)
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
      common/maxlen/maxq  
      common/blkin/gin(510),nint
      common/craypk/labs(1360)
      dimension a(*),ia(*)
      data lastb/999999/
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c
      ibl5i = lenint(ibl5)
      tmp1 = 0.0d0
      if (ipss.le.1) then
         do 20 ibuck = 1 , nbuck
            nwbuck(ibuck) = 0
            mark(ibuck) = lastb
            i = (ibuck-1)*(ibl5+ibl5i)
            ibase(ibuck) = i
            ibasen(ibuck) = lenrel(i+ibl5)
 20      continue
c
         call vclr(g,1,nsz340+nsz170)
c
         iblock = 0
      end if
c     ninb = no of elements in bucket(coreload)
      ninb = nadd*n2
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
         call unpack(gin(num2e+1),lab816,labs,numlab)
         do 50 int = 1 , nint
            n4 = int + int + int + int
            i1 = labs(n4-2)
            j1 = labs(n4-3)
            k1 = labs(n4  )
            l1 = labs(n4-1)
            if (i1.gt.nocc1) then
               if (j1.le.nocc1) then
                  if (k1.gt.nocc2) then
                     if (l1.le.nocc2) then
                        if (i1.ne.k1 .or. j1.ne.l1 .or. ipss.ne.2) then
                           if (j1.eq.1 .and. l1.eq.1) then
                           end if
                           val = gin(int)
                           tmp1 = tmp1 + dabs(val)
                           j11 = (j1-1)*noc
                           if (ipss.eq.2) j11 = (l1-1)*noc
                           i11 = (i1-nocc1-1)*nvr
                           if (ipss.eq.2) i11 = (k1-nocc2-1)*nvr
                           k11 = k1 - nocc2
                           if (ipss.eq.2) k11 = i1 - nocc1
                           l11 = l1
                           if (ipss.eq.2) l11 = j1
c----------------------------------------------------------------
c    do the contribution to k([jl],ab)
c
                           ik = k11 + i11
                           jl = j11 + l11
c
c
                           iaddr = (jl-1)*n2 + ik
c
c     iaddr is address of integral in final sequence
c
                           ibuck = (iaddr-1)/ninb
                           iaddr = iaddr - ninb*ibuck
                           ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                           nwb = nwbuck(ibuck) + 1
                           a(ibase(ibuck)+nwb) = val
                           ia(ibasen(ibuck)+nwb) = iaddr
                           nwbuck(ibuck) = nwb
                           if (nwb.eq.ibl5) then
c
c     this block full - so empty it
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
                        end if
                     end if
                  end if
               end if
            end if
c
c
c
 50      continue
         go to 30
      end if
      end
      subroutine mpsrt2(q,iq,iblki,iblki1,ifili,ifort,nrows,ncols,
     1    nocca,noccb,nvirta)
c
c -------------------------------------------------------------
c     this produces k(ij) = (ai|bj) - (aj|bi)  for all i,j occ
c     and for a strict lower triangle of ab ( antisymmetric matrix )
c--------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/craypk/labs(1360)
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
      dimension q(*),iq(*)
c
c     ibl5 = no of integrals in block of sortfile
c     only machine dependent feature should be structure of /bufa/
c
      ibl5 = nsz340
      iilen = nsz340*lenwrd()/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = nrows
      n2 = ncols
c
c       are going to sort the integrals from file ifili and
c       starting block iblki so that for ij (i.ge.j) all kl
c       integrals are available in a square on stream ifort
c
c       maxt is the number of triangles (squares for exchange ints)
c       which can be held in core (allowing n2 wkspace for reading back)
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq)/n2
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
         write (iwr,6010) maxq , maxa
         call caserr('stop')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      ipss = 1
      call rdsrt2(q(i1),iq(ii1),iblki,ifili,noccb,nocca,ipss,
     +            nocca,nvirta)
      ipss = 2
      call rdsrt2(q(i1),iq(ii1),iblki1,ifili,nocca,noccb,ipss,
     +            nocca,nvirta)
c
c       read through the sort file to give final result
c
      maxqq = maxq
      ipss = 1
      call wtsrt2(q(i1),maxqq,ifort,ipss)
      call closbf(0)
      call gmem_free(i1)
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine wtsrt1(q,maxqq,ifort,ipss)
c
c-------------------------------------------------------------
c     reads the back chained sortfile and accumulates to give
c     (ai|bj)-(aj|bi)
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension q(maxqq)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk,mkk,g(5118)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      if (ipss.le.1) call rewedz(ifort)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 50 i = 1 , nbuck
         call vclr(q,1,maxqq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     squares min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),n2,ifort)
               j = j + n2
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
            call unpack(g(nsz341),32,labs,ibl5)
            do 40 iword = 1 , nkk
               q(labs(iword)) = q(labs(iword)) + g(iword)
 40         continue
            go to 20
         end if
 50   continue
      return
      end
      subroutine wtsrt2(q,maxqq,ifort,ipss)
c
c-------------------------------------------------------------
c     reads the back chained sortfile to give (ai|jbar bbar) in
c     rectangles - all a bbar for a given i jbar
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension q(maxqq)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
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
      common/bufb/nkk,mkk,g(5118)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      if (ipss.le.1) call rewedz(ifort)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxqq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     squares min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),n2,ifort)
               j = j + n2
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
            call unpack(g(nsz341),32,labs,ibl5)
            call dsctr(nkk,g,labs,q)
            go to 20
         end if
 40   continue
      return
      end
      subroutine mpsrt0(ifort1,ifort2,ncols,nrows,q,iq,maxq)
      implicit real*8  (a-h,o-z)
c
c     out of core matrix transposition = sorting routine
c     for middle of transformation. uses the usual bucket-sort
c     algorithm
c
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      parameter (maxbuc=1500)
      dimension iq(*),q(maxq)
c
      nov = ncols
      n2 = nrows
c
c     can the matrix transposition be done in core ?
c
      if ((nov+nov*n2).le.maxq) then
         call mpdag0(ifort1,ifort2,nov,n2,q(nov+1),q(1))
         return
      end if
c
c    ibl5 = no of integrals in block of sortfile
c    only machine dependent feature should be structure of /bufb/
c
      ibl5 = nsz340
      iilen = nsz340*lenwrd()/2
      ibl52 = iilen
      ibl5i = lenint(ibl5)
c
c       are going to sort matrix which is (n2,nov) into one
c       which is (nov,n2)
c
c       maxt is the number of columns which can be held in core
c       which is the number in each bucket
c
      maxt = maxq/nov
      niqq = maxq - n2 - lenint(n2+n2)
c
c    warning - some of the areas of the mp2 code use this n2 area
c
      nword = (niqq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
      nbuck = (n2/maxt) + 1
      nadd = min(maxt,n2)
      maxa = nbuck*(ibl5+ibl5i) + n2 + lenint(n2+n2)
c
c     nbuck is the number of buckets actually needed
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxq , maxa
         call caserr('insufficient core')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q,1,maxa)
      i1 = lenrel(n2) + 1
      i2 = i1 + n2
      i4 = i2 + n2
      i3 = n2 + lenint(n2+n2) + 1
c     note that i3 and i4 are the same addresses really
      call setbfa
      call rdsrt0(q,iq(i1),iq(i2),q(i3),iq(i4),ifort1)
c
c       read through the sort file to give final result
c
      maxa = nov*nadd
      call wtsrt0(q(1),maxa,ifort2)
c
      call closbf(0)
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine rdedz0s(r,l,ifile)
      implicit real*8  (a-h,o-z)
c     this is the original routine that allows 16 bits
c     for each ij matrix element. This is insufficient
c     for large basis sets, empirically gt 170+ GTOs
      dimension r(l)
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
      parameter (m408=408,m102=102,m16=16)
      common/zbuf/z(m408),iiz(m102),b,n,nxtr
      common/bufeds/iz(m408)
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
c
      call vclr(r,1,l)
      call search(ipos(ifile),ifile)
2     call find(ifile)
      call get(z,nw)
      call unpack(iiz,m16,iz,m408)
      call upack2(b,n,nxtr)
      call dsctr(n,z,iz,r)
      if(nxtr.eq.1) go to 2
      return
      end
      subroutine rdeds0s(r,ir,nzero,ifile)
      implicit real*8  (a-h,o-z)
c     this is the original routine that allows 16 bits
c     for each ij matrix element. This is insufficient
c     for large basis sets, empirically gt 170+ GTOs
      integer b
      parameter (m408=408,m102=102,m16=16)
      common/zbuf/z(m408),iiz(m102),b,n,nxtr
      common/zlabs/iz(m408)
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
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
      dimension r(*),ir(*)
c
      call search(ipos(ifile),ifile)
2     call find(ifile)
      call get(z,nw)
      call unpack(iiz,m16,iz,m408)
      call upack2(b,n,nxtr)
      do 100 i=1,n
        nzero=nzero+1
        ir(nzero)=iz(i)
        r(nzero)=z(i)
100   continue
      if(nxtr.eq.1) go to 2
      return
      end
      subroutine rdedz0(r,l,ifile)
      implicit real*8  (a-h,o-z)
      dimension r(l)
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
      parameter (m340=340,m170=170,m32=32)
      common/zbuf/z(m340),iiz(m170),b,n,nxtr
      common/bufeds/iz(m340)
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
c
      call vclr(r,1,l)
      call search(ipos(ifile),ifile)
2     call find(ifile)
      call get(z,nw)
      call unpack(iiz,m32,iz,m340)
      call upack2(b,n,nxtr)
      call dsctr(n,z,iz,r)
      if(nxtr.eq.1) go to 2
      return
      end
      subroutine rdeds0(r,ir,nzero,ifile)
      implicit real*8  (a-h,o-z)
      integer b
      parameter (m340=340,m170=170,m32=32)
      common/zbuf/z(m340),iiz(m170),b,n,nxtr
      common/zlabs/iz(m340)
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
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
      dimension r(*),ir(*)
c
      call search(ipos(ifile),ifile)
2     call find(ifile)
      call get(z,nw)
      call unpack(iiz,m32,iz,m340)
      call upack2(b,n,nxtr)
      do 100 i=1,n
        nzero=nzero+1
        ir(nzero)=iz(i)
        r(nzero)=z(i)
100   continue
      if(nxtr.eq.1) go to 2
      return
      end
      subroutine rdsrt0(a,ia,ib,aa,iaa,ifort)
      implicit real*8  (a-h,o-z)
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
      common/bufb/nkk,mkk,g(5118)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
      dimension a(*),ia(*),ib(*),aa(*),iaa(*)
      data lastb/999999/
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c       the first n2 real words of core are used as input
c       workspace.
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
c
      m1 = 1
      call search(m1,ifort)
      iblock = 0
c     ninb = no of elements ib bucket(coreload)
      ninb = nadd*nov
c
      do 50 i = 1 , nov
c
c     read column i of original matrix
c
         nrow = n2
         call rdeds(a,ia,nrow,ifort)
c
         do 30 j = 1 , nrow
c     iaddr is address of (i,j) in transposed matrix
            iaddr = i + (ia(j)-1)*nov
            ibuck = (iaddr-1)/ninb
            ia(j) = iaddr - ninb*ibuck
            ib(j) = ibuck + 1
 30      continue
c
c
         do 40 j = 1 , nrow
c
c
            ibuck = ib(j)
c
c     element goes in bucket ibuck with modified address
c
            nwb = nwbuck(ibuck) + 1
            aa(ibase(ibuck)+nwb) = a(j)
            iaa(ibasen(ibuck)+nwb) = ia(j)
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - empty
c
               call stopbk
               mkk = mark(ibuck)
               nkk = nwb
               call dcopy(ibl5,aa(ibase(ibuck)+1),1,g,1)
               call pack(g(nsz341),32,iaa(ibasen(ibuck)+1),ibl5)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 40      continue
 50   continue
c
c      empty anything remaining in buckets
c
      do 60 ibuck = 1 , nbuck
         nwb = nwbuck(ibuck)
         if (nwb.ne.0) then
            call stopbk
            mkk = mark(ibuck)
            nkk = nwb
            call dcopy(ibl5,aa(ibase(ibuck)+1),1,g,1)
            call pack(g(nsz341),32,iaa(ibasen(ibuck)+1),ibl5)
            call sttout
            nwbuck(ibuck) = 0
            mark(ibuck) = iblock
            iblock = iblock + nsz
         end if
 60   continue
c
c
      call stopbk
      return
      end
      subroutine rdedz(r,l,ifile)
c===========================================================
c     rdedz and wtedz are sequential read and write routines
c     which take out the zeros of the array written.
c     They include buffering, attempts at vectorisation of
c     tests, gather / scatter operations, and allow for
c     different label packing schemes
c==========================================================
      implicit real*8  (a-h,o-z)
      dimension r(l)
c     maxl (1) must not exceed 65536
c          (2) must be divisable by 4
      parameter (maxl=65000)
      if(l.le.maxl)then
         call rdedz0(r,l,ifile)
      else
         npass = (l-1)/maxl + 1
         iremai=l-maxl*(npass-1)
         ioff=1
         do 10 i=1,npass
            len=maxl
            if(i.eq.npass)len=iremai
            call rdedz0(r(ioff),len,ifile)
            ioff=ioff+maxl
10       continue
      endif
      return
      end
      subroutine rdeds(r,ir,l,ifile)
      implicit real*8  (a-h,o-z)
c
c     reads a file written by wtedz, but returns
c     packed array, and its index instead
c     of inserting zeros.
c     nb the value of l supplied as argument should be
c     the maximum length of the array, but the value
c     on return will be the number of non-zero elements
c
      dimension r(l),ir(l)
      parameter (maxl=65000)
      nzero=0
      if(l.le.maxl)then
         call rdeds0(r,ir,nzero,ifile)
      else
         npass = (l-1)/maxl + 1
         ioff=1
         do 10 i=1,npass
           call rdeds0(r(ioff),ir(ioff),nzero,ifile)
           ioff=ioff+nzero
10       continue
      endif
      l=nzero
      return
      end
      subroutine rewedz(ifile)
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
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
      common/crio/last(40)
      last(ifile) = 1
      call search(1,ifile)
      return
      end
      subroutine wtedz0s(r,l,ifile)
      implicit real*8  (a-h,o-z)
c     this is the original routine that allows 16 bits
c     for each ij matrix element. This is insufficient
c     for large basis sets, empirically gt 170+ GTOs
      dimension r(l)
c
      real*8 cutomp
      common /mpcut/ cutomp
c
      parameter (maxl=65000)
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
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
      common/bufeds/map(maxl)
      integer b, pack2
      parameter (m408=408,m102=102,mword=511,m16=16)
      common/zbuf/z(m408),iiz(m102),b,num,nxtr
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
c
      call search(ipos(ifile),ifile)
c
c     set all small elements to be exactly zero
c
      do 10 i=1,l
        if(dabs(r(i)).le.cutomp)r(i)=0.0d0
10    continue
c
c     find map of non zero elements
      call dlstne(l,r,1,0.0d0,n,map)
c
c     write out the non zero elements in batchs - size
c     of batch controlled by parameter statements
c
      nblks=(n-1)/m408
      nremai=n-nblks*m408
      ioff=0
      do 20 j=1,nblks
      call dgthr(m408,r,z,map(ioff+1))
      call pack(iiz,m16,map(ioff+1),m408)
      nxtr=1
      num=m408
      b=pack2(num,nxtr)
      call put(z,mword,ifile)
      ioff=ioff+m408
20    continue
c
c     last batch
c
      call dgthr(nremai,r,z,map(ioff+1))
      call pack(iiz,m16,map(ioff+1),m408)
      nxtr=0
      num=nremai
      b=pack2(num,nxtr)
      call put(z,mword,ifile)
      return
      end
      subroutine wtedz0(r,l,ifile)
      implicit real*8  (a-h,o-z)
      dimension r(l)
c
      real*8 cutomp
      common /mpcut/ cutomp
c
      parameter (maxl=65000)
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
      common/disc/isel,npb,npbsz,irep,ichek,ipos(maxlfn),
     1  iblksz(maxlfn)
      common/bufeds/map(maxl)
      integer b, pack2
      parameter (m340=340,m170=170,mword=511,m32=32)
      common/zbuf/z(m340),iiz(m170),b,num,nxtr
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
c
      call search(ipos(ifile),ifile)
c
c     set all small elements to be exactly zero
c
      do 10 i=1,l
        if(dabs(r(i)).le.cutomp)r(i)=0.0d0
10    continue
c
c     find map of non zero elements
      call dlstne(l,r,1,0.0d0,n,map)
c
c     write out the non zero elements in batchs - size
c     of batch controlled by parameter statements
c
      nblks=(n-1)/m340
      nremai=n-nblks*m340
      ioff=0
      do 20 j=1,nblks
      call dgthr(m340,r,z,map(ioff+1))
      call pack(iiz,m32,map(ioff+1),m340)
      nxtr=1
      num=m340
      b=pack2(num,nxtr)
      call put(z,mword,ifile)
      ioff=ioff+m340
20    continue
c
c     last batch
c
      call dgthr(nremai,r,z,map(ioff+1))
      call pack(iiz,m32,map(ioff+1),m340)
      nxtr=0
      num=nremai
      b=pack2(num,nxtr)
      call put(z,mword,ifile)
      return
      end
      subroutine wtedz(r,l,ifile)
      implicit real*8  (a-h,o-z)
      dimension r(l)
c     maxl (1) must not exceed 65536
c          (2) must be divisable by 4
      parameter (maxl=65000)
      if(l.le.maxl)then
         call wtedz0(r,l,ifile)
      else
         npass = (l-1)/maxl + 1
         iremai=l-maxl*(npass-1)
         ioff=1
         do 10 i=1,npass
            len=maxl
            if(i.eq.npass)len=iremai
            call wtedz0(r(ioff),len,ifile)
            ioff=ioff+maxl
10       continue
      endif
      return
      end
      subroutine wtsrt0(q,maxq,ifort)
      implicit real*8  (a-h,o-z)
      dimension q(maxq)
      common/bufb/nkk,mkk,g(5118)
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
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/sortpk/labs(1)
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
c     ibl5i = lenint(ibl5)
c
      m1 = 1
      call search(m1,ifort)
c
      min = 1
      max = nadd
      call stopbk
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         mkk = mark(i)
         call vclr(q,1,maxq)
 20      if (mkk.eq.lastb) then
c
c     coulumns min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),nov,ifort)
               j = j + nov
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.n2) max = n2
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
            call unpack(g(nsz341),32,labs,ibl5)
            call dsctr(nkk,g,labs,q)
            go to 20
         end if
 40   continue
      return
      end
      subroutine ver_util7(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util7.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
