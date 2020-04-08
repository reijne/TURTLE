c     deck=integc
c ******************************************************
c ******************************************************
c             =   int2esp  =
c ******************************************************
c ******************************************************
      subroutine pkin70(q,iso,gout,nshels,outvv)
      implicit real*8  (a-h,p-w),integer (i-n),logical    (o)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
      character *8 fnm
      character *6 snm
c
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      dimension norgin(3)
      data fnm/'integc.m'/
      data snm/'pkin70'/
      data m25/25/
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data norgin/0,256,512/
c
c     ----- set some parameters -----
c
      inull = igmem_null()
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      ltri = ikyp(nshels)
      oskipp = .false.
      call sinset
      if(opank) osortp = .false.
c
c     ----- allocate core memory
c
      dlncutoff = dlog(cutoff)
      if (oschw)  then
        ischw = igmem_alloc_inf(ltri,fnm,snm,'schwarz',IGMEM_DEBUG)
        nschwz = 0
      endif
c
      if (osortp) then
c
c *** allocate core for symmetry matrices and buffers for p sort
c *** num*maxshellsize buffers length blocksize
c
        ngbf = num * 4
c ***
c *** need mods to clrgbf to use anything other than 340 here
c
        nav = lenwrd()
        lengbf = num2ep
        igbf = igmem_alloc_inf(ngbf*lengbf,fnm,snm,'gbf',IGMEM_DEBUG)
        itmp = (ngbf*lengbf+1)/nav
        iklbf = igmem_alloc_inf(itmp,fnm,snm,'klbf',IGMEM_DEBUG)
        itmp = (ngbf+1)/nav
        iiptbf = igmem_alloc_inf(itmp,fnm,snm,'ptbf',IGMEM_DEBUG)
c *** zero out p sort buffer counts
        call setsto(ngbf,0,q(iiptbf))
      endif
c
      if(oschw) then
c *** read in ints for schwarz inequality test
         call secget(isect(421),m25,iblk25)
         call rdedx(q(ischw),ltri,iblk25,idaf)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c     read in error function interpolation table
c
c     fill common maxc ... used to discard small (less than 10**-6)
c     integrals before they are fully evaluated
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
      call filmax
c
c
c     ----- ishell -----
c
      do 920 ii = ist0,nshels
c
c *** set ijbase to triangle of first basis function in shell
c *** and check integrity of counters
      if (osortp) then
        ijbase = (kloc(ii)*(kloc(ii)-1))/2
        if(icount.ne.1) call caserr('icount.ne.1 for new ii')
      endif
      if(kad(ii))9201,921,921
c
c     ----- print intermediate restart data -----
c
 921  dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv)write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
c
c     ----- eliminate ishell -----
c
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 9201
      mi(it) = id
  120 continue
c
c     ----- jshell -----
c
      j0 = jst0
      do 900 jj = j0,ii
      jst0 = 1
      if(kad(jj))900,901,901
 901  do 200 it = 1,nt
      id = mi(it)
      jd = iso(jj,it)
      mj(it) = jd
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id-ii) 200,180,900
  180 if (jd-jj) 200,200,900
  200 continue
c
c     ----- kshell -----
c
      k0 = kst0
      do 880 kk = k0,jj
      kst0 = 1
      if(kad(kk))880,881,881
 881  do 340 it = 1,nt
      id = mi(it)
      jd = mj(it)
      kd = iso(kk,it)
      mk(it) = kd
  240 if (id .ge. jd) go to 260
      nd = id
      id = jd
      jd = nd
  260 if (jd .ge. kd) go to 280
      nd = jd
      jd = kd
      kd = nd
      go to 240
  280 if (id-ii) 340,300,880
  300 if (jd-jj) 340,320,880
  320 if (kd-kk) 340,340,880
  340 continue
c
c     ----- lshell ----
c
      l0 = lst0
      do 860 ll = l0,kk
      lst0 = 1
      if(kad(ll))860,861,861
 861  n4 = 0
      do 540 it = 1,nt
      id = mi(it)
      jd = mj(it)
      kd = mk(it)
      ld = iso(ll,it)
  380 if (id .ge. jd) go to 400
      nd = id
      id = jd
      jd = nd
  400 if (jd .ge. kd) go to 420
      nd = jd
      jd = kd
      kd = nd
      go to 380
  420 if (kd .ge. ld) go to 440
      nd = kd
      kd = ld
      ld = nd
      go to 400
  440 if (id-ii) 540,460,860
  460 if (jd-jj) 540,480,860
  480 if (kd-kk) 540,500,860
  500 if (ld-ll) 540,520,860
  520 n4 = n4+1
      m0(n4) = it
  540 continue
c
c     ----- check for redundancies between the 3 combinations-
c           (ii,jj//kk,ll),(ii,kk//jj,ll),(ii,ll//jj,kk)
c
      oskpa = jj .eq. kk
      oskpb = (ii .eq. kk) .or. (jj .eq. ll)
      oskpc = (ii .eq. jj) .or. (kk .eq. ll)
      onpsym = .false.
      if (oskpa .or. oskpb .or. oskpc) go to 720
      onpsym = .true.
      do 640 m = 1,n4
      it = m0(m)
      ih = mi(it)
      jh = mj(it)
      if (jh .gt. ih) go to 560
      id = ih
      jd = jh
      go to 580
  560 id = jh
      jd = ih
  580 if ( .not. oskpa) oskpa = (id .eq. ii .and. jd .eq. kk) .or. (id
     +     .eq. jj .and. jd .eq. ll)
      if ( .not. oskpb) oskpb = (id .eq. ii .and. jd .eq. ll) .or. (id
     +     .eq. jj .and. jd .eq. kk)
      if (oskpa .and. oskpb) go to 660
      kh = mk(it)
      if (kh .gt. ih) go to 600
      id = ih
      kd = kh
      go to 620
  600 id = kh
      kd = ih
  620 if ( .not. oskpc) oskpc = (id .eq. ii .and. kd .eq. ll) .or. (id
     +     .eq. jj .and. kd .eq. kk)
      if (oskpa .and. oskpc) go to 680
      if (oskpb .and. oskpc) go to 700
  640 continue
      go to 720
  660 oskpc = .true.
      go to 720
  680 oskpb = .true.
      go to 720
  700 oskpa = .true.
  720 continue
      q4 =  dble(nt)/ dble(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
      iexch = 1
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
      if (oskpa .and. onpsym) qq4 = qq4+q4
      if (oskpb .and. onpsym) qq4 = qq4+q4
      go to 780
c
c     ----- (ii,kk//jj,ll) -----
c
  740 if (oskpa) go to 760
      iexch = 2
      ishell = ii
      jshell = kk
      kshell = jj
      lshell = ll
      qq4 = q4
      if (oskpc .and. onpsym) qq4 = qq4+q4
      go to 780
c
c     ----- (ii,ll//jj,kk) -----
c
  760 if (oskpb .or. oskpc) go to 840
      iexch = 3
      ishell = ii
      jshell = ll
      kshell = jj
      lshell = kk
      qq4 = q4
  780 continue
c
c     ----- initialize gout to zero -----
c
      if (oschw) then
       ijij = iky(ishell) + jshell + ischw -1
       klkl = iky(kshell) + lshell + ischw -1
       test = q(ijij) + q(klkl)
       oskipp = test.lt.dlncutoff
       if(oskipp) nschwz = nschwz + 1
      endif
      norg = norgin(iexch)
      call genr70(gout(norg+1),iexch,oskipp)
      go to (740,760,840),iexch
  840 continue
      if(opank)
     *         call pkfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,gout)
      if(.not.opank)then
        if(.not.onpsym. or. oskpa.or.oskpb.or.oskpc) then
          call pfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,gout,q)
        else
          call pfi70s(ii,jj,kk,ll,gout,q)
        endif
      endif
c
c     ----- check cpu time/ maxblock condition -----
c
      call chkout(ii,jj,kk,ll,q(inull),q)
      if(omaxb.or.tim.gt.timlim)then
         go to 940
      endif
  860 continue
  880 continue
  900 continue
      time = cpulft(1)
 9201 continue
c *** clear p sort buffers
      if (osortp) then
        call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
        do 9202 ibuf = 1,ngbf
9202        call clrgbf(ibuf,q(igbf),q(iklbf),q(iiptbf))
        call blocki
      endif
c ***
  920 continue
      call final(q,q(inull),q(inull))
      if(oschw.and.outvv) write(iwr,9010) nschwz
c
c     ----- reset core memory
c
 940  if (osortp) then
        call gmem_free_inf(iiptbf,fnm,snm,'ptbf')
        call gmem_free_inf(iklbf,fnm,snm,'klbf')
        call gmem_free_inf(igbf,fnm,snm,'gbf')
      endif
      if (oschw)  then
        call gmem_free_inf(ischw,fnm,snm,'schwarz')
      endif
c
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9010 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine genr70(gout,iexch,oskips)
c
c     gaussian two electron integral package
c     shell evaluates those integrals involving only s and p functions.
c     main loop over shells ...  accepts numbers of four shells
c     ishell  jshell  kshell  lshell
c     finds their angular quantum numbers
c     and
c     based on this orders shells in a standard manner
c     inew  jnew  knew  lnew
c     only possibilities allowed for angular quantum numbers are then
c     0000  0001  0011  0101  0111  1111
c     determines type of integral set based on the above numbers
c     calls the following routines in the order given
c     first time to preset output routines
c     filmax
c     to preset integral accurcy limits
c     sinfo
c     obtains geometrical information about the four centers
c     finds two sets of local axes
c     for centers
c     a and b  p set
c     c and d  q set
c     pinf
c     obtains information about gaussian functions connected with the p
c     set of axes
c     at this point
c     shell obtains information about the gaussian functions connected
c     with the q set of axes
c     sp0000 to sp1111
c     obtains up to 88 integrals referred to axes a b and q
c     rot2
c     rotates these integrals to up to 160 integrals on a b and q
c     tq0011 to tq1111
c     translates these integrals on a b and q to up to 256 integrals on
c     a b c and d
c     r30001 to r31111
c     rotates up to 256 integrals on a b c and d to the same number
c     referred to the fixed space axes
c
c     oskips (schwarz skip) means only ib indexing, and zeroing is done.
c
c
c
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
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
c
      integer ib
      common /flip70/ ib(4,3)
c
      common/lt/lat,lbt,lct,ldt
      common/type/itype,jtype
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
      common/miscg/mab,mcd,ngangb
      dimension gout(*)
c
      lat = ktype(ishell)-1
      lbt = ktype(jshell)-1
      lct = ktype(kshell)-1
      ldt = ktype(lshell)-1
      itype = 8*lat+4*lbt+2*lct+ldt+1
      go to (120,120,160,120,180,120,160,120,220,140,200,140,180,180,
     +     240,120),itype
c     types 0000,0001,0101,0011,0111,1111 are unaltered
  120 inew = ishell
      jnew = jshell
      knew = kshell
      lnew = lshell
      la = lat
      lb = lbt
      lc = lct
      ld = ldt
      ib(1,iexch) = 1
      ib(2,iexch) = 2
      ib(3,iexch) = 3
      ib(4,iexch) = 4
      go to 260
c     types 1001,1011 have ij switched
  140 inew = jshell
      jnew = ishell
      knew = kshell
      lnew = lshell
      la = lbt
      lb = lat
      lc = lct
      ld = ldt
      ib(1,iexch) = 2
      ib(2,iexch) = 1
      ib(3,iexch) = 3
      ib(4,iexch) = 4
      go to 260
c     types 0010,0110 have kl switched
  160 inew = ishell
      jnew = jshell
      knew = lshell
      lnew = kshell
      la = lat
      lb = lbt
      lc = ldt
      ld = lct
      ib(1,iexch) = 1
      ib(2,iexch) = 2
      ib(3,iexch) = 4
      ib(4,iexch) = 3
      go to 260
c     types 0100,1100,1101 have pairs ij and kl switched
  180 inew = kshell
      jnew = lshell
      knew = ishell
      lnew = jshell
      la = lct
      lb = ldt
      lc = lat
      ld = lbt
      ib(1,iexch) = 3
      ib(2,iexch) = 4
      ib(3,iexch) = 1
      ib(4,iexch) = 2
      go to 260
c     type 1010 has ij switched and kl switched
  200 inew = jshell
      jnew = ishell
      knew = lshell
      lnew = kshell
      la = lbt
      lb = lat
      lc = ldt
      ld = lct
      ib(1,iexch) = 2
      ib(2,iexch) = 1
      ib(3,iexch) = 4
      ib(4,iexch) = 3
      go to 260
c     type 1000  has pairs ij and kl switched followed by kl switch
  220 inew = kshell
      jnew = lshell
      knew = jshell
      lnew = ishell
      la = lct
      lb = ldt
      lc = lbt
      ld = lat
      ib(1,iexch) = 4
      ib(2,iexch) = 3
      ib(3,iexch) = 1
      ib(4,iexch) = 2
      go to 260
c     type 1110 has pairs ij and kl switched followed by ij switch
  240 inew = lshell
      jnew = kshell
      knew = ishell
      lnew = jshell
      la = ldt
      lb = lct
      lc = lat
      ld = lbt
      ib(1,iexch) = 3
      ib(2,iexch) = 4
      ib(3,iexch) = 2
      ib(4,iexch) = 1
  260 continue
c     only 6 standard types remain. 0000,0001,0011,0101,0111,1111
c     specify these by jtype
      go to (280,300,300,320,300,340,340,360,300,340,340,360,320,360,
     +     360,380),itype
  280 jtype = 1
      ngout = 1
      go to 400
  300 jtype = 2
      ngout = 4
      go to 400
  320 jtype = 3
      ngout = 16
      go to 400
  340 jtype = 4
      ngout = 64
      go to 400
  360 jtype = 5
      ngout = 64
      go to 400
  380 jtype = 6
      ngout = 256
  400 continue
c     empty common gout
      call vclr(gout,1,ngout)
      if(oskips) return
      call sinfo
      go to(440,480,500,520,540,560),jtype
 440  call sp0000(gout)
      go to 1080
 480  call sp0001(gout)
      go to 1080
 500  call sp0011(gout)
      go to 1080
 520  call sp0101(gout)
      go to 1080
 540  call sp0111(gout)
      go to 1080
 560  call sp1111(gout)
      go to 1080
1080  return
      end
      subroutine sinfo
c
c     --------------------------
c     --------------------------
c
c     obtains information about shells inew,knew,jnew,lnew
c     coordinates of abcd go into common geom
c     number of gaussians go into nga,... in common shlinf
c     shell angular quantum numbers la,... go into common shlinf
c     gaussian exponents go into arrays ag,bg,cg,dg in common shlinf
c     gaussian coefficients go into arrays csa,cpa,... in common shlinf
c
      implicit real*8  (a-h,p-w),integer (i-n),logical    (o)
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
      common/miscg/mab,mcd,ngangb
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
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
c
      real*8 p, q, r
      common /bshell/ p(mxshel),q(mxshel),r(mxshel)
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
      real*8 conp
c
c     conp = prefactors from pairs of primitives
c
      common /const/ conp(mxprms*mxprms)
c
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
      common/type/itype,jtype
      common /saveco/ uuuuuu,idadap(2),ijato,klato,ijshlo
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      data sixty/60.0d0/
      data dzero/0.0d0/,pt5/0.5d0/,pt7/0.7d0/,pt9/0.9d0/,done/1.0d0/
      data pt0001/1.0d-04/,tenm12/1.0d-12/
c
c     starting locations of shells inew jnew knew and lnew in list
c     of gaussian functions
c...   1024 changed to mxshel; seems just primitive shift only used here
c...   jvl may 2001
      ijatom=katom(inew)*mxshel+katom(jnew)
      klatom=katom(knew)*mxshel+katom(lnew)
      ijshel=inew*mxshel+jnew
      if( ijshel.ne.ijshlo ) then
        i = kstart(inew)
        j = kstart(jnew)
c     numbers of gaussian functions in shells inew jnew knew and lnew
        nga = kng(inew)
        ngb = kng(jnew)
        ngangb = nga*ngb
        mab = la+lb-1
        do 100 ni = 1,nga
        n = i-1+ni
c     maximum coefficient associated with shell
c     used to determine if any of the integrals associated with a set
c     of shells is lagge enough to warrant evaluation of the entire set
        cmaxa(ni) = cmax(n)
c     gaussian exponents
        ag(ni) = ex(n)
c     s coefficients
        csa(ni) = cs(n)
c     p coefficients
  100   cpa(ni) = cp(n)
c     repeat procedure for shells jnew knew and lnew
        do 120 nj = 1,ngb
        n = j-1+nj
        cmaxb(nj) = cmax(n)
        bg(nj) = ex(n)
        csb(nj) = cs(n)
  120   cpb(nj) = cp(n)
      endif
      k = kstart(knew)
      ngc = kng(knew)
      do 140 nk = 1,ngc
      n = k-1+nk
      cmaxc(nk) = cmax(n)
      cgg(nk) = ex(n)
      csc(nk) = cs(n)*qq4
140   cpc(nk) = cp(n)*qq4
      l = kstart(lnew)
      ngd = kng(lnew)
      do 160 nl = 1,ngd
      n = l-1+nl
      cmaxd(nl) = cmax(n)
      dg(nl) = ex(n)
      csd(nl) = cs(n)
  160 cpd(nl) = cp(n)
c     fill common misc
      mcd = lc+ld-1
c
      oijat = ijatom.eq.ijato
      oklat = klatom.eq.klato
      if(  oijat.and.oklat ) goto 10001
      if( .not. oijat) then
        ax = p(inew)
        ay = q(inew)
        az = r(inew)
        bx = p(jnew)
        by = q(jnew)
        bz = r(jnew)
        abx = bx-ax
        aby = by-ay
        abz = bz-az
        rabsq = abx**2+aby**2+abz**2
        rab = dsqrt(rabsq)
        if (rab) 1200,1400,1200
 1200   p31 = abx/rab
        p32 = aby/rab
        p33 = abz/rab
        go to 1600
 1400   p31 = dzero
        p32 = dzero
        p33 = done
1600    continue
        ijato=ijatom
      endif
      if( .not. oklat ) then
        cx = p(knew)
        cy = q(knew)
        cz = r(knew)
        dx = p(lnew)
        dy = q(lnew)
        dz = r(lnew)
        cdx = dx-cx
        cdy = dy-cy
        cdz = dz-cz
        rcdsq = cdx**2+cdy**2+cdz**2
        rcd = dsqrt(rcdsq)
        if (rcd) 1800,2000,1800
 1800   q31 = cdx/rcd
        q32 = cdy/rcd
        q33 = cdz/rcd
        go to 2200
 2000   q31 = dzero
        q32 = dzero
        q33 = done
2200    continue
        klato=klatom
      endif
      cosg = p31*q31+p32*q32+p33*q33
      cosg = dmin1(done,cosg)
      cosg = dmax1(-done,cosg)
      sing = dsqrt(done-cosg*cosg)
      if (  dabs(cosg)-pt9) 300,300,240
  240 ppq1 = p31+q31
      ppq2 = p32+q32
      ppq3 = p33+q33
      pmq1 = p31-q31
      pmq2 = p32-q32
      pmq3 = p33-q33
      p21 = pmq2*ppq3-ppq2*pmq3
      p22 = pmq3*ppq1-ppq3*pmq1
      p23 = pmq1*ppq2-ppq1*pmq2
      p2 = dsqrt(p21*p21+p22*p22+p23*p23)
      sing = pt5*p2
      if (sing-tenm12) 280,260,260
  260 temp = done/p2
      p21 = p21*temp
      p22 = p22*temp
      p23 = p23*temp
      go to 380
  280 if (  dabs(p31)-pt7) 360,360,320
  300 p21 = (p32*q33-p33*q32)/sing
      p22 = (p33*q31-p31*q33)/sing
      p23 = (p31*q32-p32*q31)/sing
      go to 380
  320 p3333 = p33*p33
      p3333 = dmin1(done,p3333)
      sinp = dsqrt(done-p3333)
      p21 = p32/sinp
      p22 = -p31/sinp
      p23 = dzero
      go to 380
  360 p3131 = p31*p31
      p3131 = dmin1(done,p3131)
      sinp = dsqrt(done-p3131)
      p21 = dzero
      p22 = p33/sinp
      p23 = -p32/sinp
  380 q21 = p21
      q22 = p22
      q23 = p23
      p11 = p22*p33-p23*p32
      p12 = p23*p31-p21*p33
      p13 = p21*p32-p22*p31
      q11 = q22*q33-q23*q32
      q12 = q23*q31-q21*q33
      q13 = q21*q32-q22*q31
      acx = (cx-ax)*p11+(cy-ay)*p12+(cz-az)*p13
      acy = (cx-ax)*p21+(cy-ay)*p22+(cz-az)*p23
      if (  dabs(acy)-pt0001) 400,400,420
  400 acy = dzero
  420 continue
      acz = (cx-ax)*p31+(cy-ay)*p32+(cz-az)*p33
      acy2 = acy*acy
c
10001 continue
      if( ijshel.ne.ijshlo ) then
        ijshlo=ijshel
        ind=0
        do 38000 i = 1,nga
        ga = ag(i)
        csai = csa(i)
        cpai = cpa(i)
        do 36000 j = 1,ngb
        ind = ind+1
        gb = bg(j)
        gab = ga+gb
        gp(ind) = gab
        eab = done/gab
        ep(ind) = eab
        gbeab = gb*eab
        app(ind) = gbeab*rab
        bpp(ind) = app(ind)-rab
        qqq = ga*gbeab*rabsq
        if (qqq-sixty) 18000,18000,10000
10000   continue
        ismlp(ind) = 2
        dp00p(ind) = dzero
        if (jtype-3) 36000,36000,12000
12000   dp01p(ind) = dzero
        conp(ind) = dzero
        if (jtype-5) 14000,14000,16000
14000   bpp(ind) = bpp(ind)*gab
        go to 36000
16000   dp10p(ind) = dzero
        dp11p(ind) = dzero
        go to 36000
18000   qq =  dexp(-qqq)*eab
        qqtest = cmaxa(i)*cmaxb(j)*qq
        if (qqtest-error1) 22000,22000,20000
20000   ismlp(ind) = 0
        go to 28000
22000   if (qqtest-error2) 26000,26000,24000
24000   ismlp(ind) = 1
        go to 28000
26000   continue
        ismlp(ind) = 2
28000   qqqq = pito52*qq
        dp00p(ind) = qqqq*csai*csb(j)
        if (jtype-3) 36000,36000,30000
30000   dp01p(ind) = qqqq*csai*cpb(j)
        if (jtype-5) 32000,32000,34000
32000   conp(ind) = dp01p(ind)*eab
        dp00p(ind) = dp00p(ind)*gab/dp01p(ind)
        bpp(ind) = bpp(ind)*gab
        go to 36000
34000   dp10p(ind) = qqqq*cpai*csb(j)
        dp11p(ind) = qqqq*cpai*cpb(j)
        conp(ind) = dp11p(ind)
        dp00p(ind) = dp00p(ind)/dp11p(ind)
        dp01p(ind) = dp01p(ind)/dp11p(ind)
        dp10p(ind) = dp10p(ind)/dp11p(ind)
36000   continue
38000   continue
      endif
      return
      end
      subroutine sinset
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common /saveco/ uuuuuu,idadap(2),ijato,klato,ijshlo
      ijato=-1
      klato=-1
      ijshlo=-1
      return
      end
      subroutine ver_integc(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integc.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
