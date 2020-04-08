c     deck=integc
c ******************************************************
c ******************************************************
c             =   int2esp  =
c ******************************************************
c ******************************************************
      subroutine pkin70(q,iso,gout,nshels,outvv)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/iofile)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/shlg70)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/picon)
INCLUDE(common/ijlab)
INCLUDE(common/timez)
INCLUDE(common/pkfil)
INCLUDE(common/sortp)
INCLUDE(common/atmblk)
INCLUDE(common/parallel)
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
_IFN1(f)        call setsto(ngbf,0,q(iiptbf))
_IF1(f)        call vclr(q(iiptbf),1,ngbf)
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
_IF(parallel)
c***   **MPP**
      next = ipg_dlbtask()
c***   **MPP**
_ENDIF
c
c     ----- ishell -----
c
_IF(parallel)
      do 920 ii = nshels, ist0, -1
_ELSE
      do 920 ii = ist0,nshels
_ENDIF
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
_IF(newints)
      oisp = (kmax(ii)-kmin(ii)+1).eq.4
_ENDIF
c
c     ----- jshell -----
c
      j0 = jst0
_IF(parallel)
      do 900 jj = ii, j0, -1
_ELSE
      do 900 jj = j0,ii
_ENDIF
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
_IF(parallel)
c***   **MPP**
      icount_dlb = icount_dlb + 1
      if(icount_dlb . eq. next) then
c***   **MPP**
_ENDIF
_IF(newints)
      oijsp = oisp.or.((kmax(jj)-kmin(jj)+1).eq.4)
_ENDIF
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
_IF(newints)
      oijksp = oijsp.or.((kmax(kk)-kmin(kk)+1).eq.4)
_ENDIF
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
      q4 =  dfloat(nt)/ dfloat(n4)
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
_IF(newints)
      oijklsp = oijksp.or.((kmax(ll)-kmin(ll)+1).eq.4)
      if (oijklsp) then
         isptype = 0
      else
         isptype = 6
      endif
_ENDIF
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
_IF(parallel)
      next = ipg_dlbtask()
      endif
_ENDIF
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
_IF(parallel)
         call pg_dlbpush
_ENDIF
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
_IF(newints)
      isptype = -1
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/flip70)
      common/lt/lat,lbt,lct,ldt
      common/type/itype,jtype
INCLUDE(common/nshel)
INCLUDE(common/miscg)
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
_IF(newints)
      if (ngangb.eq.0) return
      if (isptype.ne.0.and.isptype.ne.6) then
         call caserr('LOGICAL ERROR: isptype not initialised')
      endif
      go to(440,480,500,520,540,560,
     &      580,600,620,640,660,680),jtype+isptype
_ELSE
      go to(440,480,500,520,540,560),jtype
_ENDIF
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
_IF(newints)
 580  call sp0000(gout)
      go to 1080
 600  call sp0002(gout)
      go to 1080
 620  call sp0022(gout)
      go to 1080
 640  call sp0202(gout)
      go to 1080
 660  call sp0222(gout)
      go to 1080
 680  call sp2222(gout)
_ENDIF
1080  return
      end
_IF(newints)
      subroutine genbuild70(fock,dmat,gout,fac1,fac2,facex,iexch,ocoul,
     &                      oexch,oskips)
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/flip70)
      common/lt/lat,lbt,lct,ldt
      common/type/itype,jtype
INCLUDE(common/nshel)
INCLUDE(common/miscg)
      dimension fock(*),dmat(*),gout(*)
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
      ib(1,2) = 1
      ib(2,2) = 2
      ib(3,2) = 3
      ib(4,2) = 4
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
      ib(1,2) = 2
      ib(2,2) = 1
      ib(3,2) = 3
      ib(4,2) = 4
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
      ib(1,2) = 1
      ib(2,2) = 2
      ib(3,2) = 4
      ib(4,2) = 3
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
      ib(1,2) = 3
      ib(2,2) = 4
      ib(3,2) = 1
      ib(4,2) = 2
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
      ib(1,2) = 2
      ib(2,2) = 1
      ib(3,2) = 4
      ib(4,2) = 3
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
      ib(1,2) = 3
      ib(2,2) = 4
      ib(3,2) = 2
      ib(4,2) = 1
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
      ib(1,2) = 4
      ib(2,2) = 3
      ib(3,2) = 1
      ib(4,2) = 2
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
      if (ngangb.eq.0) return
      if (isptype.ne.0.and.isptype.ne.6) then
         call caserr('LOGICAL ERROR: isptype not initialised')
      endif
      go to(440,480,500,520,540,560,
     &      580,600,620,640,660,680),jtype+isptype
 440  call sp0000(gout)
      call build0000(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 480  call sp0001(gout)
      call build0001(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 500  call sp0011(gout)
      call build0011(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 520  call sp0101(gout)
      call build0101(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 540  call sp0111(gout)
      call build0111(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 560  call sp1111(gout)
      call build1111(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 580  call sp0000(gout)
      call build0000(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 600  call sp0002(gout)
      call build0002(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 620  call sp0022(gout)
      call build0022(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 640  call sp0202(gout)
      call build0202(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 660  call sp0222(gout)
      call build0222(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      go to 1080
 680  call sp2222(gout)
      call build2222(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
1080  return
      end
_ENDIF
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
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/shlnos)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/miscg)
INCLUDE(common/geom)
INCLUDE(common/maxc)
INCLUDE(common/nshel)
INCLUDE(common/bshell)
INCLUDE(common/qgeom)
INCLUDE(common/const)
INCLUDE(common/pgeom)
INCLUDE(common/ginf)
      common/type/itype,jtype
      common /saveco/ uuuuuu,idadap(2),ijato,klato,ijshlo
INCLUDE(common/picon)
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
_IF(newints)
        if (rabsq.ne.0.0d0) then
_ENDIF
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
_IF(newints)
        ind = ind-1
        go to 36000
_ELSE
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
_ENDIF
18000   qq =  dexp(-qqq)*eab
        qqtest = cmaxa(i)*cmaxb(j)*qq
        if (qqtest-error1) 22000,22000,20000
20000   ismlp(ind) = 0
        go to 28000
22000   if (qqtest-error2) 26000,26000,24000
24000   ismlp(ind) = 1
        go to 28000
26000   continue
_IF(newints)
        ind = ind-1
        go to 36000
_ELSE
        ismlp(ind) = 2
_ENDIF
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
_IF(newints)
        ngangb=ind
c
        else
c
        ind=0
        do i = 1,nga
        ga = ag(i)
        csai = csa(i)
        cpai = cpa(i)
        do 36010 j = 1,ngb
        ind = ind+1
        gb = bg(j)
        gab = ga+gb
        gp(ind) = gab
        eab = done/gab
        ep(ind) = eab
        gbeab = gb*eab
        app(ind) = dzero
        bpp(ind) = dzero
        qq =  eab
        ismlp(ind) = 0
        qqqq = pito52*qq
        dp00p(ind) = qqqq*csai*csb(j)
        if (jtype-3) 36010,36010,30010
30010   dp01p(ind) = qqqq*csai*cpb(j)
        if (jtype-5) 32010,32010,34010
32010   conp(ind) = dp01p(ind)*eab
        dp00p(ind) = dp00p(ind)*gab/dp01p(ind)
        bpp(ind) = bpp(ind)*gab
        go to 36010
34010   dp10p(ind) = qqqq*cpai*csb(j)
        dp11p(ind) = qqqq*cpai*cpb(j)
        conp(ind) = dp11p(ind)
        dp00p(ind) = dp00p(ind)/dp11p(ind)
        dp01p(ind) = dp01p(ind)/dp11p(ind)
        dp10p(ind) = dp10p(ind)/dp11p(ind)
36010   continue
        enddo
        endif
_ENDIF
      endif
      return
      end
      subroutine sinset
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
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
