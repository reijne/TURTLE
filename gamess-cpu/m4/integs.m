c 
c  $Author: hvd $
c  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
c  $Locker:  $
c  $Revision: 6317 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integs.m,v $
c  $State: Exp $
c  
_IFN(vector)
      subroutine jkin7a(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c
INCLUDE(common/infoa)
INCLUDE(common/ijlab)
INCLUDE(common/timez)
INCLUDE(common/scfopt)
INCLUDE(common/pkfil)
INCLUDE(common/parallel)
INCLUDE(common/morokuma)
INCLUDE(common/nshel)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
c
_IF(ccpdft)
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
_ENDIF
      if(.not.odscf) then
      else
        l2 = nx
      endif
c
      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
_IF(parallel)
c***   **MPP**
      next = ipg_dlbtask()
c***   **MPP**
_ENDIF
c
c     ----- ishell -----
c
      call filmax
_IF(parallel)
      do 400 ii = nshels, ist0 , -1
_ELSE
      do 400 ii = ist0,nshels
_ENDIF
c
c     ----- print intermediate restart data -----
c
      if(kad(ii))400,401,401
 401  dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
       endif
c
c     ----- eliminate ishell -----
c
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
_IF(newints)
      oisp = (kmax(ii)-kmin(ii)+1).eq.4
_ENDIF
c
c     ----- jshell -----
c
      j0 = jst0
_IF(parallel)
      do 380 jj = ii, j0, -1
_ELSE
      do 380 jj = j0,ii
_ENDIF
      jst0 = 1
      if(kad(jj))380,141,141
c ***
c ***
141       itrij=ikyii+jj
          tolij=dlntol+prefac(itrij)
c
c **** rejecting on basis of ij arguably too extreme
c ***  suppress
c         if(tolij .le. -3.401e0) then
c           if(odscf)intcut(1)=intcut(1)+1
c           goto 380
c         endif
c ***
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
_IF(parallel)
c***   **MPP**
      icount_dlb = icount_dlb + 1
      if(icount_dlb . eq. next) then
c***   **MPP**
_ENDIF
      ishell = ii
      jshell = jj
      mij = itrij
_IF(newints)
      oijsp = oisp.or.((kmax(jj)-kmin(jj)+1).eq.4)
_ENDIF
c
c     ----- kshell -----
c
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
 361  do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ikykk=iky(kk)
      if (odscf) then
         itrik=ikyii+kk
         mik = itrik
         itrjk=iky(max(jj,kk))+min(jj,kk)
         mjk = itrjk
         tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
_IF(newints)
      oijksp = oijsp.or.((kmax(kk)-kmin(kk)+1).eq.4)
_ENDIF
c
c     ----- lshell ----
c
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
      if(kad(ll).lt.0) goto 340
      itrkl=ikykk+ll
      tijkl=tolij+prefac(itrkl)
      if(tijkl.le.0.0d0) then
        if(odscf)intcut(2)=intcut(2)+1
        goto 340
      endif
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=tijkl+
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            intcut(3)=intcut(3)+1
            goto 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dfloat(nt)/ dfloat(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
_IF(newints)
      oijklsp = oijksp.or.((kmax(ll)-kmin(ll)+1).eq.4)
      if (oijklsp) then
         isptype = 0
      else 
         isptype = 6
      endif
_ENDIF
c
c     ----- initialize gout to zero -----
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
      if (odscf) then
        if (zscftp.eq.'uhf')then
          call genr70(gout,1,.false.)
_IF(ccpdft)
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout,
     +         fac1,fac2, facex, ocoul, oexch)
_ELSE
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout)
_ENDIF
        else if (zscftp.eq.'gvb') then
          call genr70(gout,1,.false.)
          if(nsheld.le.1) then
            call dir_build_open_70(fock,exch,dens,gout)
          else
            call dir_build_open2_70(l2,fock,exch,dens,gout)
          endif
        else
_IF(ccpdft)
_IF(cray)
          call genr70(gout,1,.false.)
          call qoutd70(fock,dens,gout,
     +          fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
          if(omorok) then
            call genr70(gout,1,.false.)
            call dbuild70_morok(fock,dens,gout)
          else
_IF(newints)
            call genbuild70(fock,dens,gout,
     &           fac1, fac2, facex, 1, ocoul, oexch, .false.)
_ELSE
            call genr70(gout,1,.false.)
            call dbuild70(fock,dens,gout,
     &           fac1, fac2, facex, ocoul, oexch)
_ENDIF
          endif
_ENDIF
_ELSE
_IF(cray)
          call genr70(gout,1,.false.)
          call qoutd70(fock,dens,gout)
_ELSE
          call genr70(gout,1,.false.)
          call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
        endif
      else
        call genr70(gout,1,.false.)
        call qout70(gout)
      endif
c
c
c     ----- check cpu time/ maxblock condition -----
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim)go to 420
c
  340 continue
  360 continue
_IF(parallel)
      next = ipg_dlbtask()
      endif
_ENDIF
  380 continue
      time = cpulft(1)
  400 continue
      ist = 1
      jst = 1
      kst = 1
      lst = 1
_IF(parallel)
      call pg_dlbpush
_ENDIF
      call final(q,fock,dens)
_IF(newints)
      isptype = -1
_ENDIF
c
c     ----- reset core memory
c
  420 continue
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
      end
      subroutine jkin7s(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels,outvv)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c
INCLUDE(common/infoa)
INCLUDE(common/ijlab)
INCLUDE(common/timez)
INCLUDE(common/scfopt)
INCLUDE(common/pkfil)
INCLUDE(common/parallel)
INCLUDE(common/morokuma)
_IF(newscf)
INCLUDE(common/newscf)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      character *8 fnm
      character *6 snm
      data fnm/'integs.m'/
      data snm/'jkin7s'/
      data m25/25/
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
_IF(ccpdft)
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
_ENDIF

      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5

      ltri = ikyp(nshels)
      l2 = nx
      oskipp = .false.

      dlncutoff = dlog(cutoff)

      nschwz = 0
      ischw = igmem_alloc_inf(ltri,fnm,snm,'schwarz',IGMEM_DEBUG)
c
      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
c *** read in ints for schwarz inequality test
_IF(newscf)
      if (onews) then
         Call get_schwartz( q( ischw ) )
      else
_ENDIF
      call secget(isect(421),m25,iblk25)
      call rdedx(q(ischw),ltri,iblk25,idaf)
_IF(newscf)
      endif
_ENDIF
      dlnmxs = -1.0d50
      do ii = 0, ltri-1
         dlnmxs = dmax1(dlnmxs,q(ischw+ii))
      enddo
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
_IF(parallel)
c***   **MPP**
      next = ipg_dlbtask()
c***   **MPP**
_ENDIF
c
c     ----- ishell -----
c
      call filmax
_IF(parallel)
      do 400 ii = nshels, ist0 ,-1
_ELSE
      do 400 ii = ist0,nshels
_ENDIF
c
c     ----- print intermediate restart data -----
c
      if(kad(ii))400,401,401
 401  dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
       endif
c
c     ----- eliminate ishell -----
c
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
_IF(newints)
      oisp = (kmax(ii)-kmin(ii)+1).eq.4
_ENDIF
c
c     ----- jshell -----
c
      j0 = jst0
_IF(parallel)
      do 380 jj = ii,j0,-1
_ELSE
      do 380 jj = j0,ii
_ENDIF
      jst0 = 1
      if(kad(jj))380,141,141
c ***
141   itrij=ikyii+jj
      mij = itrij
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
_IF(parallel)
c***   **MPP**
      icount_dlb = icount_dlb + 1
      if(icount_dlb . eq. next) then
c***   **MPP**
_ENDIF
      ishell = ii
      jshell = jj
c *** reject on basis of ij 
      ijij = itrij + ischw -1
c     call scrkl(q(ischw),iso,tolij,nshels)
c     test = q(ijij) + tolij
c     if(test.lt.dlncutoff) then
c      if(odscf)intcut(1)=intcut(1)+1
      if (odscf) then
       tijkl = dlntol + q(ijij) + dlnmxs + dlnmxd
        if (tijkl .le. 0.0d0) then
         intcut(1)=intcut(1)+1
         nschwz = nschwz + itrij
_IF(parallel)
        goto 385
_ELSE
        goto 380
_ENDIF
        endif
      else
       if (q(ijij)+dlnmxs.lt.dlncutoff) then
        nschwz = nschwz + itrij
_IF(parallel)
       goto 385
_ELSE
       goto 380
_ENDIF
       endif
      endif
_IF(newints)
      oijsp = oisp.or.((kmax(jj)-kmin(jj)+1).eq.4)
_ENDIF
c
c     ----- kshell -----
c
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
 361  do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ikykk=iky(kk)
      if (odscf) then
         itrik=ikyii+kk
         mik = itrik
         itrjk=iky(max(jj,kk))+min(jj,kk)
         mjk = itrjk
         tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
_IF(newints)
      oijksp = oijsp.or.((kmax(kk)-kmin(kk)+1).eq.4)
_ENDIF
c
c     ----- lshell ----
c
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
c ***
      if(kad(ll).lt.0) goto 340
c ***
      itrkl=ikykk+ll
      klkl = itrkl + ischw -1
      test = q(ijij) + q(klkl)
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=dlntol + test +
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            nschwz = nschwz + 1
            intcut(3)=intcut(3)+1
            goto 340
         endif
      else
         oskipp = test.lt.dlncutoff
         if(oskipp) then
           nschwz = nschwz + 1
           go to 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dfloat(nt)/ dfloat(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
_IF(newints)
      oijklsp = oijksp.or.((kmax(ll)-kmin(ll)+1).eq.4)
      if (oijklsp) then
         isptype = 0
      else
         isptype = 6
      endif
_ENDIF
      if (odscf) then
        if (zscftp.eq.'uhf')then
          call genr70(gout,1,.false.)
_IF(ccpdft)
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout,
     +         fac1,fac2, facex, ocoul, oexch)
_ELSE
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout)
_ENDIF
        else if (zscftp.eq.'gvb'.or.zscftp.eq.'grhf') then
          call genr70(gout,1,.false.)
          if(nsheld.le.1) then
            call dir_build_open_70(fock,exch,dens,gout)
          else
            call dir_build_open2_70(l2,fock,exch,dens,gout)
          endif
        else
_IF(ccpdft)
_IF(cray)
          call genr70(gout,1,.false.)
          call qoutd70(fock,dens,gout,
     +          fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
          if(omorok) then
            call genr70(gout,1,.false.)
            call dbuild70_morok(fock,dens,gout)
          else
_IF(newints)
            call genbuild70(fock,dens,gout,
     &           fac1, fac2, facex, 1, ocoul, oexch, .false.)
_ELSE
            call genr70(gout,1,.false.)
            call dbuild70(fock,dens,gout,
     &           fac1, fac2, facex, ocoul, oexch)
_ENDIF
          endif
_ENDIF
_ELSE
_IF(cray)
          call genr70(gout,1,.false.)
          call qoutd70(fock,dens,gout)
_ELSE
          call genr70(gout,1,.false.)
          call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
        endif
      else
        call genr70(gout,1,.false.)
        call qout70(gout)
      endif
c
c
c     ----- check cpu time/ maxblock condition -----
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim)go to 420
c
  340 continue
  360 continue
_IF(parallel)
  385 next = ipg_dlbtask()
      endif
_ENDIF
  380 continue
      time = cpulft(1)
  400 continue
      ist = 1
      jst = 1
      kst = 1
      lst = 1
_IF(parallel)
      call pg_dlbpush
_ENDIF
      call final(q,fock,dens)
_IF(newints)
      isptype = -1
_ENDIF
c
c     ----- reset core memory
c
*IJB Sum nschwz
c$$$  420 if(outvv) write(iwr,6030) nschwz
 420  call pg_igop( 160467, nschwz, 1, '+' )
      if(outvv) write(iwr,6030) nschwz

      call gmem_free_inf(ischw,fnm,snm,'schwarz')
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
 6030 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine scrkl(schw,iso,scrnij,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c 
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/ijlab)
      dimension schw(*),iso(nshels,*)
c
      scrnij = -1.0d50
      do 360 kk = 1,ishell
      if(kad(kk).ge.0) then
       do 220 it = 1,nt
       kd = iso(kk,it)
       if (kd .gt. ishell) go to 360
  220  maxll = kk
       if (kk .eq. ishell) maxll = jshell
       ikykk=iky(kk)
       do 340 ll = 1,maxll
       if(kad(ll).ge.0) then
        itrkl=ikykk+ll
        scrnij = max(scrnij,schw(itrkl))
       endif
  340  continue
      endif
  360 continue
      return
      end
_IF(newints)
      subroutine build0000(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for s-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
      dimension fock(*),dmat(*)
c
      i1 = kloc(inew)
      i2 = kloc(jnew)
      i3 = kloc(knew)
      i4 = kloc(lnew)
      facij = fac1
      fackl = fac2
      if (dabs(gout) .lt. cutoff) return
c
      if (ocoul) then
         itr12=iky(i1)+i2
         itr34=iky(i3)+i4
         f12 = fock(itr12) + 4.0d0*facij*gout*dmat(itr34)
         fock(itr34) = fock(itr34) + 4.0d0*fackl*gout*dmat(itr12)
         fock(itr12) = f12
      endif
c
      if (oexch) then
         itr13=iky(i1)+i3
         itr14=iky(i1)+i4
         itr23=iky(max(i2,i3))+min(i2,i3)
         itr24=iky(max(i2,i4))+min(i2,i4)
         val=gout*facex
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
         if(i2.eq.i3) val14=val+val
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
      end
      subroutine build0002(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 3 s-shell 
c     and 1 p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4)
c
      locs(1) = kloc(ishell)
      locs(2) = kloc(jshell)
      locs(3) = kloc(kshell)
      locs(4) = kloc(lshell)
      facij = fac1
      fackl = fac2
      locs(ib(4,2)) = locs(ib(4,2))-1

      do 100 l = 2,4
         val = gout(l)
         locs(ib(4,2)) = locs(ib(4,2))+1
         if (dabs(val) .lt. cutoff) goto 100
         i1 = locs(1)
         i2 = locs(2)
         i3 = locs(3)
         i4 = locs(4)
c
         if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4.0d0*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4.0d0*fackl*val*dmat(itr12)
            fock(itr12) = f12
         endif
c
         if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
         endif
c
 100  continue
      end
      subroutine build0001(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 3 s-shell 
c     and 1 l/p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4)
c
      minl    = kmin(lnew)
      locs(1) = kloc(ishell)
      locs(2) = kloc(jshell)
      locs(3) = kloc(kshell)
      locs(4) = kloc(lshell)
      facij = fac1
      fackl = fac2
      locs(ib(4,2)) = locs(ib(4,2))-1

      do 100 l = minl,4
         val = gout(l)
         locs(ib(4,2)) = locs(ib(4,2))+1
         if (dabs(val) .lt. cutoff) goto 100
         i1 = locs(1)
         i2 = locs(2)
         i3 = locs(3)
         i4 = locs(4)
c
         if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4.0d0*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4.0d0*fackl*val*dmat(itr12)
            fock(itr12) = f12
         endif
c
         if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
         endif
c
 100  continue
      end
      subroutine build0011(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 2 s-shell 
c     and 2 l/p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      mink    = kmin(knew)
      minl    = kmin(lnew)
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(3,2)) = loct(ib(3,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
c     oident = ishell .eq. kshell .and. jshell .eq. lshell
      okandl = knew.eq.lnew
c
      maxl = 4
      n1 = 1
      if (mink.eq.1) n1 = -3
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      locs(ib(3,2)) = loct(ib(3,2))
      do 110 k = mink,4
        n1 = n1 + 4
        n2 = n1
        if (minl.eq.1) n2 = n2-1
        locs(ib(3,2)) = locs(ib(3,2))+1
        locs(ib(4,2)) = loct(ib(4,2))
        if (okandl) maxl = k
        do 100 l = minl,maxl
          n2  = n2 + 1
          locs(ib(4,2)) = locs(ib(4,2))+1
          val = gout(n2)
          if (dabs(val) .lt. cutoff) goto 100
          if(locs(1).ge.locs(3)) then
            i1 = locs(1)
            i2 = locs(2)
            i3 = locs(3)
            i4 = locs(4)
            facij=fac1
            fackl=fac2
          else
            i1 = locs(3)
            i2 = locs(4)
            i3 = locs(1)
            i4 = locs(2)
            facij=fac2
            fackl=fac1
          endif
c
          if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
            fock(itr12) = f12
          endif
c
          if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
          endif
c
 100    continue
 110  continue
      end
      subroutine build0022(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 2 s-shell 
c     and 2 p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(3,2)) = loct(ib(3,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
c     oident = ishell .eq. kshell .and. jshell .eq. lshell
      okandl = knew.eq.lnew
c
      maxl = 4
      n1 = 1
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      locs(ib(3,2)) = loct(ib(3,2))
      do 110 k = 2,4
        n1 = n1 + 4
        n2 = n1
        locs(ib(3,2)) = locs(ib(3,2))+1
        locs(ib(4,2)) = loct(ib(4,2))
        if (okandl) maxl = k
        do 100 l = 2,maxl
          n2  = n2 + 1
          locs(ib(4,2)) = locs(ib(4,2))+1
          val = gout(n2)
          if (dabs(val) .lt. cutoff) goto 100
          if(locs(1).ge.locs(3)) then
            i1 = locs(1)
            i2 = locs(2)
            i3 = locs(3)
            i4 = locs(4)
            facij=fac1
            fackl=fac2
          else
            i1 = locs(3)
            i2 = locs(4)
            i3 = locs(1)
            i4 = locs(2)
            facij=fac2
            fackl=fac1
          endif
c
          if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
            fock(itr12) = f12
          endif
c
          if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
          endif
c
 100    continue
 110  continue
      end
      subroutine build0101(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 2 s-shell 
c     and 2 l/p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      minj    = kmin(jnew)
      minl    = kmin(lnew)
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(2,2)) = loct(ib(2,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
      oident = ishell .eq. kshell .and. jshell .eq. lshell
c
      maxl = 4
      n1 = 1
      if (minj.eq.1) n1 = -15
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      locs(ib(3,2)) = loct(ib(3,2))
      do 110 j = minj,4
        n1 = n1 + 16
        n2 = n1
        if (minl.eq.1) n2 = n2-1
        locs(ib(2,2)) = locs(ib(2,2))+1
        locs(ib(4,2)) = loct(ib(4,2))
        if (oident) maxl = j
        do 100 l = minl,maxl
          n2  = n2 + 1
          locs(ib(4,2)) = locs(ib(4,2))+1
          val = gout(n2)
          if (dabs(val) .lt. cutoff) goto 100
          if(locs(1).ge.locs(3)) then
            i1 = locs(1)
            i2 = locs(2)
            i3 = locs(3)
            i4 = locs(4)
            facij=fac1
            fackl=fac2
          else
            i1 = locs(3)
            i2 = locs(4)
            i3 = locs(1)
            i4 = locs(2)
            facij=fac2
            fackl=fac1
          endif
c
          if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
            fock(itr12) = f12
          endif
c
          if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
          endif
c
 100    continue
 110  continue
      end
      subroutine build0202(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 2 s-shell 
c     and 2 p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(2,2)) = loct(ib(2,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
      oident = ishell .eq. kshell .and. jshell .eq. lshell
c
      maxl = 4
      n1 = 1
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      locs(ib(3,2)) = loct(ib(3,2))
      do 110 j = 2,4
        n1 = n1 + 16
        n2 = n1
        locs(ib(2,2)) = locs(ib(2,2))+1
        locs(ib(4,2)) = loct(ib(4,2))
        if (oident) maxl = j
        do 100 l = 2,maxl
          n2  = n2 + 1
          locs(ib(4,2)) = locs(ib(4,2))+1
          val = gout(n2)
          if (dabs(val) .lt. cutoff) goto 100
          if(locs(1).ge.locs(3)) then
            i1 = locs(1)
            i2 = locs(2)
            i3 = locs(3)
            i4 = locs(4)
            facij=fac1
            fackl=fac2
          else
            i1 = locs(3)
            i2 = locs(4)
            i3 = locs(1)
            i4 = locs(2)
            facij=fac2
            fackl=fac1
          endif
c
          if (ocoul) then
            itr12=iky(i1)+i2
            itr34=iky(i3)+i4
            f12 = fock(itr12) + 4*facij*val*dmat(itr34)
            fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
            fock(itr12) = f12
          endif
c
          if (oexch) then
            itr13=iky(i1)+i3
            itr14=iky(i1)+i4
            itr23=iky(max(i2,i3))+min(i2,i3)
            itr24=iky(max(i2,i4))+min(i2,i4)
            val=val*facex
            val13=val
            val14=val
            if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
            if(i2.eq.i3) val14=val+val
            f23 = fock(itr23) - val14*dmat(itr14)
            f14 = fock(itr14) - val14*dmat(itr23)
            f13 = fock(itr13) - val13*dmat(itr24)
            fock(itr24) = fock(itr24) - val13*dmat(itr13)
            fock(itr23) = f23
            fock(itr14) = f14
            fock(itr13) = f13
          endif
c
 100    continue
 110  continue
      end
      subroutine build0111(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 1 s-shell 
c     and 3 l/p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      minj    = kmin(jnew)
      mink    = kmin(knew)
      minl    = kmin(lnew)
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(2,2)) = loct(ib(2,2))-1
      loct(ib(3,2)) = loct(ib(3,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
      okandl = knew.eq.lnew
c
      maxl = 4
      n1 = 1
      if (minj.eq.1) n1 = -15
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      do 110 j = minj,4
        n1 = n1 + 16
        n2 = n1
        if (mink.eq.1) n2 = n2-4
        locs(ib(2,2)) = locs(ib(2,2))+1
        locs(ib(3,2)) = loct(ib(3,2))
        do k = mink,4
          n2 = n2 + 4
          n3 = n2
          if (minl.eq.1) n3 = n3-1
          locs(ib(3,2)) = locs(ib(3,2))+1
          locs(ib(4,2)) = loct(ib(4,2))
          if (okandl) maxl = k
          do 100 l = minl,maxl
            n3  = n3 + 1
            locs(ib(4,2)) = locs(ib(4,2))+1
            val = gout(n3)
            if (dabs(val) .lt. cutoff) goto 100
            if(locs(1).ge.locs(3)) then
              i1 = locs(1)
              i2 = locs(2)
              i3 = locs(3)
              i4 = locs(4)
              facij=fac1
              fackl=fac2
            else
              i1 = locs(3)
              i2 = locs(4)
              i3 = locs(1)
              i4 = locs(2)
              facij=fac2
              fackl=fac1
            endif
c
            if (ocoul) then
              itr12=iky(i1)+i2
              itr34=iky(i3)+i4
              f12 = fock(itr12) + 4*facij*val*dmat(itr34)
              fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
              fock(itr12) = f12
            endif
c
            if (oexch) then
              itr13=iky(i1)+i3
              itr14=iky(i1)+i4
              itr23=iky(max(i2,i3))+min(i2,i3)
              itr24=iky(max(i2,i4))+min(i2,i4)
              val=val*facex
              val13=val
              val14=val
              if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
              if(i2.eq.i3) val14=val+val
              f23 = fock(itr23) - val14*dmat(itr14)
              f14 = fock(itr14) - val14*dmat(itr23)
              f13 = fock(itr13) - val13*dmat(itr24)
              fock(itr24) = fock(itr24) - val13*dmat(itr13)
              fock(itr23) = f23
              fock(itr14) = f14
              fock(itr13) = f13
            endif
c
 100      continue
        enddo
 110  continue
      end
      subroutine build0222(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for 1 s-shell 
c     and 3 p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
INCLUDE(common/flip70)
      dimension fock(*),dmat(*),gout(4)
      dimension locs(4),loct(4)
c
      loct(1) = kloc(ishell)
      loct(2) = kloc(jshell)
      loct(3) = kloc(kshell)
      loct(4) = kloc(lshell)
      loct(ib(2,2)) = loct(ib(2,2))-1
      loct(ib(3,2)) = loct(ib(3,2))-1
      loct(ib(4,2)) = loct(ib(4,2))-1
c
      okandl = knew.eq.lnew
c
      maxl = 4
      n1 = 1
      locs(ib(1,2)) = loct(ib(1,2))
      locs(ib(2,2)) = loct(ib(2,2))
      do 110 j = 2,4
        n1 = n1 + 16
        n2 = n1
        locs(ib(2,2)) = locs(ib(2,2))+1
        locs(ib(3,2)) = loct(ib(3,2))
        do k = 2,4
          n2 = n2 + 4
          n3 = n2
          locs(ib(3,2)) = locs(ib(3,2))+1
          locs(ib(4,2)) = loct(ib(4,2))
          if (okandl) maxl = k
          do 100 l = 2,maxl
            n3  = n3 + 1
            locs(ib(4,2)) = locs(ib(4,2))+1
            val = gout(n3)
            if (dabs(val) .lt. cutoff) goto 100
            if(locs(1).ge.locs(3)) then
              i1 = locs(1)
              i2 = locs(2)
              i3 = locs(3)
              i4 = locs(4)
              facij=fac1
              fackl=fac2
            else
              i1 = locs(3)
              i2 = locs(4)
              i3 = locs(1)
              i4 = locs(2)
              facij=fac2
              fackl=fac1
            endif
c
            if (ocoul) then
              itr12=iky(i1)+i2
              itr34=iky(i3)+i4
              f12 = fock(itr12) + 4*facij*val*dmat(itr34)
              fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
              fock(itr12) = f12
            endif
c
            if (oexch) then
              itr13=iky(i1)+i3
              itr14=iky(i1)+i4
              itr23=iky(max(i2,i3))+min(i2,i3)
              itr24=iky(max(i2,i4))+min(i2,i4)
              val=val*facex
              val13=val
              val14=val
              if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
              if(i2.eq.i3) val14=val+val
              f23 = fock(itr23) - val14*dmat(itr14)
              f14 = fock(itr14) - val14*dmat(itr23)
              f13 = fock(itr13) - val13*dmat(itr24)
              fock(itr24) = fock(itr24) - val13*dmat(itr13)
              fock(itr23) = f23
              fock(itr14) = f14
              fock(itr13) = f13
            endif
c
 100      continue
        enddo
 110  continue
      end
      subroutine build1111(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for l/p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
      dimension fock(*),dmat(*),gout(256)
c
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      iii1 = kloc(ishell)-1
      iii2 = kloc(jshell)-1
      iii3 = kloc(kshell)-1
      iii4 = kloc(lshell)-1

      oiandj = inew.eq.jnew
      oiandk = inew.eq.knew
      okandl = knew.eq.lnew
      oident = ishell .eq. kshell .and. jshell .eq. lshell

      maxj = 4
c     maxk = 4
      maxl = 4

      if (oident) then
        facij = fac1
        fackl = fac2
        i1  = iii1
        n1  = 1
        if (mini.eq.1) n1 = -63
        ijn = 0
        do i = mini, 4
          i1 = i1 + 1
          n1 = n1 + 64
          i2 = iii2
          n2 = n1
          if (minj.eq.1) n2 = n2-16
          if (oiandj) maxj = i
          do 110 j = minj, maxj
            ijn = ijn + 1
            i2  = i2  + 1
            n2  = n2  + 16
            i3  = iii3
            n3  = n2
            if (mink.eq.1) n3 = n3-4
            kln = 0
            do k = mink, 4
              i3 = i3 + 1
              n3 = n3 + 4
              i4 = iii4
              n4 = n3
              if (minl.eq.1) n4 = n4-1
              if (okandl) maxl = k
              do 100 l = minl, maxl
                kln = kln + 1
                if (kln.gt.ijn) goto 110
                i4 = i4 + 1
                n4 = n4 + 1
                val = gout(n4)
                if (dabs(val) .lt. cutoff) goto 100
c
                if (ocoul) then
                   itr12=iky(i1)+i2
                   itr34=iky(i3)+i4
                   f12 = fock(itr12) + 4*facij*val*dmat(itr34)
                   fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
                   fock(itr12) = f12
                endif
c
                if (oexch) then
                   itr13=iky(i1)+i3
                   itr14=iky(i1)+i4
                   itr23=iky(max(i2,i3))+min(i2,i3)
                   itr24=iky(max(i2,i4))+min(i2,i4)
                   val=val*facex
                   val13=val
                   val14=val
                   if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
                   if(i2.eq.i3) val14=val+val
                   f23 = fock(itr23) - val14*dmat(itr14)
                   f14 = fock(itr14) - val14*dmat(itr23)
                   f13 = fock(itr13) - val13*dmat(itr24)
                   fock(itr24) = fock(itr24) - val13*dmat(itr13)
                   fock(itr23) = f23
                   fock(itr14) = f14
                   fock(itr13) = f13
                endif
c
 100          continue
            enddo
 110      continue
        enddo
      else
        ii1 = iii1
        n1  = 1
        if (mini.eq.1) n1 = -63
        do i = mini, 4
          ii1 = ii1 + 1
          n1  = n1  + 64
          if (oiandj) maxj=i
          ii2 = iii2
          n2  = n1
          if (minj.eq.1) n2 = n2-16
          do 210 j = minj, maxj
            ii2 = ii2 + 1
            n2  = n2  + 16
            ii3 = iii3
            n3  = n2
            if (mink.eq.1) n3 = n3-4
            do k = mink, 4
              ii3 = ii3 + 1
              n3  = n3  + 4
              if (okandl) maxl=k
              ii4 = iii4
              n4  = n3
              if (minl.eq.1) n4 = n4-1
              do 200 l = minl, maxl
                ii4 = ii4 + 1
                n4  = n4  + 1
                val = gout(n4)
                if (dabs(val) .lt. cutoff) goto 200
c
                if (ii1.ge.ii3) then
                   i1 = ii1
                   i2 = ii2
                   i3 = ii3
                   i4 = ii4
                   facij = fac1
                   fackl = fac2
                else
                   i1 = ii3
                   i2 = ii4
                   i3 = ii1
                   i4 = ii2
                   facij = fac2
                   fackl = fac1
                endif
c
                if (ocoul) then
                  itr12=iky(i1)+i2
                  itr34=iky(i3)+i4
                  f12 = fock(itr12) + 4*facij*val*dmat(itr34)
                  fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
                  fock(itr12) = f12
                endif
c
                if (oexch) then
                  itr13=iky(i1)+i3
                  itr14=iky(i1)+i4
                  itr23=iky(max(i2,i3))+min(i2,i3)
                  itr24=iky(max(i2,i4))+min(i2,i4)
                  val=val*facex
                  val13=val
                  val14=val
                  if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
                  if(i2.eq.i3) val14=val+val
                  f23 = fock(itr23) - val14*dmat(itr14)
                  f14 = fock(itr14) - val14*dmat(itr23)
                  f13 = fock(itr13) - val13*dmat(itr24)
                  fock(itr24) = fock(itr24) - val13*dmat(itr13)
                  fock(itr23) = f23
                  fock(itr14) = f14
                  fock(itr13) = f13
                endif
c
 200          continue
            enddo
 210      continue
        enddo
      endif
      end
      subroutine build2222(fock,dmat,gout,fac1,fac2,facex,ocoul,oexch)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     This is subroutine dbuild70 specialised for p-shell quartets
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/nshel)
INCLUDE(common/shlt)
      dimension fock(*),dmat(*),gout(256)
c
      iii1 = kloc(ishell)-1
      iii2 = kloc(jshell)-1
      iii3 = kloc(kshell)-1
      iii4 = kloc(lshell)-1

      oiandj = inew.eq.jnew
      oiandk = inew.eq.knew
      okandl = knew.eq.lnew
      oident = ishell .eq. kshell .and. jshell .eq. lshell

      maxj = 4
c     maxk = 4
      maxl = 4

      if (oident) then
        facij = fac1
        fackl = fac2
        i1  = iii1
        n1  = 1
        ijn = 0
        do i = 2, 4
          i1 = i1 + 1
          n1 = n1 + 64
          i2 = iii2
          n2 = n1
          if (oiandj) maxj = i
          do 110 j = 2, maxj
            ijn = ijn + 1
            i2  = i2  + 1
            n2  = n2  + 16
            i3  = iii3
            n3  = n2
            kln = 0
            do k = 2, 4
              i3 = i3 + 1
              n3 = n3 + 4
              i4 = iii4
              n4 = n3
              if (okandl) maxl = k
              do 100 l = 2, maxl
                kln = kln + 1
                if (kln.gt.ijn) goto 110
                i4 = i4 + 1
                n4 = n4 + 1
                val = gout(n4)
                if (dabs(val) .lt. cutoff) goto 100
c
                if (ocoul) then
                   itr12=iky(i1)+i2
                   itr34=iky(i3)+i4
                   f12 = fock(itr12) + 4*facij*val*dmat(itr34)
                   fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
                   fock(itr12) = f12
                endif
c
                if (oexch) then
                   itr13=iky(i1)+i3
                   itr14=iky(i1)+i4
                   itr23=iky(max(i2,i3))+min(i2,i3)
                   itr24=iky(max(i2,i4))+min(i2,i4)
                   val=val*facex
                   val13=val
                   val14=val
                   if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
                   if(i2.eq.i3) val14=val+val
                   f23 = fock(itr23) - val14*dmat(itr14)
                   f14 = fock(itr14) - val14*dmat(itr23)
                   f13 = fock(itr13) - val13*dmat(itr24)
                   fock(itr24) = fock(itr24) - val13*dmat(itr13)
                   fock(itr23) = f23
                   fock(itr14) = f14
                   fock(itr13) = f13
                endif
c
 100          continue
            enddo
 110      continue
        enddo
      else
        ii1 = iii1
        n1  = 1
        do i = 2, 4
          ii1 = ii1 + 1
          n1  = n1  + 64
          if (oiandj) maxj=i
          ii2 = iii2
          n2  = n1
          do 210 j = 2, maxj
            ii2 = ii2 + 1
            n2  = n2  + 16
            ii3 = iii3
            n3  = n2
            do k = 2, 4
              ii3 = ii3 + 1
              n3  = n3  + 4
              if (okandl) maxl=k
              ii4 = iii4
              n4  = n3
              do 200 l = 2, maxl
                ii4 = ii4 + 1
                n4  = n4  + 1
                val = gout(n4)
                if (dabs(val) .lt. cutoff) goto 200
c
                if (ii1.ge.ii3) then
                   i1 = ii1
                   i2 = ii2
                   i3 = ii3
                   i4 = ii4
                   facij = fac1
                   fackl = fac2
                else
                   i1 = ii3
                   i2 = ii4
                   i3 = ii1
                   i4 = ii2
                   facij = fac2
                   fackl = fac1
                endif
c
                if (ocoul) then
                  itr12=iky(i1)+i2
                  itr34=iky(i3)+i4
                  f12 = fock(itr12) + 4*facij*val*dmat(itr34)
                  fock(itr34) = fock(itr34) + 4*fackl*val*dmat(itr12)
                  fock(itr12) = f12
                endif
c
                if (oexch) then
                  itr13=iky(i1)+i3
                  itr14=iky(i1)+i4
                  itr23=iky(max(i2,i3))+min(i2,i3)
                  itr24=iky(max(i2,i4))+min(i2,i4)
                  val=val*facex
                  val13=val
                  val14=val
                  if(i1.eq.i3 .or. i2.eq.i4) val13=val+val
                  if(i2.eq.i3) val14=val+val
                  f23 = fock(itr23) - val14*dmat(itr14)
                  f14 = fock(itr14) - val14*dmat(itr23)
                  f13 = fock(itr13) - val13*dmat(itr24)
                  fock(itr24) = fock(itr24) - val13*dmat(itr13)
                  fock(itr23) = f23
                  fock(itr14) = f14
                  fock(itr13) = f13
                endif
c
 200          continue
            enddo
 210      continue
        enddo
      endif
      end
_ENDIF
      subroutine sp0000(gout)
c        *****  special fast routine for -p- loop for 0000 ****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc(5001)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/ginf)
INCLUDE(common/pqgeom)
INCLUDE(common/pgeom)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
c
      dimension gout(*)
c
      data sixty/60.0d0/
c
      g0000 = 0.d0
      do 220 k = 1,ngc
      gc = cgg(k)
      csck = csc(k)
      gcrcds = gc*rcdsq
      do 220 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.d0/gcd
      gdecd = gd*ecd
      ppp = -gdecd*gcrcds
      if(ppp+sixty)220,500,500
  500 pp = ecd* dexp(ppp)
  520 qqtest = pp*cmaxc(k)*cmaxd(l)
      if (qqtest .le. error1) go to 100
      ismlq = 0
      go to 120
  100 if (qqtest .le. error2) go to 220
      ismlq = 1
  120 cq = gdecd*rcd
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      h0000 = 0.d0
      do 200 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 200
      auxvar = var(isml+1)
      p = ((aqz-app(i))**2+qperp2)/(ep(i)+ecd)
      if (p .le. auxvar) go to 180
      h0000 = h0000+dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      go to 200
  180 continue
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2)
      h0000 = h0000+ dp00p(i)*f0/dsqrt(gp(i)+gcd)
  200 continue
      g0000 = g0000+h0000*csck*csd(l)*pp
  220 continue
      gout(1) = g0000
      return
      end
      subroutine filmax
c
c     --------------------------
c     --------------------------
c
c     finds maximum value of s and p coefficients
c     also sets limits determining how accurately a set of integrals
c     need be evaluated in order to guarantee an overall integral
c     accuracy of 10**-6
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/auxvar)
INCLUDE(common/nshel)
INCLUDE(common/maxc)
INCLUDE(common/shlt)
c
      data five/5.0d0/,fiften/15.0d0/
      data ten,twenty/10.0d0,20.0d0/
c
      do 100 i = 1,nshell
      l = kstart(i)
      n = l+kng(i)-1
      do 100 j = l,n
      a1 =  dabs(cs(j))
      a2 =  dabs(cp(j))
      cmax(j) = dmax1(a1,a2)
  100 continue
c
      if (ifasti.eq.1) then
c     settings in gamess-us
       var(1) = twenty
       var(2) = ten
       error1 = dmin1(cutoff,1.0d-10)
       error2 = dmin1(cutoff,1.0d-10)
      else if (ifasti.eq.2) then
c     tighten thresholds
       var(1) = fiften
       var(2) = five
       error1 = dsqrt(cutoff) * 0.01d0
       error2 = cutoff * 0.001d0
      else
_IF1()c      error1 = sqrt(cutoff)
_IF1()c      error2 = cutoff
       error1 = dsqrt(cutoff) * 0.1d0
       var(1) = fiften
       var(2) = five
       error2 = cutoff * 0.01d0
      endif
      return
      end
      subroutine qout70(gout)
c
c ... qout70 now handles conventional integrals only
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IFN1(iv)      common/craypk/integ(1)
_IF1(iv)      common/craypk/integ(340),intkl(340)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/nshel)
INCLUDE(common/misc)
      common/blkin/goutx(510),nword
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/shlg70)
      dimension gout(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c     data pt5 /0.5/
c
c     ----- pack the 4 indices of integral into one word
c     ----- write label + integral on mainfile
c
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      ijn = 0
      jmax = maxj
        do 2600 i = 1,maxi
          if (oianj) jmax = i
          i1 = loci + i
          ipack = i4096(i1)
          do 2400 j = 1,jmax
            ijn = ijn+1
            n1 = ib(ib1,i)+ib(jb1,j)+1
            lmax = maxl
            i2 = locj + j
            if(i1-i2)10, 20, 20
  10        lab1=i4096(i2) + i1
            go to 30
  20        lab1=ipack + i2
  30        kln = 0
            do 2200 k = 1,maxk
              if (okanl) lmax = k
              i3 = lock+k
              kpack = i4096(i3)
              do 2000 l = 1,lmax
                kln = kln+1
                if (oident .and. kln .gt. ijn) go to 2400
                nn = n1+ib(kb1,k)+ib(lb1,l)
                val = gout(nn)
                if ( dabs(val) .lt. cutoff) go to 2000
                i4 = locl+l
                if(i3-i4) 40, 50, 50
 40             lab2=i4096(i4) + i3
                go to 60
 50             lab2=kpack + i4
 60             goutx(icount) = val
_IFN1(iv)                integ(ic4  )= max(lab1,lab2)
_IFN1(iv)                integ(ic4+1)= min(lab1,lab2)
_IFN1(iv)                ic4 = ic4 + 2
_IF1(iv)                integ(icount)= max(lab1,lab2)
_IF1(iv)                intkl(icount)= min(lab1,lab2)
                icount = icount+1
                if (icount .le. nintmx) go to 2000
                call blocki
                if(omaxb)go to 261
2000          continue
2200        continue
2400      continue
2600    continue
  261 return
      end
_IF(unicos)
_IF(ccpdft)
      subroutine qoutd70(fock,dmat,g,
     +    fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine qoutd70(fock,dmat,g)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/cslosc)
_IFN1(iv)      common/craypk/integ(1)
_IF1(iv)      common/craypk/integ(340),intkl(340)
INCLUDE(common/ijlab)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/nshel)
INCLUDE(common/misc)
      common/blkin/goutx(510),nword
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/shlg70)
      dimension fock(*),dmat(*),g(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
c     data pt5 /0.5/
c
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      ijn = 0
      jmax = maxj
_IF(ccpdft)
      if (.not. odft) then
_ENDIF
c     ----- pack the 4 indices of integral into one word
c     ----- write label + integral on mainfile to use CAL dbuild
c
      ifac = 4096
        do 260 i = 1,maxi
        i1=loci+i
        ipack = ifac*i1
          if (oianj) jmax = i
          do 240 j = 1,jmax
            ijn = ijn+1
            n1 = ib(ib1,i)+ib(jb1,j)+1
            lmax = maxl
            i2=locj+j
            if(i1-i2) 1,2,2
    1       lab1 = ifac*i2 + i1
            go to 3
   2        lab1=ipack+i2
   3        continue
            kln=0
            do 220 k = 1,maxk
              if (okanl) lmax = k
              i3=lock+k
              kpack = ifac*i3
              do 200 l = 1,lmax
                kln = kln+1
                if (oident .and. kln .gt. ijn) go to 240
                nn = n1+ib(kb1,k)+ib(lb1,l)
                val = g(nn)
                if ( dabs(val) .lt. cutoff) go to 200
                goutx(icount) = val
                i4=locl+l
                if(i3-i4)4,5,5
    4           lab2 = ifac*i4 + i3
                go to 6
    5           lab2=kpack+i4
    5           integ(ic4  )=max(lab1,lab2)
                integ(ic4+1)=min(lab1,lab2)
                ic4=ic4+2
                icount=icount+1
                if(icount.le.nintmx) go to 200
                ochek=.false.
                nrec=nrec+1
                call dbuild(fock,dmat)
  200         continue
  220       continue
  240     continue
  260   continue
_IF(ccpdft)
      else
c
c     non conventional fock build involved .. cannot use CAL
c
      do 360 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 340 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 320 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 300 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 340
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 300
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
c
c Coulomb terms
c
      if(ocoul)then
         itr12=iky(i1)+i2
         itr34=iky(i3)+i4
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c    
c Full exchange term for HF or weighted exchange term for b3lyp etc
c    
         itr13=iky(i1)+i3
         itr14=iky(i1)+i4
         itr23=iky(max(i2,i3))+min(i2,i3)
         itr24=iky(max(i2,i4))+min(i2,i4)
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
  300 continue
  320 continue
  340 continue
  360 continue
c
      endif
_ENDIF
      return
      end
_ELSE
_IF(ccpdft)
      subroutine dbuild70(fock,dmat,gout,
     &     fac1,fac2, facex, ocoul, oexch)
_ELSE
      subroutine dbuild70(fock,dmat,gout)
_ENDIF

      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/shlnos)
INCLUDE(common/nshel)
INCLUDE(common/misc)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/shlg70)

      dimension ib(4,4)
      dimension fock(*),dmat(*),gout(*)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
_IF(newints)
      ii1 = int1
      ii2 = int2
_ELSE
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
_ENDIF
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
_IF(newints)
        ii3 = int3
        ii4 = int4
_ELSE
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
_ENDIF
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
_IF(ccpdft)
           facij=fac1
           fackl=fac2
_ENDIF
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
_IF(ccpdft)
           facij=fac2
           fackl=fac1
_ENDIF
        endif
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
_IF(ccpdft)
c
c Coulomb terms
c
      if(ocoul)then
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c     
c Full exchange term for HF or weighted exchange term for b3lyp etc
c     
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
_ELSE
        val2=val+val
        val4=val2+val2
        val13=val
        val14=val
        if(i1.eq.i3 .or. i2.eq.i4) val13=val2
        if(i2.eq.i3) val14=val2
        f12 = val4*dmat(itr34) + fock(itr12)
        fock(itr34) = val4*dmat(itr12) + fock(itr34)
        fock(itr12) = f12
        f23 = fock(itr23) - val14*dmat(itr14)
        f14 = fock(itr14) - val14*dmat(itr23)
        f13 = fock(itr13) - val13*dmat(itr24)
        fock(itr24) = fock(itr24) - val13*dmat(itr13)
        fock(itr23) = f23
        fock(itr14) = f14
        fock(itr13) = f13
_ENDIF
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_ENDIF

_IF(ccpdft)
      subroutine dir_build_uhf70(fock,ak,p,q,gout,
     &   fac1, fac2, facex, ocoul, oexch)
_ELSE
      subroutine dir_build_uhf70(fock,ak,p,q,gout)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/shlg70)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),fock(*),ak(*),p(*),q(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell UHF Fock builder 
c
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
_IF(ccpdft)
           facij=fac1
           fackl=fac2
_ENDIF
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
_IF(ccpdft)
           facij=fac2
           fackl=fac1
_ENDIF
        endif
      gik=val
      val2=gik+gik
      val4=val2+val2
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
_IF(ccpdft)
c
c Coulomb term
c
      if(ocoul)then
         gik=val
         val2=gik+gik
         val4=val2+val2
         fock(itr12) = facij*val4*p(itr34) + fock(itr12)
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
         endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
         gik=val*facex
         val2=gik+gik
         val4=val2+val2
         gil=gik
         if(i1.eq.i3.or.i2.eq.i4)gik=val2
         if(i2.eq.i3)gil=val2
         if(i2.ge.i3)goto 1
         itr23=ikyk+i2
         if(i2.ge.i4)goto 1
         itr24=iky(i4)+i2
 1       ajk=fock(itr23)-gil*p(itr14)
         bjk=ak(itr23)+gil*q(itr14)
         ail=fock(itr14)-gil*p(itr23)
         bil=ak(itr14)+gil*q(itr23)
         aik=fock(itr13)-gik*p(itr24)
         bik=ak(itr13)+gik*q(itr24)
         fock(itr24)=fock(itr24)-gik*p(itr13)
         ak(itr24)=ak(itr24)+gik*q(itr13)
         fock(itr23)=ajk
         ak(itr23)=bjk
         fock(itr14)=ail
         ak(itr14)=bil
         fock(itr13)=aik
         ak(itr13)=bik
      endif
_ELSE
c
c  non-dft UHF version
c
      aij=val4*p(itr34)+fock(itr12)
      fock(itr34)=val4*p(itr12)+fock(itr34)
      fock(itr12)=aij
c... exchange
      gil=gik
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)goto 1
      itr23=ikyk+i2
      if(i2.ge.i4)goto 1
      itr24=iky(i4)+i2
1     ajk=fock(itr23)-gil*p(itr14)
      bjk=ak(itr23)+gil*q(itr14)
      ail=fock(itr14)-gil*p(itr23)
      bil=ak(itr14)+gil*q(itr23)
      aik=fock(itr13)-gik*p(itr24)
      bik=ak(itr13)+gik*q(itr24)
      fock(itr24)=fock(itr24)-gik*p(itr13)
      ak(itr24)=ak(itr24)+gik*q(itr13)
      fock(itr23)=ajk
      ak(itr23)=bjk
      fock(itr14)=ail
      ak(itr14)=bil
      fock(itr13)=aik
      ak(itr13)=bik
_ENDIF
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_IF()
      subroutine dir_build_uhf70_orig(fock,ak,p,q,gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/shlg70)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),fock(*),ak(*),p(*),q(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell UHF Fock builder 
c
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=gik+gik
      val4=val2+val2
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      aij=val4*p(itr34)+fock(itr12)
      fock(itr34)=val4*p(itr12)+fock(itr34)
      fock(itr12)=aij
c... exchange
      gil=gik
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)goto 1
      itr23=ikyk+i2
      if(i2.ge.i4)goto 1
      itr24=iky(i4)+i2
1     ajk=fock(itr23)-gil*p(itr14)
      bjk=ak(itr23)+gil*q(itr14)
      ail=fock(itr14)-gil*p(itr23)
      bil=ak(itr14)+gil*q(itr23)
      aik=fock(itr13)-gik*p(itr24)
      bik=ak(itr13)+gik*q(itr24)
      fock(itr24)=fock(itr24)-gik*p(itr13)
      ak(itr24)=ak(itr24)+gik*q(itr13)
      fock(itr23)=ajk
      ak(itr23)=bjk
      fock(itr14)=ail
      ak(itr14)=bil
      fock(itr13)=aik
      ak(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_ENDIF
      subroutine dir_build_open_70(coul,exch,dens,gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/shlg70)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell SCF J & K builder (nshell=1)
c
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        gik=val
        val2=val+val
        ikyi=iky(i1)
        ikyj=iky(i2)
        ikyk=iky(i3)
        itr13=ikyi+i3
        itr14=ikyi+i4
        itr12=ikyi+i2
        itr23=ikyj+i3
        itr24=ikyj+i4
        itr34=ikyk+i4
        gil=val
        if(i1.eq.i3.or.i2.eq.i4)gik=val2
        if(i2.eq.i3)gil=val2
        if(i2.ge.i3)goto 280
        itr23=ikyk+i2
        if(i2.ge.i4)goto 280
        itr24=iky(i4)+i2
  280   bij=val2*dens(itr34)+coul(itr12)
        coul(itr34)=val2*dens(itr12)+coul(itr34)
        coul(itr12)=bij
        bjk=exch(itr23)+gil*dens(itr14)
        bil=exch(itr14)+gil*dens(itr23)
        bik=exch(itr13)+gik*dens(itr24)
        exch(itr24)=exch(itr24)+gik*dens(itr13)
        exch(itr23)=bjk
        exch(itr14)=bil
        exch(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open2_70(l2,coul,exch,dens,gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/shlg70)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell SCF J & K builder (nshell > 1)
c
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=val+val
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      gil=val
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)go to 280
      itr23=ikyk+i2
      if(i2.ge.i4)go to 280
      itr24=iky(i4)+i2
  280 continue
_IF1(x)c$dir scalar
_IF1(c)cdir$ novector
      do 300 iiii=1,nsheld
      bij=val2*dens(itr34)+coul(itr12)
      coul(itr34)=val2*dens(itr12)+coul(itr34)
      coul(itr12)=bij
      bjk=exch(itr23)+gil*dens(itr14)
      bil=exch(itr14)+gil*dens(itr23)
      bik=exch(itr13)+gik*dens(itr24)
      exch(itr24)=exch(itr24)+gik*dens(itr13)
      exch(itr23)=bjk
      exch(itr14)=bil
      exch(itr13)=bik
      itr12=itr12+l2
      itr34=itr34+l2
      itr13=itr13+l2
      itr14=itr14+l2
      itr23=itr23+l2
      itr24=itr24+l2
  300 continue
_IF1(c)cdir$ vector
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine sp0001(gout)
c        *****  special fast routine for -p- loop for 0001 ****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/pqgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(334),a1(333),b1(333),c1(334),
     +  a2(4000)
      dimension gout(*)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
      data sixty,tenm12/60.0d0,1.0d-12/
c
      gout1 = 0.0d0
      gout2 = 0.0d0
      gout3 = 0.0d0
      gout4 = 0.0d0
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.0d0/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = 1.0d0
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      g = 1.d0/(ep(i)+ecd)
      p = (pqab*pqab+qperp2)*g
      if (p .le. auxvar) go to 140
      q0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      q1 = 0.5d0*q0/p
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      q0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      q1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*q
  160 u = g*q1
      h0000 = h0000+q0
      h0001 = h0001+u
      h0003 = h0003-u*pqab
  180 continue
      h0001 = h0001*ecd*qperp
      h0003 = h0003*ecd
      p = dq*h0000
      g0001 = h0001*cosp+p*sing
      g0002 = h0001*sinp
      g0003 = h0003+p*cosg
      gout1 = gout1+dq00*h0000
      gout2 = gout2+dq01*g0001
      gout3 = gout3+dq01*g0002
      gout4 = gout4+dq01*g0003
 940  continue
      t1 = gout2
      t2 = gout3
      t3 = gout4
      gout(1) = gout1
      gout(2) = p11*t1+p21*t2+p31*t3
      gout(3) = p12*t1+p22*t2+p32*t3
      gout(4) = p13*t1+p23*t2+p33*t3
      return
      end
_EXTRACT(sp0011,pclinux)
      subroutine sp0011(gout)
c        *****  special fast routine for -p- loop for 0011 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
c
      dimension gout(*)
c
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
 660  h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0011 = 0.d0
      h0013 = 0.d0
      h0033 = 0.d0
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = g*(pqab2+qperp2)
      g = g*ecd
      if (p .le. auxvar) go to 140
      f0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 h0000 = h0000+f0
      h0001 = h0001+f1
      h0003 = h0003-f1*pqab
      h0011 = h0011+f2
      h0013 = h0013-f2*pqab
      h0033 = h0033+f2*pqab2
  180 continue
      h0022 = 0.5d0*ecd*(h0000-h0001)
      h0001 = h0001*qperp
      h0011 = h0011*qperp2+h0022
      h0013 = h0013*qperp
      h0033 = h0033+h0022
      if(sinp)120,100,120
 100  if(cosp)1000,120,920
 120  v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      g0011 = v44*h0011+v47*h0022
      g0012 = v54*h0011+v57*h0022
      g0022 = v74*h0011+v77*h0022
      g0013 = cosp*h0013
      g0023 = sinp*h0013
      g0033 = h0033
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0000 = h0000
      go to 2000
  920 g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
      go to 2000
1000  g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = -h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      g0010 = g0001
      g0020 = g0002
      g0021 = g0012
      g0030 = g0003
      g0031 = g0013
      g0032 = g0023
      if (rcdsq) 220,220,200
 200  g0010 = g0010+r13*g0000
      g0011 = g0011+r13*g0001
      g0012 = g0012+r13*g0002
      g0013 = g0013+r13*g0003
      g0030 = g0030+r33*g0000
      g0031 = g0031+r33*g0001
      g0032 = g0032+r33*g0002
      g0033 = g0033+r33*g0003
      g0001 = g0001+r14*g0000
      g0011 = g0011+r14*g0010
      g0021 = g0021+r14*g0020
      g0031 = g0031+r14*g0030
      g0003 = g0003+r34*g0000
      g0013 = g0013+r34*g0010
      g0023 = g0023+r34*g0020
      g0033 = g0033+r34*g0030
220   gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 5) = gout( 5)+g0010*dq10
      gout( 6) = gout( 6)+g0011*dq11
      gout( 7) = gout( 7)+g0012*dq11
      gout( 8) = gout( 8)+g0013*dq11
      gout( 9) = gout( 9)+g0020*dq10
      gout( 10) = gout( 10)+g0021*dq11
      gout( 11) = gout( 11)+g0022*dq11
      gout( 12) = gout( 12)+g0023*dq11
      gout( 13) = gout( 13)+g0030*dq10
      gout( 14) = gout( 14)+g0031*dq11
      gout( 15) = gout( 15)+g0032*dq11
      gout( 16) = gout( 16)+g0033*dq11
 940  continue
      ind = 0
      do 700 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
  700 continue
      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
      do 720 k = 1,4
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i3 ) = p13*t1+p23*t2+p33*t3
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
  720 continue
      return
      end
_ENDEXTRACT
      subroutine sp0101(gout)
c        *****  special fast routine for -p- loop for 0101 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/astore)
INCLUDE(common/shllfo)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      dimension gout(*)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/const)
INCLUDE(common/geom)
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0100 = 0.d0
      h0101 = 0.d0
      h0103 = 0.d0
      h0300 = 0.d0
      h0301 = 0.d0
      h0303 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = (pqab2+qperp2)*g
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 continue
      g03 = -pqab*f1
      h0000 = h0000+f0 *dp00
      h0001 = h0001+f1 *dp00
      h0003 = h0003+g03*dp00
      h0100 = h0100-f1
      h0101 = h0101-f2
      h0103 = h0103+pqab*f2
      h0300 = h0300-g03+bp*f0
      h0301 = h0301+bp*f1
      h0303 = h0303-pqab2*f2+bp*g03
  180 continue
      p = qperp*ecd
      h0001 = h0001*p
      h0003 = h0003*ecd
      h0202 = -0.5d0*ecd*h0100
      h0100 = h0100*qperp
      h0101 = h0101*qperp2*ecd
      h0103 = h0103*p
      h0301 = h0301*p
      h0303 = h0303*ecd
      h0301 = h0301+h0103
      h0101 = h0101+h0202
      h0303 = h0303+h0202
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      g0101 = cosp*h0101
      g0102 = sinp*h0101
      g0201 = u12*h0202
      g0202 = cosp*h0202
      g0301 = cosp*h0301
      g0302 = sinp*h0301
      g0303 = h0303
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0300 = h0300
      g0000 = h0000
      h0101 = g0101
      h0102 = g0102
      h0201 = g0201
      h0202 = g0202
      g0101 = cosp*h0101+u12*h0201
      g0102 = cosp*h0102+u12*h0202
      g0103 = cosp*h0103
      g0201 = sinp*h0101+cosp*h0201
      g0202 = sinp*h0102+cosp*h0202
      g0203 = sinp*h0103
      g0100 = cosp*h0100
      g0200 = sinp*h0100
      go to 2000
  920 g0100 = h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      go to 2000
 1000 g0100 = -h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = -h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = -h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
2000  r14 = dq*sing
      r34 = dq*cosg
      if (rcdsq) 720,720,700
  700 g0001 = g0001+r14*g0000
      g0101 = g0101+r14*g0100
      g0201 = g0201+r14*g0200
      g0301 = g0301+r14*g0300
      g0003 = g0003+r34*g0000
      g0103 = g0103+r34*g0100
      g0203 = g0203+r34*g0200
      g0303 = g0303+r34*g0300
  720 gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 17) = gout( 17)+g0100*dq00
      gout( 18) = gout( 18)+g0101*dq01
      gout( 19) = gout( 19)+g0102*dq01
      gout( 20) = gout( 20)+g0103*dq01
      gout( 33) = gout( 33)+g0200*dq00
      gout( 34) = gout( 34)+g0201*dq01
      gout( 35) = gout( 35)+g0202*dq01
      gout( 36) = gout( 36)+g0203*dq01
      gout( 49) = gout( 49)+g0300*dq00
      gout( 50) = gout( 50)+g0301*dq01
      gout( 51) = gout( 51)+g0302*dq01
      gout( 52) = gout( 52)+g0303*dq01
 940  continue
c ***
c ***
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 0101
c
c
c
_IFN1(f)      ind = 0
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 1 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   1  continue
_IFN1(f)      ind = -15
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 j = 1,4
_IFN1(f)      ind = ind+16
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      call mvml3(p11,1,gout(2),1,16,gout(2),1,16,4)
      return
      end
_EXTRACT(sp0111,ultra)
      subroutine sp0111(gout)
c        *****  special fast routine for -p- loop for 0111 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension g(64),h(64)
INCLUDE(common/shllfo)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/geom)
INCLUDE(common/const)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
c
      dimension gout(*)
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      eab = ep(i)
      dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(pqab2+qperp2)
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
      f3 = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)+theta4*c3(n+2))*gggy
  160 continue
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f2pqa2 = f2*pqab2
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      q1 = q1+f0 *bp
      q2 = q2+f1 *bp
      q3 = q3+f2 *bp
      q4 = q4+f1pqab*bp
      q5 = q5+f2pqab*bp
      q6 = q6+f2pqa2*bp
      r1 = r1+f1
      r2 = r2+f2
      r3 = r3+f3
      r4 = r4+f1pqab
      r5 = r5+f2pqab
      r6 = r6+f3pqab
      r7 = r7+f2pqa2
      r8 = r8+f3*pqab2
      r9 = r9+f3pqab*pqab2
  180 continue
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      qecd = qperp*ecd
      qecd2 = qperp*ecd2
      q2ecd = qperp2*ecd
      q2ecd2 = qperp2*ecd2
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = 0.5d0*qecd2*r2
      h( 44) = -0.5d0*ecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q2ecd2*qperp*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      g( 22) = v44*h( 22)+v47*h( 27)
      g( 23) = v54*h( 22)+v57*h( 27)
      g( 27) = v74*h( 22)+v77*h( 27)
      g( 24) = cosp*h( 24)
      g( 28) = sinp*h( 24)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g( 64) = h( 64)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 16) = h( 16)
      g( 18) = cosp*h( 18)
      g( 19) = sinp*h( 18)
      g( 20) = h( 20)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g( 52) = h( 52)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g(  4) = h(  4)
      g( 49) = h( 49)
      g(  1) = h(  1)
      h( 22) = g( 22)
      h( 23) = g( 23)
      h( 24) = g( 24)
      h( 27) = g( 27)
      h( 28) = g( 28)
      h( 38) = g( 38)
      h( 39) = g( 39)
      h( 40) = g( 40)
      h( 43) = g( 43)
      h( 44) = g( 44)
      h( 18) = g( 18)
      h( 19) = g( 19)
      h( 34) = g( 34)
      h( 35) = g( 35)
      g( 22) = cosp*h( 22)+u12*h( 38)
      g( 23) = cosp*h( 23)+u12*h( 39)
      g( 24) = cosp*h( 24)+u12*h( 40)
      g( 27) = cosp*h( 27)+u12*h( 43)
      g( 28) = cosp*h( 28)+u12*h( 44)
      g( 32) = cosp*h( 32)
      g( 38) = sinp*h( 22)+cosp*h( 38)
      g( 39) = sinp*h( 23)+cosp*h( 39)
      g( 40) = sinp*h( 24)+cosp*h( 40)
      g( 43) = sinp*h( 27)+cosp*h( 43)
      g( 44) = sinp*h( 28)+cosp*h( 44)
      g( 48) = sinp*h( 32)
      g( 18) = cosp*h( 18)+u12*h( 34)
      g( 19) = cosp*h( 19)+u12*h( 35)
      g( 20) = cosp*h( 20)
      g( 34) = sinp*h( 18)+cosp*h( 34)
      g( 35) = sinp*h( 19)+cosp*h( 35)
      g( 36) = sinp*h( 20)
      g( 17) = cosp*h( 17)
      g( 33) = sinp*h( 17)
      go to 2000
  920 g( 17) = h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = h( 20)
      g( 22) = h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = h( 27)
      g( 28) = dzero
      g( 32) = h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
      go to 2000
 1000 g( 17) = -h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = -h( 20)
      g( 22) = -h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = -h( 27)
      g( 28) = dzero
      g( 32) = -h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = -h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = -h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = -h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = -h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = -h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
2000  r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,50,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 720,720,700
700   do 701 kq1=1,49,16
          t1=g(kq1)
          t5=g(kq1+4)
          t13=g(kq1+12)
          t5=t5+r13*t1
          t13=t13+r33*t1
          t2=g(kq1+1)
          g(kq1+5)=g(kq1+5)+r13*t2+r14*t5
          g(kq1+13)=g(kq1+13)+r33*t2+r14*t13
          g(kq1+1)=t2+r14*t1
          t3=g(kq1+2)
          g(kq1+6)=g(kq1+6)+t3*r13
          g(kq1+14)=g(kq1+14)+t3*r33
          t9=g(kq1+8)
          g(kq1+9)=g(kq1+9)+r14*t9
          g(kq1+11)=g(kq1+11)+r34*t9
          t4=g(kq1+3)
          g(kq1+7)=g(kq1+7)+t4*r13+t5*r34
          g(kq1+4)=t5
          g(kq1+15)=g(kq1+15)+t4*r33+t13*r34
          g(kq1+12)=t13
701       g(kq1+3)=t4+t1*r34
720   do 721 kq1=2,62,4
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      dq01dd=dq01-dq11
      do 722 kq1=1,49,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq00
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq01dd
          gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq01dd
722       gout(kq1+3)=gout(kq1+3)+g(kq1+3)*dq01dd
      do 723 kq1=5,53,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq10
          gout(kq1+4)=gout(kq1+4)+g(kq1+4)*dq10
723       gout(kq1+8)=gout(kq1+8)+g(kq1+8)*dq10
  940 continue
_IFN1(f)      ind = 0
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 1 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   1  continue
_IFN1(f)      ind = -12
_IFN1(f)      do 2 j = 1,4
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 l = 1,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IFN1(f)      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 3 loopjk = 1,16
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   3  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      ind=5
_IF1(f)      do 4 j=1,4
_IF1(f)      call  mvml3(p11,1,gout(ind),4,1,gout(ind),4,1,4)
_IF1(f)  4   ind=ind+16
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,16)
      return
      end
_ENDEXTRACT
_EXTRACT(sp1111,ultra)
c ******************************************************
c ******************************************************
c             =   sp1111  =
c ******************************************************
c ******************************************************
      subroutine sp1111(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF1(cuf)      parameter (dzero=0.0e0)
_IFN1(cuf)      parameter (dzero=0.0d0)
      dimension g(256),h(256)
      dimension kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
INCLUDE(common/miscg)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
INCLUDE(common/pgeom)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/qgeom)
INCLUDE(common/const)
INCLUDE(common/maxc)
INCLUDE(common/astore)
      common/inttab/
     +  aa(333),ba(333),ca(333),abc1,
     +  ab(333),bb(333),cb(333),abc2,
     +  ac(333),bc(333),cc(333),abc3,
     +  ad(333),bd(333),cd(333),abc4,
     +  ae(333),be(333),ce(333),abc5,
     +  af(333),bf(333),cf(333),abc6
INCLUDE(common/auxvar)
c
      dimension gout(*)
c      data dzero/0.0e0/,done/1.0e0/
      data done/1.0d0/
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c ***
c *** this is the fps version of sp1111.
c ***
c *** as much code as possible reduced to loops (>=4)
c *** to avoid ps cache misses and to enhance compiler
c *** optimisation. will probably run like a drain on
c *** the cray-1s]
c ***
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
660   p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v4 = 0.d0
      v5 = 0.d0
      v6 = 0.d0
      w1 = 0.d0
      w2 = 0.d0
      w3 = 0.d0
      w4 = 0.d0
      w5 = 0.d0
      w6 = 0.d0
      w7 = 0.d0
      w8 = 0.d0
      w9 = 0.d0
      s1 = 0.d0
      s2 = 0.d0
      s3 = 0.d0
      s4 = 0.d0
      s6 = 0.d0
      s7 = 0.d0
      s8 = 0.d0
      s9 = 0.d0
      s10 = 0.d0
      s11 = 0.d0
      s12 = 0.d0
      s13 = 0.d0
      s14 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      t4 = 0.d0
      t5 = 0.d0
      t6 = 0.d0
      t7 = 0.d0
      t8 = 0.d0
      t9 = 0.d0
      t10 = 0.d0
      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
      t14 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0
      c5 = 0.d0
      c6 = 0.d0
_IF1(x)c$dir scalar
      do 180 ind = 1,ngangb
      isml = ismlq+ismlp(ind)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
       eab = ep(ind)
      dp00 = dp00p(ind)
      dp01 = dp01p(ind)
      dp10 = dp10p(ind)
      ap = app(ind)
      bp = bpp(ind)
      pqab = aqz-ap
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(qperp2+pqab2)
      if (p .le. auxvar) go to 140
      f0 = dsqrt(0.7853981625d0/(p*(gp(ind)+gcd)))*conp(ind)
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      f4 = 3.5d0*f3*gtx
      go to 160
  140 q = conp(ind)/dsqrt(gp(ind)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)+theta4*ca(n+2))*q
      f1 = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)+theta4*cb(n+2))*gy
      f2 = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)+theta4*cc(n+2))*ggy
      f3 = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)+theta4*cd(n+2))*gggy
      f4 = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)+theta4*ce(n+2))*gggy*
     &     gabcd
  160 apbp = ap*bp
      eab2 = eab*eab
      bpdp01 = bp*dp01
      apdp10 = ap*dp10
      edp01 = eab*dp01
      edp10 = eab*dp10
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f4pqab = f4*pqab
      f1pqa2 = f1*pqab2
      f2pqa2 = f2*pqab2
      f3pqa2 = f3*pqab2
      f4pqa2 = f4*pqab2
      f2pqa3 = f2pqa2*pqab
      f3pqa3 = f3pqa2*pqab
      f4pqa3 = f4pqa2*pqab
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      r1 = r1+f1 *edp01
      r2 = r2+f2 *edp01
      r3 = r3+f3 *edp01
      r4 = r4+f1pqab *edp01
      r5 = r5+f2pqab *edp01
      r6 = r6+f3pqab *edp01
      r7 = r7+f2pqa2 *edp01
      r8 = r8+f3pqa2 *edp01
      r9 = r9+f3pqa3 *edp01
      w1 = w1+f1 *edp10
      w2 = w2+f2 *edp10
      w3 = w3+f3 *edp10
      w4 = w4+f1pqab *edp10
      w5 = w5+f2pqab *edp10
      w6 = w6+f3pqab *edp10
      w7 = w7+f2pqa2 *edp10
      w8 = w8+f3pqa2 *edp10
      w9 = w9+f3pqa3 *edp10
      s1 = s1+f0 *eab
      s2 = s2+f1 *eab
      s3 = s3+f2 *eab
      s4 = s4+f3 *eab
      s6 = s6+f1pqab*eab
      s7 = s7+f2pqab*eab
      s8 = s8+f3pqab*eab
      s9 = s9+f1pqa2*eab
      s10 = s10+f2pqa2*eab
      s11 = s11+f3pqa2*eab
      s12 = s12+f2pqa3*eab
      s13 = s13+f3pqa3*eab
      s14 = s14+f3pqa3*pqab*eab
      t1 = t1+f0 *eab2
      t2 = t2+f1 *eab2
      t3 = t3+f2 *eab2
      t4 = t4+f3 *eab2
      t5 = t5+f4 *eab2
      t6 = t6+f2pqab*eab2
      t7 = t7+f3pqab*eab2
      t8 = t8+f4pqab*eab2
      t9 = t9+f2pqa2*eab2
      t10 = t10+f3pqa2*eab2
      t11 = t11+f4pqa2*eab2
      t12 = t12+f3pqa3*eab2
      t13 = t13+f4pqa3*eab2
      t14 = t14+f4pqa3*pqab*eab2
      if (rabsq .eq. 0.0d0) go to 180
      q1 = q1+f0 *bpdp01
      q2 = q2+f1 *bpdp01
      q3 = q3+f2 *bpdp01
      q4 = q4+f1pqab*bpdp01
      q5 = q5+f2pqab*bpdp01
      q6 = q6+f2pqa2*bpdp01
      v1 = v1+f0 *apdp10
      v2 = v2+f1 *apdp10
      v3 = v3+f2 *apdp10
      v4 = v4+f1pqab*apdp10
      v5 = v5+f2pqab*apdp10
      v6 = v6+f2pqa2*apdp10
      c1 = c1+f0 *apbp
      c2 = c2+f1 *apbp
      c3 = c3+f2 *apbp
      c4 = c4+f1pqab*apbp
      c5 = c5+f2pqab*apbp
      c6 = c6+f2pqa2*apbp
  180 continue
      a1 = aqz*s2-s6
      a2 = aqz*s3-s7
      a3 = aqz*s4-s8
      a4 = aqz*s6-s9
      a5 = aqz*s7-s10
      a6 = aqz*s8-s11
      a8 = aqz*s10-s12
      a9 = aqz*s11-s13
      a10 = aqz*s13-s14
      bqz = aqz-rab
      b1 = bqz*s2-s6
      b2 = bqz*s3-s7
      b3 = bqz*s4-s8
      b4 = bqz*s6-s9
      b5 = bqz*s7-s10
      b6 = bqz*s8-s11
      b8 = bqz*s10-s12
      b9 = bqz*s11-s13
      b10 = bqz*s13-s14
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      hecd2 = 0.5d0*ecd2
      qecd = qperp*ecd
      hqecd = 0.5d0*qecd
      qecd2 = qperp*ecd2
      hqecd2 = 0.5d0*qecd2
      q2ecd = qperp2*ecd
      q3ecd = qperp*q2ecd
      q2ecd2 = qperp2*ecd2
      q3ecd2 = q2ecd2*qperp
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = hqecd2*r2
      h( 44) = -hecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q3ecd2*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      h( 65) = -qperp*w1
      h(193) = w4+v1
      h(131) = hecd*w1
      h( 66) = h(131)-q2ecd*w2
      h( 68) = qecd*w5
      h(194) = h( 68)+qecd*v2
      h(196) = h(131)-ecd*w7-ecd*v4
      h(135) = hqecd2*w2
      h(140) = -hecd2*w5
      h( 75) = h(135)-qperp*h(131)
      h(203) = h(140)+hecd*(h(193)-ecd*v2)
      h( 72) = h(140)+q2ecd2*w6
      h(200) = h(135)-qecd2*(w8+v5)
      h( 70) = h( 75)+h(135)+h(135)-q3ecd2*w3
      h( 80) = h( 75)-qecd2*w8
      h(198) = h(203)+q2ecd2*(w6+v3)
      h(208) = h(203)+h(140)+h(140)+ecd2*(w9+v6)
      h(161) = 0.5d0*(s1-t2)
      h( 81) = h(161)+qperp2*t3
      h(113) = -qperp*(t6+b1)
      h(209) = -qperp*(t6+a1)
      h(241) = h(161)+t9+a4+b4+c1
      h(162) = hqecd*(s2-t3)
      h( 82) = h(162)-qecd*t3+q3ecd*t4
      temp = hecd*t6-q2ecd*t7
      h(114) = temp+hecd*b1-q2ecd*b2
      h(210) = temp+hecd*a1-q2ecd*a2
      h(242) = h(162)+qecd*(t10+a5+b5+c2)
      h( 99) = -hqecd*t3
      h(147) = h( 99)
      h(179) = hecd*(t6+b1)
      h(227) = hecd*(t6+a1)
      h(164) = hecd*(t6-s6)
      h( 84) = h(164)-q2ecd*t7
      temp = -hqecd*t3+qecd*t10
      h(116) = temp+qecd*b5
      h(212) = temp+qecd*a5
      h(244) = h(164)+ecd*(t6-t12-a8-b8-c4)+hecd*(a1+b1)
      h(103) = 0.25d0*ecd2*t3-0.5d0*q2ecd2*t4
      h(151) = h(103)
      h(183) = hqecd2*(t7+b2)
      h(231) = hqecd2*(t7+a2)
      h(108) = hqecd2*t7
      h(156) = h(108)
      h(188) = hecd2*(0.5d0*t3-t10-b5)
      h(236) = hecd2*(0.5d0*t3-t10-a5)
      hxxyy = 0.25d0*(ecd*(s1-t2)-ecd2*(s2-t3))
      h(171) = hxxyy+hecd2*t3
      h( 91) = hxxyy+0.5d0*(q2ecd*t3-q2ecd2*t4)
      temp = hqecd*(ecd*t7-t6)
      h(123) = temp+hqecd*(ecd*b2-b1)
      h(219) = temp+hqecd*(ecd*a2-a1)
      h(251) = hxxyy+hecd*(t9+a4+b4+c1)-hecd2*(t10+a5+b5+c2)
      h(166) = hxxyy+0.5d0*q2ecd2*(s3-t4)
      h( 86) = hxxyy+(hecd2+0.5d0*q2ecd)*t3+q2ecd2*(-3.d0*t4+
     +    0.5d0*s3+qperp2*t5)
      h(118) = 1.5d0*qecd2*(t7+b2)-hqecd*(t6+b1)-q3ecd2*(b3+t8)
      h(214) = 1.5d0*qecd2*(t7+a2)-hqecd*(t6+a1)-q3ecd2*(a3+t8)
      h(246) = hxxyy-hecd2*(qperp2*t4+t10+a5+b5)+hecd*(t9+a4+b4+c1-ecd*
     +    c2)+q2ecd2*(t11+0.5d0*s3+a6+b6+c3)
      h(168) = hqecd2*(t7-s7)
      h( 88) = 1.5d0*qecd2*t7-hqecd2*s7-q3ecd2*t8
      temp = hecd2*(0.5d0*t3-t10)+q2ecd2*(t11-0.5d0*t4)
      h(120) = temp-hecd2*b5+q2ecd2*b6
      h(216) = temp-hecd2*a5+q2ecd2*a6
      h(248) = qecd2*(1.5d0*t7-t13-a9-b9-c5)-hqecd2*(s7-a2-b2)
      h(176) = hxxyy+hecd2*(s10-t10)
      h( 96) = hxxyy-hecd2*(qperp2*t4+t10-s10)+0.5d0*q2ecd*t3+q2ecd2*
     +     t11
      h(128) = qecd2*(1.5d0*t7-t13-b9)-hqecd*(t6+b1)+hqecd2*b2
      h(224) = qecd2*(1.5d0*t7-t13-a9)-hqecd*(t6+a1)+hqecd2*a2
      h(256) = hxxyy+hecd2*(-3.d0*(a5+b5)+t3+s10-c2)+ecd2*(-3.d0*t10+
     +     t14+a10+b10+c6)+hecd*(t9+a4+b4+c1)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
_IF1(a)cvd$  shortloop
      do 103 kq1=22,214,48
          g(kq1  ) = v44*h(kq1) + v47*h(kq1+5)
          g(kq1+1) = v54*h(kq1) + v57*h(kq1+5)
103       g(kq1+5) = v74*h(kq1) + v77*h(kq1+5)
_IF1(a)cvd$  shortloop
      do 101 kq1=24,216,48
          g(kq1  ) = cosp*h(kq1)
101       g(kq1+4) = sinp*h(kq1)
_IF1(a)cvd$  shortloop
      do 102 kq1=18,210,48
          g(kq1  ) = cosp*h(kq1)
102       g(kq1+1) = sinp*h(kq1)
      g( 80) = h( 80)
      g( 86) = v44*h( 86)+v47*h( 91)
      g( 87) = v54*h( 86)+v57*h( 91)
      g( 91) = v74*h( 86)+v77*h( 91)
      g( 88) = cosp*h( 88)
      g( 92) = sinp*h( 88)
      g( 96) = h( 96)
      g(102) = v45*h(103)
      g(103) = v55*h(103)
      g(107) = -g(102)
      g(104) = u12*h(108)
      g(108) = cosp*h(108)
      g(112) = dzero
      g(128) = h(128)
      g(134) = v45*h(135)
      g(135) = v55*h(135)
      g(139) = -g(134)
      g(136) = u12*h(140)
      g(140) = cosp*h(140)
      g(144) = dzero
      g(150) = v45*h(151)
      g(151) = v55*h(151)
      g(155) = -g(150)
      g(152) = u12*h(156)
      g(156) = cosp*h(156)
      g(160) = dzero
      g(176) = h(176)
      g(182) = v45*h(183)
      g(183) = v55*h(183)
      g(187) = -g(182)
      g(184) = u12*h(188)
      g(188) = cosp*h(188)
      g(192) = dzero
      g(198) = v44*h(198)+v47*h(203)
      g(199) = v54*h(198)+v57*h(203)
      g(203) = v74*h(198)+v77*h(203)
      g(200) = cosp*h(200)
      g(204) = sinp*h(200)
      g(230) = v45*h(231)
      g(231) = v55*h(231)
      g(235) = -g(230)
      g(232) = u12*h(236)
      g(236) = cosp*h(236)
      g(240) = dzero
      g(246) = v44*h(246)+v47*h(251)
      g(247) = v54*h(246)+v57*h(251)
      g(251) = v74*h(246)+v77*h(251)
      g(248) = cosp*h(248)
      g(252) = sinp*h(248)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 48) = dzero
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 68) = h( 68)
      g( 82) = cosp*h( 82)
      g( 83) = sinp*h( 82)
      g( 84) = h( 84)
      g( 98) = u12*h( 99)
      g( 99) = cosp*h( 99)
      g(100) = dzero
      g(116) = h(116)
      g(130) = u12*h(131)
      g(131) = cosp*h(131)
      g(132) = dzero
      g(146) = u12*h(147)
      g(147) = cosp*h(147)
      g(148) = dzero
      g(164) = h(164)
      g(178) = u12*h(179)
      g(179) = cosp*h(179)
      g(180) = dzero
      g(194) = cosp*h(194)
      g(195) = sinp*h(194)
      g(226) = u12*h(227)
      g(227) = cosp*h(227)
      g(228) = dzero
      g(242) = cosp*h(242)
      g(243) = sinp*h(242)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 36) = dzero
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g( 65) = h( 65)
      g( 81) = h( 81)
      g( 97) = dzero
      g(113) = h(113)
      g(129) = dzero
      g(145) = dzero
      g(161) = h(161)
      g(177) = dzero
      g(225) = dzero
      g( 33) = dzero
      h( 80) = cosp*g( 80)
      h( 96) = cosp*g( 96)
      h(112) =           u12*g(176)
      h(128) = cosp*g(128)
      h(144) = sinp*g( 80)
      h(160) = sinp*g( 96)
      h(176) =           cosp*g(176)
      h(192) = sinp*g(128)
_IF1(a)cvd$  shortloop
      do 121 kq1=70,118,16
          h(kq1   ) = cosp*g(kq1   ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1   ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+ 1) + u12*g(kq1+65)
          h(kq1+65) = sinp*g(kq1+ 1) + cosp*g(kq1+65)
          h(kq1+ 2) = cosp*g(kq1+ 2) + u12*g(kq1+66)
          h(kq1+66) = sinp*g(kq1+ 2) + cosp*g(kq1+66)
          h(kq1+ 5) = cosp*g(kq1+ 5) + u12*g(kq1+69)
          h(kq1+69) = sinp*g(kq1+ 5) + cosp*g(kq1+69)
          h(kq1+ 6) = cosp*g(kq1+ 6) + u12*g(kq1+70)
          h(kq1+70) = sinp*g(kq1+ 6) + cosp*g(kq1+70)
121   continue
      h( 68) = cosp*g( 68)
      h( 84) = cosp*g( 84)
      h(100) =           u12*g(164)
      h(116) = cosp*g(116)
      h(132) = sinp*g( 68)
      h(148) = sinp*g( 84)
      h(164) =           cosp*g(164)
      h(180) = sinp*g(116)
_IF1(a)cvd$  shortloop
      do 122 kq1=66,114,16
          h(kq1   ) = cosp*g(kq1  ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1  ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+1) + u12*g(kq1+65)
122       h(kq1+65) = sinp*g(kq1+1) + cosp*g(kq1+65)
_IF1(a)cvd$  shortloop
      do 1221 kq1=2,50,16
          h(kq1   ) = g(kq1   )
          h(kq1+ 1) = g(kq1+ 1)
          h(kq1+ 4) = g(kq1+ 4)
          h(kq1+ 5) = g(kq1+ 5)
          h(kq1+ 6) = g(kq1+ 6)
1221      h(kq1+ 9) = g(kq1+ 9)
_IF1(a)cvd$  shortloop
      do 1222 kq1=12,60,16
          h(kq1    ) = g(kq1    )
          h(kq1+182) = g(kq1+182)
          h(kq1+183) = g(kq1+183)
          h(kq1+186) = g(kq1+186)
          h(kq1+187) = g(kq1+187)
1222      h(kq1+188) = g(kq1+188)
_IF1(a)cvd$  shortloop
      do 1223 kq1=203,251,16
          h(kq1  ) = g(kq1  )
1223      h(kq1+1) = g(kq1+1)
      h( 65) = cosp*g( 65)
      h( 81) = cosp*g( 81)
      h( 97) =           u12*g(161)
      h(113) = cosp*g(113)
      h(129) = sinp*g( 65)
      h(145) = sinp*g( 81)
      h(161) =           cosp*g(161)
      h(177) = sinp*g(113)
      h( 48) = g( 48)
      h( 36) = g( 36)
      h(228) = g(228)
      h(240) = g(240)
      h(225) = g(225)
      h( 33) = g( 33)
_IF1(a)cvd$  shortloop
      do 123 kq1=22,214,64
          g(kq1   ) = cosp*h(kq1   ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1   ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+ 1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+ 1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+ 2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+ 2) + cosp*h(kq1+18)
          g(kq1+ 5) = cosp*h(kq1+ 5) + u12* h(kq1+21)
          g(kq1+21) = sinp*h(kq1+ 5) + cosp*h(kq1+21)
          g(kq1+ 6) = cosp*h(kq1+ 6) + u12* h(kq1+22)
          g(kq1+22) = sinp*h(kq1+ 6) + cosp*h(kq1+22)
          g(kq1+10) = cosp*h(kq1+10) + u12* h(kq1+26)
123       g(kq1+26) = sinp*h(kq1+10) + cosp* h(kq1+26)
_IF1(a)cvd$  shortloop
      do 124 kq1=17,209,64
          g(kq1   ) = cosp*h(kq1  ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1  ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+2) + cosp*h(kq1+18)
          g(kq1+ 3) = cosp*h(kq1+3) + u12* h(kq1+19)
124       g(kq1+19) = sinp*h(kq1+3) + cosp* h(kq1+19)
_IF1(a)cvd$  concur
      do 125 kq1=49,177,64
          kkq1=kq1
_IF1(a)cvd$  shortloop
          do 126 kkkq1=1,32
              g(kkq1)=h(kkq1)
              kkq1=kkq1+1
126       continue
125   continue
_IF1(a)cvd$  shortloop
      do 127 kq1=1,16
          g(kq1)=h(kq1)
127       g(kq1+240)=h(kq1+240)
      goto 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
920   do 921 kkq1=1,10
          kq1=kq1off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) = h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) = h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
921       g(kq1+15) = h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
      do 922 kkq1=1,6
          kq1=kq2off(kkq1)
          g(kq1   ) = dzero
          g(kq1+ 1) = dzero
          g(kq1+ 2) = h(kq1+ 2)
          g(kq1+ 3) = dzero
          g(kq1+ 5) = dzero
          g(kq1+ 6) = h(kq1+ 6)
          g(kq1+ 7) = dzero
          g(kq1+10) = dzero
          g(kq1+11) = h(kq1+11)
922       g(kq1+15) = dzero
      go to 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
1000  do 1001 kkq1=1,6
          kq1=kq3off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) =-h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) =-h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
1001      g(kq1+15) = h(kq1+15)
      do 1002 kkq1=1,4
          kq1=kq4off(kkq1)
          g(kq1   ) = -h(kq1   )
          g(kq1+ 1) =  h(kq1+ 1)
          g(kq1+ 2) =  dzero
          g(kq1+ 3) = -h(kq1+ 3)
          g(kq1+ 5) = -h(kq1+ 5)
          g(kq1+ 6) =  dzero
          g(kq1+ 7) =  h(kq1+ 7)
          g(kq1+10) = -h(kq1+10)
          g(kq1+11) =  dzero
1002      g(kq1+15) = -h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
      do 1003 kkq1=1,6
          kq1=kq5off(kkq1)
          g(kq1   ) =  dzero
          g(kq1+ 1) =  dzero
          g(kq1+ 2) =  h(kq1+ 2)
          g(kq1+ 3) =  dzero
          g(kq1+ 5) =  dzero
          g(kq1+ 6) = -h(kq1+ 6)
          g(kq1+ 7) =  dzero
          g(kq1+10) =  dzero
          g(kq1+11) =  h(kq1+11)
1003      g(kq1+15) =  dzero
          g(99)=-g(99)
          g(108)=-g(108)
          g(147)=-g(147)
          g(156)=-g(156)
          g(103)=-g(103)
          g(151)=-g(151)
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
_IF1(a)cvd$  shortloop
      do 2001 kq1=2,242,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 1200,1200,1300
_IF1(a)cvd$  concur
1300  do 1301 kq1=1,4
          kkq1=kq1
_IF1(a)cvd$  shortloop
          do 1302 jq1=1,16
              g(kkq1+4) = r13*g(kkq1) + g(kkq1+4)
              g(kkq1+12)= r33*g(kkq1) + g(kkq1+12)
1302          kkq1=kkq1+16
1301  continue
c ***
      do 1303 kq1=1,253,4
          g(kq1+1) = r14*g(kq1) + g(kq1+1)
1303      g(kq1+3) = r34*g(kq1) + g(kq1+3)
c ***
c1200  do 1201 kq1=1,256
c1201      gout(kq1) = dq11*g(kq1) + gout(kq1)
 1200  continue
       call daxpy(256,dq11,g,1,gout,1)
      dq01x=dq01-dq11
_IF1(a)cvd$  vector
      do 1202 kq1=2,242,16
          gout(kq1  ) = dq01x*g(kq1  ) + gout(kq1  )
          gout(kq1+1) = dq01x*g(kq1+1) + gout(kq1+1)
1202      gout(kq1+2) = dq01x*g(kq1+2) + gout(kq1+2)
      dq10x=dq10-dq11
      dq00x=dq00-dq11
_IF1(a)cvd$  vector
      do 1203 kq1=1,241,16
          gout(kq1   ) = dq00x*g(kq1   ) + gout(kq1   )
          gout(kq1+ 4) = dq10x*g(kq1+ 4) + gout(kq1+ 4)
          gout(kq1+ 8) = dq10x*g(kq1+ 8) + gout(kq1+ 8)
1203      gout(kq1+12) = dq10x*g(kq1+12) + gout(kq1+12)
940   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 1111
c
c
_IFN1(f)      i1 = 64
_IFN1(f)      i2 = 128
_IFN1(f)      i3 = 192
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)      do 1 jkl = 1,64
_IFN1(f)      i1 = i1+1
_IFN1(f)      i2 = i2+1
_IFN1(f)      i3 = i3+1
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    1 continue
_IFN1(f)      ind = -48
_IFN1(f)      do 2 i = 1,4
_IFN1(f)      ind = ind+48
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    2 continue
_IFN1(f)      ind = -12
_IFN1(f)      do 3 loopij = 1,16
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 3 l = 1,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    3 continue
_IFN1(f)      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)      do 4 ijk = 1,64
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    4 continue
_IF1(f)      call mvml3(p11,1,gout(65),64,1,gout(65),64,1,64)
_IF1(f)      ind=17
_IF1(f)      do 5 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),16,1,gout(ind),16,1,16)
_IF1(f)   5  ind=ind+64
_IF1(f)      ind=5
_IF1(f)      do 6 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),4,16,gout(ind),4,16,16)
_IF1(f)   6  ind=ind+1
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,64)
c
      return
      end
_ENDEXTRACT
_IF(newints)
      subroutine sp0002(gout)
c        *****  special fast routine for -p- loop for 0001 ****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/pqgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(334),a1(333),b1(333),c1(334),
     +  a2(4000)
      dimension gout(*)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
      data sixty,tenm12/60.0d0,1.0d-12/
c
      gout1 = 0.0d0
      gout2 = 0.0d0
      gout3 = 0.0d0
      gout4 = 0.0d0
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.0d0/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
c     sd = csd(l)
c     pc = cpc(k)
      pd = cpd(l)
c     dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = 1.0d0
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      g = 1.d0/(ep(i)+ecd)
      p = (pqab*pqab+qperp2)*g
      if (p .le. auxvar) go to 140
      q0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      q1 = 0.5d0*q0/p
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      q0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      q1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*q
  160 u = g*q1
      h0000 = h0000+q0
      h0001 = h0001+u
      h0003 = h0003-u*pqab
  180 continue
      h0001 = h0001*ecd*qperp
      h0003 = h0003*ecd
      p = dq*h0000
      g0001 = h0001*cosp+p*sing
      g0002 = h0001*sinp
      g0003 = h0003+p*cosg
c     gout1 = gout1+dq00*h0000
      gout2 = gout2+dq01*g0001
      gout3 = gout3+dq01*g0002
      gout4 = gout4+dq01*g0003
 940  continue
      t1 = gout2
      t2 = gout3
      t3 = gout4
c     gout(1) = gout1
      gout(2) = p11*t1+p21*t2+p31*t3
      gout(3) = p12*t1+p22*t2+p32*t3
      gout(4) = p13*t1+p23*t2+p33*t3
      return
      end
      subroutine sp0022(gout)
c        *****  special fast routine for -p- loop for 0011 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
c
      dimension gout(*)
c
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 continue
c     sc = csc(k)
c     sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
c     dq00 = sc*sd*v
c     dq01 = sc*pd*v
c     dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
 660  h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0011 = 0.d0
      h0013 = 0.d0
      h0033 = 0.d0
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = g*(pqab2+qperp2)
      g = g*ecd
      if (p .le. auxvar) go to 140
      f0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 h0000 = h0000+f0
      h0001 = h0001+f1
      h0003 = h0003-f1*pqab
      h0011 = h0011+f2
      h0013 = h0013-f2*pqab
      h0033 = h0033+f2*pqab2
  180 continue
      h0022 = 0.5d0*ecd*(h0000-h0001)
      h0001 = h0001*qperp
      h0011 = h0011*qperp2+h0022
      h0013 = h0013*qperp
      h0033 = h0033+h0022
      if(sinp)120,100,120
 100  if(cosp)1000,120,920
 120  v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      g0011 = v44*h0011+v47*h0022
      g0012 = v54*h0011+v57*h0022
      g0022 = v74*h0011+v77*h0022
      g0013 = cosp*h0013
      g0023 = sinp*h0013
      g0033 = h0033
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0000 = h0000
      go to 2000
  920 g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
      go to 2000
1000  continue
      g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = -h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      g0010 = g0001
      g0020 = g0002
      g0021 = g0012
      g0030 = g0003
      g0031 = g0013
      g0032 = g0023
      if (rcdsq) 220,220,200
 200  continue
      g0010 = g0010+r13*g0000
      g0011 = g0011+r13*g0001
      g0012 = g0012+r13*g0002
      g0013 = g0013+r13*g0003
      g0030 = g0030+r33*g0000
      g0031 = g0031+r33*g0001
      g0032 = g0032+r33*g0002
      g0033 = g0033+r33*g0003
c     g0001 = g0001+r14*g0000
      g0011 = g0011+r14*g0010
      g0021 = g0021+r14*g0020
      g0031 = g0031+r14*g0030
c     g0003 = g0003+r34*g0000
      g0013 = g0013+r34*g0010
      g0023 = g0023+r34*g0020
      g0033 = g0033+r34*g0030
220   continue
c     gout( 1) = gout( 1)+g0000*dq00
c     gout( 2) = gout( 2)+g0001*dq01
c     gout( 3) = gout( 3)+g0002*dq01
c     gout( 4) = gout( 4)+g0003*dq01
c     gout( 5) = gout( 5)+g0010*dq10
      gout( 6) = gout( 6)+g0011*dq11
      gout( 7) = gout( 7)+g0012*dq11
      gout( 8) = gout( 8)+g0013*dq11
c     gout( 9) = gout( 9)+g0020*dq10
      gout( 10) = gout( 10)+g0021*dq11
      gout( 11) = gout( 11)+g0022*dq11
      gout( 12) = gout( 12)+g0023*dq11
c     gout( 13) = gout( 13)+g0030*dq10
      gout( 14) = gout( 14)+g0031*dq11
      gout( 15) = gout( 15)+g0032*dq11
      gout( 16) = gout( 16)+g0033*dq11
 940  continue
      ind = 1
      do 700 l = 2,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
  700 continue
      ind = 1
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
      do 720 k = 2,4
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i3 ) = p13*t1+p23*t2+p33*t3
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
  720 continue
      return
      end
      subroutine sp0202(gout)
c        *****  special fast routine for -p- loop for 0101 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/astore)
INCLUDE(common/shllfo)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      dimension gout(*)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/const)
INCLUDE(common/geom)
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 continue
      sc = csc(k)
c     sd = csd(l)
c     pc = cpc(k)
      pd = cpd(l)
c     dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 continue
c     h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0100 = 0.d0
      h0101 = 0.d0
      h0103 = 0.d0
      h0300 = 0.d0
      h0301 = 0.d0
      h0303 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
c     dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = (pqab2+qperp2)*g
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 continue
      g03 = -pqab*f1
c     h0000 = h0000+f0 *dp00
c     h0001 = h0001+f1 *dp00
c     h0003 = h0003+g03*dp00
      h0100 = h0100-f1
      h0101 = h0101-f2
      h0103 = h0103+pqab*f2
      h0300 = h0300-g03+bp*f0
      h0301 = h0301+bp*f1
      h0303 = h0303-pqab2*f2+bp*g03
  180 continue
      p = qperp*ecd
      h0001 = h0001*p
      h0003 = h0003*ecd
      h0202 = -0.5d0*ecd*h0100
      h0100 = h0100*qperp
      h0101 = h0101*qperp2*ecd
      h0103 = h0103*p
      h0301 = h0301*p
      h0303 = h0303*ecd
      h0301 = h0301+h0103
      h0101 = h0101+h0202
      h0303 = h0303+h0202
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      g0101 = cosp*h0101
      g0102 = sinp*h0101
      g0201 = u12*h0202
      g0202 = cosp*h0202
      g0301 = cosp*h0301
      g0302 = sinp*h0301
      g0303 = h0303
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0300 = h0300
c     g0000 = h0000
      h0101 = g0101
      h0102 = g0102
      h0201 = g0201
      h0202 = g0202
      g0101 = cosp*h0101+u12*h0201
      g0102 = cosp*h0102+u12*h0202
      g0103 = cosp*h0103
      g0201 = sinp*h0101+cosp*h0201
      g0202 = sinp*h0102+cosp*h0202
      g0203 = sinp*h0103
      g0100 = cosp*h0100
      g0200 = sinp*h0100
      go to 2000
  920 g0100 = h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = h0301
      g0302 = dzero
      g0303 = h0303
c     g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      go to 2000
 1000 g0100 = -h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = -h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = -h0301
      g0302 = dzero
      g0303 = h0303
c     g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
2000  r14 = dq*sing
      r34 = dq*cosg
      if (rcdsq) 720,720,700
  700 continue
c     g0001 = g0001+r14*g0000
      g0101 = g0101+r14*g0100
      g0201 = g0201+r14*g0200
      g0301 = g0301+r14*g0300
c     g0003 = g0003+r34*g0000
      g0103 = g0103+r34*g0100
      g0203 = g0203+r34*g0200
      g0303 = g0303+r34*g0300
  720 continue
c     gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
c     gout( 17) = gout( 17)+g0100*dq00
      gout( 18) = gout( 18)+g0101*dq01
      gout( 19) = gout( 19)+g0102*dq01
      gout( 20) = gout( 20)+g0103*dq01
c     gout( 33) = gout( 33)+g0200*dq00
      gout( 34) = gout( 34)+g0201*dq01
      gout( 35) = gout( 35)+g0202*dq01
      gout( 36) = gout( 36)+g0203*dq01
c     gout( 49) = gout( 49)+g0300*dq00
      gout( 50) = gout( 50)+g0301*dq01
      gout( 51) = gout( 51)+g0302*dq01
      gout( 52) = gout( 52)+g0303*dq01
 940  continue
c ***
c ***
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 0101
c
c
c
_IFN1(f)      ind = 1
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 1 loopkl = 2,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   1  continue
_IFN1(f)      ind = 1
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 j = 2,4
_IFN1(f)      ind = ind+16
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      call mvml3(p11,1,gout(2),1,16,gout(2),1,16,4)
      return
      end
_EXTRACT(sp0111,ultra)
      subroutine sp0222(gout)
c        *****  special fast routine for -p- loop for 0111 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension g(64),h(64)
INCLUDE(common/shllfo)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/geom)
INCLUDE(common/const)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
c
      dimension gout(*)
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 continue
c     sc = csc(k)
c     sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
c     dq00 = sc*sd*v
c     dq01 = sc*pd*v
c     dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 continue
c     p1 = 0.d0
c     p2 = 0.d0
c     p3 = 0.d0
c     p4 = 0.d0
c     p5 = 0.d0
c     p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      eab = ep(i)
c     dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(pqab2+qperp2)
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
      f3 = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)+theta4*c3(n+2))*gggy
  160 continue
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f2pqa2 = f2*pqab2
c     p1 = p1+f0 *dp00
c     p2 = p2+f1 *dp00
c     p3 = p3+f2 *dp00
c     p4 = p4+f1pqab*dp00
c     p5 = p5+f2pqab*dp00
c     p6 = p6+f2pqa2*dp00
      q1 = q1+f0 *bp
      q2 = q2+f1 *bp
      q3 = q3+f2 *bp
      q4 = q4+f1pqab*bp
      q5 = q5+f2pqab*bp
      q6 = q6+f2pqa2*bp
      r1 = r1+f1
      r2 = r2+f2
      r3 = r3+f3
      r4 = r4+f1pqab
      r5 = r5+f2pqab
      r6 = r6+f3pqab
      r7 = r7+f2pqa2
      r8 = r8+f3*pqab2
      r9 = r9+f3pqab*pqab2
  180 continue
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      qecd = qperp*ecd
      qecd2 = qperp*ecd2
      q2ecd = qperp2*ecd
      q2ecd2 = qperp2*ecd2
c     h(  1) = p1
c     h(  2) = qecd*p2
c     h(  4) = -ecd*p4
c     h( 11) = hecd*(p1-ecd*p2)
c     h(  6) = h( 11)+q2ecd2*p3
c     h(  8) = -qecd2*p5
c     h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = 0.5d0*qecd2*r2
      h( 44) = -0.5d0*ecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q2ecd2*qperp*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      g( 22) = v44*h( 22)+v47*h( 27)
      g( 23) = v54*h( 22)+v57*h( 27)
      g( 27) = v74*h( 22)+v77*h( 27)
      g( 24) = cosp*h( 24)
      g( 28) = sinp*h( 24)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g( 64) = h( 64)
c     g(  6) = v44*h(  6)+v47*h( 11)
c     g(  7) = v54*h(  6)+v57*h( 11)
c     g( 11) = v74*h(  6)+v77*h( 11)
c     g(  8) = cosp*h(  8)
c     g( 12) = sinp*h(  8)
c     g( 16) = h( 16)
      g( 18) = cosp*h( 18)
      g( 19) = sinp*h( 18)
      g( 20) = h( 20)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g( 52) = h( 52)
c     g(  2) = cosp*h(  2)
c     g(  3) = sinp*h(  2)
c     g(  4) = h(  4)
      g( 49) = h( 49)
c     g(  1) = h(  1)
      h( 22) = g( 22)
      h( 23) = g( 23)
      h( 24) = g( 24)
      h( 27) = g( 27)
      h( 28) = g( 28)
      h( 38) = g( 38)
      h( 39) = g( 39)
      h( 40) = g( 40)
      h( 43) = g( 43)
      h( 44) = g( 44)
      h( 18) = g( 18)
      h( 19) = g( 19)
      h( 34) = g( 34)
      h( 35) = g( 35)
      g( 22) = cosp*h( 22)+u12*h( 38)
      g( 23) = cosp*h( 23)+u12*h( 39)
      g( 24) = cosp*h( 24)+u12*h( 40)
      g( 27) = cosp*h( 27)+u12*h( 43)
      g( 28) = cosp*h( 28)+u12*h( 44)
      g( 32) = cosp*h( 32)
      g( 38) = sinp*h( 22)+cosp*h( 38)
      g( 39) = sinp*h( 23)+cosp*h( 39)
      g( 40) = sinp*h( 24)+cosp*h( 40)
      g( 43) = sinp*h( 27)+cosp*h( 43)
      g( 44) = sinp*h( 28)+cosp*h( 44)
      g( 48) = sinp*h( 32)
      g( 18) = cosp*h( 18)+u12*h( 34)
      g( 19) = cosp*h( 19)+u12*h( 35)
      g( 20) = cosp*h( 20)
      g( 34) = sinp*h( 18)+cosp*h( 34)
      g( 35) = sinp*h( 19)+cosp*h( 35)
      g( 36) = sinp*h( 20)
      g( 17) = cosp*h( 17)
      g( 33) = sinp*h( 17)
      go to 2000
  920 continue
      g( 17) = h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = h( 20)
      g( 22) = h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = h( 27)
      g( 28) = dzero
      g( 32) = h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
c     g(  1) = h(  1)
c     g(  2) = h(  2)
c     g(  3) = dzero
c     g(  4) = h(  4)
c     g(  6) = h(  6)
c     g(  7) = dzero
c     g(  8) = h(  8)
c     g( 11) = h( 11)
c     g( 12) = dzero
c     g( 16) = h( 16)
      go to 2000
 1000 g( 17) = -h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = -h( 20)
      g( 22) = -h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = -h( 27)
      g( 28) = dzero
      g( 32) = -h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = -h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = -h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = -h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
c     g(  1) = h(  1)
c     g(  2) = -h(  2)
c     g(  3) = dzero
c     g(  4) = h(  4)
c     g(  6) = h(  6)
c     g(  7) = dzero
c     g(  8) = -h(  8)
c     g( 11) = h( 11)
c     g( 12) = dzero
c     g( 16) = h( 16)
2000  r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
c     do 2001 kq1=2,50,16
      do 2001 kq1=18,50,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 720,720,700
700   do 701 kq1=17,49,16
          t1=g(kq1)
          t5=g(kq1+4)
          t13=g(kq1+12)
          t5=t5+r13*t1
          t13=t13+r33*t1
          t2=g(kq1+1)
          g(kq1+5)=g(kq1+5)+r13*t2+r14*t5
          g(kq1+13)=g(kq1+13)+r33*t2+r14*t13
c         g(kq1+1)=t2+r14*t1
          t3=g(kq1+2)
          g(kq1+6)=g(kq1+6)+t3*r13
          g(kq1+14)=g(kq1+14)+t3*r33
          t9=g(kq1+8)
          g(kq1+9)=g(kq1+9)+r14*t9
          g(kq1+11)=g(kq1+11)+r34*t9
          t4=g(kq1+3)
          g(kq1+7)=g(kq1+7)+t4*r13+t5*r34
c         g(kq1+4)=t5
          g(kq1+15)=g(kq1+15)+t4*r33+t13*r34
c         g(kq1+12)=t13
c         g(kq1+3)=t4+t1*r34
701   continue
c720   do 721 kq1=2,62,4
720   continue
c      do 721 kq1=22,62,4
c          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
c          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
c721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      do ii=0,2
      do 721 kq1=22+ii*16,30+ii*16,4
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      enddo
c     dq01dd=dq01-dq11
c     dq01dd=-dq11
c      do 722 kq1=1,49,16
c          gout(kq1  )=gout(kq1  )+g(kq1  )*dq00
c          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq01dd
c          gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq01dd
c722       gout(kq1+3)=gout(kq1+3)+g(kq1+3)*dq01dd
c      do 723 kq1=5,53,16
c          gout(kq1  )=gout(kq1  )+g(kq1  )*dq10
c          gout(kq1+4)=gout(kq1+4)+g(kq1+4)*dq10
c723       gout(kq1+8)=gout(kq1+8)+g(kq1+8)*dq10
  940 continue
_IFN1(f)c     ind = 0
_IFN1(f)      ind = 4
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)c     do 1 loopkl = 1,16,4
_IFN1(f)      do 1 loopkl = 5,16,4
_IFN1(f)      ind = ind+1
_IFN1(f)      do ikl=2,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)      enddo
_IFN1(f)   1  continue
_IFN1(f)c     ind = -12
_IFN1(f)c     do 2 j = 1,4
_IFN1(f)      ind = 4
_IFN1(f)      do 2 j = 2,4
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      ind = ind+1
_IFN1(f)      do 2 l = 2,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)c     ind = -3
_IFN1(f)c     do 3 loopjk = 1,16,4
_IFN1(f)      ind = 13
_IFN1(f)      do 3 loopjk = 5,16,4
_IFN1(f)      ind = ind+4
_IFN1(f)      do ijk = 2,4
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)      enddo
_IFN1(f)   3  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      ind=5
_IF1(f)      do 4 j=1,4
_IF1(f)      call  mvml3(p11,1,gout(ind),4,1,gout(ind),4,1,4)
_IF1(f)  4   ind=ind+16
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,16)
      return
      end
_ENDEXTRACT
_EXTRACT(sp1111,ultra)
c ******************************************************
c ******************************************************
c             =   sp1111  =
c ******************************************************
c ******************************************************
      subroutine sp2222(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF1(cuf)      parameter (dzero=0.0e0)
_IFN1(cuf)      parameter (dzero=0.0d0)
      dimension g(256),h(256)
      dimension kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
INCLUDE(common/miscg)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
INCLUDE(common/pgeom)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/qgeom)
INCLUDE(common/const)
INCLUDE(common/maxc)
INCLUDE(common/astore)
      common/inttab/
     +  aa(333),ba(333),ca(333),abc1,
     +  ab(333),bb(333),cb(333),abc2,
     +  ac(333),bc(333),cc(333),abc3,
     +  ad(333),bd(333),cd(333),abc4,
     +  ae(333),be(333),ce(333),abc5,
     +  af(333),bf(333),cf(333),abc6
INCLUDE(common/auxvar)
c
      dimension gout(*)
c      data dzero/0.0e0/,done/1.0e0/
      data done/1.0d0/
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c ***
c *** this is the fps version of sp1111.
c ***
c *** as much code as possible reduced to loops (>=4)
c *** to avoid ps cache misses and to enhance compiler
c *** optimisation. will probably run like a drain on
c *** the cray-1s]
c ***
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 continue
c     sc = csc(k)
c     sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
c     dq00 = sc*sd*v
c     dq01 = sc*pd*v
c     dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
660   continue
c     p1 = 0.d0
c     p2 = 0.d0
c     p3 = 0.d0
c     p4 = 0.d0
c     p5 = 0.d0
c     p6 = 0.d0
c     q1 = 0.d0
c     q2 = 0.d0
c     q3 = 0.d0
c     q4 = 0.d0
c     q5 = 0.d0
c     q6 = 0.d0
c     r1 = 0.d0
c     r2 = 0.d0
c     r3 = 0.d0
c     r4 = 0.d0
c     r5 = 0.d0
c     r6 = 0.d0
c     r7 = 0.d0
c     r8 = 0.d0
c     r9 = 0.d0
c     v1 = 0.d0
c     v2 = 0.d0
c     v3 = 0.d0
c     v4 = 0.d0
c     v5 = 0.d0
c     v6 = 0.d0
c     w1 = 0.d0
c     w2 = 0.d0
c     w3 = 0.d0
c     w4 = 0.d0
c     w5 = 0.d0
c     w6 = 0.d0
c     w7 = 0.d0
c     w8 = 0.d0
c     w9 = 0.d0
      s1 = 0.d0
      s2 = 0.d0
      s3 = 0.d0
      s4 = 0.d0
      s6 = 0.d0
      s7 = 0.d0
      s8 = 0.d0
      s9 = 0.d0
      s10 = 0.d0
      s11 = 0.d0
      s12 = 0.d0
      s13 = 0.d0
      s14 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      t4 = 0.d0
      t5 = 0.d0
      t6 = 0.d0
      t7 = 0.d0
      t8 = 0.d0
      t9 = 0.d0
      t10 = 0.d0
      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
      t14 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0
      c5 = 0.d0
      c6 = 0.d0
_IF1(x)c$dir scalar
      do 180 ind = 1,ngangb
      isml = ismlq+ismlp(ind)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      eab  = ep(ind)
c     dp00 = dp00p(ind)
c     dp01 = dp01p(ind)
c     dp10 = dp10p(ind)
      ap = app(ind)
      bp = bpp(ind)
      pqab = aqz-ap
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(qperp2+pqab2)
      if (p .le. auxvar) go to 140
      f0 = dsqrt(0.7853981625d0/(p*(gp(ind)+gcd)))*conp(ind)
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      f4 = 3.5d0*f3*gtx
      go to 160
  140 q = conp(ind)/dsqrt(gp(ind)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dfloat(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)+theta4*ca(n+2))*q
      f1 = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)+theta4*cb(n+2))*gy
      f2 = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)+theta4*cc(n+2))*ggy
      f3 = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)+theta4*cd(n+2))*gggy
      f4 = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)+theta4*ce(n+2))*gggy*
     &     gabcd
  160 apbp = ap*bp
      eab2 = eab*eab
c     bpdp01 = bp*dp01
c     apdp10 = ap*dp10
c     edp01 = eab*dp01
c     edp10 = eab*dp10
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f4pqab = f4*pqab
      f1pqa2 = f1*pqab2
      f2pqa2 = f2*pqab2
      f3pqa2 = f3*pqab2
      f4pqa2 = f4*pqab2
      f2pqa3 = f2pqa2*pqab
      f3pqa3 = f3pqa2*pqab
      f4pqa3 = f4pqa2*pqab
c     p1 = p1+f0 *dp00
c     p2 = p2+f1 *dp00
c     p3 = p3+f2 *dp00
c     p4 = p4+f1pqab*dp00
c     p5 = p5+f2pqab*dp00
c     p6 = p6+f2pqa2*dp00
c     r1 = r1+f1 *edp01
c     r2 = r2+f2 *edp01
c     r3 = r3+f3 *edp01
c     r4 = r4+f1pqab *edp01
c     r5 = r5+f2pqab *edp01
c     r6 = r6+f3pqab *edp01
c     r7 = r7+f2pqa2 *edp01
c     r8 = r8+f3pqa2 *edp01
c     r9 = r9+f3pqa3 *edp01
c     w1 = w1+f1 *edp10
c     w2 = w2+f2 *edp10
c     w3 = w3+f3 *edp10
c     w4 = w4+f1pqab *edp10
c     w5 = w5+f2pqab *edp10
c     w6 = w6+f3pqab *edp10
c     w7 = w7+f2pqa2 *edp10
c     w8 = w8+f3pqa2 *edp10
c     w9 = w9+f3pqa3 *edp10
      s1 = s1+f0 *eab
      s2 = s2+f1 *eab
      s3 = s3+f2 *eab
      s4 = s4+f3 *eab
      s6 = s6+f1pqab*eab
      s7 = s7+f2pqab*eab
      s8 = s8+f3pqab*eab
      s9 = s9+f1pqa2*eab
      s10 = s10+f2pqa2*eab
      s11 = s11+f3pqa2*eab
      s12 = s12+f2pqa3*eab
      s13 = s13+f3pqa3*eab
      s14 = s14+f3pqa3*pqab*eab
      t1 = t1+f0 *eab2
      t2 = t2+f1 *eab2
      t3 = t3+f2 *eab2
      t4 = t4+f3 *eab2
      t5 = t5+f4 *eab2
      t6 = t6+f2pqab*eab2
      t7 = t7+f3pqab*eab2
      t8 = t8+f4pqab*eab2
      t9 = t9+f2pqa2*eab2
      t10 = t10+f3pqa2*eab2
      t11 = t11+f4pqa2*eab2
      t12 = t12+f3pqa3*eab2
      t13 = t13+f4pqa3*eab2
      t14 = t14+f4pqa3*pqab*eab2
      if (rabsq .eq. 0.0d0) go to 180
c     q1 = q1+f0 *bpdp01
c     q2 = q2+f1 *bpdp01
c     q3 = q3+f2 *bpdp01
c     q4 = q4+f1pqab*bpdp01
c     q5 = q5+f2pqab*bpdp01
c     q6 = q6+f2pqa2*bpdp01
c     v1 = v1+f0 *apdp10
c     v2 = v2+f1 *apdp10
c     v3 = v3+f2 *apdp10
c     v4 = v4+f1pqab*apdp10
c     v5 = v5+f2pqab*apdp10
c     v6 = v6+f2pqa2*apdp10
      c1 = c1+f0 *apbp
      c2 = c2+f1 *apbp
      c3 = c3+f2 *apbp
      c4 = c4+f1pqab*apbp
      c5 = c5+f2pqab*apbp
      c6 = c6+f2pqa2*apbp
  180 continue
      a1 = aqz*s2-s6
      a2 = aqz*s3-s7
      a3 = aqz*s4-s8
      a4 = aqz*s6-s9
      a5 = aqz*s7-s10
      a6 = aqz*s8-s11
      a8 = aqz*s10-s12
      a9 = aqz*s11-s13
      a10 = aqz*s13-s14
      bqz = aqz-rab
      b1 = bqz*s2-s6
      b2 = bqz*s3-s7
      b3 = bqz*s4-s8
      b4 = bqz*s6-s9
      b5 = bqz*s7-s10
      b6 = bqz*s8-s11
      b8 = bqz*s10-s12
      b9 = bqz*s11-s13
      b10 = bqz*s13-s14
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      hecd2 = 0.5d0*ecd2
      qecd = qperp*ecd
      hqecd = 0.5d0*qecd
      qecd2 = qperp*ecd2
      hqecd2 = 0.5d0*qecd2
      q2ecd = qperp2*ecd
      q3ecd = qperp*q2ecd
      q2ecd2 = qperp2*ecd2
      q3ecd2 = q2ecd2*qperp
c     h(  1) = p1
c     h(  2) = qecd*p2
c     h(  4) = -ecd*p4
c     h( 11) = hecd*(p1-ecd*p2)
c     h(  6) = h( 11)+q2ecd2*p3
c     h(  8) = -qecd2*p5
c     h( 16) = h( 11)+ecd2*p6
c     h( 17) = -qperp*r1
c     h( 49) = r4+q1
c     h( 35) = hecd*r1
c     h( 18) = h( 35)-q2ecd*r2
c     h( 20) = qecd*r5
c     h( 50) = h( 20)+qecd*q2
c     h( 52) = h( 35)-ecd*r7-ecd*q4
c     h( 39) = hqecd2*r2
c     h( 44) = -hecd2*r5
c     h( 27) = h( 39)-qperp*h( 35)
c     h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
c     h( 24) = h( 44)+q2ecd2*r6
c     h( 56) = h( 39)-qecd2*(r8+q5)
c     h( 22) = h( 27)+h( 39)+h( 39)-q3ecd2*r3
c     h( 32) = h( 27)-qecd2*r8
c     h( 54) = h( 59)+q2ecd2*(r6+q3)
c     h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
c     h( 65) = -qperp*w1
c     h(193) = w4+v1
c     h(131) = hecd*w1
c     h( 66) = h(131)-q2ecd*w2
c     h( 68) = qecd*w5
c     h(194) = h( 68)+qecd*v2
c     h(196) = h(131)-ecd*w7-ecd*v4
c     h(135) = hqecd2*w2
c     h(140) = -hecd2*w5
c     h( 75) = h(135)-qperp*h(131)
c     h(203) = h(140)+hecd*(h(193)-ecd*v2)
c     h( 72) = h(140)+q2ecd2*w6
c     h(200) = h(135)-qecd2*(w8+v5)
c     h( 70) = h( 75)+h(135)+h(135)-q3ecd2*w3
c     h( 80) = h( 75)-qecd2*w8
c     h(198) = h(203)+q2ecd2*(w6+v3)
c     h(208) = h(203)+h(140)+h(140)+ecd2*(w9+v6)
      h(161) = 0.5d0*(s1-t2)
      h( 81) = h(161)+qperp2*t3
      h(113) = -qperp*(t6+b1)
      h(209) = -qperp*(t6+a1)
      h(241) = h(161)+t9+a4+b4+c1
      h(162) = hqecd*(s2-t3)
      h( 82) = h(162)-qecd*t3+q3ecd*t4
      temp = hecd*t6-q2ecd*t7
      h(114) = temp+hecd*b1-q2ecd*b2
      h(210) = temp+hecd*a1-q2ecd*a2
      h(242) = h(162)+qecd*(t10+a5+b5+c2)
      h( 99) = -hqecd*t3
      h(147) = h( 99)
      h(179) = hecd*(t6+b1)
      h(227) = hecd*(t6+a1)
      h(164) = hecd*(t6-s6)
      h( 84) = h(164)-q2ecd*t7
      temp = -hqecd*t3+qecd*t10
      h(116) = temp+qecd*b5
      h(212) = temp+qecd*a5
      h(244) = h(164)+ecd*(t6-t12-a8-b8-c4)+hecd*(a1+b1)
      h(103) = 0.25d0*ecd2*t3-0.5d0*q2ecd2*t4
      h(151) = h(103)
      h(183) = hqecd2*(t7+b2)
      h(231) = hqecd2*(t7+a2)
      h(108) = hqecd2*t7
      h(156) = h(108)
      h(188) = hecd2*(0.5d0*t3-t10-b5)
      h(236) = hecd2*(0.5d0*t3-t10-a5)
      hxxyy = 0.25d0*(ecd*(s1-t2)-ecd2*(s2-t3))
      h(171) = hxxyy+hecd2*t3
      h( 91) = hxxyy+0.5d0*(q2ecd*t3-q2ecd2*t4)
      temp = hqecd*(ecd*t7-t6)
      h(123) = temp+hqecd*(ecd*b2-b1)
      h(219) = temp+hqecd*(ecd*a2-a1)
      h(251) = hxxyy+hecd*(t9+a4+b4+c1)-hecd2*(t10+a5+b5+c2)
      h(166) = hxxyy+0.5d0*q2ecd2*(s3-t4)
      h( 86) = hxxyy+(hecd2+0.5d0*q2ecd)*t3+q2ecd2*(-3.d0*t4+
     +    0.5d0*s3+qperp2*t5)
      h(118) = 1.5d0*qecd2*(t7+b2)-hqecd*(t6+b1)-q3ecd2*(b3+t8)
      h(214) = 1.5d0*qecd2*(t7+a2)-hqecd*(t6+a1)-q3ecd2*(a3+t8)
      h(246) = hxxyy-hecd2*(qperp2*t4+t10+a5+b5)+hecd*(t9+a4+b4+c1-ecd*
     +    c2)+q2ecd2*(t11+0.5d0*s3+a6+b6+c3)
      h(168) = hqecd2*(t7-s7)
      h( 88) = 1.5d0*qecd2*t7-hqecd2*s7-q3ecd2*t8
      temp = hecd2*(0.5d0*t3-t10)+q2ecd2*(t11-0.5d0*t4)
      h(120) = temp-hecd2*b5+q2ecd2*b6
      h(216) = temp-hecd2*a5+q2ecd2*a6
      h(248) = qecd2*(1.5d0*t7-t13-a9-b9-c5)-hqecd2*(s7-a2-b2)
      h(176) = hxxyy+hecd2*(s10-t10)
      h( 96) = hxxyy-hecd2*(qperp2*t4+t10-s10)+0.5d0*q2ecd*t3+q2ecd2*
     +     t11
      h(128) = qecd2*(1.5d0*t7-t13-b9)-hqecd*(t6+b1)+hqecd2*b2
      h(224) = qecd2*(1.5d0*t7-t13-a9)-hqecd*(t6+a1)+hqecd2*a2
      h(256) = hxxyy+hecd2*(-3.d0*(a5+b5)+t3+s10-c2)+ecd2*(-3.d0*t10+
     +     t14+a10+b10+c6)+hecd*(t9+a4+b4+c1)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
_IF1(a)cvd$  shortloop
c     do 103 kq1=22,214,48
      do 103 kq1=118,214,48
          g(kq1  ) = v44*h(kq1) + v47*h(kq1+5)
          g(kq1+1) = v54*h(kq1) + v57*h(kq1+5)
103       g(kq1+5) = v74*h(kq1) + v77*h(kq1+5)
_IF1(a)cvd$  shortloop
c     do 101 kq1=24,216,48
      do 101 kq1=120,216,48
          g(kq1  ) = cosp*h(kq1)
101       g(kq1+4) = sinp*h(kq1)
_IF1(a)cvd$  shortloop
c     do 102 kq1=18,210,48
      do 102 kq1=114,210,48
          g(kq1  ) = cosp*h(kq1)
102       g(kq1+1) = sinp*h(kq1)
c     g( 80) = h( 80)
      g( 86) = v44*h( 86)+v47*h( 91)
      g( 87) = v54*h( 86)+v57*h( 91)
      g( 91) = v74*h( 86)+v77*h( 91)
      g( 88) = cosp*h( 88)
      g( 92) = sinp*h( 88)
      g( 96) = h( 96)
      g(102) = v45*h(103)
      g(103) = v55*h(103)
      g(107) = -g(102)
      g(104) = u12*h(108)
      g(108) = cosp*h(108)
      g(112) = dzero
      g(128) = h(128)
c     g(134) = v45*h(135)
c     g(135) = v55*h(135)
c     g(139) = -g(134)
c     g(136) = u12*h(140)
c     g(140) = cosp*h(140)
      g(144) = dzero
      g(150) = v45*h(151)
      g(151) = v55*h(151)
      g(155) = -g(150)
      g(152) = u12*h(156)
      g(156) = cosp*h(156)
      g(160) = dzero
      g(176) = h(176)
      g(182) = v45*h(183)
      g(183) = v55*h(183)
      g(187) = -g(182)
      g(184) = u12*h(188)
      g(188) = cosp*h(188)
      g(192) = dzero
c     g(198) = v44*h(198)+v47*h(203)
c     g(199) = v54*h(198)+v57*h(203)
c     g(203) = v74*h(198)+v77*h(203)
c     g(200) = cosp*h(200)
c     g(204) = sinp*h(200)
      g(230) = v45*h(231)
      g(231) = v55*h(231)
      g(235) = -g(230)
      g(232) = u12*h(236)
      g(236) = cosp*h(236)
      g(240) = dzero
      g(246) = v44*h(246)+v47*h(251)
      g(247) = v54*h(246)+v57*h(251)
      g(251) = v74*h(246)+v77*h(251)
      g(248) = cosp*h(248)
      g(252) = sinp*h(248)
c     g( 38) = v45*h( 39)
c     g( 39) = v55*h( 39)
c     g( 43) = -g( 38)
c     g( 40) = u12*h( 44)
c     g( 44) = cosp*h( 44)
c2    g( 48) = dzero
c     g( 54) = v44*h( 54)+v47*h( 59)
c     g( 55) = v54*h( 54)+v57*h( 59)
c     g( 59) = v74*h( 54)+v77*h( 59)
c     g( 56) = cosp*h( 56)
c     g( 60) = sinp*h( 56)
c     g(  6) = v44*h(  6)+v47*h( 11)
c     g(  7) = v54*h(  6)+v57*h( 11)
c     g( 11) = v74*h(  6)+v77*h( 11)
c     g(  8) = cosp*h(  8)
c     g( 12) = sinp*h(  8)
c     g( 68) = h( 68)
      g( 82) = cosp*h( 82)
      g( 83) = sinp*h( 82)
      g( 84) = h( 84)
      g( 98) = u12*h( 99)
      g( 99) = cosp*h( 99)
      g(100) = dzero
      g(116) = h(116)
c     g(130) = u12*h(131)
c     g(131) = cosp*h(131)
      g(132) = dzero
      g(146) = u12*h(147)
      g(147) = cosp*h(147)
      g(148) = dzero
      g(164) = h(164)
      g(178) = u12*h(179)
      g(179) = cosp*h(179)
      g(180) = dzero
c     g(194) = cosp*h(194)
c     g(195) = sinp*h(194)
      g(226) = u12*h(227)
      g(227) = cosp*h(227)
      g(228) = dzero
      g(242) = cosp*h(242)
      g(243) = sinp*h(242)
c     g( 34) = u12*h( 35)
c     g( 35) = cosp*h( 35)
c2    g( 36) = dzero
c     g( 50) = cosp*h( 50)
c     g( 51) = sinp*h( 50)
c     g(  2) = cosp*h(  2)
c     g(  3) = sinp*h(  2)
c     g( 65) = h( 65)
      g( 81) = h( 81)
      g( 97) = dzero
      g(113) = h(113)
      g(129) = dzero
      g(145) = dzero
      g(161) = h(161)
      g(177) = dzero
      g(225) = dzero
c2    g( 33) = dzero
c     h( 80) = cosp*g( 80)
      h( 96) = cosp*g( 96)
      h(112) =           u12*g(176)
      h(128) = cosp*g(128)
c     h(144) = sinp*g( 80)
      h(160) = sinp*g( 96)
      h(176) =           cosp*g(176)
      h(192) = sinp*g(128)
_IF1(a)cvd$  shortloop
c     do 121 kq1=70,118,16
      do 121 kq1=86,118,16
          h(kq1   ) = cosp*g(kq1   ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1   ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+ 1) + u12*g(kq1+65)
          h(kq1+65) = sinp*g(kq1+ 1) + cosp*g(kq1+65)
          h(kq1+ 2) = cosp*g(kq1+ 2) + u12*g(kq1+66)
          h(kq1+66) = sinp*g(kq1+ 2) + cosp*g(kq1+66)
          h(kq1+ 5) = cosp*g(kq1+ 5) + u12*g(kq1+69)
          h(kq1+69) = sinp*g(kq1+ 5) + cosp*g(kq1+69)
          h(kq1+ 6) = cosp*g(kq1+ 6) + u12*g(kq1+70)
          h(kq1+70) = sinp*g(kq1+ 6) + cosp*g(kq1+70)
121   continue
c     h( 68) = cosp*g( 68)
      h( 84) = cosp*g( 84)
      h(100) =           u12*g(164)
      h(116) = cosp*g(116)
c     h(132) = sinp*g( 68)
      h(148) = sinp*g( 84)
      h(164) =           cosp*g(164)
      h(180) = sinp*g(116)
_IF1(a)cvd$  shortloop
c     do 122 kq1=66,114,16
      do 122 kq1=82,114,16
          h(kq1   ) = cosp*g(kq1  ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1  ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+1) + u12*g(kq1+65)
122       h(kq1+65) = sinp*g(kq1+1) + cosp*g(kq1+65)
_IF1(a)cvd$  shortloop
c      do 1221 kq1=2,50,16
c          h(kq1   ) = g(kq1   )
c          h(kq1+ 1) = g(kq1+ 1)
c          h(kq1+ 4) = g(kq1+ 4)
c          h(kq1+ 5) = g(kq1+ 5)
c          h(kq1+ 6) = g(kq1+ 6)
c1221      h(kq1+ 9) = g(kq1+ 9)
_IF1(a)cvd$  shortloop
c     do 1222 kq1=12,60,16
      do 1222 kq1=28,60,16
c         h(kq1    ) = g(kq1    )
          h(kq1+182) = g(kq1+182)
          h(kq1+183) = g(kq1+183)
          h(kq1+186) = g(kq1+186)
          h(kq1+187) = g(kq1+187)
1222      h(kq1+188) = g(kq1+188)
_IF1(a)cvd$  shortloop
c     do 1223 kq1=203,251,16
      do 1223 kq1=219,251,16
          h(kq1  ) = g(kq1  )
1223      h(kq1+1) = g(kq1+1)
c     h( 65) = cosp*g( 65)
      h( 81) = cosp*g( 81)
      h( 97) =           u12*g(161)
      h(113) = cosp*g(113)
c     h(129) = sinp*g( 65)
      h(145) = sinp*g( 81)
      h(161) =           cosp*g(161)
      h(177) = sinp*g(113)
      h( 48) = g( 48)
      h( 36) = g( 36)
      h(228) = g(228)
      h(240) = g(240)
      h(225) = g(225)
      h( 33) = g( 33)
_IF1(a)cvd$  shortloop
c     do 123 kq1=22,214,64
      do 123 kq1=86,214,64
          g(kq1   ) = cosp*h(kq1   ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1   ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+ 1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+ 1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+ 2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+ 2) + cosp*h(kq1+18)
          g(kq1+ 5) = cosp*h(kq1+ 5) + u12* h(kq1+21)
          g(kq1+21) = sinp*h(kq1+ 5) + cosp*h(kq1+21)
          g(kq1+ 6) = cosp*h(kq1+ 6) + u12* h(kq1+22)
          g(kq1+22) = sinp*h(kq1+ 6) + cosp*h(kq1+22)
          g(kq1+10) = cosp*h(kq1+10) + u12* h(kq1+26)
123       g(kq1+26) = sinp*h(kq1+10) + cosp* h(kq1+26)
_IF1(a)cvd$  shortloop
c     do 124 kq1=17,209,64
      do 124 kq1=81,209,64
          g(kq1   ) = cosp*h(kq1  ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1  ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+2) + cosp*h(kq1+18)
          g(kq1+ 3) = cosp*h(kq1+3) + u12* h(kq1+19)
124       g(kq1+19) = sinp*h(kq1+3) + cosp* h(kq1+19)
_IF1(a)cvd$  concur
c     do 125 kq1=49,177,64
      do 125 kq1=113,177,64
          kkq1=kq1
_IF1(a)cvd$  shortloop
          do 126 kkkq1=1,32
              g(kkq1)=h(kkq1)
              kkq1=kkq1+1
126       continue
125   continue
_IF1(a)cvd$  shortloop
      do 127 kq1=1,16
c         g(kq1)=h(kq1)
127       g(kq1+240)=h(kq1+240)
      goto 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
c920   do 921 kkq1=1,10
920   do 921 kkq1=4,10
          kq1=kq1off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) = h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) = h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
921       g(kq1+15) = h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
c     do 922 kkq1=1,6
      do 922 kkq1=2,6
          kq1=kq2off(kkq1)
          g(kq1   ) = dzero
          g(kq1+ 1) = dzero
          g(kq1+ 2) = h(kq1+ 2)
          g(kq1+ 3) = dzero
          g(kq1+ 5) = dzero
          g(kq1+ 6) = h(kq1+ 6)
          g(kq1+ 7) = dzero
          g(kq1+10) = dzero
          g(kq1+11) = h(kq1+11)
922       g(kq1+15) = dzero
      go to 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
c1000  do 1001 kkq1=1,6
1000  do 1001 kkq1=3,6
          kq1=kq3off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) =-h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) =-h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
1001      g(kq1+15) = h(kq1+15)
c     do 1002 kkq1=1,4
      do 1002 kkq1=2,4
          kq1=kq4off(kkq1)
          g(kq1   ) = -h(kq1   )
          g(kq1+ 1) =  h(kq1+ 1)
          g(kq1+ 2) =  dzero
          g(kq1+ 3) = -h(kq1+ 3)
          g(kq1+ 5) = -h(kq1+ 5)
          g(kq1+ 6) =  dzero
          g(kq1+ 7) =  h(kq1+ 7)
          g(kq1+10) = -h(kq1+10)
          g(kq1+11) =  dzero
1002      g(kq1+15) = -h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
c     do 1003 kkq1=1,6
      do 1003 kkq1=2,6
          kq1=kq5off(kkq1)
          g(kq1   ) =  dzero
          g(kq1+ 1) =  dzero
          g(kq1+ 2) =  h(kq1+ 2)
          g(kq1+ 3) =  dzero
          g(kq1+ 5) =  dzero
          g(kq1+ 6) = -h(kq1+ 6)
          g(kq1+ 7) =  dzero
          g(kq1+10) =  dzero
          g(kq1+11) =  h(kq1+11)
1003      g(kq1+15) =  dzero
          g(99)=-g(99)
          g(108)=-g(108)
          g(147)=-g(147)
          g(156)=-g(156)
          g(103)=-g(103)
          g(151)=-g(151)
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
_IF1(a)cvd$  shortloop
c     do 2001 kq1=2,242,16
      do 2001 kq1=66,242,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+12) = g(kq1+ 6)
          g(kq1+13) = g(kq1+10)
2001  continue
      if (rcdsq) 1200,1200,1300
_IF1(a)cvd$  concur
1300  do 1301 kq1=1,4
c         kkq1=kq1
          kkq1=kq1+5*16
_IF1(a)cvd$  shortloop
c         do 1302 jq1=1,16
          do 1302 jq1=6,16
              g(kkq1+4) = r13*g(kkq1) + g(kkq1+4)
              g(kkq1+12)= r33*g(kkq1) + g(kkq1+12)
1302          kkq1=kkq1+16
1301  continue
c ***
c      do 1303 kq1=1,253,4
c          g(kq1+1) = r14*g(kq1) + g(kq1+1)
c1303      g(kq1+3) = r34*g(kq1) + g(kq1+3)
c ***
      do 1303 ii=86,246,64
         kq2=ii
         do jj=1,3
            kq1=kq2
            do kk=1,3
               g(kq1  )    = r14 *g(kq1-1) + g(kq1  )
               g(kq1+2)    = r34 *g(kq1-1) + g(kq1+2)
               kq1=kq1+4
            enddo 
            kq2=kq2+16
         enddo
1303  continue
1200  continue
c     call daxpy(256,dq11,g,1,gout,1)
      do ii=86,246,64
         kq2=ii
         do jj=1,3
            kq1=kq2
            do kk=1,3
               gout(kq1  ) = dq11*g(kq1  ) + gout(kq1  )
               gout(kq1+1) = dq11*g(kq1+1) + gout(kq1+1)
               gout(kq1+2) = dq11*g(kq1+2) + gout(kq1+2)
               kq1=kq1+4
            enddo 
            kq2=kq2+16
         enddo
      enddo
c     dq01x=dq01-dq11
c     dq01x=-dq11
_IF1(a)cvd$  vector
c      do 1202 kq1=2,242,16
c          gout(kq1  ) = dq01x*g(kq1  ) + gout(kq1  )
c          gout(kq1+1) = dq01x*g(kq1+1) + gout(kq1+1)
c1202      gout(kq1+2) = dq01x*g(kq1+2) + gout(kq1+2)
c     dq10x=dq10-dq11
c     dq00x=dq00-dq11
c     dq10x=-dq11
c     dq00x=-dq11
_IF1(a)cvd$  vector
c      do 1203 kq1=1,241,16
c          gout(kq1   ) = dq00x*g(kq1   ) + gout(kq1   )
c          gout(kq1+ 4) = dq10x*g(kq1+ 4) + gout(kq1+ 4)
c          gout(kq1+ 8) = dq10x*g(kq1+ 8) + gout(kq1+ 8)
c1203      gout(kq1+12) = dq10x*g(kq1+12) + gout(kq1+12)
940   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 1111
c
c
_IFN1(f)c     i1 = 64
_IFN1(f)c     i2 = 128
_IFN1(f)c     i3 = 192
_IFN1(f)      i1 = 84
_IFN1(f)      i2 = 148
_IFN1(f)      i3 = 212
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)c     do 1 jkl = 1,64,4
_IFN1(f)      do 1 jkl = 21,64,16
_IFN1(f)      do ii=1,3
_IFN1(f)      i1 = i1+1
_IFN1(f)      i2 = i2+1
_IFN1(f)      i3 = i3+1
_IFN1(f)      do njkl = 2,4
_IFN1(f)      i1 = i1+1
_IFN1(f)      i2 = i2+1
_IFN1(f)      i3 = i3+1
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)      enddo
_IFN1(f)      enddo
_IFN1(f)      i1 = i1+4
_IFN1(f)      i2 = i2+4
_IFN1(f)      i3 = i3+4
_IFN1(f)    1 continue
_IFN1(f)c     ind = -48
_IFN1(f)c     do 2 i = 1,4
_IFN1(f)      ind = 16
_IFN1(f)      do 2 i = 2,4
_IFN1(f)      ind = ind+48
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)c     do 2 loopkl = 1,16,4
_IFN1(f)      ind = ind+4
_IFN1(f)      do 2 loopkl = 5,16,4
_IFN1(f)      ind = ind+1
_IFN1(f)      do nkl = 2,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)      enddo
_IFN1(f)    2 continue
_IFN1(f)c     ind = -12
_IFN1(f)c     do 3 loopij = 1,16
_IFN1(f)c     ind = ind+13
_IFN1(f)      ind = 52
_IFN1(f)      do 3 loopij = 2,4
_IFN1(f)      ind = ind+16
_IFN1(f)      do 3 ij2 = 2,4
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      ind = ind+1
_IFN1(f)      do 3 l = 2,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    3 continue
_IFN1(f)c     ind = -3
_IFN1(f)      ind = 61
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)c     do 4 ijk = 1,64,4
_IFN1(f)      do 4 ijk = 21,64,16
_IFN1(f)      ind = ind+16
_IFN1(f)      do 4 ijk2 = 2,4
_IFN1(f)      ind = ind+4
_IFN1(f)      do nijk = 2,4
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)      enddo
_IFN1(f)    4 continue
_IF1(f)      call mvml3(p11,1,gout(65),64,1,gout(65),64,1,64)
_IF1(f)      ind=17
_IF1(f)      do 5 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),16,1,gout(ind),16,1,16)
_IF1(f)   5  ind=ind+64
_IF1(f)      ind=5
_IF1(f)      do 6 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),4,16,gout(ind),4,16,16)
_IF1(f)   6  ind=ind+1
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,64)
c
      return
      end
_ENDEXTRACT
_ENDIF
_ELSE
      subroutine jkin7a(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels)
_IF1(a)cvd$r noconcur
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
c ***
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
INCLUDE(common/ijlab)
INCLUDE(common/timez)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
INCLUDE(common/pkfil)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
INCLUDE(common/morokuma)
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
c     data m22,m25,mword1/22,25,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c ***
c *** ** preliminary ** vectorised integrals code.
c *** should provide full accuracy integrals (12-13 sig. fig.)
c *** only (00/11),(01/11) and (11/11) are currently vectorised.
c ***
_IF(parallel)
c***   **MPP**
      iflop = iipsci()
c***   **MPP**
_ENDIF
c
_IF(ccpdft)
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
_ENDIF
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      if(.not.odscf) then
      else
        l2 = ikyp(numorb)
      endif
c
      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
      call filmax
      if(kst0.le.0. or. lst0.le.0) go to 4010
      do 400 ii = ist0,nshels
      if(kad(ii))400,401,401
 401  isht=ktype(ii)*10
      dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
      endif
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
      j0 = jst0
      do 380 jj = j0,ii
      jst0 = 1
      if(kad(jj))380,141,141
c ***
  141  ijsht=isht + ktype(jj)
          if(ijsht.eq.22) goto 380
c ***
          itrij=ikyii+jj
          mij = itrij
          tolij=dlntol+prefac(itrij)
_IF1()c
_IF1()c **** rejecting on basis of ij arguably too extreme
_IF1()c ***  suppress
_IF1()c         if(tolij .le. -3.401e0) then
_IF1()c           if(odscf)intcut(1)=intcut(1)+1
_IF1()c           goto 380
_IF1()c         endif
c ***
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
      ishell = ii
      jshell = jj
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
  361 do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ksht = ktype(kk)*10
      ikykk=iky(kk)
      if (odscf) then
        itrik=ikyii+kk
        mik = itrik
        itrjk=iky(max(jj,kk))+min(jj,kk)
        mjk = itrjk
        tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
_IF(parallel)
c***   **MPP**
      if (oipsci()) go to 340
c***   **MPP**
_ENDIF
c ***
      klsht = ksht + ktype(ll)
      if(klsht .eq. 22.or. kad(ll).lt.0) goto 340
c ***
      itrkl=ikykk+ll
      tijkl=tolij+prefac(itrkl)
      if(tijkl.le.0.0d0) then
       if(odscf)intcut(2)=intcut(2)+1
       goto 340
      endif
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=tijkl+
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            intcut(3)=intcut(3)+1
            goto 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dfloat(nt)/ dfloat(n4)
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
      call genr70(gout,1,.false.)
      if (odscf) then
        if (zscftp.eq.'uhf')then
_IF(ccpdft)
        call dir_build_uhf70(fock,fockb,
     +  dens,densb,gout,
     +  fac1,fac2, facex, ocoul, oexch)
_ELSE
        call dir_build_uhf70(fock,fockb,
     +                dens,densb,gout)
_ENDIF
        else if (zscftp.eq.'gvb'.or.zscftp.eq.'grhf') then
         if(nsheld.le.1) then
          call dir_build_open_70(fock,exch,dens,gout)
         else
          call dir_build_open2_70(l2,fock,exch,dens,gout)
         endif
        else
_IF(ccpdft)
_IF(cray)
        call qoutd70(fock,dens,gout,
     +        fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
       if(omorok) then
        call dbuild70_morok(fock,dens,gout)
        else
        call dbuild70(fock,dens,gout,
     +        fac1, fac2, facex, ocoul, oexch)
        endif
_ENDIF
_ELSE
_IF(cray)
        call qoutd70(fock,dens,gout)
_ELSE
        call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
        endif
      else
        call qout70(gout)
      endif
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim) then
          go to 420
      endif
  340 continue
  360 continue
  380 continue
      time = cpulft(1)
  400 continue
c ***
      ist = 1
      jst = 1
      kst = 1
      lst = 1
_IF(ccpdft)
4010  call jk70va(zscftp,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +            iso,gout,nshels,fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
4010  call jk70va(zscftp,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +            iso,gout,nshels)
_ENDIF
c
      if(omaxb.or.tim.gt.timlim)go to 420
c ***
      call final(q,fock,dens)
  420 continue

      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
      end
      subroutine jkin7s(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels,outvv)
_IF1(a)cvd$r noconcur
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
c ***
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/symtry)
INCLUDE(common/shlg70)
INCLUDE(common/picon)
INCLUDE(common/ijlab)
INCLUDE(common/timez)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
INCLUDE(common/pkfil)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
INCLUDE(common/morokuma)
c
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      data m25/25/
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c ***
c *** ** preliminary ** vectorised integrals code.
c *** should provide full accuracy integrals (12-13 sig. fig.)
c *** only (00/11),(01/11) and (11/11) are currently vectorised.
c ***
_IF(parallel)
c***   **MPP**
      iflop = iipsci()
c***   **MPP**
_ENDIF
_IF(ccpdft)
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
_ENDIF
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      ltri = ikyp(nshels)
      oskipp = .false.

      dlncutoff = dlog(cutoff)

      ischw = igmem_alloc(ltri)
      nschwz = 0
      if(.not.odscf) then
      else
         l2 = ikyp(numorb)
      endif

      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c *** read in ints for schwarz inequality test
      call secget(isect(421),m25,iblk25)
      call rdedx(q(ischw),ltri,iblk25,idaf)
      dlnmxs = -1.0d50
      do ii = 0, ltri-1
         dlnmxs = dmax1(dlnmxs,q(ischw+ii))
      enddo
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
      call filmax
      if(kst0.le.0. or. lst0.le.0) go to 4010
      do 400 ii = ist0,nshels
      if(kad(ii))400,401,401
 401  isht=ktype(ii)*10
      dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
      endif
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
      j0 = jst0
      do 380 jj = j0,ii
      jst0 = 1
      if(kad(jj))380,141,141
c ***
  141  ijsht=isht + ktype(jj)
          if(ijsht.eq.22) goto 380
c ***
          itrij=ikyii+jj
          mij = itrij
_IF1()c         tolij=dlntol+prefac(itrij)
_IF1()c
_IF1()c **** rejecting on basis of ij arguably too extreme
_IF1()c ***  suppress
_IF1()c         if(tolij .le. -3.401e0) then
_IF1()c           if(odscf)intcut(1)=intcut(1)+1
_IF1()c           goto 380
_IF1()c         endif
c ***
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
      ishell = ii
      jshell = jj
c *** reject on basis of ij 
      call scrkl(q(ischw),iso,tolij,nshels)
      ijij = itrij + ischw -1
      test = q(ijij) + tolij
      if(test.lt.dlncutoff) then
       if(odscf)intcut(1)=intcut(1)+1
       goto 380
      endif
c
c     ----- kshell -----
c
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
  361 do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ksht = ktype(kk)*10
      ikykk=iky(kk)
      if (odscf) then
         itrik=ikyii+kk
         mik = itrik
         itrjk=iky(max(jj,kk))+min(jj,kk)
         mjk = itrjk
         tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
_IF(parallel)
c***   **MPP**
      if (oipsci()) go to 340
c***   **MPP**
_ENDIF
c ***
      klsht = ksht + ktype(ll)
      if(klsht .eq. 22.or. kad(ll).lt.0) goto 340
c ***
      itrkl=ikykk+ll
      ijij = itrij + ischw -1
      klkl = itrkl + ischw -1
      test = q(ijij) + q(klkl)
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=dlntol+ test +
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            nschwz = nschwz + 1
            intcut(3)=intcut(3)+1
            goto 340
         endif
       else
         oskipp = test.lt.dlncutoff
         if(oskipp) then
           nschwz = nschwz + 1
           go to 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dfloat(nt)/ dfloat(n4)
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
      call genr70(gout,1,.false.)
      if (odscf) then
        if (zscftp.eq.'uhf')then
_IF(ccpdft)
        call dir_build_uhf70(fock,fockb,
     +  dens,densb,gout,
     +  fac1,fac2, facex, ocoul, oexch)
_ELSE
        call dir_build_uhf70(fock,fockb,
     +                dens,densb,gout)
_ENDIF
        else if (zscftp.eq.'gvb'.or.zscftp.eq.'grhf') then
         if(nsheld.le.1) then
          call dir_build_open_70(fock,exch,dens,gout)
         else
          call dir_build_open2_70(l2,fock,exch,dens,gout)
         endif
        else
_IF(ccpdft)
_IF(cray)
        call qoutd70(fock,dens,gout,
     +        fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
       if(omorok) then
        call dbuild70_morok(fock,dens,gout)
        else
        call dbuild70(fock,dens,gout,
     +        fac1, fac2, facex, ocoul, oexch)
        endif
_ENDIF
_ELSE
_IF(cray)
        call qoutd70(fock,dens,gout)
_ELSE
        call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
        endif
      else
        call qout70(gout)
      endif
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim)then
         go to 420
      endif
  340 continue
  360 continue
  380 continue
      time = cpulft(1)
  400 continue
c ***
      ist = 1
      jst = 1
      kst = 1
      lst = 1
_IF(ccpdft)
4010  call jk70vs(zscftp,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +     q(ischw),iso,gout,nshels,nschwz,
     +     fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
4010  call jk70vs(zscftp,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +     q(ischw),iso,gout,nshels,nschwz)
_ENDIF
      if(omaxb.or.tim.gt.timlim)go to 420
c ***
      call final(q,fock,dens)
  420 if(outvv) write(iwr,6030) nschwz

      call gmem_free(ischw)
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
 6030 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine filmax
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/auxvar)
INCLUDE(common/nshel)
INCLUDE(common/maxc)
INCLUDE(common/shlt)
      data five/5.0d0/,fiften/15.0d0/
      data ten,twenty/10.0d0,20.0d0/
      do 100 i = 1,nshell
      l = kstart(i)
      n = l+kng(i)-1
      do 100 j = l,n
      a1 =  dabs(cs(j))
      a2 =  dabs(cp(j))
      cmax(j) = dmax1(a1,a2)
  100 continue
c
      if (ifasti.eq.1) then
c     settings in gamess-us
       var(1) = twenty
       var(2) = ten
       error1 = dmin1(cutoff,1.0d-10)
       error2 = dmin1(cutoff,1.0d-10)
      else if (ifasti.eq.2) then
c     tighten thresholds
       var(1) = fiften
       var(2) = five
       error1 = dsqrt(cutoff) * 0.1d0
       error2 = cutoff * 0.01d0
      else
_IF1()c      error1 = sqrt(cutoff)
_IF1()c      error2 = cutoff
       error1 = dsqrt(cutoff) * 0.1d0
       var(1) = fiften
       var(2) = five
       error2 = cutoff * 0.01d0
      endif
c ***
c *** set up /values/ for accurate evaluation of incomplete gamma
c ***
      call setfm
c ***
      return
      end
      subroutine sp0000(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc(5001)
      common/junk/p(mxprms2),f0(mxprms2),ind(mxprms2)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/ginf)
INCLUDE(common/pqgeom)
INCLUDE(common/pgeom)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
      dimension q(mxprms2),gout(*)
      g0000 = 0.d0
      do 220 k = 1,ngc
      gc = cgg(k)
      csck = csc(k)
      gcrcds = gc*rcdsq
      do 220 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.d0/gcd
      gdecd = gd*ecd
      pp = ecd* dexp(-gdecd*gcrcds)
      qqtest = pp*cmaxc(k)*cmaxd(l)
      if (qqtest .le. error1) go to 100
      ismlq = 0
      go to 120
  100 if (qqtest .le. error2) go to 220
      ismlq = 1
  120 cq = gdecd*rcd
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      h0000 = 0.d0
_IF1(a)c scalar concurrent wins for all ngangb !!
_IF1(a)cvd$  novector
      do 201 i=1,ngangb
         p(i) = ((aqz-app(i))**2+qperp2)/(ep(i)+ecd)
         q(i) = dp00p(i)/dsqrt(gp(i)+gcd)
 201  continue
      call compfs(ngangb,ngangb,p,f0,ind,0)
_IF1(a)cvd$  novector
      do 202 i=1,ngangb
202       h0000 = h0000 + f0(i)*q(i)
c ***
      g0000 = g0000+h0000*csck*csd(l)*pp
  220 continue
      gout(1) = g0000
      return
      end
      subroutine sp0001(gout)
c        *****  special fast routine for -p- loop for 0001 ****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c
INCLUDE(common/pqgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(334),a1(333),b1(333),c1(334),
     +  a2(4000)
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
c
c mxprms2 = square of max no. of primitives in a contraction
      dimension p(mxprms2),fs(mxprms2,2),g(mxprms2)
      dimension pqab(mxprms2),q(mxprms2),iind(mxprms2)
c
      dimension gout(*)
c
      data sixty,tenm12/60.0d0,1.0d-12/
c
      gout1 = 0.0d0
      gout2 = 0.0d0
      gout3 = 0.0d0
      gout4 = 0.0d0
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.0d0/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = 1.0d0
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
c modifications to bring up to full accuracy
_IF1(a)c scalar concurrent preferable to vector
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 178 i = 1,ngangb
         pqab(i) = aqz-app(i)
         g(i) = 1.d0/(ep(i)+ecd)
         p(i) = (pqab(i)*pqab(i)+qperp2)*g(i)
         q(i) = dp00p(i) / dsqrt(gp(i)+gcd)
 178  continue
      call compfs(ngangb,mxprms2,p,fs,iind,1)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 179 i = 1,ngangb
         q0 = q(i) * fs(i,1)
         q1 = q(i) * fs(i,2)
         u = g(i)*q1
         h0000 = h0000+q0
         h0001 = h0001+u
         h0003 = h0003-u*pqab(i)
 179  continue
_IF1()c      do 180 i = 1,ngangb
_IF1()c      isml = ismlq+ismlp(i)
_IF1()c      if (isml .ge. 2) go to 180
_IF1()c      auxvar = var(isml+1)
_IF1()c      pqab = aqz-app(i)
_IF1()c      g = 1.e0/(ep(i)+ecd)
_IF1()c      p = (pqab*pqab+qperp2)*g
_IF1()c      if (p .le. auxvar) go to 140
_IF1()c      q0 = dp00p(i)*sqrt (0.7853981625e0/(p*(gp(i)+gcd)))
_IF1()c      q1 = 0.5e0*q0/p
_IF1()c      go to 160
_IF1()c  140 q = dp00p(i)/sqrt (gp(i)+gcd)
_IF1()c      qq = p*12.5e0
_IF1()c      n =  int (qq)
_IF1()c      theta = qq- float(n)
_IF1()c      theta2 = theta*(theta-1.e0)
_IF1()c      theta3 = theta2*(theta-2.e0)
_IF1()c      theta4 = theta2*(theta+1.e0)
_IF1()c      q0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
_IF1()c      q1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*q
_IF1()c  160 continue
_IF1()c      u = g*q1
_IF1()c      h0000 = h0000+q0
_IF1()c      h0001 = h0001+u
_IF1()c      h0003 = h0003-u*pqab
_IF1()c  180 continue
      h0001 = h0001*ecd*qperp
      h0003 = h0003*ecd
      pp = dq*h0000
      g0001 = h0001*cosp+pp*sing
      g0002 = h0001*sinp
      g0003 = h0003+pp*cosg
      gout1 = gout1+dq00*h0000
      gout2 = gout2+dq01*g0001
      gout3 = gout3+dq01*g0002
      gout4 = gout4+dq01*g0003
 940  continue
      t1 = gout2
      t2 = gout3
      t3 = gout4
      gout(1) = gout1
      gout(2) = p11*t1+p21*t2+p31*t3
      gout(3) = p12*t1+p22*t2+p32*t3
      gout(4) = p13*t1+p23*t2+p33*t3
      return
      end
      subroutine sp0011(gout)
c        *****  special fast routine for -p- loop for 0011 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
INCLUDE(common/maxc)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
c
      dimension g(mxprms2),p(mxprms2),pqab(mxprms2)
      dimension q(mxprms2),iind(mxprms2),fs(mxprms2,3)
      dimension gout(*)
      equivalence (v12,u12)
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
 660  h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0011 = 0.d0
      h0013 = 0.d0
      h0033 = 0.d0
c bring up to full accuracy
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 178 i = 1,ngangb
         pqab(i) = aqz-app(i)
         gg = 1.d0/(ep(i)+ecd)
         p(i) = gg*(pqab(i)*pqab(i)+qperp2)
         g(i) = gg*ecd
         q(i) = dp00p(i)/dsqrt(gp(i)+gcd)
 178  continue
      call compfs(ngangb,mxprms2,p,fs,iind,2)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 179 i = 1,ngangb
         gy = g(i)*q(i)
         ggy = g(i)*gy
         f0 = fs(i,1)*q(i)
         f1 = fs(i,2)*gy
         f2 = fs(i,3)*ggy
         h0000 = h0000+f0
         h0001 = h0001+f1
         h0003 = h0003-f1*pqab(i)
         h0011 = h0011+f2
         h0013 = h0013-f2*pqab(i)
         h0033 = h0033+f2*pqab(i)*pqab(i)
 179  continue
_IF1()c      do 180 i = 1,ngangb
_IF1()c      isml = ismlq+ismlp(i)
_IF1()c      if (isml .ge. 2) go to 180
_IF1()c      auxvar = var(isml+1)
_IF1()c      pqab = aqz-app(i)
_IF1()c      pqab2 = pqab*pqab
_IF1()c      g = 1.e0/(ep(i)+ecd)
_IF1()c      p = g*(pqab2+qperp2)
_IF1()c      g = g*ecd
_IF1()c      if (p .le. auxvar) go to 140
_IF1()c      f0 = dp00p(i)*sqrt (0.7853981625e0/(p*(gp(i)+gcd)))
_IF1()c      gtx = g/p
_IF1()c      f1 = 0.5e0*f0*gtx
_IF1()c      f2 = 1.5e0*f1*gtx
_IF1()c      go to 160
_IF1()c  140 q = dp00p(i)/sqrt (gp(i)+gcd)
_IF1()c      gy = g*q
_IF1()c      ggy = g*gy
_IF1()c      qq = p*12.5e0
_IF1()c      n =  int (qq)
_IF1()c      theta = qq- float(n)
_IF1()c      theta2 = theta*(theta-1.e0)
_IF1()c      theta3 = theta2*(theta-2.e0)
_IF1()c      theta4 = theta2*(theta+1.e0)
_IF1()c      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
_IF1()c      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
_IF1()c      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
_IF1()c  160 h0000 = h0000+f0
_IF1()c      h0001 = h0001+f1
_IF1()c      h0003 = h0003-f1*pqab
_IF1()c      h0011 = h0011+f2
_IF1()c      h0013 = h0013-f2*pqab
_IF1()c      h0033 = h0033+f2*pqab2
_IF1()c  180 continue
      h0022 = 0.5d0*ecd*(h0000-h0001)
      h0001 = h0001*qperp
      h0011 = h0011*qperp2+h0022
      h0013 = h0013*qperp
      h0033 = h0033+h0022
      if(sinp)120,100,120
 100  if(cosp)1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      g0011 = v44*h0011+v47*h0022
      g0012 = v54*h0011+v57*h0022
      g0022 = v74*h0011+v77*h0022
      g0013 = cosp*h0013
      g0023 = sinp*h0013
      g0033 = h0033
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0000 = h0000
      go to 2000
  920 g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
      go to 2000
1000  g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = -h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      g0010 = g0001
      g0020 = g0002
      g0021 = g0012
      g0030 = g0003
      g0031 = g0013
      g0032 = g0023
      if (rcdsq) 220,220,200
 200  g0010 = g0010+r13*g0000
      g0011 = g0011+r13*g0001
      g0012 = g0012+r13*g0002
      g0013 = g0013+r13*g0003
      g0030 = g0030+r33*g0000
      g0031 = g0031+r33*g0001
      g0032 = g0032+r33*g0002
      g0033 = g0033+r33*g0003
      g0001 = g0001+r14*g0000
      g0011 = g0011+r14*g0010
      g0021 = g0021+r14*g0020
      g0031 = g0031+r14*g0030
      g0003 = g0003+r34*g0000
      g0013 = g0013+r34*g0010
      g0023 = g0023+r34*g0020
      g0033 = g0033+r34*g0030
220   gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 5) = gout( 5)+g0010*dq10
      gout( 6) = gout( 6)+g0011*dq11
      gout( 7) = gout( 7)+g0012*dq11
      gout( 8) = gout( 8)+g0013*dq11
      gout( 9) = gout( 9)+g0020*dq10
      gout( 10) = gout( 10)+g0021*dq11
      gout( 11) = gout( 11)+g0022*dq11
      gout( 12) = gout( 12)+g0023*dq11
      gout( 13) = gout( 13)+g0030*dq10
      gout( 14) = gout( 14)+g0031*dq11
      gout( 15) = gout( 15)+g0032*dq11
      gout( 16) = gout( 16)+g0033*dq11
 940  continue
      ind = 0
      do 700 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
  700 continue
      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
      do 720 k = 1,4
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i3 ) = p13*t1+p23*t2+p33*t3
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
  720 continue
      return
      end
      subroutine sp0101(gout)
c        *****  special fast routine for -p- loop for 0101 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c
INCLUDE(common/astore)
INCLUDE(common/shllfo)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/const)
INCLUDE(common/geom)
c
      dimension p(mxprms2),pqab(mxprms2),q(mxprms2)
      dimension g(mxprms2),fs(mxprms2,3),iind(mxprms2)
      dimension gout(*)
      equivalence(v12,u12)
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0100 = 0.d0
      h0101 = 0.d0
      h0103 = 0.d0
      h0300 = 0.d0
      h0301 = 0.d0
      h0303 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 178 i = 1,ngangb
         pqab(i) = aqz-app(i)
         g(i) = 1.d0/(ep(i)+ecd)
         p(i) = (pqab(i)*pqab(i)+qperp2)*g(i)
         q(i) = conp(i)/dsqrt(gp(i)+gcd)
 178     continue
      call compfs(ngangb,mxprms2,p,fs,iind,2)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 179 i = 1,ngangb
         gy = g(i)*q(i)
         ggy = g(i)*gy
         f0 = fs(i,1)*q(i)
         f1 = fs(i,2)*gy
         f2 = fs(i,3)*ggy
         g03 = -pqab(i)*f1
         h0000 = h0000+f0 *dp00p(i)
         h0001 = h0001+f1 *dp00p(i)
         h0003 = h0003+g03*dp00p(i)
         h0100 = h0100-f1
         h0101 = h0101-f2
         h0103 = h0103+pqab(i)*f2
         h0300 = h0300-g03+bpp(i)*f0
         h0301 = h0301+bpp(i)*f1
         h0303 = h0303-pqab(i)*pqab(i)*f2+bpp(i)*g03
 179  continue
_IF1()c      do 180 i = 1,ngangb
_IF1()c      isml = ismlq+ismlp(i)
_IF1()c      if (isml .ge. 2) go to 180
_IF1()c      auxvar = var(isml+1)
_IF1()c      eab = ep(i)
_IF1()c      dp00 = dp00p(i)
_IF1()c      bp = bpp(i)
_IF1()c      pqab = aqz-app(i)
_IF1()c      pqab2 = pqab*pqab
_IF1()c      g = 1.e0/(ep(i)+ecd)
_IF1()c      p = (pqab2+qperp2)*g
_IF1()c      if (p .le. auxvar) go to 140
_IF1()c      f0 = conp(i)*sqrt (0.7853981625e0/(p*(gp(i)+gcd)))
_IF1()c      gtx = g/p
_IF1()c      f1 = 0.5e0*f0*gtx
_IF1()c      f2 = 1.5e0*f1*gtx
_IF1()c      go to 160
_IF1()c  140 q = conp(i)/sqrt (gp(i)+gcd)
_IF1()c      gy = g*q
_IF1()c      ggy = g*gy
_IF1()c      qq = p*12.5e0
_IF1()c      n =  int (qq)
_IF1()c      theta = qq- float(n)
_IF1()c      theta2 = theta*(theta-1.e0)
_IF1()c      theta3 = theta2*(theta-2.e0)
_IF1()c      theta4 = theta2*(theta+1.e0)
_IF1()c      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
_IF1()c      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
_IF1()c      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
_IF1()c  160 continue
_IF1()c      g03 = -pqab*f1
_IF1()c      h0000 = h0000+f0 *dp00
_IF1()c      h0001 = h0001+f1 *dp00
_IF1()c      h0003 = h0003+g03*dp00
_IF1()c      h0100 = h0100-f1
_IF1()c      h0101 = h0101-f2
_IF1()c      h0103 = h0103+pqab*f2
_IF1()c      h0300 = h0300-g03+bp*f0
_IF1()c      h0301 = h0301+bp*f1
_IF1()c      h0303 = h0303-pqab2*f2+bp*g03
_IF1()c  180 continue
      pp = qperp*ecd
      h0001 = h0001*pp
      h0003 = h0003*ecd
      h0202 = -0.5d0*ecd*h0100
      h0100 = h0100*qperp
      h0101 = h0101*qperp2*ecd
      h0103 = h0103*pp
      h0301 = h0301*pp
      h0303 = h0303*ecd
      h0301 = h0301+h0103
      h0101 = h0101+h0202
      h0303 = h0303+h0202
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      g0101 = cosp*h0101
      g0102 = sinp*h0101
      g0201 = v12*h0202
      g0202 = cosp*h0202
      g0301 = cosp*h0301
      g0302 = sinp*h0301
      g0303 = h0303
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0300 = h0300
      g0000 = h0000
      h0101 = g0101
      h0102 = g0102
      h0201 = g0201
      h0202 = g0202
      g0101 = cosp*h0101+u12*h0201
      g0102 = cosp*h0102+u12*h0202
      g0103 = cosp*h0103
      g0201 = sinp*h0101+cosp*h0201
      g0202 = sinp*h0102+cosp*h0202
      g0203 = sinp*h0103
      g0100 = cosp*h0100
      g0200 = sinp*h0100
      go to 2000
  920 g0100 = h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      go to 2000
 1000 g0100 = -h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = -h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = -h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
2000  r14 = dq*sing
      r34 = dq*cosg
      if (rcdsq) 720,720,700
  700 g0001 = g0001+r14*g0000
      g0101 = g0101+r14*g0100
      g0201 = g0201+r14*g0200
      g0301 = g0301+r14*g0300
      g0003 = g0003+r34*g0000
      g0103 = g0103+r34*g0100
      g0203 = g0203+r34*g0200
      g0303 = g0303+r34*g0300
  720 gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 17) = gout( 17)+g0100*dq00
      gout( 18) = gout( 18)+g0101*dq01
      gout( 19) = gout( 19)+g0102*dq01
      gout( 20) = gout( 20)+g0103*dq01
      gout( 33) = gout( 33)+g0200*dq00
      gout( 34) = gout( 34)+g0201*dq01
      gout( 35) = gout( 35)+g0202*dq01
      gout( 36) = gout( 36)+g0203*dq01
      gout( 49) = gout( 49)+g0300*dq00
      gout( 50) = gout( 50)+g0301*dq01
      gout( 51) = gout( 51)+g0302*dq01
      gout( 52) = gout( 52)+g0303*dq01
 940  continue
c ***
c ***
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 0101
c
c
c
_IFN1(f)      ind = 0
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 1 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   1  continue
_IFN1(f)      ind = -15
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 j = 1,4
_IFN1(f)      ind = ind+16
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      call mvml3(p11,1,gout(2),1,16,gout(2),1,16,4)
      return
      end
_EXTRACT(sp0111,ultra)
      subroutine sp0111(gout)
c        *****  special fast routine for -p- loop for 0111 *****
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c
      dimension g(64),h(64)
INCLUDE(common/shllfo)
INCLUDE(common/miscg)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/pgeom)
INCLUDE(common/qgeom)
INCLUDE(common/maxc)
INCLUDE(common/geom)
INCLUDE(common/const)
INCLUDE(common/astore)
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
INCLUDE(common/auxvar)
      dimension gout(*)
      dimension p(mxprms2),gabcd(mxprms2)
      dimension q(mxprms2),fs(mxprms2,4),iind(mxprms2)
      equivalence(v12,u12)
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
c        *****  begin -p- loop                   *****
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 179 i = 1,ngangb
         pqab = aqz-app(i)
         gabcd(i) = 1.d0/(ep(i)+ecd)
         p(i) = gabcd(i)*(pqab*pqab+qperp2)
         q(i) = conp(i)/dsqrt(gp(i)+gcd)
 179  continue
      call compfs(ngangb,mxprms2,p,fs,iind,3)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 180 i = 1,ngangb
         gy = gabcd(i)*q(i)
         ggy = gabcd(i)*gy
         gggy = gabcd(i)*ggy
         f0 = fs(i,1) * q(i)
         f1 = fs(i,2) * gy
         f2 = fs(i,3) * ggy
         f3 = fs(i,4) * gggy
         dp00 = dp00p(i)
         bp = bpp(i)
         pqab = aqz-app(i)
         pqab2 = pqab*pqab
_IF1()c      do 180 i = 1,ngangb
_IF1()c      isml = ismlq+ismlp(i)
_IF1()c      if (isml .ge. 2) go to 180
_IF1()c      auxvar = var(isml+1)
_IF1()c      eab = ep(i)
_IF1()c      dp00 = dp00p(i)
_IF1()c      bp = bpp(i)
_IF1()c      pqab = aqz-app(i)
_IF1()c      pqab2 = pqab*pqab
_IF1()c      gabcd = 1.e0/(eab+ecd)
_IF1()c      p = gabcd*(pqab2+qperp2)
_IF1()c     if (p .le. auxvar) go to 140
_IF1()c      f0 = conp(i)*sqrt (0.7853981625e0/(p*(gp(i)+gcd)))
_IF1()c      gtx = gabcd/p
_IF1()c      f1 = 0.5e0*f0*gtx
_IF1()c      f2 = 1.5e0*f1*gtx
_IF1()c      f3 = 2.5e0*f2*gtx
_IF1()c      go to 160
_IF1()c  140 q = conp(i)/sqrt (gp(i)+gcd)
_IF1()c      gy = gabcd*q
_IF1()c      ggy = gabcd*gy
_IF1()c      gggy = gabcd*ggy
_IF1()c      qq = p*12.5e0
_IF1()c      n =  int (qq)
_IF1()c      theta = qq- float(n)
_IF1()c      theta2 = theta*(theta-1.e0)
_IF1()c      theta3 = theta2*(theta-2.e0)
_IF1()c      theta4 = theta2*(theta+1.e0)
_IF1()c      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
_IF1()c      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
_IF1()c      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
_IF1()c      f3 = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)+theta4*c3(n+2))*gggy
_IF1()c  160 continue
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f2pqa2 = f2*pqab2
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      q1 = q1+f0 *bp
      q2 = q2+f1 *bp
      q3 = q3+f2 *bp
      q4 = q4+f1pqab*bp
      q5 = q5+f2pqab*bp
      q6 = q6+f2pqa2*bp
      r1 = r1+f1
      r2 = r2+f2
      r3 = r3+f3
      r4 = r4+f1pqab
      r5 = r5+f2pqab
      r6 = r6+f3pqab
      r7 = r7+f2pqa2
      r8 = r8+f3*pqab2
      r9 = r9+f3pqab*pqab2
  180 continue
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      qecd = qperp*ecd
      qecd2 = qperp*ecd2
      q2ecd = qperp2*ecd
      q2ecd2 = qperp2*ecd2
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = 0.5d0*qecd2*r2
      h( 44) = -0.5d0*ecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q2ecd2*qperp*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      g( 22) = v44*h( 22)+v47*h( 27)
      g( 23) = v54*h( 22)+v57*h( 27)
      g( 27) = v74*h( 22)+v77*h( 27)
      g( 24) = cosp*h( 24)
      g( 28) = sinp*h( 24)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = v12*h( 44)
      g( 44) = cosp*h( 44)
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g( 64) = h( 64)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 16) = h( 16)
      g( 18) = cosp*h( 18)
      g( 19) = sinp*h( 18)
      g( 20) = h( 20)
      g( 34) = v12*h( 35)
      g( 35) = cosp*h( 35)
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g( 52) = h( 52)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g(  4) = h(  4)
      g( 49) = h( 49)
      g(  1) = h(  1)
      h( 22) = g( 22)
      h( 23) = g( 23)
      h( 24) = g( 24)
      h( 27) = g( 27)
      h( 28) = g( 28)
      h( 38) = g( 38)
      h( 39) = g( 39)
      h( 40) = g( 40)
      h( 43) = g( 43)
      h( 44) = g( 44)
      h( 18) = g( 18)
      h( 19) = g( 19)
      h( 34) = g( 34)
      h( 35) = g( 35)
      g( 22) = cosp*h( 22)+u12*h( 38)
      g( 23) = cosp*h( 23)+u12*h( 39)
      g( 24) = cosp*h( 24)+u12*h( 40)
      g( 27) = cosp*h( 27)+u12*h( 43)
      g( 28) = cosp*h( 28)+u12*h( 44)
      g( 32) = cosp*h( 32)
      g( 38) = sinp*h( 22)+cosp*h( 38)
      g( 39) = sinp*h( 23)+cosp*h( 39)
      g( 40) = sinp*h( 24)+cosp*h( 40)
      g( 43) = sinp*h( 27)+cosp*h( 43)
      g( 44) = sinp*h( 28)+cosp*h( 44)
      g( 48) = sinp*h( 32)
      g( 18) = cosp*h( 18)+u12*h( 34)
      g( 19) = cosp*h( 19)+u12*h( 35)
      g( 20) = cosp*h( 20)
      g( 34) = sinp*h( 18)+cosp*h( 34)
      g( 35) = sinp*h( 19)+cosp*h( 35)
      g( 36) = sinp*h( 20)
      g( 17) = cosp*h( 17)
      g( 33) = sinp*h( 17)
      go to 2000
  920 g( 17) = h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = h( 20)
      g( 22) = h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = h( 27)
      g( 28) = dzero
      g( 32) = h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
      go to 2000
 1000 g( 17) = -h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = -h( 20)
      g( 22) = -h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = -h( 27)
      g( 28) = dzero
      g( 32) = -h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = -h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = -h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = -h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = -h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = -h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
2000  r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,50,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 720,720,700
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
700   do 701 kq1=1,49,16
          t1=g(kq1)
          t5=g(kq1+4)
          t13=g(kq1+12)
          t5=t5+r13*t1
          t13=t13+r33*t1
          t2=g(kq1+1)
          g(kq1+5)=g(kq1+5)+r13*t2+r14*t5
          g(kq1+13)=g(kq1+13)+r33*t2+r14*t13
          g(kq1+1)=t2+r14*t1
          t3=g(kq1+2)
          g(kq1+6)=g(kq1+6)+t3*r13
          g(kq1+14)=g(kq1+14)+t3*r33
          t9=g(kq1+8)
          g(kq1+9)=g(kq1+9)+r14*t9
          g(kq1+11)=g(kq1+11)+r34*t9
          t4=g(kq1+3)
          g(kq1+7)=g(kq1+7)+t4*r13+t5*r34
          g(kq1+4)=t5
          g(kq1+15)=g(kq1+15)+t4*r33+t13*r34
          g(kq1+12)=t13
701       g(kq1+3)=t4+t1*r34
720   do 721 kq1=2,62,4
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      dq01dd=dq01-dq11
      do 722 kq1=1,49,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq00
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq01dd
          gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq01dd
722       gout(kq1+3)=gout(kq1+3)+g(kq1+3)*dq01dd
      do 723 kq1=5,53,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq10
          gout(kq1+4)=gout(kq1+4)+g(kq1+4)*dq10
723       gout(kq1+8)=gout(kq1+8)+g(kq1+8)*dq10
  940 continue
_IFN1(f)      ind = 0
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 1 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   1  continue
_IFN1(f)      ind = -12
_IFN1(f)      do 2 j = 1,4
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 l = 1,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   2  continue
_IFN1(f)      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 3 loopjk = 1,16
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)   3  continue
_IF1(f)      call mvml3(p11,1,gout(17),16,1,gout(17),16,1,16)
_IF1(f)      ind=5
_IF1(f)      do 4 j=1,4
_IF1(f)      call  mvml3(p11,1,gout(ind),4,1,gout(ind),4,1,4)
_IF1(f)  4   ind=ind+16
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,16)
      return
      end
_ENDEXTRACT
_EXTRACT(sp1111,ultra)
c ******************************************************
c ******************************************************
c             =   sp1111  =
c ******************************************************
c ******************************************************
      subroutine sp1111(gout)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
_IF1(cfu)      parameter (dzero=0.0d0)
_IFN1(cfu)      parameter (dzero=0.0d0)
      dimension g(256),h(256)
      dimension kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
INCLUDE(common/miscg)
INCLUDE(common/shllfo)
INCLUDE(common/geom)
INCLUDE(common/pgeom)
INCLUDE(common/pqgeom)
INCLUDE(common/ginf)
INCLUDE(common/qgeom)
INCLUDE(common/const)
INCLUDE(common/maxc)
INCLUDE(common/astore)
      common/inttab/
     +  aa(333),ba(333),ca(333),abc1,
     +  ab(333),bb(333),cb(333),abc2,
     +  ac(333),bc(333),cc(333),abc3,
     +  ad(333),bd(333),cd(333),abc4,
     +  ae(333),be(333),ce(333),abc5,
     +  af(333),bf(333),cf(333),abc6
INCLUDE(common/auxvar)
      dimension gout(*)
      dimension p(mxprms2),q(mxprms2),gabcd(mxprms2)
      dimension fs(mxprms2,5),iind(mxprms2)
      equivalence(v12,u12)
      data done/1.0d0/
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c ***
c *** this is the fps version of sp1111.
c ***
c *** as much code as possible reduced to loops (>=4)
c *** to avoid ps cache misses and to enhance compiler
c *** optimisation. will probably run like a drain on
c *** the cray-1s
c ***
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
660   p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v4 = 0.d0
      v5 = 0.d0
      v6 = 0.d0
      w1 = 0.d0
      w2 = 0.d0
      w3 = 0.d0
      w4 = 0.d0
      w5 = 0.d0
      w6 = 0.d0
      w7 = 0.d0
      w8 = 0.d0
      w9 = 0.d0
      s1 = 0.d0
      s2 = 0.d0
      s3 = 0.d0
      s4 = 0.d0
      s6 = 0.d0
      s7 = 0.d0
      s8 = 0.d0
      s9 = 0.d0
      s10 = 0.d0
      s11 = 0.d0
      s12 = 0.d0
      s13 = 0.d0
      s14 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      t4 = 0.d0
      t5 = 0.d0
      t6 = 0.d0
      t7 = 0.d0
      t8 = 0.d0
      t9 = 0.d0
      t10 = 0.d0
      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
      t14 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0
      c5 = 0.d0
      c6 = 0.d0
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 179 ind = 1,ngangb
         pqab = aqz-app(ind)
         gabcd(ind) = 1.d0/(ep(ind)+ecd)
         p(ind) = gabcd(ind)*(qperp2+pqab*pqab)
         q(ind) = conp(ind)/dsqrt(gp(ind)+gcd)
 179  continue
      call compfs(ngangb,mxprms2,p,fs,iind,4)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 180 ind = 1,ngangb
         gy = gabcd(ind)*q(ind)
         ggy = gabcd(ind)*gy
         gggy = gabcd(ind)*ggy
         ggggy = gabcd(ind)*gggy
         f0 = fs(ind,1) * q(ind)
         f1 = fs(ind,2) * gy
         f2 = fs(ind,3) * ggy
         f3 = fs(ind,4) * gggy
         f4 = fs(ind,5) * ggggy
         eab = ep(ind)
         dp00 = dp00p(ind)
         dp01 = dp01p(ind)
         dp10 = dp10p(ind)
         ap = app(ind)
         bp = bpp(ind)
         pqab = aqz-ap
         pqab2 = pqab*pqab
_IF1()c      do 180 ind = 1,ngangb
_IF1()c      isml = ismlq+ismlp(ind)
_IF1()c      if (isml .ge. 2) go to 180
_IF1()c      auxvar = var(isml+1)
_IF1()c       eab = ep(ind)
_IF1()c      dp00 = dp00p(ind)
_IF1()c      dp01 = dp01p(ind)
_IF1()c      dp10 = dp10p(ind)
_IF1()c      ap = app(ind)
_IF1()c      bp = bpp(ind)
_IF1()c      pqab = aqz-ap
_IF1()c      pqab2 = pqab*pqab
_IF1()c      gabcd = 1.e0/(eab+ecd)
_IF1()c      p = gabcd*(qperp2+pqab2)
_IF1()c      if (p .le. auxvar) go to 140
_IF1()c      f0 = sqrt (0.7853981625e0/(p*(gp(ind)+gcd)))*conp(ind)
_IF1()c      gtx = gabcd/p
_IF1()c      f1 = 0.5e0*f0*gtx
_IF1()c      f2 = 1.5e0*f1*gtx
_IF1()c      f3 = 2.5e0*f2*gtx
_IF1()c      f4 = 3.5e0*f3*gtx
_IF1()c      go to 160
_IF1()c  140 q = conp(ind)/sqrt (gp(ind)+gcd)
_IF1()c      gy = gabcd*q
_IF1()c      ggy = gabcd*gy
_IF1()c      gggy = gabcd*ggy
_IF1()c      qq = p*12.5e0
_IF1()c      n =  int (qq)
_IF1()c      theta = qq- float(n)
_IF1()c      theta2 = theta*(theta-1.e0)
_IF1()c      theta3 = theta2*(theta-2.e0)
_IF1()c      theta4 = theta2*(theta+1.e0)
_IF1()c      f0 = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)+theta4*ca(n+2))*q
_IF1()c      f1 = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)+theta4*cb(n+2))*gy
_IF1()c      f2 = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)+theta4*cc(n+2))*ggy
_IF1()c      f3 = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)+theta4*cd(n+2))*gggy
_IF1()c      f4 = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)+theta4*ce(n+2))*gggy*
_IF1()c     &     gabcd
_IF1()c  160 continue
      apbp = ap*bp
      eab2 = eab*eab
      bpdp01 = bp*dp01
      apdp10 = ap*dp10
      edp01 = eab*dp01
      edp10 = eab*dp10
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f4pqab = f4*pqab
      f1pqa2 = f1*pqab2
      f2pqa2 = f2*pqab2
      f3pqa2 = f3*pqab2
      f4pqa2 = f4*pqab2
      f2pqa3 = f2pqa2*pqab
      f3pqa3 = f3pqa2*pqab
      f4pqa3 = f4pqa2*pqab
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      r1 = r1+f1 *edp01
      r2 = r2+f2 *edp01
      r3 = r3+f3 *edp01
      r4 = r4+f1pqab *edp01
      r5 = r5+f2pqab *edp01
      r6 = r6+f3pqab *edp01
      r7 = r7+f2pqa2 *edp01
      r8 = r8+f3pqa2 *edp01
      r9 = r9+f3pqa3 *edp01
      w1 = w1+f1 *edp10
      w2 = w2+f2 *edp10
      w3 = w3+f3 *edp10
      w4 = w4+f1pqab *edp10
      w5 = w5+f2pqab *edp10
      w6 = w6+f3pqab *edp10
      w7 = w7+f2pqa2 *edp10
      w8 = w8+f3pqa2 *edp10
      w9 = w9+f3pqa3 *edp10
      s1 = s1+f0 *eab
      s2 = s2+f1 *eab
      s3 = s3+f2 *eab
      s4 = s4+f3 *eab
      s6 = s6+f1pqab*eab
      s7 = s7+f2pqab*eab
      s8 = s8+f3pqab*eab
      s9 = s9+f1pqa2*eab
      s10 = s10+f2pqa2*eab
      s11 = s11+f3pqa2*eab
      s12 = s12+f2pqa3*eab
      s13 = s13+f3pqa3*eab
      s14 = s14+f3pqa3*pqab*eab
      t1 = t1+f0 *eab2
      t2 = t2+f1 *eab2
      t3 = t3+f2 *eab2
      t4 = t4+f3 *eab2
      t5 = t5+f4 *eab2
      t6 = t6+f2pqab*eab2
      t7 = t7+f3pqab*eab2
      t8 = t8+f4pqab*eab2
      t9 = t9+f2pqa2*eab2
      t10 = t10+f3pqa2*eab2
      t11 = t11+f4pqa2*eab2
      t12 = t12+f3pqa3*eab2
      t13 = t13+f4pqa3*eab2
      t14 = t14+f4pqa3*pqab*eab2
      if (rabsq .eq. 0.0d0) go to 180
      q1 = q1+f0 *bpdp01
      q2 = q2+f1 *bpdp01
      q3 = q3+f2 *bpdp01
      q4 = q4+f1pqab*bpdp01
      q5 = q5+f2pqab*bpdp01
      q6 = q6+f2pqa2*bpdp01
      v1 = v1+f0 *apdp10
      v2 = v2+f1 *apdp10
      v3 = v3+f2 *apdp10
      v4 = v4+f1pqab*apdp10
      v5 = v5+f2pqab*apdp10
      v6 = v6+f2pqa2*apdp10
      c1 = c1+f0 *apbp
      c2 = c2+f1 *apbp
      c3 = c3+f2 *apbp
      c4 = c4+f1pqab*apbp
      c5 = c5+f2pqab*apbp
      c6 = c6+f2pqa2*apbp
  180 continue
      a1 = aqz*s2-s6
      a2 = aqz*s3-s7
      a3 = aqz*s4-s8
      a4 = aqz*s6-s9
      a5 = aqz*s7-s10
      a6 = aqz*s8-s11
      a8 = aqz*s10-s12
      a9 = aqz*s11-s13
      a10 = aqz*s13-s14
      bqz = aqz-rab
      b1 = bqz*s2-s6
      b2 = bqz*s3-s7
      b3 = bqz*s4-s8
      b4 = bqz*s6-s9
      b5 = bqz*s7-s10
      b6 = bqz*s8-s11
      b8 = bqz*s10-s12
      b9 = bqz*s11-s13
      b10 = bqz*s13-s14
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      hecd2 = 0.5d0*ecd2
      qecd = qperp*ecd
      hqecd = 0.5d0*qecd
      qecd2 = qperp*ecd2
      hqecd2 = 0.5d0*qecd2
      q2ecd = qperp2*ecd
      q3ecd = qperp*q2ecd
      q2ecd2 = qperp2*ecd2
      q3ecd2 = q2ecd2*qperp
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = hqecd2*r2
      h( 44) = -hecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q3ecd2*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      h( 65) = -qperp*w1
      h(193) = w4+v1
      h(131) = hecd*w1
      h( 66) = h(131)-q2ecd*w2
      h( 68) = qecd*w5
      h(194) = h( 68)+qecd*v2
      h(196) = h(131)-ecd*w7-ecd*v4
      h(135) = hqecd2*w2
      h(140) = -hecd2*w5
      h( 75) = h(135)-qperp*h(131)
      h(203) = h(140)+hecd*(h(193)-ecd*v2)
      h( 72) = h(140)+q2ecd2*w6
      h(200) = h(135)-qecd2*(w8+v5)
      h( 70) = h( 75)+h(135)+h(135)-q3ecd2*w3
      h( 80) = h( 75)-qecd2*w8
      h(198) = h(203)+q2ecd2*(w6+v3)
      h(208) = h(203)+h(140)+h(140)+ecd2*(w9+v6)
      h(161) = 0.5d0*(s1-t2)
      h( 81) = h(161)+qperp2*t3
      h(113) = -qperp*(t6+b1)
      h(209) = -qperp*(t6+a1)
      h(241) = h(161)+t9+a4+b4+c1
      h(162) = hqecd*(s2-t3)
      h( 82) = h(162)-qecd*t3+q3ecd*t4
      temp = hecd*t6-q2ecd*t7
      h(114) = temp+hecd*b1-q2ecd*b2
      h(210) = temp+hecd*a1-q2ecd*a2
      h(242) = h(162)+qecd*(t10+a5+b5+c2)
      h( 99) = -hqecd*t3
      h(147) = h( 99)
      h(179) = hecd*(t6+b1)
      h(227) = hecd*(t6+a1)
      h(164) = hecd*(t6-s6)
      h( 84) = h(164)-q2ecd*t7
      temp = -hqecd*t3+qecd*t10
      h(116) = temp+qecd*b5
      h(212) = temp+qecd*a5
      h(244) = h(164)+ecd*(t6-t12-a8-b8-c4)+hecd*(a1+b1)
      h(103) = 0.25d0*ecd2*t3-0.5d0*q2ecd2*t4
      h(151) = h(103)
      h(183) = hqecd2*(t7+b2)
      h(231) = hqecd2*(t7+a2)
      h(108) = hqecd2*t7
      h(156) = h(108)
      h(188) = hecd2*(0.5d0*t3-t10-b5)
      h(236) = hecd2*(0.5d0*t3-t10-a5)
      hxxyy = 0.25d0*(ecd*(s1-t2)-ecd2*(s2-t3))
      h(171) = hxxyy+hecd2*t3
      h( 91) = hxxyy+0.5d0*(q2ecd*t3-q2ecd2*t4)
      temp = hqecd*(ecd*t7-t6)
      h(123) = temp+hqecd*(ecd*b2-b1)
      h(219) = temp+hqecd*(ecd*a2-a1)
      h(251) = hxxyy+hecd*(t9+a4+b4+c1)-hecd2*(t10+a5+b5+c2)
      h(166) = hxxyy+0.5d0*q2ecd2*(s3-t4)
      h( 86) = hxxyy+(hecd2+0.5d0*q2ecd)*t3+q2ecd2*(-3.d0*t4+
     +     0.5d0*s3+qperp2*t5)
      h(118) = 1.5d0*qecd2*(t7+b2)-hqecd*(t6+b1)-q3ecd2*(b3+t8)
      h(214) = 1.5d0*qecd2*(t7+a2)-hqecd*(t6+a1)-q3ecd2*(a3+t8)
      h(246) = hxxyy-hecd2*(qperp2*t4+t10+a5+b5)+hecd*(t9+a4+b4+c1-ecd*c
     +     2)+q2ecd2*(t11+0.5d0*s3+a6+b6+c3)
      h(168) = hqecd2*(t7-s7)
      h( 88) = 1.5d0*qecd2*t7-hqecd2*s7-q3ecd2*t8
      temp = hecd2*(0.5d0*t3-t10)+q2ecd2*(t11-0.5d0*t4)
      h(120) = temp-hecd2*b5+q2ecd2*b6
      h(216) = temp-hecd2*a5+q2ecd2*a6
      h(248) = qecd2*(1.5d0*t7-t13-a9-b9-c5)-hqecd2*(s7-a2-b2)
      h(176) = hxxyy+hecd2*(s10-t10)
      h( 96) = hxxyy-hecd2*(qperp2*t4+t10-s10)+0.5d0*q2ecd*t3+q2ecd2*
     +     t11
      h(128) = qecd2*(1.5d0*t7-t13-b9)-hqecd*(t6+b1)+hqecd2*b2
      h(224) = qecd2*(1.5d0*t7-t13-a9)-hqecd*(t6+a1)+hqecd2*a2
      h(256) = hxxyy+hecd2*(-3.d0*(a5+b5)+t3+s10-c2)+ecd2*(-3.d0*t10+
     +     t14+a10+b10+c6)+hecd*(t9+a4+b4+c1)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
_IF1(a)cvd$  shortloop
      do 103 kq1=22,214,48
          g(kq1  ) = v44*h(kq1) + v47*h(kq1+5)
          g(kq1+1) = v54*h(kq1) + v57*h(kq1+5)
103       g(kq1+5) = v74*h(kq1) + v77*h(kq1+5)
_IF1(a)cvd$  shortloop
      do 101 kq1=24,216,48
          g(kq1  ) = cosp*h(kq1)
101       g(kq1+4) = sinp*h(kq1)
_IF1(a)cvd$  shortloop
      do 102 kq1=18,210,48
          g(kq1  ) = cosp*h(kq1)
102       g(kq1+1) = sinp*h(kq1)
      g( 80) = h( 80)
      g( 86) = v44*h( 86)+v47*h( 91)
      g( 87) = v54*h( 86)+v57*h( 91)
      g( 91) = v74*h( 86)+v77*h( 91)
      g( 88) = cosp*h( 88)
      g( 92) = sinp*h( 88)
      g( 96) = h( 96)
      g(102) = v45*h(103)
      g(103) = v55*h(103)
      g(107) = -g(102)
      g(104) = v12*h(108)
      g(108) = cosp*h(108)
      g(112) = dzero
      g(128) = h(128)
      g(134) = v45*h(135)
      g(135) = v55*h(135)
      g(139) = -g(134)
      g(136) = v12*h(140)
      g(140) = cosp*h(140)
      g(144) = dzero
      g(150) = v45*h(151)
      g(151) = v55*h(151)
      g(155) = -g(150)
      g(152) = v12*h(156)
      g(156) = cosp*h(156)
      g(160) = dzero
      g(176) = h(176)
      g(182) = v45*h(183)
      g(183) = v55*h(183)
      g(187) = -g(182)
      g(184) = v12*h(188)
      g(188) = cosp*h(188)
      g(192) = dzero
      g(198) = v44*h(198)+v47*h(203)
      g(199) = v54*h(198)+v57*h(203)
      g(203) = v74*h(198)+v77*h(203)
      g(200) = cosp*h(200)
      g(204) = sinp*h(200)
      g(230) = v45*h(231)
      g(231) = v55*h(231)
      g(235) = -g(230)
      g(232) = v12*h(236)
      g(236) = cosp*h(236)
      g(240) = dzero
      g(246) = v44*h(246)+v47*h(251)
      g(247) = v54*h(246)+v57*h(251)
      g(251) = v74*h(246)+v77*h(251)
      g(248) = cosp*h(248)
      g(252) = sinp*h(248)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = v12*h( 44)
      g( 44) = cosp*h( 44)
      g( 48) = dzero
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 68) = h( 68)
      g( 82) = cosp*h( 82)
      g( 83) = sinp*h( 82)
      g( 84) = h( 84)
      g( 98) = v12*h( 99)
      g( 99) = cosp*h( 99)
      g(100) = dzero
      g(116) = h(116)
      g(130) = v12*h(131)
      g(131) = cosp*h(131)
      g(132) = dzero
      g(146) = v12*h(147)
      g(147) = cosp*h(147)
      g(148) = dzero
      g(164) = h(164)
      g(178) = v12*h(179)
      g(179) = cosp*h(179)
      g(180) = dzero
      g(194) = cosp*h(194)
      g(195) = sinp*h(194)
      g(226) = v12*h(227)
      g(227) = cosp*h(227)
      g(228) = dzero
      g(242) = cosp*h(242)
      g(243) = sinp*h(242)
      g( 34) = v12*h( 35)
      g( 35) = cosp*h( 35)
      g( 36) = dzero
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g( 65) = h( 65)
      g( 81) = h( 81)
      g( 97) = dzero
      g(113) = h(113)
      g(129) = dzero
      g(145) = dzero
      g(161) = h(161)
      g(177) = dzero
      g(225) = dzero
      g( 33) = dzero
      h( 80) = cosp*g( 80)
      h( 96) = cosp*g( 96)
      h(112) =           u12*g(176)
      h(128) = cosp*g(128)
      h(144) = sinp*g( 80)
      h(160) = sinp*g( 96)
      h(176) =           cosp*g(176)
      h(192) = sinp*g(128)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 121 kq1=70,118,16
          h(kq1   ) = cosp*g(kq1   ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1   ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+ 1) + u12*g(kq1+65)
          h(kq1+65) = sinp*g(kq1+ 1) + cosp*g(kq1+65)
          h(kq1+ 2) = cosp*g(kq1+ 2) + u12*g(kq1+66)
          h(kq1+66) = sinp*g(kq1+ 2) + cosp*g(kq1+66)
          h(kq1+ 5) = cosp*g(kq1+ 5) + u12*g(kq1+69)
          h(kq1+69) = sinp*g(kq1+ 5) + cosp*g(kq1+69)
          h(kq1+ 6) = cosp*g(kq1+ 6) + u12*g(kq1+70)
          h(kq1+70) = sinp*g(kq1+ 6) + cosp*g(kq1+70)
121   continue
      h( 68) = cosp*g( 68)
      h( 84) = cosp*g( 84)
      h(100) =           u12*g(164)
      h(116) = cosp*g(116)
      h(132) = sinp*g( 68)
      h(148) = sinp*g( 84)
      h(164) =           cosp*g(164)
      h(180) = sinp*g(116)
_IF1(a)cvd$  shortloop
      do 122 kq1=66,114,16
          h(kq1   ) = cosp*g(kq1  ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1  ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+1) + u12*g(kq1+65)
122       h(kq1+65) = sinp*g(kq1+1) + cosp*g(kq1+65)
_IF1(a)cvd$  shortloop
      do 1221 kq1=2,50,16
          h(kq1   ) = g(kq1   )
          h(kq1+ 1) = g(kq1+ 1)
          h(kq1+ 4) = g(kq1+ 4)
          h(kq1+ 5) = g(kq1+ 5)
          h(kq1+ 6) = g(kq1+ 6)
1221      h(kq1+ 9) = g(kq1+ 9)
_IF1(a)cvd$  shortloop
      do 1222 kq1=12,60,16
          h(kq1    ) = g(kq1    )
          h(kq1+182) = g(kq1+182)
          h(kq1+183) = g(kq1+183)
          h(kq1+186) = g(kq1+186)
          h(kq1+187) = g(kq1+187)
1222      h(kq1+188) = g(kq1+188)
_IF1(a)cvd$  shortloop
      do 1223 kq1=203,251,16
          h(kq1  ) = g(kq1  )
1223      h(kq1+1) = g(kq1+1)
      h( 65) = cosp*g( 65)
      h( 81) = cosp*g( 81)
      h( 97) =           u12*g(161)
      h(113) = cosp*g(113)
      h(129) = sinp*g( 65)
      h(145) = sinp*g( 81)
      h(161) =           cosp*g(161)
      h(177) = sinp*g(113)
      h( 48) = g( 48)
      h( 36) = g( 36)
      h(228) = g(228)
      h(240) = g(240)
      h(225) = g(225)
      h( 33) = g( 33)
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 123 kq1=22,214,64
          g(kq1   ) = cosp*h(kq1   ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1   ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+ 1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+ 1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+ 2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+ 2) + cosp*h(kq1+18)
          g(kq1+ 5) = cosp*h(kq1+ 5) + u12* h(kq1+21)
          g(kq1+21) = sinp*h(kq1+ 5) + cosp*h(kq1+21)
          g(kq1+ 6) = cosp*h(kq1+ 6) + u12* h(kq1+22)
          g(kq1+22) = sinp*h(kq1+ 6) + cosp*h(kq1+22)
          g(kq1+10) = cosp*h(kq1+10) + u12* h(kq1+26)
123       g(kq1+26) = sinp*h(kq1+10) + cosp* h(kq1+26)
_IF1(a)cvd$  shortloop
      do 124 kq1=17,209,64
          g(kq1   ) = cosp*h(kq1  ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1  ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+2) + cosp*h(kq1+18)
          g(kq1+ 3) = cosp*h(kq1+3) + u12* h(kq1+19)
124       g(kq1+19) = sinp*h(kq1+3) + cosp* h(kq1+19)
_IF1(a)cvd$  concur
      do 125 kq1=49,177,64
          kkq1=kq1
_IF1(a)cvd$  shortloop
          do 126 kkkq1=1,32
              g(kkq1)=h(kkq1)
              kkq1=kkq1+1
126       continue
125   continue
_IF1(a)cvd$  shortloop
      do 127 kq1=1,16
          g(kq1)=h(kq1)
127       g(kq1+240)=h(kq1+240)
      goto 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
920   do 921 kkq1=1,10
          kq1=kq1off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) = h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) = h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
921       g(kq1+15) = h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
      do 922 kkq1=1,6
          kq1=kq2off(kkq1)
          g(kq1   ) = dzero
          g(kq1+ 1) = dzero
          g(kq1+ 2) = h(kq1+ 2)
          g(kq1+ 3) = dzero
          g(kq1+ 5) = dzero
          g(kq1+ 6) = h(kq1+ 6)
          g(kq1+ 7) = dzero
          g(kq1+10) = dzero
          g(kq1+11) = h(kq1+11)
922       g(kq1+15) = dzero
      go to 2000
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
1000  do 1001 kkq1=1,6
          kq1=kq3off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) =-h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) =-h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
1001      g(kq1+15) = h(kq1+15)
      do 1002 kkq1=1,4
          kq1=kq4off(kkq1)
          g(kq1   ) = -h(kq1   )
          g(kq1+ 1) =  h(kq1+ 1)
          g(kq1+ 2) =  dzero
          g(kq1+ 3) = -h(kq1+ 3)
          g(kq1+ 5) = -h(kq1+ 5)
          g(kq1+ 6) =  dzero
          g(kq1+ 7) =  h(kq1+ 7)
          g(kq1+10) = -h(kq1+10)
          g(kq1+11) =  dzero
1002      g(kq1+15) = -h(kq1+15)
_IF1(a)cvd$  shortloop
_IF1(a)cvd$  nodepchk
      do 1003 kkq1=1,6
          kq1=kq5off(kkq1)
          g(kq1   ) =  dzero
          g(kq1+ 1) =  dzero
          g(kq1+ 2) =  h(kq1+ 2)
          g(kq1+ 3) =  dzero
          g(kq1+ 5) =  dzero
          g(kq1+ 6) = -h(kq1+ 6)
          g(kq1+ 7) =  dzero
          g(kq1+10) =  dzero
          g(kq1+11) =  h(kq1+11)
1003      g(kq1+15) =  dzero
          g(99)=-g(99)
          g(108)=-g(108)
          g(147)=-g(147)
          g(156)=-g(156)
          g(103)=-g(103)
          g(151)=-g(151)
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
_IF1(a)cvd$  shortloop
      do 2001 kq1=2,242,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 1200,1200,1300
_IF1(a)cvd$  concur
1300  do 1301 kq1=1,4
          kkq1=kq1
_IF1(a)cvd$  shortloop
          do 1302 jq1=1,16
              g(kkq1+4) = r13*g(kkq1) + g(kkq1+4)
              g(kkq1+12)= r33*g(kkq1) + g(kkq1+12)
1302          kkq1=kkq1+16
1301  continue
c ***
      do 1303 kq1=1,253,4
          g(kq1+1) = r14*g(kq1) + g(kq1+1)
1303      g(kq1+3) = r34*g(kq1) + g(kq1+3)
c ***
c1200  do 1201 kq1=1,256
c1201      gout(kq1) = dq11*g(kq1) + gout(kq1)
 1200  continue
       call daxpy(256,dq11,g,1,gout,1)
      dq01x=dq01-dq11
_IF1(a)cvd$  vector
      do 1202 kq1=2,242,16
          gout(kq1  ) = dq01x*g(kq1  ) + gout(kq1  )
          gout(kq1+1) = dq01x*g(kq1+1) + gout(kq1+1)
1202      gout(kq1+2) = dq01x*g(kq1+2) + gout(kq1+2)
      dq10x=dq10-dq11
      dq00x=dq00-dq11
_IF1(a)cvd$  vector
      do 1203 kq1=1,241,16
          gout(kq1   ) = dq00x*g(kq1   ) + gout(kq1   )
          gout(kq1+ 4) = dq10x*g(kq1+ 4) + gout(kq1+ 4)
          gout(kq1+ 8) = dq10x*g(kq1+ 8) + gout(kq1+ 8)
1203      gout(kq1+12) = dq10x*g(kq1+12) + gout(kq1+12)
940   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 1111
c
c
_IFN1(f)      i1 = 64
_IFN1(f)      i2 = 128
_IFN1(f)      i3 = 192
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)      do 1 jkl = 1,64
_IFN1(f)      i1 = i1+1
_IFN1(f)      i2 = i2+1
_IFN1(f)      i3 = i3+1
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    1 continue
_IFN1(f)      ind = -48
_IFN1(f)      do 2 i = 1,4
_IFN1(f)      ind = ind+48
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 2 loopkl = 1,16
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 16+ind
_IFN1(f)      i2 = 32+ind
_IFN1(f)      i3 = 48+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    2 continue
_IFN1(f)      ind = -12
_IFN1(f)      do 3 loopij = 1,16
_IFN1(f)      ind = ind+12
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(a)cvd$  shortloop
_IFN1(f)      do 3 l = 1,4
_IFN1(f)      ind = ind+1
_IFN1(f)      i1 = 4+ind
_IFN1(f)      i2 = 8+ind
_IFN1(f)      i3 = 12+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    3 continue
_IFN1(f)      ind = -3
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IFN1(f)      do 4 ijk = 1,64
_IFN1(f)      ind = ind+4
_IFN1(f)      i1 = 1+ind
_IFN1(f)      i2 = 2+ind
_IFN1(f)      i3 = 3+ind
_IFN1(f)      t1 = gout(i1)
_IFN1(f)      t2 = gout(i2)
_IFN1(f)      t3 = gout(i3)
_IFN1(f)      gout(i1 ) = p11*t1+p21*t2+p31*t3
_IFN1(f)      gout(i2 ) = p12*t1+p22*t2+p32*t3
_IFN1(f)      gout(i3 ) = p13*t1+p23*t2+p33*t3
_IFN1(f)    4 continue
_IF1(f)      call mvml3(p11,1,gout(65),64,1,gout(65),64,1,64)
_IF1(f)      ind=17
_IF1(f)      do 5 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),16,1,gout(ind),16,1,16)
_IF1(f)   5  ind=ind+64
_IF1(f)      ind=5
_IF1(f)      do 6 i=1,4
_IF1(f)      call mvml3(p11,1,gout(ind),4,16,gout(ind),4,16,16)
_IF1(f)   6  ind=ind+1
_IF1(f)      call mvml3(p11,1,gout(2),1,4,gout(2),1,4,64)
c
      return
_ENDEXTRACT
      end
_IF(ccpdft)
      subroutine iv0011(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,
     +           gout,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +           nfree,ifree,
     +           fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine iv0011(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,
     +           gout,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +           nfree,ifree)
_ENDIF
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
      dimension q(*),kv(*),lv(*),qq4v(*),gout(*)
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
c ***
INCLUDE(common/confus)
INCLUDE(common/cslosc)
INCLUDE(common/flip70)
INCLUDE(common/geom)
INCLUDE(common/iofile)
INCLUDE(common/maxc)
INCLUDE(common/pkfil)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/shlnos)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
INCLUDE(common/morokuma)
      common/type  /itype,jtype
      jtype=3
      l2 = nx
c *** i,j sp; k,l s.
      if(i.ge.j) then
         iish=i
         jjsh=j
         knew=i
         lnew=j
      else
         iish=j
         jjsh=i
         knew=j
         lnew=i
      endif
      ijsh=iish*4096+jjsh
      if(nkl.lt.7) then
         call sinset
         do 10 iab=1,nkl
            inew=kv(iab)
            jnew=lv(iab)
            kshell=kv(iab)
            lshell=lv(iab)
            klsh=kshell*4096+lshell
            if(ijsh.ge.klsh) then
               ishell=iish
               jshell=jjsh
               ib(1,1)=3
               ib(2,1)=4
               ib(3,1)=1
               ib(4,1)=2
            else
               ishell=kshell
               jshell=lshell
               kshell=iish
               lshell=jjsh
               ib(1,1)=1
               ib(2,1)=2
               ib(3,1)=3
               ib(4,1)=4
            endif
            call vclr(gout,1,16)
            qq4=qq4v(iab)
            call sinfo
            call sp0011(gout)
            if (odscf) then
              if (zscftp.eq.'uhf')then
_IF(ccpdft)
              call dir_build_uhf70(fock,fockb,
     +        dens,densb,gout,
     +        fac1,fac2, facex, ocoul, oexch)
_ELSE
              call dir_build_uhf70(fock,fockb,
     +        dens,densb,gout)
_ENDIF
              else if (zscftp.eq.'gvb') then
               if(nsheld.le.1) then
               call dir_build_open_70(fock,exch,dens,gout)
              else
               call dir_build_open2_70(l2,fock,exch,
     +                                 dens,gout)
              endif
              else
_IF(ccpdft)
_IF(cray)
              call qoutd70(fock,dens,gout,
     +        fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
              if(omorok) then
               call dbuild70_morok(fock,dens,gout)
              else
               call dbuild70(fock,dens,gout,
     +         fac1, fac2, facex, ocoul, oexch)
              endif
_ENDIF
_ELSE
_IF(cray)
              call qoutd70(fock,dens,gout)
_ELSE
              call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
              endif
            else
              call qout70(gout)
            endif
10       continue
         nkl=0
         return
      endif
c ***
      nnkl=nkl
      if((nnkl/2)*2.eq.nnkl) nnkl=nnkl+1
      nnkl=min(nnkl,nklmax)
      nnklng=nnkl*ngkngl
      igout=ifree
      ipmat=igout+nnkl*16
      irab=ipmat+nnkl*9
      igp=irab+nnkl
      iep=igp+nnklng
      idp00p=iep+nnklng
      iapp=idp00p+nnklng
      ibpp=iapp+nnklng
      iismlp=ibpp+nnklng
      iacx=iismlp+nnklng
      iacy=iacx+nnkl
      iacz=iacy+nnkl
      iacy2=iacz+nnkl
      icosg=iacy2+nnkl
      ising=icosg+nnkl
      iaqz=ising+nnkl
      iw1=iaqz+nnkl
      isinp=iw1+nnkl
      icosp=isinp+nnkl
      impaab=icosp+nnkl
      iisave=impaab+nnkl
      iqperp=iisave+nnkl
      ih=iqperp+nnkl
      ig=ih+nnkl*16
c *** arrays below overlap g(*,1...)
      ieab=ig
      idp00=ieab+nnkl
      iap=idp00+nnkl
      igpw=iap+nnkl
      iismlw=igpw+nnkl
      ipqab=iismlw+nnkl
      igabcd=ipqab+nnkl
      ip=igabcd+nnkl
      if0=ip+nnkl
      if1=if0+nnkl
      if2=if1+nnkl
c ***
      iused=ig+nnkl*16-ifree
      if(iused.gt.nfree) call caserr('core error in iv0011')
c ***
      c1=cpulft(1)
      call sinfov(jtype,knew,lnew,kv,lv,qq4v,ngkngl, ngc,cgg,csc,cpc,
     -ngd,dg,csd,cpd, q(irab),rcd,q(ipmat),dummy, q(igp),q(iep),q(idp00p
     -),dummy,dummy,q(iapp),q(ibpp), q(iacx),q(iacy),q(iacz),q(iacy2)
     -,q(icosg),q(ising), dummy, cmaxc,cmaxd,q(iismlp),error1,error2,
     - nkl,nnkl)
      c2=cpulft(1)
      cpus(17)=cpus(17)+c2-c1
c ***
      call vclr(q(igout),1,16*nnkl)
      call v0011(ngkngl,ngc,cgg,csc,cpc,ngd,dg,csd,cpd, rcd,q(ipmat)
     -,q(igp),q(iep),q(idp00p),q(iapp), q(iacx),q(iacy),q(iacz),q(iacy2)
     - ,q(icosg),q(ising), cmaxc,cmaxd,q(iismlp),error1,error2, nkl,
     -nnkl,q(ig),q(ih),q(impaab),q(iisave),q(iqperp),q(iaqz), q(iw1)
     -,q(isinp),q(icosp),q(ieab),q(idp00),q(iap),q(igpw), q(iismlw)
     -, q(ipqab),q(igabcd),q(ip),q(if0),q(if1),q(if2),q(igout))
c ***
      c3=cpulft(1)
      cpus(18)=cpus(18)+c3-c2
      do 30 iab=1,nkl
         inew=kv(iab)
         jnew=lv(iab)
         kshell=kv(iab)
         lshell=lv(iab)
         klsh=kshell*4096+lshell
         if(ijsh.ge.klsh) then
            ishell=iish
            jshell=jjsh
            ib(1,1)=3
            ib(2,1)=4
            ib(3,1)=1
            ib(4,1)=2
         else
            ishell=kshell
            jshell=lshell
            kshell=iish
            lshell=jjsh
            ib(1,1)=1
            ib(2,1)=2
            ib(3,1)=3
            ib(4,1)=4
         endif
         iipq=igout+iab-1
         do 20 ipq=1,16
            gout(ipq)=q(iipq)
20       iipq=iipq+nnkl
         if (odscf) then
           if (zscftp.eq.'uhf')then
_IF(ccpdft)
            call dir_build_uhf70(fock,fockb,
     +      dens,densb,gout,
     +      fac1,fac2, facex, ocoul, oexch)
_ELSE
            call dir_build_uhf70(fock,fockb,
     +      dens,dens,gout)
_ENDIF

           else if (zscftp.eq.'gvb') then
             if(nsheld.le.1) then
             call dir_build_open_70(fock,exch,dens,
     +                              gout)
             else
             call dir_build_open2_70(l2,fock,exch,
     +                              dens,gout)
             endif
           else
_IF(ccpdft)
_IF(cray)
             call qoutd70(fock,dens,gout,
     +       fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
             if(omorok) then
              call dbuild70_morok(fock,dens,gout)
             else
              call dbuild70(fock,dens,gout,
     +        fac1, fac2, facex, ocoul, oexch)
             endif
_ENDIF
_ELSE
_IF(cray)
             call qoutd70(fock,dens,gout)
_ELSE
             call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
           endif
         else
           call qout70(gout)
         endif
30    continue
      c4=cpulft(1)
      cpus(19)=cpus(19)+c4-c3
      nkl=0
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine v0011(
     +ngangb,
     +ngc,cg,csc,cpc,ngd,dg,csd,cpd,
     +rcd,pmat,
     +gp,ep,dp00p,app,
     +acx,acy,acz,acy2,cosg,sing,
     +cmaxc,cmaxd,ismlp,error1,error2,
c ***
     +nab,nnab,g,h,mpaab,isave,qperp,aqz,w1,
     +sinp,cosp,eab,dp00,ap,gpw,ismlpw,
     +pqab,gabcd,p,f0,f1,f2,gout)
c ***
_IF1(a)cvd$r altcode(32) vector
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c ***
      dimension cg(ngc),csc(ngc),cpc(ngc),dg(ngd),csd(ngd),cpd(ngd),
     +          pmat(9,nab),gp(nnab,ngangb),ep(nnab,ngangb),
     +          dp00p(nnab,ngangb),app(nnab,ngangb),
     +          acx(nab),acy(nab),acz(nab),acy2(nab),cosg(nab),
     +          sing(nab),cmaxc(*),cmaxd(*),ismlp(nnab,ngangb),
     +          g(nnab,16),h(nnab,16),mpaab(nab),isave(nab),qperp(nab),
     +          aqz(nab),w1(nab),sinp(nab),cosp(nab),eab(nab),dp00(nab),
     +          ap(nab),gpw(nab),ismlpw(nab),pqab(nab),gabcd(nab),
     +          p(nab),f0(nab),f1(nab),f2(nab),gout(nnab,16)
c ***
INCLUDE(common/auxvar)
      common/inttab/a0(333),b0(333),c0(333),abc1,a1(333),b1(333)
     -              ,c1(333),abc2,a2(333),b2(333),c2(333),abc3,a3(333)
     -              ,b3(333),c3(333),abc4,a4(333),b4(333),c4(333)
     -              ,abc5,a5(333),b5(333),c5(333),abc6
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
      do 400 k = 1,ngc
         gc = cg(k)
         do 401 l = 1,ngd
            gd = dg(l)
            gcd = gc+gd
            ecd = done/gcd
            cq = gd*ecd*rcd
            dq = cq-rcd
            qqq = cq*dq*gcd
            if (qqq+sixty) 10,20,20
10          v = 0.0d0
            go to 30
20          v = dexp(qqq)*ecd
30          qqtest = cmaxc(k)*cmaxd(l)*v
            if (qqtest-error1) 50,50,40
40          ismlq = 0
            go to 70
50          if (qqtest-error2) 401,401,60
60          ismlq = 1
70          sc = csc(k)
            sd = csd(l)
            pc = cpc(k)
            pd = cpd(l)
            dq00 = sc*sd*v
            dq01 = sc*pd*v
            dq10 = pc*sd*v
            dq11 = pc*pd*v
c ***
            n120 = 0
            n920 = 0
            n1000 = 0
_IF1(a)cvd$  novector
            do 80 iab=1,nab
               aqx = acx(iab)+sing(iab)*cq
               aqz(iab) = acz(iab)+cosg(iab)*cq
               qperp2 = aqx*aqx+acy2(iab)
               qperp(iab) = dsqrt(qperp2)
               w1(iab)=aqx
80          continue
            do 160 iab=1,nab
               if (qperp(iab)-tenm12) 100,100,90
90             cospp = -w1(iab)/qperp(iab)
               sinpp = -acy(iab)/qperp(iab)
               go to 110
100            cospp = done
               sinpp = 0.0d0
110            continue
c also form vector map around 120,920,1000 test.
               if (sinpp) 130,120,130
120            if (cospp) 150,130,140
130            n120=n120+1
               isave(iab)=n120
               mpaab(n120)=iab
               sinp(n120)=sinpp
               cosp(n120)=cospp
               goto 160
140            isave(iab)=nab-n920
               mpaab(nab-n920)=iab
               n920=n920+1
               goto 160
150            isave(iab)=1000
               n1000=n1000+1
160         continue
c ***
            m1000 = n120
            do 170 iab=1,nab
               if(isave(iab).eq.1000) then
                  m1000=m1000+1
                  mpaab(m1000)=iab
               endif
170         continue
c
_IF1(cu)            call gather(nab,w1,qperp,mpaab)
_IFN1(cfu)            call dgthr(nab,qperp,w1,mpaab)
_IF1(f)            call viindx(qperp,mpaab,1,w1,1,nab)
            do 180 iab=1,nab
               qperp(iab)=w1(iab)
180         continue
_IF1(cu)            call gather(nab,w1,aqz,mpaab)
_IFN1(cfu)            call dgthr(nab,aqz,w1,mpaab)
_IF1(f)            call viindx(aqz,mpaab,1,w1,1,nab)
            do 190 iab=1,nab
               aqz(iab)=w1(iab)
190         continue
c ***
            do 200 iab=1,nab
               h(iab,1) = 0.d0
               h(iab,2) = 0.d0
               h(iab,4) = 0.d0
               h(iab,6) = 0.d0
               h(iab,8) = 0.d0
               h(iab,16) = 0.d0
200         continue
c ***
_IF1(c)            call sensefi(mode)
_IF1(c)            call clearfi
c ***fortran gather scatter vectorises on xmp
            do 260 kkll=1,ngangb
               do 210 iab=1,nab
                  iiab=mpaab(iab)
                  eab(iab)=ep(iiab,kkll)
                  dp00(iab)=dp00p(iiab,kkll)
                  ap(iab)=app(iiab,kkll)
                  gpw(iab)=gp(iiab,kkll)
c                 ismlpw(iab)=ismlp(iiab,kkll)
210            continue
_IF1(a)cvd$  novector
               do 220 iab=1,nab
                  pqab(iab) = aqz(iab)-ap(iab)
                  pqab2 = pqab(iab)*pqab(iab)
                  gg = 1.d0/(eab(iab)+ecd)
                  p(iab) = gg*(pqab2+qperp(iab)*qperp(iab))
                  gabcd(iab) = gg*ecd
220            continue
_IF1()c               do 230 iab=1,nab
_IF1()cc *** temp fix to ensure that get same no.s as gaussian
_IF1()c                  temp=cvmgt(0.0,dp00(iab),((ismlpw(iab)+ismlq).ge.2))
_IF1()c                  f0(iab) = temp*sqrt(0.7853981625e0/(p(iab)*(gpw(iab)
_IF1()c     -            +gcd)))
_IF1()c                  gtx = gabcd(iab)/p(iab)
_IF1()c                  f1(iab) = 0.5e0*f0(iab)*gtx
_IF1()c                  f2(iab) = 1.5e0*f1(iab)*gtx
_IF1()c230            continue
_IF1()c               do 240 iab=1,nab
_IF1()c                  isml=ismlq+ismlpw(iab)
_IF1()c                  if(isml.ge.2) goto 240
_IF1()c                  auxvar = var(isml+1)
_IF1()c                  if(p(iab).gt.auxvar) goto 240
_IF1()c                  q = dp00(iab)/sqrt (gpw(iab)+gcd)
_IF1()c                  gy = gabcd(iab)*q
_IF1()c                  ggy = gabcd(iab)*gy
_IF1()c                  qq = p(iab)*12.5e0
_IF1()c                  n = int (qq)
_IF1()c                  theta = qq- float(n)
_IF1()c                  theta2 = theta*(theta-1.e0)
_IF1()c                  theta3 = theta2*(theta-2.e0)
_IF1()c                  theta4 = theta2*(theta+1.e0)
_IF1()c                  f0(iab)=(a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4
_IF1()c     -            *c0(n+2)) *q
_IF1()c                  f1(iab)=(a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4
_IF1()c     -            *c1(n+2)) *gy
_IF1()c                  f2(iab)=(a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4
_IF1()c     -            *c2(n+2)) *ggy
_IF1()c240            continue
_IF1()c *** assume here that f0...f2 are sequential in core
               call compfs(nab,nnab,p,f0,ismlpw,2)
_IF1(a)cvd$  novector
               do 231 iab=1,nab
                  q = dp00(iab)/dsqrt(gpw(iab)+gcd)
                  gy = gabcd(iab)*q
                  ggy = gabcd(iab)*gy
                  f0(iab)=f0(iab)*q
                  f1(iab)=f1(iab)*gy
                  f2(iab)=f2(iab)*ggy
231            continue
c
               do 250 iab=1,nab
                  h(iab,1) = h(iab,1)+f0(iab)
                  h(iab,2) = h(iab,2)+f1(iab)
                  h(iab,4) = h(iab,4)-f1(iab)*pqab(iab)
                  h(iab,6) = h(iab,6)+f2(iab)
                  h(iab,8) = h(iab,8)-(f2(iab)*pqab(iab))
                  h(iab,16) = h(iab,16)+(f2(iab)*pqab(iab))*pqab(iab)
250            continue
260         continue
c ***
_IF1(c)            if(mode.ne.0) call setfi
c ***
            do 280 iab=1,nab
               h(iab,11) = 0.5d0*ecd*(h(iab,1)-h(iab,2))
               h(iab,2) = h(iab,2)*qperp(iab)
               h(iab,6) = h(iab,6)*qperp(iab)*qperp(iab)+h(iab,11)
               h(iab,8) = h(iab,8)*qperp(iab)
               h(iab,16) = h(iab,16)+h(iab,11)
280         continue
c ***
            if(n120.eq.0) goto 300
            iablo=1
            iabhi=n120
            do 290 iab=iablo,iabhi
               v44 = cosp(iab)*cosp(iab)
               v77 = v44
               v47 = done-v44
               v74 = v47
               v54 = cosp(iab)*sinp(iab)
               v57 = -v54
               g(iab,6) = v44*h(iab,6)+v47*h(iab,11)
               g(iab,7) = v54*h(iab,6)+v57*h(iab,11)
               g(iab,11) = v74*h(iab,6)+v77*h(iab,11)
               g(iab,8) = cosp(iab)*h(iab,8)
               g(iab,12) = sinp(iab)*h(iab,8)
               g(iab,16) = h(iab,16)
               g(iab,2) = cosp(iab)*h(iab,2)
               g(iab,3) = sinp(iab)*h(iab,2)
               g(iab,4) = h(iab,4)
               g(iab,1) = h(iab,1)
290         continue
c ***
300         if(n920.eq.0) goto 320
            iablo=nab-n920+1
            iabhi=nab
            do 310 iab=iablo,iabhi
               g(iab,1) = h(iab,1)
               g(iab,2) = h(iab,2)
               g(iab,3) = dzero
               g(iab,4) = h(iab,4)
               g(iab,6) = h(iab,6)
               g(iab,7) = dzero
               g(iab,8) = h(iab,8)
               g(iab,11) = h(iab,11)
               g(iab,12) = dzero
               g(iab,16) = h(iab,16)
310         continue
c ***
320         if(n1000.eq.0) goto 340
            iablo=n120+1
            iabhi=n120+n1000
            do 330 iab=iablo,iabhi
               g(iab,1) = h(iab,1)
               g(iab,2) = -h(iab,2)
               g(iab,3) = dzero
               g(iab,4) = h(iab,4)
               g(iab,6) = h(iab,6)
               g(iab,7) = dzero
               g(iab,8) = -h(iab,8)
               g(iab,11) = h(iab,11)
               g(iab,12) = dzero
               g(iab,16) = h(iab,16)
330         continue
c ***
340         continue
c ***
c *** now recover from scatter,gather
c ***
            call dcopy(nnab*16,g,1,h,1)
_IF1(a)cvd$  nodepchk
            do 350 iab=1,nab
               iiab=mpaab(iab)
               do 351 ipq=1,16
                  g(iiab,ipq)=h(iab,ipq)
351            continue
350         continue
c ***
            do 360 iab=1,nab
               g(iab,5) = g(iab,2)
               g(iab,9) = g(iab,3)
               g(iab,10) = g(iab,7)
               g(iab,13) = g(iab,4)
               g(iab,14) = g(iab,8)
               g(iab,15) = g(iab,12)
360         continue
c ***
            if(rcd.le.0) goto 380
            do 370 iab=1,nab
               r13 = cq*sing(iab)
               r33 = cq*cosg(iab)
               r14 = dq*sing(iab)
               r34 = dq*cosg(iab)
               g(iab,5) = g(iab,5)+r13*g(iab,1)
               g(iab,6) = g(iab,6)+r13*g(iab,2)
               g(iab,7) = g(iab,7)+r13*g(iab,3)
               g(iab,8) = g(iab,8)+r13*g(iab,4)
               g(iab,13) = g(iab,13)+r33*g(iab,1)
               g(iab,14) = g(iab,14)+r33*g(iab,2)
               g(iab,15) = g(iab,15)+r33*g(iab,3)
               g(iab,16) = g(iab,16)+r33*g(iab,4)
               g(iab,2) = g(iab,2)+r14*g(iab,1)
               g(iab,6) = g(iab,6)+r14*g(iab,5)
               g(iab,10) = g(iab,10)+r14*g(iab,9)
               g(iab,14) = g(iab,14)+r14*g(iab,13)
               g(iab,4) = g(iab,4)+r34*g(iab,1)
               g(iab,8) = g(iab,8)+r34*g(iab,5)
               g(iab,12) = g(iab,12)+r34*g(iab,9)
               g(iab,16) = g(iab,16)+r34*g(iab,13)
370         continue
c ***
380         do 390 iab=1,nab
               gout(iab, 1) = gout(iab, 1)+g(iab,1)*dq00
               gout(iab, 2) = gout(iab, 2)+g(iab,2)*dq01
               gout(iab, 3) = gout(iab, 3)+g(iab,3)*dq01
               gout(iab, 4) = gout(iab, 4)+g(iab,4)*dq01
               gout(iab, 5) = gout(iab, 5)+g(iab,5)*dq10
               gout(iab, 6) = gout(iab, 6)+g(iab,6)*dq11
               gout(iab, 7) = gout(iab, 7)+g(iab,7)*dq11
               gout(iab, 8) = gout(iab, 8)+g(iab,8)*dq11
               gout(iab, 9) = gout(iab, 9)+g(iab,9)*dq10
               gout(iab, 10) = gout(iab, 10)+g(iab,10)*dq11
               gout(iab, 11) = gout(iab, 11)+g(iab,11)*dq11
               gout(iab, 12) = gout(iab, 12)+g(iab,12)*dq11
               gout(iab, 13) = gout(iab, 13)+g(iab,13)*dq10
               gout(iab, 14) = gout(iab, 14)+g(iab,14)*dq11
               gout(iab, 15) = gout(iab, 15)+g(iab,15)*dq11
               gout(iab, 16) = gout(iab, 16)+g(iab,16)*dq11
390         continue
401      continue
400   continue
c ***
      ind = 0
      do 410 l = 1,4
         ind = ind+1
         i1 = 4+ind
         i2 = 8+ind
         i3 = 12+ind
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodpchk
_IF1(i)c*vdir: ignore recrdeps
         do 411 iab=1,nab
            t1 = gout(iab,i1)
            t2 = gout(iab,i2)
            t3 = gout(iab,i3)
            gout(iab,i1 ) = pmat(1,iab)*t1+pmat(4,iab)*t2+pmat(7,iab)
     -      *t3
            gout(iab,i2 ) = pmat(2,iab)*t1+pmat(5,iab)*t2+pmat(8,iab)
     -      *t3
            gout(iab,i3 ) = pmat(3,iab)*t1+pmat(6,iab)*t2+pmat(9,iab)
     -      *t3
411      continue
410   continue
      ind = -3
      do 420 k = 1,4
         ind = ind+4
         i1 = 1+ind
         i2 = 2+ind
         i3 = 3+ind
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
         do 421 iab=1,nab
            t1 = gout(iab,i1)
            t2 = gout(iab,i2)
            t3 = gout(iab,i3)
            gout(iab,i3 ) = pmat(3,iab)*t1+pmat(6,iab)*t2+pmat(9,iab)
     -      *t3
            gout(iab,i1 ) = pmat(1,iab)*t1+pmat(4,iab)*t2+pmat(7,iab)
     -      *t3
            gout(iab,i2 ) = pmat(2,iab)*t1+pmat(5,iab)*t2+pmat(8,iab)
     -      *t3
421      continue
420   continue
      return
      end
_IF(ccpdft)
      subroutine iv0111(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,gout,
     +                  q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                  nfree,ifree,
     _                  fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine iv0111(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,gout,
     +                  q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                  nfree,ifree)
_ENDIF
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
      dimension q(*)
      dimension kv(*),lv(*),qq4v(*),gout(*)
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
c ***
c ***
INCLUDE(common/confus)
INCLUDE(common/cslosc)
INCLUDE(common/flip70)
INCLUDE(common/geom)
INCLUDE(common/iofile)
INCLUDE(common/maxc)
INCLUDE(common/pkfil)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/shlnos)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
INCLUDE(common/morokuma)
      common/type  /itype,jtype
      jtype=5
      l2 = nx
c **** i,j,l sp; k s.
      if(i.ge.j) then
         iish=i
         jjsh=j
         knew=i
         lnew=j
      else
         iish=j
         jjsh=i
         knew=j
         lnew=i
      endif
      ijsh=iish*4096+jjsh
      if(nkl.lt.7) then
c *** this is just to unscrew saved variables in sinfo
         call sinset
         do 10 iab=1,nkl
         c1=cpulft(1)
            inew=kv(iab)
            jnew=lv(iab)
            kshell=kv(iab)
            lshell=lv(iab)
            ib(3,1)=1
            ib(4,1)=2
            if(kshell.lt.lshell) then
               kshell=lv(iab)
               lshell=kv(iab)
               ib(3,1)=2
               ib(4,1)=1
            endif
            klsh=kshell*4096+lshell
            if(ijsh.ge.klsh) then
               ishell=iish
               jshell=jjsh
               ib(1,1)=3
               ib(2,1)=4
            else
               ishell=kshell
               jshell=lshell
               kshell=iish
               lshell=jjsh
               ib(1,1)=ib(3,1)
               ib(2,1)=ib(4,1)
               ib(3,1)=3
               ib(4,1)=4
            endif
            call vclr(gout,1,64)
            qq4=qq4v(iab)
            call sinfo
            c2=cpulft(1)
            cpus(14)=cpus(14)+c2-c1
            call sp0111(gout)
            c3=cpulft(1)
            cpus(15)=cpus(15)+c3-c2
            if (odscf) then
              if (zscftp.eq.'uhf')then
_IF(ccpdft)
               call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout,
     +         fac1,fac2, facex, ocoul, oexch)
_ELSE
               call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout)
_ENDIF
              else if (zscftp.eq.'gvb') then
               if(nsheld.le.1) then
               call dir_build_open_70(fock,exch,dens,
     +                                gout)
               else
               call dir_build_open2_70(l2,fock,exch,
     +                                dens,gout)
               endif
              else
_IF(ccpdft)
_IF(cray)
              call qoutd70(fock,dens,gout,
     +        fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
              if(omorok) then
               call dbuild70_morok(fock,dens,gout)
              else
               call dbuild70(fock,dens,gout,
     +         fac1, fac2, facex, ocoul, oexch)
              endif
_ENDIF
_ELSE
_IF(cray)
              call qoutd70(fock,dens,gout)
_ELSE
              call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
              endif
            else
              call qout70(gout)
            endif
            c4=cpulft(1)
            cpus(16)=cpus(16)+c4-c3
10       continue
         nkl=0
         return
      endif
c ***
      nnkl=nkl
      if( (nnkl/2)*2.eq.nnkl ) nnkl=nnkl+1
      nnkl=min(nnkl,nklmax)
      nnklng=nnkl*ngkngl
      igout=ifree
      irab=igout+nnkl*64
      ipmat=irab+nnkl
      iqmat=ipmat+nnkl*9
      igp=iqmat+nnkl*9
      iep=igp+nnklng
      idp00p=iep+nnklng
      idp01p=idp00p+nnklng
      iapp=idp01p+nnklng
      ibpp=iapp+nnklng
      iconp=ibpp+nnklng
      iismlp=iconp+nnklng
      iacx=iismlp+nnklng
      iacy=iacx+nnkl
      iacz=iacy+nnkl
      iacy2=iacz+nnkl
      icosg=iacy2+nnkl
      ising=icosg+nnkl
      iaqz=ising+nnkl
      isinp=iaqz+nnkl
      icosp=isinp+nnkl
      impaab=icosp+nnkl
      iisave=impaab+nnkl
      iqperp=iisave+nnkl
      ir13=iqperp+nnkl
      ir14=ir13+nnkl
      ir33=ir14+nnkl
      ir34=ir33+nnkl
      iaqx=ir34+nnkl
      ih=iaqx+nnkl
      ig=ih+nnkl*64
c *** arrays below overlap g(*,1...)
      ipn=ig
      iqn=ipn+nnkl*6
      irn=iqn+nnkl*6
      ieab=irn+nnkl*9
      idp00=ieab+nnkl
      iap=idp00+nnkl
      ibp=iap+nnkl
      igpw=ibp+nnkl
      iconpw=igpw+nnkl
      iismlw=iconpw+nnkl
      ipqab=iismlw+nnkl
      igabcd=ipqab+nnkl
      ip=igabcd+nnkl
      if0=ip+nnkl
      if1=if0+nnkl
      if2=if1+nnkl
      if3=if2+nnkl
      iqprpn=if3+nnkl
      iaqzn=iqprpn+nnkl
      iqecd2=iaqzn+nnkl
      iq2cd2=iqecd2+nnkl
      iq2ecd=iq2cd2+nnkl
c ***
      iused= ig+nnkl*64-ifree
      if(iused.gt.nfree) call caserr(' core error in iv0111.')
c ***
      c1=cpulft(1)
      call sinfov(jtype,knew,lnew,kv,lv,qq4v,ngkngl,ngc,cgg,csc,cpc,
     -ngd,dg,csd,cpd, q(irab),rcd,q(ipmat),q(iqmat), q(igp),q(iep)
     -,q(idp00p),q(idp01p),dummy,q(iapp),q(ibpp), q(iacx),q(iacy),
     -q(iacz),q(iacy2),q(icosg),q(ising), q(iconp), cmaxc,cmaxd,q(iismlp
     -),error1,error2, nkl,nnkl)
      c2=cpulft(1)
      cpus(14)=cpus(14)+c2-c1
c ***
      call vclr(q(igout),1,64*nnkl)
      call v0111(ngkngl,ngc,cgg,csc,cpc,ngd,dg,csd,cpd, rcd,q(ipmat)
     -, q(igp),q(iep),q(idp00p),q(iapp),q(ibpp), q(iacx),q(iacy),q(iacz)
     - ,q(iacy2),q(icosg),q(ising), q(iconp), cmaxc,cmaxd,q(iismlp)
     -,error1,error2, nkl,nnkl, q(ig),q(ih),q(impaab),q(iisave),q(iqperp
     -),q(iaqz),q(iaqx), q(isinp),q(icosp),q(ieab),q(idp00), q(iap)
     -,q(ibp),q(igpw),q(iconpw),q(iismlw),q(ipqab), q(igabcd),q(ip)
     -, q(if0),q(if1),q(if2),q(if3),q(ipn),q(iqn),q(irn), q(iqprpn)
     -,q(iaqzn),q(iqecd2),q(iq2cd2),q(iq2ecd), q(ir13),q(ir14),q(ir33)
     -,q(ir34),q(igout))
      c3=cpulft(1)
      cpus(15)=cpus(15)+c3-c2
      c3=cpulft(1)
      do 20 iab=1,nkl
         inew=kv(iab)
         jnew=lv(iab)
         kshell=kv(iab)
         lshell=lv(iab)
         ib(3,1)=1
         ib(4,1)=2
         if(kshell.lt.lshell) then
            kshell=lv(iab)
            lshell=kv(iab)
            ib(3,1)=2
            ib(4,1)=1
         endif
         klsh=kshell*4096+lshell
         if(ijsh.ge.klsh) then
            ishell=iish
            jshell=jjsh
            ib(1,1)=3
            ib(2,1)=4
         else
            ishell=kshell
            jshell=lshell
            kshell=iish
            lshell=jjsh
            ib(1,1)=ib(3,1)
            ib(2,1)=ib(4,1)
            ib(3,1)=3
            ib(4,1)=4
         endif
         call dcopy(64,q(igout+iab-1),nnkl,gout,1)
         if (odscf) then
           if (zscftp.eq.'uhf')then
_IF(ccpdft)
            call dir_build_uhf70(fock,fockb,
     +      dens,densb,gout,
     +      fac1,fac2, facex, ocoul, oexch)
_ELSE
            call dir_build_uhf70(fock,fockb,
     +      dens,densb,gout)
_ENDIF
           else if (zscftp.eq.'gvb') then
            if(nsheld.le.1) then
            call dir_build_open_70(fock,exch,dens,
     +                             gout)
            else
            call dir_build_open2_70(l2,fock,exch,
     +                             dens,gout)
            endif
           else
_IF(ccpdft)
_IF(cray)
           call qoutd70(fock,dens,gout,
     +     fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
           if(omorok) then
            call dbuild70_morok(fock,dens,gout)
           else
            call dbuild70(fock,dens,gout,
     +        fac1, fac2, facex, ocoul, oexch)
           endif
_ENDIF
_ELSE
_IF(cray)
           call qoutd70(fock,dens,gout)
_ELSE
           call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
           endif
         else
           call qout70(gout)
         endif
20    continue
      c4=cpulft(1)
      cpus(16)=cpus(16)+c4-c3
      nkl=0
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine v0111(
     + ngangb,
     + ngc,cg,csc,cpc,ngd,dg,csd,cpd,
     + rcd,pmat,
     + gp,ep,dp00p,app,bpp,
     + acx,acy,acz,acy2,cosg,sing,
     + conp,
     + cmaxc,cmaxd,ismlp,error1,error2,
     + nab,nnab,
     + g,h,mpaab,isave,qperp,aqz,aqx,
     + sinp,cosp,eab,dp00,ap,bp,
     + gpw,conpw,ismlpw,pqab,gabcd,p,f0,f1,f2,f3,
     + pn,qn,rn,qperpn,
     + aqzn,qecd2,q2ecd2,q2ecd,r13,r14,r33,r34,gout)
c ***
_IF1(a)cvd$r altcode(32) vector
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cg(ngc),csc(ngc),cpc(ngc),dg(ngd),csd(ngd),cpd(ngd),
     +          pmat(9,nab),gp(nnab,*),ep(nnab,*),dp00p(nnab,*),
     +          app(nnab,*),bpp(nnab,*),
     +          acx(nab),acy(nab),acz(nab),acy2(nab),cosg(nab),
     +          sing(nab),conp(nnab,*),cmaxc(ngc),cmaxd(ngd),
     +          ismlp(nnab,*),g(nnab,64),h(nnab,64),mpaab(nab),
     +          isave(nab),qperp(nab),aqz(nab),sinp(nab),cosp(nab),
     +          aqx(nab),eab(nab),dp00(nab),ap(nab),bp(nab),gpw(nab),
     +          conpw(nab),ismlpw(nab),pqab(nab),gabcd(nab),
     +          p(nab),f0(nab),f1(nab),f2(nab),f3(nab),
     +          pn(nnab,6),qn(nnab,6),rn(nnab,9),qperpn(nab),
     +          aqzn(nab),qecd2(nab),q2ecd2(nab),q2ecd(nab),
     +          r13(nab),r14(nab),r33(nab),r34(nab),
     +          gout(nnab,64)
c ***
_IF1()      common/auxvar/auxnnn,var(2)
      common/inttab/a0(333),b0(333),c0(333),abc1,a1(333),b1(333)
     -              ,c1(333),abc2,a2(333),b2(333),c2(333),abc3,a3(333)
     -              ,b3(333),c3(333),abc4,a4(333),b4(333),c4(333)
     -              ,abc5,a5(333),b5(333),c5(333),abc6
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      rcdsq=rcd*rcd
      do 640 k = 1,ngc
         gc = cg(k)
c
         do 630 l = 1,ngd
            gd = dg(l)
            gcd = gc+gd
            ecd = done/gcd
            cq = gd*ecd*rcd
            dq = cq-rcd
            qqq = cq*dq*gcd
            if (qqq+sixty) 10,20,20
10          v = 0.0d0
            go to 30
20          v =  dexp(qqq)*ecd
30          qqtest = cmaxc(k)*cmaxd(l)*v
            if (qqtest-error1) 50,50,40
40          ismlq = 0
            go to 70
50          if (qqtest-error2) 630,630,60
60          ismlq = 1
70          sc = csc(k)
            sd = csd(l)
            pc = cpc(k)
            pd = cpd(l)
            dq00 = sc*sd*v
            dq01 = sc*pd*v
            dq10 = pc*sd*v
            dq11 = pc*pd*v
c
_IF1(a)cvd$  novector
            do 80 iab=1,nab
               aqx(iab) = acx(iab)+sing(iab)*cq
               aqz(iab) = acz(iab)+cosg(iab)*cq
               qperp(iab) =  dsqrt(aqx(iab)*aqx(iab)+acy2(iab))
80          continue
c
c while making sinp also form mapping vectors to circumvent
c 120,1000,920 tests later.
c
            n120=0
            n920=0
            n1000=0
            do 160 iab=1,nab
               if (qperp(iab)-tenm12) 100,100,90
90             cospp = -aqx(iab)/qperp(iab)
               sinpp = -acy(iab)/qperp(iab)
               go to 110
100            cospp = done
               sinpp = 0.0d0
110            continue
c
c sort all 120's at front, 920's at end. 1000's in middle.
c
               if (sinpp) 130,120,130
120            if (cospp) 150,130,140
130            n120=n120+1
               isave(iab)=n120
               mpaab(n120)=iab
               sinp(n120)=sinpp
               cosp(n120)=cospp
               goto 160
140            isave(iab)=nab-n920
               mpaab(nab-n920)=iab
               n920=n920+1
               goto 160
150            isave(iab)=1000
               n1000=n1000+1
160         continue
c
            m1000=n120
            do 170 iab=1,nab
               if(isave(iab).eq.1000) then
                  m1000=m1000+1
                  mpaab(m1000)=iab
               endif
170         continue
c have to reorder qperp, aqz
_IF1(cu)            call gather(nab,qperpn,qperp,mpaab)
_IFN1(cuf)            call dgthr(nab,qperp,qperpn,mpaab)
_IF1(f)            call viindx(qperp,mpaab,1,qperpn,1,nab)
_IF1(cu)            call gather(nab,aqzn,aqz,mpaab)
_IFN1(cuf)            call dgthr(nab,aqz,aqzn,mpaab)
_IF1(f)            call viindx(aqz,mpaab,1,aqzn,1,nab)
c
            call vclr(pn,1,6*nnab)
            call vclr(qn,1,6*nnab)
            call vclr(rn,1,9*nnab)
c
_IF1(c)            call sensefi(mode)
_IF1(c)            call clearfi
            do 230 kkll = 1,ngangb
               do 180 iab=1,nab
                  iiab=mpaab(iab)
                  eab(iab)=ep(iiab,kkll)
                  dp00(iab)=dp00p(iiab,kkll)
                  ap(iab)=app(iiab,kkll)
                  bp(iab)=bpp(iiab,kkll)
                  gpw(iab)=gp(iiab,kkll)
                  conpw(iab)=conp(iiab,kkll)
c                 ismlpw(iab)=ismlp(iiab,kkll)
180            continue
c ***
_IF1(a)cvd$  novector
               do 190 iab=1,nab
                  pqab(iab)=aqzn(iab)-ap(iab)
                  pqab2=pqab(iab)*pqab(iab)
                  gabcd(iab)=1.0d0/(eab(iab)+ecd)
                  qperp2=qperpn(iab)*qperpn(iab)
190            p(iab)=gabcd(iab)*(pqab2+qperp2)
_IF1()c               do 200 iab=1,nab
_IF1()c                  temp=cvmgt(0.0,conpw(iab),((ismlpw(iab)+ismlq).ge.2)
_IF1()c     -            )
_IF1()c                  f0(iab)= sqrt(0.7853981625e0/(p(iab)*(gpw(iab)+gcd)
_IF1()c     -            ))* temp
_IF1()c                  gtx=gabcd(iab)/p(iab)
_IF1()c                  f1(iab) = 0.5e0*f0(iab)*gtx
_IF1()c                  f2(iab) = 1.5e0*f1(iab)*gtx
_IF1()c                  f3(iab) = 2.5e0*f2(iab)*gtx
_IF1()c200            continue
_IF1()c ***
_IF1()c               do 210 iab=1,nab
_IF1()c                  isml=ismlq+ismlpw(iab)
_IF1()c                  if(isml.ge.2) goto 210
_IF1()c                  auxvar=var(isml+1)
_IF1()c                  if(p(iab).gt.auxvar) goto 210
_IF1()c                  q=conpw(iab)/ sqrt(gpw(iab)+gcd)
_IF1()c                  gy = gabcd(iab)*q
_IF1()c                  ggy = gabcd(iab)*gy
_IF1()c                  gggy = gabcd(iab)*ggy
_IF1()c                  qq = p(iab)*12.5e0
_IF1()c                  n =   int(qq)
_IF1()c                  theta = qq- float(n)
_IF1()c                  theta2 = theta*(theta-1.0e0)
_IF1()c                  theta3 = theta2*(theta-2.0e0)
_IF1()c                  theta4 = theta2*(theta+1.0e0)
_IF1()c                  f0(iab) = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)
_IF1()c     -            + theta4*c0(n+2))*q
_IF1()c                  f1(iab) = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)
_IF1()c     -            + theta4*c1(n+2))*gy
_IF1()c                  f2(iab) = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)
_IF1()c     -            + theta4*c2(n+2))*ggy
_IF1()c                  f3(iab) = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)
_IF1()c     -            + theta4*c3(n+2))*gggy
_IF1()c210            continue
c *** assume that f0...f3 sequential in core
               call compfs(nab,nnab,p,f0,ismlpw,3)
_IF1(a)cvd$  novector
               do 231 iab=1,nab
                  q = conpw(iab)/dsqrt(gpw(iab)+gcd)
                  gy = gabcd(iab)*q
                  ggy = gabcd(iab)*gy
                  gggy = gabcd(iab)*ggy
                  f0(iab)=f0(iab)*q
                  f1(iab)=f1(iab)*gy
                  f2(iab)=f2(iab)*ggy
                  f3(iab)=f3(iab)*gggy
231            continue
c ***
               do 220 iab=1,nab
                  f1pqab = f1(iab)*pqab(iab)
                  f2pqab = f2(iab)*pqab(iab)
                  f3pqab = f3(iab)*pqab(iab)
                  f2pqa2 = f2pqab *pqab(iab)
                  pqab2 = pqab(iab)*pqab(iab)
                  pn(iab,1) = pn(iab,1)+f0(iab)*dp00(iab)
                  pn(iab,2) = pn(iab,2)+f1(iab)*dp00(iab)
                  pn(iab,3) = pn(iab,3)+f2(iab)*dp00(iab)
                  pn(iab,4) = pn(iab,4)+f1pqab *dp00(iab)
                  pn(iab,5) = pn(iab,5)+f2pqab *dp00(iab)
                  pn(iab,6) = pn(iab,6)+f2pqa2 *dp00(iab)
                  qn(iab,1) = qn(iab,1)+f0(iab)*bp(iab)
                  qn(iab,2) = qn(iab,2)+f1(iab)*bp(iab)
                  qn(iab,3) = qn(iab,3)+f2(iab)*bp(iab)
                  qn(iab,4) = qn(iab,4)+f1pqab *bp(iab)
                  qn(iab,5) = qn(iab,5)+f2pqab *bp(iab)
                  qn(iab,6) = qn(iab,6)+f2pqa2 *bp(iab)
                  rn(iab,1) = rn(iab,1)+f1(iab)
                  rn(iab,2) = rn(iab,2)+f2(iab)
                  rn(iab,3) = rn(iab,3)+f3(iab)
                  rn(iab,4) = rn(iab,4)+f1pqab
                  rn(iab,5) = rn(iab,5)+f2pqab
                  rn(iab,6) = rn(iab,6)+f3pqab
                  rn(iab,7) = rn(iab,7)+f2pqa2
                  rn(iab,8) = rn(iab,8)+f3(iab)*pqab2
                  rn(iab,9) = rn(iab,9)+f3pqab *pqab2
220            continue
230         continue
_IF1(c)            if(mode.ne.0) call setfi
c ***
            hecd = 0.5d0*ecd
            ecd2 = ecd*ecd
c
            do 240 iab=1,nab
               t1=qperpn(iab)
               t2=t1*t1
c qperp contains qecd
               qperp(iab) = t1*ecd
               qecd2(iab) = t1*ecd2
               q2ecd(iab) = t2*ecd
240         q2ecd2(iab) = t2*ecd2
c
            do 250 iab=1,nab
               h(iab, 1) = pn(iab,1)
               h(iab, 2) = qperp(iab)*pn(iab,2)
               temp = hecd*(pn(iab,1)-ecd*pn(iab,2))
               h(iab, 11)=temp
               h(iab, 16) = temp+ecd2*pn(iab,6)
250         h(iab, 6) = temp+q2ecd2(iab)*pn(iab,3)
c
            do 260 iab=1,nab
               h(iab, 4) = -ecd*pn(iab,4)
260         h(iab, 8) = -qecd2(iab)*pn(iab,5)
c
            do 270 iab=1,nab
               h(iab, 17) = -qperpn(iab)*rn(iab,1)
               temp = hecd*rn(iab,1)
               h(iab, 35)=temp
               h(iab, 18) = temp-q2ecd(iab)*rn(iab,2)
               h(iab, 52) = temp-ecd*rn(iab,7)-ecd*qn(iab,4)
               temp1 = 0.5d0*qecd2(iab)*rn(iab,2)
               h(iab, 39) = temp1
               temp2 = temp1-qperpn(iab)*temp
               h(iab, 27)=temp2
               h(iab, 56) = temp1-qecd2(iab)*(rn(iab,8)+qn(iab,5))
               h(iab, 22) = temp2+temp1+temp1- q2ecd2(iab)*qperpn(iab)
     -         *rn(iab,3)
               h(iab, 32) = temp2-qecd2(iab)*rn(iab,8)
270         continue
c
            do 280 iab=1,nab
               temp=rn(iab,4)+qn(iab,1)
               h(iab, 49) = temp
               temp1 = -0.5d0*ecd2*rn(iab,5)
               h(iab, 44)=temp1
               temp2 = temp1+hecd*(temp-ecd*qn(iab,2))
               h(iab, 59)=temp2
               temp3 = qperp(iab)*rn(iab,5)
               h(iab, 20)=temp3
               h(iab, 50) = temp3+qperp(iab)*qn(iab,2)
               h(iab, 24) = temp1+q2ecd2(iab)*rn(iab,6)
               h(iab, 54) = temp2+q2ecd2(iab)*(rn(iab,6)+qn(iab,3)
     -         )
               h(iab, 64) = temp2+temp1+temp1+ ecd2*(rn(iab,9)+qn(iab,
     -         6))
280         continue
c
c      if (sinp) 120,100,120
c  100 if (cosp) 1000,120,920
c  120 continue
c
            if(n120.eq.0) goto 390
c
            do 290 iab=1,n120
               v44 = cosp(iab)*cosp(iab)
               v77 = v44
               v47 = done-v44
               v74 = v47
               v54 = cosp(iab)*sinp(iab)
               v57 = -v54
               v45 = v57+v57
               v55 = v44-v47
               g(iab, 38) = v45*h(iab, 39)
               g(iab, 39) = v55*h(iab, 39)
               g(iab, 6) = v44*h(iab, 6)+v47*h(iab, 11)
               g(iab, 7) = v54*h(iab, 6)+v57*h(iab, 11)
               g(iab, 11) = v74*h(iab, 6)+v77*h(iab, 11)
               g(iab, 22) = v44*h(iab, 22)+v47*h(iab, 27)
               g(iab, 23) = v54*h(iab, 22)+v57*h(iab, 27)
               g(iab, 27) = v74*h(iab, 22)+v77*h(iab, 27)
               g(iab, 54) = v44*h(iab, 54)+v47*h(iab, 59)
               g(iab, 55) = v54*h(iab, 54)+v57*h(iab, 59)
               g(iab, 59) = v74*h(iab, 54)+v77*h(iab, 59)
290         continue
c
            do 300 iab=1,n120
               g(iab, 24) = cosp(iab)*h(iab, 24)
               g(iab, 28) = sinp(iab)*h(iab, 24)
               g(iab, 56) = cosp(iab)*h(iab, 56)
               g(iab, 60) = sinp(iab)*h(iab, 56)
               g(iab, 8) = cosp(iab)*h(iab, 8)
               g(iab, 12) = sinp(iab)*h(iab, 8)
300         continue
            do 310 iab=1,n120
               g(iab, 18) = cosp(iab)*h(iab, 18)
               g(iab, 19) = sinp(iab)*h(iab, 18)
               g(iab, 50) = cosp(iab)*h(iab, 50)
               g(iab, 51) = sinp(iab)*h(iab, 50)
               g(iab, 2) = cosp(iab)*h(iab, 2)
               g(iab, 3) = sinp(iab)*h(iab, 2)
310         continue
            do 320 iab=1,n120
               g(iab, 44) = cosp(iab)*h(iab, 44)
               g(iab, 40) =-sinp(iab)*h(iab, 44)
               g(iab, 35) = cosp(iab)*h(iab, 35)
               g(iab, 34) =-sinp(iab)*h(iab, 35)
320         continue
            do 330 iab=1,n120
               g(iab, 43) = -g(iab, 38)
               g(iab, 43) = -g(iab, 38)
               g(iab, 64) = h(iab, 64)
               g(iab, 16) = h(iab, 16)
               g(iab, 20) = h(iab, 20)
               g(iab, 52) = h(iab, 52)
               g(iab, 4) = h(iab, 4)
               g(iab, 49) = h(iab, 49)
               g(iab, 1) = h(iab, 1)
330         continue
c
            do 340 iab=1,n120
               h(iab, 22) = g(iab, 22)
               h(iab, 23) = g(iab, 23)
               h(iab, 24) = g(iab, 24)
               h(iab, 27) = g(iab, 27)
               h(iab, 28) = g(iab, 28)
               h(iab, 38) = g(iab, 38)
               h(iab, 39) = g(iab, 39)
               h(iab, 40) = g(iab, 40)
               h(iab, 43) = g(iab, 43)
340         continue
c
            do 350 iab=1,n120
               h(iab, 44) = g(iab, 44)
               h(iab, 18) = g(iab, 18)
               h(iab, 19) = g(iab, 19)
               h(iab, 34) = g(iab, 34)
               h(iab, 35) = g(iab, 35)
350         continue
c
            do 360 iab=1,n120
               g(iab, 18) = cosp(iab)*h(iab, 18)-sinp(iab)*h(iab,
     -         34)
               g(iab, 34) = sinp(iab)*h(iab, 18)+cosp(iab)*h(iab,
     -         34)
               g(iab, 19) = cosp(iab)*h(iab, 19)-sinp(iab)*h(iab,
     -         35)
               g(iab, 35) = sinp(iab)*h(iab, 19)+cosp(iab)*h(iab,
     -         35)
               g(iab, 22) = cosp(iab)*h(iab, 22)-sinp(iab)*h(iab,
     -         38)
               g(iab, 38) = sinp(iab)*h(iab, 22)+cosp(iab)*h(iab,
     -         38)
360         continue
            do 370 iab=1,n120
               g(iab, 23) = cosp(iab)*h(iab, 23)-sinp(iab)*h(iab,
     -         39)
               g(iab, 39) = sinp(iab)*h(iab, 23)+cosp(iab)*h(iab,
     -         39)
               g(iab, 24) = cosp(iab)*h(iab, 24)-sinp(iab)*h(iab,
     -         40)
               g(iab, 40) = sinp(iab)*h(iab, 24)+cosp(iab)*h(iab,
     -         40)
               g(iab, 27) = cosp(iab)*h(iab, 27)-sinp(iab)*h(iab,
     -         43)
               g(iab, 43) = sinp(iab)*h(iab, 27)+cosp(iab)*h(iab,
     -         43)
               g(iab, 28) = cosp(iab)*h(iab, 28)-sinp(iab)*h(iab,
     -         44)
               g(iab, 44) = sinp(iab)*h(iab, 28)+cosp(iab)*h(iab,
     -         44)
370         continue
            do 380 iab=1,n120
               g(iab, 32) = cosp(iab)*h(iab, 32)
               g(iab, 48) = sinp(iab)*h(iab, 32)
               g(iab, 20) = cosp(iab)*h(iab, 20)
               g(iab, 36) = sinp(iab)*h(iab, 20)
               g(iab, 17) = cosp(iab)*h(iab, 17)
               g(iab, 33) = sinp(iab)*h(iab, 17)
380         continue
c
c      go to 2000
c
390         if(n920.eq.0) goto 440
c
            iablo=nab-n920+1
            do 400 iab=iablo,nab
               g(iab, 19) = dzero
               g(iab, 23) = dzero
               g(iab, 28) = dzero
               g(iab, 33) = dzero
               g(iab, 34) = dzero
               g(iab, 36) = dzero
               g(iab, 38) = dzero
               g(iab, 40) = dzero
               g(iab, 43) = dzero
               g(iab, 48) = dzero
               g(iab, 51) = dzero
               g(iab, 55) = dzero
               g(iab, 60) = dzero
               g(iab, 3) = dzero
               g(iab, 7) = dzero
               g(iab, 12) = dzero
400         continue
            do 410 iab=iablo,nab
               g(iab, 17) = h(iab, 17)
               g(iab, 18) = h(iab, 18)
               g(iab, 20) = h(iab, 20)
               g(iab, 22) = h(iab, 22)
               g(iab, 24) = h(iab, 24)
               g(iab, 27) = h(iab, 27)
               g(iab, 32) = h(iab, 32)
               g(iab, 35) = h(iab, 35)
               g(iab, 39) = h(iab, 39)
410         continue
            do 420 iab=iablo,nab
               g(iab, 44) = h(iab, 44)
               g(iab, 49) = h(iab, 49)
               g(iab, 50) = h(iab, 50)
               g(iab, 52) = h(iab, 52)
               g(iab, 54) = h(iab, 54)
               g(iab, 56) = h(iab, 56)
               g(iab, 59) = h(iab, 59)
               g(iab, 64) = h(iab, 64)
               g(iab, 1) = h(iab, 1)
420         continue
            do 430 iab=iablo,nab
               g(iab, 2) = h(iab, 2)
               g(iab, 4) = h(iab, 4)
               g(iab, 6) = h(iab, 6)
               g(iab, 8) = h(iab, 8)
               g(iab, 11) = h(iab, 11)
               g(iab, 16) = h(iab, 16)
430         continue
c
c      go to 2000
c
440         if(n1000.eq.0) goto 490
            iablo=n120+1
            iabhi=n120+n1000
            do 450 iab=iablo,iabhi
               g(iab, 19) = dzero
               g(iab, 23) = dzero
               g(iab, 28) = dzero
               g(iab, 33) = dzero
               g(iab, 34) = dzero
               g(iab, 36) = dzero
               g(iab, 38) = dzero
               g(iab, 40) = dzero
               g(iab, 43) = dzero
               g(iab, 48) = dzero
               g(iab, 51) = dzero
               g(iab, 55) = dzero
               g(iab, 60) = dzero
               g(iab, 3) = dzero
               g(iab, 7) = dzero
               g(iab, 12) = dzero
450         continue
            do 460 iab=iablo,iabhi
               g(iab, 17) = -h(iab, 17)
               g(iab, 20) = -h(iab, 20)
               g(iab, 22) = -h(iab, 22)
               g(iab, 27) = -h(iab, 27)
               g(iab, 32) = -h(iab, 32)
               g(iab, 39) = -h(iab, 39)
               g(iab, 50) = -h(iab, 50)
               g(iab, 56) = -h(iab, 56)
               g(iab, 2) = -h(iab, 2)
               g(iab, 8) = -h(iab, 8)
460         continue
            do 470 iab=iablo,iabhi
               g(iab, 18) = h(iab, 18)
               g(iab, 24) = h(iab, 24)
               g(iab, 35) = h(iab, 35)
               g(iab, 44) = h(iab, 44)
               g(iab, 49) = h(iab, 49)
               g(iab, 52) = h(iab, 52)
               g(iab, 54) = h(iab, 54)
               g(iab, 59) = h(iab, 59)
               g(iab, 64) = h(iab, 64)
470         continue
            do 480 iab=iablo,iabhi
               g(iab, 1) = h(iab, 1)
               g(iab, 4) = h(iab, 4)
               g(iab, 6) = h(iab, 6)
               g(iab, 11) = h(iab, 11)
               g(iab, 16) = h(iab, 16)
480         continue
c
c
490         continue
c
c undo the gathering forced by 120/920/1000
c
            call dcopy(nnab*64,g,1,h,1)
_IF1(a)cvd$  nodepchk
            do 510 iab=1,nab
               iiab=mpaab(iab)
               do 500 i=1,64
500            g(iiab,i)=h(iab,i)
510         continue
c
            do 520 iab=1,nab
               r13(iab) = cq*sing(iab)
               r14(iab) = dq*sing(iab)
               r33(iab) = cq*cosg(iab)
               r34(iab) = dq*cosg(iab)
520         continue
c
            do 540 kq1=2,50,16
               do 530 iab=1,nab
                  g(iab,kq1+ 3) = g(iab,kq1 )
                  g(iab,kq1+ 7) = g(iab,kq1+ 1)
                  g(iab,kq1+ 8) = g(iab,kq1+ 5)
                  g(iab,kq1+11) = g(iab,kq1+ 2)
                  g(iab,kq1+12) = g(iab,kq1+ 6)
                  g(iab,kq1+13) = g(iab,kq1+10)
530            continue
540         continue
c
            if (rcdsq) 580,580,550
c
550         do 570 kq1=1,49,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
               do 560 iab=1,nab
                  t1=g(iab,kq1)
                  t5=g(iab,kq1+4)
                  t13=g(iab,kq1+12)
                  t5=t5+r13(iab)*t1
                  t13=t13+r33(iab)*t1
                  t2=g(iab,kq1+1)
                  g(iab,kq1+5)=g(iab,kq1+5)+r13(iab)*t2+r14(iab)*t5
                  g(iab,kq1+13)=g(iab,kq1+13)+r33(iab)*t2+r14(iab)
     -            *t13
                  g(iab,kq1+1)=t2+r14(iab)*t1
                  t3=g(iab,kq1+2)
                  g(iab,kq1+6)=g(iab,kq1+6)+t3*r13(iab)
                  g(iab,kq1+14)=g(iab,kq1+14)+t3*r33(iab)
                  t9=g(iab,kq1+8)
                  g(iab,kq1+9)=g(iab,kq1+9)+r14(iab)*t9
                  g(iab,kq1+11)=g(iab,kq1+11)+r34(iab)*t9
                  t4=g(iab,kq1+3)
                  g(iab,kq1+7)=g(iab,kq1+7)+t4*r13(iab)+t5*r34(iab)
                  g(iab,kq1+4)=t5
                  g(iab,kq1+15)=g(iab,kq1+15)+t4*r33(iab)+t13*r34(iab)
                  g(iab,kq1+12)=t13
                  g(iab,kq1+3)=t4+t1*r34(iab)
560            continue
570         continue
c
580         continue
            do 620 kq1=1,49,16
               do 590 iab=1,nab
590            gout(iab,kq1 )=gout(iab,kq1 ) + g(iab,kq1 )*dq00
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
               do 600 iab=1,nab
                  gout(iab,kq1+4 )=gout(iab,kq1+4 ) + g(iab,kq1+4
     -            )*dq10
                  gout(iab,kq1+8 )=gout(iab,kq1+8 ) + g(iab,kq1+8
     -            )*dq10
600            gout(iab,kq1+12)=gout(iab,kq1+12) + g(iab,kq1+12)*dq10
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
               do 610 iab=1,nab
                  gout(iab,kq1+1 )=gout(iab,kq1+1 ) + g(iab,kq1+1
     -            )*dq01
                  gout(iab,kq1+2 )=gout(iab,kq1+2 ) + g(iab,kq1+2
     -            )*dq01
                  gout(iab,kq1+3 )=gout(iab,kq1+3 ) + g(iab,kq1+3
     -            )*dq01
                  gout(iab,kq1+5 )=gout(iab,kq1+5 ) + g(iab,kq1+5
     -            )*dq11
                  gout(iab,kq1+6 )=gout(iab,kq1+6 ) + g(iab,kq1+6
     -            )*dq11
                  gout(iab,kq1+7 )=gout(iab,kq1+7 ) + g(iab,kq1+7
     -            )*dq11
                  gout(iab,kq1+9 )=gout(iab,kq1+9 ) + g(iab,kq1+9
     -            )*dq11
                  gout(iab,kq1+10)=gout(iab,kq1+10) + g(iab,kq1+10)
     -            *dq11
                  gout(iab,kq1+11)=gout(iab,kq1+11) + g(iab,kq1+11)
     -            *dq11
                  gout(iab,kq1+13)=gout(iab,kq1+13) + g(iab,kq1+13)
     -            *dq11
                  gout(iab,kq1+14)=gout(iab,kq1+14) + g(iab,kq1+14)
     -            *dq11
610            gout(iab,kq1+15)=gout(iab,kq1+15) + g(iab,kq1+15)*dq11
620         continue
c
630      continue
640   continue
c
      do 690 iab=1,nab
         p11=pmat(1,iab)
         p12=pmat(2,iab)
         p13=pmat(3,iab)
         p21=pmat(4,iab)
         p22=pmat(5,iab)
         p23=pmat(6,iab)
         p31=pmat(7,iab)
         p32=pmat(8,iab)
         p33=pmat(9,iab)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
         do 650 ig=17,32
            t1 = gout(iab,ig)
            t2 = gout(iab,ig+16)
            t3 = gout(iab,ig+32)
            gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
            gout(iab,ig+16) = p12*t1 + p22*t2 + p32*t3
650      gout(iab,ig+32) = p13*t1 + p23*t2 + p33*t3
cfps      call mvml3(pmat(1,iab),1,gout(iab,17),16*nnab,nnab,
cfps +                             gout(iab,17),16*nnab,nnab,16)
c          call calmv3(pmat(1,iab),gout(iab,17),16*nnab,nnab,16)
         ind=5
         do 670 j=1,4
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
            do 660 ig=ind,ind+3
               t1 = gout(iab,ig)
               t2 = gout(iab,ig+4)
               t3 = gout(iab,ig+8)
               gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
               gout(iab,ig+4 ) = p12*t1 + p22*t2 + p32*t3
660         gout(iab,ig+8 ) = p13*t1 + p23*t2 + p33*t3
cfps          call  mvml3(pmat(1,iab),1,gout(iab,ind),4*nnab,nnab,
cfps +                                  gout(iab,ind),4*nnab,nnab,4)
c              call calmv3(pmat(1,iab),gout(iab,ind),4*nnab,nnab,4)
670      ind=ind+16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
         do 680 ig=2,62,4
            t1 = gout(iab,ig)
            t2 = gout(iab,ig+1)
            t3 = gout(iab,ig+2)
            gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
            gout(iab,ig+1 ) = p12*t1 + p22*t2 + p32*t3
680      gout(iab,ig+2 ) = p13*t1 + p23*t2 + p33*t3
cfps      call mvml3(pmat(1,iab),1,gout(iab,2),nnab,4*nnab,
cfps +                             gout(iab,2),nnab,4*nnab,16)
c          call calmv3(pmat(1,iab),gout(iab,2),nnab,4*nnab,16)
690   continue
      return
      end
      subroutine sinfov(jtype,ish,jsh,kshv,lshv,qq4,
c *** /miscg/
     +ngangb,
c *** /shllfo/
     +ngc,cgg,csc,cpc,ngd,dg,csd,cpd,
c *** /geom/
     +rab,rcd,pmat,qmat,
c *** /pgeom/
     +gp,ep,dp00p,dp01p,dp10p,app,bpp,
c *** /qgeom/
     +acx,acy,acz,acy2,cosg,sing,
c *** /const/
     +conp,
c *** /maxc/
     +cmaxc,cmaxd,ismlp,error1,error2,
c ***
c ***
     +nab,nnab)
c ***
_IF1(a)cvd$r novector
_IF1(a)cvd$r noconcur
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/bshell)
      common/maxc  /cmax(mxprim)
INCLUDE(common/nshel)
INCLUDE(common/picon)
c ***
      dimension kshv(nab),lshv(nab),qq4(nab),
     +          cgg(*),csc(*),cpc(*),dg(*),csd(*),cpd(*),
     +          rab(nab),pmat(9,nab),qmat(9,nab),
     +          gp(nnab,1),ep(nnab,1),dp00p(nnab,1),dp01p(nnab,1),
     +          dp10p(nnab,1),app(nnab,1),bpp(nnab,1),
     +          acx(nab),acy(nab),acz(nab),acy2(nab),cosg(nab),
     +          conp(nnab,ngangb),sing(nab),
     +          cmaxc(*),cmaxd(*),ismlp(nnab,1)
c ***
      dimension cmaxa(mxprms),csa(mxprms),cpa(mxprms),ag(mxprms)
      dimension cmaxb(mxprms),csb(mxprms),cpb(mxprms),bg(mxprms)
      data sixty/60.0d0/
      data dzero/0.0d0/,pt5/0.5d0/,pt7/0.7d0/,pt9/0.9d0/,done/1.0d0/
      data pt0001/1.0d-04/,tenm12/1.0d-12/
c ***
      knew=ish
      lnew=jsh
      k = kstart(knew)
      ngc = kng(knew)
      do 10 nk = 1,ngc
         n = k-1+nk
         cmaxc(nk) = cmax(n)
         cgg(nk) = ex(n)
         csc(nk) = cs(n)
10    cpc(nk) = cp(n)
      l = kstart(lnew)
      ngd = kng(lnew)
      do 20 nl = 1,ngd
         n = l-1+nl
         cmaxd(nl) = cmax(n)
         dg(nl) = ex(n)
         csd(nl) = cs(n)
20    cpd(nl) = cp(n)
c ***
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
      rcd =  dsqrt(rcdsq)
      if (rcd) 30,40,30
30    q31 = cdx/rcd
      q32 = cdy/rcd
      q33 = cdz/rcd
      go to 50
40    q31 = dzero
      q32 = dzero
      q33 = done
50    continue
c ***
      ipt=1
      do 350iab=1,nab
         inew=kshv(iab)
         jnew=lshv(iab)
         i = kstart(inew)
         j = kstart(jnew)
         nga = kng(inew)
         ngb = kng(jnew)
         do 60 ni = 1,nga
            n = i-1+ni
            cmaxa(ni) = cmax(n)
            ag(ni) = ex(n)
            csa(ni) = cs(n)*qq4(iab)
60       cpa(ni) = cp(n)*qq4(iab)
         do 70 nj = 1,ngb
            n = j-1+nj
            cmaxb(nj) = cmax(n)
            bg(nj) = ex(n)
            csb(nj) = cs(n)
70       cpb(nj) = cp(n)
c ***
c
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
         rab(iab) =  dsqrt(rabsq)
         if (rab(iab)) 80,90,80
80       p31 = abx/rab(iab)
         p32 = aby/rab(iab)
         p33 = abz/rab(iab)
         go to 100
90       p31 = dzero
         p32 = dzero
         p33 = done
100      continue
         cosgg = p31*q31+p32*q32+p33*q33
         cosgg = dmin1(done,cosgg)
         cosgg = dmax1(-done,cosgg)
         singg =  dsqrt(done-cosgg*cosgg)
         cosg(iab)=cosgg
         sing(iab)=singg
         if (  dabs(cosgg)-pt9) 140,140,110
110      ppq1 = p31+q31
         ppq2 = p32+q32
         ppq3 = p33+q33
         pmq1 = p31-q31
         pmq2 = p32-q32
         pmq3 = p33-q33
         p21 = pmq2*ppq3-ppq2*pmq3
         p22 = pmq3*ppq1-ppq3*pmq1
         p23 = pmq1*ppq2-ppq1*pmq2
         p2 =  dsqrt(p21*p21+p22*p22+p23*p23)
         singg = pt5*p2
         sing(iab)=singg
         if (singg-tenm12) 130,120,120
120      temp = done/p2
         p21 = p21*temp
         p22 = p22*temp
         p23 = p23*temp
         go to 170
130      if (  dabs(p31)-pt7) 160,160,150
140      p21 = (p32*q33-p33*q32)/singg
         p22 = (p33*q31-p31*q33)/singg
         p23 = (p31*q32-p32*q31)/singg
         go to 170
150      p3333 = p33*p33
         p3333 = dmin1(done,p3333)
         sinp =  dsqrt(done-p3333)
         p21 = p32/sinp
         p22 = -p31/sinp
         p23 = dzero
         go to 170
160      p3131 = p31*p31
         p3131 = dmin1(done,p3131)
         sinp =  dsqrt(done-p3131)
         p21 = dzero
         p22 = p33/sinp
         p23 = -p32/sinp
170      q21 = p21
         q22 = p22
         q23 = p23
         p11 = p22*p33-p23*p32
         p12 = p23*p31-p21*p33
         p13 = p21*p32-p22*p31
         q11 = q22*q33-q23*q32
         q12 = q23*q31-q21*q33
         q13 = q21*q32-q22*q31
         acx(iab) = (cx-ax)*p11+(cy-ay)*p12+(cz-az)*p13
         acy(iab) = (cx-ax)*p21+(cy-ay)*p22+(cz-az)*p23
         if (  dabs(acy(iab))-pt0001) 180,180,190
180      acy(iab) = dzero
190      continue
         acz(iab) = (cx-ax)*p31+(cy-ay)*p32+(cz-az)*p33
         acy2(iab) = acy(iab)*acy(iab)
c
         pmat(1,iab)=p11
         pmat(2,iab)=p12
         pmat(3,iab)=p13
         pmat(4,iab)=p21
         pmat(5,iab)=p22
         pmat(6,iab)=p23
         pmat(7,iab)=p31
         pmat(8,iab)=p32
         pmat(9,iab)=p33
c        qmat(1,iab)=q11
c        qmat(2,iab)=q12
c        qmat(3,iab)=q13
c        qmat(4,iab)=q21
c        qmat(5,iab)=q22
c        qmat(6,iab)=q23
c        qmat(7,iab)=q31
c        qmat(8,iab)=q32
c        qmat(9,iab)=q33
c
         ind=0
         do 340 i = 1,nga
            ga = ag(i)
            csai = csa(i)
            cpai = cpa(i)
            do 330 j = 1,ngb
               ind = ind+1
               gb = bg(j)
               gab = ga+gb
               gp(iab,ind) = gab
               eab = done/gab
               ep(iab,ind) = eab
               gbeab = gb*eab
               app(iab,ind) = gbeab*rab(iab)
               bpp(iab,ind) = app(iab,ind)-rab(iab)
               qqq = ga*gbeab*rabsq
               if (qqq-sixty) 240,240,200
200            ismlp(iab,ind) = 2
               dp00p(iab,ind) = dzero
               if (jtype-3) 330,330,210
210            dp01p(iab,ind) = dzero
               conp(iab,ind) = dzero
               if (jtype-5) 220,220,230
220            bpp(iab,ind) = bpp(iab,ind)*gab
               go to 330
230            dp10p(iab,ind) = dzero
               dp11p = dzero
               go to 330
240            qq =  dexp(-qqq)*eab
               qqtest = cmaxa(i)*cmaxb(j)*qq
               if (qqtest-error1) 260,260,250
250            ismlp(iab,ind) = 0
               go to 290
260            if (qqtest-error2) 280,280,270
270            ismlp(iab,ind) = 1
               go to 290
280            ismlp(iab,ind) = 2
290            qqqq = pito52*qq
               dp00p(iab,ind) = qqqq*csai*csb(j)
               if (jtype-3) 330,330,300
300            dp01p(iab,ind) = qqqq*csai*cpb(j)
               if (jtype-5) 310,310,320
310            conp(iab,ind) = dp01p(iab,ind)*eab
               dp00p(iab,ind) = dp00p(iab,ind)*gab/dp01p(iab,ind)
               bpp(iab,ind) = bpp(iab,ind)*gab
               go to 330
320            dp10p(iab,ind) = qqqq*cpai*csb(j)
               dp11p = qqqq*cpai*cpb(j)
               conp(iab,ind) = dp11p
               dp00p(iab,ind) = dp00p(iab,ind)/dp11p
               dp01p(iab,ind) = dp01p(iab,ind)/dp11p
               dp10p(iab,ind) = dp10p(iab,ind)/dp11p
330         continue
340      continue
c ***
350   continue
c      write(6,*) ' qq4,ngangb,rab,gp,ep,dp00p,dp01p,dp10p'
c      write(6,1)  qq4,ngangb,rab,gp(1,1),ep(1,1),dp00p(1,1),
c     +dp01p(1,1),dp10p(1,1)
c      write(6,*) ' app,bpp,acx,acy,acz,acy2,cosg,sing,conp'
c      write(6,1)   app(1,1),bpp(1,1),acx,acy,acz,acy2,cosg,sing
c     +,conp(1,1)
c      write(6,*) ' ismlp,iptab'
c      write(6,1)   ismlp(1,1)
c1     format(1h ,5f7.3)
      return
      end
      subroutine stempl
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
      common/junk/pp(2*mxprms2),ind(mxprms2),indsp(256,4)
      do 10index=1,64
         indsp(index,1)=1
         indsp(index+64,1)=2
         indsp(index+128,1)=3
10    indsp(index+192,1)=4
      do 20ipt=1,193,64
         do 20index=ipt,ipt+15
            indsp(index,2)=1
            indsp(index+16,2)=2
            indsp(index+32,2)=3
20    indsp(index+48,2)=4
      do 30ipt=1,4
         do 30index=ipt,ipt+240,16
            indsp(index,3)=1
            indsp(index+4,3)=2
            indsp(index+8,3)=3
30    indsp(index+12,3)=4
      do 40index=1,253,4
         indsp(index,4)=1
         indsp(index+1,4)=2
         indsp(index+2,4)=3
40    indsp(index+3,4)=4
c ***
      return
      end
      subroutine mkngkl(nshell,ktype,kng,nsp,ngklmx,
     + mapsp,mapinv,ngvals,istng,lstng)
      dimension ktype(nshell),kng(nshell),mapsp(nshell),
     +          mapinv(nshell),ngvals(*),istng(*),lstng(*)
c ***
      nsp=0
      do 10 ii=1,nshell
         kit=ktype(ii)
         if(kit.eq.1 .or. kit.eq.2) then
            nsp=nsp+1
            mapsp(nsp)=ii
         endif
10    continue
      if(nsp.eq.0) return
      do 30 ii=1,nsp
         no=kng(mapsp(ii))
         do 20 jj=ii+1,nsp
            nn=kng(mapsp(jj))
            if(nn.gt.no) then
               mm=mapsp(jj)
               mapsp(jj)=mapsp(ii)
               mapsp(ii)=mm
               no=nn
            endif
20       continue
30    continue
      do 40 ii=1,nshell
40    mapinv(ii)=0
      do 50 ii=1,nsp
         mapinv(mapsp(ii))=ii
50    continue
c *** mapsp maps ordered s & sp shells into usual list
c *** mapinv maps usual list into ordered s & sp shell list
      do 60 ii=1,100
         ngvals(ii)=0
         istng(ii)=0
60    lstng(ii)=0
      ngp=kng(mapsp(1))
      istng(ngp)=1
      ngklmx=0
      do 70 ii=1,nsp
         ngi=kng(mapsp(ii))
         if(ngp.ne.ngi) then
            lstng(ngp)=ii-1
            istng(ngi)=ii
            ngp=ngi
         endif
         do 70 jj=1,nsp
            ngingj=ngi*kng(mapsp(jj))
            ngklmx=max(ngingj,ngklmx)
70    ngvals(ngingj)=ngvals(ngingj)+1
      lstng(ngi)=nsp
      return
      end
c ******************************************************
c ******************************************************
c             =   jk70v  =
c ******************************************************
c ******************************************************
_IF(ccpdft)
      subroutine jk70vs(zscftp,q,fock,fockb,exch,dens,densb,
     &     prefac,rdmat,coul,iso,gout,nshels,nschwz,
     &     fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine jk70vs(zscftp,q,fock,fockb,exch,dens,densb,
     &      prefac,rdmat,coul,iso,gout,nshels,nschwz)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c ***
c ***
_IF1(a)cvd$r noconcur
      dimension iso(nshels,*),q(*),coul(*)
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
c ***
INCLUDE(common/sortp)
      common/bufbb/ mapsp(mxshel),mapinv(mxshel),
     &          kvpppp(255),lvpppp(255),q4pppp(255),
     &          kvsppp(255),lvsppp(255),q4sppp(255),
     &          kvsspp(255),lvsspp(255),q4sspp(255),
     &          istng(mxprms2),lstng(mxprms2),ngvals(mxprms2)
c ***
      dimension mi(48),mj(48),mk(48),m0(48)
      dimension gout(*)
INCLUDE(common/confus)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/mapper)
INCLUDE(common/nshel)
INCLUDE(common/picon)
INCLUDE(common/pkfil)
INCLUDE(common/restar)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
INCLUDE(common/segm)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/timez)
INCLUDE(common/symtry)
      common/ttqout/nok,nnok
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c ***
c *** code for (sp sp / sp sp), (s sp / sp sp) and (s s / sp sp) /blocks
c *** code for other vectorised and also scalar evaluations
c *** can eventually be included here when i get a grip on things].
c ***
c *** allocate remainder of core.
c ***
      oskipp = .false.
      dlncutoff = dlog(cutoff)
      ifree = igmem_alloc_all(nfree)
c
      if(outv)
     *write(iwr,*) ' no. of words in jk70v = ',nfree
      call vclr(cpus,1,20)
c ***
c *** set up data for vectorisation of qout70 (sp sp/sp sp)
c ***
_IF(parallel)
c***   **MPP**
      iflop = iipsci()
c***   **MPP**
_ENDIF
      nok=0
      nnok=0
      call stempl
c ***
c *** form list of s & sp shell indices, order according to
c *** contraction level and form some other offset arrays
c *** for addressing purposes
c ***
      call mkngkl(nshels,ktype,kng,nsp,ngklmx,
     + mapsp,mapinv,ngvals,istng,lstng)
c ***
c *** istng(ng),lstng(ng) = range of shell indices with ng prims
c *** ngvals(ng) > 0 if can form ngk*ngl=ng, 0 otherwise
c *** ngklmx = maximum value of ngk*ngl
c *** mapsp maps from ordered shell list into usual one
c *** mapinv reverses this
c ***
      time = cpulft(1)
      tim0 = time
      tim1 = time
c ***
      ist0 = ist
      jst0 = jst
      kst0 = 1
      lst0 = 1
c ***
      do 380 iii = ist0,nsp
         ii=mapsp(iii)
c *** i and j must be sp shells.
         if(ktype(ii).ne.2) goto 380
         if(kad(ii))380,80,80
80       dt0 = time-tim0
         dt1 = time-tim1
         tim1 = time
      if(outv ) then
       if(odscf) then
            if(iter.le.0)
     *      write(iwr,1002)ii,jst0,kst0,lst0,nrec,dt1,dt0
            if(ologf) then
            write(text,1002)ii,jst0,kst0,lst0,nrec,dt1,dt0
            call sptchk(text)
            endif
       else
              write(iwr,1001)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
       endif
      endif
         do 130 it = 1,nt
            id = iso(ii,it)
            id=mapinv(id)
            if (id .gt. iii) go to 380
            m0(it) = id
130      continue
         j0 = jst0
         do 370 jjj = j0,iii
            jst0=1
            jj=mapsp(jjj)
            if(ktype(jj).ne.2.or.kad(jj).lt.0) goto 370
c ***
            itrij=iky(max(ii,jj)) + min(ii,jj)
            tolij=dlntol+prefac(itrij)
_IF1()c           if(tolij .le. -3.401e0) then
_IF1()c              if(odscf)intcut(1)=intcut(1)+1
_IF1()c              goto 370
_IF1()c           endif
c ***
            do 170 it = 1,nt
               jd = iso(jj,it)
               jd = mapinv(jd)
               if (jd .gt. iii) goto 370
               id = m0(it)
               if (id .ge. jd) go to 160
               nd = id
               id = jd
               jd = nd
160            if (id .eq. iii .and. jd .gt. jjj) go to 370
               mi(it) = id
               mj(it) = jd
170         continue
            k0 = kst0
c ***
c *** loop thru blocks of kl with same product contraction
c *** level.
c ***
            do 360 ngkngl=1,ngklmx
               if(ngvals(ngkngl).eq.0) goto 360
               mxpppp=min(255,nfree/(256*3+32+9*ngkngl))
               mxsppp=min(255,nfree/(3*64+8*ngkngl+36))
               mxsspp=min(255,nfree/(71+6*ngkngl))
c *** adjust mx.... to avoid bank conflicts
               if((mxpppp/2)*2.eq.mxpppp) mxpppp=mxpppp-1
          if(mxpppp.le.0) call caserr('insufficient  core for v1111')
               if((mxsppp/2)*2.eq.mxsppp) mxsppp=mxsppp-1
          if(mxsppp.le.0) call caserr('insufficient  core for v0111')
               if((mxsspp/2)*2.eq.mxsspp) mxsspp=mxsspp-1
          if(mxsspp.le.0) call caserr('insufficient  core for v0011')
               npppp=0
               nsppp=0
               nsspp=0
c *** loop thru all k,l that can form ngk*ngl=ngkngl
c *** this gives longest possible loop lengths
c *** note that need different symmetry constraints for
c *** (sp sp/sp sp) and (sp sp/s sp) shell blocks ...
c *** didn't in standard loop structure 'cos drove over
c *** fully canonical indices rather than constraining shell
c *** types for given indices.
               do 350 kkk=k0,nsp
                  kst0=1
                  kk=mapsp(kkk)
                  kkt=ktype(kk)
c *** restrict range if (sp sp/sp sp).
                  if(kkt.eq.2 .and. kkk.gt.iii) goto 350
                  ngk=kng(kk)
                  ngl=ngkngl/ngk
                  if(ngk*ngl.ne.ngkngl) goto 350
                  llo=istng(ngl)
                  if(llo.eq.0) goto 350
                  itrik = iky(max(ii,kk)) + min(ii,kk)
                  itrjk = iky(max(jj,kk)) + min(jj,kk)
                  if(kad(kk))350,180,180
180                  do 200 it = 1,nt
                     kd = iso(kk,it)
                     kd = mapinv(kd)
                     if(kkt.eq.2) then
                        if(kd.gt.iii) goto 350
                     else
                        id=mi(it)
                        jd=mj(it)
                        if(id.eq.iii .and. jd.eq.jjj .and. kd.gt.kkk)
     -                   goto 350
                     endif
200               mk(it) = kd
                  maxll = kkk
                  if (kkk .eq. iii) maxll = jjj
                  l0 = lst0
                  llo=max(llo,l0)
c ***
                  if(kkt.eq.2) then
c ***
                     if(llo.gt.maxll) goto 350
                     lhi=min(maxll,lstng(ngl))
                     do 260 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.2.or.kad(ll).lt.0) goto 260
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 260
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                        test = coul(itrij) + coul(itrkl)
                        if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=dlntol+ test +
     +                dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                      rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              nschwz = nschwz + 1
                              goto 260
                           endif
                        else
                         oskipp = test.lt.cutoff
                         if(oskipp) then
                           nschwz = nschwz + 1
                           go to 260
                         endif
                        endif
                        n4 = 0
                        do 250 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           if (ld .gt. iii) go to 260
                           kd = mk(it)
                           if(kd .gt. ld) go to 230
                           nd = kd
                           kd = ld
                           ld = nd
230                        id = mi(it)
                           jd = mj(it)
                           if (id.ne.iii .and. kd.ne.iii) goto 250
                           if (kd .lt. id) go to 240
                           if (kd.eq.id .and. ld.le.jd) goto 240
                           nd = id
                           id = kd
                           kd = nd
                           nd = jd
                           jd = ld
                           ld = nd
240                        if (jd .lt. jjj) go to 250
                           if (jd .gt. jjj) go to 260
                           if (kd .lt. kkk) go to 250
                           if (kd .gt. kkk) go to 260
                           if (ld .lt. lll) go to 250
                           if (ld .gt. lll) go to 260
                           n4 = n4+1
250                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        npppp=npppp+1
                        kvpppp(npppp)=max(kk,ll)
                        lvpppp(npppp)=min(kk,ll)
                        q4pppp(npppp)=qq4
_IF(ccpdft)
                        if(npppp.eq.mxpppp)call iv1111(npppp,mxpppp,
     +                   ngkngl,ii,jj,kvpppp,lvpppp,q4pppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(npppp.eq.mxpppp)call iv1111(npppp,mxpppp,
     +                   ngkngl,ii,jj,kvpppp,lvpppp,q4pppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
260                  continue
c ***
                  else
c ***
c *** loop for (sp sp / s sp)
                     lhi=lstng(ngl)
                     do 300 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.2.or.kad(ll).lt.0) goto 300
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 300
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                         test = coul(itrij) + coul(itrkl)
                         if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=dlntol+ test + 
     +                 dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                       rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              intcut(3)=intcut(3)+1
                             nschwz = nschwz + 1
                             goto 300
                           endif
                        else
                         oskipp = test.lt.cutoff
                         if(oskipp) then
                           nschwz = nschwz + 1
                           go to 300
                         endif
                        endif
                        n4 = 0
                        do 290 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           kd = mk(it)
                           id = mi(it)
                           jd = mj(it)
                           if(id.ne.iii) goto 290
                           if(jd.ne.jjj) goto 290
                           if(kd.ne.kkk) goto 290
                           if(ld.gt.lll) goto 300
                           if(ld.lt.lll) goto 290
                           n4 = n4+1
290                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        nsppp=nsppp+1
                        kvsppp(nsppp)=kk
                        lvsppp(nsppp)=ll
                        q4sppp(nsppp)=qq4
_IF(ccpdft)
                        if(nsppp.eq.mxsppp)call iv0111(nsppp,mxsppp,
     +                   ngkngl,ii,jj,kvsppp,lvsppp,q4sppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(nsppp.eq.mxsppp)call iv0111(nsppp,mxsppp,
     +                   ngkngl,ii,jj,kvsppp,lvsppp,q4sppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
300                  continue
c *** loop for (sp sp / s s)
                     lhi=min(lstng(ngl),kkk)
                     if(lhi.lt.llo) goto 350
                     do 340 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.1.or.kad(ll).lt.0) goto 340
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 340
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                        test = coul(itrij) + coul(itrkl)
                        if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=dlntol+ test + 
     +                 dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                       rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              intcut(3)=intcut(3)+1
                              nschwz = nschwz + 1
                              goto 340
                           endif
                         else
                          oskipp = test.lt.cutoff
                          if(oskipp) then
                            nschwz = nschwz + 1
                            go to 340
                          endif
                         endif
                        n4 = 0
                        do 330 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           kd = mk(it)
                           if(ld.gt.kd) then
                              nd=kd
                              kd=ld
                              ld=nd
                           endif
                           id = mi(it)
                           jd = mj(it)
                           if(id.ne.iii) goto 330
                           if(jd.ne.jjj) goto 330
                           if(kd.gt.kkk) goto 340
                           if(kd.lt.kkk) goto 330
                           if(ld.gt.lll) goto 340
                           if(ld.lt.lll) goto 330
                           n4 = n4+1
330                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        nsspp=nsspp+1
                        kvsspp(nsspp)=max(kk,ll)
                        lvsspp(nsspp)=min(kk,ll)
                        q4sspp(nsspp)=qq4
_IF(ccpdft)
                        if(nsspp.eq.mxsspp)call iv0011(nsspp,mxsspp,
     -                   ngkngl,ii,jj,kvsspp,lvsspp,q4sspp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(nsspp.eq.mxsspp)call iv0011(nsspp,mxsspp,
     -                   ngkngl,ii,jj,kvsspp,lvsspp,q4sspp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
340                  continue
                  endif
350            continue
_IF(ccpdft)
               if(npppp.ne.0) call iv1111(npppp,mxpppp,ngkngl, ii,
     -         jj,kvpppp,lvpppp,q4pppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
               if(nsppp.ne.0) call iv0111(nsppp,mxsppp,ngkngl, ii,
     -         jj,kvsppp,lvsppp,q4sppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
               if(nsspp.ne.0) call iv0011(nsspp,mxsspp,ngkngl, ii,
     -         jj,kvsspp,lvsspp,q4sspp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
               if(npppp.ne.0) call iv1111(npppp,mxpppp,ngkngl, ii,
     -         jj,kvpppp,lvpppp,q4pppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
               if(nsppp.ne.0) call iv0111(nsppp,mxsppp,ngkngl, ii,
     -         jj,kvsppp,lvsppp,q4sppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
               if(nsspp.ne.0) call iv0011(nsspp,mxsspp,ngkngl, ii,
     -         jj,kvsspp,lvsspp,q4sspp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
_ENDIF
360         continue
               call chkov(iii,jjj,q,fock,dens)
               if(omaxb.or.tim.gt.timlim)go to 390
370      continue
         time = cpulft(1)
380   continue
c ***
c ***
390   continue
      call gmem_free(ifree)
      if(outv) then
      write(iwr,*) ' v1111 timings '
      write(iwr,1000) (cpus(iqiq),iqiq=10,19)
c1000  format(1h ,'to 1808',f8.3,'  do 1801',f8.3,'  do 1901',f8.3/
c     +       1h ,'do 12..',f8.3,'  do 92..',f8.3,'  do 10..',f8.3/
c     +       1h ,'do 20..',f8.3,'  do 13..',f8.3,'  do 120.',f8.3/
1000  format(1h ,'mvml3  ',f8.3,'  sinfov ',f8.3,'  v1111  ',f8.3/
     +       1h ,'qout70 ',f8.3,'  spsinf ',f8.3,'  sp0111 ',f8.3/
     +       1h ,'spqout ',f8.3,'  sssinf ',f8.3,'  sp0011 ',f8.3/
     +       1h ,'ssqout ',f8.3)
      write(iwr,*) ' nok ',nok,' nnok ',nnok
      endif
      return
1001  format(i4,3i5,1x,i10,i9,f11.2,f9.2)
1002  format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
      end
_IF(ccpdft)
      subroutine jk70va(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels,
     +                  fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine jk70va(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c ***
c ***
_IF1(a)cvd$r noconcur
      dimension iso(nshels,*),q(*),gout(*)
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
c ***
      common/bufbb/ mapsp(mxshel),mapinv(mxshel),
     &          kvpppp(255),lvpppp(255),q4pppp(255),
     &          kvsppp(255),lvsppp(255),q4sppp(255),
     &          kvsspp(255),lvsspp(255),q4sspp(255),
     &          istng(mxprms2),lstng(mxprms2),ngvals(mxprms2)
c ***
      dimension mi(48),mj(48),mk(48),m0(48)
INCLUDE(common/confus)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/mapper)
INCLUDE(common/nshel)
INCLUDE(common/picon)
INCLUDE(common/pkfil)
INCLUDE(common/restar)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
INCLUDE(common/segm)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
INCLUDE(common/timez)
INCLUDE(common/symtry)
      common/ttqout/nok,nnok
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c ***
c *** code for (sp sp / sp sp), (s sp / sp sp) and (s s / sp sp) /blocks
c *** code for other vectorised and also scalar evaluations
c *** can eventually be included here when i get a grip on things].
c ***
c *** allocate remainder of core.
c ***
      oskipp = .false.
      ifree = igmem_alloc_all(nfree)
c
      if(outv)
     *write(iwr,*) ' no. of words in jk70v = ',nfree
      call vclr(cpus,1,20)
c ***
c *** set up data for vectorisation of qout70 (sp sp/sp sp)
c ***
_IF(parallel)
c***   **MPP**
      iflop = iipsci()
c***   **MPP**
_ENDIF
      nok=0
      nnok=0
      call stempl
c ***
c *** form list of s & sp shell indices, order according to
c *** contraction level and form some other offset arrays
c *** for addressing purposes
c ***
      call mkngkl(nshels,ktype,kng,nsp,ngklmx,
     + mapsp,mapinv,ngvals,istng,lstng)
c ***
c *** istng(ng),lstng(ng) = range of shell indices with ng prims
c *** ngvals(ng) > 0 if can form ngk*ngl=ng, 0 otherwise
c *** ngklmx = maximum value of ngk*ngl
c *** mapsp maps from ordered shell list into usual one
c *** mapinv reverses this
c ***
      time = cpulft(1)
      tim0 = time
      tim1 = time
c ***
      ist0 = ist
      jst0 = jst
      kst0 = 1
      lst0 = 1
c ***
      do 380 iii = ist0,nsp
         ii=mapsp(iii)
c *** i and j must be sp shells.
         if(ktype(ii).ne.2) goto 380
         if(kad(ii))380,80,80
80       dt0 = time-tim0
         dt1 = time-tim1
         tim1 = time
      if(outv ) then
       if(odscf) then
            if(iter.le.0)
     *      write(iwr,1002)ii,jst0,kst0,lst0,nrec,dt1,dt0
            if(ologf) then
            write(text,1002)ii,jst0,kst0,lst0,nrec,dt1,dt0
            call sptchk(text)
            endif
       else
              write(iwr,1001)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
       endif
      endif
         do 130 it = 1,nt
            id = iso(ii,it)
            id=mapinv(id)
            if (id .gt. iii) go to 380
            m0(it) = id
130      continue
         j0 = jst0
         do 370 jjj = j0,iii
            jst0=1
            jj=mapsp(jjj)
            if(ktype(jj).ne.2.or.kad(jj).lt.0) goto 370
c ***
            itrij=iky(max(ii,jj)) + min(ii,jj)
            tolij=dlntol+prefac(itrij)
_IF1()c           if(tolij .le. -3.401e0) then
_IF1()c              if(odscf)intcut(1)=intcut(1)+1
_IF1()c              goto 370
_IF1()c           endif
c ***
            do 170 it = 1,nt
               jd = iso(jj,it)
               jd = mapinv(jd)
               if (jd .gt. iii) goto 370
               id = m0(it)
               if (id .ge. jd) go to 160
               nd = id
               id = jd
               jd = nd
160            if (id .eq. iii .and. jd .gt. jjj) go to 370
               mi(it) = id
               mj(it) = jd
170         continue
            k0 = kst0
c ***
c *** loop thru blocks of kl with same product contraction
c *** level.
c ***
            do 360 ngkngl=1,ngklmx
               if(ngvals(ngkngl).eq.0) goto 360
               mxpppp=min(255,nfree/(256*3+32+9*ngkngl))
               mxsppp=min(255,nfree/(3*64+8*ngkngl+36))
               mxsspp=min(255,nfree/(71+6*ngkngl))
c *** adjust mx.... to avoid bank conflicts
               if((mxpppp/2)*2.eq.mxpppp) mxpppp=mxpppp-1
          if(mxpppp.le.0) call caserr('insufficient  core for v1111')
               if((mxsppp/2)*2.eq.mxsppp) mxsppp=mxsppp-1
          if(mxsppp.le.0) call caserr('insufficient  core for v0111')
               if((mxsspp/2)*2.eq.mxsspp) mxsspp=mxsspp-1
          if(mxsspp.le.0) call caserr('insufficient  core for v0011')
               npppp=0
               nsppp=0
               nsspp=0
c *** loop thru all k,l that can form ngk*ngl=ngkngl
c *** this gives longest possible loop lengths
c *** note that need different symmetry constraints for
c *** (sp sp/sp sp) and (sp sp/s sp) shell blocks ...
c *** didn't in standard loop structure 'cos drove over
c *** fully canonical indices rather than constraining shell
c *** types for given indices.
               do 350 kkk=k0,nsp
                  kst0=1
                  kk=mapsp(kkk)
                  kkt=ktype(kk)
c *** restrict range if (sp sp/sp sp).
                  if(kkt.eq.2 .and. kkk.gt.iii) goto 350
                  ngk=kng(kk)
                  ngl=ngkngl/ngk
                  if(ngk*ngl.ne.ngkngl) goto 350
                  llo=istng(ngl)
                  if(llo.eq.0) goto 350
                  itrik = iky(max(ii,kk)) + min(ii,kk)
                  itrjk = iky(max(jj,kk)) + min(jj,kk)
                  if(kad(kk))350,180,180
180                  do 200 it = 1,nt
                     kd = iso(kk,it)
                     kd = mapinv(kd)
                     if(kkt.eq.2) then
                        if(kd.gt.iii) goto 350
                     else
                        id=mi(it)
                        jd=mj(it)
                        if(id.eq.iii .and. jd.eq.jjj .and. kd.gt.kkk)
     -                   goto 350
                     endif
200               mk(it) = kd
                  maxll = kkk
                  if (kkk .eq. iii) maxll = jjj
                  l0 = lst0
                  llo=max(llo,l0)
c ***
                  if(kkt.eq.2) then
c ***
                     if(llo.gt.maxll) goto 350
                     lhi=min(maxll,lstng(ngl))
                     do 260 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.2.or.kad(ll).lt.0) goto 260
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 260
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                        tijkl=tolij+prefac(itrkl)
                        if(tijkl.le.0.0d0) then
                           if(odscf) intcut(2)=intcut(2)+1
                           goto 260
                        endif
                        if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=tijkl +
     +                 dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                       rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              intcut(3)=intcut(3)+1
                              goto 260
                           endif
                        endif
                        n4 = 0
                        do 250 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           if (ld .gt. iii) go to 260
                           kd = mk(it)
                           if(kd .gt. ld) go to 230
                           nd = kd
                           kd = ld
                           ld = nd
230                        id = mi(it)
                           jd = mj(it)
                           if (id.ne.iii .and. kd.ne.iii) goto 250
                           if (kd .lt. id) go to 240
                           if (kd.eq.id .and. ld.le.jd) goto 240
                           nd = id
                           id = kd
                           kd = nd
                           nd = jd
                           jd = ld
                           ld = nd
240                        if (jd .lt. jjj) go to 250
                           if (jd .gt. jjj) go to 260
                           if (kd .lt. kkk) go to 250
                           if (kd .gt. kkk) go to 260
                           if (ld .lt. lll) go to 250
                           if (ld .gt. lll) go to 260
                           n4 = n4+1
250                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        npppp=npppp+1
                        kvpppp(npppp)=max(kk,ll)
                        lvpppp(npppp)=min(kk,ll)
                        q4pppp(npppp)=qq4
_IF(ccpdft)
                        if(npppp.eq.mxpppp)call iv1111(npppp,mxpppp,
     -                   ngkngl,ii,jj,kvpppp,lvpppp,q4pppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(npppp.eq.mxpppp)call iv1111(npppp,mxpppp,
     -                   ngkngl,ii,jj,kvpppp,lvpppp,q4pppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
260                  continue
c ***
                  else
c ***
c *** loop for (sp sp / s sp)
                     lhi=lstng(ngl)
                     do 300 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.2.or.kad(ll).lt.0) goto 300
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 300
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                        tijkl=tolij+prefac(itrkl)
                        if(tijkl.le.0.0d0) then
                           if(odscf) intcut(2)=intcut(2)+1
                           goto 300
                        endif
                        if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=tijkl+
     +                 dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                       rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              intcut(3)=intcut(3)+1
                              goto 300
                           endif
                        endif
                        n4 = 0
                        do 290 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           kd = mk(it)
                           id = mi(it)
                           jd = mj(it)
                           if(id.ne.iii) goto 290
                           if(jd.ne.jjj) goto 290
                           if(kd.ne.kkk) goto 290
                           if(ld.gt.lll) goto 300
                           if(ld.lt.lll) goto 290
                           n4 = n4+1
290                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        nsppp=nsppp+1
                        kvsppp(nsppp)=kk
                        lvsppp(nsppp)=ll
                        q4sppp(nsppp)=qq4
_IF(ccpdft)
                        if(nsppp.eq.mxsppp)call iv0111(nsppp,mxsppp,
     +                   ngkngl,ii,jj,kvsppp,lvsppp,q4sppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(nsppp.eq.mxsppp)call iv0111(nsppp,mxsppp,
     -                   ngkngl,ii,jj,kvsppp,lvsppp,q4sppp, gout,
     +                   q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
300                  continue
c *** loop for (sp sp / s s)
                     lhi=min(lstng(ngl),kkk)
                     if(lhi.lt.llo) goto 350
                     do 340 lll=llo,lhi
                        lst0=1
                        ll=mapsp(lll)
                        if(ktype(ll).ne.1.or.kad(ll).lt.0) goto 340
_IF(parallel)
c***   **MPP**
                        if (oipsci()) go to 340
c***   **MPP**
_ENDIF
                        itrkl = iky(max(kk,ll)) + min(kk,ll)
                        tijkl=tolij+prefac(itrkl)
                        if(tijkl.le.0.0d0) then
                           if(odscf) intcut(2)=intcut(2)+1
                           goto 340
                        endif
                        if (odscf) then
                           mij = itrij
                           mik = itrik
                           mil = iky(max(ii,ll)) + min(ii,ll)
                           mjk = itrjk
                           mjl = iky(max(jj,ll)) + min(jj,ll)
                           mkl = itrkl
                           tijkl=tijkl + 
     +                dmax1(rdmat(mij),rdmat(mik),rdmat(mil),
     +                      rdmat(mjk),rdmat(mjl),rdmat(mkl))
                           if(tijkl.le.0.0d0) then
                              intcut(3)=intcut(3)+1
                              goto 340
                           endif
                        endif
                        n4 = 0
                        do 330 it = 1,nt
                           ld = iso(ll,it)
                           ld=mapinv(ld)
                           kd = mk(it)
                           if(ld.gt.kd) then
                              nd=kd
                              kd=ld
                              ld=nd
                           endif
                           id = mi(it)
                           jd = mj(it)
                           if(id.ne.iii) goto 330
                           if(jd.ne.jjj) goto 330
                           if(kd.gt.kkk) goto 340
                           if(kd.lt.kkk) goto 330
                           if(ld.gt.lll) goto 340
                           if(ld.lt.lll) goto 330
                           n4 = n4+1
330                     continue
                        q4 = dfloat(nt)/ dfloat(n4)
                        qq4 = q4
                        nsspp=nsspp+1
                        kvsspp(nsspp)=max(kk,ll)
                        lvsspp(nsspp)=min(kk,ll)
                        q4sspp(nsspp)=qq4
_IF(ccpdft)
                        if(nsspp.eq.mxsspp)call iv0011(nsspp,mxsspp,
     +                   ngkngl,ii,jj,kvsspp,lvsspp,q4sspp, 
     +                   gout,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree,
     +                   fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
                        if(nsspp.eq.mxsspp)call iv0011(nsspp,mxsspp,
     +                   ngkngl,ii,jj,kvsspp,lvsspp,q4sspp, 
     +                   gout,q,fock,fockb,exch,dens,densb,prefac,rdmat,
     +                   nfree,ifree)
_ENDIF
340                  continue
                  endif
350            continue
_IF(ccpdft)
               if(npppp.ne.0) call iv1111(npppp,mxpppp,ngkngl, ii,
     +         jj,kvpppp,lvpppp,q4pppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
               if(nsppp.ne.0) call iv0111(nsppp,mxsppp,ngkngl, ii,
     +         jj,kvsppp,lvsppp,q4sppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
               if(nsspp.ne.0) call iv0011(nsspp,mxsspp,ngkngl, ii,
     +         jj,kvsspp,lvsspp,q4sspp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree,
     +         fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
               if(npppp.ne.0) call iv1111(npppp,mxpppp,ngkngl, ii,
     -         jj,kvpppp,lvpppp,q4pppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
               if(nsppp.ne.0) call iv0111(nsppp,mxsppp,ngkngl, ii,
     -         jj,kvsppp,lvsppp,q4sppp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
               if(nsspp.ne.0) call iv0011(nsspp,mxsspp,ngkngl, ii,
     -         jj,kvsspp,lvsspp,q4sspp,gout,q,fock,fockb,exch,
     +         dens,densb,prefac,rdmat,nfree,ifree)
_ENDIF
360         continue
               call chkov(iii,jjj,q,fock,dens)
               if(omaxb.or.tim.gt.timlim)go to 390
370      continue
         time = cpulft(1)
380   continue
c ***
c ***
390   continue
      call gmem_free(ifree)
c
      if(outv) then
      write(iwr,*) ' v1111 timings '
      write(iwr,1000) (cpus(iqiq),iqiq=10,19)
c1000  format(1h ,'to 1808',f8.3,'  do 1801',f8.3,'  do 1901',f8.3/
c     +       1h ,'do 12..',f8.3,'  do 92..',f8.3,'  do 10..',f8.3/
c     +       1h ,'do 20..',f8.3,'  do 13..',f8.3,'  do 120.',f8.3/
1000  format(1h ,'mvml3  ',f8.3,'  sinfov ',f8.3,'  v1111  ',f8.3/
     +       1h ,'qout70 ',f8.3,'  spsinf ',f8.3,'  sp0111 ',f8.3/
     +       1h ,'spqout ',f8.3,'  sssinf ',f8.3,'  sp0011 ',f8.3/
     +       1h ,'ssqout ',f8.3)
      write(iwr,*) ' nok ',nok,' nnok ',nnok
      endif
      return
1001  format(i4,3i5,1x,i10,i9,f11.2,f9.2)
1002  format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
      end
      subroutine chkov(ii,jj,core,fock,dens)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/sortp)
INCLUDE(common/cslosc)
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/machin)
INCLUDE(common/ijlab)
INCLUDE(common/shlt)
INCLUDE(common/iofile)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpcut(7),iter
      character*10 charwall
c ***
      dimension core(*)
c ***
      if(omaxb) then
       call maxout
       return
      endif
      if(ochek) return
      ochek = .true.
      t4=cpulft(0)-safety
      if((timlst-t4)*2.5d0.lt.t4) then
        timlst = t4
      else
       osafe = safety.le.safe
       irest = 1
       ista = ii
       jsta = jj
       ksta = -1
       lsta = -1
       jsta = jj+1
       if (jsta .le. ii) go to 3
       jsta = 1
       ista = ii+1
       if (ista .gt. nshell) return
3      if(odscf) then
         if(icount.gt.1)  then
         ochek=.false.
         nrec=nrec+1
         endif
c        call dbuild(fock,dens)
       else
         call blocki
       endif
c
       isti = ista
       jsti = jsta
       ksti = ksta
       lsti = lsta
       lastb = iblkmp
       lastu = mfilep
       m2file = mfilep
       m2last(m2file)=iblkmp
       ist = ista
       jst = jsta
       kst = ksta
       lst = lsta
       ss   = cpulft(1)
       if(osafe.or.outv) write(iwr,5)ss,charwall(),ist,jst,kst,lst
        if(odscf) then
         call wrt3(fock,lentri,ibl171+iof171(3),idaf)
        else
         if(osafe.or.outv) then
           write(iwr,51)
           call filprn(m2file,m2blk,m2last,m2tape)
         endif
        endif
        itask(mtask)=irest
         call revise
         if(omaxb) then
         call maxout
         return
         endif
       t4=dumtim*1.5d0
       if(osafe) then
          write(iwr,400)
          write(iwr,9068)
          tim=timlim+0.5d0
       else
         safety=cpulft(0)
         if(safety.gt.t4)then
         safety=safety-dumtim
         else
         safety=safe
         endif
        timlst=cpulft(0)-safety
        if(ista.lt.nshell) then
         if(outv) then
          if(odscf) then
          if(iter.le.0) write(iwr,401)
          else
          write(iwr,401)
          endif
         endif
        endif
       endif
       call clredx
       endif
      return
c
 400  format(//
     *' *** insufficient time to complete integral evaluation'//)
 401  format(//
     * 1x,58('-')/
     * 1x,'ist',2x,'jst',2x,'kst',2x,'lst',7x,'nrec',
     * 3x,'intloc',5x,'del(t)',5x,'time'/
     * 1x,58('-')/)
 5     format(/1x,58('-')//
     * ' job dumped at ',f8.2,' seconds',a10,' wall'//
     * ' next batch ',4i5/)
 51   format(/
     * ' status of mainfile'/1x,18('-'))
 9068 format(/10x,27('*')/
     *10x,'*** warning ***'/
     *10x,'this job must be restarted'/
     *10x,27('*')/)
c
      end
c ******************************************************
c ******************************************************
c             =   iv1111  =
c ******************************************************
c ******************************************************
_IF(ccpdft)
      subroutine iv1111(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,gout,
     +                  q,fock,fockb,exch,dens,densb,prefac,
     +                  rdmat,nfree,ifree,
     +                  fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine iv1111(nkl,nklmax,ngkngl,i,j,kv,lv,qq4v,gout,
     +                  q,fock,fockb,exch,dens,densb,prefac,
     +                  rdmat,nfree,ifree)
_ENDIF
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character    (z),character    (x)
      implicit character    (y)
INCLUDE(common/sizes)
c ***
      dimension q(*)
      dimension kv(*),lv(*),qq4v(*),gout(*)
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
c ***
c ***
c *** data initialisation of pointers into core for the
c *** various arrays ... assume integer words same length as reals.
c ***
INCLUDE(common/confus)
INCLUDE(common/cslosc)
INCLUDE(common/flip70)
INCLUDE(common/geom)
INCLUDE(common/iofile)
INCLUDE(common/maxc)
INCLUDE(common/pkfil)
INCLUDE(common/shlg70)
INCLUDE(common/shllfo)
INCLUDE(common/shlnos)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
      common/type  /itype,jtype
      l2 = nx
      jtype=6
      itype=16
c *** this is just to unscrew saved variables in sinfo
      call sinset
c ***
      if(i.ge.j) then
         iish=i
         jjsh=j
         knew=i
         lnew=j
      else
         iish=j
         jjsh=i
         knew=j
         lnew=i
      endif
      ijsh=iish*4096+jjsh
      if(nkl.lt.7) then
         do 10 iab=1,nkl
            inew=kv(iab)
            jnew=lv(iab)
            kshell=kv(iab)
            lshell=lv(iab)
            klsh=kshell*4096+lshell
            if(ijsh.ge.klsh) then
               ishell=iish
               jshell=jjsh
               ib(1,1)=3
               ib(2,1)=4
               ib(3,1)=1
               ib(4,1)=2
            else
               ishell=kshell
               jshell=lshell
               kshell=iish
               lshell=jjsh
               ib(1,1)=1
               ib(2,1)=2
               ib(3,1)=3
               ib(4,1)=4
            endif
            call vclr(gout,1,256)
            qq4=qq4v(iab)
            call sinfo
            call sp1111(gout)
            if (odscf) then
              if (zscftp.eq.'uhf')then
_IF(ccpdft)
              call dir_build_uhf70(fock,fockb,
     +        dens,densb,gout,
     +        fac1,fac2, facex, ocoul, oexch)
_ELSE
              call dir_build_uhf70(fock,fockb,
     +        dens,densb,gout)
_ENDIF
              else if (zscftp.eq.'gvb') then
              if(nsheld.le.1) then
              call dir_build_open_70(fock,exch,dens,
     +                               gout)
              else
              call dir_build_open2_70(l2,fock,exch,
     +                               dens,gout)
              endif
              else
_IF(ccpdft)
_IF(cray)
              call qoutd70(fock,dens,gout,
     +        fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
              call dbuild70(fock,dens,gout,
     +        fac1, fac2, facex, ocoul, oexch)
_ENDIF
_ELSE
_IF(cray)
              call qoutd70(fock,dens,gout)
_ELSE
              call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
              endif
            else
              call qout70(gout)
            endif
10       continue
         nkl=0
         return
      endif
c ***
      nnkl=nkl
      if((nnkl/2)*2.eq.nnkl) nnkl=nnkl+1
      nnkl=min(nnkl,nklmax)
      nnklng=nnkl*ngkngl
c      write(6,*) ' nkl,nnkl,nnklng ',nkl,nnkl,nnklng
      igout=ifree
      irab=igout+nnkl*256
      ipmat=irab+nnkl
      iqmat=ipmat+nnkl*9
      igp=iqmat+nnkl*9
      iep=igp+nnklng
      idp00p=iep+nnklng
      idp01p=idp00p+nnklng
      idp10p=idp01p+nnklng
      iapp=idp10p+nnklng
      ibpp=iapp+nnklng
      iconp=ibpp+nnklng
      iismlp=iconp+nnklng
      iacx=iismlp+nnklng
      iacy=iacx+nnkl
      iacz=iacy+nnkl
      iacy2=iacz+nnkl
      icosg=iacy2+nnkl
      ising=icosg+nnkl
      iw1=ising+nnkl
      iaqz=iw1+nnkl
      isinp=iaqz+nnkl
      icosp=isinp+nnkl
      impaab=icosp+nnkl
      iisave=impaab+nnkl
      iqperp=iisave+nnkl
      ih=iqperp+nnkl
      ig=ih+nnkl*256
c *** arrays below overlap g(*,71...)
      ieab=ig+nnkl*70
      idp00=ieab+nnkl
      idp01=idp00+nnkl
      idp10=idp01+nnkl
      iap=idp10+nnkl
      ibp=iap+nnkl
      igpw=ibp+nnkl
      iconpw=igpw+nnkl
      iismlw=iconpw+nnkl
      ipqab=iismlw+nnkl
      igabcd=ipqab+nnkl
      ip=igabcd+nnkl
      if0=ip+nnkl
      if1=if0+nnkl
      if2=if1+nnkl
      if3=if2+nnkl
      if4=if3+nnkl
c ***
c *** note that g(*,100 ...) are used internally
c ***
      iused= ig+nnkl*256-ifree
      if(iused.gt.nfree) call caserr(' core error in iv1111.')
c ***
      c1=cpulft(1)
      call sinfov(jtype,knew,lnew,kv,lv,qq4v,ngkngl,ngc,cgg,csc,cpc,
     - ngd,dg,csd,cpd, q(irab),rcd,q(ipmat),q(iqmat), q(igp),q(iep)
     - ,q(idp00p),q(idp01p),q(idp10p),q(iapp),q(ibpp), q(iacx),q(iacy)
     - ,q(iacz),q(iacy2),q(icosg),q(ising), q(iconp), cmaxc,cmaxd,
     - q(iismlp),error1,error2, nkl,nnkl)
      c2=cpulft(1)
      cpus(11)=cpus(11)+c2-c1
c ***
      call vclr(q(igout),1,256*nnkl)
      call v1111(ngkngl,ngc,cgg,csc,cpc,ngd,dg,csd,cpd, q(irab),rcd,
     - q(ipmat),q(iqmat), q(igp),q(iep),q(idp00p),q(idp01p),q(idp10p)
     - ,q(iapp),q(ibpp), q(iacx),q(iacy),q(iacz),q(iacy2),q(icosg)
     - , q(ising), q(iconp), cmaxc,cmaxd,q(iismlp),error1,error2,
     -nkl, nnkl, q(ig),q(ih),q(impaab),q(iisave),q(iqperp),q(iaqz)
     -,q(iw1) , q(isinp),q(icosp),q(ieab),q(idp00),q(idp01),q(idp10)
     -, q(iap) ,q(ibp),q(igpw),q(iconpw),q(iismlw),q(ipqab), q(igabcd)
     -,q(ip) , q(if0),q(if1),q(if2),q(if3),q(if4),q(igout))
      c3=cpulft(1)
      cpus(12)=cpus(12)+c3-c2
      do 20 iab=1,nkl
         kshell=kv(iab)
         lshell=lv(iab)
         klsh=kshell*4096+lshell
         if(ijsh.ge.klsh) then
            ishell=iish
            jshell=jjsh
            ib(1,1)=3
            ib(2,1)=4
            ib(3,1)=1
            ib(4,1)=2
         else
            ishell=kshell
            jshell=lshell
            kshell=iish
            lshell=jjsh
            ib(1,1)=1
            ib(2,1)=2
            ib(3,1)=3
            ib(4,1)=4
         endif
         call dcopy(256,q(igout-1+iab),nnkl,gout,1)
         if (odscf) then
           if (zscftp.eq.'uhf')then
_IF(ccpdft)
            call dir_build_uhf70(fock,fockb,
     +      dens,densb,gout,
     +      fac1,fac2, facex, ocoul, oexch)
_ELSE
            call dir_build_uhf70(fock,fockb,
     +      dens,densb,gout)
_ENDIF
           else if (zscftp.eq.'gvb') then
           if(nsheld.le.1) then
           call dir_build_open_70(fock,exch,dens,
     +                            gout)
           else
           call dir_build_open2_70(l2,fock,exch,
     +                            dens,gout)
           endif
           else
_IF(ccpdft)
_IF(cray)
           call qoutd70(fock,dens,gout,
     +     fac1,fac2,facex,ocoul,oexch,odft)
_ELSE
           call dbuild70(fock,dens,gout,
     +     fac1, fac2, facex, ocoul, oexch)
_ENDIF
_ELSE
_IF(cray)
           call qoutd70(fock,dens,gout)
_ELSE
           call dbuild70(fock,dens,gout)
_ENDIF
_ENDIF
           endif
         else
           call qout70(gout)
         endif
20    continue
      c4=cpulft(1)
      cpus(13)=cpus(13)+c4-c3
      nkl=0
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine v1111(
c *** /miscg/
     +ngangb,
c *** /shllfo/
     +ngc,cg,csc,cpc,ngd,dg,csd,cpd,
c *** /geom/
     +rab,rcd,pmat,qmat,
c *** /pgeom/
     +gp,ep,dp00p,dp01p,dp10p,app,bpp,
c *** /qgeom/
     +acx,acy,acz,acy2,cosg,sing,
c *** /const/
     +conp,
c *** /maxc/
     +cmaxc,cmaxd,ismlp,error1,error2,
c ***
c ***
     +nab,nnab,g,h,mpaab,isave,qperp,aqz,w1,
     +sinp,cosp,eab,dp00,dp01,dp10,ap,bp,gpw,conpw,ismlpw,
     +pqab,gabcd,p,f0,f1,f2,f3,f4,gout)
c ***
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character    (z),character    (x)
      implicit character    (y)
_IF1(a)cvd$r altcode(32) vector
c ***
_IF1(cfu)      parameter (dzero=0.0,done=1.0)
_IFN1(cfu)      parameter (dzero=0.0d0,done=1.0d0)
      dimension cg(*),csc(*),cpc(*),dg(*),csd(*),
     +          cpd(*),rab(nab),pmat(9,nab),qmat(9,nab),
     +          gp(nnab,ngangb),ep(nnab,ngangb),
     +          dp00p(nnab,ngangb),dp01p(nnab,ngangb),
     +          dp10p(nnab,ngangb),app(nnab,ngangb),
     +          bpp(nnab,ngangb),acx(nab),
     +          acy(nab),acz(nab),acy2(nab),cosg(nab),sing(nab),
     +          conp(nnab,ngangb),cmaxc(*),cmaxd(*),
     +          ismlp(nnab,ngangb)
      dimension eab(nab),dp00(nab),dp01(nab),dp10(nab),ap(nab),
     +          bp(nab),gpw(nab),conpw(nab),ismlpw(nab),pqab(nab),
     +          p(nab),gabcd(nab),
     +          f0(nab),f1(nab),f2(nab),f3(nab),f4(nab)
c ***
      dimension g(nnab,256),h(nnab,256),mpaab(nab),
     +          isave(nab),qperp(nab),w1(nab),aqz(nab),
     +          sinp(nab),cosp(nab),gout(nnab,256)
c ***
      dimension kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
INCLUDE(common/auxvar)
_IF1()      common/ccount/nall,ngt
INCLUDE(common/confus)
      common/inttab/aa(333),ba(333),ca(333),abc1,ab(333),bb(333)
     -              ,cb(333),abc2,ac(333),bc(333),cc(333),abc3,ad(333)
     -              ,bd(333),cd(333),abc4,ae(333),be(333),ce(333)
     -              ,abc5,af(333),bf(333),cf(333),abc6
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c ***
c *** extrinsically vectorised version of sp1111
c ***
      do 860 k = 1,ngc
         gc = cg(k)
         do 850 l = 1,ngd
_IF1()            nall=nall+nab
            gd = dg(l)
            gcd = gc+gd
            ecd = done/gcd
            cq = gd*ecd*rcd
            dq = cq-rcd
            qqq = cq*dq*gcd
            if (qqq+sixty) 10,20,20
10          v = dzero
            go to 30
20          v =  dexp(qqq)*ecd
30          qqtest = cmaxc(k)*cmaxd(l)*v
            if (qqtest-error1) 50,50,40
40          ismlq = 0
            go to 70
50          if (qqtest-error2) 850,850,60
60          ismlq = 1
70          sc = csc(k)
            sd = csd(l)
            pc = cpc(k)
            pd = cpd(l)
            dq00 = sc*sd*v
            dq01 = sc*pd*v
            dq10 = pc*sd*v
            dq11 = pc*pd*v
            n120=0
            n920=0
            n1000=0
c
            do 150 iab=1,nab
               aqx = acx(iab)+sing(iab)*cq
               aqz(iab) = acz(iab)+cosg(iab)*cq
               qperp2 = aqx*aqx+acy2(iab)
               qperp(iab) = dsqrt(qperp2)
               if (qperp(iab)-tenm12) 90,90,80
80             cospp = -aqx/qperp(iab)
               sinpp = -acy(iab)/qperp(iab)
               go to 100
90             cospp = done
               sinpp = dzero
100            continue
c also form vector map around 120,920,1000 test.
               if (sinpp) 120,110,120
110            if (cospp) 140,120,130
120            n120=n120+1
               isave(iab)=n120
               mpaab(n120)=iab
               sinp(n120)=sinpp
               cosp(n120)=cospp
               goto 150
130            isave(iab)=nab-n920
               mpaab(nab-n920)=iab
               n920=n920+1
               goto 150
140            isave(iab)=1000
               n1000=n1000+1
150         continue
c
            m1000=n120
            do 160 iab=1,nab
               if(isave(iab).eq.1000) then
                  m1000=m1000+1
                  mpaab(m1000)=iab
               endif
160         continue
c have to reorder qperp,aqz,rab
_IF1(cfu)            call gather(nab,w1,qperp,mpaab)
_IFN1(cfu)            call dgthr(nab,qperp,w1,mpaab)
            do 170 i=1,nab
170         qperp(i)=w1(i)
_IF1(cfu)            call gather(nab,w1,aqz,mpaab)
_IFN1(cfu)            call dgthr(nab,aqz,w1,mpaab)
            do 180 i=1,nab
180         aqz(i)=w1(i)
_IF1(cfu)            call gather(nab,w1,rab,mpaab)
_IFN1(cfu)            call dgthr(nab,rab,w1,mpaab)
c ***
c *** p1-6,r1-9,w1-9,s1-14,t1-14,q1-6,v1-6,c1-6
c ***   1    7   16   25    39    53   59   65
c *** now overlay the first 70*nab elements of g
            call vclr(g,1,70*nnab)
c ***
c *** from here on use order split over 120,920,1000 test
c ***
c *** loop over primitives of k,l shells. have to gather
c *** some stuff first.
c ***
c *** some floating point exceptions are generated in the vectorised
c *** code. switch them off for all of do 1801
c ***
_IF1(c)            call sensefi(mode)
_IF1(c)            call clearfi
c ***
            do 310 kkll=1,ngangb
c *** fortran gather/scatter vectorises on xmp
               do 190 iab=1,nab
                  iiab=mpaab(iab)
                  eab(iab)=ep(iiab,kkll)
                  dp00(iab)=dp00p(iiab,kkll)
                  dp01(iab)=dp01p(iiab,kkll)
                  dp10(iab)=dp10p(iiab,kkll)
                  ap(iab)=app(iiab,kkll)
                  bp(iab)=bpp(iiab,kkll)
                  gpw(iab)=gp(iiab,kkll)
                  conpw(iab)=conp(iiab,kkll)
c                 ismlpw(iab)=ismlp(iiab,kkll)
190            continue
c              call gather(nab,eab,ep(1,kkll),mpaab)
c              call gather(nab,dp00,dp00p(1,kkll),mpaab)
c              call gather(nab,dp01,dp01p(1,kkll),mpaab)
c              call gather(nab,dp10,dp10p(1,kkll),mpaab)
c              call gather(nab,ap,app(1,kkll),mpaab)
c              call gather(nab,bp,bpp(1,kkll),mpaab)
c              call gather(nab,gpw,gp(1,kkll),mpaab)
c              call gather(nab,conpw,conp(1,kkll),mpaab)
c              call gather(nab,ismlpw,ismlp(1,kkll),mpaab)
c
_IF1(a)cvd$  novector
               do 200 iab=1,nab
                  pqab(iab) = aqz(iab)-ap(iab)
                  pqab2 = pqab(iab)*pqab(iab)
                  gabcd(iab) = done/(eab(iab)+ecd)
                  qperp2=qperp(iab)*qperp(iab)
                  p(iab) = gabcd(iab)*(qperp2+pqab2)
200            continue
_IF1()c              do 210 iab=1,nab
_IF1()c *** temporary zap for isml>=2
_IF1()c                 temp=cvmgt(0.0,conpw(iab),((ismlpw(iab)+ismlq).ge.2)
_IF1()c    -             )
_IF1()c                 f0(iab)= sqrt(0.7853981625/(p(iab)*(gpw(iab) +gcd)
_IF1()c    -             ))* temp
_IF1()c                 gtx = gabcd(iab)/p(iab)
_IF1()c                 f1(iab) = 0.5*f0(iab)*gtx
_IF1()c                 f2(iab) = 1.5*f1(iab)*gtx
_IF1()c                 f3(iab) = 2.5*f2(iab)*gtx
_IF1()c                 f4(iab) = 3.5*f3(iab)*gtx
_IF1()c210           continue
_IF1()c *** this interpolation formula is used less frequently
_IF1()c *** for large molecules
_IF1()c               do 230 iab=1,nab
_IF1()c                  isml=ismlq+ismlpw(iab)
_IF1()c                  if(isml.ge.2) goto 230
_IF1()c                  auxvar = var(isml+1)
_IF1()c                  if(p(iab).gt.auxvar) goto 230
_IF1()c                  ngt=ngt+1
_IF1()c                  q = conpw(iab)/sqrt(gpw(iab)+gcd)
_IF1()c                  gy = gabcd(iab)*q
_IF1()c                  ggy = gabcd(iab)*gy
_IF1()c                  gggy = gabcd(iab)*ggy
_IF1()c                  qq = p(iab)*12.5
_IF1()c                  n =   int(qq)
_IF1()c                  theta = qq-dfloat(n)
_IF1()c                  theta2 = theta*(theta-done)
_IF1()c                  theta3 = theta2*(theta-2.0)
_IF1()c                  theta4 = theta2*(theta+done)
_IF1()c                  f0(iab) = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)
_IF1()c     -             + theta4*ca(n+2))*q
_IF1()c                  f1(iab) = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)
_IF1()c     -             + theta4*cb(n+2))*gy
_IF1()c                  f2(iab) = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)
_IF1()c     -             + theta4*cc(n+2))*ggy
_IF1()c                  f3(iab) = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)
_IF1()c     -             + theta4*cd(n+2))*gggy
_IF1()c                  f4(iab) = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)
_IF1()c     -             + theta4*ce(n+2))*gggy* gabcd(iab)
_IF1()c230            continue
c ***
c *** use here fact that f0..f4 are allocated sequentially in core
c ***
               call compfs(nab,nnab,p,f0,ismlpw,4)
_IF1(a)cvd$  novector
               do 231 iab=1,nab
                  q = conpw(iab)/dsqrt(gpw(iab)+gcd)
                  gy = gabcd(iab)*q
                  ggy = gabcd(iab)*gy
                  gggy = gabcd(iab)*ggy
                  f0(iab)=f0(iab)*q
                  f1(iab)=f1(iab)*gy
                  f2(iab)=f2(iab)*ggy
                  f3(iab)=f3(iab)*gggy
                  f4(iab)=f4(iab)*gggy*gabcd(iab)
231            continue
c *** these are the loops that take all the time for high
c *** contraction levels .... now fully vectorised.
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 240 iab=1,nab
                  pqab2=pqab(iab)*pqab(iab)
c *** 100 = f1pqab, 101 = f2pqab, 102 = f3pqab, 103 = f1pqa2
c *** 104 = f2pqa2, 105 = f3pqa2, 106 = f3pqa3
                  g(iab,100) = f1(iab)*pqab(iab)
                  g(iab,101) = f2(iab)*pqab(iab)
                  g(iab,102) = f3(iab)*pqab(iab)
                  g(iab,103) = f1(iab)*pqab2
                  g(iab,104) = f2(iab)*pqab2
                  g(iab,105) = f3(iab)*pqab2
                  g(iab,106) = g(iab,105)*pqab(iab)
240            continue
c *** p1-6,r1-9,w1-9,s1-14,t1-14,q1-6,v1-6,c1-6
c ***   1    7   16   25    39    53   59   65
c *** now overlay the first 70*nab elements of g
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
               do 250 iab=1,nab
                  g(iab, 1) = g(iab, 1)+f0(iab) *dp00(iab)
                  g(iab, 2) = g(iab, 2)+f1(iab) *dp00(iab)
                  g(iab, 3) = g(iab, 3)+f2(iab) *dp00(iab)
                  g(iab, 4) = g(iab, 4)+g(iab,100)*dp00(iab)
                  g(iab, 5) = g(iab, 5)+g(iab,101)*dp00(iab)
                  g(iab, 6) = g(iab, 6)+g(iab,104)*dp00(iab)
250            continue
c ***
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 260 iab=1,nab
                  edp01 = eab(iab)*dp01(iab)
                  g(iab, 7) = g(iab, 7)+f1(iab) *edp01
                  g(iab, 8) = g(iab, 8)+f2(iab) *edp01
                  g(iab, 9) = g(iab, 9)+f3(iab) *edp01
                  g(iab,10) = g(iab,10)+g(iab,100) *edp01
                  g(iab,11) = g(iab,11)+g(iab,101) *edp01
                  g(iab,12) = g(iab,12)+g(iab,102) *edp01
                  g(iab,13) = g(iab,13)+g(iab,104) *edp01
                  g(iab,14) = g(iab,14)+g(iab,105) *edp01
                  g(iab,15) = g(iab,15)+g(iab,106) *edp01
260            continue
c ***
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 270 iab=1,nab
                  edp10 = eab(iab)*dp10(iab)
                  g(iab,16) = g(iab,16)+f1(iab) *edp10
                  g(iab,17) = g(iab,17)+f2(iab) *edp10
                  g(iab,18) = g(iab,18)+f3(iab) *edp10
                  g(iab,19) = g(iab,19)+g(iab,100) *edp10
                  g(iab,20) = g(iab,20)+g(iab,101) *edp10
                  g(iab,21) = g(iab,21)+g(iab,102) *edp10
                  g(iab,22) = g(iab,22)+g(iab,104) *edp10
                  g(iab,23) = g(iab,23)+g(iab,105) *edp10
                  g(iab,24) = g(iab,24)+g(iab,106) *edp10
270            continue
c ***
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 280 iab=1,nab
                  g(iab,25) = g(iab,25)+f0(iab) *eab(iab)
                  g(iab,26) = g(iab,26)+f1(iab) *eab(iab)
                  g(iab,27) = g(iab,27)+f2(iab) *eab(iab)
                  g(iab,28) = g(iab,28)+f3(iab) *eab(iab)
c *** no s5
                  g(iab,30) = g(iab,30)+g(iab,100)*eab(iab)
                  g(iab,31) = g(iab,31)+g(iab,101)*eab(iab)
                  g(iab,32) = g(iab,32)+g(iab,102)*eab(iab)
                  g(iab,33) = g(iab,33)+g(iab,103)*eab(iab)
                  g(iab,34) = g(iab,34)+g(iab,104)*eab(iab)
                  g(iab,35) = g(iab,35)+g(iab,105)*eab(iab)
                  g(iab,36) = g(iab,36)+g(iab,104)*pqab(iab)*eab(iab)
                  g(iab,37) = g(iab,37)+g(iab,106)*eab(iab)
                  g(iab,38) = g(iab,38)+g(iab,106)*pqab(iab)*eab(iab)
280            continue
c *** t1 not used
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 290 iab=1,nab
                  eab2 = eab(iab)*eab(iab)
                  f4pqab = f4(iab)*pqab(iab)
                  f4pqa2 = f4pqab*pqab(iab)
                  f4pqa3 = f4pqa2*pqab(iab)
                  g(iab,40) = g(iab,40)+f1(iab) *eab2
                  g(iab,41) = g(iab,41)+f2(iab) *eab2
                  g(iab,42) = g(iab,42)+f3(iab) *eab2
                  g(iab,43) = g(iab,43)+f4(iab) *eab2
                  g(iab,44) = g(iab,44)+g(iab,101)*eab2
                  g(iab,45) = g(iab,45)+g(iab,102)*eab2
                  g(iab,46) = g(iab,46)+f4pqab*eab2
                  g(iab,47) = g(iab,47)+g(iab,104)*eab2
                  g(iab,48) = g(iab,48)+g(iab,105)*eab2
                  g(iab,49) = g(iab,49)+f4pqa2*eab2
                  g(iab,50) = g(iab,50)+g(iab,106)*eab2
                  g(iab,51) = g(iab,51)+f4pqa3*eab2
                  g(iab,52) = g(iab,52)+f4pqa3*pqab(iab)*eab2
290            continue
c              if (w1(iab) .eq. dzero) go to end of loop
c ***
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
               do 300 iab=1,nab
                  bpdp01 = bp(iab)*dp01(iab)
                  g(iab,53) = g(iab,53)+f0(iab) *bpdp01
                  g(iab,54) = g(iab,54)+f1(iab) *bpdp01
                  g(iab,55) = g(iab,55)+f2(iab) *bpdp01
                  g(iab,56) = g(iab,56)+g(iab,100)*bpdp01
                  g(iab,57) = g(iab,57)+g(iab,101)*bpdp01
                  g(iab,58) = g(iab,58)+g(iab,104)*bpdp01
c ***
                  apdp10 = ap(iab)*dp10(iab)
                  g(iab,59) = g(iab,59)+f0(iab) *apdp10
                  g(iab,60) = g(iab,60)+f1(iab) *apdp10
                  g(iab,61) = g(iab,61)+f2(iab) *apdp10
                  g(iab,62) = g(iab,62)+g(iab,100)*apdp10
                  g(iab,63) = g(iab,63)+g(iab,101)*apdp10
                  g(iab,64) = g(iab,64)+g(iab,104)*apdp10
c ***
                  apbp = ap(iab)*bp(iab)
                  g(iab,65) = g(iab,65)+f0(iab) *apbp
                  g(iab,66) = g(iab,66)+f1(iab) *apbp
                  g(iab,67) = g(iab,67)+f2(iab) *apbp
                  g(iab,68) = g(iab,68)+g(iab,100)*apbp
                  g(iab,69) = g(iab,69)+g(iab,101)*apbp
                  g(iab,70) = g(iab,70)+g(iab,104)*apbp
300            continue
310         continue
c ***
_IF1(c)            if(mode.ne.0) call setfi
c ***
c      write(6,*) ' p r w s t q v c'
c      write(6,7181) (g(1,iio),iio=1,70)
c
c ***
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 320 iab=1,nab
               aqzz=aqz(iab)
               a1 = aqzz*g(iab,26)-g(iab,30)
               a2 = aqzz*g(iab,27)-g(iab,31)
               a3 = aqzz*g(iab,28)-g(iab,32)
               a4 = aqzz*g(iab,30)-g(iab,33)
               a5 = aqzz*g(iab,31)-g(iab,34)
               a6 = aqzz*g(iab,32)-g(iab,35)
               a8 = aqzz*g(iab,34)-g(iab,36)
               a9 = aqzz*g(iab,35)-g(iab,37)
               a10 = aqzz*g(iab,37)-g(iab,38)
               bqz = aqzz-w1(iab)
               b1 = bqz*g(iab,26)-g(iab,30)
               b2 = bqz*g(iab,27)-g(iab,31)
               b3 = bqz*g(iab,28)-g(iab,32)
               b4 = bqz*g(iab,30)-g(iab,33)
               b5 = bqz*g(iab,31)-g(iab,34)
               b6 = bqz*g(iab,32)-g(iab,35)
               b8 = bqz*g(iab,34)-g(iab,36)
               b9 = bqz*g(iab,35)-g(iab,37)
               b10 = bqz*g(iab,37)-g(iab,38)
               hecd = 0.5d0*ecd
               ecd2 = ecd*ecd
               hecd2 = 0.5d0*ecd2
               qecd = qperp(iab)*ecd
               hqecd = 0.5d0*qecd
               qecd2 = qperp(iab)*ecd2
               hqecd2 = 0.5d0*qecd2
               qperp2=qperp(iab)*qperp(iab)
               q2ecd = qperp2*ecd
               q3ecd = qperp(iab)*q2ecd
               q2ecd2 = qperp2*ecd2
               q3ecd2 = q2ecd2*qperp(iab)
               h(iab, 1) = g(iab, 1)
               h(iab, 2) = qecd*g(iab, 2)
               h(iab, 4) = -ecd*g(iab, 4)
               h(iab, 11) = hecd*(g(iab, 1)-ecd*g(iab, 2))
               h(iab, 6) = h(iab, 11)+q2ecd2*g(iab, 3)
               h(iab, 8) = -qecd2*g(iab, 5)
               h(iab, 16) = h(iab, 11)+ecd2*g(iab, 6)
               h(iab, 17) = -qperp(iab)*g(iab, 7)
               h(iab, 49) = g(iab,10)+g(iab,53)
               h(iab, 35) = hecd*g(iab, 7)
               h(iab, 18) = h(iab, 35)-q2ecd*g(iab, 8)
               h(iab, 20) = qecd*g(iab,11)
               h(iab, 50) = h(iab, 20)+qecd*g(iab,54)
               h(iab, 52) = h(iab, 35)-ecd*g(iab,13)-ecd*g(iab,56)
               h(iab, 39) = hqecd2*g(iab, 8)
               h(iab, 44) = -hecd2*g(iab,11)
               h(iab, 27) = h(iab, 39)-qperp(iab)*h(iab, 35)
               h(iab, 59) = h(iab, 44)+hecd*(h(iab, 49)-ecd*g(iab,
     -          54))
               h(iab, 24) = h(iab, 44)+q2ecd2*g(iab,12)
               h(iab, 56) = h(iab, 39)-qecd2*(g(iab,14)+g(iab,57))
               h(iab, 22) = h(iab, 27)+h(iab, 39)+h(iab, 39)-q3ecd2
     -          *g(iab, 9)
               h(iab, 32) = h(iab, 27)-qecd2*g(iab,14)
               h(iab, 54) = h(iab, 59)+q2ecd2*(g(iab,12)+g(iab,55)
     -          )
               h(iab, 64) = h(iab, 59)+h(iab, 44)+h(iab, 44)+ ecd2
     -          *(g(iab,15)+g(iab,58))
               h(iab, 65) = -qperp(iab)*g(iab,16)
               h(iab,193) = g(iab,19)+g(iab,59)
               h(iab,131) = hecd*g(iab,16)
               h(iab, 66) = h(iab,131)-q2ecd*g(iab,17)
               h(iab, 68) = qecd*g(iab,20)
               h(iab,194) = h(iab, 68)+qecd*g(iab,60)
               h(iab,196) = h(iab,131)-ecd*g(iab,22)-ecd*g(iab,62)
               h(iab,135) = hqecd2*g(iab,17)
               h(iab,140) = -hecd2*g(iab,20)
               h(iab, 75) = h(iab,135)-qperp(iab)*h(iab,131)
               h(iab,203) = h(iab,140)+hecd*(h(iab,193)-ecd*g(iab,
     -          60))
               h(iab, 72) = h(iab,140)+q2ecd2*g(iab,21)
               h(iab,200) = h(iab,135)-qecd2*(g(iab,23)+g(iab,63))
               h(iab, 70) = h(iab, 75)+h(iab,135)+h(iab,135)-q3ecd2
     -          *g(iab,18)
               h(iab, 80) = h(iab, 75)-qecd2*g(iab,23)
               h(iab,198) = h(iab,203)+q2ecd2*(g(iab,21)+g(iab,61)
     -          )
               h(iab,208) = h(iab,203)+h(iab,140)+h(iab,140)+ ecd2
     -          *(g(iab,24)+g(iab,64))
               h(iab,161) = 0.5d0*(g(iab,25)-g(iab,40))
               h(iab, 81) = h(iab,161)+qperp2*g(iab,41)
               h(iab,113) = -qperp(iab)*(g(iab,44)+b1)
               h(iab,209) = -qperp(iab)*(g(iab,44)+a1)
               h(iab,241) = h(iab,161)+g(iab,47)+a4+b4+g(iab,65)
               h(iab,162) = hqecd*(g(iab,26)-g(iab,41))
               h(iab, 82) = h(iab,162)-qecd*g(iab,41)+q3ecd*g(iab,
     -          42)
               temp = hecd*g(iab,44)-q2ecd*g(iab,45)
               h(iab,114) = temp+hecd*b1-q2ecd*b2
               h(iab,210) = temp+hecd*a1-q2ecd*a2
               h(iab,242) = h(iab,162)+qecd*(g(iab,48)+a5+b5+g(iab,
     -          66))
               h(iab, 99) = -hqecd*g(iab,41)
               h(iab,147) = h(iab, 99)
               h(iab,179) = hecd*(g(iab,44)+b1)
               h(iab,227) = hecd*(g(iab,44)+a1)
               h(iab,164) = hecd*(g(iab,44)-g(iab,30))
               h(iab, 84) = h(iab,164)-q2ecd*g(iab,45)
               temp = -hqecd*g(iab,41)+qecd*g(iab,48)
               h(iab,116) = temp+qecd*b5
               h(iab,212) = temp+qecd*a5
               h(iab,244) = h(iab,164)+ ecd*(g(iab,44)-g(iab,50)-a8
     -          -b8-g(iab,68))+hecd*(a1+b1)
               h(iab,103) = 0.25d0*ecd2*g(iab,41)-0.5d0*q2ecd2*g(iab,
     -          42)
               h(iab,151) = h(iab,103)
               h(iab,183) = hqecd2*(g(iab,45)+b2)
               h(iab,231) = hqecd2*(g(iab,45)+a2)
               h(iab,108) = hqecd2*g(iab,45)
               h(iab,156) = h(iab,108)
               h(iab,188) = hecd2*(0.5d0*g(iab,41)-g(iab,48)-b5)
               h(iab,236) = hecd2*(0.5d0*g(iab,41)-g(iab,48)-a5)
               hxxyy = 0.25d0*(ecd*(g(iab,25)-g(iab,40))- ecd2*(g(iab,
     -          26)-g(iab,41)))
               h(iab,171) = hxxyy+hecd2*g(iab,41)
               h(iab, 91) = hxxyy+0.5d0*(q2ecd*g(iab,41)-q2ecd2*g(iab,
     -          42))
               temp = hqecd*(ecd*g(iab,45)-g(iab,44))
               h(iab,123) = temp+hqecd*(ecd*b2-b1)
               h(iab,219) = temp+hqecd*(ecd*a2-a1)
               h(iab,251) = hxxyy+hecd*(g(iab,47)+a4+b4+g(iab,65))
     -          - hecd2*(g(iab,48)+a5+b5+g(iab,66))
               h(iab,166) = hxxyy+0.5d0*q2ecd2*(g(iab,27)-g(iab,42)
     -          )
               h(iab, 86) = hxxyy+(hecd2+0.5d0*q2ecd)*g(iab,41)+ q2ecd2
     -          *(-3.0d0*g(iab,42)+0.5d0*g(iab,27)+qperp2*g(iab, 43))
               h(iab,118) = 1.5d0*qecd2*(g(iab,45)+b2)- hqecd*(g(iab,
     -          44)+b1)-q3ecd2*(b3+g(iab,46))
               h(iab,214) = 1.5d0*qecd2*(g(iab,45)+a2)-hqecd*(g(iab,
     -          44)+a1)- q3ecd2*(a3+g(iab,46))
               h(iab,246) = hxxyy-hecd2*(qperp2*g(iab,42)+g(iab,48)
     -          +a5+b5)+ hecd*(g(iab,47)+a4+b4+g(iab,65)-ecd*g(iab,
     -          66))+ q2ecd2*(g(iab,49)+0.5d0*g(iab,27)+a6+b6+g(iab,
     -          67))
               h(iab,168) = hqecd2*(g(iab,45)-g(iab,31))
               h(iab, 88) = 1.5d0*qecd2*g(iab,45)-hqecd2*g(iab,31)
     -         - q3ecd2*g(iab,46)
               temp = hecd2*(0.5d0*g(iab,41)-g(iab,48))+q2ecd2*(g(iab,
     -          49)- 0.5d0*g(iab,42))
               h(iab,120) = temp-hecd2*b5+q2ecd2*b6
               h(iab,216) = temp-hecd2*a5+q2ecd2*a6
               h(iab,248) = qecd2*(1.5d0*g(iab,45)-g(iab,51)-a9-b9
     -         -g(iab,69))- hqecd2*(g(iab,31)-a2-b2)
               h(iab,176) = hxxyy+hecd2*(g(iab,34)-g(iab,48))
               h(iab, 96) = hxxyy-hecd2*(qperp2*g(iab,42)+g(iab,48)
     -          -g(iab,34))+ 0.5d0*q2ecd*g(iab,41)+q2ecd2*g(iab,49)
               h(iab,128) = qecd2*(1.5d0*g(iab,45)-g(iab,51)-b9)- hqecd
     -          *(g(iab,44)+b1)+hqecd2*b2
               h(iab,224) = qecd2*(1.5d0*g(iab,45)-g(iab,51)-a9)- hqecd
     -          *(g(iab,44)+a1)+hqecd2*a2
               h(iab,256) = hxxyy+ hecd2*(-3.0d0*(a5+b5)+g(iab,41)
     -         +g(iab,34)-g(iab,66))+ ecd2*(-3.0d0*g(iab,48)+g(iab,
     -          52)+a10+b10+g(iab,70))+ hecd*(g(iab,47)+a4+b4+g(iab,
     -          65))
320         continue
c ***
            if(n120.eq.0) goto 550
c ***
            iablo=1
            iabhi=n120
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 330 iab=iablo,iabhi
               v44 = cosp(iab)*cosp(iab)
               v77 = v44
               v47 = done-v44
               v74 = v47
               v54 = cosp(iab)*sinp(iab)
               v57 = -v54
               g(iab, 22 ) = v44*h(iab, 22) + v47*h(iab, 22+5)
               g(iab, 22+1) = v54*h(iab, 22) + v57*h(iab, 22+5)
               g(iab, 22+5) = v74*h(iab, 22) + v77*h(iab, 22+5)
               g(iab, 70 ) = v44*h(iab, 70) + v47*h(iab, 70+5)
               g(iab, 70+1) = v54*h(iab, 70) + v57*h(iab, 70+5)
               g(iab, 70+5) = v74*h(iab, 70) + v77*h(iab, 70+5)
               g(iab,118 ) = v44*h(iab,118) + v47*h(iab,118+5)
               g(iab,118+1) = v54*h(iab,118) + v57*h(iab,118+5)
               g(iab,118+5) = v74*h(iab,118) + v77*h(iab,118+5)
               g(iab,166 ) = v44*h(iab,166) + v47*h(iab,166+5)
               g(iab,166+1) = v54*h(iab,166) + v57*h(iab,166+5)
               g(iab,166+5) = v74*h(iab,166) + v77*h(iab,166+5)
               g(iab,214 ) = v44*h(iab,214) + v47*h(iab,214+5)
               g(iab,214+1) = v54*h(iab,214) + v57*h(iab,214+5)
               g(iab,214+5) = v74*h(iab,214) + v77*h(iab,214+5)
               g(iab, 86) = v44*h(iab, 86)+v47*h(iab, 91)
               g(iab, 87) = v54*h(iab, 86)+v57*h(iab, 91)
               g(iab, 91) = v74*h(iab, 86)+v77*h(iab, 91)
               g(iab,198) = v44*h(iab,198)+v47*h(iab,203)
               g(iab,199) = v54*h(iab,198)+v57*h(iab,203)
               g(iab,203) = v74*h(iab,198)+v77*h(iab,203)
               g(iab,246) = v44*h(iab,246)+v47*h(iab,251)
               g(iab,247) = v54*h(iab,246)+v57*h(iab,251)
               g(iab,251) = v74*h(iab,246)+v77*h(iab,251)
               g(iab, 54) = v44*h(iab, 54)+v47*h(iab, 59)
               g(iab, 55) = v54*h(iab, 54)+v57*h(iab, 59)
               g(iab, 59) = v74*h(iab, 54)+v77*h(iab, 59)
               g(iab, 6) = v44*h(iab, 6)+v47*h(iab, 11)
               g(iab, 7) = v54*h(iab, 6)+v57*h(iab, 11)
               g(iab, 11) = v74*h(iab, 6)+v77*h(iab, 11)
330         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 340 iab=iablo,iabhi
               v44 = cosp(iab)*cosp(iab)
               v47 = done-v44
               v54 = cosp(iab)*sinp(iab)
               v57 = -v54
               v45 = v57+v57
               v55 = v44-v47
               g(iab,102) = v45*h(iab,103)
               g(iab,103) = v55*h(iab,103)
               g(iab,134) = v45*h(iab,135)
               g(iab,135) = v55*h(iab,135)
               g(iab,150) = v45*h(iab,151)
               g(iab,151) = v55*h(iab,151)
               g(iab,182) = v45*h(iab,183)
               g(iab,183) = v55*h(iab,183)
               g(iab,230) = v45*h(iab,231)
               g(iab,231) = v55*h(iab,231)
               g(iab, 38) = v45*h(iab, 39)
               g(iab, 39) = v55*h(iab, 39)
340         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 350 iab=iablo,iabhi
               g(iab, 24 ) = cosp(iab)*h(iab, 24)
               g(iab, 24+4) = sinp(iab)*h(iab, 24)
               g(iab, 72 ) = cosp(iab)*h(iab, 72)
               g(iab, 72+4) = sinp(iab)*h(iab, 72)
               g(iab,120 ) = cosp(iab)*h(iab,120)
               g(iab,120+4) = sinp(iab)*h(iab,120)
               g(iab,168 ) = cosp(iab)*h(iab,168)
               g(iab,168+4) = sinp(iab)*h(iab,168)
               g(iab,216 ) = cosp(iab)*h(iab,216)
               g(iab,216+4) = sinp(iab)*h(iab,216)
               g(iab, 18 ) = cosp(iab)*h(iab, 18)
               g(iab, 18+1) = sinp(iab)*h(iab, 18)
               g(iab, 66 ) = cosp(iab)*h(iab, 66)
               g(iab, 66+1) = sinp(iab)*h(iab, 66)
               g(iab,114 ) = cosp(iab)*h(iab,114)
               g(iab,114+1) = sinp(iab)*h(iab,114)
               g(iab,162 ) = cosp(iab)*h(iab,162)
               g(iab,162+1) = sinp(iab)*h(iab,162)
               g(iab,210 ) = cosp(iab)*h(iab,210)
               g(iab,210+1) = sinp(iab)*h(iab,210)
               g(iab, 88) = cosp(iab)*h(iab, 88)
               g(iab, 92) = sinp(iab)*h(iab, 88)
               g(iab,200) = cosp(iab)*h(iab,200)
               g(iab,204) = sinp(iab)*h(iab,200)
               g(iab,248) = cosp(iab)*h(iab,248)
               g(iab,252) = sinp(iab)*h(iab,248)
               g(iab, 56) = cosp(iab)*h(iab, 56)
               g(iab, 60) = sinp(iab)*h(iab, 56)
               g(iab, 8) = cosp(iab)*h(iab, 8)
               g(iab, 12) = sinp(iab)*h(iab, 8)
               g(iab, 82) = cosp(iab)*h(iab, 82)
               g(iab, 83) = sinp(iab)*h(iab, 82)
               g(iab,194) = cosp(iab)*h(iab,194)
               g(iab,195) = sinp(iab)*h(iab,194)
               g(iab,242) = cosp(iab)*h(iab,242)
               g(iab,243) = sinp(iab)*h(iab,242)
               g(iab, 50) = cosp(iab)*h(iab, 50)
               g(iab, 51) = sinp(iab)*h(iab, 50)
               g(iab, 2) = cosp(iab)*h(iab, 2)
               g(iab, 3) = sinp(iab)*h(iab, 2)
350         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 360 iab=iablo,iabhi
               g(iab,104) = -sinp(iab)*h(iab,108)
               g(iab,108) = cosp(iab)*h(iab,108)
               g(iab,136) = -sinp(iab)*h(iab,140)
               g(iab,140) = cosp(iab)*h(iab,140)
               g(iab,152) = -sinp(iab)*h(iab,156)
               g(iab,156) = cosp(iab)*h(iab,156)
               g(iab,184) = -sinp(iab)*h(iab,188)
               g(iab,188) = cosp(iab)*h(iab,188)
               g(iab,232) = -sinp(iab)*h(iab,236)
               g(iab,236) = cosp(iab)*h(iab,236)
               g(iab, 40) = -sinp(iab)*h(iab, 44)
               g(iab, 44) = cosp(iab)*h(iab, 44)
               g(iab, 98) = -sinp(iab)*h(iab, 99)
               g(iab, 99) = cosp(iab)*h(iab, 99)
               g(iab,130) = -sinp(iab)*h(iab,131)
               g(iab,131) = cosp(iab)*h(iab,131)
               g(iab,146) = -sinp(iab)*h(iab,147)
               g(iab,147) = cosp(iab)*h(iab,147)
               g(iab,178) = -sinp(iab)*h(iab,179)
               g(iab,179) = cosp(iab)*h(iab,179)
               g(iab,226) = -sinp(iab)*h(iab,227)
               g(iab,227) = cosp(iab)*h(iab,227)
               g(iab, 34) = -sinp(iab)*h(iab, 35)
               g(iab, 35) = cosp(iab)*h(iab, 35)
360         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 370 iab=iablo,iabhi
               g(iab, 80) = h(iab, 80)
               g(iab, 96) = h(iab, 96)
               g(iab,107) = -g(iab,102)
               g(iab,112) = dzero
               g(iab,128) = h(iab,128)
               g(iab,139) = -g(iab,134)
               g(iab,144) = dzero
               g(iab,155) = -g(iab,150)
               g(iab,160) = dzero
               g(iab,176) = h(iab,176)
               g(iab,187) = -g(iab,182)
               g(iab,192) = dzero
               g(iab,235) = -g(iab,230)
               g(iab,240) = dzero
               g(iab, 43) = -g(iab, 38)
               g(iab, 48) = dzero
               g(iab, 68) = h(iab, 68)
               g(iab, 84) = h(iab, 84)
               g(iab,100) = dzero
               g(iab,116) = h(iab,116)
               g(iab,132) = dzero
               g(iab,148) = dzero
               g(iab,164) = h(iab,164)
               g(iab,180) = dzero
               g(iab,228) = dzero
               g(iab, 36) = dzero
               g(iab, 65) = h(iab, 65)
               g(iab, 81) = h(iab, 81)
               g(iab, 97) = dzero
               g(iab,113) = h(iab,113)
               g(iab,129) = dzero
               g(iab,145) = dzero
               g(iab,161) = h(iab,161)
               g(iab,177) = dzero
               g(iab,225) = dzero
               g(iab, 33) = dzero
370         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 380 iab=iablo,iabhi
               h(iab, 80) = cosp(iab)*g(iab, 80)
               h(iab, 96) = cosp(iab)*g(iab, 96)
               h(iab,112) = -sinp(iab)*g(iab,176)
               h(iab,128) = cosp(iab)*g(iab,128)
               h(iab,144) = sinp(iab)*g(iab, 80)
               h(iab,160) = sinp(iab)*g(iab, 96)
               h(iab,176) = cosp(iab)*g(iab,176)
               h(iab,192) = sinp(iab)*g(iab,128)
               h(iab, 68) = cosp(iab)*g(iab, 68)
               h(iab, 84) = cosp(iab)*g(iab, 84)
               h(iab,100) = -sinp(iab)*g(iab,164)
               h(iab,116) = cosp(iab)*g(iab,116)
               h(iab,132) = sinp(iab)*g(iab, 68)
               h(iab,148) = sinp(iab)*g(iab, 84)
               h(iab,164) = cosp(iab)*g(iab,164)
               h(iab,180) = sinp(iab)*g(iab,116)
               h(iab, 65) = cosp(iab)*g(iab, 65)
               h(iab, 81) = cosp(iab)*g(iab, 81)
               h(iab, 97) = -sinp(iab)*g(iab,161)
               h(iab,113) = cosp(iab)*g(iab,113)
               h(iab,129) = sinp(iab)*g(iab, 65)
               h(iab,145) = sinp(iab)*g(iab, 81)
               h(iab,161) = cosp(iab)*g(iab,161)
               h(iab,177) = sinp(iab)*g(iab,113)
380         continue
c ***
            do 400 kq1=70,118,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 390 iab=iablo,iabhi
                  h(iab,kq1 ) = cosp(iab)*g(iab,kq1 ) - sinp(iab)
     -             *g(iab,kq1+64)
                  h(iab,kq1+64) = sinp(iab)*g(iab,kq1 ) + cosp(iab)
     -             *g(iab,kq1+64)
                  h(iab,kq1+ 1) = cosp(iab)*g(iab,kq1+ 1) - sinp(iab)
     -             *g(iab,kq1+65)
                  h(iab,kq1+65) = sinp(iab)*g(iab,kq1+ 1) + cosp(iab)
     -             *g(iab,kq1+65)
                  h(iab,kq1+ 2) = cosp(iab)*g(iab,kq1+ 2) - sinp(iab)
     -             *g(iab,kq1+66)
                  h(iab,kq1+66) = sinp(iab)*g(iab,kq1+ 2) + cosp(iab)
     -             *g(iab,kq1+66)
                  h(iab,kq1+ 5) = cosp(iab)*g(iab,kq1+ 5) - sinp(iab)
     -             *g(iab,kq1+69)
                  h(iab,kq1+69) = sinp(iab)*g(iab,kq1+ 5) + cosp(iab)
     -             *g(iab,kq1+69)
                  h(iab,kq1+ 6) = cosp(iab)*g(iab,kq1+ 6) - sinp(iab)
     -             *g(iab,kq1+70)
                  h(iab,kq1+70) = sinp(iab)*g(iab,kq1+ 6) + cosp(iab)
     -             *g(iab,kq1+70)
c ***
                  h(iab,kq1- 4) = cosp(iab)*g(iab,kq1-4) - sinp(iab)
     -             *g(iab,kq1+60)
                  h(iab,kq1+60) = sinp(iab)*g(iab,kq1-4) + cosp(iab)
     -             *g(iab,kq1+60)
                  h(iab,kq1- 3) = cosp(iab)*g(iab,kq1-3) - sinp(iab)
     -             *g(iab,kq1+61)
                  h(iab,kq1+61) = sinp(iab)*g(iab,kq1-3) + cosp(iab)
     -             *g(iab,kq1+61)
390            continue
400         continue
            do 420 kq1=2,50,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 410 iab=iablo,iabhi
                  h(iab,kq1 ) = g(iab,kq1 )
                  h(iab,kq1+ 1) = g(iab,kq1+ 1)
                  h(iab,kq1+ 4) = g(iab,kq1+ 4)
                  h(iab,kq1+ 5) = g(iab,kq1+ 5)
                  h(iab,kq1+ 6) = g(iab,kq1+ 6)
                  h(iab,kq1+ 9) = g(iab,kq1+ 9)
410            continue
420         continue
            do 440 kq1=12,60,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 430 iab=iablo,iabhi
                  h(iab,kq1 ) = g(iab,kq1 )
                  h(iab,kq1+182) = g(iab,kq1+182)
                  h(iab,kq1+183) = g(iab,kq1+183)
                  h(iab,kq1+186) = g(iab,kq1+186)
                  h(iab,kq1+187) = g(iab,kq1+187)
                  h(iab,kq1+188) = g(iab,kq1+188)
430            continue
440         continue
            do 460 kq1=203,251,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 450 iab=iablo,iabhi
                  h(iab,kq1 ) = g(iab,kq1 )
                  h(iab,kq1+1) = g(iab,kq1+1)
450            continue
460         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 470 iab=iablo,iabhi
               h(iab, 48) = g(iab, 48)
               h(iab, 36) = g(iab, 36)
               h(iab,228) = g(iab,228)
               h(iab,240) = g(iab,240)
               h(iab,225) = g(iab,225)
               h(iab, 33) = g(iab, 33)
470         continue
            do 490 kq1=22,214,64
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 480 iab=iablo,iabhi
                  g(iab,kq1 ) = cosp(iab)*h(iab,kq1 ) - sinp(iab)
     -             * h(iab,kq1+16)
                  g(iab,kq1+16) = sinp(iab)*h(iab,kq1 ) + cosp(iab)
     -             * h(iab,kq1+16)
                  g(iab,kq1+ 1) = cosp(iab)*h(iab,kq1+ 1) - sinp(iab)
     -             * h(iab,kq1+17)
                  g(iab,kq1+17) = sinp(iab)*h(iab,kq1+ 1) + cosp(iab)
     -             * h(iab,kq1+17)
                  g(iab,kq1+ 2) = cosp(iab)*h(iab,kq1+ 2) - sinp(iab)
     -             * h(iab,kq1+18)
                  g(iab,kq1+18) = sinp(iab)*h(iab,kq1+ 2) + cosp(iab)
     -             * h(iab,kq1+18)
                  g(iab,kq1+ 5) = cosp(iab)*h(iab,kq1+ 5) - sinp(iab)
     -             * h(iab,kq1+21)
                  g(iab,kq1+21) = sinp(iab)*h(iab,kq1+ 5) + cosp(iab)
     -             * h(iab,kq1+21)
                  g(iab,kq1+ 6) = cosp(iab)*h(iab,kq1+ 6) - sinp(iab)
     -             * h(iab,kq1+22)
                  g(iab,kq1+22) = sinp(iab)*h(iab,kq1+ 6) + cosp(iab)
     -             * h(iab,kq1+22)
                  g(iab,kq1+10) = cosp(iab)*h(iab,kq1+10) - sinp(iab)
     -             * h(iab,kq1+26)
                  g(iab,kq1+26) = sinp(iab)*h(iab,kq1+10) + cosp(iab)
     -             * h(iab,kq1+26)
                  g(iab,kq1- 5) = cosp(iab)*h(iab,kq1-5) - sinp(iab)
     -             * h(iab,kq1+11)
                  g(iab,kq1+11) = sinp(iab)*h(iab,kq1-5) + cosp(iab)
     -             * h(iab,kq1+11)
                  g(iab,kq1- 4) = cosp(iab)*h(iab,kq1-4) - sinp(iab)
     -             * h(iab,kq1+12)
                  g(iab,kq1+12) = sinp(iab)*h(iab,kq1-4) + cosp(iab)
     -             * h(iab,kq1+12)
                  g(iab,kq1- 3) = cosp(iab)*h(iab,kq1-3) - sinp(iab)
     -             * h(iab,kq1+13)
                  g(iab,kq1+13) = sinp(iab)*h(iab,kq1-3) + cosp(iab)
     -             * h(iab,kq1+13)
                  g(iab,kq1- 2) = cosp(iab)*h(iab,kq1-2) - sinp(iab)
     -             * h(iab,kq1+14)
                  g(iab,kq1+14) = sinp(iab)*h(iab,kq1-2) + cosp(iab)
     -             * h(iab,kq1+14)
480            continue
490         continue
            do 520 kq1=49,177,64
               kkq1=kq1
               do 510 kkkq1=1,32
                  do 500 iab=iablo,iabhi
                     g(iab,kkq1)=h(iab,kkq1)
500               continue
                  kkq1=kkq1+1
510            continue
520         continue
            do 540 kq1=1,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 530 iab=iablo,iabhi
                  g(iab,kq1)=h(iab,kq1)
                  g(iab,kq1+240)=h(iab,kq1+240)
530            continue
540         continue
c
550         if(n920.eq.0) goto 600
c
            iablo=nab-n920+1
            iabhi=nab
            do 570 kkq1=1,10
               kq1=kq1off(kkq1)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 560 iab=iablo,iabhi
                  g(iab,kq1 ) = h(iab,kq1 )
                  g(iab,kq1+ 1) = h(iab,kq1+ 1)
                  g(iab,kq1+ 2) = dzero
                  g(iab,kq1+ 3) = h(iab,kq1+ 3)
                  g(iab,kq1+ 5) = h(iab,kq1+ 5)
                  g(iab,kq1+ 6) = dzero
                  g(iab,kq1+ 7) = h(iab,kq1+ 7)
                  g(iab,kq1+10) = h(iab,kq1+10)
                  g(iab,kq1+11) = dzero
                  g(iab,kq1+15) = h(iab,kq1+15)
560            continue
570         continue
c ***
            do 590 kkq1=1,6
               kq1=kq2off(kkq1)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 580 iab=iablo,iabhi
                  g(iab,kq1 ) = dzero
                  g(iab,kq1+ 1) = dzero
                  g(iab,kq1+ 2) = h(iab,kq1+ 2)
                  g(iab,kq1+ 3) = dzero
                  g(iab,kq1+ 5) = dzero
                  g(iab,kq1+ 6) = h(iab,kq1+ 6)
                  g(iab,kq1+ 7) = dzero
                  g(iab,kq1+10) = dzero
                  g(iab,kq1+11) = h(iab,kq1+11)
                  g(iab,kq1+15) = dzero
580            continue
590         continue
c
600         if(n1000.eq.0) goto 680
c
            iablo=n120+1
            iabhi=n120+n1000
            do 620 kkq1=1,6
               kq1=kq3off(kkq1)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 610 iab=iablo,iabhi
                  g(iab,kq1 ) = h(iab,kq1 )
                  g(iab,kq1+ 1) =-h(iab,kq1+ 1)
                  g(iab,kq1+ 2) = dzero
                  g(iab,kq1+ 3) = h(iab,kq1+ 3)
                  g(iab,kq1+ 5) = h(iab,kq1+ 5)
                  g(iab,kq1+ 6) = dzero
                  g(iab,kq1+ 7) =-h(iab,kq1+ 7)
                  g(iab,kq1+10) = h(iab,kq1+10)
                  g(iab,kq1+11) = dzero
                  g(iab,kq1+15) = h(iab,kq1+15)
610            continue
620         continue
            do 640 kkq1=1,4
               kq1=kq4off(kkq1)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 630 iab=iablo,iabhi
                  g(iab,kq1 ) = -h(iab,kq1 )
                  g(iab,kq1+ 1) = h(iab,kq1+ 1)
                  g(iab,kq1+ 2) = dzero
                  g(iab,kq1+ 3) = -h(iab,kq1+ 3)
                  g(iab,kq1+ 5) = -h(iab,kq1+ 5)
                  g(iab,kq1+ 6) = dzero
                  g(iab,kq1+ 7) = h(iab,kq1+ 7)
                  g(iab,kq1+10) = -h(iab,kq1+10)
                  g(iab,kq1+11) = dzero
                  g(iab,kq1+15) = -h(iab,kq1+15)
630            continue
640         continue
            do 660 kkq1=1,6
               kq1=kq5off(kkq1)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 650 iab=iablo,iabhi
                  g(iab,kq1 ) = dzero
                  g(iab,kq1+ 1) = dzero
                  g(iab,kq1+ 2) = h(iab,kq1+ 2)
                  g(iab,kq1+ 3) = dzero
                  g(iab,kq1+ 5) = dzero
                  g(iab,kq1+ 6) = -h(iab,kq1+ 6)
                  g(iab,kq1+ 7) = dzero
                  g(iab,kq1+10) = dzero
                  g(iab,kq1+11) = h(iab,kq1+11)
                  g(iab,kq1+15) = dzero
650            continue
660         continue
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 670 iab=iablo,iabhi
               g(iab,99)=-g(iab,99)
               g(iab,108)=-g(iab,108)
               g(iab,147)=-g(iab,147)
               g(iab,156)=-g(iab,156)
               g(iab,103)=-g(iab,103)
               g(iab,151)=-g(iab,151)
670         continue
c
680         continue
c ***
c *** now have to recover from the 120,920,1000 gather
c *** nnab should be set to avoid bank conflicts
c ***
            call dcopy(nnab*256,g,1,h,1)
            do 690 iab=1,nab
               iiab=mpaab(iab)
               call dcopy(256,h(iab,1),nnab,g(iiab,1),nnab)
690         continue
c ***
            do 710 kq1=2,242,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 700 iab=1,nab
                  g(iab,kq1+ 3) = g(iab,kq1 )
                  g(iab,kq1+ 7) = g(iab,kq1+ 1)
                  g(iab,kq1+11) = g(iab,kq1+ 2)
                  g(iab,kq1+ 8) = g(iab,kq1+ 5)
                  g(iab,kq1+12) = g(iab,kq1+ 6)
                  g(iab,kq1+13) = g(iab,kq1+10)
700            continue
710         continue
c
            if (rcd.eq.0.0d0) goto 780
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 720 iab=1,nab
c *** r13
               qperp(iab) = cq*sing(iab)
c *** r33
               w1(iab) = cq*cosg(iab)
720         continue
c ***
            kkq1=1
            do 740 jq1=1,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 730 iab=1,nab
                  g(iab,kkq1+4) = qperp(iab)*g(iab,kkq1) + g(iab,kkq1
     -             +4)
                  g(iab,kkq1+12)= w1(iab)*g(iab,kkq1) + g(iab,kkq1
     -             +12)
                  g(iab,kkq1+5) = qperp(iab)*g(iab,kkq1+1) + g(iab,
     -             kkq1+5)
                  g(iab,kkq1+13)= w1(iab)*g(iab,kkq1+1) + g(iab,kkq1
     -             +13)
                  g(iab,kkq1+6) = qperp(iab)*g(iab,kkq1+2) + g(iab,
     -             kkq1+6)
                  g(iab,kkq1+14)= w1(iab)*g(iab,kkq1+2) + g(iab,kkq1
     -             +14)
                  g(iab,kkq1+7) = qperp(iab)*g(iab,kkq1+3) + g(iab,
     -             kkq1+7)
                  g(iab,kkq1+15)= w1(iab)*g(iab,kkq1+3) + g(iab,kkq1
     -             +15)
730            continue
740         kkq1=kkq1+16
            do 750 iab=1,nab
c *** r14
               qperp(iab) = dq*sing(iab)
c *** r34
               w1(iab) = dq*cosg(iab)
750         continue
            do 770 kq1=1,241,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 760 iab=1,nab
                  g(iab,kq1+1) = qperp(iab)*g(iab,kq1) + g(iab,kq1
     -             +1)
                  g(iab,kq1+3) = w1(iab)*g(iab,kq1) + g(iab,kq1+3)
                  g(iab,kq1+5) = qperp(iab)*g(iab,kq1+4) + g(iab,kq1
     -             +5)
                  g(iab,kq1+7) = w1(iab)*g(iab,kq1+4) + g(iab,kq1
     -             +7)
                  g(iab,kq1+9) = qperp(iab)*g(iab,kq1+8) + g(iab,kq1
     -             +9)
                  g(iab,kq1+11) = w1(iab)*g(iab,kq1+8) + g(iab,kq1
     -             +11)
                  g(iab,kq1+13) = qperp(iab)*g(iab,kq1+12) + g(iab,
     -             kq1+13)
                  g(iab,kq1+15) = w1(iab)*g(iab,kq1+12) + g(iab,kq1
     -            +15)
760            continue
770         continue
c
c *** extra work here can be avoided in this version
c
_IF1(a)cvd$  nodepchk
780         do 800 kq1=1,256
               do 790 iab=1,nab
                  gout(iab,kq1)=gout(iab,kq1)+dq11*g(iab,kq1)
790            continue
800         continue
            dq01x=dq01-dq11
_IF1(a)cvd$  nodepchk
            do 820 kq1=2,242,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 810 iab=1,nab
                  gout(iab,kq1 ) = dq01x*g(iab,kq1 ) + gout(iab,
     -            kq1)
                  gout(iab,kq1+1) = dq01x*g(iab,kq1+1) + gout(iab,
     -             kq1+1)
                  gout(iab,kq1+2) = dq01x*g(iab,kq1+2) + gout(iab,
     -             kq1+2)
810            continue
820         continue
            dq10x=dq10-dq11
            dq00x=dq00-dq11
_IF1(a)cvd$  nodepchk
            do 840 kq1=1,241,16
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
               do 830 iab=1,nab
                  gout(iab,kq1 ) = dq00x*g(iab,kq1 ) + gout(iab,
     -            kq1)
                  gout(iab,kq1+ 4) = dq10x*g(iab,kq1+ 4) + gout(iab,
     -             kq1+ 4)
                  gout(iab,kq1+ 8) = dq10x*g(iab,kq1+ 8) + gout(iab,
     -             kq1+ 8)
                  gout(iab,kq1+12) = dq10x*g(iab,kq1+12) + gout(iab,
     -             kq1+12)
830            continue
840         continue
c
850      continue
c
860   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order (  1),(  2),(  3),...(  5),(  6),...( 17),( 18),.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case ( 86)
c
c
c
_IF1(a)cvd$  nodepchk
      do 930 iab=1,nab
         p11=pmat(1,iab)
         p12=pmat(2,iab)
         p13=pmat(3,iab)
         p21=pmat(4,iab)
         p22=pmat(5,iab)
         p23=pmat(6,iab)
         p31=pmat(7,iab)
         p32=pmat(8,iab)
         p33=pmat(9,iab)
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
         do 870 ig=65,128
            t1=gout(iab,ig )
            t2=gout(iab,ig+ 64)
            t3=gout(iab,ig+128)
            gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
            gout(iab,ig+ 64) = p12*t1 + p22*t2 + p32*t3
            gout(iab,ig+128) = p13*t1 + p23*t2 + p33*t3
870      continue
         ind=17
         do 890 i=1,4
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 880 ig=ind,ind+15
               t1=gout(iab,ig )
               t2=gout(iab,ig+ 16)
               t3=gout(iab,ig+ 32)
               gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
               gout(iab,ig+ 16) = p12*t1 + p22*t2 + p32*t3
               gout(iab,ig+ 32) = p13*t1 + p23*t2 + p33*t3
880         continue
890      ind=ind+64
         ind=5
         do 910 i=1,4
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
            do 900 ig=ind,ind+240,16
               t1=gout(iab,ig )
               t2=gout(iab,ig+ 4)
               t3=gout(iab,ig+ 8)
               gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
               gout(iab,ig+ 4) = p12*t1 + p22*t2 + p32*t3
               gout(iab,ig+ 8) = p13*t1 + p23*t2 + p33*t3
900         continue
910      ind=ind+1
_IF1(ct)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
         do 920 ig=2,254,4
            t1=gout(iab,ig )
            t2=gout(iab,ig+ 1)
            t3=gout(iab,ig+ 2)
            gout(iab,ig ) = p11*t1 + p21*t2 + p31*t3
            gout(iab,ig+ 1) = p12*t1 + p22*t2 + p32*t3
            gout(iab,ig+ 2) = p13*t1 + p23*t2 + p33*t3
920      continue
930   continue
c
      return
      end
      subroutine setfm
      implicit REAL  (a-h,o-z)
      common/values/fm(2001,9),rdelta,delta,delo2
      dimension t(2001),et(2001)
      delta=40.0d0/2000.0d0
      delo2=delta*0.5d0
      rdelta=1.0d0/delta
      maxm=8
      do 10 i=1,2001
          tt=delta*dfloat(i-1)
          et(i)=dexp(-tt)
          t(i)=2.0d0*tt
10        fm(i,maxm+1)=0.0d0
      do 20 i=200,maxm,-1
          rr=1.0d0/dfloat(2*i+1)
          do 20 ii=1,2001
20            fm(ii,maxm+1)=(et(ii)+t(ii)*fm(ii,maxm+1))*rr
      do 30 i=maxm,1,-1
          rr=1.0d0/dfloat(2*i-1)
          do 30 ii=1,2001
30            fm(ii,i)=(et(ii)+t(ii)*fm(ii,i+1))*rr
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine compfs(nval,nnval,t,f,ind,mhigh)
_IF1(a)cvd$r novector
      implicit REAL  (a-h,o-z)
      dimension t(nval),f(nnval,5),ind(nval)
      common/values/fm(2001,9),rdelta,delta,delo2
_IF1(cfu)      parameter(fac4=5.81586419828372,
_IF1(cfu)     1          fac3=1.66167548522392,
_IF1(cfu)     2          fac2=0.66467019408957,
_IF1(cfu)     3          fac1=0.44311346272638,
_IF1(cfu)     4          fac0=0.88622692545276,
_IF1(cfu)     5          rhalf=0.5,rthird=0.333333333333333,rquart=0.25,
_IF1(cfu)     6          r2o7=0.285714285714286,r2o5=0.4,
_IF1(cfu)     7          r2o3=0.666666666666667)
_IFN1(cfu)      parameter(fac4=5.81586419828372d0,
_IFN1(cfu)     1          fac3=1.66167548522392d0,
_IFN1(cfu)     2          fac2=0.66467019408957d0,
_IFN1(cfu)     3          fac1=0.44311346272638d0,
_IFN1(cfu)     4          fac0=0.88622692545276d0,
_IFN1(cfu)     5          rhalf=0.5d0,rthird=0.333333333333333d0,rquart=0.25d0,
_IFN1(cfu)     6          r2o7=0.285714285714286d0,r2o5=0.4d0,
_IFN1(cfu)     7          r2o3=0.666666666666667d0)
      data t0,t1,t2,t3,t4/28.0d0,32.0d0,35.0d0,38.0d0,40.0d0/
c
c     r.j.harrison 6/3/87.
c
c     computes f0,f1,f2,f3,f4 for all values of the argument
c     to a relative accuracy of better then 4.e-13 for all.
c     uses 4th order taylor expansion on grid out to t=40.0
c     asymptotic expansion accurate for t greater than
c     f0(28),f1(32),f2(35),f3(38),f4(40).
c
c     branch to desired max value of m
c
      if(nval.eq.0) return
      goto (100,200,300,400,500),mhigh+1
c
100   continue
c     all loops vectorise, cray-xmp/48 cft1.15
c     on alliant 8CE (not super CE) always better for these loops
c     in scalar concurrent (even asymptotically !!)
      do 110 i=1,nval
          if(t(i).lt.t0) goto 110
          f(i,1)=fac0/dsqrt(t(i))
110   continue
c
_IF1()      call whenflt(nval,t,1,t0,ind,nlt)
_IF1(x)      call dlstlt(nval,t,1,t0,nlt,ind)
_IFN1(x)      nlt=0
_IFN1(x)      do 115 i=1,nval
_IFN1(x)          if(t(i).ge.t0) goto 115
_IFN1(x)          nlt=nlt+1
_IFN1(x)          ind(nlt)=i
_IFN1(x)115   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 120 ii=1,nlt
          i=ind(ii)
          ts=t(i)
          n=idint((ts+delo2)*rdelta)
          x=delta*dfloat(n)-ts
          n=n+1
          f(i,1)=fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     &       rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
120   continue
      return
c
200   continue
_IF1()      call whenfge(nval,t,1,t1,ind,nge)
_IF1(x)      call dlstge(nval,t,1,t1,nge,ind)
_IFN1(x)      nge=0
_IFN1(x)      do 205 i=1,nval
_IFN1(x)          if(t(i).lt.t1) goto 205
_IFN1(x)          nge=nge+1
_IFN1(x)          ind(nge)=i
_IFN1(x)205   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 210 ii=1,nge
          i=ind(ii)
          ts=t(i)
          f1=fac1/(ts*dsqrt(ts))
          f(i,1)=ts*(f1+f1)
          f(i,2)=f1
210   continue
c
_IF1()      call whenflt(nval,t,1,t1,ind,nlt)
_IF1(x)      call dlstlt(nval,t,1,t1,nlt,ind)
_IFN1(x)      nlt=0
_IFN1(x)      do 215 i=1,nval
_IFN1(x)          if(t(i).ge.t1) goto 215
_IFN1(x)          nlt=nlt+1
_IFN1(x)          ind(nlt)=i
_IFN1(x)215   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 220 ii=1,nlt
          i=ind(ii)
          ts=t(i)
          n=idint((ts+delo2)*rdelta)
          x=delta*dfloat(n)-ts
          n=n+1
          f(i,1)=fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     &       rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
          f(i,2)=fm(n,2)+x*(fm(n,3)+rhalf*x*(fm(n,4)+
     &       rthird*x*(fm(n,5)+rquart*x*fm(n,6))))
220   continue
      return
c
300   continue
_IF1()      call whenfge(nval,t,1,t2,ind,nge)
_IF1(x)      call dlstge(nval,t,1,t2,nge,ind)
_IFN1(x)      nge=0
_IFN1(x)      do 305 i=1,nval
_IFN1(x)          if(t(i).lt.t2) goto 305
_IFN1(x)          nge=nge+1
_IFN1(x)          ind(nge)=i
_IFN1(x)305   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 310 ii=1,nge
          i=ind(ii)
          ts=t(i)
          f2=fac2/(ts*ts*dsqrt(ts))
          f1=r2o3*ts*f2
          f(i,1)=ts*(f1+f1)
          f(i,2)=f1
          f(i,3)=f2
310   continue
c
_IF1()      call whenflt(nval,t,1,t2,ind,nlt)
_IF1(x)      call dlstlt(nval,t,1,t2,nlt,ind)
_IFN1(x)      nlt=0
_IFN1(x)      do 315 i=1,nval
_IFN1(x)          if(t(i).ge.t2) goto 315
_IFN1(x)          nlt=nlt+1
_IFN1(x)          ind(nlt)=i
_IFN1(x)315   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 320 ii=1,nlt
          i=ind(ii)
          ts=t(i)
          n=idint((ts+delo2)*rdelta)
          x=delta*dfloat(n)-ts
          n=n+1
          f(i,1)=fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     &       rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
          f(i,2)=fm(n,2)+x*(fm(n,3)+rhalf*x*(fm(n,4)+
     &       rthird*x*(fm(n,5)+rquart*x*fm(n,6))))
          f(i,3)=fm(n,3)+x*(fm(n,4)+rhalf*x*(fm(n,5)+
     &       rthird*x*(fm(n,6)+rquart*x*fm(n,7))))
320   continue
      return
c
400   continue
_IF1()      call whenfge(nval,t,1,t3,ind,nge)
_IF1(x)      call dlstge(nval,t,1,t3,nge,ind)
_IFN1(x)      nge=0
_IFN1(x)      do 405 i=1,nval
_IFN1(x)          if(t(i).lt.t3) goto 405
_IFN1(x)          nge=nge+1
_IFN1(x)          ind(nge)=i
_IFN1(x)405   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 410 ii=1,nge
          i=ind(ii)
          ts=t(i)
          tt=ts*ts
          f3=fac3/(ts*tt*dsqrt(ts))
          f2=r2o5*ts*f3
          f1=r2o3*ts*f2
          f0=ts*(f1+f1)
          f(i,1)=f0
          f(i,2)=f1
          f(i,3)=f2
          f(i,4)=f3
410   continue
c
_IF1()      call whenflt(nval,t,1,t3,ind,nlt)
_IF1(x)      call dlstlt(nval,t,1,t3,nlt,ind)
_IFN1(x)      nlt=0
_IFN1(x)      do 415 i=1,nval
_IFN1(x)          if(t(i).ge.t3) goto 415
_IFN1(x)          nlt=nlt+1
_IFN1(x)          ind(nlt)=i
_IFN1(x)415   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 420 ii=1,nlt
          i=ind(ii)
          ts=t(i)
          n=idint((ts+delo2)*rdelta)
          x=delta*dfloat(n)-ts
          n=n+1
          f(i,1)=fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     &       rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
          f(i,2)=fm(n,2)+x*(fm(n,3)+rhalf*x*(fm(n,4)+
     &       rthird*x*(fm(n,5)+rquart*x*fm(n,6))))
          f(i,3)=fm(n,3)+x*(fm(n,4)+rhalf*x*(fm(n,5)+
     &       rthird*x*(fm(n,6)+rquart*x*fm(n,7))))
          f(i,4)=fm(n,4)+x*(fm(n,5)+rhalf*x*(fm(n,6)+
     &       rthird*x*(fm(n,7)+rquart*x*fm(n,8))))
420   continue
      return
c
500   continue
_IF1()      call whenfge(nval,t,1,t4,ind,nge)
_IF1(x)      call dlstge(nval,t,1,t4,nge,ind)
_IFN1(x)      nge=0
_IFN1(x)      do 505 i=1,nval
_IFN1(x)          if(t(i).lt.t4) goto 505
_IFN1(x)          nge=nge+1
_IFN1(x)          ind(nge)=i
_IFN1(x)505   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 510 ii=1,nge
          i=ind(ii)
          ts=t(i)
          tt=ts*ts
          f4=fac4/(tt*tt*dsqrt(ts))
          f3=r2o7*ts*f4
          f2=r2o5*ts*f3
          f1=r2o3*ts*f2
          f0=ts*(f1+f1)
          f(i,1)=f0
          f(i,2)=f1
          f(i,3)=f2
          f(i,4)=f3
          f(i,5)=f4
510   continue
c
_IF1()      call whenflt(nval,t,1,t4,ind,nlt)
_IF1(x)      call dlstlt(nval,t,1,t4,nlt,ind)
_IFN1(x)      nlt=0
_IFN1(x)      do 515 i=1,nval
_IFN1(x)          if(t(i).ge.t4) goto 515
_IFN1(x)          nlt=nlt+1
_IFN1(x)          ind(nlt)=i
_IFN1(x)515   continue
_IF1(x)c$dir no_recurrence
_IF1(i)c*vdir: ignore recrdeps
_IF1(a)cvd$  nodepchk
      do 520 ii=1,nlt
          i=ind(ii)
          ts=t(i)
          n=idint((ts+delo2)*rdelta)
          x=delta*dfloat(n)-ts
          n=n+1
          f(i,1)=fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     &       rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
          f(i,2)=fm(n,2)+x*(fm(n,3)+rhalf*x*(fm(n,4)+
     &       rthird*x*(fm(n,5)+rquart*x*fm(n,6))))
          f(i,3)=fm(n,3)+x*(fm(n,4)+rhalf*x*(fm(n,5)+
     &       rthird*x*(fm(n,6)+rquart*x*fm(n,7))))
          f(i,4)=fm(n,4)+x*(fm(n,5)+rhalf*x*(fm(n,6)+
     &       rthird*x*(fm(n,7)+rquart*x*fm(n,8))))
          f(i,5)=fm(n,5)+x*(fm(n,6)+rhalf*x*(fm(n,7)+
     &       rthird*x*(fm(n,8)+rquart*x*fm(n,9))))
520   continue
      return
      end
      subroutine qout70(gout)
_IF1(a)cvd$r novector
c
c ... qout70 now handles conventional integrals only
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
      dimension gout(*)
      dimension ib(4,4)
      common/blkin/goutx(510),nword
_IFN1(iv)      common/craypk/integ(680)
_IF1(iv)      common/craypk/integ(340),intkl(340)
INCLUDE(common/cslosc)
      common/flip70 /ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -            lab1q(256),lab2q(256),iptge(256)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
      ifac=256
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 210 i = 1,maxi
            if (oianj) jmax = i
            i1 = loci + i
            ipack = i4096(i1)
            do 200 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2 = locj + j
               if(i1-i2)120, 130, 130
120            lab1=i4096(i2) + i1
               go to 140
130            lab1=ipack + i2
140            kln = 0
               do 190 k = 1,maxk
                  if (okanl) lmax = k
                  i3 = lock+k
                  kpack = i4096(i3)
                  do 180 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 200
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = gout(nn)
                     if ( dabs(val) .lt. cutoff) go to 180
                     i4 = locl+l
                     if(i3-i4) 150, 160, 160
150                  lab2=i4096(i4) + i3
                     go to 170
160                  lab2=kpack + i4
170                  goutx(icount) = val
_IFN1(iv)                     integ(ic4 )= max(lab1,lab2)
_IFN1(iv)                     integ(ic4+1)= min(lab1,lab2)
_IFN1(iv)                     ic4 = ic4 + 2
_IF1(iv)                integ(icount)= max(lab1,lab2)
_IF1(iv)                intkl(icount)= min(lab1,lab2)
                     icount = icount+1
                     if (icount .le. nintmx) go to 180
                     call blocki
                     if(omaxb)go to 220
180               continue
190            continue
200         continue
210      continue
220   return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in qout70.')
      endif
c mainfile output
         if(.not.oianj) then
            do 440 iq=1,ngout
               i=indsp(iq,ib1)+loci
               j=indsp(iq,jb1)+locj
440         lab1q(iq)=i*ifac+j
         else
            do 450 iq=1,ngout
               i=indsp(iq,ib1)+loci
               j=indsp(iq,jb1)+locj
_IF1(c)450         lab1q(iq)=cvmgt(i*ifac+j,0,i.ge.j)
_IFN1(c)         if(i.ge.j) then
_IFN1(c)          lab1q(iq)=i*ifac+j
_IFN1(c)          else
_IFN1(c)          lab1q(iq)=0
_IFN1(c)          endif
_IFN1(c) 450      continue
         endif
         if(.not.okanl) then
            do 460 iq=1,ngout
               k=indsp(iq,kb1)+lock
               l=indsp(iq,lb1)+locl
460         lab2q(iq)=k*ifac+l
         else
            do 470 iq=1,ngout
               k=indsp(iq,kb1)+lock
               l=indsp(iq,lb1)+locl
_IF1(c)470         lab2q(iq)=cvmgt(k*ifac+l,0,k.ge.l)
_IFN1(c)         if(k.ge.l) then
_IFN1(c)         lab2q(iq)=k*ifac+l
_IFN1(c)         else
_IFN1(c)         lab2q(iq)=0
_IFN1(c)         endif
_IFN1(c) 470     continue
         endif
         if(oident) then
            do 480 iq=1,ngout
_IF1(c)480         lab1q(iq)=cvmgt(lab1q(iq),0,lab1q(iq).ge.lab2q(iq))
_IFN1(c)         if(lab1q(iq).lt.lab2q(iq)) then
_IFN1(c)         lab1q(iq) = 0
_IFN1(c)         endif
_IFN1(c) 480     continue
         elseif(ishell.eq.kshell) then
               do 490 iq=1,ngout
                  lab1=lab1q(iq)
                  lab2=lab2q(iq)
                  lab1q(iq)=max(lab1,lab2)
490            lab2q(iq)=min(lab1,lab2)
         endif
_IF1()c     do 500 iq=1,ngout
_IF1()c500  gout(iq)=abs(g(iq))
_IF1()c     call whenfge(ngout,gout,1,cutoff,iptge,nval)
_IF1()c do 501 vectorises on the xmp
         nval=0
         do 501 iq=1,ngout
            if(dabs(gout(iq)).lt.cutoff) goto 501
            nval=nval+1
            iptge(nval)=iq
501      continue
         iptq=1
510      nleft=min(nval,nintmx-icount+1)
         if(.not.(oident.or.oianj.or.okanl)) then
            nok=nok+1
            do 520 iiq=iptq,iptq+nleft-1
               iq=iptge(iiq)
_IFN1(iv)               integ(ic4 )=lab1q(iq)
_IFN1(iv)               integ(ic4+1)=lab2q(iq)
_IFN1(iv)               ic4=ic4+2
_IF1(iv)               integ(icount )=lab1q(iq)
_IF1(iv)               intkl(icount)=lab2q(iq)
               goutx(icount)=gout(iq)
520         icount=icount+1
         else
            nnok=nnok+1
            do 530 iiq=iptq,iptq+nleft-1
               iq=iptge(iiq)
               if(lab1q(iq)*lab2q(iq).eq.0) goto 530
_IFN1(iv)               integ(ic4 )=lab1q(iq)
_IFN1(iv)               integ(ic4+1)=lab2q(iq)
_IFN1(iv)               ic4=ic4+2
_IF1(iv)               integ(icount )=lab1q(iq)
_IF1(iv)               intkl(icount)=lab2q(iq)
               goutx(icount)=gout(iq)
               icount=icount+1
530         continue
         endif
         if(icount.gt.nintmx) then
               call blocki
            if(omaxb) goto 220
         endif
         nval=nval-nleft
         if(nval.le.0) goto 220
         iptq=iptq+nleft
         goto 510
c
      end
_IF(unicos)
_IF(ccpdft)
      subroutine qoutd70(fock,dmat,g,
     +    fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine qoutd70(fock,dmat,g)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c ***
c ***
      dimension fock(*),dmat(*),g(*)
      dimension ib(4,4)
      common/blkin/goutx(510),nword
      common/craypk/integ(680)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -                 lab1q(256),lab2q(256),iptge(256)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
_IF(ccpdft)
      if (.not. odft) then
      ifac = 4096
_ENDIF
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 110 i = 1,maxi
            i1=loci+i
            ipack=ifac*i1
            if (oianj) jmax = i
            do 100 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2=locj+j
               if(i1-i2) 20,30,30
20             lab1=ifac*i2+i1
               go to 40
30             lab1=ipack+i2
40             continue
               kln=0
               do 90 k = 1,maxk
                  if (okanl) lmax = k
                  i3=lock+k
                  kpack=ifac*i3
                  do 80 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 100
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = g(nn)
                     if ( dabs(val) .lt. cutoff) go to 80
                     goutx(icount) = val
                     i4=locl+l
                     if(i3-i4)50,60,60
50                   lab2=ifac*i4+i3
                     go to 70
60                   lab2=kpack+i4
70                   integ(ic4 )=max(lab1,lab2)
                     integ(ic4+1)=min(lab1,lab2)
                     ic4=ic4+2
                     icount=icount+1
                     if(icount.le.nintmx) go to 80
                     ochek=.false.
                     nrec=nrec+1
                     call dbuild(fock,dmat)
80                continue
90             continue
100         continue
110      continue
220   return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in qout70.')
      endif
         if(.not.oianj) then
            do 240 iq=1,ngout
               i=indsp(iq,ib1)+loci
               j=indsp(iq,jb1)+locj
240         lab1q(iq)=i*ifac+j
         else
            do 250 iq=1,ngout
               i=indsp(iq,ib1)+loci
               j=indsp(iq,jb1)+locj
250         lab1q(iq)=cvmgt(i*ifac+j,0,i.ge.j)
         endif
         if(.not.okanl) then
            do 260 iq=1,ngout
               k=indsp(iq,kb1)+lock
               l=indsp(iq,lb1)+locl
260         lab2q(iq)=k*ifac+l
         else
            do 270 iq=1,ngout
               k=indsp(iq,kb1)+lock
               l=indsp(iq,lb1)+locl
270         lab2q(iq)=cvmgt(k*ifac+l,0,k.ge.l)
            endif
         if(oident) then
            do 280 iq=1,ngout
280         lab1q(iq)=cvmgt(lab1q(iq),0,lab1q(iq).ge.lab2q(iq))
         elseif(ishell.eq.kshell) then
               do 290 iq=1,ngout
                  lab1=lab1q(iq)
                  lab2=lab2q(iq)
                  lab1q(iq)=max(lab1,lab2)
290            lab2q(iq)=min(lab1,lab2)
         endif
            nval=0
         do 301 iq=1,ngout
            if(dabs(g(iq)).lt.cutoff) goto 301
            nval=nval+1
            iptge(nval)=iq
301      continue
         iptq=1
310      nleft=min(nval,nintmx-icount+1)
         if(.not.(oident.or.oianj.or.okanl)) then
            nok=nok+1
            do 320 iiq=iptq,iptq+nleft-1
               iq=iptge(iiq)
               integ(ic4 )=lab1q(iq)
               integ(ic4+1)=lab2q(iq)
               ic4=ic4+2
               goutx(icount)=g(iq)
320         icount=icount+1
         else
            nnok=nnok+1
            do 330 iiq=iptq,iptq+nleft-1
               iq=iptge(iiq)
               if(lab1q(iq)*lab2q(iq).eq.0) goto 330
               integ(ic4 )=lab1q(iq)
               integ(ic4+1)=lab2q(iq)
               ic4=ic4+2
               goutx(icount)=g(iq)
               icount=icount+1
330         continue
         endif
         if(icount.gt.nintmx) then
               ochek=.false.
               nrec=nrec+1
               call dbuild(fock,dmat)
         endif
         nval=nval-nleft
         if(nval.le.0) goto 220
         iptq=iptq+nleft
         goto 310
c
_IF(ccpdft)
      else
c
c     non conventional fock build involved .. cannot use CAL
c
      if(jtype.ge.5 .or. jtype.eq.3) goto 2300
c ***
      ijn = 0
      jmax = maxj
      do 460 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 440 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 420 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 400 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 440
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 400
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
c
c Coulomb terms
c
      if(ocoul)then
         itr12=iky(i1)+i2
         itr34=iky(i3)+i4
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c   
c Full exchange term for HF or weighted exchange term for b3lyp etc
c   
         itr13=iky(i1)+i3
         itr14=iky(i1)+i4
         itr23=iky(max(i2,i3))+min(i2,i3)
         itr24=iky(max(i2,i4))+min(i2,i4)
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
  400 continue
  420 continue
  440 continue
  460 continue
c
      return
2300  if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
         ngout=64
      elseif(jtype.eq.3) then
         ngout=16
      else
         call caserr('confusion in qoutd70')
      endif
          if(.not.(oident.or.oianj.or.okanl)) then
            do 2400 iq=1,ngout
            val = g(iq)
            if( dabs(val).lt.cutoff) goto 2400
               i1 = indsp(iq,ib1)+loci
               i2 = indsp(iq,jb1)+locj
               i3 = indsp(iq,kb1)+lock
               i4 = indsp(iq,lb1)+locl
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               ii3 = max(i3,i4)
               ii4 = min(i3,i4)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
                  facij=fac1
                  fackl=fac2
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
                  facij=fac2
                  fackl=fac1
               endif
c
c Coulomb terms
c
               if(ocoul)then
                  itr12=iky(i1)+i2
                  itr34=iky(i3)+i4
                  val2=val+val
                  val4=val2+val2
                  f12 = facij*val4*dmat(itr34) + fock(itr12)
                  fock(itr12) = f12
                  if(itr12 .ne. itr34)then
                     fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
                  endif
               endif
c
               if(oexch)then
c Full exchange term for HF or weighted exchange term for b3lyp etc
c
                  itr13=iky(i1)+i3
                  itr14=iky(i1)+i4
                  itr23=iky(max(i2,i3))+min(i2,i3)
                  itr24=iky(max(i2,i4))+min(i2,i4)
                  val=val*facex
                  val2=val+val
                  val13=val
                  val14=val
                  if(i1.eq.i3 .or. i2.eq.i4) val13=val2
                  if(i2.eq.i3) val14=val2
                  f23 = fock(itr23) - val14*dmat(itr14)
                  f14 = fock(itr14) - val14*dmat(itr23)
                  f13 = fock(itr13) - val13*dmat(itr24)
                  fock(itr24) = fock(itr24) - val13*dmat(itr13)
                  fock(itr23) = f23
                  fock(itr14) = f14
                  fock(itr13) = f13
               endif
2400        continue
           else
            do 2500 iq=1,ngout
            val =g(iq)
            if( dabs(val).lt.cutoff) goto 2500
               i = indsp(iq,ib1)+loci
               j = indsp(iq,jb1)+locj
                if(oianj.and.i.lt.j)  go to 2500
               lab1q(iq) = iky(i) + j
               k = indsp(iq,kb1)+lock
               l = indsp(iq,lb1)+locl
                if(okanl.and.k.lt.l)  go to 2500
               lab2q(iq) = iky(k) + l
                if(oident.and.lab1q(iq).lt.lab2q(iq)) go to 2500
               ii1 = max(i,j)
               ii2 = min(i,j)
               ii3 = max(k,l)
               ii4 = min(k,l)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
                  facij=fac1
                  fackl=fac2
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
                  facij=fac2
                  fackl=fac1
               endif
c
c Coulomb terms
c
               if(ocoul)then
                 itr12=iky(i1)+i2
                 itr34=iky(i3)+i4
                 val2=val+val
                 val4=val2+val2
                 f12 = facij*val4*dmat(itr34) + fock(itr12)
                 fock(itr12) = f12
                 if(itr12 .ne. itr34)then
                    fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
                 endif
               endif
c
               if (oexch)then
c    
c Full exchange term for HF or weighted exchange term for b3lyp etc
c    
                 itr13=iky(i1)+i3
                 itr14=iky(i1)+i4
                 itr23=iky(max(i2,i3))+min(i2,i3)
                 itr24=iky(max(i2,i4))+min(i2,i4)
                 val=val*facex
                 val2=val+val
                 val13=val
                 val14=val
                 if(i1.eq.i3 .or. i2.eq.i4) val13=val2
                 if(i2.eq.i3) val14=val2
                 f23 = fock(itr23) - val14*dmat(itr14)
                 f14 = fock(itr14) - val14*dmat(itr23)
                 f13 = fock(itr13) - val13*dmat(itr24)
                 fock(itr24) = fock(itr24) - val13*dmat(itr13)
                 fock(itr23) = f23
                 fock(itr14) = f14
                 fock(itr13) = f13
               endif
2500        continue
           endif
c
      endif
_ENDIF
      return
      end
_ELSE
_IF(ccpdft)
      subroutine dbuild70(fock,dmat,gout,
     +     fac1,fac2, facex, ocoul, oexch)
_ELSE
      subroutine dbuild70(fock,dmat,gout)
_ENDIF
_IF1(a)cvd$r novector
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
c ***
c ***
      dimension fock(*),dmat(*),gout(*)
      dimension ib(4,4)
      common/blkin/goutx(510),nword
_IFN1(iv)      common/craypk/integ(680)
_IF1(iv)      common/craypk/integ(340),intkl(340)
INCLUDE(common/cslosc)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -            lab1q(256),lab2q(256),iptge(256)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
      ifac=256
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 110 i = 1,maxi
            i1=loci+i
            if (oianj) jmax = i
            do 100 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2=locj+j
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               kln=0
               do 90 k = 1,maxk
                  if (okanl) lmax = k
                  i3=lock+k
                  do 80 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 100
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = gout(nn)
                     if ( dabs(val) .lt. cutoff) go to 80
                     i4=locl+l
                     ii3 = max(i3,i4)
                     ii4 = min(i3,i4)
                     if(ii1.ge.ii3) then
                        i1 = ii1
                        i2 = ii2
                        i3 = ii3
                        i4 = ii4
_IF(ccpdft)
                        facij=fac1
                        fackl=fac2
_ENDIF
                     else
                        i1 = ii3
                        i2 = ii4
                        i3 = ii1
                        i4 = ii2
_IF(ccpdft)
                        facij=fac2
                        fackl=fac1
_ENDIF
                     endif
                     itr12=iky(i1)+i2
                     itr13=iky(i1)+i3
                     itr14=iky(i1)+i4
                     itr34=iky(i3)+i4
                     itr23=iky(max(i2,i3))+min(i2,i3)
                     itr24=iky(max(i2,i4))+min(i2,i4)
_IF(ccpdft)
c
c Coulomb terms
c
      if(ocoul)then
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c     
c Full exchange term for HF or weighted exchange term for b3lyp etc
c     
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
_ELSE
                     val2=val+val
                     val4=val2+val2
                     val13=val
                     val14=val
                     if(i1.eq.i3 .or. i2.eq.i4) val13=val2
                     if(i2.eq.i3) val14=val2
                     f12 = val4*dmat(itr34) + fock(itr12)
                     fock(itr34) = val4*dmat(itr12) + fock(itr34)
                     fock(itr12) = f12
                     f23 = fock(itr23) - val14*dmat(itr14)
                     f14 = fock(itr14) - val14*dmat(itr23)
                     f13 = fock(itr13) - val13*dmat(itr24)
                     fock(itr24) = fock(itr24) - val13*dmat(itr13)
                     fock(itr23) = f23
                     fock(itr14) = f14
                     fock(itr13) = f13
_ENDIF
 80               continue
 90            continue
100         continue
110      continue
      return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in qout70.')
      endif
          if(.not.(oident.or.oianj.or.okanl)) then
            do 240 iq=1,ngout
            val = gout(iq)
            if( dabs(val).lt.cutoff) goto 240
               i1 = indsp(iq,ib1)+loci
               i2 = indsp(iq,jb1)+locj
               i3 = indsp(iq,kb1)+lock
               i4 = indsp(iq,lb1)+locl
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               ii3 = max(i3,i4)
               ii4 = min(i3,i4)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
_IF(ccpdft)
                  facij=fac1
                  fackl=fac2
_ENDIF
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
_IF(ccpdft)
                  facij=fac2
                  fackl=fac1
_ENDIF
               endif
               itr12=iky(i1)+i2
               itr13=iky(i1)+i3
               itr14=iky(i1)+i4
               itr34=iky(i3)+i4
               itr23=iky(max(i2,i3))+min(i2,i3)
               itr24=iky(max(i2,i4))+min(i2,i4)
_IF(ccpdft)
c
c Coulomb terms
c
      if(ocoul)then
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c     
c Full exchange term for HF or weighted exchange term for b3lyp etc
c     
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
_ELSE
               val2=val+val
               val4=val2+val2
               val13=val
               val14=val
               if(i1.eq.i3 .or. i2.eq.i4) val13=val2
               if(i2.eq.i3) val14=val2
               f12 = val4*dmat(itr34) + fock(itr12)
               fock(itr34) = val4*dmat(itr12) + fock(itr34)
               fock(itr12) = f12
               f23 = fock(itr23) - val14*dmat(itr14)
               f14 = fock(itr14) - val14*dmat(itr23)
               f13 = fock(itr13) - val13*dmat(itr24)
               fock(itr24) = fock(itr24) - val13*dmat(itr13)
               fock(itr23) = f23
               fock(itr14) = f14
               fock(itr13) = f13
_ENDIF
240         continue
           else
            do 250 iq=1,ngout
            val =gout(iq)
            if( dabs(val).lt.cutoff) goto 250
               i = indsp(iq,ib1)+loci
               j = indsp(iq,jb1)+locj
                if(oianj.and.i.lt.j)  go to 250
               lab1q(iq) = iky(i) + j
               k = indsp(iq,kb1)+lock
               l = indsp(iq,lb1)+locl
                if(okanl.and.k.lt.l)  go to 250
               lab2q(iq) = iky(k) + l
                if(oident.and.lab1q(iq).lt.lab2q(iq)) go to 250
               ii1 = max(i,j)
               ii2 = min(i,j)
               ii3 = max(k,l)
               ii4 = min(k,l)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
_IF(ccpdft)
                  facij=fac1
                  fackl=fac2
_ENDIF
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
_IF(ccpdft)
                  facij=fac2
                  fackl=fac1
_ENDIF
               endif
               itr12=iky(i1)+i2
               itr13=iky(i1)+i3
               itr14=iky(i1)+i4
               itr34=iky(i3)+i4
               itr23=iky(max(i2,i3))+min(i2,i3)
               itr24=iky(max(i2,i4))+min(i2,i4)
_IF(ccpdft)
c
c Coulomb terms
c
      if(ocoul)then
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c     
c Full exchange term for HF or weighted exchange term for b3lyp etc
c     
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
_ELSE
               val2=val+val
               val4=val2+val2
               val13=val
               val14=val
               if(i1.eq.i3 .or. i2.eq.i4) val13=val2
               if(i2.eq.i3) val14=val2
               f12 = val4*dmat(itr34) + fock(itr12)
               fock(itr34) = val4*dmat(itr12) + fock(itr34)
               fock(itr12) = f12
               f23 = fock(itr23) - val14*dmat(itr14)
               f14 = fock(itr14) - val14*dmat(itr23)
               f13 = fock(itr13) - val13*dmat(itr24)
               fock(itr24) = fock(itr24) - val13*dmat(itr13)
               fock(itr23) = f23
               fock(itr14) = f14
               fock(itr13) = f13
_ENDIF
250         continue
           endif
      end
_ENDIF
_IF(ccpdft)
      subroutine dir_build_uhf70(fock,ak,p,q,gout,
     + fac1, fac2, facex, ocoul, oexch)
_ELSE
      subroutine dir_build_uhf70(fock,ak,p,q,gout)
_ENDIF
_IF1(a)cvd$r novector
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
      dimension gout(*),fock(*),ak(*),p(*),q(*)
      dimension ib(4,4)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -            lab1q(256),lab2q(256),iptge(256)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
c
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell UHF Fock builder 
c
_IF1(f)      magint(v)=max(extract(v,0,11)-990,1)
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
      ifac=256
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 110 i = 1,maxi
            i1=loci+i
            if (oianj) jmax = i
            do 100 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2=locj+j
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               kln=0
               do 90 k = 1,maxk
                  if (okanl) lmax = k
                  i3=lock+k
                  do 80 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 100
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = gout(nn)
                     if ( dabs(val) .lt. cutoff) go to 80
                     i4=locl+l
                     ii3 = max(i3,i4)
                     ii4 = min(i3,i4)
                     if(ii1.ge.ii3) then
                        i1 = ii1
                        i2 = ii2
                        i3 = ii3
                        i4 = ii4
_IF(ccpdft)
                        facij=fac1
                        fackl=fac2
_ENDIF
                     else
                        i1 = ii3
                        i2 = ii4
                        i3 = ii1
                        i4 = ii2
_IF(ccpdft)
                        facij=fac2
                        fackl=fac1
_ENDIF
                     endif
                    gik=val
                    val2=gik+gik
                    val4=val2+val2
                    ikyi=iky(i1)
                    ikyj=iky(i2)
                    ikyk=iky(i3)
                    itr13=ikyi+i3
                    itr14=ikyi+i4
                    itr12=ikyi+i2
                    itr23=ikyj+i3
                    itr24=ikyj+i4
                    itr34=ikyk+i4
_IF(ccpdft)
c
c Coulomb term
c
                    if(ocoul)then
                      gik=val
                      val2=gik+gik
                      val4=val2+val2
                      fock(itr12) = facij*val4*p(itr34) + fock(itr12)
                      if(itr12 .ne. itr34)then
                         fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
                      endif
                   endif
c
c exchange
c
                   if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
                      gik=val*facex
                      val2=gik+gik
                      val4=val2+val2
                      gil=gik
                      if(i1.eq.i3.or.i2.eq.i4)gik=val2
                      if(i2.eq.i3)gil=val2
                      if(i2.ge.i3)goto 1
                      itr23=ikyk+i2
                      if(i2.ge.i4)goto 1
                      itr24=iky(i4)+i2
1                     ajk=fock(itr23)-gil*p(itr14)
                      bjk=ak(itr23)+gil*q(itr14)
                      ail=fock(itr14)-gil*p(itr23)
                      bil=ak(itr14)+gil*q(itr23)
                      aik=fock(itr13)-gik*p(itr24)
                      bik=ak(itr13)+gik*q(itr24)
                      fock(itr24)=fock(itr24)-gik*p(itr13)
                      ak(itr24)=ak(itr24)+gik*q(itr13)
                      fock(itr23)=ajk
                      ak(itr23)=bjk
                      fock(itr14)=ail
                      ak(itr14)=bil
                      fock(itr13)=aik
                      ak(itr13)=bik
                   endif
_ELSE
c
c  non-dft UHF version
c
                    aij=val4*p(itr34)+fock(itr12)
                    fock(itr34)=val4*p(itr12)+fock(itr34)
                    fock(itr12)=aij
c... exchange
                    gil=gik
                    if(i1.eq.i3.or.i2.eq.i4)gik=val2
                    if(i2.eq.i3)gil=val2
                    if(i2.ge.i3)goto 10
                    itr23=ikyk+i2
                    if(i2.ge.i4)goto 10
                    itr24=iky(i4)+i2
10                  ajk=fock(itr23)-gil*p(itr14)
                    bjk=ak(itr23)+gil*q(itr14)
                    ail=fock(itr14)-gil*p(itr23)
                    bil=ak(itr14)+gil*q(itr23)
                    aik=fock(itr13)-gik*p(itr24)
                    bik=ak(itr13)+gik*q(itr24)
                    fock(itr24)=fock(itr24)-gik*p(itr13)
                    ak(itr24)=ak(itr24)+gik*q(itr13)
                    fock(itr23)=ajk
                    ak(itr23)=bjk
                    fock(itr14)=ail
                    ak(itr14)=bil
                    fock(itr13)=aik
                    ak(itr13)=bik
_ENDIF
 80               continue
 90            continue
100         continue
110      continue
      return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in qout70.')
      endif
          if(.not.(oident.or.oianj.or.okanl)) then
            do 240 iq=1,ngout
            val = gout(iq)
            if( dabs(val).lt.cutoff) goto 240
               i1 = indsp(iq,ib1)+loci
               i2 = indsp(iq,jb1)+locj
               i3 = indsp(iq,kb1)+lock
               i4 = indsp(iq,lb1)+locl
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               ii3 = max(i3,i4)
               ii4 = min(i3,i4)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
_IF(ccpdft)
                        facij=fac1
                        fackl=fac2
_ENDIF
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
_IF(ccpdft)
                        facij=fac2
                        fackl=fac1
_ENDIF
               endif
               gik=val
               val2=gik+gik
               val4=val2+val2
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
_IF(ccpdft)
c
c Coulomb term
c
                if(ocoul)then
                  gik=val
                  val2=gik+gik
                  val4=val2+val2
                  fock(itr12) = facij*val4*p(itr34) + fock(itr12)
                  if(itr12 .ne. itr34)then
                     fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
                  endif
               endif
c
c exchange
c
               if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
                  gik=val*facex
                  val2=gik+gik
                  val4=val2+val2
                  gil=gik
                  if(i1.eq.i3.or.i2.eq.i4)gik=val2
                  if(i2.eq.i3)gil=val2
                  if(i2.ge.i3)goto 20
                  itr23=ikyk+i2
                  if(i2.ge.i4)goto 20
                  itr24=iky(i4)+i2
20                ajk=fock(itr23)-gil*p(itr14)
                  bjk=ak(itr23)+gil*q(itr14)
                  ail=fock(itr14)-gil*p(itr23)
                  bil=ak(itr14)+gil*q(itr23)
                  aik=fock(itr13)-gik*p(itr24)
                  bik=ak(itr13)+gik*q(itr24)
                  fock(itr24)=fock(itr24)-gik*p(itr13)
                  ak(itr24)=ak(itr24)+gik*q(itr13)
                  fock(itr23)=ajk
                  ak(itr23)=bjk
                  fock(itr14)=ail
                  ak(itr14)=bil
                  fock(itr13)=aik
                  ak(itr13)=bik
               endif
_ELSE
c
c  non-dft UHF version
c
               aij=val4*p(itr34)+fock(itr12)
               fock(itr34)=val4*p(itr12)+fock(itr34)
               fock(itr12)=aij
c... exchange
               gil=gik
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto 20
               itr23=ikyk+i2
               if(i2.ge.i4)goto 20
               itr24=iky(i4)+i2
20             ajk=fock(itr23)-gil*p(itr14)
               bjk=ak(itr23)+gil*q(itr14)
               ail=fock(itr14)-gil*p(itr23)
               bil=ak(itr14)+gil*q(itr23)
               aik=fock(itr13)-gik*p(itr24)
               bik=ak(itr13)+gik*q(itr24)
               fock(itr24)=fock(itr24)-gik*p(itr13)
               ak(itr24)=ak(itr24)+gik*q(itr13)
               fock(itr23)=ajk
               ak(itr23)=bjk
               fock(itr14)=ail
               ak(itr14)=bil
               fock(itr13)=aik
               ak(itr13)=bik
_ENDIF
240         continue
           else
            do 250 iq=1,ngout
            val =gout(iq)
            if( dabs(val).lt.cutoff) goto 250
               i = indsp(iq,ib1)+loci
               j = indsp(iq,jb1)+locj
                if(oianj.and.i.lt.j)  go to 250
               lab1q(iq) = iky(i) + j
               k = indsp(iq,kb1)+lock
               l = indsp(iq,lb1)+locl
                if(okanl.and.k.lt.l)  go to 250
               lab2q(iq) = iky(k) + l
                if(oident.and.lab1q(iq).lt.lab2q(iq)) go to 250
               ii1 = max(i,j)
               ii2 = min(i,j)
               ii3 = max(k,l)
               ii4 = min(k,l)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
_IF(ccpdft)
                  facij=fac1
                  fackl=fac2
_ENDIF
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
_IF(ccpdft)
                  facij=fac2
                  fackl=fac1
_ENDIF
               endif
               gik=val
               val2=gik+gik
               val4=val2+val2
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
_IF(ccpdft)
c
c Coulomb term
c
               if(ocoul)then
                 gik=val
                 val2=gik+gik
                 val4=val2+val2
                 fock(itr12) = facij*val4*p(itr34) + fock(itr12)
                 if(itr12 .ne. itr34)then
                    fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
                 endif
               endif
c
c exchange
c
               if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
                 gik=val*facex
                 val2=gik+gik
                 val4=val2+val2
                 gil=gik
                 if(i1.eq.i3.or.i2.eq.i4)gik=val2
                 if(i2.eq.i3)gil=val2
                 if(i2.ge.i3)goto 30
                 itr23=ikyk+i2
                 if(i2.ge.i4)goto 30
                 itr24=iky(i4)+i2
30               ajk=fock(itr23)-gil*p(itr14)
                 bjk=ak(itr23)+gil*q(itr14)
                 ail=fock(itr14)-gil*p(itr23)
                 bil=ak(itr14)+gil*q(itr23)
                 aik=fock(itr13)-gik*p(itr24)
                 bik=ak(itr13)+gik*q(itr24)
                 fock(itr24)=fock(itr24)-gik*p(itr13)
                 ak(itr24)=ak(itr24)+gik*q(itr13)
                 fock(itr23)=ajk
                 ak(itr23)=bjk
                 fock(itr14)=ail
                 ak(itr14)=bil
                 fock(itr13)=aik
                 ak(itr13)=bik
               endif
_ELSE
c
c  non-dft UHF version
c
               aij=val4*p(itr34)+fock(itr12)
               fock(itr34)=val4*p(itr12)+fock(itr34)
               fock(itr12)=aij
c... exchange
               gil=gik
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto 30
               itr23=ikyk+i2
               if(i2.ge.i4)goto 30
               itr24=iky(i4)+i2
30             ajk=fock(itr23)-gil*p(itr14)
               bjk=ak(itr23)+gil*q(itr14)
               ail=fock(itr14)-gil*p(itr23)
               bil=ak(itr14)+gil*q(itr23)
               aik=fock(itr13)-gik*p(itr24)
               bik=ak(itr13)+gik*q(itr24)
               fock(itr24)=fock(itr24)-gik*p(itr13)
               ak(itr24)=ak(itr24)+gik*q(itr13)
               fock(itr23)=ajk
               ak(itr23)=bjk
               fock(itr14)=ail
               ak(itr14)=bil
               fock(itr13)=aik
               ak(itr13)=bik
_ENDIF
250         continue
           endif
      end
      subroutine dir_build_open_70(coul,exch,dens,gout)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -            lab1q(256),lab2q(256),iptge(256)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
c
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell SCF J & K builder (nshell=1)
c
c
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
cdir$ list
cdir$ novector
      ifac=256
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 110 i = 1,maxi
            i1=loci+i
            if (oianj) jmax = i
            do 100 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2=locj+j
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               kln=0
               do 90 k = 1,maxk
                  if (okanl) lmax = k
                  i3=lock+k
                  do 80 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 100
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = gout(nn)
                     if ( dabs(val) .lt. cutoff) go to 80
                     i4=locl+l
                     ii3 = max(i3,i4)
                     ii4 = min(i3,i4)
                     if(ii1.ge.ii3) then
                        i1 = ii1
                        i2 = ii2
                        i3 = ii3
                        i4 = ii4
                     else
                        i1 = ii3
                        i2 = ii4
                        i3 = ii1
                        i4 = ii2
                     endif
                     gik=val
                     val2=val+val
                     ikyi=iky(i1)
                     ikyj=iky(i2)
                     ikyk=iky(i3)
                     itr13=ikyi+i3
                     itr14=ikyi+i4
                     itr12=ikyi+i2
                     itr23=ikyj+i3
                     itr24=ikyj+i4
                     itr34=ikyk+i4
                     gil=val
                     if(i1.eq.i3.or.i2.eq.i4)gik=val2
                     if(i2.eq.i3)gil=val2
                     if(i2.ge.i3)goto 280
                     itr23=ikyk+i2
                     if(i2.ge.i4)goto 280
                     itr24=iky(i4)+i2
  280                bij=val2*dens(itr34)+coul(itr12)
                     coul(itr34)=val2*dens(itr12)+coul(itr34)
                     coul(itr12)=bij
                     bjk=exch(itr23)+gil*dens(itr14)
                     bil=exch(itr14)+gil*dens(itr23)
                     bik=exch(itr13)+gik*dens(itr24)
                     exch(itr24)=exch(itr24)+gik*dens(itr13)
                     exch(itr23)=bjk
                     exch(itr14)=bil
                     exch(itr13)=bik
 80               continue
 90            continue
100         continue
110      continue
      return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in qout70.')
      endif
          if(.not.(oident.or.oianj.or.okanl)) then
            do 240 iq=1,ngout
            val = gout(iq)
            if( dabs(val).lt.cutoff) goto 240
               i1 = indsp(iq,ib1)+loci
               i2 = indsp(iq,jb1)+locj
               i3 = indsp(iq,kb1)+lock
               i4 = indsp(iq,lb1)+locl
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               ii3 = max(i3,i4)
               ii4 = min(i3,i4)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
               endif
               gik=val
               val2=val+val
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
               gil=val
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto  20
               itr23=ikyk+i2
               if(i2.ge.i4)goto  20
               itr24=iky(i4)+i2
   20          bij=val2*dens(itr34)+coul(itr12)
               coul(itr34)=val2*dens(itr12)+coul(itr34)
               coul(itr12)=bij
               bjk=exch(itr23)+gil*dens(itr14)
               bil=exch(itr14)+gil*dens(itr23)
               bik=exch(itr13)+gik*dens(itr24)
               exch(itr24)=exch(itr24)+gik*dens(itr13)
               exch(itr23)=bjk
               exch(itr14)=bil
               exch(itr13)=bik
240         continue
           else
            do 250 iq=1,ngout
            val =gout(iq)
            if( dabs(val).lt.cutoff) goto 250
               i = indsp(iq,ib1)+loci
               j = indsp(iq,jb1)+locj
                if(oianj.and.i.lt.j)  go to 250
               lab1q(iq) = iky(i) + j
               k = indsp(iq,kb1)+lock
               l = indsp(iq,lb1)+locl
                if(okanl.and.k.lt.l)  go to 250
               lab2q(iq) = iky(k) + l
                if(oident.and.lab1q(iq).lt.lab2q(iq)) go to 250
               ii1 = max(i,j)
               ii2 = min(i,j)
               ii3 = max(k,l)
               ii4 = min(k,l)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
               endif
               gik=val
               val2=val+val
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
               gil=val
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto  30
               itr23=ikyk+i2
               if(i2.ge.i4)goto  30
               itr24=iky(i4)+i2
   30          bij=val2*dens(itr34)+coul(itr12)
               coul(itr34)=val2*dens(itr12)+coul(itr34)
               coul(itr12)=bij
               bjk=exch(itr23)+gil*dens(itr14)
               bil=exch(itr14)+gil*dens(itr23)
               bik=exch(itr13)+gik*dens(itr24)
               exch(itr24)=exch(itr24)+gik*dens(itr13)
               exch(itr23)=bjk
               exch(itr14)=bil
               exch(itr13)=bik
250         continue
           endif
cdir$ vector
cdir$ nolist
      return
      end
      subroutine dir_build_open2_70(l2,coul,exch,dens,gout)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxprms2 = mxprms * mxprms)
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
INCLUDE(common/cslosc)
INCLUDE(common/ijlab)
INCLUDE(common/indez)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/misc)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/shlg70)
INCLUDE(common/shlnos)
INCLUDE(common/shlt)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      common/junk/pppp(2*mxprms2),inddd(mxprms2),indsp(256,4),
     -            lab1q(256),lab2q(256),iptge(256)
      common/ttqout/nok,nnok
      common/type  /itype,jtype
c
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/

c
c     ----- Direct open-shell SCF J & K builder (nshell > 1)
c
c
c ***
c *** vectorised for (sp sp/sp sp) blocks.
c *** assume that shells are in cannonical order and that
c *** that also implies cannonical order of basis functions.
c ***
c ***
cdir$ list
cdir$ novector
      ifac=256
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      if(jtype.ge.5 .or. jtype.eq.3) goto 230
c ***
      ijn = 0
      jmax = maxj
         do 110 i = 1,maxi
            i1=loci+i
            if (oianj) jmax = i
            do 100 j = 1,jmax
               ijn = ijn+1
               n1 = ib(ib1,i)+ib(jb1,j)+1
               lmax = maxl
               i2=locj+j
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               kln=0
               do 90 k = 1,maxk
                  if (okanl) lmax = k
                  i3=lock+k
                  do 80 l = 1,lmax
                     kln = kln+1
                     if (oident .and. kln .gt. ijn) go to 100
                     nn = n1+ib(kb1,k)+ib(lb1,l)
                     val = gout(nn)
                     if ( dabs(val) .lt. cutoff) go to 80
                     i4=locl+l
                     ii3 = max(i3,i4)
                     ii4 = min(i3,i4)
                     if(ii1.ge.ii3) then
                        i1 = ii1
                        i2 = ii2
                        i3 = ii3
                        i4 = ii4
                     else
                        i1 = ii3
                        i2 = ii4
                        i3 = ii1
                        i4 = ii2
                     endif
                     gik=val
                     val2=val+val
                     ikyi=iky(i1)
                     ikyj=iky(i2)
                     ikyk=iky(i3)
                     itr13=ikyi+i3
                     itr14=ikyi+i4
                     itr12=ikyi+i2
                     itr23=ikyj+i3
                     itr24=ikyj+i4
                     itr34=ikyk+i4
                     gil=val
                     if(i1.eq.i3.or.i2.eq.i4)gik=val2
                     if(i2.eq.i3)gil=val2
                     if(i2.ge.i3)goto 280
                     itr23=ikyk+i2
                     if(i2.ge.i4)goto 280
                     itr24=iky(i4)+i2
  280                continue
                     do 150 iiii=1,nsheld
                      bij=val2*dens(itr34)+coul(itr12)
                      coul(itr34)=val2*dens(itr12)+coul(itr34)
                      coul(itr12)=bij
                      bjk=exch(itr23)+gil*dens(itr14)
                      bil=exch(itr14)+gil*dens(itr23)
                      bik=exch(itr13)+gik*dens(itr24)
                      exch(itr24)=exch(itr24)+gik*dens(itr13)
                      exch(itr23)=bjk
                      exch(itr14)=bil
                      exch(itr13)=bik
                      itr12=itr12+l2
                      itr34=itr34+l2
                      itr13=itr13+l2
                      itr14=itr14+l2
                      itr23=itr23+l2
                      itr24=itr24+l2
150                   continue
 80               continue
 90            continue
100         continue
110      continue
      return
c *** for (sp sp / sp sp) only
230   if(jtype.eq.6) then
         ngout=256
      elseif(jtype.eq.5) then
            ngout=64
         elseif(jtype.eq.3) then
               ngout=16
            else
               call caserr(' confusion in dir_build_open2_70')
      endif
          if(.not.(oident.or.oianj.or.okanl)) then
            do 240 iq=1,ngout
            val = gout(iq)
            if( dabs(val).lt.cutoff) goto 240
               i1 = indsp(iq,ib1)+loci
               i2 = indsp(iq,jb1)+locj
               i3 = indsp(iq,kb1)+lock
               i4 = indsp(iq,lb1)+locl
               ii1 = max(i1,i2)
               ii2 = min(i1,i2)
               ii3 = max(i3,i4)
               ii4 = min(i3,i4)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
               endif
               gik=val
               val2=val+val
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
               gil=val
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto  20
               itr23=ikyk+i2
               if(i2.ge.i4)goto  20
               itr24=iky(i4)+i2
   20          continue
               do 200 iiii=1,nsheld
                bij=val2*dens(itr34)+coul(itr12)
                coul(itr34)=val2*dens(itr12)+coul(itr34)
                coul(itr12)=bij
                bjk=exch(itr23)+gil*dens(itr14)
                bil=exch(itr14)+gil*dens(itr23)
                bik=exch(itr13)+gik*dens(itr24)
                exch(itr24)=exch(itr24)+gik*dens(itr13)
                exch(itr23)=bjk
                exch(itr14)=bil
                exch(itr13)=bik
                itr12=itr12+l2
                itr34=itr34+l2
                itr13=itr13+l2
                itr14=itr14+l2
                itr23=itr23+l2
                itr24=itr24+l2
200            continue
240         continue
           else
            do 250 iq=1,ngout
            val =gout(iq)
            if( dabs(val).lt.cutoff) goto 250
               i = indsp(iq,ib1)+loci
               j = indsp(iq,jb1)+locj
                if(oianj.and.i.lt.j)  go to 250
               lab1q(iq) = iky(i) + j
               k = indsp(iq,kb1)+lock
               l = indsp(iq,lb1)+locl
                if(okanl.and.k.lt.l)  go to 250
               lab2q(iq) = iky(k) + l
                if(oident.and.lab1q(iq).lt.lab2q(iq)) go to 250
               ii1 = max(i,j)
               ii2 = min(i,j)
               ii3 = max(k,l)
               ii4 = min(k,l)
               if(ii1.ge.ii3) then
                  i1 = ii1
                  i2 = ii2
                  i3 = ii3
                  i4 = ii4
               else
                  i1 = ii3
                  i2 = ii4
                  i3 = ii1
                  i4 = ii2
               endif
               gik=val
               val2=val+val
               ikyi=iky(i1)
               ikyj=iky(i2)
               ikyk=iky(i3)
               itr13=ikyi+i3
               itr14=ikyi+i4
               itr12=ikyi+i2
               itr23=ikyj+i3
               itr24=ikyj+i4
               itr34=ikyk+i4
               gil=val
               if(i1.eq.i3.or.i2.eq.i4)gik=val2
               if(i2.eq.i3)gil=val2
               if(i2.ge.i3)goto  30
               itr23=ikyk+i2
               if(i2.ge.i4)goto  30
               itr24=iky(i4)+i2
   30          continue
               do 300 iiii=1,nsheld
                bij=val2*dens(itr34)+coul(itr12)
                coul(itr34)=val2*dens(itr12)+coul(itr34)
                coul(itr12)=bij
                bjk=exch(itr23)+gil*dens(itr14)
                bil=exch(itr14)+gil*dens(itr23)
                bik=exch(itr13)+gik*dens(itr24)
                exch(itr24)=exch(itr24)+gik*dens(itr13)
                exch(itr23)=bjk
                exch(itr14)=bil
                exch(itr13)=bik
                itr12=itr12+l2
                itr34=itr34+l2
                itr13=itr13+l2
                itr14=itr14+l2
                itr23=itr23+l2
                itr24=itr24+l2
300            continue
250         continue
           endif
cdir$ vector
cdir$ nolist
      return
      end
_ENDIF
      subroutine ver_integs(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integs.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
