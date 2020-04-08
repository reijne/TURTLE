c     deck=scf
c ******************************************************
c ******************************************************
c             =   scfgvb  =
c ******************************************************
c ******************************************************
      subroutine rhfgvb(q,ogrhf)
c
c     ----- grhf and gvb program driver -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/tran)
INCLUDE(common/scra7)
INCLUDE(common/timez)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/scfwfn)
INCLUDE(common/field)
INCLUDE(common/drfopt)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     +ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     +nsytyp,nsos
     +,tti(25),e1i(25),asym
     +,cilow(12),jpair,kone,ktwo,kcorb(2,12)
     +,iojk(49),ioham(25),iojkao(49)
      common/blkcore/corev(512),array(10)
_IF(drf)
cdrf
INCLUDE(../drf/comdrf/drfpar)
cdrf
_ENDIF
      character*10 charwall
      dimension q(*),o1e(6)
c
c     ----- set up length values -----
c
      l1 = num
      l2 = l1*(l1+1)/2
      l3 = l1*l1
c
c     ----- set up grhf-gvb data from /scfwfn/ -----
c
      call gvbset(nl2max,l1)
c
c     ----- allocate dynamic storage space -----
c
c     -e-     at  q(i10)
c     -trans- at  q(i20)
c     -q   -  at  q(i50)
c     -dens-   at  q(i60)
c     -xx,eig- at  q(i40)
c     -coul-  at  q(i70)
c     -exch-  at  q(i80)
c
c
      mc = 10*l2
c
c     first determine total memory available, and use
c     this in deriving data storage, no of passes of
c     fock build etc. Free allocated storage, and reasssign
c     with amount actually needed
c
      nmaxly  = igmem_max_memory()
      left = nmaxly - mc
      if (left .lt. 1) left = 0
      nwscmu = left + 6*l2
      nwscmu = (nwscmu/l2)*l2
      nwscmu = 2*((nwscmu/l2)/3) + ((nwscmu/l2)/3)
      if (nwscmu .ge. nl2max) then
        nwscmu = nl2max
        nl2max = 0
      endif
      nl2max = nl2max - nwscmu
      nl2max = 3*((nl2max+2)/3)
      nwscmu = nwscmu*l2
      nwscmu = max(nwscmu,6*l2)
      nwscm4 = ((nwscmu/l2)/3)*l2
      nwscm4 = max(nwscm4,l1+l3)
c
c    determine total requirements
c
      i10  = 0
      i20  = i10+l1
      i30  = i10+l2
      i40  = i30+l2
      i50  = i40+l1
      i60  = i50+l3
      i70  = i60+nwscm4
      i80  = i70+nwscm4
      last = i80+nwscm4
c
      if (ofield) then
         i51  = i40+l2
         i61  = i51+l2
         i71  = i61+l2
         lastf  = i71+l2
         last   = max(last,lastf)
      endif
c  
c     ----- get core memory -----
c
      i10 = igmem_alloc(last)
      i20  = i10+l1
      i30  = i10+l2
      i40  = i30+l2
      i50  = i40+l1
      i60  = i50+l3
      i70  = i60+nwscm4
      i80  = i70+nwscm4
      last = i80+nwscm4
c
      if (ofield) then
         i51  = i40+l2
         i61  = i51+l2
         i71  = i61+l2
         lastf  = i71+l2
         last   = max(last,lastf)
      endif
      if (nprint .eq. 5) write (iwr,9008)i10,i20,i30,i40,i50 ,i60 ,i70
     +,i80, last
c
c ----- first restore s,t,f from section 192 and store on ed7
c
c ----- are field values to be introduced
c
      if(ofield) then
c
      do loop = 1,6
       o1e(loop) = .true.
      enddo
      call getmat(q(i10),q(i30 ),q(i40),q(i51),q(i61),q(i71),
     +            array,num,o1e,ionsec)
      do 220 loop=1,l2
      q(i40+loop-1) = q(i40+loop-1) - fieldx * q(i51+loop-1)
     *                              - fieldy * q(i61+loop-1)
     *                              - fieldz * q(i71+loop-1)
 220  continue
      else
c
      do loop =1,3
       o1e(loop) = .true.
       o1e(loop+3) = .false.
      enddo
      call getmat(q(i10),q(i30 ),q(i40),q(i40),q(i40),q(i40),
     +            array,num,o1e,ionsec)
      endif
c
_IF(drf)
      if (oreact. and. field .eq. ' ') then
_ENDIF
      call wrt3(q(i10),l2,ibl7s,num8)
      call wrt3(q(i30),l2,ibl7t,num8)
      call wrt3(q(i40),l2,ibl7f,num8)
c
c ----- transform s-matrix
c
      call tranp(q(i10),q(i30 ))
      call wrt3(q(i30),l2,ibl7st,num8)
_IF(drf)
      endif
_ENDIF
      call gvbitr(q,q(i20),q(i50),q(i60),q(i70),q(i80),q(i10),
     +            q(i40),nwscmu,l1,l2,l3)
c
c     ----- reset memory -----
c
      call gmem_free(i10)
c
      if(nprint.eq.-5)return
      cpu=cpulft(1)
      write(iwr,8888)cpu,charwall()
 8888 format(//' end of grhf-gvb scf at ',f10.2,' seconds',a10,' wall'
     *,//1x,104('-')/)
c
      return
 9008 format(' core assignment'/'    i10,   i20,   i30,   i40,   i50,'
     +     , '   i60,   i70,   i80'      /8i7 /' last = ',
     +     i7)
      end
      subroutine gvbset(nl2max,l1)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/atmol3)
INCLUDE(common/scfwfn)
INCLUDE(common/infoa)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
c
      if (nconv .le. 0) nconv = 5
      norb = 0
      nsos = 0
c
c     ----- set up some option switches -----
c
      nsytyp = 1
      icoupl = 0
      if (ifolow .ne. 1) ifolow = 0
      irotb = 1
c
c     ----- set up orbital description -----
c
      if(gapa1.eq.0.0d0.and.gapa2.eq.0.0d0) irotb=-1
      norb = nco
      ncores=1
      if(nco.eq.0) then
        ncores=0
      endif
      nopen = 0
      if (nseto .ne. 0) then
         do 100 i = 1,nseto
         nop = no(i)
         nopen = nopen+nop
         norb = norb+nop
  100    continue
      endif
      if (npair .ne. 0) then
         norb = norb+npair+npair
      endif
c
c     ----- calculate the number of fock operators nham -----
c
      nham = ncores+npair+npair+nseto
c
c     ----- calculate the maximum number of triangular matrices
c     which fill be needed if all is to fit into core -----
c
      nl2max = 3*ncores + 3*(npair+npair+nseto)
      nl2max = max(6,nl2max)
c
c     ----- set up symmetry - only one type currently -----
c
      jorb = 0
      ns(1) = l1
      melsym(1) = norb
      do 180 i = 1,nsytyp
      ilo = nsos+1
      jhi = nsos+melsym(i)
      nsos = nsos+ns(i)
      do 160 j = ilo,jhi
      jorb = jorb+1
  160 msympr(jorb) = j
  180 continue
c
c     ----- generate nconf -----
c
      if (nco .ne. 0) then
      call setsto(nco,1,nconf)
      endif
      if (nseto .ne. 0) then
         nbase = ncores
         ic = 0
         do 260 i = 1,nseto
         nop = no(i)
         call setsto(nop,nbase+1,nconf(ic+nco+1))
         ic = ic+nop
  260    nbase = nbase+1
      endif
      if (npair .ne. 0) then
         np2 = npair+npair
         do 300 i = 1,np2
 300     nconf(i+nco+nopen) = ncores+nseto+i
      endif
      nvirt=nsos-norb
      if(nvirt.ne.0) then
         call setsto(nvirt,nham+1,nconf(norb+1))
      endif
c
c     ----- generate kcorb -----
c
      nbase = nco+nopen
      if (npair .ne. 0) then
         do 340 kpair = 1,npair
         kcorb(1,kpair) = nbase+1
         kcorb(2,kpair) = nbase+2
         nbase = nbase+2
  340    continue
      endif
      return
      end
      subroutine gvbout(trans,gi,e,l1,nprint)
c
c     ----- grhf-gvb output program -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/prints)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/scfopt)
INCLUDE(common/scfwfn)
INCLUDE(common/dm)
INCLUDE(common/harmon)
      character*10 charwall
      dimension trans(l1,*),gi(l1,*),e(*)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      data dzero,done/0.0d0,1.0d0/
      data m20/20/
      data idzero/0/
c
      mpunch = 1
      lprnt = l1
      if(.not.oprint(20)) lprnt = min(norb+5,l1)
      lprnt=min(lprnt,newbas0)
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1
      if (nprint .ne. 7 .and. mpunch .ne. 1) go to 220
      if (nprint .ne. 7) go to 200
c
c     ----- punch out restart data -----
c
      ifab = 3
      write (ipu,9008) nco,nseto,npair,ifab
      if (nseto .gt. 0) then
       do 100 i = 1,nseto
       write (ipu,9008) no(i)
  100  continue
      endif
      if (npair .gt. 0) then
       do 140 kpair = 1,npair
       write (ipu,9028) (cicoef(j,kpair),j = 1,2)
  140  continue
      endif
      write (ipu,9028) (f(i),i = 1,nham)
      ij = 0
      do 180 i = 1,nham
      write (ipu,9028) (alpha(ij+j),j = 1,i)
      write (ipu,9028) (beta(ij+j),j = 1,i)
      ij = ij + i
  180 continue
      write (ipu,9008) norb,idzero,idzero
  200 continue
      call pusql(trans,l1,l1,l1)
  220 continue
c
c     ----- print final results -----
c
      cpu=cpulft(1)
      write (iwr,9168) iter,cpu,charwall(),ek,ehf,en,etot,vir
      if (nprint .eq. -5) go to 420
      if(.not.oprint(50).or.oprint(25))go to 430
      write (iwr,9048)
      do 240 i = 1,norb
      write (iwr,9068) i,e(i)
  240 continue
      write (iwr,9088)
      call tdown(trans,ilifq,trans,ilifq,newbas0)
      call prev(trans,e,lprnt,l1,l1)
 430  if (npair .eq. 0) go to 280
      write (iwr,9128)
      write (iwr,9108)
      do 260 i = 1,npair
      sab = (cicoef(1,i)+cicoef(2,i))/(cicoef(1,i)-cicoef(2,i))
      sab =  dabs(sab)
      write (iwr,9148) i,(kcorb(j,i),j = 1,2),(cicoef(j,i),j = 1,2),
     +     cilow(i),sab
  260 continue
  280 if(.not.oprint(50).or.oprint(25))go to 420
      if (npair .eq. 0) go to 380
      write (iwr,9188)
      do 360 i = 1,npair
      c1 = cicoef(1,i)
      c2 = cicoef(2,i)
      m1 = kcorb(1,i)
      m2 = kcorb(2,i)
      m1 = msympr(m1)
      m2 = msympr(m2)
      alph = dsqrt(-c1/c2)
      bet = done/dsqrt(done+alph**2)
      alph = alph*bet
      rmax1 = dzero
      rmax2 = dzero
      do 320 k = 1,l1
      gi(k,1) = alph*trans(k,m1) + bet*trans(k,m2)
      gi(k,2) = -alph*trans(k,m1)+bet*trans(k,m2)
      if (  dabs(gi(k,1)) .lt. rmax1) go to 300
      iph1 = 1
      rmax1 =  dabs(gi(k,1))
      if (gi(k,1) .lt. dzero) iph1 = -1
  300 continue
      if ( dabs(gi(k,2)) .lt. rmax2) go to 320
      iph2 = 1
      rmax2 =  dabs(gi(k,2))
      if (gi(k,2) .lt. dzero) iph2 = -1
  320 continue
      do 340 k = 1,l1
      if (iph1 .gt. 0) go to 340
      gi(k,1) = -gi(k,1)
      if (iph2 .gt. 0) go to 340
      gi(k,2) = -gi(k,2)
  340 continue
      write (iwr,9208) i
      call prsql(gi,2,l1,l1)
  360 continue
c
c     ----- print out the final f coefficient matrix -----
c
  380 write (iwr,9228)
      do 400 i = 1,nham
      write (iwr,9248) i,f(i)
  400 continue
c
c     ----- get the lagrangian multipliers  -----
c
      l2norb = norb*norb
      call secget(iseclg,m20,iblk20)
      call rdedx(gi,l2norb,iblk20,idaf)
      write (iwr,9268)
      call prsq(gi,norb,norb,norb)
  420 write (iwr,9288) asym
      return
 9008 format(14i5)
 9028 format(5e15.8)
 9048 format(/10x,11('-')/10x,'eigenvalues'/10x,11('-')/)
 9068 format(1x,i5,f20.12,4x,f20.12)
 9088 format(/10x,12('-')/10x,'eigenvectors'/10x,12('-')/)
 9108 format(10x,' pair  orbital 1  orbital 2   ci coef 1   ci coef 2',
     +3x,'ci lowering  overlap'/)
 9128 format(/10x,16('-')/10x,'pair information'/10x,16('-')/)
 9148 format(10x,i4,i8,i11,5x,f10.7,2x,f10.7,2x,f10.7,1x,f10.7)
 9168 format(/10x,14('-')/10x,'final energies after',i4,' cycles at '
     *,f12.2,' seconds ',a10,' wall'/
     *10x,14('-')/
     +10x,'kinetic energy    ',f18.10/ 
     +10x,'electronic energy ',f18.10/ 
     +10x,'nuclear energy    ',f18.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy      ',f18.10,8x,f10.7/)
 9188 format(/10x,11('-')/10x,'gi orbitals'/10x,11('-')/)
 9208 format(//20x,'pair ',i3)
 9228 format(/10x,20('-')/10x,'f coefficient matrix'/10x,20('-')/)
 9248 format(3x,5x,i5,f20.12)
 9268 format(/10x,22('-')/10x,'lagrangian multipliers'/10x,22('-')/)
 9288 format(/10x,'max lagrangian asymmetry',2x,f15.10)
      end
      subroutine redden(trans,da,db,nconf,l1,l2)
c
c     ----- calculate the reduced density matrix -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scfwfn)
      dimension trans(l1,*),da(l2),db(l2)
      dimension nconf(*)
      data dzero,two/0.0d0,-2.0d0/
c
c     ----- set up some control data needed to generate the data -----
c
      if (nco .ne. 0) call setsto(nco,1,nconf)
      nopen = 0
      iconf = ncores+1
      iorb = nco
      if (nseto .ne. 0) then
      do 160 i = 1,nseto
      nopen = nopen+no(i)
      nop = no(i)
      call setsto(nop,iconf,nconf(iorb+1))
      iorb = iorb + nop
  160 iconf = iconf + 1
      endif
      if (npair .ne. 0) then
      npair2 = npair + npair
      do 200 i = 1,npair2
      nconf(i+iorb) = iconf
  200 iconf = iconf + 1
      endif
      norb = nco + nopen + npair + npair
c
      call vclr(da,1,l2)
      call vclr(db,1,l2)
c
c
c     ----- get alpha part first -----
c
c     the open part of alpha will be off by 0.5
c
      iadd=1
      do 280 i=1,l1
      do 260 k = 1,norb
      dum=f(nconf(k))*trans(i,k)
      if(dum.ne.dzero) then
      call daxpy(i,dum,trans(1,k),1,da(iadd),1)
         endif
  260 continue
c
c     ----- now get beta part - first get only the open part -----
c
  280 iadd=iadd+i
      ilo = nco+1
      ihi = nco+nopen
      if (nopen .gt. 0) then
      iadd=1
      do 300 i=1,l1
      do 320 k = ilo,ihi
      dum=f(nconf(k))*trans(i,k)
      if(dum.ne.dzero) then
      call daxpy(i,dum,trans(1,k),1,db(iadd),1)
         endif
  320 continue
 300  iadd=iadd+i
      endif
c
c     ----- now add this into alpha to get the correct alpha
c           density -----
c
      call vadd(da,1,db,1,da,1,l2)
c
c     ----- now subtract out two times what's in db to give the
c           correct beta -----
c
      call vsma(db,1,two,da,1,db,1,l2)
c
c     ----- save the reduced density matrices -----
c
      call wrt3(da,l2,ibl3pa,idaf)
      call wrt3(db,l2,ibl3pb,idaf)
      return
      end
      subroutine gvbitr(core,vect,temp,dens,coul,exch,e,pp,
     +                  nwscmu,l1,l2,l3)
c
c     ----- main grhf-gvb driver -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/drfopt)
_IF(drf)
cafc
INCLUDE(../drf/comdrf/sizesrf)
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfbem)
      integer idafh, navh
      common/hdafile/idafh,navh,ioda(2,1000)
      common/enrghf/enhh,etothh,ehfhh
c
_ENDIF
      dimension vect(l1,*),temp(l1,*),dens(*),coul(*),exch(*)
      dimension core(*),e(*)
      dimension pp(*)
INCLUDE(common/scra7)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scfwfn)
INCLUDE(common/dm)
INCLUDE(common/timez)
INCLUDE(common/restar)
INCLUDE(common/scfopt)
INCLUDE(common/infoa)
INCLUDE(common/harmon)
INCLUDE(common/zorac)
INCLUDE(common/atmol3)
INCLUDE(common/machin)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/runlab)
INCLUDE(common/timeperiods)
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
      common/diisd/
     +        st(210),ccc(20),rrr(19),derror,scale(20),iposit(20),
     +        nstore,mp,ondiis,nspaca,
     +        stb(210),cccb(20),rrrb(19),derrb,scaleb(20),iposb(24),
     +        ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     +        nsytyp,nsos,
     +        tti(25),e1i(25),asym,
     +        cilow(12),jpair,kone,ktwo,kcorb(2,12),
     +        iojk(49),ioham(25),iojkao(49)
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c           run through all code section
c
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
INCLUDE(common/fsymas)
c ... for dummy symass section
      parameter(isymtp=99)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF

      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      data zscf/'scf'/
      data done,two,ten/1.0d0,2.0d0,1.0d+01/
      data pt2,twopt2/0.2d0,2.2d0/
      data m23,m20/23,20/
      data dmptlc/1.0d-02/
      data oredo / .false./
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      call check_feature('gvbitr')
      out = nprint .eq. 5
      outon = nprint .ne. -5
_IF(parallel)
c
c - this should be default anyway -
c
      nav = lenwrd()
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
*     write(6,*)'master',ipg_nodeid(),omaster
_ENDIF
      call start_time_period(TP_TEST1)
      skale=2.0d0
      lprnt = l1
      if(.not.oprint(20))lprnt = min(norb+5,l1)
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9008)
_IF(ccpdft)
      if(CD_active())then
       call caserr('DFT not available for GRHF-GVB, use UHF')
      endif
_ENDIF
      tim=0.0d0
      en=enucf(nat,czan,c)
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(vect,l2,ibl7s,num8)
         call rdedx(dens,l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,vect,dens)
      endif
_ELSE
      call rdedx(dens,l2,ibl7st,num8)
_ENDIF
      call rdedx(vect,l3,ibl3qa,idaf)
      call vclr(e,1,norb)
      if(out) then
         call prev(vect,e,lprnt,l1,l1)
         call prtri(dens,l1)
      endif
c
c
c     ----- initialize individual timers for the five sections
c           of grhf-gvb-----
c
      timed = 0.0d0
      timef = 0.0d0
      timit = 0.0d0
      timeg = 0.0d0
      timem = 0.0d0
      timeo = 0.0d0
c
c     ---- initialize some constants for convergence control -----
c
c
      iprint=nprint
      if (nconv .le. 0) nconv = 5
      acurcy = ten**(-nconv)
      damp = 0.0d0
      damp0 = 0.0d0
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0
c     sqcdf = 0.0d0
      iter = 0
      iterv = 0
      kcount = 0
      rshift = 0.0d0
      diff = 0.0d0
      diffo = diff
      de = 0.0d0
      dep = 0.0d0
      deavg = 0.0d0
      ehf = 0.0d0
      lockt = lock
      ehf0 = 0.0d0
      ek = 0.0d0
      vir = 0.0d0
c
c     ----- initialize disk for extrapolation -----
c
      ndaf = ibl7la
      lnoc = l1*norb
      if (lnoc .lt. l2) lnoc = l2
      lene=lensec(lnoc)
      lene=lene*3
      lenq=lensec(l3)
      ibl7qa=ndaf+lene
      iblko=ibl7qa+lenq
      lenp=lensec(l2)
      ibl7la=iblko+lenp
c
c ----- allocate space on ed7 for j,k,ham
c
      lenh=lensec(l2)
      iojk(1)=ibl7la
      ij=ibl7la+lenh
      nham1=nham-ncores
      if(nham1.lt.1)go to 3001
      do 3000 i=1,nham1
      ii=i+i
      iojk(ii)=ij
      iojk(ii+1)=ij+lenh
 3000 ij=ij+lenh+lenh
 3001 iojkao(1)=ij
      ij=ij+lenh
      if(nham1.lt.1) goto 3004
      do 3003 i=1,nham1
      ii=i+i
      iojkao(ii)=ij
      iojkao(ii+1)=ij+lenh
 3003 ij=ij+lenh+lenh
 3004 continue
      do 3002 i=1,nham
      ioham(i)=ij
 3002 ij=ij+lenh
      ndafd=ij
      ocvged = .false.
c     orotdm = .false.
c     if (mconv .gt. 7) orotdm = .true.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if(odiis) go to 500
      if (odamph .or. oshift) damp = dmax1(done,dmpcut)
500   call wdisk(vect,vect,vect,ndaf,num8,lnoc)
c
c     ----- compute a set of orthonormal orbitals from the
c           overlap matrix -----
c
      l0 = newbas0
      call qmat(core,dens,coul,e,exch,iky,l0,l1,l3,l1,out)
      lprnt = min(lprnt,l0)
      l0 = l1
      ltran = l0
c
c....  the use of l0 in this routine is confusing
c....  it basically is l1; newbas0 has taken it's role
c....  see e.g. ortho1
c....  (hvd): this is explained in more detail in subroutine preharm
c
c     ----- get initial time remaining for time check purposes -----
c
      call timrem(tlefts)
      call end_time_period(TP_TEST1)
_IF(drf)
cafc
c  write start_density to da31
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        call dawrit(idafdrf,iodadrf,dens,l2*2,70,navdrf)
      endif
_ENDIF
c
c     ----- if iter .eq. 0 , skip the convergence routines -----
c
      timit=cpulft(1)
 100  continue
      call start_time_period(TP_TEST2)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
      if (iter .gt. 0) then
        if(iter.eq.maxcyc-1) nprint=15
c
c     ----- damp and extrapolate the eigenvectors if necessary -----
c
        if (iter .eq. 1) deavg = 0.0d0
        if (iter .eq. 2) deavg =  dabs(de)
        if (iter .ge. 3) deavg = ( dabs(de)+ dabs(dep)+pt2*deavg)/
     +                             twopt2
        dmptst = acurcy
        if ((diffo-diff) .lt. 0.0d0) dmptst = diffo-diff
        if (iter .gt. 2) call dampd(de,dep,deavg,damp,dmptst,diff,diffp,
     +       dmptlc)
        if (damp .lt. dmpcut) damp = dmpcut
        call extrpd(de,damp,damp0,vect,dens,coul,exch,
     +       l1,lnoc,ndaf,num8,iterv,1,2)
        ltran = newbas0
      else
        ltran = norb
      endif
c
c     ----- reorthonormalize if damping, extrapolation or
c           level shifting was done -----
c
      if(iter.eq.0)go to 142
      if(mod(iter,5).eq.0) goto 142
      if(.not.odiis.and. kcount.eq.0) goto 142
      if(damp.ne.0.0d0) goto 142
      goto 160
 142  call ortho1(dens,dens,vect,pp,iky,ltran,l0,l1,l2,l3,l1)
      ltran=l0
      call rdedx(dens,l3,ibl3qs,idaf)
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,vect,l1)
         call load_ga_from_square(ih_vec,dens,l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(vect,ih_scr2, l1)

      else
         call tfsqc(vect,dens,coul,l0,l1,l1)
      endif
_ELSE
      call tfsqc(vect,dens,coul,l0,l1,l1)
_ENDIF
 160  if (out) call prev(vect,e,lprnt,l1,l1)
c
c     ----- update the contents of the transformation vector
c           file -----
c
      call wrt3(vect,l3,ibl3qa,idaf)
      time0=cpulft(1)
      call tdown(vect,ilifq,vect,ilifq,ltran)
c
c     ----- set up the data for the virial calculation -----
c
      call virset(vect,dens,dens(l2+1),l1,l2,num8,core)
_IF(parallel)
      else
c
c...  other roots ** idle **
c...
      endif
c...
      if(ipiomode .eq. IO_NZ_S)then
        call pg_brdcst(7124,ltran,8,0)
        call pg_brdcst(7123,vect,l3*8,0)
        call pg_synch(4444)
        oswed3(4)=.true.
        oswed3(8)=.true.
      endif
      call pg_synch(5555)
c...
_ENDIF
      call end_time_period(TP_TEST2)
c
c calculate zora corrections if need be;
c here it is a dummy add to core(iscr+l2)
c they are added to one electron fock in jkform
c
      if (ozora .and. iter .ne. 0) then
         iscr = igmem_alloc(l2*2)
         call    rdedx(core(iscr),l2,iblko,num8)
         itert = iter
         if (ocvged) itert = 999999999
         if ((mod(iter,niter_z).ne.0.and.itert.ne.999999999)
     1       .or.nat_z.ne.0)  then
            call zora(core,core(iscr),core(iscr+l2),'calc')
         else
            call zora(core,core(iscr),core(iscr+l2),'force')
         end if
         call gmem_free(iscr)
      endif
c
c     ----- form the j and k fock matrices -----
c
      call jkform(vect,temp,dens,coul,exch,nwscmu,nopk,
     +     npass,l1,l1,l2,nprint,core)
      timdum=cpulft(1)
      timef = timef +timdum - timit
      timit = timdum
c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
c
c     ----- determine the optimal ci coefficients -----
c
      call optci(coul,exch,pp,iky,l1,l2,nprint)
      timdum=cpulft(1)
      timeg = timeg + timdum - timit
      timit = timdum
c
c ------ determine gvb orbital gradients
c
      diffpp=diffp
      diffp=diff
      diff=0.0d0
      call gvbgrd(temp,coul,coul(l2+1),iky,l2)
c
c ------ monitor convergence here
c
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
c
cgdf  not sure if needed: convergence check here seems to be for 
c     skipping some work, when the real check comes later
c
      if (osrso) ocvged = .true.
c
      if (ocvged)  then
       oredo = .false.
       if (ozora) then
c...         see if our coulomb matrix is up to date
          if (mod(iter-1,niter_z).ne.0.and.itert.ne.999999999.and.
     1        .not.(is_z .ne. 0.or.icoul_z .eq. 3.or.nat_z.ne.0)) then
             if (opzora) write(iwr,*) ' ** restart coulomb calc **'
c...          reset diis   ...
             nstore = 0
             oredo = .true.
             go to 667
          end if
        endif
       end if 
667    continue
c
c ---- additional code for diis2 method
c
      if(.not.odiis) goto 180
c
c ------ take precautions if tester is increasing again 
c
      if(diff.gt.diffp.and.diffp.gt.diffpp) then
*   the line below which is intended to ignore the impact
*   of rising tester on diis does NOT yield the correct
*   energies .. to be resolved .. (1/10/95)
*      if(diff.ge.acurcy*1000.0d0.or.iter.eq.0) then
        nstore=0
        goto 180
*      endif
      endif
      if (ocvged .and. .not. oredo) go to 270
      call start_time_period(TP_DIIS)
      call diiso(temp,dens,coul,exch,vect,ndafd,lockt,num8)
      call end_time_period(TP_DIIS)
c
c ------ read in symmetry adapted vectors
c
      call rdedx(vect,l3,ibl3qa,idaf)
      if (ondiis) then
        call start_time_period(TP_ORFOG)
        njkmat=2*(nham-ncores)+1
        call rdedx(dens,l2,ibl7st,num8)
_IF(ga)
        if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
           call porth(vect,dens, num, newbas0)
        else
           call mult2(vect,coul,dens,l0,l0,l1)
           call orfog(vect,vect,coul,coul(l2+1),iky,ilifq,newbas0,l1,1)
        endif
_ELSE
        call mult2(vect,coul,dens,l0,l0,l1)
        call orfog(vect,vect,coul,coul(l2+1),iky,ilifq,newbas0,l1,1)
_ENDIF
        call wrt3(vect,l3,ibl3qa,idaf)
        call tdown(vect,ilifq,vect,ilifq,newbas0)
        call end_time_period(TP_ORFOG)
c
c ------ transform diis generated jk matrices to mo basis
c
        call start_time_period(TP_TEST3)
        do 163 jk=1,njkmat
        call rdedx(coul,l2,iojkao(jk),num8)
_IF(ga)
        if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
           call load_ga_from_square(ih_vec,vect,l1)
           call load_ga_from_triangle(ih_scr2,coul,l1)
           call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
           call load_triangle_from_ga(dens,ih_scr2, l1)
        else
           call mult2(vect,dens,coul,l1,l1,l1)
        endif
_ELSE
        call mult2(vect,dens,coul,l1,l1,l1)
_ENDIF
 163    call wrt3(dens,l2,iojk(jk)  ,num8)
c
c     ----- determine the optimal ci coefficients -----
c
        call optci(coul,exch,pp,iky,l1,l2,nprint)
        call end_time_period(TP_TEST3)
      else
        call tdown(vect,ilifq,vect,ilifq,l1)
      endif
 180  timdum=cpulft(1)
      timed = timed + timdum - timit
      timit = timdum
c
c     ----- perform orbital mixing among occupied orbitals -----
c
c     ----- do the two by two rotations among the occupied orbitals ----
c
      call start_time_period(TP_TEST4)
      call rdedx(vect,l3,ibl3qa,idaf)
      call rotb(dens,vect,coul(l2+1),coul,iky,l1,l2,gapa1,ibrk,gapa2,
     +          nprint)
      timdum=cpulft(1)
      timem = timem + timdum - timit
      timit = timdum
c
c      ----- compute the fock matrices in the sab basis
c
      call gvbham(dens,dens(l2+1),coul,nwscmu,l2)
_IF(drf)
cafc
c  add drf contributions to fock matrix: here ???
c  last number in argument-list has got to be the variable
c  concerning number of shells: which ???
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        print *, 'calling drfhamo'
caleko
c        call drfhamo(q(i200),q(i210),q(i110),q(i120),q(i130),
c     +         q(i140),q(i170),q(i150),q(i180),q(i160),q(i190),2)
c        do 6666 i = 1, l2*2
c          coul(i) = coul(i) + scffact*q(i200-1+i)
c NOTE
c  the following CANNOT work, as none of the pointers (i200 etc)
c  have been set .. caserr for now
c
c       call drfhamo(core(i200),core(i210),core(i110),core(i120),
c    +         core(i130),core(i140),core(i170),core(i150),core(i180),
c    +         core(i160),core(i190),2)
        call caserr(' integer pointers must be set in gvbitr')
c       do 6666 i = 1, l2*2
c         coul(i) = coul(i) + scffact*core(i200-1+i)
c6666    continue
      endif
_ENDIF
      call end_time_period(TP_TEST4)
c
c     ----- perform the ocbse step - mix in virtuals -----
c
      call search(ioham(1),num8)
      ncall=1
      do 260 iham = 1,nham
      call reads(coul,l2,num8)
      if (out) call prtri(coul,l1)
      call start_time_period(TP_DIAG)
      call ocbse(iham,vect,dens,coul,temp,pp,e,de,dep,iterv,ncall,
     +           l1,l2)
      call end_time_period(TP_DIAG)
      ncall=0
  260 continue
      call wrt3(vect,l3,ibl3qa,idaf)
      call tdown(vect,ilifq,vect,ilifq,l1)
 270  continue
      timdum=cpulft(1)
      timeo = timeo + timdum - timit
      timit = timdum
c
c     ----- calculate the reduced density matrices -----
c
      call redden(vect,dens,dens(l2+1),nconf,l1,l2)
c
c     ----- combine the alpha and beta parts into one density
c           matrix. then get the old density matrix and save the
c           new density matrix -----
c
      call vadd(dens,1,dens(l2+1),1,dens,1,l2)
      if (iter .gt. 0) call rdedx(dens(l2+1),l2,iblko,num8)
      call wrt3(dens,l2,iblko,num8)
_IF(drf)
cahv - In the analysis we need the total density (a+b)
c      this is calculated here and needs to be written to dafile
c      The write in hfscf (master.m) is removed because it wrote
c      only the a density
cahv 0206
       if (oreact) then
          call dawrit(idafh,ioda,dens,l2,16,navh)
       endif
cafc
c  diff dens calculation for drf
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        diff2 = cvgden(dens,dens(l2+1),l2)
        call dawrit(idafdrf,iodadrf,dens(l2+1),l2,70,navdrf)
      endif
c
_ENDIF
      call rdedx(vect,l3,ibl3qa,idaf)
      if (out) call prev(vect,e,norb,l1,l1)
      time1=cpulft(1)
      shift2=rshift
      shift1=gapa1
      if(iter.gt.ibrk) shift1=gapa2
c     if (iter .gt. maxcyc) go to 340
      ek = ek+ek
      iter = iter + 1
      dep = de
      de = ehf - ehf0
      ehf0 = ehf
      etot = ehf+en
      delt = time1 - time0
      if (iter.eq.1.and.outon) then
        write(iwr,9048) maxcyc,mconv,nconv,npunch,en
        if (odebug(31)) then
         write (iwr,9028)
        else
         write (iwr,9029)
        endif
      endif
      diffo = diff
      vir = (etot-ek)/(two*etot)
      if(outon.or.(maxcyc-iter.lt.10)) then
       if (odebug(31)) then
        write (iwr,9068) iter,kcount,etot,ehf,de,diff,shift1,shift2,
     +   damp,derror,delt,time1
       else
        write (iwr,9069) iter,kcount,etot,ehf,de,diff,shift1,shift2,
     +   damp,derror
       endif
      endif
_IF(parallel)
      else
          iter = iter + 1
      endif
      if(ipiomode .eq. IO_NZ_S)then
        lscf = 20+12/nav
        call pg_brdcst(7125,maxcyc,lscf*8,0)
        call pg_brdcst(7126,vect,l3*8,0)
        call pg_synch(4444)
        oswed3(4)=.true.
        oswed3(8)=.true.
      endif
      call pg_synch(5555)
_ENDIF
c
c ------ monitor convergence here
c
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
cgdf
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if
c
      if (ocvged)  then
       oredo = .false.
       if (ozora) then
c...         see if our coulomb matrix is up to date
          if (mod(iter-1,niter_z).ne.0.and.itert.ne.999999999.and.
     1        .not.(is_z .ne. 0.or.icoul_z .eq. 3.or.nat_z.ne.0)) then
             if (opzora) write(iwr,*) ' ** restart coulomb calc **'
c...          reset diis ...
             nstore = 0
             oredo = .true.
             go to 666
          else if (oscalz) then
c            print *,'WARNING scaling not yet implemented'
             iscal_z = igmem_alloc(l2)
c sf compute scaling factors
             call zora(core,dens,core(iscal_z),'scale')
c sf perform scaling
             id_z = igmem_alloc(l2)
             call scale_z(core,e,core(id_z),vect,core(iscal_z),
     1                    ehf,etot,ne,num,0)
             call gmem_free(id_z)
             call gmem_free(iscal_z)
             goto 666
          end if
        endif
       end if 
666    continue
      if (.not.ocvged .or.oredo) then
c
c     ----- check for time remaining -----
c
      call timrem(tlefti)
      if (tlefti .gt. skale*(tlefts-tlefti)/iter) then
c       sufficient time to usefully continue
          if (iter .lt. maxcyc) go to 100
          write (iwr,9108)
          etot = 0.0d0
          irest=3
          ehf0 =ehf
          ehf = -en
        else
c     insufficient time to usefully continue
          write (iwr,9088)
          irest=3
          call texit(0,irest)
          etot = 0.0d0
          ehf = -en
        endif
      else
c
c     ----- energy converged -----
c
      if(outon)write (iwr,9128)
c
      endif
c
c     ----- calculate the five time steps per iteration -----
c
  360 continue
_IF(drf)
      if (field .ne. ' ') then
        ehfhh = ehf
        enhh = en
        etothh = etot
        call dawrit(idafh,ioda,dens,l2,16,navh)
      endif
_ENDIF
      timed = timed/ dfloat(iter)
      timef = timef/ dfloat(iter)
      timeg = timeg/ dfloat(iter)
      timem = timem/ dfloat(iter)
      timeo = timeo/ dfloat(iter)
      if(outon)write(iwr,9148)npass,timef,timeg,timem,
     +                        timeo,timed
c
c     ----- save the final /scfwfn/ -----
c
      m4=lensec(mach(5))
      i =m4+m4
      call secput(isect(504),m23,i,iblk23)
      iblk23 = iblk23 + m4
      call wrt3(cicoef,mach(5),iblk23,idaf)
c
c     ----- get the true eigenvalues -----
c
      do 400 i = 1,nham
      do 380 j = 1,norb
      if (nconf(j) .ne. i) go to 380
      e(j) = two*e(j)*f(i)
  380 continue
  400 continue
c
      l2norb = norb*norb
      len2=lensec(l2norb)
      call secput(iseclg,m20,len2,iblk20)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
c
c     ----- save the final eigenvalues and orbitals -----
c
      call start_time_period(TP_TEST5)
      call gvbsav(vect,e,mouta,l1,ibl3ea)
      call rdedx(vect,l3,ibl3qa,idaf)
      call lagrng(l1,norb,l2,coul,dens,vect,iky)
c
c     ----- save the lagrangian multpliers in section iseclg -----
c
      call wrt3(dens,l2norb,iblk20,idaf)
      if(ocvged)irest=0
      nprint=iprint
      call rdedx(vect,l3,ibl3qa,idaf)
c
_IF(drf)
cafc
      if (field .ne. ' ') then
        call dawrit(idafh,ioda,vect,l3,15,navh)
      endif
cafc
_ENDIF
      call gvbout(vect,temp,e,l1,nprint)
c
c     canonicalise gvb mos
c
      call canon(core,vect,dens,coul,e,l0,l1,l2,l3,nprint)
      call end_time_period(TP_TEST5)
_IF(parallel)
      else
c...  other roots
c...  allow for sections created in canon (moutb) and symass
c...  save results of the symmetry assignment into ed3
c...  assuming this has been called
      len1=lensec(mach(8))
      len2=lensec(mach(9))
      j=len2+lenq+len1+1
      call secput(moutb,3,j,iblnum)
       if (otsym) then
        isymsc = isect(499)
        call secput(isymsc,isymtp,1+lensec((num-1)/nav+1)
     +                             +lensec(num),iblnum)
       endif
      endif
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
c
      call copq_again(mouta)
      return
c
 9008 format(//40x,24('*')/ 40x,
     *'grhf-gvb scf calculation'/40x,24('*'))
 9028 format(/1x,127('=')/
     *3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',4x,'level shift',2x,'damping',11x,'diis',
     * 3x,'del(t)',5x,'time'/
     *18x,'energy',10x,'energy'/1x,127('='))
 9029 format(/1x,109('=')/
     *3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',4x,'level shift',2x,'damping',11x,'diis'/
     *18x,'energy',10x,'energy'/1x,109('='))
 9068 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,3f8.3,f15.8,2f9.2)
 9069 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,3f8.3,f15.8)
 9048 format(/15x,'convergence data'/15x,16('=')//
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criteria         =1.0e-',i2/
     +     ' punch out option             = ',i6//
     +    ' ----- nuclear energy ----- = ',f20.12/)
 9088 format(//10x,26('*')//
     *10x,'*** warning ***'/
     *10x,'scf has not converged yet'/
     *10x,'this job must be restarted'/
     *10x,26('*')//)
 9108 format(/10x,30('-')/10x,'excessive number of iterations' /
     +        10x,30('-'))
 9128 format(/10x,16('-')/10x,'energy converged'/10x,16('-')/)
 9148 format(/
     +  10x,28('-')/
     +  10x,'scf statistics per iteration'/
     +  10x,28('-')/
     +  10x,'number of integral passes   ',i5,/
     +  10x,'j+k  formation time    ',f10.3/
     +  10x,'geminal opt time       ',f10.3/
     +  10x,'mixorb  opt time       ',f10.3/
     +  10x,'ocbse   opt time       ',f10.3/
     +  10x,'orbital grad. and diis ',f10.3)
      end
      subroutine gvbsav(q,e,ndaf,l1,iblk)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),e(*)
INCLUDE(common/sizes)
INCLUDE(common/runlab)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
      data m0,m1/0,1/
c
      l0 = newbas0
      call putq(zcom,ztitle,e,e,l1,l1,l0,
     *m1,m0,q,ndaf,iblkf)
      call wrt3(e,l1,iblk,idaf)
      return
c
      end
      subroutine denshl(trans,d,l1,iham)
c
c     ----- calculate the density matrix for the iham'th
c           fock operator -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/scfwfn)
      dimension trans(l1,*),d(*)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     * tti(25),e1i(25),asym,
     * cilow(12),jpair,kone,ktwo,kcorb(2,12),
     * iojk(49),ioham(25),iojkao(49)
      data dzero/0.0d0/
c
      call vclr(d,1,ikyp(l1))
      do 120 kx = 1,norb
      kbf = msympr(kx)
      nx = nconf(kx)
      if (iham .gt. nham) go to 100
      if (nx .ne. iham) go to 120
      ij=1
      do 140 i=1,l1
      if(trans(i,kbf).ne.dzero) then
      call daxpy(i,trans(i,kbf),trans(1,kbf),1,d(ij),1)
          endif
140   ij=ij+i
      go to 120
  100 fact=f(nx)
      ij=1
      do 150 i=1,l1
      dum=fact*trans(i,kbf)
      if(dum.ne.dzero) then
      call daxpy(i,dum,trans(1,kbf),1,d(ij),1)
         endif
 150  ij=ij+i
  120 continue
      return
c
      end
      subroutine ocbse(iham,trans,dj,h,temp,pp,eig,dee,dep,
     1                 itervv,ncall,l1,l2)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/mapper)
INCLUDE(common/atmol3)
      dimension trans(l1,*),dj(l1,*),h(*),pp(*),eig(*),temp(*)
      common/diisd/st(210),ccc(20),rrr(20),scale(20),iposit(20),
     + nstore,mp,ondiis,nspaca,stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos
      dimension istrt(maxorb)
c
c ------- perform ocbse step
c
       ioc=0
      do 900 isy=1,nsytyp
      nop=melsym(isy)
      mo=0
c
c ------ if any occupied orbital belongs to the iham'th shell set istrt
c
      do 10 i=1,nop
      ioc=ioc+1
      if(nconf(ioc).ne.iham) goto 10
      mo=mo+1
      istrt(mo)=msympr(ioc)
 10   continue
c
c ------ move all virtuals to istrt
c
      if(mo.eq.0) goto 900
      nsize=ns(isy)-nop+mo
      jj=msympr(ioc)
      mop=mo+1
      do 20 i=mop,nsize
      jj=jj+1
 20   istrt(i)=jj
c
c ------ create the transformation matrix for this shell
c
      ij=1
      do 30 i=1,nsize
      call dcopy(l1,trans(1,istrt(i)),1,temp(ij),1)
 30   ij=ij+l1
c
c ------ transform the fock matrix to a shell fock matrix
c
      call mult2(temp,dj,h,nsize,nsize,l1)
      call dcopy(l2,dj,1,h,1)
      diaacc=diff*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc=5.0d-5
      if(diaacc.lt.1.0d-11)diaacc=1.0d-11
c
c ------ level shift mo fock matrix
c
      rshift=0.0d0
      call shiftq(h,mo,0,nsize,dee,dep,itervv,ncall,gapb1,ibrk,gapb2)
c
c ------ diagonalise this symmetry block
c
      m2=2
      call jacobi(h,iky,nsize,dj,ilifq,l1,pp,m2,lock,diaacc)
c
c ------ store relevant eigen values
c
      do 41 i=1,nsize
      sub=0.0d0
      if(i.gt.mo) sub=rshift
      eig(istrt(i))=pp(i)-sub
      pp(i)=0.0d0
 41   continue
c
c ------ ensure that phases agree
c
      do 50 i=1,nsize
      phase=1.0d0
      if(dj(i,i).lt.0.0d0) phase=-1.0d0
      call dscal(nsize,phase,dj(1,i),1)
 50   continue
c
c ------ apply transformation to relavent shells of vector
c
      do 90 i=1,l1
      do 70 k=1,nsize
      val=trans(i,istrt(k))
      if(val.eq.0.0d0) goto 70
c
c ------ accumulate changes to occupied orbitals
c
_IF1(civu)      do 60 j=1,nsize
_IF1(civu)   60 pp(j)=pp(j)+val*dj(k,j)
_IFN1(civu)      call daxpy(nsize,val,dj(k,1),l1,pp,1)
 70   continue
      do 80 k=1,nsize
      trans(i,istrt(k))=pp(k)
 80   pp(k)=0.0d0
 90   continue
 900  continue
c
      return
      end
      subroutine schort(u,nrow,ncol,ndim)
c
c     ----- schmidt orthonormalize the eigenvectors in the
c           transformed representation -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension u(ndim,*)
INCLUDE(common/iofile)
      data thresh/1.0d-10/
      do 160 j = 1,ncol
      do 120 k = 1,j
      sum = ddot(nrow,u(1,k),1,u(1,j),1)
      if (k .eq. j) go to 140
      summ = -sum
      call daxpy(nrow,summ,u(1,k),1,u(1,j),1)
 120  continue
  140 if (sum .lt. thresh) go to 180
      sum = 1.0d0/dsqrt(sum)
      call dscal(nrow,sum,u(1,j),1)
 160  continue
      return
  180 write (iwr,9028) j
      write (iwr,9008) sum
      call caserr('eigenvector with zero norm detected')
 9008 format(1x,'norm is ',e13.6)
 9028 format(1h0/5x,'eigenvector',i4,5x,'has dzero norm' /)
      return
      end
      subroutine optci(h,dj,b,ia,l1,l2,nprint)
c
c     ----- loop over all pairs and optimize the ci coefficients -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/scfwfn)
INCLUDE(common/scfopt)
      dimension dj(*),h(*)
      dimension hci(3),eig(2,2),eci(2),b(*)
      dimension ia(*)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      data dzero,pt01/0.0d0,0.01d0/
      if (npair .eq. 0) return
      test = acurcy*pt01
      nterm = 2
      iaa = l2+1
c
c     ib = l2+l2
c
c     ----- get j(i,j) and k(i,j) integrals in place -----
c
      call rdedx(h,l2,iojk(1),num8)
      kcold=0
      do 120 i = 1,norb
      kc=nconf(i)
      if (kc .eq. ncores) go to 120
      if (kc.eq.kcold) go to 150
      call reads(h,l2,num8)
      call reads(h(iaa),l2,num8)
  150 do 100 j = 1,norb
      if (nconf(j) .eq. ncores) go to 100
      jbf = msympr(j)
      jj = ia(jbf)+jbf
      kk = jj+l2
      ij = ia(max(i,j))+min(i,j)
      ijk = ij+l2
      dj(ij) = h(jj)
      dj(ijk) = h(kk)
  100 continue
  120 kcold=kc
      call rdedx(h,l2,iojk(1),num8)
c
c     ----- start ci iteration over pairs -----
c
      iterci = 0
  140 continue
      iterci = iterci+1
      conci = dzero
      do 160 jpair = 1,npair
      kone = kcorb(1,jpair)
      ktwo = kcorb(2,jpair)
c
c     ----- optimize the jpair'th pair -----
c
      call newci(h,dj,dj(iaa),conci,hci,eig,eci,b,ia,nterm,
     +           nprint)
  160 continue
      if (npair .eq. 1) conci = dzero
c
c     ----- test for convergence of the ci optimization -----
c
      if (conci .gt. test) go to 140
c
c     ----- normalize the ci coefficients and update the f, alpha,
c           and beta matrices -----
c
      call ciexpr(kcorb,nconf,nham,ia,l1)
      if (nprint .ne. 5) go to 200
      write (iwr,9008) iterci
      write (iwr,9028)
      do 180 j = 1,nham
      write (iwr,9048) j,f(j)
  180 continue
      call prtri(alpha,nham)
      call prtri(beta,nham)
  200 continue
c
      return
 9008 format(/5x,'ci part converged in ',i3,2x,'iterations' )
 9028 format(/5x,'energy coefs' /)
 9048 format(5x,'shell ',i2,2x,'hf= ',e15.8)
      end
      subroutine newci(h,dj,dk,conci,hci,eig,eci,b,ia,
     + nterm,nprint)
c
c     ----- calculate a new set of ci coefficients for the
c           jpair'th pair -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/scfwfn)
      dimension dj(*),dk(*),h(*)
      dimension b(*)
      dimension ia(*)
      dimension hci(*),eig(nterm,*),eci(*)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      data dzero,two,thresh/0.0d0,2.0d0,1.0d-08/
c     itest = 1
      ipick = 1
      ibf = msympr(kone)
      nn = (nterm*(nterm+1))/2
      jbf = msympr(ktwo)
      i1 = ia(ibf)+ibf
c     i2 = ia(max(ibf,jbf)) + min(ibf,jbf)
      i3 = ia(jbf)+jbf
      call vclr(hci,1,nn)
c
c     add hcore to ci matrix
c
      hci(1) = two*h(i1)
      hci(3) = two*h(i3)
      hci(2) = dzero
c
c      ----- process j and k matrices -----
c
      do 120 ii = 1,norb
      nx = nconf(ii)
      ix = ia(max(kone,ii))+min(kone,ii)
      iy = ia(max(ktwo,ii))+min(ktwo,ii)
      if (nx .eq. ncores) go to 120
      if (ii .eq. kone .or. ii .eq. ktwo) go to 120
c
c     ----- add potential terms from other shells -----
c
      cnst = f(nx)
      hci(1) = hci(1)+cnst*two*(two*dj(ix)-dk(ix))
      hci(3) = hci(3)+cnst*two*(two*dj(iy)-dk(iy))
  120 continue
      ix = ia(kone)+kone
      iy = ia(ktwo)+ktwo
      iz = ia(max(kone,ktwo))+min(kone,ktwo)
      hci(1) = hci(1)+dj(ix)
      hci(3) = hci(3)+dj(iy)
      hci(2) = dk(iz)
      if (nprint .ne. 5) go to 160
      ib = 0
      write (iwr,9008)
      do 140 j = 1,2
      iaa = ib+1
      ib = ib+j
  140 write (iwr,9028) (hci(i),i = iaa,ib)
  160 continue
      hstr = hci(1)
      if (hci(3) .lt. hstr) hstr = hci(3)
c
c     ----- get new ci coefficients -----
c
      call gldiag(nterm,nterm,nterm,hci,b,eci,eig,ia,2)
      signsw = eig(1,ipick)
      do 180 k = 1,nterm
      if (signsw .lt. dzero) eig(k,ipick) = -eig(k,ipick)
      conci = conci+(cicoef(k,jpair)-eig(k,ipick))**2
      cicoef(k,jpair) = eig(k,ipick)
  180 continue
      cilow(jpair) = eci(1)-hstr
      if (nprint .ne. 5) go to 220
      write (iwr,9048)
      do 200 j = 1,nterm
      write (iwr,9068) j,eci(j),(eig(i,j),i = 1,nterm)
  200 continue
  220 continue
      lone = nconf(kone)
      ltwo = nconf(ktwo)
      f(lone) = eig(1,ipick)**2
      f(ltwo) = eig(2,ipick)**2
      if (f(lone) .lt. thresh) go to 240
      if (f(ltwo) .lt. thresh) go to 240
      if ((eig(1,ipick)/eig(2,ipick)) .gt. dzero) go to 260
      return
  240 write (iwr,9088) f(nx),ii
      call caserr('invalid occupation number detected in gvb pair')
  260 continue
      write (iwr,9108) (eig(k,ipick),k = 1,nterm)
      call caserr('gvb pair orbitals are ill-defined')
 9008 format(/5x,'ci matrix '/)
 9028 format(/5x,5e15.8/(5x,5e15.8))
 9048 format(//5x,'ci solutions' /)
 9068 format(5x,i5,5x,'energy = ',f20.10/(5x,5e15.8))
 9088 format(/5x,'occ no = ',e15.8,5x,'orb ',i5)
 9108 format(//5x,'ci coef same sign -goodby',2e15.8//)
c
      return
      end
      subroutine gvbham(dj,dk,h,nwscmu,l2)
c
c     ----- calculate a fock matrix for each fock operator -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension dj(*),dk(*),h(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/scfwfn)
INCLUDE(common/scfopt)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      data done/1.0d0/
c
c ------ calculate as many fock operators as possible
c ------ in the given core.
c
      nfock=nwscmu/l2/3
      nfock=nfock+nfock
      nstart=1
 100  nend=nstart+nfock-1
      if(nend.gt.nham) nend=nham
c
c ------ deal with 1-e part
c
      call rdedx(dj,l2,iojkao(1),num8)
      ioff=1
      do 120 iham=nstart,nend
      call dcopy(l2,dj,1,h(ioff),1)
 120  ioff=ioff+l2
c
c ----- deal with j/k matrices nb working in ao basis
c
      do 150 jham=ncores+1,nham
      call reads(dj,l2,num8)
      call reads(dk,l2,num8)
      ioff=1
      do 150 iham=nstart,nend
      fact=done/f(iham)
      n=iky(max(iham,jham))+min(iham,jham)
      a=alpha(n)*fact
      call daxpy(l2,a,dj,1,h(ioff),1)
      b=beta(n)*fact
      call daxpy(l2,b,dk,1,h(ioff),1)
 150  ioff=ioff+l2
      ioff=1
      do 210 iham=nstart,nend
      call tranp(h(ioff),dj)
      call wrt3(dj,l2,ioham(iham),num8)
 210  ioff=ioff+l2
c
c ------ now transform hams
c
      nstart=nstart+nfock
      if(nstart.le.nham) goto 100
      return
      end
_EXTRACT(lagrng,ultra)
      subroutine lagrng(l1,noc,l2,h,dj,trans,ia)
c
c     ----- calculate the lagrangian multiplier matrix -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/scfwfn)
      common/diisd/ sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      dimension h(l2),dj(noc,*),trans(l1,*)
      dimension ia(*)
      data dzero,two/0.0d0,2.0d0/
c
c     ----- skip past the j and k matrices to the fock matrices -----
c
      call search(ioham(1),num8)
c
c     ----- loop over the fock operators -----
c
      do 200 iham = 1,nham
      call reads(h,l2,num8)
      do 180 i = 1,noc
      ii = msympr(i)
      if (nconf(i) .ne. iham) go to 180
      do 160 j = 1,noc
      jj = msympr(j)
      dj(jj,i) = dzero
      do 140 k = 1,l1
      do 140 l = 1,l1
      kl = ia(max(k,l))+min(k,l)
      dj(jj,i) = dj(jj,i)+trans(k,j)*h(kl)*trans(l,ii)
  140 continue
  160 continue
  180 continue
  200 continue
c
c     ----- convert to the generalized fock operator -----
c
      do 260 i = 1,noc
      sum=f(nconf(i))*two
      call dscal(noc,sum,dj(1,i),1)
 260  continue
c
c     ----- find the largest asymmetry -----
c
      asym = dzero
      do 280 i = 1,noc
      do 280 j = 1,i
      diff =  dabs(dj(i,j)-dj(j,i))
      if (diff .gt. asym) asym = diff
  280 continue
      return
      end
_ENDEXTRACT
      subroutine gvbgrd(hold,h,h2,ia,l2)
c
c     ----- rotate the occupied orbitals among themselves -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension hold(*),h(*),h2(*)
      dimension ia(*)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/scfwfn)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
c
      ek=ddot(nham,tti(1),1,f(1),1)
      ehf=ddot(nham,e1i(1),1,f(1),1)
      call vclr(hold,1,ikyp(nsos))
      call rdedx(h,l2,iojk(1),num8)
c
c ------ deal with closed shell fock matrix
c
      nhi=0
      do 80 is=1,nsytyp
      maxs=nhi+melsym(is)
      nlo=nhi+1
      nhi=nhi+ns(is)
      if(nlo.gt.nhi) goto 80
      ibf=msympr(nlo)-1
      do 20 i=nlo,nhi
      ibf=ibf+1
      kc=nconf(i)
      ii=ia(ibf)+ibf
      jbf=msympr(nlo)-1
      if(i.le.maxs)ehf=ehf+f(kc)*h(ii)
      do 20 j=nlo,nhi
      jbf=jbf+1
      lc=nconf(j)
      jj=ia(jbf)+jbf
      if(kc.eq.lc) goto 20
      if(kc.lt.lc) goto 15
      ij=ia(ibf)+jbf
c
c ------ 1st derivative
c
      hold(ij)=hold(ij)-f(kc)*h(ij)
      goto 20
 15   ij=ia(jbf)+ibf
c
c ---- 1st derivative
c
      hold(ij)=hold(ij)+f(kc)*h(ij)
 20   continue
 80   continue
      nhi = 0
      do 240 is = 1,nsytyp
      maxs=nhi+melsym(is)
      nlo = nhi+1
      nhi = nhi+ns(is)
      if (nlo .gt. nhi) go to 240
      kcold = 0
      ibf=msympr(nlo)-1
      do 220 i = nlo,nhi
      ibf=ibf+1
      kc = nconf(i)
      if (kc .le. ncores .or. kc .eq. (nham+1)) go to 220
      if (kc .eq. kcold) goto 220
      call reads(h,l2,num8)
      call reads(h2,l2,num8)
c
c ------ fock matrix contributions dealt with here
c ------ first loop over possible fock indices
c
      jbf=msympr(nlo)-1
      do 180 j=nlo,nhi
      jbf=jbf+1
      lc=nconf(j)
      jj=ia(jbf)+jbf
      kclc=ia(max(kc,lc))+min(kc,lc)
      aj=alpha(kclc)
      ak=beta(kclc)
      if(j.le.maxs) ehf=ehf+aj*h(jj)+ak*h2(jj)
c
c ------ loop over orbital index
c
      kbf=msympr(nlo)-1
      do 170 k=nlo,nhi
      kbf=kbf+1
      mc=nconf(k)
      if(mc.eq.lc) goto 170
      if(mc.gt.lc) goto 165
      jk=ia(jbf)+kbf
c
c ------ 1st deriv
c
      hold(jk)=hold(jk)-aj*h(jk)-ak*h2(jk)
      goto 170
 165   jk=ia(kbf)+jbf
c
c ------ 1st deriv
c
      hold(jk)=hold(jk)+aj*h(jk)+ak*h2(jk)
 170  continue
 180  continue
      kcold = kc
  220 continue
  240 continue
c
c ------ determine largest error
c
      loop=idamax(ikyp(nsos),hold,1)
      diff=dabs(hold(loop))
      return
      end
      subroutine rotb(hold,trans,h,h2,ia,l1,l2,gapa,ibrk,gapb,
     +                nprint)
c
c     ----- rotate the occupied orbitals among themselves -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/scfwfn)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
      dimension trans(l1,*),hold(l1,*),h(*),h2(*)
      dimension ia(*)
      common/diisd/st(210),ccc(20),rrr(20),scale(20),ipoosa(20),
     + nstore,mp,ondiis,nspaca, stb(270), iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      data dzero,pt5,two,ten/0.0d0,0.5d0,2.0d0,10.0d0/
c
      shift=gapa
      if(iter.gt.ibrk)shift=gapb
      rottol = diff*5.0d-3
      if(rottol .gt. 5.0d-5) rottol = 5.0d-5
      if(rottol .lt. acurcy*0.1d0) rottol = acurcy*0.1d0
      atest = rottol
      dtest = atest/ten
      low = 1
      if (irotb .lt. 0) low = 0
      rotmax = pt5
      do 120 j = 1,nsos
      call vclr(hold(1,j),1,nsos)
  120 continue
      call rdedx(h,l2,iojk(1),num8)
c
c ------ deal with closed shell fock matrix
c
      nhi=0
      do 80 is=1,nsytyp
      nlo=nhi+1
      nhi=nhi+melsym(is)
      if(nlo.gt.nhi) goto 80
      do 20 i=nlo,nhi
      kc=nconf(i)
      ibf=msympr(i)
      ii=ia(ibf)+ibf
      do 20 j=nlo,nhi
      lc=nconf(j)
      jbf=msympr(j)
      jj=ia(jbf)+jbf
      if(kc.eq.lc) goto 20
      if(kc.lt.lc) goto 15
      ij=ia(ibf)+jbf
c
c ------ 2nd derivative
c
      hold(jbf,ibf)=hold(jbf,ibf)+f(kc)*(h(jj)-h(ii))
c
c ------ 1st derivative
c
      hold(ibf,jbf)=hold(ibf,jbf)-f(kc)*h(ij)
      goto 20
c
c ---- 2nd derivative
c
 15   ij=ia(jbf)+ibf
      hold(ibf,jbf)=hold(ibf,jbf)+f(kc)*(h(jj)-h(ii))
c
c ---- 1st derivative
c
      hold(jbf,ibf)=hold(jbf,ibf)+f(kc)*h(ij)
 20   continue
 80   continue
      nhi = 0
      do 240 is = 1,nsytyp
      nlo = nhi+1
      nhi = nhi+melsym(is)
      if (nlo .gt. nhi) go to 240
      kcold = 0
      do 220 i = nlo,nhi
      kc = nconf(i)
      ibf = msympr(i)
      kckc=ia(kc)+kc
      ii=ia(ibf)+ibf
      if (kc .eq. ncores) go to 220
      if (kc .eq. kcold) goto 220
      call reads(h,l2,num8)
      call reads(h2,l2,num8)
      do 161 k=nlo,nhi
      mc=nconf(k)
      if(mc.ne.kc) goto 161
      kbf=msympr(k)
      do 160 j = nlo,k
      lc = nconf(j)
      kclc = ia(kc)+lc
      lclc = ia(lc)+lc
      jbf = msympr(j)
      if (kc .eq. lc) go to 160
      cj = beta(kckc)-beta(kclc)+beta(lclc)-beta(kclc)
      ck = beta(kckc)-beta(kclc)+two*alpha(kckc)-two*alpha(kclc)+
     +     beta(lclc)-beta(kclc)+two*alpha(lclc)-two*alpha(kclc)
      ij = ia(jbf)+jbf
      hold(jbf,kbf) = hold(jbf,kbf)+cj*h(ij)+ck*h2(ij)
  160 continue
 161  continue
c
c ------ fock matrix contributions dealt with here
c ------ first loop over possible fock indices
c
      do 180 j=nlo,nhi
      lc=nconf(j)
      lclc=ia(lc)+lc
      jbf=msympr(j)
      jj=ia(jbf)+jbf
      kclc=ia(max(kc,lc))+min(kc,lc)
      aj=alpha(kclc)
      ak=beta(kclc)
c
c ------ loop over orbital index
c
      do 170 k=nlo,nhi
      kbf=msympr(k)
      kk=ia(kbf)+kbf
      mc=nconf(k)
      if(mc.eq.lc) goto 170
      if(mc.gt.lc) goto 165
      jk=ia(j)+k
c
c ------ 2nd deriv
c
      hold(kbf,jbf)=hold(kbf,jbf)+aj*(h(kk)-h(jj))+ak*(h2(kk)-h2(jj))
c
c ------ 1st deriv
c
      hold(jbf,kbf)=hold(jbf,kbf)-aj*h(jk)-ak*h2(jk)
      goto 170
c
c ------ 2nd deriv
c
 165   jk=ia(k)+j
      hold(jbf,kbf)=hold(jbf,kbf)+aj*(h(kk)-h(jj))+ak*(h2(kk)-h2(jj))
c
c ------ 1st deriv
c
      hold(kbf,jbf)=hold(kbf,jbf)+aj*h(jk)+ak*h2(jk)
 170  continue
 180  continue
      kcold = kc
  220 continue
  240 continue
      if(nprint.ne.15.and.nprint.ne.5)go to 602
      write(iwr,603)
 603  format(//1x,39('-')/
     *1x,'elements of derivative matrix -- b(ij,ij)'/1x,39('-')//)
      call printd(hold,nsos,iwr)
c
c
c     ----- perform the actual rotation -----
c
602   ihi = 0
      nhi = 0
      thmax=0.0d0
      do 440 is = 1,nsytyp
      nlo = nhi+1
      nhi = nhi+melsym(is)
c     ncob = melsym(is)
c     jhi = ihi+ncob
c     ilo = ihi+1
      nsize = ns(is)
      ihi = ihi+nsize
      if (nhi .lt. nlo) go to 420
      do 400 j = nlo,nhi
      jbf = msympr(j)
      jc = nconf(j)
      do 380 i = nlo,j
      ibf = msympr(i)
      ic = nconf(i)
      if (jc .eq. ic) go to 380
      anum = hold(jbf,ibf)
      denom = hold(ibf,jbf) + shift
      if (low .ne. 0) denom =  dabs(denom)
      if ( dabs(denom) .lt. dtest) go to 360
      if ( dabs(anum) .lt. atest) go to 360
      theta = -anum/denom
      if( dabs(theta).gt.thmax)thmax= dabs(theta)
      hold(jbf,ibf) = theta
      hold(ibf,jbf) = -hold(jbf,ibf)
      go to 380
  360 continue
      hold(jbf,ibf) = dzero
      hold(ibf,jbf) = dzero
  380 continue
  400 continue
  420 continue
  440 continue
c
c ------ scale the transformation matrix
c
      scl=1.0d0
      if(thmax.gt.rotmax) scl=rotmax/thmax
      do 470 i=1,nsos
      call dscal(nsos,scl,hold(1,i),1)
 470  hold(i,i)=1.0d0
c
c ------ orthogonalise transformation matrix
c
      lt=l1*l1
      call dcopy(lt,hold,1,h2,1)
      call shalf(hold,nsos,nsos,l1,h,ifail)
      if(ifail.ne.0) then
      call dcopy(lt,h2,1,hold,1)
          endif
      call schort(hold,nsos,nsos,l1)
c
c ------ revise current set of vectors
c
      call tfsqc(hold,trans,h2,l1,l1,l1)
      call dcopy(lt,hold,1,trans,1)
      return
      end
      subroutine printd(der,nnn,iw)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension der(*),dd(6),ii(6),jj(6)
      data thr/1.0d-4/
c
      n=iabs(nnn)
      olow=.false.
      if(nnn.lt.0) olow=.true.
      k=0
      icount=0
      do 1 i=1,n
      max=n
      if(olow)max=i
      do 1 j=1,max
      icount=icount+1
      d=der(icount)
      if( dabs(d).lt.thr)go to 1
      k=k+1
      ii(k)=i
      jj(k)=j
      dd(k)=d
      if(k.lt.5)goto 1
      write(iw,2)(ii(k),jj(k),dd(k),k=1,5)
   2  format(5(2x,2i4,f10.4))
      k=0
   1  continue
      if(k.eq.0)return
      write(iw,2)(ii(j),jj(j),dd(j),j=1,k)
      return
      end
_IFN1(f)      subroutine jkform(trans,temp,dens,coul,exch,nwscmu,
_IF1(f)      subroutine jkform(trans,temp,dens,exch,coul,nwscmu,
     + nopk,npass,ltran,l1,l2,nprint,core)
c
c -- jkform drives the two electron integral part of the fock matrix
c    formation.  multiple passes over the integrals are done
c    only if needed. -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/scfwfn)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/zorac)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      dimension trans(l1,*)
      dimension dens(*),coul(*),exch(*),temp(*)
      dimension core(*)
c
      nscmu = nwscmu/l2
      kham = 1
      npass = 0
      iscmf = 1+ncores
      nleft = nham
  100 npass = npass + 1
      nclf = 0
      nscmf = nscmu/3
      if (npass .gt. 1) go to 120
      nclf = ncores
      nscmf = (nscmu-ncores*3)/3
      nleft = nleft - ncores
 120  if (nleft .gt. nscmf) go to 140
      nscmf = nleft
 140  nleft = nleft - nscmf
c
c     ----- symmetrize the fock matrix and transform to molecular orbita
c     operators -----
c
      call hstarg(trans,coul,exch,dens,nopk,
     +     nclf,nscmf,iscmf,l1,l2,nprint)
c
c     ----- closed shell part -----
c
      ioff = 1
      if (npass .gt. 1 ) go to 260
      call rdedx(dens,l2,ibl7f,num8)
c
c...  add zora corrections if required
c
      if (ozora) call zora(core,core,dens,'read')
c
      if (nco .le. 0) then
      call dcopy(l2,dens,1,temp,1)
      call wrt3(temp,l2,iojkao(kham),num8)
      call mult2(trans,dens,temp,ltran,ltran,l1)
      call wrt3(dens,l2,iojk(kham),num8)
      else
      if (nprint.eq.5) then
       write(iwr,*)' 2J-K closed shell'
       call prtri(coul,l1)
      endif
      call symh(coul,temp,iky,1,npair)
      call vadd(temp,1,dens,1,temp,1,l2)
      call wrt3(temp,l2,iojkao(kham),num8)
      call mult2(trans,dens,temp,ltran,ltran,l1)
      call wrt3(dens,l2,iojk(kham),num8)
      ioff=ioff+l2
      endif
 260  if (nscmf .ne. 0) then
c
c     ----- open shell part -----
c
      do 280 ifo = 1,nscmf
      if (nprint.eq.5) then
       write(iwr,*)' coulomb matrix, ifo = ',ifo
       call prtri(coul(ioff),l1)
      endif
      call symh(coul( ioff),temp,iky,1,npair)
      call wrt3(temp,l2,iojkao(kham+1),num8)
      call mult2(trans,dens,temp,ltran,ltran,l1)
      call wrt3(dens,l2,iojk(kham+1),num8)
      if (nprint.eq.5) then
       write(iwr,*)' exchange matrix, ifo = ',ifo
       call prtri(exch(ioff),l1)
      endif
      call symh(exch(ioff),temp,iky,1,npair)
      call wrt3(temp,l2,iojkao(kham+2),num8)
      call mult2(trans,dens,temp,ltran,ltran,l1)
      call wrt3(dens,l2,iojk(kham+2),num8)
      ioff = ioff + l2
  280 kham = kham+2
c
c     ----- get the original trans back -----
c
      endif
      if (nleft .lt. 1) go to 400
      iscmf = iscmf + nscmf
      go to 100
c
 400  return
      end
      subroutine virset(trans,h,d,l1,l2,numscr,core)
c
c     ----- virset sets up the data needed for calculating the
c                  virial ratio -----
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/zorac)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     *ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos
     *,tti(25),e1i(25),asym
     *,cilow(12),jpair,kone,ktwo,kcorb(2,12)
     *,iojk(49),ioham(25),iojkao(49)
      dimension ehold(25,2)
      equivalence (ehold(1,1),tti(1))
      dimension trans(l1,*),h(*),d(*)
      dimension core(*)
c
c     ----- read t matrix -----
c
      call rdedx(h,l2,ibl7t,numscr)
      do 100 i = 1,nham
      call denshl(trans,d,l1,i)
      ehold(i,1) = tracep(h,d,l1)
  100 continue
c
c     ----- read t + v matrix -----
c
      call rdedx(h,l2,ibl7f,numscr)
c
c...  add zora corrections if required
c
      if (ozora) call zora(core,core,h,'read')
c
      do 120 i = 1,nham
      call denshl(trans,d,l1,i)
      ehold(i,2) = tracep(h,d,l1)
  120 continue
      return
      end
      subroutine rhfvir(dens,ehf,enuc,ek,vir,numscr,l1,l2)
c
c  compute virial (vir = v/2e)
c
      implicit none
      REAL dens(*), ehf,enuc,ek,vir
      integer numscr,l1,l2

INCLUDE(common/vcore)
INCLUDE(common/scra7)

      external igmem_alloc
      integer igmem_alloc

      external tracep
      REAL tracep

      REAL etot

      integer itmp
      itmp = igmem_alloc(l2)
      call rdedx(Q(itmp),l2,ibl7t,numscr)
      ek = tracep(Q(itmp),dens,l1)
      call gmem_free(itmp)
      etot = enuc + ehf
      vir = (etot-ek)/(2.0d0*etot)
      end
c
      subroutine uhfvir(adens,bdens,ehf,enuc,ek,vir,numscr,l1,l2)
      implicit none
      REAL adens(*), bdens(*), ehf,enuc,ek,vir
      integer numscr,l1,l2

INCLUDE(common/vcore)
INCLUDE(common/scra7)

      external igmem_alloc
      integer igmem_alloc

      external tracep
      REAL tracep

      REAL etot

      integer itmp
      itmp = igmem_alloc(l2)
      call rdedx(Q(itmp),l2,ibl7t,numscr)
      ek = tracep(Q(itmp),adens,l1) + 
     &     tracep(Q(itmp),bdens,l1) 
      call gmem_free(itmp)
      etot = enuc + ehf
      vir = (etot-ek)/(2.0d0*etot)
      end
      subroutine shalf(c,nbas,norb,ndim,work,ifail)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension c(ndim,*),work(norb,norb)
      dimension delk(4),temp(maxorb)
      ifail=0
c
c --- symmetric orthogonalisation routine
c --- coded by j.kendrick july 1981 ici runcorn england
c
      itmax=20
      do 100 iter=1,itmax
      del=0.0d0
      do 40 i=1,norb
      do 40 j=1,norb
      sij=-ddot(nbas,c(1,i),1,c(1,j),1)
      if(sij.gt.1.0d+7) goto 140
      if(i.ne.j)del=del+sij*sij
 40   work(j,i)=sij
      if(iter.le.3) goto 48
c
c ------ check on convergence
c
      do 44 i=1,3
      if(del.lt.delk(i)) goto 48
   44 continue
      goto 140
  48  ipos=mod(iter-1,3)+1
      delk(ipos)=del
      do 50 k=1,norb
 50   work(k,k)=work(k,k)+3.0d0
      do 80 i=1,norb
      do 70 j=1,norb
      temp(j)=0.5d0*ddot(norb,c(i,1),ndim,work(1,j),1)
 70   continue
      do 75 k=1,norb
 75   c(i,k)=temp(k)
 80   continue
      if(del.lt.1.0d-12) return
 100  continue
 140  continue
      ifail=1
       write(6,120)delk
 120   format(1x,'***** symmetric orthogonalisation failure ****'/
     1 1x,4(1x,e10.4))
      return
      end
      subroutine canon(core,trans,dj,h,e,l0,l1,l2,l3,nprint)
c
c     ----- canonicalize gvb/grhf mos -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension core(*)
      dimension trans(*),dj(*),h(*),e(*)
INCLUDE(common/sizes)
INCLUDE(common/fsymas)
INCLUDE(common/runlab)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scfwfn)
INCLUDE(common/atmol3)
INCLUDE(common/harmon)
INCLUDE(common/dnfnw)
      common/junk/space(maxorb)
      common/diisd/ sta(270),iposa(24),stb(270),iposb(24),
     + ns(30),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      dimension ib(maxorb)
      data m1,m2/1,2/
c
      call rdedx(trans,l3,ibl3qa,idaf)
      call rdedx(h,l2,ioham(1),num8)
      do 6000 i=1,l0
 6000 ib(i)=(i-1)*l0
      nopen=0
      if(nseto.gt.0) then
       do 6002 i=1,nseto
6002   nopen=nopen+no(i)
      endif
      no1=nco+1
      no2=nco+nopen+npair+npair
      nvirt=l0-no2
      nv1=no2+1
      lprnt = newbas0
      if(.not.oprint(20)) lprnt = min(no2+5,newbas0)
c
c ----- first canon. the domos
c
      if(nco.le.1)go to 6003
      call mult1a(dj,trans,ib,nco,h,iky,l0)
      call jacobi(dj,iky,nco,trans,ib,l0,e,m1,m2,1.0d-8)
c
c ----- now canon. the pomos in shells
c
 6003 if(nopen.le.1)go to 6004
       do 6007 i=1,nseto
       nopen=no(i)
       if(nopen.le.1)go to 6007
      call mult1a(dj,trans,ib(no1),nopen,h,iky,l0)
      call jacobi(dj,iky,nopen,trans,ib(no1),l0,e(no1),
     * m1,m2,1.0d-8)
 6007 no1=no1+nopen
c
c ----- finally vmos
c
 6004 if(nvirt.le.1)go to 6005
      call mult1a(dj,trans,ib(nv1),nvirt,h,iky,l0)
c...   if harmonic reduction last vectors are 0.0 => 
c...   last elements of dj are 0.0 => shift them away
      if (newbas0.lt.l0) then
         do i=nvirt-l0+newbas0+1,nvirt
           if (dj(iky(i+1)).ne.0.0d0)call caserr('harmonic fails canon')
           dj(iky(i+1)) = 999999.99d99
         end do
      end if
c
      call jacobi(dj,iky,nvirt,trans,ib(nv1),l0,e(nv1),
     * m1,m2,1.0d-8)
c
c ----- back transform to ao basis
c
 6005 continue
      if(moutb.le.0)go to 6010
      if(nprint.ne.-5)write(iwr,6011)moutb
 6011 format(/
     *' output canonicalised m.o to section ',i3,
     *' of dumpfile')
c
c ----- setup occupation numbers -
c
      call vclr(space,1,l1)
      do 6009 i=1,norb
      no1=nconf(i)
 6009 space(i)=f(no1)+f(no1)
      call putq(zcom,ztitle,e,space,l1,l1,newbas0,
     *m1,m1,trans,moutb,i)
6010  if(nprint.eq.-5)return
      call analmo(trans,e,space,ilifq,newbas0,l1)
      call tdown(trans,ilifq,trans,ilifq,newbas0)
      if(otsym) call symass(trans,e,space,core)
      if(oprint(25))return
      write(iwr,6006)
 6006 format(/1x,104('-')//
     *40x,36('=')/
     *40x,'gvb natural orbitals (canonicalised)'/
     *40x,36('=')/)
      call prev(trans,e,lprnt,l1,l1)
c
c     evaluate the correlation energy
c     using density functional formula
c     only HF wavefunction
c
      if(npair.le.0) then
       if(odenfu) call becke (trans,2)
      endif
      return
      end
      subroutine scf(q)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/cslosc)
INCLUDE(common/field)
INCLUDE(common/scrf)
INCLUDE(common/gvalue)
INCLUDE(common/tran)
INCLUDE(common/runlab)
INCLUDE(common/scra7)
INCLUDE(common/symtry)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/funct)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/restrj)
INCLUDE(common/restrl)
INCLUDE(common/filel)
INCLUDE(common/phycon)
INCLUDE(common/zorac)
INCLUDE(common/statis)
INCLUDE(common/timeperiods)
INCLUDE(common/segm)
INCLUDE(common/dnfnw)
INCLUDE(common/xfield)
INCLUDE(common/datgue)
INCLUDE(common/psscrf)
      common/craypk/ijkl(1360)
      common/blkcore/corev(512),array(10)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/diisd/ sta(210),cca(20),ra(20),scalea(20),iposa(20),
     + nstora,mpa,odiisa,
     + nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     + mpb,odiisb,nspacb,nsss(30),igvbo(maxorb),igvb1(maxorb),igsp(4),
     + ekk(63),intci(150)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/scra/iso(mxshel,48)
c
_IF(drf)
INCLUDE(../drf/comdrf/sizesrf)
_ENDIF
_IF(parallel)
INCLUDE(common/nodeio)
INCLUDE(common/parcntl)
_ENDIF
      character *8 title,guess
      common/restrz/title(12),guess
      dimension zcas(3),q(*),o1e(6)
c
INCLUDE(common/drfopt)
_IF(drf)
cdrf
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfbem)
c
_ENDIF
      character*5 fnm
      character*3 snm
      data fnm,snm/"scf.m","scf"/
      data zrhf /'rhf'/,zcas   /'casscf','mcscf','vb'/
      data zuhf /'uhf'/
      data m1,m10,m16/1,10,16/
      data zgvb,zgrhf /'gvb','grhf'/
      data dzero /0.0d0/
c
      call cpuwal(begin,ebegin)
      call start_time_period(TP_SCF)
      if(nprint.ne.-5)write(iwr,180)
 180  format(/1x,104('-'))
c
c     ----- set convergence criteria -----
c
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
      nstora= 0
      nstorb= 0
      accdi2=acurcy*acurcy
c
c     clear buffer for integral output
c
      call setsto(1360,0,ijkl)
      odiisa=.false.
      odiisb=.false.
      if(zscftp.eq.zcas(1).or.zscftp.eq.zcas(2)
     &   .or. zscftp.eq.zcas(3)  ) then
        go to 130
      endif
c
c     ----- read in transformation matrices for s,p,d,f and g basis 
c           functions.
c
      nav = lenwrd()
      if (nt.gt.1) then
        call rdedx(ptr,nw196(1),ibl196(1),idaf)
        if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
        if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
        if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
        call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      endif
      l2=num*(num+1)/2
      len2=l2+l2
      ilen = 3*l2
      if (ofield) then
        ilen = ilen + 3*l2
      endif
      i10 = igmem_alloc_inf(ilen,fnm,snm,"i10",IGMEM_DEBUG)      
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
      if (ofield) then
       i50  = i40 + l2
       i60  = i50 + l2
       last = i60 + l2
      else
       last=i40
      endif
      length=last-i10

      if(length .ne. ilen)then
         write(6,*)length,ilen
         call caserr('mem size error')
      endif
c
c ----- if scrf then restore dipole integrals
c       and write to ed7                       ----
c
      if(oscrf) call dipmat(zscftp,q)
c
c ----- first restore s,t,f from section 192 and store on ed7
c
c ----- are field values to be introduced
c
      if(ofield) then
c
       do loop =1,6
        o1e(loop) = .true.
       enddo
       call getmat(q(i10),q(i20 ),q(i30),q(i40),q(i50),q(i60),
     * array,num,o1e,ionsec)
c ----- field calculation
       do 220 loop=1,l2
       q(i30+loop-1) = q(i30+loop-1) - fieldx * q(i40+loop-1)
     *                               - fieldy * q(i50+loop-1)
     *                               - fieldz * q(i60+loop-1)
 220   continue
       if(nprint.ne.-5) then
        write(iwr,230)fieldx,fieldy,fieldz
 230    format(/5x,38('*')/
     *          5x,'one-electron hamiltonian modified with'/
     *          5x,'x *',f11.7/
     *          5x,'y *',f11.7/
     *          5x,'z *',f11.7/
     *          5x,38('*') /)
       endif
      else
c ------ scrf calculation
       if(oscrf) then
c ------ convert angstrom to a.u.in g aradix,y,z
        aradix=aradix/toang(1)
        aradiy=aradiy/toang(1)
        aradiz=aradiz/toang(1)
c ------ calculate g factor (if sphere then only need use x)
        gx = (2.0d0*(dielec-1.0d0))/
     +                ((2.0d0*dielec + 1.0d0)*(aradix**3))
c ----- check if charged, and alter g accordingly
c ----- +ve ions
        if(ich.gt.0)then
                gx = gx*ne
                gx = gx/(ne+ich)
        endif
c ----- -ve ions
        if(ich.lt.0)then
                gx = gx*(ne+ich)
                gx = gx/ne
        endif
c
        gy = gx
        gz = gx
c
        if(nprint.ne.-5) then
        write(iwr,*)'    **********************************************'
        write(iwr,231)gx
        write(iwr,232)aradix
        write(iwr,*)'    **********************************************'
 231    format(/5x,'*** f factor set at            ',f11.7,' ***')
 232    format(/5x,'*** solvent cavity set at      ',f7.3,'au   ***'/)
        endif
       endif
c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S) then
       oswed3(4) = .false.
       oswed3(8) = .false.
      endif

      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
       if(omaster) then
_ENDIF
        do loop =1,3
         o1e(loop) = .true.
         o1e(loop+3) = .false.
        enddo
        call getmat(q(i10),q(i20 ),q(i30),q(i10),q(i10),q(i10),
     *  array,num,o1e,ionsec)
_IF(parallel)
       endif
_ENDIF
      endif
_IF(drf)
cdrf
c ------ drf addition to one-el ham.--------------------drf extension
c
      if (oreact) then
c
      if (field(:4) .ne. 'scf' .and. field(5:) .ne. 'scf') then
        if (idrfout .ge. 1) write (iwr,*) ' no scf reactionfield'
      else
        if (idrfout .ge. 1)
     +    write (iwr,*) '     **** reaction field additions ****'
c
c * * for mcscf or ci wfn's 2-el. drf-contributions should be included in
c     the 2e-integrals, for scf: wfnflg=0
c
c  -----  read internal one-electron hamiltonian
c
         if (idrfout .eq. 4)
     +   call hatout(q(i30),num,num,3,'h1-int')
c
c  -----  add potential external charges
c
	 i80 = igmem_alloc_inf(nchd*ngran,fnm,snm,"i80",IGMEM_DEBUG)
	 i90 = igmem_alloc_inf(nchd*ngran,fnm,snm,"i90",IGMEM_DEBUG)
         if (field(:4) .eq. 'scf') then
           call daread(idafdrf,iodadrf,q(i80),nchd*ngran,12)
          if (idrfout .ge. 1)
     +     write (iwr,229) 'potential external charges'
           do 25, i = 1, ngran
             call addup(q(i30),q(i80+(i-1)*nchd),q(i30),nchd)
   25      continue
c
           if (idrfout .eq. 4)
     +      call hatout(q(i30),num,num,3,'h1-vext')
c
cafc  add neqsta potential to be implemented (see
c     hfmcci on wfnnew.f)
c
           if (neqrf .eq. 1) then
c      -----  add potential non-equilibrium rf
c
             call daread(idafdrf,iodadrf,q(i80),nchd,14)
             if (idrfout .ge. 1)
     +       write (iwr,229) 'potential non-equilibrium rf'
             do 29, i = 1, nchd
              q(i30+i-1) = q(i30+i-1) + q(i80+i-1)
   29        continue
c
             if (idrfout .eq. 4)
     +         call hatout(q(i30),num,num,3,'h1-vneq')
           endif
         endif
c
         if ((field(5:) .eq. 'scf') .and. (iarfcal .eq. 0)) then
c    -----  add potential electronic induced dipoles at external charges
c
           call daread(idafdrf,iodadrf,q(i80),nchd*ngran,52)
           if (idrfout .ge. 1)
     1     write (iwr,229) 'pot. electr. ind. dip.s at ext. charg'
           do 30, i = 1, ngran
             call addup(q(i30),q(i80+(i-1)*nchd),q(i30),nchd)
   30      continue
c
           if (idrfout .eq. 4)
     1       call hatout(q(i30),num,num,3,'hexdip')
c
          if (iextdip .eq. 0) then
c     3-----
c      -----  add potential external induced dipoles at expansion centra
c
            call daread(idafdrf,iodadrf,q(i80),nchd*ngran,17)
            do 31, i = 1, ngran
              call addup(q(i30),q(i80+(i-1)*nchd),q(i30),nchd)
  31        continue
c
            if (idrfout .eq. 4)
     1      call hatout(q(i30),num,num,3,'hexdip2')
c     3-----
          endif
c
c    -----  add one-electron self energy and "reverse" dispersion
c
           if (gamdrf .ne. dzero) then
             call daread(idafdrf,iodadrf,q(i80),nchd,58)
             if (idrfout .ge. 1)
     1       write (iwr,229) 'one-el  self-energy'
             do 33, i = 1, nchd
              q(i30+i-1) = q(i30+i-1) + q(i80+i-1)
   33        continue
c
             if (irevdis .eq. 1) then
               call daread(idafdrf,iodadrf,q(i80),nchd,13)
             if (idrfout .ge. 1)
     1       write (iwr,229) 'reverse disp'
               do 36, i = 1, nchd
                q(i30+i-1) = q(i30+i-1) + q(i80+i-1)
   36          continue
             endif
           endif
c
           if (idrfout .eq. 4)
     1     call hatout(q(i30),num,num,3,'hddis1')
c
c    -----  add screening of nuclear attraction
c
           call daread(idafdrf,iodadrf,q(i80),nchd,54)
           call daread(idafdrf,iodadrf,q(i90),nchd,15)
           if (idrfout .ge. 1)
     1     write (iwr,229) 'screening of nuc. attr'
           do 37, i = 1, nchd
            q(i30+i-1) = q(i30+i-1) + scffact*
     1                     (q(i80+i-1) + q(i90+i-1))
   37      continue
         endif
c
cafc     coupling ground-state potential to excited state
c        to be implemented (see hfmcci on wfnnew.f)
c
         call dawrit(idafh,ioda,q(i30),nchd,11,navh)
         if (idrfout .eq. 4)
     1    call hatout(q(i30),num,num,3,'hdrf-1')
  102  continue
  229      format(/5x,38('*')/
     +             5x,'one-electron hamiltonian modified with'/
     +             5x,'reactionfield contributions of type'/
     +             5x,a/
     +             5x,38('*') /)

	call gmem_free_inf(i90,fnm,snm,"i90")
	call gmem_free_inf(i80,fnm,snm,"i80")
      endif
c
      endif
c
cdrf-hmat--------------------------------------------------end drf
c
_ENDIF
_IF(parallel)
      if(omaster) then
_ENDIF
       call wrt3(q(i10),l2,ibl7s,num8)
       call wrt3(q(i20),l2,ibl7t,num8)
       call wrt3(q(i30),l2,ibl7f,num8)
_IF(drf)
      if(zscftp.eq.zcas(1).or.zscftp.eq.zcas(2)  ) then
         goto 130
      else
       if (oreact) call wrt3(q(i30),l2,ibl3f,idaf)
      endif
_ENDIF
c
c ----- transform the s-matrix
c
       call tranp(q(i10),q(i20))
       call wrt3(q(i20),l2,ibl7st,num8)
_IF(parallel)
      endif
      oswed3(4) = .true.
      oswed3(8) = .true.
_ENDIF
c
c ----- reset core
c
      call gmem_free_inf(i10,fnm,snm,"i10")
c
c     ----- establish buffers for the two electron integral file(s)
c
 130  continue
c     opank = .false.
c     if (nopk .eq. (-1)) opank = .true.
c     if ((zscftp .ne.zrhf .or. na .ne. nb) .and. nopk .ne. 1) opank =
c    +     .true.
      ogrhf = .false.
c
c     ----- execute scf procedure -----
c
      if (zscftp .eq. zgrhf) ogrhf = .true.
      sz = dzero
      s2 = dzero
c
      if(.not.odscf) then
       lfile=m2file
       do 141 i=1,lfile
       lotape(i)=m2tape(i)
       liblk(i)  =m2blk(i)
 141   llblk(i)  =liblk(i)-m2last(i)
      endif
c
c        atomic density startup for open and closed shells
c        modified for correct restart
c

      if (guess.eq.'atoms') then
c
c         do one cycle scf with the density matrix from denat
c
_IF(ga)
         if(zscftp.eq.zrhf.and. odscf.and.
     +         num.ge.idporth .and.
     +         num.ge.idpmult2.and.
     +         num.ge.idpdiis .and.
     +         num.ge.idpdiag .and.
     +         ipiomode .ne. IO_NZ_S)then
             call denscf_ga(q,zscftp)
         else
             call denscf(q,zscftp)
         endif
_ELSE
         call denscf(q,zscftp)
_ENDIF
         if (irest.ne.0) go to 160
c
c        set zguess to 'anything' so denscf will be called but once
c        but obviously not if atdens will be redone
c
         if (.not.oalway) then
            zguess = 'anything'
            guess  = 'anything'
         end if
      end if
c
      if(zscftp . eq. zcas(1) ) go to 150
      if(zscftp . eq. zcas(2) ) go to 250
      if(zscftp . eq. zcas(3) ) go to 251
      if (zscftp .eq. zgvb ) go to 140
      if (zscftp .eq. zgrhf) go to 145
      if (zscftp .eq. zuhf )go to 120
      if (zscftp .ne. zuhf .and. na .ne. nb) go to 100
c
c     ----- closed shell restricted hartree fock -----
c

      if (zscftp.eq.zrhf) then
        if (odscf) then
              if (oscrf) then
                 call dscrf(q)
              else
                 if(odnew) then
                  call drcdrv(q,'scf')
                  goto 160
                 else
_IF(ga)
                 if(num.ge.idporth .and.
     +              num.ge.idpmult2.and.
     +              num.ge.idpdiis .and.
     +              num.ge.idpdiag .and.
     +              ipiomode .ne. IO_NZ_S.and..not.oreact)then
                  call drhfcl_ga(q,c,czan)
                 else
                  call drhfcl(q,c,czan)
                 endif
_ELSE
                 call drhfcl(q,c,czan)
_ENDIF
                 endif
              endif
        else
c
c  first determine total amount of free memory available
         lwor = igmem_max_memory()
c  now determine core required given that lwor is available
c
cdft
           lscdf = 0
           if(osc)lscdf = l2
c  memory for grid - @ make this switchable
           if (osc.or.odenfu) then
            lscdf = lscdf + 50*194*nat
            write(6,*)'lscdf',lscdf
           endif

            nex = 0
            if(ocryst)nex = 1
            if(lwor.gt.((35+nex)*l2+num+lscdf)) then
               len=(31+nex)*l2+num+lscdf
               omem=.true.
            else
               len=len2+2*num+lscdf + nex*l2
               omem=.false.
            endif

_IF1()c         if(lwor.gt.(35*l2+num)) then
_IF1()c         len=31*l2+num
_IF1()c         omem=.true.
_IF1()c         else
_IF1()c         len=len2+2*num
_IF1()c         omem=.false.
_IF1()c         endif
_IF(unicos)
         nss=max(2,min(10,(lwor-len)/len2))
_ENDIF
_IF(convex)
         nss=max(2,min(18,(lwor-len)/len2))
_ENDIF
_IF(titan)
         nss=max(2,min(5,(lwor-len)/len2))
_ENDIF
_IFN(unicos,convex,titan)
      nss=2
_ENDIF
         lword=nss*len2+len
_IF(drf)
cdrf-----------------------create memory for drf
      if ((field(5:) .eq. 'scf') .and.
     1     (intdrf .eq. 0)) then
        if (iarfcal .eq. 0) then
        lword = lword + 5 + 9*l2 + nomga
        if (itwoeps .eq. 1) 
     1    lword = lword + nomga
        else
          lword = lword + 5 + 7*l2 + ndim*ndim
     1  + nwtr*nwtc + 4*ndim
        endif
      endif
cdrf--------------------------------------------
_ENDIF
c 
c sf ....   spin orbit zora calculation
c
         if (oso .and. ozora)  then
_IF(zora)
           print *,'Invoking Spin Orbit Zora Driver'
c ....  'closed' shells (aa density = bb density) and open shells
            call scf_so(q)
            return
_ELSE
            call caserr('Zora not available')
_ENDIF
         end if 
c
c        now allocate required core 
c 
c ------ force scrf calc through 'd' routine not 'm'
c
         odvsg=.true.
         if(oscrf) odvsg=.false.
c
c ------ carry on original code
c
         ibase = igmem_alloc_inf(lword,fnm,snm,'ibase',IGMEM_DEBUG)
         if(omem.and.odvsg.and..not.odisc) then
            call rhfclm(q,q(ibase),c,czan,lword,nss)
         else
            call rhfcld(q,q(ibase),c,czan,lword,nss)
         endif
         call gmem_free_inf(ibase,fnm,snm,'ibase')
        endif
      endif
      go to 160
c
c     ----- open shell restricted hartree fock -----
c           e.r.davidson rohf not available in this
c           version.  ocbse type rohf in gvb.
c
  100 call caserr('davidson rohf not available')
c
c     ----- open shell unrestricted hartree fock -----
c
  120 if (odscf) then
         if(odnew) then
             call caserr('direct-UHF not available')
         else
             call duhfop(sz,s2,q)
         endif
      else
         if (oso .and. ozora)  then
_IF(zora)
            print *,'Invoking Spin Orbit Zora Driver'
c ....  'closed' shells (aa density = bb density) and open shells
            call scf_so(q)
            return
_ELSE
             call caserr('Zora not available')
_ENDIF
         end if

         call uhfop(sz,s2,q)
      endif
      go to 160
c
c     ----- generalized valence bond calculations -----
c
  140 if (odscf) then
         if(odnew) then
             call caserr('direct-GVB not available')
         else
             call drhfgvb(q,ogrhf)
         endif
      else
         call rhfgvb(q,ogrhf)
      endif
      go to 160
c
c     ----- generalized hartree fock
c
  145 if (odscf) then
         if(odnew) then
             call caserr('direct-GVB not available')
         else
             call drhfgvb(q,ogrhf)
         endif
      else
         call rhfgvb(q,ogrhf)
      endif
      go to 160
c
c     ----- complete active space scf calculation
c
  150 call casscf(q)
      call timana(9)
      go to 170
c
c     ----- 2nd order mcscf calculation
c
  250 call multi(q)
      call timana(10)
      go to 170
c
c     ----- vbscf calculation
c
  251 if (opssc) then
         call mod1ed(c,czan,1)
      endif
      call vbstart(q,2)
      call timana(11)
      go to 170
c
  160 call timana(5)
  170 enrgy = etot
c
c     ----- save data
c
      array(1) = en
      array(2) = ehf
      array(3) = etot
      array(4) = sz
      array(5) = s2
      array(6) = ek
      array(7) = vir
      do loop  = 8, 10
       array(loop) = 0.0d0
      enddo
c
      call secput(isect(494),m16,m1,iblk16)
      call wrt3(array,m10,iblk16,idaf)
c
      call end_time_period(TP_SCF)
c
      return
      end
      subroutine scf2(q)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/cslosc)
INCLUDE(common/tran)
INCLUDE(common/runlab)
INCLUDE(common/scra7)
INCLUDE(common/symtry)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/funct)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/filel)
INCLUDE(common/statis)
INCLUDE(common/segm)
INCLUDE(common/timez)
      common/blkcore/corev(512),array(10)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/sta(210),cca(20),ra(20),scalea(20),iposa(20),nstora,
     + mpa,odiisa,
     + nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     + mpb,odiisb,nspacb,nsss(30),igvbo(maxorb),igvb1(maxorb),igsp(4),
     + ekk(63),intci(150)
      common/scra/iso(mxshel,48)
      common/craypk/ijkl(1360)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
c
c  ----- common block for psscrf calc.
c
INCLUDE(common/psscrf)
INCLUDE(common/xfield)
c
      character*5 fnm
      character*4 snm
      data fnm,snm/'scf.m','scf2'/
c
      character*10 charwall
      dimension zcas(3),q(*),o1e(6)
      dimension dump(10)
      data zrhf /'rhf'/,zcas   /'casscf','mcscf','vb'/
      data zuhf /'uhf'/
      data m1,m10,m16/1,10,16/
      data zgvb,zgrhf /'gvb','grhf'/
      data dzero /0.0d0/
c
      ktemp = ibl7la
      mwor16 = 3 + 6/lenwrd()
c
      if(irest.eq.0) then
        call scf(q)
      else
       call secget(msecp,m16,iblkin)
       call rdedx(ptspac,mwor16,iblkin,idaf)
       if(iscf.ne.0) go to 500
        call scf(q)
      endif
      if(irest.ne.0) then
       call mod1ed(c,czan,1)
       write(iwr,9367)
       go to 400
      endif
c
500   ibl7la= ktemp
c
c     ----- find the time remaining at the start -----
c
      call timrem(tlefts)
      call cpuwal(begin,ebegin)
      write(iwr,181)
 181  format(//20x,52('*')/20x,12('*'),' Pisa solvation calculation ',
     * 12('*')/20x
     *,3('*'),'J Tomasi et al, Chem. Phys., 55(1981), 117-129',3('*')
     */20x,6('*'),'GAMESS version by Darren Green Feb. 1991 ',6('*'),
     * /20x,52('*')/)
c
      skale=2.0d0
      thresh = 10.0d0**(-iconv)
       max10 = maxat*20
c
c      determine resident allocation of core
c
       i100 = 0
       i100z=i100 + max10*3
       i100d=i100z+ max10
       i100c=i100d+ max10*2
       i100e=i100c+ max10*3
       i100f=i100e+ max10*2
       last =i100f+ max10
c      now allocate required core
       i100 = igmem_alloc_inf(last,fnm,snm,'i100',IGMEM_DEBUG)
c
c      now allocate pointers
c
       i100z=i100 + max10*3
       i100d=i100z+ max10
       i100c=i100d+ max10*2
       i100e=i100c+ max10*3
       i100f=i100e+ max10*2
c
       if(irest.eq.0) iscf = iscf + 1
       call surgen(q,q(i100),q(i100z),q(i100d),q(i100c),
     +             q(i100f),iscf,iwr)
c
c     ----- set convergence criteria -----
c
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
      nstora= 0
      nstorb= 0
      accdi2=acurcy*acurcy
c
c ----- loop for psscrf
c       1) get surface charges from gas vectors
c       2) do scf calc.
c       3) get revised surface charges
c       4) do scf
c       5) if energy not stationary go to 3)
c
           call psanin
           call potini(q,q(i100),q(i100z),q(i100d),q(i100c),q(i100e),
     +      q(i100f))
           call clrefg(q(i100e))
           ibl7la= ktemp
           call mod1e(q,q(i100),q(i100z))
           call mod1ed(q(i100),q(i100z),1)
c
      etold=0.0d0
c
      do 1 loop=1,maxcyc
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
      nstora= 0
      nstorb= 0
      accdi2=acurcy*acurcy
c
c     clear buffer for integral output
c
           call setsto(1360,0,ijkl)
           odiisa=.false.
           odiisb=.false.
           if(zscftp.eq.zcas(1).or.zscftp.eq.zcas(2)
     &        .or.zscftp.eq.zcas(3) )go to 130
c
c     ----- read in transformation matrices for s,p,d,f and g basis functions.
c
           if (nt.gt.1) then
           nav = lenwrd()
             call rdedx(ptr,nw196(1),ibl196(1),idaf)
             if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
             if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
             if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
             call readi(iso,nw196(5)*nav,ibl196(5),idaf)
           endif
c     allocate memory for transforming 1e-integrals
           l2=num*(num+1)/2
           len2=l2+l2
           i10 = igmem_alloc_inf(len2+l2,fnm,snm,'i100',IGMEM_DEBUG)
           i20=i10+l2
           i30=i20+l2
           last=i30+l2
c
           do i = 1,3
            o1e(i) = .true.
            o1e(i+3) = .false.
           enddo
           call getmat(q(i10),q(i20 ),q(i30),q(i10),q(i10),q(i10),
     *      array,num,o1e,ionsec)
           call wrt3(q(i10),l2,ibl7s,num8)
           call wrt3(q(i20),l2,ibl7t,num8)
           call wrt3(q(i30),l2,ibl7f,num8)
c
c
c ----- transform the s-matrix
c
           call tranp(q(i10),q(i20))
           call wrt3(q(i20),l2,ibl7st,num8)
c
c ----- reset core
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
c     ----- establish buffers for the two electron integral file(s)
c
 130  continue
c     opank = .false.
c     if (nopk .eq. (-1)) opank = .true.
c     if ((zscftp .ne.zrhf .or. na .ne. nb) .and. nopk .ne. 1) opank =
c    +     .true.
      ogrhf = .false.
c
c     ----- execute scf procedure -----
c
      if (zscftp .eq. zgrhf) ogrhf = .true.
      sz = dzero
      s2 = dzero
c
      if(zscftp . eq. zcas(1) ) go to 150
      if(zscftp . eq. zcas(2) ) go to 250
      if(zscftp . eq. zcas(3) ) go to 251
      if(.not.odscf) then
       lfile=m2file
       do 141 i=1,lfile
       lotape(i)=m2tape(i)
       liblk(i)  =m2blk(i)
 141   llblk(i)  =liblk(i)-m2last(i)
      end if
c
      if (zscftp .eq. zgvb ) go to 140
      if (zscftp .eq. zgrhf) go to 145
      if (zscftp .eq. zuhf )go to 120
      if (zscftp .ne. zuhf .and. na .ne. nb) go to 100
c
c     ----- closed shell restricted hartree fock -----
c
      if (zscftp.eq.zrhf) then
        if (odscf) then
         if(odnew) then
          call drcdrv(q,'scf')
          goto 160
         else
          call drhfcl(q,q(i100),q(i100z))
         endif
        else
c
c  first determine total amount of free memory available, lwor
         lwor  = igmem_max_memory()
c  now determine core required given that lwor is available

         nex = 0
         if(ocryst)nex = 1
         if(lwor.gt.((35+nex)*l2+num)) then
            len=(31+nex)*l2+num
            omem=.true.
         else
            len=len2+2*num
            omem=.false.
         endif

_IF(unicos)
         nss=max(2,min(10,(lwor-len)/len2))
_ENDIF
_IF(convex)
         nss=max(2,min(18,(lwor-len)/len2))
_ENDIF
_IF(titan)
         nss=max(2,min(5,(lwor-len)/len2))
_ENDIF
_IFN(unicos,convex,titan)
      nss=2
_ENDIF
         lword=nss*len2+len
c now allocate memory
         ibase = igmem_alloc_inf(lword,fnm,snm,'ibase',IGMEM_DEBUG)
         if(omem) then
           call rhfclm(q,q(ibase),q(i100),q(i100z),lword,nss)
         else
           call rhfcld(q,q(ibase),q(i100),q(i100z),lword,nss)
         endif
         call gmem_free_inf(ibase,fnm,snm,'ibase')
        endif
      endif
c
c     ----- open shell restricted hartree fock -----
c           e.r.davidson rohf not available in this
c           version.  ocbse type rohf in gvb.
c
      go to 160
  100 call caserr('davidson rohf not available')
c
c     ----- open shell unrestricted hartree fock -----
c
  120 call uhfop(sz,s2,q)
      go to 160
c
c     ----- generalized valence bond calculations -----
c
  140 call rhfgvb(q,ogrhf)
      go to 160
c
c     ----- generalized hartree fock
c
  145 call rhfgvb(q,ogrhf)
      go to 160
c
c     ----- complete active space scf calculation
c
  150 call casscf(q)
      call timana(9)
      go to 170
  250 call multi(q)
      call timana(10)
      go to 170
  251 call mod1ed(q(i100),q(i100z),1)
      call vbstart(q,2)
      call secget(isect(494),16,iblkrem)
      call rdedx(dump,10,iblkrem,idaf)
      etot=dump(3)
      call mod1ed(q(i100),q(i100z),0)
      call timana(11)
      go to 170
  160 call timana(5)
  170 enrgy = etot
      if(irest.ne.0) then
      irest=3
      write(iwr,9367)
      call texit(0,irest)
      go to 400
      endif
c
c     ----- now check for convergence
c           only energy at moment dvsg -----
           contst=(etold-etot)
           if (contst.le.thresh) goto 2
c
c     ----- exit in case of time limit -----
c     ----- punch out restart data -----
c
      call timrem(tlefti)
      if (tlefti .gt. skale*(tlefts-tlefti)/loop) go to 380
      irest=3
      write(iwr,9367)
 9367 format(/10x,31('*')/
     *10x,'*** warning ***'/
     *10x,'solvation scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,31('*')/)
      call texit(0,irest)
      go to 400
c
c     ----- save data
c
 380       etold=etot
           array(1) = en
           array(2) = ehf
           array(3) = etot
           array(4) = sz
           array(5) = s2
           array(6) = ek
           array(7) = vir
           do loop2  = 8, 10
             array(loop2) = 0.0d0
           enddo
c
           call psanin
           call potini(q,q(i100),q(i100z),q(i100d),q(i100c),q(i100e),
     +                 q(i100f))
           call clrefg(q(i100e))
           ibl7la= ktemp
           call mod1e(q,q(i100),q(i100z))
           call mod1ed(q(i100),q(i100z),1)
c
1     continue
400   write(iwr,4)etot
4     format(1x,'***** psscrf not converged   *****'/1x,
     * '*****  energy= ',f13.6,' *****'/)
      goto 5
c
2     array(1) = en
      array(2) = ehf
      array(3) = etot
      array(4) = sz
      array(5) = s2
      array(6) = ek
      array(7) = vir
      do loop2  = 8, 10
       array(loop2) = 0.0d0
      enddo
c
      write(iwr,3)loop,etot
3     format(/1x,42('*')/
     +1x,'***** psscrf calculation converged   *****'/
     *1x,'*****       after',i4,'  cycles        *****'/
     *1x,'***** to an energy of ',f13.6,' H *****'/1x,42('*')//)
c
      irest = 0
      call tofren(etot,etot1,q(i100z),q(i100f))
c
      write(iwr,31)etot1
31    format(/1x,42('*')/
     + 1x,'*****  corrected for free energy     *****'/
     + 1x,'***** psscrf energy = ',f13.6,' H   *****'/1x,42('*')//)
      call secput(isect(494),m16,m1,iblk16)
      call wrt3(array,m10,iblk16,idaf)
c
5     continue
c
      call gmem_free_inf(i100,fnm,snm,'i100')
c
      cpu=cpulft(1)
      write(iwr,7777)cpu,charwall()
7777  format(/
     *' end of solvation scf at ',f12.2,' seconds',a10,' wall',/)
      return
      end
      subroutine extrpd(dee,dampp,damp00,h0,h1,h2,h3,l1,l2,ndaf,
     +  numscr,iterl,ncall,ityp)
c
c     -----  use either davidson's damping and extrapolation
c            or pople's extrapolation -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/runlab)
INCLUDE(common/extsav)
      dimension h0(*),h1(*),h2(*),h3(*)
      data dzero,done,two,four/0.0d0,1.0d0,2.0d0,4.0d0/
      data pt5/0.5d0/
      data tol,tol1,tol2/1.0d-07,1.9d0,0.99d0/
      data shrnkf,dmpmin,rsmin/100.0d0,0.01d0,0.8d0/
      data zgvb,zgrhf /'gvb','grhf'/
c
c     ----- decode convergence method -----
c
      oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (ncall .eq. 1) kcount = kcount+1
      if (iter .gt. 1) go to 140
c
c     ----- iter = 1 -----
c
      call dcopy(l2,h0,1,h1,1)
      call dcopy(l2,h0,1,h2,1)
      go to 660
c
c     ----- current fock matrix is in -h0-
c           get previous fock matrices (h1,h2,h3) -----
c
  140 call rdisk(h1,h2,h3,ndaf,numscr,l2)
      if (iter .gt. 2) go to 160
c
c     ----- iter = 2 -----
c
      if (odamph) go to 320
      if ( .not. oshift) go to 420
      if ((zscftp.eq.zgvb.or.zscftp.eq.zgrhf).and. odiis) go to 420
      dampp = done
      go to 320
c
c     ----- iter > 2 -----
c
  160 if (oshift) go to 220
      if (odampr) go to 200
      if ( .not. odamph) go to 420
      if ( dabs(dee) .gt. exttol) go to 180
      if (dee .gt. dzero .and. kcount .gt. iextin) go to 180
      go to 420
  180 dmptol = dmptol/shrnkf
      exttol = exttol/shrnkf
      go to 300
  200 if ( dabs(dee) .gt. dmptol) go to 300
      if (dee .gt. dzero) go to 300
      if (dampp .gt. dmpmin) go to 300
      go to 420
  220 if ( .not. oextpr) go to 260
      if ( dabs(dee) .gt. exttol) go to 240
      if (dee .gt. dzero .and. kcount .gt. iextin) go to 240
      go to 420
  240 exttol = exttol/shrnkf
      dmptol = dmptol/shrnkf
      rshift = rsmin
      iterl = 0
      if (odamph) go to 300
      go to 620
  260 if ( dabs(dee) .gt. dmptol) go to 280
      if (dee .gt. dzero) go to 280
      if (odamph .and. dampp .gt. dmpmin) go to 280
      if ( .not. oextra) go to 280
      if (rshift .ge. vshtol) go to 280
      if (iterl .eq. 0) go to 280
      go to 420
  280 if (odamph) go to 300
      if (iterl .lt. 2) go to 620
      if ( .not. oextra) go to 620
      go to 440
c
c     ----- davidson's damping -----
c
  300 if (kcount .lt. iextin) go to 320
c
      if ( .not. oshift .or. iterl .ge. 2) go to 360
 320  temp=done/(done+dampp)
      do 340 i = 1,l2
 340  h0(i) = (h0(i)+dampp*h1(i))*temp
      go to 400
c
c     ----- iter > 2 , damping -----
c
  360 cutoff = pt5-damp00
      temp=done/(done+dampp)
      do 380 i = 1,l2
380   h0(i) = (h0(i)+dampp*h1(i))*temp
      if ( .not. oextra .or. cutoff .lt. dzero) go to 400
      odampr = .true.
      oextpr = .true.
c
      go to 460
  400 if (ncall .ne. 1) go to 660
      odampr = .true.
      oextpr = .false.
      go to 660
  420 if ( .not. oextra) go to 620
c
c     ----- pople's extrapolation procedure
c
      if (ncall .ne. 1) go to 460
      odampr = .false.
      oextpr = .true.
  440 dampp = dzero
c
c     ----- skip to end if first cycle or after extrapolation -----
c
  460 call wdisk(h0,h1,h2,ndaf,numscr,l2)
      call vsub(h2,1,h3,1,h3,1,l2)
      call vsub(h1,1,h2,1,h2,1,l2)
      call vsub(h0,1,h1,1,h1,1,l2)
      if (kcount .lt. iextin .or. iter .lt. 4) go to 680
c
c     ----- find displacement dp1,dp2,dp3 -----
c
      if (ityp .eq. 2) go to 500
      sp11 = tracep(h1,h1,l1)
      sp12 = tracep(h2,h1,l1)
      sp13 = tracep(h3,h1,l1)
      sp22 = tracep(h2,h2,l1)
      sp23 = tracep(h3,h2,l1)
      sp33 = tracep(h3,h3,l1)
      go to 520
  500 continue
      sp11 = ddot(l2,h1,1,h1,1)
      sp12 = ddot(l2,h2,1,h1,1)
      sp13 = ddot(l2,h3,1,h1,1)
      sp22 = ddot(l2,h2,1,h2,1)
      sp23 = ddot(l2,h3,1,h2,1)
      sp33 = ddot(l2,h3,1,h3,1)
  520 continue
      dp1 = dsqrt(sp11)
      dp2 = dsqrt(sp22)
c
c     ----- find cosine of angle between successive displacements -----
c
      dp3 = dsqrt(sp33)
c
c     ----- find cosine of angle between -dp(3)- and
c           plane of =dp(1)- and -dp(2)-.
c
      cosphi = sp12/(dp1*dp2)
      r = sp11*sp22-sp12*sp12
      p = (sp13*sp22-sp12*sp23)/r
      q = (sp23*sp11-sp12*sp13)/r
c
c     ----- do not extrapolate unless -4- consecutive points are
c           nearly coplanar -----
c
      cospsi = dsqrt(p*p*sp11+q*q*sp22+two*p*q*sp12)/dp3
      if (cospsi .le. tol) go to 680
c
c     ----- express -dp(1)- as x*dp(3)(projected)+y*dp(2) -----
c
      if (dampp .gt. dmpmin) go to 680
      q = -q/p
c
c     ----- test if 2*2 matrix has real eigenvalues
c           between -tol/2 and +tol/2 -----
c
      p = done/p
      pq = q*q+four*p
      if (pq .lt. dzero) go to 680
      pq =  dabs(q)+dsqrt(pq)
c
c     ----- if -4- point extrapolation is not possible,
c           try -3- point
c
      if (pq .le. tol1) go to 560
      if ( dabs(cosphi) .le. tol2) go to 680
      p = dp1/(dp2*cosphi-dp1)
      call daxpy(l2,p,h1,1,h0,1)
      go to 600
  560 ppp = p/(done-p-q)
      qqq = (p+q)/(done-p-q)
      do 580 i = 1,l2
      h0(i) = h0(i)+ppp*h2(i)+qqq*h1(i)
  580 continue
c
  600 kcount = 0
      call wrt3(h0,l2,ndaf,numscr)
      go to 680
c
c     ----- no damping or extrapolation -----
c
  620 if (ncall .ne. 1) go to 640
      odampr = .false.
      oextpr = .false.
c
c     ----- save new (modified) fock matrix -----
c
  640 dampp = dzero
c
  660 call wdisk(h0,h1,h2,ndaf,numscr,l2)
  680 continue
      return
      end
      subroutine rhfcld(q0,q,coord,atmchg,lword,nss)
c ===
c === 6-triangle version (with parallel extensions)
c ===
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
_IF(drf)
cafc ##################################################################
INCLUDE(../drf/comdrf/sizesrf)
cafc ##################################################################
_ENDIF
c     ----- closed shell hf-scf calculation -----
c           c.c.j. roothaan, rev.mod.phys. 23, 69 (1951)
c
c     irest = 3 ..... restart scf
c     nprint = 5   ... mo*s + convergence data printed for each
c                      iteration.
c              0   ... normal printing after convergenge or at exit
c
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/atmol3)
      common/blkin/scr(1)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/disc)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
INCLUDE(common/runlab)
INCLUDE(common/dnfnw)
INCLUDE(common/psscrf)
INCLUDE(common/xfield)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/symtry)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     + iposit(20),nsti(2),ondiis,junkj
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c           run through all code section
INCLUDE(common/nodeio)
INCLUDE(common/parcntl)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/drive_dft)
INCLUDE(common/statis)
_ENDIF
      Logical use_symmetry, force_serial
      Integer n_irreps
c
       integer domo, vmo
       logical logtest
       common/testerp/tester(200), idomo(200), ivmo(200),
     +               testerd(200),idomod(200), ivmod(200),
     +               testdum(2,200),idomodum(2,200),
     +               logtest(200)
c
INCLUDE(common/zorac)
INCLUDE(common/fermidirac)
c
c -----    common blocks for scrf         -----
c
INCLUDE(common/scrf)
INCLUDE(common/gvalue)
c
      dimension q(*),q0(*),coord(3,*),atmchg(*)
      dimension tmol(3),tele(3),tnuc(3)
c     dimension osign(2)
      character*10 charwall
c
INCLUDE(common/drfopt)
_IF(drf)
cafc ##################################################################
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfbem)
cafc
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
      common/enrghf/enhh,etothh,ehfhh
      logical odrf, oarf
cafc ##################################################################
_ENDIF
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data bigd/1.0d4/
      data small/1.0d-16/
      data m3/3/
            call check_feature('rhfcld')
_IF(drf)
c
cahv
        odrf = .false.
        oarf = .false.
cahv
_ENDIF
_IF(ccpdft)
      idum = CD_update_geom(c)
_ENDIF
c
c     variables for monitoring tester + convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmn = 0.0d0
      esmear = esmear_start
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      dft_accu = 1.0d0
      out = nprint .eq. 5
      outon = nprint .ne. -5
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)
_IF(parallel)
      lscf = 20+12/nav
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
_ENDIF
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = l2+l2
      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
c
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
      i50 = 1
      i60 = i50+l2
      i80 = i60+l2
      i90 = i80+l1
      ix  = i90+l1
      isiz=0
      if(ocryst)isiz=l2
      i10 = ix + isiz
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      last = i10+nss*l4-1
c
_IF(drf)
cahv --- DRF extension ---                                                   
      if (oreact) then
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
cafc
c-----  set acur for use in 2-electron rf routines drfhamx
c
        acur = acurcy
        if (iarfcal .eq. 1) then
          oarf = .true.
        else
          odrf = .true.
        endif
c
c      difference density at i70
c      ovlap at i110  from da10, 12, l2
c      dipx  at i120  from da10, 53, l2
c      dipy  at i130  from da10, 54, l2
c      dipz  at i140  from da10, 55, l2
c      omega(s) at i150    da31, 50, nomega
c      omega(op) at i160   da31, 51, nomega, if itwoeps=1
c      iexp  at i170       da31, 2,  l2
c      ijbit(s) at i180    da31, 56, l2
c      ijbit(op) at i190   da31, 57, l2, if itwoeps=1 else clearit
c      fdrf      at i200   calc in drfhamc, l2
cafc
        i110 = last + 1 
        i120 = i110 + l2
        i130 = i120 + l2
        i140 = i130 + l2
        i150 = i140 + l2
        if (odrf) then
          i160 = i150 + nomga
          if (itwoeps .eq. 1) then
            i170 = i160 + nomga
          else
            i170 = i160 + 1
          endif
        else
c
c      wt and vr at i150
c      relay at i160
c
          i160 = i150 + nwtc*nwtr
          i170 = i160 + ndim*ndim
        endif
        i180 = i170 + l2
        i190 = i180 + l2
        i200 = i190 + l2
        i210 = i200 + l2
        last = i210 + l2
c
c       i110 = igmem_alloc(l2)
c       i120 = igmem_alloc(l2)
c       i130 = igmem_alloc(l2)
c       i140 = igmem_alloc(l2)
c       i150 = igmem_alloc(nomga)
c       i200 = igmem_alloc(l2)
c       i210 = igmem_alloc(l2)
c       if (itwoeps .eq. 1) i160 = igmem_alloc(nomga)
c       i170 = igmem_alloc(l2)
c       i180 = igmem_alloc(l2)
c       i190 = igmem_alloc(l2)
      endif
cahv
      endif
_ENDIF
c
      if(last.gt.lword) then

      write(iwr,9309)last,lword
      call caserr('insufficient memory available')
      endif
      if (out) write (iwr,9308) i50,i60,i80,i90,i10,i20,i21,
     +     i30,i31,i40,last
c
c     ----- occupation numbers -----
c
c     -o- at q(i80)
c
c     ----- initialize variables -----
c
_IF(drf)
cahv --- DRF extension ---                                                   
      if (odrf) then
cafc
c  -----  read necessary data from da10
        call daread(idafh,ioda,q(i110),l2,12)
        call daread(idafh,ioda,q(i120),l2,53)
        call daread(idafh,ioda,q(i130),l2,54)
        call daread(idafh,ioda,q(i140),l2,55)
cafc
c  -----  read necessary data from da31
        call daread(idafdrf,iodadrf,q(i170),l2,2)
        call daread(idafdrf,iodadrf,q(i150),nomga,50)
        call daread(idafdrf,iodadrf,q(i180),l2,56)
cafc
        if (itwoeps .eq. 1) then
          call daread(idafdrf,iodadrf,q(i160),nomga,51)
          call daread(idafdrf,iodadrf,q(i190),l2,57)
        else
          call clear(q(i190),l2)
        endif
        call clear(q(i200),l2)
cafc
      endif
cahv
_ENDIF
      if (maxcyc .le. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
       tester(loop) = 0.0d0
       testerd(loop) = 0.0d0
       logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1
      i=lensec(l2)
      iblkh0 = ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)

      if(numdis.eq.0)then
         numdis=num8
         ndafd=iblkh+i
      else
         ndafd=ibldis+2
      endif
      ndafdi=ndafd+3*i
      if(irest.ne.3.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         kcount = 0
         iterv = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
c
      diffdp = bigd
      diffdpp = bigd
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3

      ovir = .false.
      if(oaimpac())ovir = .true.
c
      if (odamph .or. oshift) damp = done
c
      lockt=lock
      tim0=cpulft(1)
_IF(drf)
cdrf1
c    ---- drf: 2-el. reactionfield contributions
cdrf1
c     only if really 'scf'type and if intdrf=0 
c     (intdrf=2 2-electron DRF integrals are added to the 
c       vacuum 2-electron integrals on disk)
c     change core partitioning to likes of drf.
c     to store are:
c       ovlap (=i110)
c       dipx dipy dipz (i120,130,140)
c       omega(s)  omega(op) (i150,160)
c       iexp (i170)
c       ijbit(s)  ijbit(op) (i180,i190)
c      check if itwoeps .eq. 1
c    itwoeps: flag for calculation of two relay-matrices with subsequent
c            solnam for surching parameters) only if itwoeps .eq. 0
c            create space for them in partitioning
c            endif
c    read from da10 (idafh=10) records 12,53,54,55
c    in gamess= call sec192(ss,h-core,ekinet,l2,ncall) intega.f:3476+
c
c       rec.11    h-core matrix      inttel.f:376
c     with gamess= q(i20)
c       rec.12    ss (overlap mat)   inttel.f:375
c     with gamess= q(i10)
c       rec.16    density
c       rec.41    ekinet matrix      inttel.f:377
c     with gamess= q(i30)
c       rec.53-55 dipole moments     intnew.f:175
c       rec.45    xexp \ expansion centra (zmat ?)
c       rec.46    iexpc/
c       rec.97    blocked overlap s mat inttel:379
c       rec.210-213 (rfcal.f:740)
c    read from da31 (idafdrf=31) records 2,50,56 itwoeps:51,57
c       rec.1     assignation
c       rec.2     expansion centra
c       rec.50    ilomg,  omega matrix ieps=0 kappa=kappa1 (static)
c       rec.51    ilomg,  omega matrix ieps=1 kappa=kappa2 (optic )
c       rec.56    ilij,  (ij/ij)  ieps=0 (static)
c       rec.57    ilij,  (ij/ij)  ieps=1 (optic)
c     created by drf routines: (in da31)
c       rec.12    interaction integrals with external charges
c       rec.13    reverse dispersion contrib. to interac. en.
c       rec.14    non-equil potential
c       rec.58-   interaction self to be contracted with dens.mat
c       rec.61-66 second moment integrals
c       rec.70    (start)density matrix
c    if itwoeps .eq. 1 then read some more and clear some
cdrf1(i.e. skip something)
_ENDIF
c
c     ----- nuclear energy
c     ----- psscrf dvsg
c
      if (opssc) then
         en = denuc(nat,atmchg,coord,itotbq,iwr)
      else if (ocryst) then
         en = crnuc(nat,czan,coord,ecrn)
      else
         en = enucf(nat,atmchg,coord)
      endif
c
c     l0 = number of canonical vectors kept.
c
      l0=newbas0
c
      n_irreps = 0
      Do i = 1, l1
         isymmo( i ) = isymao( i )
         n_irreps = Max( n_irreps, isymao( i ) )
      End Do
      use_symmetry = n_irreps .GT. 1  .and . symm_diag
      if (outon) then
       write (iwr,9008) en
       if(use_symmetry) write(iwr,9009)
       if(odebug(31)) then
         write (iwr,9028) maxcyc,mconv,nconv,mpunch
_IF1(tcx)     +                    ,nss
       else
         write (iwr,9029) maxcyc,mconv,nconv,mpunch
_IF1(tcx)     +                    ,nss
       endif
      endif
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i10),l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,q(i20),q(i10))
      endif
_ENDIF
c
c     ----- compute canonical orthonormal vectors -----
c
c     -s - at q(i10)
c     -sv- at q(i20)
c     -se- at q(i31)
c
       call rdedx(q(i10),l2,ibl7st,num8)
       call qmat_symm(q,q(i10),q(i20),q(i31),q(i40),iky,l0,l1,l3,l1,out,
     +               isymmo, use_symmetry )
c
c     ----- compute initial guess density matrix -----
c
c     -d- at q(i30)
c
c     ----- orthogonalise vectors
c
c load guess vectors..
      call rdedx(q(i30),l3,ibl3qa,idaf)
      if (.not.onoor) then
c load overlap
      call rdedx(q(i50),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i50), num, l0)
       else
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     +    ilifq,l0,l1,1)
       endif
_ELSE
      call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
      call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     +     ilifq,l0,l1,1)
_ENDIF
      else
        write(iwr,'(2a)') '** warning start-orbitals not orthogonal **',
     1                    '** result may not be variational         **'
c...    if need be one can pick up the occupation numbers here to
c...    calculate coulomb energies for correlated monomers (jvl,2003)
      end if
      call wrt3(q(i30),l3,ibl3qa,idaf)
      call tdown(q(i30),ilifq,q(i30),ilifq,l1)
      call rdedx(q(i90),l1,ibl3ea,idaf)
      ndoc = na
      if (osmear) then
        call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                   iwr,oprint(60))
        emrmn = 2*emrmn
      else
         call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
      call wrt3(q(i50),l2,ibl3pa,idaf)
c...
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      lprnt = l1
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      lprnt=min(lprnt,l0)
_IF(drf)
cahv --- DRF extension ---                                                   
cafc  write initial density to da31#70 :
cafc  drf expects the diffdens to be written on this spot, but that's
cafc  hardly possible at this stage, isn't it ;-)
cafc  to be used in calculating contributions to 2e-part of Fock in SCF
cafc
      if (odrf .and. igetden.eq.1) 
     1 call daread(idafh,ioda,q(i50),l2,16)
      if (odrf) call dcopy(l2,q(i50),1,q(i210),1)
cahv
_ENDIF
_IF(ccpdft)
      if (CD_active()) then
         call retrieve_spare(imemspare)
         imemfree  = igmem_max_memory() - imemspare 
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q0,q0,iwr)
         else if (ks_bas.eq.KS_MO) then
            imemreq = CD_memreq_energy_mo(l1,na,0,q0,q0,iwr)
         else if (ks_bas.eq.KS_AOMO) then
            imemreq = CD_memreq_energy(l1,na,0,q0,q0,iwr)
         endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            write(iwr,600)ierror
            call caserr('Out of memory in incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
c
c     ----- find the time remaining at the start -----
c
      call timrem(tlefts)
c
c     ----- start scf procedure -----
c     ----- construct a skeleton fock matrix -----
c
c     -d- at q(i50)
c     -h- at q(i10)
c  -exch- at q(i20)
c
 140  continue
_IF(parallel)
c***   ***node-MPP***
      call pg_synch(4444)
_ENDIF
c
_IF(ccpdft)
c
c Modify fock builder options
c
      idum = CD_set_2e()
_ENDIF
c
_IF(cray,convex,titan)
      call hstar(q(i50),q(i10),q(i20),nopk,nss)
_ELSE
      call hstar(q(i50),q(i10),q(i20),nopk)
_ENDIF
_IF(ccpdft)
c
c restore fock builder options
c
      idum = CD_reset_2e()
_ENDIF

_IF(parallel)
c***   ***node-MPP***
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
c
      if (out) then
      write (iwr,9068)
      call prtril(q(i10),l1)
      endif
_IF(drf)
cahv --- DRF extension ---                                                   
cafc  it's about time to add the drf-contributions to the fock-matrix
cafc  this is done by the routine drfhamc, which uses the diff dens
cafc  by reading it from da31 and places the results to be added at i200
cafc
      if (odrf) then
        call drfhamc(q(i200),q(i210),q(i110),q(i120),q(i130),
     1    q(i140),q(i170),q(i150),q(i180),q(i160),q(i190))
cafc
c --- add to fock matrix   (scffact * fdrf(i) )
        do 1111, i = 1, l2
          q(i10-1+i) = q(i10-1+i) + scffact*q(i200-1+i)
 1111   continue
c
        call dcopy(l2,q(i50),1,q(i210),1)
      endif
cahv
_ENDIF
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at q(i10)
c     scratch area at q(i20)
c
      call symh(q(i10),q(i20),iky,0,0)
      if (out) then
          write (iwr,9088)
          call prtril(q(i10),l1)
      endif
      if (.not.oscrf) then
c
c        ----- read in core hamiltonian matrix
c              and calculate hf energy -----
c
c        -h0- at q(i20)
c        - h- at q(i10)
c        - e- at q(i30)
c
         if (ozora.and..not.ocvged) then
c
c           sf compute zora contributions to one-electron integrals
c              needs density matrix at q(i50)
c
            call rdedx(q(i20),l2,ibl7f,num8)
            if ((mod(iter,niter_z).ne.0.and.iter.ne.999999999)
     1          .or.nat_z.ne.0) then
               call zora(q0,q(i50),q(i20),'calc')
         else
            call zora(q0,q(i50),q(i20),'force')
         end if
c
      else
         call rdedx(q(i20),l2,ibl7f,num8)
c...     add zora corrections (only ato,mic or restore)
         if (ozora) call zora(q0,q(i50),q(i20),'read')
      endif
c
      call vadd(q(i10),1,q(i20),1,q(i10),1,nx)
c     save previous total
      ehf0 = ehf
c
c     ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i50),q(i20),l1)
c
_IF(ccpdft)
      if(CD_active())then
c
c        ccpdft: evaluate Kohn Sham energy expression
c
         if(CD_HF_exchange() .or. CD_HF_coulomb())then
c
c           Coulomb or exchange operator is in q(i10) augmented by H_1,
c           compute energy using HF expression (maybe without exchange)
c
            ehf2 = tracep(q(i50),q(i10),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c           Coulomb operator has not yet been evaluated
c           (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1
         endif
c
c        Update Kohn-Sham matrix and compute fitted/integrated
c        energy terms
c
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(q0(isma),1,l2)
_ENDIF
         if (ks_bas.eq.KS_AO) then
            idum = CD_energy_ao(c,q(i10),dum,q(i50),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_MO) then
            idum = CD_energy_mo(l1,na,0,c,q(i10),dum,q(i30),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_AOMO) then
            idum = CD_energy(l1,na,0,c,q(i10),dum,q(i30),dum,q(i50),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else
            call caserr("rhfcld: illegal ks_bas")
         endif
_IF(debug_S)
         call compare_S(q0(isma),q0(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF
         call symm_op(q(i10))
 
         ehf = etmp+edft
 
         if (out) then
          if(opg_root()) then
            call CD_print_dftresults(.true.,.false.,iwr)
            write(iwr,9407)ehf1, ehf2
          endif
         endif
         call timana(31)
         call cpuwal(begin,ebegin)
 
      else
c
c        Hartree-Fock energy expression
c        E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                      ehf1            ehf2
c
         ehf2 = tracep(q(i50),q(i10),l1)
         ehf = (ehf1+ehf2)*pt5
      endif
_ELSE
        ehf2 = tracep(q(i50),q(i10),l1)
        ehf = (ehf1+ehf2)*pt5

_ENDIF

        if(ocryst)then
           m53=53
           call secget(isecfx,m53,iblk)
c           write(6,*)'secget',isecfx,m53,iblk
           call rdedx(q(ix),l2,iblk,idaf)
           ehf3 = tracep(q(i50),q(ix),l1)
c           write(6,*)'ehf3',ehf3
           esave = ehf
           ehf = (ehf1+ehf2-ehf3)*pt5
           ecre = esave - ehf
c           write(6,*)'new ehf, ecre',ehf,ecre
        endif
c
      else
c
c ----- start of scrf procedure
c ----- write 2e matrix to ed7
c
        call wrt3(q(i10),l2,iblkh,num8)
c
c ----- calculate dipole moment
c ----- read in dipole integrals
c
        call rdedx(q(i10),l2,ibl7x,num8)
        call rdedx(q(i20),l2,ibl7y,num8)
        call rdedx(q(i30),l2,ibl7z,num8)
c
c ----- construct electronic part
c
        tele(1) = -tracep(q(i50),q(i10),l1)
        tele(2) = -tracep(q(i50),q(i20),l1)
        tele(3) = -tracep(q(i50),q(i30),l1)
c
c ----- and nuclear part
        tnuc(1) = ddot(nat,czan,1,c(1,1),3)
        tnuc(2) = ddot(nat,czan,1,c(2,1),3)
        tnuc(3) = ddot(nat,czan,1,c(3,1),3)
c ----- and totals x,y,z
        do iq=1,3
          tmol(iq)=tnuc(iq)+tele(iq)
        enddo
c
c ----- total dipole moment
c
        dtot =ddot(m3,tmol,1,tmol,1)
        if(dtot.gt.small)dtot=dsqrt(dtot)
c ----- scale dipole integrals by g's and dipole moments
        do loop=1,l2
           q(i10+loop-1)=q(i10+loop-1)*tmol(1)*gx
           q(i20+loop-1)=q(i20+loop-1)*tmol(2)*gy
           q(i30+loop-1)=q(i30+loop-1)*tmol(3)*gz
        enddo
c
c ----- now add x,y,z components into i30
c ----- leaves i10,i20 for fock matrices
        do loop=1,l2
         q(i30+loop-1)=q(i30+loop-1)+q(i20+loop-1)+q(i10+loop-1)
        enddo
c
c ----- read back fock matrices
        call rdedx(q(i10),l2,iblkh,num8)
        call rdedx(q(i20),l2,ibl7f,num8)
c...    add zora corrections (only atomic / restore)
        if (ozora) call zora(q0,q0,q(i20),'read')
c
c ----- add dipole terms to 1e matrix
        do loop=1,l2
         q(i20+loop-1)=q(i20+loop-1)+q(i30+loop-1)
        enddo
c
c ----- do hf energy calculation
        call vadd(q(i10),1,q(i20),1,q(i10),1,l2)
        ehf1 = tracep(q(i50),q(i20),l1)
        ehf2 = tracep(q(i50),q(i10),l1)
        ehf0 = ehf
        ehf = (ehf1+ehf2)*pt5
c ----- energy expression for scrf
        ehf = ehf + gx*(dtot**2)*pt5
c ----- now allow for charged systems
        ehf = ehf + ((1.0d0-dielec)/(dielec*aradix))*(ich**2)
c
      endif
c
c compute kinetic energy and virial
c
      ehf = ehf + 0.5d0*emrmn
      if(ovir)call rhfvir(q(i50),ehf,en,ek,vir,num8,l1,l2)
c
c      ----- save fock matrix
c
c     -h- at q(i10)
c
      call wrt3(q(i10),l2,iblkh,num8)

      if (out) then
          write (iwr,9048)
          call prtril(q(i10),l1)
      endif
c
      iter = iter+1
      etot0 = etot
      etot = ehf+en
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
c
c     ----- damp and extrapolate hamiltonian matrix -----
c
c     -h - at q(i10)     hamiltonian matrix (n th iteration)
c     -ho- at q(i20)     old h matrix (n-1 th iteration)
c     -ha- at q(i30)     ancient h matrix (n-2 th iteration)
c     -hp- at q(i40)     prehistoric h matrix (n-3 th iteration)
c
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
c
      if (damp .lt. dmpcut) damp = dmpcut
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafd,
     +numdis,iterv,1,1)
      diffpp=diffp
      diffp=diff
      if(odiis) then
_IF(ga)
      if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
         call diis_ga(q(i10),q(i30),q(i50),ibl3qa,diffd,domo,vmo)
      else
         call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     +   ndafdi,numdis,diffd,lockt,domo,vmo)
      endif
_ELSE
       call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     + ndafdi,numdis,diffd,lockt,domo,vmo)
_ENDIF
      if (iter.le.maxcycp) then
       idomod(iter) = domo
       ivmod(iter) =  vmo
       testerd(iter) = diffd
      endif
      endif
c
c ------ take precautions if tester is increasing again
c
      if(iter.eq.1 .or. .not.odiis) goto 223
      if(ondiis) then
        if (oprind) then
         write(iwr,9824)diffd,diffdp
         write(iwr,9823)diffd, acurcy*diisf
         write(iwr,*) 'de = ', de
        endif
        if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.acurcy*diisf
     1  .or.(otestd.and.de.lt.1.0d-10)) then
*        if (de.lt.1.0d-5) then
          if(iter.le.maxcycp) logtest(iter) = .true.
          diffdpp = diffdp
          diffdp = diffd
          go to 223
*        else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
*         odynamic = .true.
*         accdi1 = 0.0d0
*         itdiis = 0
*         if(oprind) write(iwr,9828) iter
*        endif
        else
          if(oprind) write(iwr,9829) iter
        endif
        nsti(1)=0
        ondiis=.false.
        call rdedx(q(i10),l2,iblkh,num8)
        diffdp = bigd
        diffdpp = bigdp
      else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c       if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c        iwild = iwild + 1
c        osign(iwild) = de.lt.0.0d0
c        if (iwild.eq.2.and..not.oreset.and.
c    +       .not.(osign(1).and.osign(2)) ) then
c         gapa1 = 5.0d0
c         gapa2 = gapa1
c         oreset = .true.
c         iwild = 0
c         write(iwr,9831) gapa1
c         maxcyc = maxcyc + maxcyc
c         maxcycp = min(200,maxcyc)
c        endif
c       endif
c
        if (odynamic) then
         if(diff.le.diffp.or.diffp.le.diffpp) then
          itdiis = itdiis + 1
          if (itdiis.eq.3) then
            accdi1 = diff
            if (oprind) write(iwr,9826) accdi1
          endif
         else
          if (oprind) write(iwr,9827)
          itdiis = 0
         endif
        endif
      endif
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -q- at q(i30) orthonormalizing transformation
c     -h- at q(i10)
c    -h'- at q(i50) transformed h matrix
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
 223  call rdedx(q(i30),l3,ibl3qa,idaf)
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao, 
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diff=0.0d0
      ii=ndoc+1
      if(ii.le.l0)then
        do 250 i=ii,l0
        ij=iky(i)+i50
        loop=idamax(min(i-1,nocmx),q(ij),1)
        if(loop.gt.0) then
         dum = dabs(q(ij+loop-1))
         if (dum.gt.diff) then
           diff = dum
           if(iter.le.maxcycp) then
            tester(iter) = diff
            idomo(iter) = loop
            ivmo(iter) = i
           endif
         endif
        endif
 250    continue
      endif
      dmplim = dmax1(dmpcut,2.0d0)
      if(ondiis.and..not.osmear)diff=diffd
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy) .and. (iter.gt.1)
      if (ocvged) then
c
c sf perform one extra scf-cycle after convergence;
c    to compute 'effective scaled' fock operator.
c
       if ( ozora .and. oscalz) then
          id_z =  igmem_alloc(l2*2)
c sf compute scaling factors 
          call rdedx(q0(id_z),l2,ibl3pa,idaf) 
          call zora(q0,q0(id_z),q0(id_z+l2),'scale')
c sf perform scaling 
          call scale_z(q0,q(i50),q0(id_z),q(i30),
     1                 q0(id_z+l2),ehf,etot,ne,num,0) 
c 
          call gmem_free(id_z) 
          goto 666 
        endif 
c 
        go to 321 
       endif 
666   continue 
c 
c            shift the diagonal of the transformed 
c           h matrix by rshift for the virtual part.  
c 
      rshift=0.0d0 
      if(.not.ondiis.or.olevd) 
     * call shiftq(q(i50),nocmx,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2) 
c 
c     ----- diagonalize new hamiltonian matrix ----- 
c 
c     -h- at q(i50) 
c     -v- at q(i30) 
c     -d- at q(i90) 
c 
      m2=2 
      diaacc = diff*5.0d-3 
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5 
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11 

      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
c
c     ----- back-transform the eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i41)
c     -q- at q(i10)
c     scratch area at q(i50)
c
      call rdedx(q(i10),l3,ibl3qa,idaf)
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
c
c     ----- if vshift is true, reorthogonalize the vectors and
c           thereby the transformation matrix every igs th iteration.
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
      ig = mod(iter,igs)
      if (ig .eq. 0) then
         call rdedx(q(i10),l2,ibl7st,num8)
c
_IF(ga)
         if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
            call porth(q(i30),q(i10), num, l0)
         else
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
         endif
_ELSE
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
         call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
      endif
      call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
c     ----- form density matrix -----
c
c     -v- at q(i30)
c     -o- at q(i80)
c     -d- at q(i50)
c
      ig=nocmx
      ndoc = na
      if (osmear) then
         esmear = max(esmear_final,
     &                min(esmear,esmear_scale*(diff-acurcy)))
         call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                    iwr,oprint(60))
         emrmn = 2*emrmn
      endif
      if(ig.lt.l0)
_IFN1(civu)     * call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civu)     * call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      if (.not.osmear) then
         call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
c
c     print vectors
c
      if (out) then
         write (iwr,9148)
         call prev(q(i30),q(i90),lprnt,l1,l1)
      endif
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      if (out) then
         write (iwr,9168)
         call prtril(q(i50),l1)
      endif
c
c     compute RMS convergence on density
c
      if (iter.ne.1) then
       call rdedx(q(i10),l2,ibl3pa,idaf)
       dsq = cvgdens(q(i50),q(i10),l2)
       dsq = dsq / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
c
c     ----- save mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i50)
c     -d- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i50),q(i90),q(i80),ndaf,l1,l2,
     *ibl3pa,ibl3ea)
      if(ozora) then
c he print dens
      print *,'dens in rhf'
      call rdedx(q(i50),l2,ibl3pa,ndaf)
      call prtri(q(i50),l1)
      end if
c
 321  if(numdis.ne.num8) then
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
      endif

      tim1=cpulft(1)
      delt = tim1-tim0
      tim0 = tim1
_IF(drf)
cafc ##################################################################
cafc  calculate the diff dens (old at i70, new at i50)
cafc  save this to da31#70
cafc
      if((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
         diff2=cvgden(q(i50),q(i210),l2)
c        call dawrit(idafdrf,iodadrf,q(i70),l2,70,navdrf)
      endif
cafc ##################################################################
_ENDIF
c
c     monitor static tester with diis on - at the moment remove diis
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if (oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif
 
c
c     ----- check and print convergence behavior -----
c
      otrigp = maxcyc-iter.lt.10
      if(outon.or.otrigp) then
       if(odebug(31)) then
         write (iwr,9188) iter,kcount,etot,
     +    ehf,de,diff,rshift,damp,derror,delt,tim1,yavr
       else
         write (iwr,9189) iter,kcount,etot,
     +    ehf,de,diff,rshift,damp,derror,yavr
       endif
      endif
_IF(parallel)
      else
          iter = iter + 1
      endif
******
      if(ipiomode .eq. IO_NZ_S)then
         call pg_brdcst(7125,maxcyc,lscf*8,0)
         call pg_brdcst(7123,q(i50),l2*8,0)
         call pg_synch(4444)
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif
_ENDIF
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy) .and. (iter.gt.1)
c
c     now flag covergence if STATIC tester has been encountered
c
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
c
c  update files in core or GAs
c
      if (ioupd.gt.1) call ioupdate
cgdf  
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if
c
      if (ocvged) go to 400
c
c     ----- exit in case of time limit -----
c     ----- punch out restart data -----
c
      call timit(0)
      call timrem(tlefti)
      if (tlefti .gt. skale*(tlefts-tlefti)/iter) go to 380
      irest=3
      write(iwr,9367)
      call texit(0,irest)
      go to 400
  380 if (iter .lt. maxcyc) go to 140
       if(omaxcyc) then
         write (iwr,9289)
         ocvged = .true.
       else
         write (iwr,9288)
c        write out error description to xml/punchfile
         call blkerror('excessive number of SCF iterations',0)
         etot = dzero
         irest=3
         ehf0=  ehf
         ehf = -en
       endif
  400 continue
c
c     ----- print actual value of energy -----
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
      lock = lockt
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr(
     +         'Memory failure in rhfcld:CD_jfit_clean2')
         endif
      endif
_ENDIF
      if (ocvged.and.outon) write (iwr,9208)
_IF(drf)
cafc ##################################################################
cafc  supposedly the hh-variables are being used in analyzing drf-
cafc  contributions later on (no good this way !)
cafc
      if (field.ne.' ') then
        ehfhh=ehf
        enhh=en
        etothh=etot
      endif
cafc ##################################################################
_ENDIF
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if(omodel) then
       write (iwr,9371) iter,tim1,charwall(),ehf,en,etot,diffdens
      else if (ocryst) then
       write (iwr,9372) iter,tim1,charwall(),ecre,ecrn,ecre+ecrn,ehf,en,etot,
     +                  diffdens
      else if (ovir) then
       write (iwr,9365) iter,tim1,charwall(),ek,ehf,en,etot,vir,diffdens
      else
       write (iwr,9368) iter,tim1,charwall(),ehf,en,etot,diffdens
      endif
c     write out energies to xml file      
      call blkscfe(ehf,en,etot)
      if(yavr.ne.yblnk)write(iwr,9369)
c
_IF(drf)
       if(.not.oreact) then
_ENDIF
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       write(iwr,9199)
       idum = min(iter,200)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       else
       write(iwr,9198) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       endif
       enddo
       write(iwr,9197)
       endif
_IF(drf)
       endif
_ENDIF
c
      if(numdis.eq.num8)numdis=0
      if (nprint .eq. 5 .or. nprint .eq. -5) go to 440
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call rdedx(q(i30),l2,ibl3pa,idaf)
      call rdedx(q(i21),l1,ibl3ea,idaf)
_IF(drf)
cafc ##################################################################
cafc  final density at i30 is saved to da10#16 (why ?!)
cahv  FOR ENERGY ANALYSIS!
cafc
      if (field .ne. ' ') call dawrit(idafh,ioda,q(i30),l2,16,navh)
cafc ##################################################################
_ENDIF
      call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
      call tdown(q(i10),ilifq,q(i10),ilifq,l0)
      if(otsym) 
     +   call symass(q(i10),q(i21),q(i80),q0)
      if(oprint(25))go to 440
      write (iwr,9148)
      call prev(q(i10),q(i21),lprnt,l1,l1)
c
c     evaluate the correlation energy using density functional
      if (odenfu) call becke(q(i10),1)
c
      if(.not.outon) then
        write (iwr,9168)
        call prtril(q(i30),l1)
      endif
  440 if(oprint(49))then
         write(iwr,9169)
         call rdedx(q(i10),l2,iblkh,num8)
         ij=i10
         do 452 i=1,l1
         scr(i)=q(ij)
 452     ij=ij+i+1
         write(iwr,9170)(scr(i),i=1,l1)
      endif
_IF(nbo)
c
c save a copy of the fock matrix on the dumpfile
c (section 242) for possible nbo analysis
c
      len42 = lensec(nx)
      m = 42
      call secput(isect(42),m,len42,iblk42)
      call rdedx(q(i10),l2,iblkh,num8)
      call wrt3(q(i10),nx,iblk42,idaf)
      lds(isect(42)) = nx
_ENDIF
_IF(ga)
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.false.)
      endif
_ENDIF
      if (mpunch .ne. 0) then
c
c     ----- punch the occupied orbitals -----
c
c     -v- at q(i10)
c
         call rdedx(q(i10),l3,ibl3qa,idaf)
         call pusql(q(i10),na,l1,l1)
      endif
      if (ocvged) irest = 0
      if(outon)then
      cpu=cpulft(1)
      write(iwr,7777)cpu,charwall()
      endif
      accdi1 = accdin
      call copq_again(mouta)
c
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapa2 = gapa1
c      maxcyc = maxcyc / 2
c     endif
c
      return
c
c9831 format(/15x,'**** increase level shifter to ',f8.2,' ****'/)
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9827 format(1x,'rising diff - de = ', f20.10)
 9825 format(/' *** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffd, diffdp = ',2f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9197 format(10x,64('+'))
 9198 format(10x,i5,f11.7,' (',2i4,') (*)',6x,f11.7,' (',2i4,')' )
 9298 format(10x,i5,f11.7,' (',2i4,')',10x,f11.7,' (',2i4,') (*)' )
 9199 format(/10x,64('+')/
     + 10x,'CONVERGENCE / TESTER ANALYSIS',12x,'From DIIS'/
     + 10x,64('+')/
     + 10x, 'Iter.', '   Tester   (domo/vmo)',
     + 10x,          '   Tester   (domo/vmo)'/)
 9008 format(/1x,'----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9028 format(/15x,'convergence data'/15x,16('=')///
     +     1x,'maximum number of iterations = ',i6/
     +     1x,'method of convergence        = ',i6/
     +     1x,'convergence criterion        =1.0e-',i2/
_IFN1(tcx)     +     1x,'punch out option             = ',i6//1x,126('=')/
_IF1(tcx)     +     1x,'punch out option             = ',i6/
_IF1(tcx)     +     1x,'vectorisation factor in h-build',i6//1x,126('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis',
     * 4x,'del(t)',6x,'time'/
     * 18x,'energy',10x,'energy',35x,'shift'/1x,126('='))
 9029 format(/15x,'convergence data'/15x,16('=')///
     +     1x,'maximum number of iterations = ',i6/
     +     1x,'method of convergence        = ',i6/
     +     1x,'convergence criterion        =1.0e-',i2/
_IFN1(tcx)     +     1x,'punch out option             = ',i6//1x,106('=')/
_IF1(tcx)     +     1x,'punch out option             = ',i6/
_IF1(tcx)     +     1x,'vectorisation factor in h-build',i6//1x,106('=')/
     * 3x,'cycle',10x,'total',5x,'electronic',8x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     * 18x,'energy',10x,'energy',35x,'shift'/1x,106('='))
 9188 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,
     1       2f10.3,a3)
 9189 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,a3)
 9048 format(20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23(
     +     '-'))
 9148 format(//1x,120('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9168 format(/20x,14('-')/20x,'density matrix'/20x,14('-'))
 9169 format(/40x,32('=')/
     *40x,'diagonal elements of fock matrix'/40x,32('=')/)
 9170 format(/5x,7f15.7)
 9208 format(/20x,16('-')/20x,'energy converged'/20x,16('-'))
 9288 format(/20x,30('-')/20x,'excessive number of iterations'/
     +     20x,30('-'))
 9289 format(/20x,30('-')/
     +        20x,'excessive number of iterations'/
     +        20x,'but flag SCF convergence'/
     +        20x,30('-'))
 9308 format(' core assignment'/' i50, i60, i80,',
     +     ' i90, i10, i20, i21, i30, i31, i40  = '/ 10i8/
     +     ' last = ', i8)
 9348 format(//40x,32('*')/40x,'closed-shell rhf scf calculation'
     +       /40x,32('*'))
 9365 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     *,f12.2,' seconds',a10,' wall',/
     *10x,14('-')/
     +10x,'kinetic energy         ',f18.10/ 
     +10x,'electronic energy      ',f18.10/ 
     +10x,'nuclear energy         ',f18.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy           ',f18.10,8x,f10.7/
     +10x,'convergence on density ',f18.10)
 9368 format(  /20x,14('-')/
     *20x,'final energies   after',i4,' cycles at ',f10.2,' seconds',
     *a10,' wall',/20x,14('-')//
     *20x,'electronic energy      ',f18.10/
     *20x,'nuclear energy         ',f18.10/
     *20x,'total energy           ',f18.10/
     +20x,'convergence on density ',f18.10)
 9371 format(  /20x,14('-')/
     *20x,'final energies   after',i4,' cycles at ',f10.2,' seconds',
     *a10,' wall',/20x,14('-')//
     *20x,'electronic energy      ',f18.10/
     *20x,'nuclear energy         ',f18.10/
     *20x,'total qm energy        ',f18.10/
     +20x,'convergence on density ',f18.10)
 9372 format(  /20x,14('-')/
     *20x,'final energies   after',i4,' cycles at ',f10.2,' seconds',
     *a10,' wall',/20x,14('-')//
     *20x,'electron-crystal energy',f18.10/
     *20x,'nuclear-crystal energy ',f18.10/
     *20x,'total-crystal energy   ',f18.10/
     *20x,'total electronic energy',f18.10/
     *20x,'total nuclear energy   ',f18.10/
     *20x,'total energy           ',f18.10/
     +20x,'convergence on density ',f18.10)
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9369 format(//20x,90('*')/
     *20x,'*',88x,'*'/20x,'*',88x,'*'/
     *20x,'*',31x,'warning  state is averaged',31x,'*'/
     =20x,'*',88x,'*'/20x,90('*')//)
7777  format(//
     *' end of closed shell scf at ',f12.2,
     *' seconds',a10,' wall',//1x,104('-')//)
9309  format(//1x,'words requested ',i8/
     *         1x,'words available ',i8/)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'=============================='/)
      end
_IFN(parallel)
      subroutine diiscd(h0,h1,h2,h3,iblkh0,iblkqa,ndaf,numscr,
     * diff0,lockt,idomo,ivmo)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension h0(*),h1(*),h2(*),h3(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis
c
c
c diisc is specifically for the closed shell  cases
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written out to either:
c (i)  ed7 - scratch (num8)
c (ii) numdis - a stream for maintaining diis information across runs
c the parameter numscr is used to specify this stream
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting block for the vectors .. need space for 20
c                                           .. 40 for uhf
c ndafd   the block which holds the beginning of the extrapolated vecs
c accdi1  diis comes in only when diff is less than this
c accdi2  the diis solution is only used when the residuals
c          are less than thia
c
c num8   the stream for ed7
c numscr the stream for saving diis information (defaults to num8)
c num    the dimension of the arrays
c nx     the size of the arrays
c
INCLUDE(common/harmon)
INCLUDE(common/prints)
c ===
c ===  6-triangle variant of diis - h0 overlays h2
c ===
      derror=0.0d0
      l3=num*num
      num0 = newbas0
      idomo = 0
      ivomo = 0
      nmin=3
      iblk=ndaf
      nmax=8
      length=lensec(nx)
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1) goto 300
      call wrt3(h0,nx,iblkh0,num8)
      call rdedx(h3,l3,iblkqa,idaf)
c
c ----- sort out the vectors
c
      call tdown(h3,ilifq,h3,ilifq,num0)
c
c ----- calculate the error vector
c
      call mult2(h3,h1,h2,num0,num0,num)
      diff0=0.0d0
      ij=0
      do 5 i=1,num0
      do 4 j=1,i
      ij=ij+1
      if(i.le.na.and.j.le.na) h1(ij)=0.0d0
      if(i.gt.na.and.j.gt.na) h1(ij)=0.0d0
      if( dabs(h1(ij)).gt.diff0) then
       diff0= dabs(h1(ij))
        idomo = j
        ivmo =  i
      endif
 4    continue
 5    h1(ij)=0.0d0
      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) goto 310
c
c ------ transform back to the ao basis
c
      do 10 i=1,num
      do 10 j=1,i
      dum=h3(ilifq(i)+j)
      h3(ilifq(i)+j)=h3(ilifq(j)+i)
 10   h3(ilifq(j)+i)=dum
      call mult2(h3,h2,h1,num,num,num)
      call rdedx(h1,nx,ibl7s,num8)
      call square(h3,h1,num,num)
      call mult2(h3,h1,h2,num,num,num)
      call rdedx(h0,nx,iblkh0,num8)
      call search(iblk,numscr)
      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)

      if(ipos1.ge.1) then
       do 50 i=1,ipos1
       ibl=iposit(i)+length
       call rdedx(h3,nx,ibl,numscr)
       mp=mp+1
       st(mp)=tracep(h1,h3,num)
 50    continue
      endif
      mp=mp+1
      st(mp)=tracep(h1,h1,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
      iposit(ipos)=iposun(numscr)
      call wrt3s(h0,nx,numscr)
      call wrt3s(h1,nx,numscr)
      ipos1=ipos+1
      if(ipos1.le.itot) then
       do 110 i=ipos1,itot
       ibl=iposit(i)+length
       call rdedx(h3,nx,ibl,numscr)
       mp=mp+i-1
       st(mp)=tracep(h1,h3,num)
 110   continue
      endif
c
c --- now solve the diis equations
c
      if(nstore.le.nmin)go to 400
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h3,h3(ndim+1),ifail)
c     call fixnag(ifail,'f04atf-diiscd')
      if(ifail.ne.0.and.oprint(47))write(iwr,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 210
      do 114 i=1,ndim
  114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
c     del=0.0d0
      call vclr(h0,1,nx)
      do 150 k=1,itot
      ibl=iposit(k)
      call rdedx(h1,nx,ibl,numscr)
      call daxpy(nx,ct(k),h1,1,h0,1)
 150  continue
      if(.not.ondiis) kcount=0
      ondiis=.true.
      return
 210  ibl=iposit(ipos)
      call rdedx(h0,nx,ibl,numscr)
      call reads(h1,nx,numscr)
      call wrt3(h0,nx,iblk,numscr)
      call wrt3s(h1,nx,numscr)
      nstore=1
      itot=1
      mp=iky(ipos)+ipos
      st(1)=st(mp)
      go to 410
 310  call rdedx(h0,nx,iblkh0,num8)
 300  nstore=0
      itot=0
      mp=0
      go to 410
 400  call rdedx(h0,nx,iblkh0,num8)
 410  ondiis=.false.
      return
      end
_ENDIF
_IF(parallel)
      subroutine diiscd(h0,h1,h2,h3,iblkh0,iblkqa,ndaf,numscr,
     * diff0,lockt,idomo,ivmo)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension h0(*),h1(*),h2(*),h3(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis
      common/scftim/tdiag(4),tdiis,tmult(5)
c
c
c diisc is specifically for the closed shell  cases
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and wriiten out to ed7
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting block for the vectors .. need space for 20
c                                           .. 40 for uhf
c ndafd   the block which holds the beginnig of the extrapolated vecs
c accdi1  diis comes in only when diff is less than this
c accdi2  the diis solution is only used when the residuals
c          are less than thia
c
c num8   the stream for ed7
c num     the dimension of the arrays
c nx     the size of the arrays
c
INCLUDE(common/prints)
c ===
c ===  6-triangle variant of diis - h0 overlays h2
c ===
      dumtim=dclock()
      derror=0.0d0
      num0 = newbas0
      l3=num*num
      l22=nx+nx
      nmin=3
      iblk=ndaf
      nmax=8
      length=lensec(nx+nx)
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1) goto 300
      call rdedx(h3,l3,iblkqa,idaf)
c
c ----- sort out the vectors
c
      call tdown(h3,ilifq,h3,ilifq,num0)
c
c ----- calculate the error vector
c
c     call scopy(nx,h0,1,h0(nx+1),1)
      call wrt3(h0,nx,iblkh0,num8)
      call mult2(h3,h1,h2,num0,num0,num)
      diff0=0.0d0
      ij=0
      do 5 i=1,num0
      do 4 j=1,i
      ij=ij+1
      if(i.le.na.and.j.le.na) h1(ij)=0.0d0
      if(i.gt.na.and.j.gt.na) h1(ij)=0.0d0
      if( dabs(h1(ij)).gt.diff0) then
       diff0= dabs(h1(ij))
       idomo = j
       ivmo = i
      endif
 4    continue
 5    h1(ij)=0.0d0
      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) goto 310
c
c ------ transform back to the ao basis
c
      do 10 i=1,num
      do 10 j=1,i
      dum=h3(ilifq(i)+j)
      h3(ilifq(i)+j)=h3(ilifq(j)+i)
 10   h3(ilifq(j)+i)=dum
      call mult2(h3,h2,h1,num,num,num)
      call rdedx(h1,nx,ibl7s,num8)
      call square(h3,h1,num,num)
      call mult2(h3,h1,h2,num,num,num)
c     call scopy(nx,h0(nx+1),1,h0,1)
      call rdedx(h0,nx,iblkh0,num8)
      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.lt.1) goto 100
      do 50 i=1,ipos1
      ibl=iposit(i)
      call rdedx(h3,l22,ibl,numscr)
      mp=mp+1
      st(mp)=tracep(h1,h3(nx+1),num)
 50   continue
 100  mp=mp+1
      st(mp)=tracep(h1,h1,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
      iposit(ipos)=(ipos-1)*length + iblk
      ibl=iposit(ipos)
      call dcopy(nx,h0,1,h3,1)
      call dcopy(nx,h1,1,h3(nx+1),1)
      call wrt3(h3,l22,ibl,numscr)
      ipos1=ipos+1
      if(ipos1.gt.itot) goto 115
      do 110 i=ipos1,itot
      ibl=iposit(i)
      call rdedx(h3,l22,ibl,numscr)
      mp=mp+i-1
      st(mp)=tracep(h1,h3(nx+1),num)
 110  continue
 115  continue
c
c --- now solve the diis equations
c
      if(nstore.le.nmin)goto 400
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h3,h3(ndim+1),ifail)
      if(ifail.ne.0.and.oprint(47))write(iwr,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 210
      do 114 i=1,ndim
  114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
      del=0.0d0
      call vclr(h0,1,nx)
      do 150 k=1,itot
      ibl=iposit(k)
      call rdedx(h1,l22,ibl,numscr)
      call daxpy(nx,ct(k),h1,1,h0,1)
 150  continue
      if(.not.ondiis) kcount=0
      ondiis=.true.
      go to 9999
c
c     **** see comments above
c
 210  ibl=iposit(ipos)
      call rdedx(h0,l22,ibl,numscr)
      call wrt3(h0,l22,iblk,numscr)
      nstore=1
      itot=1
      mp=iky(ipos)+ipos
      st(1)=st(mp)
      go to 410
 310  continue
c     call scopy(nx,h0(nx+1),1,h0,1)
      call rdedx(h0,nx,iblkh0,num8)
 300  nstore=0
      itot=0
      mp=0
      go to 410
400    continue
c     call scopy(nx,h0(nx+1),1,h0,1)
 410  ondiis=.false.
9999  continue
      tdiis=tdiis+(dclock()-dumtim)
      return
      end
_ENDIF
      subroutine extrpm(dee,dampp,damp00,h0,h1,h2,h3,l1,l2,ndaf1,
     +  iterl,ncall,ityp,core)
c
c     -----  use either davidson's damping and extrapolation
c            or pople's extrapolation -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/runlab)
INCLUDE(common/extsav)
      dimension core(*),h0(*),h1(*),h2(*),h3(*)
      data dzero,done,two,four/0.0d0,1.0d0,2.0d0,4.0d0/
      data pt5/0.5d0/
      data tol,tol1,tol2/1.0d-07,1.9d0,0.99d0/
      data shrnkf,dmpmin,rsmin/100.0d0,0.01d0,0.8d0/
      data zgvb,zgrhf /'gvb','grhf'/
c
      ndaf2 = ndaf1+l2
      ndaf3 = ndaf2+l2
c
c     ----- decode convergence method -----
c
      oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (ncall .eq. 1) kcount = kcount+1
      if (iter .gt. 1) go to 140
c
c     ----- iter = 1 -----
c
      call dcopy(l2,h0,1,core(ndaf1),1)
      call dcopy(l2,h0,1,core(ndaf2),1)
      call dcopy(l2,h0,1,core(ndaf3),1)
      return
c
c     ----- current fock matrix is in -h0-
c           get previous fock matrices (h1,h2,h3) -----
c
  140 call dcopy(l2,core(ndaf1),1,h1,1)
      call dcopy(l2,core(ndaf2),1,h2,1)
      call dcopy(l2,core(ndaf3),1,h3,1)
      if (iter .gt. 2) go to 160
c
c     ----- iter = 2 -----
c
      if (odamph) go to 320
      if ( .not. oshift) go to 420
      if ((zscftp.eq.zgvb.or.zscftp.eq.zgrhf).and. odiis) go to 420
      dampp = done
      go to 320
c
c     ----- iter > 2 -----
c
  160 if (oshift) go to 220
      if (odampr) go to 200
      if ( .not. odamph) go to 420
      if ( dabs(dee) .gt. exttol) go to 180
      if (dee .gt. dzero .and. kcount .gt. iextin) go to 180
      go to 420
  180 dmptol = dmptol/shrnkf
      exttol = exttol/shrnkf
      go to 300
  200 if ( dabs(dee) .gt. dmptol) go to 300
      if (dee .gt. dzero) go to 300
      if (dampp .gt. dmpmin) go to 300
      go to 420
  220 if ( .not. oextpr) go to 260
      if ( dabs(dee) .gt. exttol) go to 240
      if (dee .gt. dzero .and. kcount .gt. iextin) go to 240
      go to 420
  240 exttol = exttol/shrnkf
      dmptol = dmptol/shrnkf
      rshift = rsmin
      iterl = 0
      if (odamph) go to 300
      go to 620
  260 if ( dabs(dee) .gt. dmptol) go to 280
      if (dee .gt. dzero) go to 280
      if (odamph .and. dampp .gt. dmpmin) go to 280
      if ( .not. oextra) go to 280
      if (rshift .ge. vshtol) go to 280
      if (iterl .eq. 0) go to 280
      go to 420
  280 if (odamph) go to 300
      if (iterl .lt. 2) go to 620
      if ( .not. oextra) go to 620
      go to 440
c
c     ----- davidson's damping -----
c
  300 if (kcount .lt. iextin) go to 320
c
      if ( .not. oshift .or. iterl .ge. 2) go to 360
 320  temp=done/(done+dampp)
      do 340 i = 1,l2
 340  h0(i) = (h0(i)+dampp*h1(i))*temp
      go to 400
c
c     ----- iter > 2 , damping -----
c
  360 cutoff = pt5-damp00
      temp=done/(done+dampp)
      do 380 i = 1,l2
380   h0(i) = (h0(i)+dampp*h1(i))*temp
      if ( .not. oextra .or. cutoff .lt. dzero) go to 400
      odampr = .true.
      oextpr = .true.
c
      go to 460
  400 if (ncall .ne. 1) go to 660
      odampr = .true.
      oextpr = .false.
      go to 660
  420 if ( .not. oextra) go to 620
c
c     ----- pople's extrapolation procedure
c
      if (ncall .ne. 1) go to 460
      odampr = .false.
      oextpr = .true.
  440 dampp = dzero
c
c     ----- skip to end if first cycle or after extrapolation -----
c
  460 call dcopy(l2,h0,1,core(ndaf1),1)
      call dcopy(l2,h1,1,core(ndaf2),1)
      call dcopy(l2,h2,1,core(ndaf3),1)
      call vsub(h2,1,h3,1,h3,1,l2)
      call vsub(h1,1,h2,1,h2,1,l2)
      call vsub(h0,1,h1,1,h1,1,l2)
      if (kcount .lt. iextin .or. iter .lt. 4) go to 680
c
c     ----- find displacement dp1,dp2,dp3 -----
c
      if (ityp .eq. 2) go to 500
      sp11 = tracep(h1,h1,l1)
      sp12 = tracep(h2,h1,l1)
      sp13 = tracep(h3,h1,l1)
      sp22 = tracep(h2,h2,l1)
      sp23 = tracep(h3,h2,l1)
      sp33 = tracep(h3,h3,l1)
      go to 520
  500 continue
      sp11 = ddot(l2,h1,1,h1,1)
      sp12 = ddot(l2,h2,1,h1,1)
      sp13 = ddot(l2,h3,1,h1,1)
      sp22 = ddot(l2,h2,1,h2,1)
      sp23 = ddot(l2,h3,1,h2,1)
      sp33 = ddot(l2,h3,1,h3,1)
  520 continue
      dp1 = dsqrt(sp11)
      dp2 = dsqrt(sp22)
c
c     ----- find cosine of angle between successive displacements -----
c
      dp3 = dsqrt(sp33)
c
c     ----- find cosine of angle between -dp(3)- and
c           plane of =dp(1)- and -dp(2)-.
c
      cosphi = sp12/(dp1*dp2)
      r = sp11*sp22-sp12*sp12
      p = (sp13*sp22-sp12*sp23)/r
      q = (sp23*sp11-sp12*sp13)/r
c
c     ----- do not extrapolate unless -4- consecutive points are
c           nearly coplanar -----
c
      cospsi = dsqrt(p*p*sp11+q*q*sp22+two*p*q*sp12)/dp3
      if (cospsi .le. tol) go to 680
c
c     ----- express -dp(1)- as x*dp(3)(projected)+y*dp(2) -----
c
      if (dampp .gt. dmpmin) go to 680
      q = -q/p
c
c     ----- test if 2*2 matrix has real eigenvalues
c           between -tol/2 and +tol/2 -----
c
      p = done/p
      pq = q*q+four*p
      if (pq .lt. dzero) go to 680
      pq =  dabs(q)+dsqrt(pq)
c
c     ----- if -4- point extrapolation is not possible,
c           try -3- point
c
      if (pq .le. tol1) go to 560
      if ( dabs(cosphi) .le. tol2) go to 680
      p = dp1/(dp2*cosphi-dp1)
      call daxpy(l2,p,h1,1,h0,1)
      go to 600
  560 ppp = p/(done-p-q)
      qqq = (p+q)/(done-p-q)
      do 580 i = 1,l2
      h0(i) = h0(i)+ppp*h2(i)+qqq*h1(i)
  580 continue
c
  600 kcount = 0
      call dcopy(l2,h0,1,core(ndaf1),1)
      go to 680
c
c     ----- no damping or extrapolation -----
c
  620 if (ncall .ne. 1) go to 640
      odampr = .false.
      oextpr = .false.
c
c     ----- save new (modified) fock matrix -----
c
  640 dampp = dzero
c
  660 call dcopy(l2,h0,1,core(ndaf1),1)
      call dcopy(l2,h1,1,core(ndaf2),1)
      call dcopy(l2,h2,1,core(ndaf3),1)
  680 return
      end
      subroutine rhfclm(q0,q,coord,atmchg,lword,nss)
c ===
c === parallel in-store version (ed3/ed7)
c ===
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c     ----- closed shell hf-scf calculation -----
c           c.c.j. roothaan, rev.mod.phys. 23, 69 (1951)
c
c     irest = 3 ..... restart scf
c
c     nprint = 5   ... mo*s + convergence data printed for each
c                      iteration.
c              6   ... energy data only printed for telex
c              0   ... normal printing after convergenge or at exit
c
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/atmol3)
      common/blkin/scr(1)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/disc)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
INCLUDE(common/runlab)
INCLUDE(common/dnfnw)
INCLUDE(common/prnprn)
INCLUDE(common/psscrf)
INCLUDE(common/xfield)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/symtry)
INCLUDE(common/fermidirac)
INCLUDE(common/harmon)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     +             iposit(20),
     +             nsti(2),ondiis,junkj
c
       integer domo, vmo
       logical logtest
       common/testerp/tester(200), idomo(200), ivmo(200),
     +               testerd(200),idomod(200), ivmod(200),
     +               testdum(2,200),idomodum(2,200),
     +               logtest(200)
c
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c           run through all code section
INCLUDE(common/nodeio)
INCLUDE(common/parcntl)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/drive_dft)
INCLUDE(common/statis)
_ENDIF
INCLUDE(common/blur)
_IF(charmm)
INCLUDE(common/chmgms)
_ENDIF
      logical use_symmetry, force_serial
      integer n_irreps
INCLUDE(common/zorac)
      character *10 charwall
      dimension q(*),q0(*),coord(3,*),atmchg(*)
c     dimension osign(2)
INCLUDE(common/drfopt)
_IF(ccpdft)
c     now allocated dynamically prior to gden_energy call
c     REAL vblur(20000)
_ENDIF

_IF(drf)
cdrf2-------------------------------------------------start drf addition
INCLUDE(../drf/comdrf/sizesrf)
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfbem)
cmw this is hdafile because dafil collides with gamess  include 'da/dafil'
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
cmw this common cannot be included because gamess uses the varnames
      common/enrghf/enhh,etothh,ehfhh
      logical odrf, oarf
cdrf2-------------------------------------------------end
_ENDIF
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data bigd/1.0d4/
            call check_feature('rhfclm')
      dft_accu = 1.0d0
      oredo = .false.
      odft  = .false.
_IF(drf)
cahv --- DRF extension ---
      odrf = .false.
      oarf = .false.
      enadd = dzero
cahv
_ENDIF
_IF(ccpdft)
      idum = CD_update_geom(c)
      odft = CD_active()
_ENDIF
c
c     variables for monitoring tester + convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmn = 0.0d0
      esmear = esmear_start
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      out = nprint .eq. 5
      outon = nprint .ne. -5
      timeit = cpulft(1)
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)
_IF(parallel)
      lscf = 20+12/nav
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
_ENDIF
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = l2+l2
      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
c
c     ----- set pointers for partitioning of core -----
c
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
      i50 = 1
      i60 = i50+l2
      i80 = i60+l2
      i90 = i80+l1
      i10 = i90+l1
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      jblks = i10+nss*l4
c
      jblkf = jblks + l2
c  extra triangle for crystal energy correction
      jblkfx = jblkf+l2
      isiz=0
      if(ocryst)isiz=l2
      jblkst= jblkfx+isiz
      jblkqa= jblkst+l2
      jblkpa= jblkqa+l3
      jblkea= jblkpa+l2
      jblkqs= jblkea+l1
      jblkh0= jblkqs+l3
      jblkh = jblkh0+l2
      jblkd = jblkh+l2
      ii    = jblkd+3*l2
      ndafdi = ii
c
c allocate extra space for scaled density matrix
c
      if (ozora) then
         iscald = igmem_alloc(l2)
         iz_1 = igmem_alloc(4*num*num)
         iz_2 = iz_1 + num*num
         iz_3 = iz_2 + num*num
         iz_4 = iz_3 + num*num
         call vclr(q0(iz_1),1,4*num*num)
      endif
c 
      do 180 loop=1,8
      iposit(loop)=ii
180   ii=ii+l2+l2
      last = ii-1
_IF(drf)
cdrf2----------------------------drf: 2-el. reaction f. contrib
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
cafc
c----  set acur for use in 2-electron rf routines drfhamx
c
        acur = acurcy
        if (iarfcal .eq. 0) then
          odrf = .true.
        else
          oarf = .true.
        endif
c
c      ovlap at i110  from da10, 12, l2
c      dipx  at i120  from da10, 53, l2
c      dipy  at i130  from da10, 54, l2
c      dipz  at i140  from da10, 55, l2
c      omega(s) at i150    da31, 50, nomega
c      omega(op) at i160   da31, 51, nomega, if itwoeps=1
c      iexp  at i170       da31, 2,  l2
c      ijbit(s) at i180    da31, 56, l2
c      ijbit(op) at i190   da31, 57, l2, if itwoeps=1 else clearit
c      fdrf      at i200   calc in drfhamc, l2
c      ddrf      at i210   l2
        i110 = last + 1       
        i120 = i110 + l2      
        i130 = i120 + l2      
        i140 = i130 + l2      
        i150 = i140 + l2         
        if (odrf) then
          i200 = i150 + nomga   
        else
c
c         relay at i150
c
          i200 = i150 + ndim*ndim
        endif
        i210 = i200 + l2      
        i170 = i210 + l2        
        i180 = i170 + l2        
        if (odrf) then
        if (itwoeps .eq. 1) then
          i160 = i180 + l2        
          i190 = i160 + nomga
        else
          i160 = i180 + 1       
          i190 = i180 + l2
        endif
        last = i190 + l2
        else
c       wt/vr at i180
c       nuclear + external source at i190
c       indx at i160
c

          i190 = i180 + nwtr*nwtc
          i160 = i190 + ndim
          i193 = i160 + ndim
          i196 = i193 + ndim
          last = i196 + ndim
        endif
      endif
cdrf2--------------------------------------------------end
_ENDIF
c
      if(last.gt.lword)then

      write(iwr,9309)last,lword
      call caserr('insufficient memory available')
      endif
      if (out) write (iwr,9308) i50,i60,i80,i90,i10,i20,i21,
     +     i30,i31,i40,last
c
c     ----- set timing variables
c
      timedf= 0.0d0
      timeh = 0.0d0
      times = 0.0d0
      timedi= 0.0d0
      timed = 0.0d0
      timet = 0.0d0
      timev = 0.0d0
      timeo = 0.0d0
      timetc= 0.0d0

      ovir = .false.
      if(oaimpac())ovir = .true.
c
c     ----- occupation numbers -----
c
_IF(drf)
cdrf2--------------------------------------------------------drf extension
      if (odrf .or. oarf) then
c
c  -----  read necessary data from da10
c
c     ----- occupation numbers -----
        call daread(idafh,ioda,q(i110),l2,12)
        call daread(idafh,ioda,q(i120),l2,53)
        call daread(idafh,ioda,q(i130),l2,54)
        call daread(idafh,ioda,q(i140),l2,55)
        call daread(idafdrf,iodadrf,q(i170),l2,2)
      endif
      if (odrf) then
c
c  -----  read necessary data from da31
c
        call daread(idafdrf,iodadrf,q(i150),nomga,50)
        call daread(idafdrf,iodadrf,q(i180),l2,56)
c
        if (itwoeps .eq. 1) then
          call daread(idafdrf,iodadrf,q(i160),nomga,51)
          call daread(idafdrf,iodadrf,q(i190),l2,57)
        else
          call clear(q(i190),l2)
        endif
c
c       iodadrf(1,72) = 0
        call clear(q(i200),l2)
      else if (oarf) then
        call daread(idafdrf,iodadrf,q(i150),ndim*ndim,43)
        call daread(idafdrf,iodadrf,q(i160),ndim,44)
        call daread(idafdrf,iodadrf,q(i190),ndim,101)
      endif
cdrf2-------------------------------------------------------end drf extension
_ENDIF
c     -o- at q(i80)
c
c
c     ----- initialize variables -----
c
      if (maxcyc .le. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
       tester(loop) = 0.0d0
       testerd(loop) = 0.0d0
       logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1
c
c     i=lensec(l2)
c     iblkh0 = ibl7la
c     iblkqq = iblkh0 + i
c     iblkh=iblkqq+lensec(l3)
c
      if(numdis.eq.0) numdis=num8
      if(irest.ne.3.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         kcount = 0
         iterv = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      diffdp = bigd
      diffdpp = bigd
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
c
      ovir = .false.
      if(oaimpac())ovir = .true.
c
      if (odamph .or. oshift) damp = done

c
      lockt=lock
c
      if(odenfu)then
      write(iwr,*)'denfun',odenfu,osymm,osc,idenfu,oprint(25)
c
c - put aside memory for dft contribution
c
         if(osc)then
            idf = last+1
            last = last + nx
         endif
         write(iwr,*)'mem chk',last, lword, lword - last
         if(last.gt.lword)call caserr('memory ran out')
c
c  --- set up dft integration grid
c
         write(iwr,*)'start set up df grid',cpulft(1)
         call scdfgr(q(last + 1),lword - last,npt)
         write(iwr,*)'end set up df grid',cpulft(1)
         write(iwr,*)'there are',npt,'points'
         iwei = last + 1
         last = last + npt
         timdum=cpulft(1)
         timegw =  timdum - timeit
         timeit = timdum
      endif
c
c     ----- nuclear energy
c     ----- psscrf dvsg
c
      if (opssc) then
         en = denuc(nat,atmchg,coord,itotbq,iwr)
      else if (ocryst) then
         en = crnuc(nat,czan,coord,ecrn)
      else
         en = enucf(nat,atmchg,coord)
      endif
c
c     l0 = number of canonical vectors kept.
c
      l0=newbas0
      n_irreps = 0
      Do i = 1, l1
         isymmo( i ) = isymao( i )
         n_irreps = Max( n_irreps, isymao( i ) )
      End Do
      use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
      if(outon) then
        write (iwr,9008) en
        if(use_symmetry) write(iwr,9009)
_IFN1(tcx)        write (iwr,9028) maxcyc,mconv,nconv,mpunch
_IF1(tcx)        write (iwr,9028) maxcyc,mconv,nconv,mpunch,nss
      endif
c     Write out the header for the SCF module that contains
c     the parameters used for the following cycles
      call blkscfhead(nconv, mconv, maxcyc,dmpcut,accdi1,accdi2)
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i10),l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,q(i20),q(i10))
      endif
_ENDIF
c
c     ----- compute canonical orthonormal vectors -----
c
c     -s - at q(i10)
c     -sv- at q(i20)
c     -se- at q(i31)
c
      call rdedx(q(i10),l2,ibl7st,num8)
      call dcopy(l2,q(i10),1,q(jblkst),1)
      call qmat_symm(q,q(i10),q(jblkqs),q(i31),q(i40),iky,l0,l1,l3,
     * l1,out,isymmo,use_symmetry)
c
c     ----- compute initial guess density matrix -----
c
c     -d- at q(i30)
c
c     ----- orthogonalise vectors
c
c load guess vectors..
      call rdedx(q(jblkqa),l3,ibl3qa,idaf)
c load overlap
      call dcopy(l2,q(jblkst),1,q(i50),1)
      if (.not.onoor) then
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(jblkqa),q(i50), num, l0)
       else
          call mult2(q(jblkqa),q(i10),q(i50),l0,l0,l1)
          call orfog(q(jblkqa),q(jblkqa),q(i10),q(i20),iky,
     *         ilifq,l0,l1,1)
       endif
_ELSE
      call mult2(q(jblkqa),q(i10),q(i50),l0,l0,l1)
      call orfog(q(jblkqa),q(jblkqa),q(i10),q(i20),iky,
     *     ilifq,l0,l1,1)
_ENDIF
      else
        write(iwr,'(2a)') '** warning start-orbitals not orthogonal **',
     1                    '** result may not be variational         **'
      end if
c
      call wrt3(q(jblkqa),l3,ibl3qa,idaf)
      call tdown(q(i30),ilifq,q(jblkqa),ilifq,l0)
      call rdedx(q(i90),l1,ibl3ea,idaf)
      call dcopy(l1,q(i90),1,q(jblkea),1)
      ndoc = na
      if (osmear) then
        call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                   iwr,oprint(60))
        emrmn = 2*emrmn
      else
        call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
      call dcopy(l2,q(i50),1,q(jblkpa),1)
      call rdedx(q(jblks),l2,ibl7s,num8)
      call rdedx(q(jblkf),l2,ibl7f,num8)
c...  add zora corrections (atom/restore)
      if (ozora) call zora(q0,q0,q(jblkf),'read')
c
_IF(ccpdft)
      if(oblur)then
c
c add contribution from gaussian blurred charges
c to one electron hamiltonian
c
         eblur = 0.0d0

         idum = gden_energy(c,q(jblkf),dum,q(i50),dum,
     +        e_blur_elec, e_blur_nuc,q0,q0,outon,ochmdbg,iwr)

         write(iwr,*)'Initial electron-blur energy',e_blur_elec
         write(iwr,*)'Nuclear-blur energy',e_blur_nuc
         en = en + e_blur_nuc

      endif
_ENDIF

c...
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      lprnt = l1
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
_IF(drf)
cdrf2-------------------------------------------------------start drf add
      if (oreact) then
c      write start-density mat to da31(idafdrf) too
c     (which is a copy of density matrix)
c
       if (igetden.eq.1) 
     + call daread(idafh,ioda,q(i50),l2,16)
       if (odrf .or. oarf) call dcopy(l2,q(i50),1,q(i210),1)
       call dawrit(idafh,ioda,q(i50),l2,16,navh)
      else
       lprnt=min(lprnt,l0)
      endif
c
cdrf2----------------------------------------------------------------end
_ELSE
      lprnt=min(lprnt,l0)
_ENDIF
c
c     ----- find the time remaining at the start -----
c
      call timrem(tlefts)
c
      timdum=cpulft(1)
      timei =  timdum - timeit
      timeit = timdum
_IF(ccpdft)
      if (CD_active()) then
         call retrieve_spare(imemspare)
         imemfree  = igmem_max_memory() - imemspare 
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q0,q0,iwr)
         else if (ks_bas.eq.KS_MO) then
            imemreq = CD_memreq_energy_mo(l1,na,0,q0,q0,iwr)
         else if (ks_bas.eq.KS_AOMO) then
            imemreq = CD_memreq_energy(l1,na,0,q0,q0,iwr)
         endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            write(iwr,600)ierror
            call caserr('Out of memory in incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
c     ----- start scf procedure -----
c
c
c     ----- construct a skeleton fock matrix -----
c
c     -d- at q(i50)
c     -h- at q(i10)
c  -exch- at q(i20)
c
 140  continue
c
_IF(parallel)
c***   ***node-MPP***
      call pg_synch(4444)
_ENDIF
c
_IF(ccpdft)
c
c Modify fock builder options
c
      idum = CD_set_2e()
_ENDIF
_IF(cray,convex,titan)
      call hstar(q(i50),q(i10),q(i20),nopk,nss)
_ELSE
      call hstar(q(i50),q(i10),q(i20),nopk)
_ENDIF
c

_IF(ccpdft)
c
c restore fock builder options
c
      idum = CD_reset_2e()
_ENDIF

_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
c***   ***node-MPP***
_ENDIF
      if (ozora) then
c
c sf compute zora contributions to one-electron integrals
c
         itert = iter
         if (ocvged) itert = 999999999
         call rdedx(q(jblkf),l2,ibl7f,num8)
         if ((mod(iter,niter_z).ne.0.and.itert.ne.999999999)
     1       .or.nat_z.ne.0) then
            call zora(q0,q(i50),q(jblkf),'calc')
         else
            call zora(q0,q(i50),q(jblkf),'force')
         end if
         if (oioraz) then
            call zora(q0,q0(iz_1),q0(iz_2),'scale')
c
c q0(iz_2) now contains -X.X-dagger
c
c transform to adapted basis and complete metric
c
            call tranp(q0(iz_2),q0(iz_1))
            call dscal(l2,-1.0d0,q0(iz_1),1)
c
            call vadd(q0(iz_1),1,q(jblkst),1,q0(iz_1),1,l2)
c at iz_1 the adapted 1+X.X-dagger being the iora metric
         endif
      endif
      if (out) then
      write (iwr,9068)
      call prtril(q(i10),l1)
      endif
      timdum=cpulft(1)
      timeh = timeh + timdum - timeit
      timeit = timdum
_IF(drf)
c
      if (oreact) then
c
      call dawrit(idafh,ioda,q(i10),nx,103,navh)
cdrf2-------------------------------------------------------------add
c     at this point, add the drf contributions
      if (odrf) then
c
c --- update drf 2_el. hamiltonian
c     fdrf    at i200
c     difdens at i210 read in from dafile in drfhamc
c
      call drfhamc(q(i200),q(i210),q(i110),q(i120),q(i130),
     1    q(i140),q(i170),q(i150),q(i180),q(i160),q(i190))
c
c --- add to fock matrix   (scffact * fdrf(i) )
c
        do 1111, i = 1, l2
          q(i10-1+i) = q(i10-1+i) + scffact*q(i200-1+i)
 1111   continue
        call dcopy(l2,q(i50),1,q(i210),1)
      else if (oarf) then
        call arfone(q(i200),q(i50),
     1  q(i110),q(i120),q(i130),q(i140),
     2  q(i150),q(i160),q(i180),q(i170),
     3  q(i190),q(i193),q(i196),enadd)
        call dcopy(nx,q(jblkf),1,q(i210),1)
        call vadd(q(jblkf),1,q(i200),1,q(jblkf),1,nx)
        call wrt3(q(jblkf),nx,ibl3f,idaf)
      endif
c
      endif
cdrf2-------------------------------------------------------------end
_ENDIF
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at q(i10)
c     scratch area at q(i20)
c
      call symh(q(i10),q(i20),iky,0,0)
      if (out) then
         write (iwr,9088)
         call prtril(q(i10),l1)
      endif
cdft+
      if(osc)then
c
c-compute fock matrix contributions
c-set orbital occupancies
c
         call scdfoc(q(i80),nocmx)
         call scdft(q(i30),q(idf),q(iwei),edft,1)
         if(out)then
            write(iwr,*)'density functional fock contribution'
            call prtril(q(idf),l1)
         endif
         timdum = cpulft(1)
         timedf = timedf + timdum - timeit
         timeit = timdum
      endif
c
c     ----- read in core hamiltonian matrix
c           and calculate hf energy -----
c
c     -h0- at q(i20)
c     - h- at q(i10)
c     - e- at q(i30)
c

      call vadd(q(i10),1,q(jblkf),1,q(i10),1,nx)

c
c  save previous total 
c
      ehf0 = ehf
c
c  ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i50),q(jblkf),l1)

_IF(drf)
      if (oarf)  call dcopy(nx,q(i210),1,q(jblkf),1)
_ENDIF

_IF(ccpdft)
      if(oblur)then
c
c add contribution from gaussian blurred charges
c to one electron hamiltonian
c allocate vblur array (originally hard-wired?)
c
         nblur = 20000
         iblur = igmem_alloc(nblur)
         idum = gden_energy(c,q(iblur),dum,q(i50),dum,
     +        e_blur_elec, e_blur_nuc,q0,q0,outon,ochmdbg,iwr)
         call gmem_free(iblur)
         write(iwr,*)'Electron-blur energy',e_blur_elec
      endif

      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(CD_HF_exchange() .or. CD_HF_coulomb())then
c
c Coulomb or exchange operator is in q(i10) augmented by H_1, 
c compute energy using HF expression (maybe without exchange)
c
            ehf2 = tracep(q(i50),q(i10),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1
         endif
c
c Update Kohn-Sham matrix and compute fitted/integrated 
c energy terms
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(q0(isma),1,l2)
_ENDIF
         if (ks_bas.eq.KS_AO) then
            idum = CD_energy_ao(c,q(i10),dum,q(i50),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_MO) then
            idum = CD_energy_mo(l1,na,0,c,q(i10),dum,q(i30),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_AOMO) then
            idum = CD_energy(l1,na,0,c,q(i10),dum,q(i30),dum,q(i50),dum,
     +           edft,q0,q0,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q0(isma)
_ENDIF
     +           )
         else
            call caserr("rhfclm: illegal ks_bas")
         endif
_IF(debug_S)
         call compare_S(q0(isma),q0(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF
         call symm_op(q(i10))

         ehf = etmp+edft

         if (out) then
          if(opg_root()) then
            call CD_print_dftresults(.true.,.false.,iwr)
            write(iwr,9407)ehf1, ehf2
          endif
         endif
         call timana(31)
         call cpuwal(begin,ebegin)

      else
c
c  Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
         ehf2 = tracep(q(i50),q(i10),l1)

c         
         ehf = (ehf1+ehf2)*pt5

      endif
_ELSE
      ehf2 = tracep(q(i50),q(i10),l1)
c
c probably time to delete this code?!
c
      if(osc)call vadd(q(i10),1,q(idf),1,q(i10),1,nx)
      ehf = (ehf1+ehf2)*pt5
_ENDIF
      ehf = ehf + 0.5d0*emrmn


c for xtal field calculations, the core hamiltonian
c includes the crystal potential, but for the energy 
c all terms modelling molecule-molecule interactions 
c must be divided by two 

      if(ocryst)then
         m53=53   !!!!!!!!
         write(iwr,*)'secget',isecfx,m53
         call secget(isecfx,m53,iblk)
         call rdedx(q(jblkfx),l2,iblk,idaf)
         ehf3 = tracep(q(i50),q(jblkfx),l1)
         write(iwr,*)'ehf3',ehf3
         esave = ehf
         ehf = (ehf1+ehf2-ehf3)*pt5
         ecre = esave - ehf
         write(iwr,*)'new ehf, ecre',ehf,ecre
      endif

cdft+
      if(osc)then
         ehf = ehf+edft
         if(out)write(iwr,*)' dft correction to scf energy ',edft
      endif

c
c compute kinetic energy and virial
c
      if(ovir)call rhfvir(q(i50),ehf,en,ek,vir,num8,l1,l2)
c
c      ----- save fock matrix
c
c     -h- at q(i10)
c
      call dcopy(l2,q(i10),1,q(jblkh),1)
      if (out) then
         write (iwr,9048)
         call prtril(q(i10),l1)
      endif
c
      iter = iter+1
      etot0 = etot
      etot = ehf+en
_IF(drf)
      etot = etot+enadd
_ENDIF
      dep = de
      de = ehf-ehf0

      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
c
c     ----- damp and extrapolate hamiltonian matrix -----
c
c     -h - at q(i10)     hamiltonian matrix (n th iteration)
c     -ho- at q(i20)     old h matrix (n-1 th iteration)
c     -ha- at q(i30)     ancient h matrix (n-2 th iteration)
c     -hp- at q(i40)     prehistoric h matrix (n-3 th iteration)
c
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
c
      if (damp .lt. dmpcut) damp = dmpcut
      call extrpm(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,jblkd,
     +iterv,1,1,q)

      timdum=cpulft(1)
      times = times + timdum - timeit
      timeit = timdum
      diffpp=diffp
      diffp=diff
c
      if(odiis) then
_IF(ga)
      if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
         call diis_ga(q(i10),q(i30),q(i50),ibl3qa,diffd,domo,vmo)
      else
        call diiscm(q,q(i10),q(i30),q(i10),q(i50),jblkh0,jblkqa,
     *   jblks,ndafdi,diffd,domo,vmo)
      endif
_ELSE
        call diiscm(q,q(i10),q(i30),q(i10),q(i50),jblkh0,jblkqa,
     *   jblks,ndafdi,diffd,domo,vmo)
_ENDIF
        if (iter.le.maxcycp) then
          idomod(iter) = domo
          ivmod(iter) =  vmo
          testerd(iter) = diffd
        endif
        timdum=cpulft(1)
        timedi = timedi + timdum - timeit
        timeit = timdum
      endif
c
c ------ take precautions if tester is increasing again
c
      if (iter.eq.1 .or. .not.odiis) goto 223
      if (ondiis) then
        if (oprind) then
         write(iwr,9824)diffd,diffdp
         write(iwr,9823)diffd, acurcy*diisf
         write(iwr,*) 'de = ', de
        endif
        if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.acurcy*diisf
     1  .or.(otestd.and.de.lt.1.0d-10)) then
*        if (de.lt.1.0d-5) then
          if(iter.le.maxcycp) logtest(iter) = .true.
          diffdpp = diffdp
          diffdp = diffd
          go to 223
*        else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
*         odynamic = .true.
*         accdi1 = 0.0d0
*         itdiis = 0
*         if(oprind) write(iwr,9828) iter
*        endif
        else
          if(oprind) write(iwr,9829) iter
        endif
        nsti(1)=0
        ondiis=.false.
        call dcopy(l2,q(jblkh),1,q(i10),1)
        diffdp = bigd
        diffdpp = bigd
      else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c
c       if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c        iwild = iwild + 1
c        osign(iwild) = de.lt.0.0d0
c        if (iwild.eq.2.and..not.oreset.and.
c    +       .not.(osign(1).and.osign(2)) ) then
c         gapa1 = 5.0d0
c         gapa2 = gapa1
c         oreset = .true.
c         iwild = 0
c         write(iwr,9831) gapa1
c         maxcyc = maxcyc + maxcyc
c         maxcycp = min(200,maxcyc)
c        endif
c       endif
        if (odynamic) then
         if(diff.le.diffp.or.diffp.le.diffpp) then
          itdiis = itdiis + 1
          if (itdiis.eq.3) then
            accdi1 = diff
            if (oprind) write(iwr,9826) accdi1
          endif
         else
          if (oprind) write(iwr,9827)
          itdiis = 0
         endif
        endif
      endif
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     - q- at q(i30) orthonormalizing transformation
c     - h- at q(i10)
c     -h'- at q(i50) transformed h matrix
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
 223  continue
      if (symm_diag) then
         call characterize_mo( l1, l0, q( jblkqa ), isymao, 
     +                         n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      if (ozora .and. oioraz) then
         if (iter .eq. 1) print *,'IORA option invoked'
         jblkq = jblkqs
         if (oshift) jblkq = jblkqa
c
c re-orthonormalize transformation for new metric
c
         call mult2(q(jblkq),q0(iz_2),q0(iz_1),l0,l0,l1)   
         call orfog(q(jblkq),q(jblkq),q0(iz_2),q0(iz_1),iky,
     *              ilifq,l0,l1,1)
c        
         if (op2zora) then
             print *,'original ortho trafo'
             call prsq(q(jblkq),l0,l0,l1)
         end if
c new metric constructed 
      end if
c
c     transform the Fock matrix to MO-basis
c
      call tdown(q(i30),ilifq,q(jblkqa),ilifq,l0)
c
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diff=0.0d0
      ii=ndoc+1
      if(ii.le.l0)then
        do 250 i=ii,l0
          ij=iky(i)+i50
          loop=idamax(min(i-1,nocmx),q(ij),1)
          if(loop.gt.0) then
           dum = dabs(q(ij+loop-1))
           if (dum.gt.diff) then
             diff = dum
             if(iter.le.maxcycp) then
              tester(iter) = diff
              idomo(iter) = loop
              ivmo(iter) = i
             endif
           endif
          endif
 250    continue
      endif
c
      dmplim = dmax1(dmpcut,2.0d0)
      if(ondiis.and..not.osmear)diff=diffd
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy) .and. (iter.gt.1)
      if (ocvged)  then
       oredo = .false.
       if (ozora) then
c...         see if our coulomb matrix is up to date
          if (mod(iter-1,niter_z).ne.0.and.itert.ne.999999999.and.
     1        .not.(is_z .ne. 0.or.icoul_z .eq. 3.or.nat_z.ne.0)) then
             if (opzora) write(iwr,*) ' ** restart coulomb calc **'
c...          reset diis ...
             nsti(1) = 0
             oredo = .true.
             go to 666
          else if (oscalz) then
c sf compute scaling factors
             call zora(q0,q(jblkpa),q(i80),'scale')
c sf perform scaling
             id_zora1 = igmem_alloc(l2)
             call scale_z(q0,q(i50),q0(id_zora1),q(i30),q(i80),
     1                    ehf,etot,ne,num,0)
             call gmem_free(id_zora1)
             goto 666
          end if
       endif
c
       timdum=cpulft(1)
       timet = timet + timdum - timeit
       timeit = timdum
       go to 321
      endif
c
c            shift the diagonal of the transformed
c           h matrix by rshift for the virtual part.
c
666   continue
      rshift=0.0d0
      if(.not.ondiis.or.olevd)
     *  call shiftq(q(i50),nocmx,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)
      timdum=cpulft(1)
      timet = timet + timdum - timeit
      timeit = timdum
c
c     ----- diagonalize new hamiltonian matrix -----
c
c     -h- at q(i50)
c     -v- at q(i30)
c     -d- at q(i90)
c
      m2=2
      diaacc = diff*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11

*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
      timdum=cpulft(1)
      timed = timed + timdum - timeit
      timeit = timdum
c
c     ----- back-transform the eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i41)
c     -q- at q(i10)
c     scratch area at q(i50)
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(jblkqa),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(jblkqa),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(jblkqa),q(i50),l0,l1,l1)
_ENDIF
      timdum=cpulft(1)
      timev = timev + timdum - timeit
      timeit = timdum
c
c     ----- if vshift is true, reorthogonalize the vectors and
c           thereby the transformation matrix every igs th iteration.
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
      ig = mod(iter,igs)
      if (ig .eq. 0 .and. .not. oioraz) then
         call dcopy(l2,q(jblkst),1,q(i10),1)
c
_IF(ga)
         if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
            call porth(q(i30),q(i10), num, l0)
         else
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
         endif
_ELSE
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
         call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
c
         timdum=cpulft(1)
         timeo = timeo + timdum - timeit
         timeit = timdum
      endif
      call dcopy(l3,q(i30),1,q(jblkqa),1)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
c     ----- form density matrix -----
c
c     -v- at q(i30)
c     -o- at q(i80)
c     -d- at q(i50)
c
      ndoc = na
      ig=nocmx
      if (osmear) then 
         esmear = max(esmear_final,
     &                min(esmear,esmear_scale*(diff-acurcy)))
         call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                    iwr,oprint(60))
         emrmn = 2*emrmn
      endif
      if(ig.lt.l0)
_IFN1(civu)     * call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civu)     * call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      if (.not.osmear) then
         call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      if (out) then
         write (iwr,9168)
         call prtril(q(i50),l1)
      endif
c
c     print vectors
c
      if (out) then
         write (iwr,9148)
         call prev(q(i30),q(i90),lprnt,l1,l1)
      endif
c
c     compute RMS convergence on density
c
         if (iter.ne.1) then
          dsq = cvgdens(q(i50),q(jblkpa),l2)
          dsq = dsq / dfloat(l2)
          diffdens = dsqrt(dsq)
          if (oprind) write(iwr,9876) diffdens
         endif
c
c     ----- save mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i50)
c     -d- at q(i90)
c
         call dcopy(l2,q(i50),1,q(jblkpa),1)
         call dcopy(l1,q(i90),1,q(jblkea),1)
c
 321  if(numdis.ne.num8)then
         call wrt3(en,nw1d,ibldis,numdis)
         call wrt3s(st,nw2d,numdis)
      endif
      tim1=cpulft(1)
_IF(drf)
cdrf2
      if(odrf) then
         diff2=cvgden(q(i50),q(i210),l2)
      endif
cdrf2
_ENDIF
c
c     monitor static tester with diis on - at the moment remove diis
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if(oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif

c
c     ----- check and print convergence behavior -----
c
      otrigp = maxcyc-iter.lt.10
      if(outon.or.otrigp) write (iwr,9188) iter,kcount,etot,
     +                ehf,de,diff,rshift,damp,derror,yavr
c     Write on cycle energy to xml/punchfile
      call blkscfc(iter,etot,ehf,de,diff)
_IF(parallel)
_IF(ga)
c
c     must update vectors on dumpfile for diis_ga
c
      call wrt3(q(jblkqa),l3,ibl3qa,idaf)
c
_ENDIF
      else
          iter = iter + 1
      endif
c 
      if(ipiomode .eq. IO_NZ_S)then
         call pg_brdcst(7125,maxcyc,lscf*8,0)
         call pg_brdcst(7123,q(i50),l2*8,0)
         call pg_synch(4444)
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif
_ENDIF
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy) .and. (iter.gt.1)
c
c     now flag covergence if STATIC tester has been encountered
c
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
c
c  update files in core or GAs
c
      if (ioupd.gt.1) call ioupdate
      timdum=cpulft(1)
      timetc = timetc + timdum - timeit
      timeit = timdum
cgdf  
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if
c
      if (ocvged.and..not.oredo) go to 400
c
c     ----- exit in case of time limit -----
c     ----- punch out restart data -----
c
      call timit(0)
      call timrem(tlefti)
      if (tlefti .gt. skale*(tlefts-tlefti)/iter) go to 380
      irest=3
      write(iwr,9367)
      call texit(0,irest)
      go to 400
  380 if (iter .lt. maxcyc) go to 140
       if(omaxcyc) then
         write (iwr,9289)
         ocvged = .true.
       else
         write (iwr,9288)
c        write out error description to xml/punchfile
         call blkerror('excessive number of SCF iterations',0)
         etot = dzero
         ehf0 = ehf
         irest=3
         ehf = -en
       endif
  400 continue
c
c     ----- print actual value of energy -----
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
      lock = lockt
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr(
     +         'Memory failure in rhfclm:CD_jfit_clean2')
         endif
      endif
_ENDIF
      if (ocvged.and.outon) write (iwr,9208)
cdft
      if(osc)then
         write(iwr,9366)'HF elect. energy ',ehf - edft
         write(iwr,9366)'DF corrln energy ',edft
         write(iwr,9366)'total elec energy',ehf
      endif
c
_IF(drf)
cdrf----------
      if (field.ne.' ') then
       ehfhh=ehf
       enhh=en
       etothh=etot
      endif
cdrf----------
_ENDIF
      write(iwr,9373) iter,tim1,charwall()
_IF(ccpdft)
      if (odft) then
         call CD_get_dftresults(npts,aelec,belec)
         write(iwr,9374) npts,aelec,(aelec-ne)/ne,edft
      endif
_ENDIF
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if (omodel) then
         write (iwr,9371) ehf,en,etot,diffdens
      else if (ocryst) then
         write (iwr,9372) ecre,ecrn,ecre+ecrn,ehf,en,etot,
     +                    diffdens
      else if (ovir)then
         write (iwr,9365) ek,ehf,en,etot,vir,diffdens
      else
         write (iwr,9368) ehf,en,etot,diffdens
      endif
      if(yavr.ne.yblnk)write(iwr,9369)
c     write out energies to xml file      
      call blkscfe(ehf,en,etot)
_IF(drf)
       if(.not.oreact) then
_ENDIF
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       write (iwr,9199)
       idum = min(iter,200)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       else
       write(iwr,9198) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       endif
       enddo
       write (iwr,9197)
       endif
_IF(drf)
       endif
_ENDIF
c
      if(outon.and.odebug(31))  then
cdft  
         total = timeh + timed + timedi + timet +
     +           timeo + times + timetc + timei + timedf
cdft
         write(iwr,9178) timeh,timed,timedi,timet,timeo,times,
     *                   timetc,timei,timedf,total
cdft+ 
         if(odenfu)then
            write (iwr,9177)' density functional grid weights     ',
     +         timegw
            if(osc)then
               write (iwr,9177)' density functional fock formation   ',
     +         timedf
            else
               write (iwr,9177)' density functional energy calcn.    ',
     +         timedf
            endif
         endif
      endif
_IF(parallel)
c
c *** now send in-core q,p and e to other nodes
c
      if(ipiomode .eq. IO_NZ_S)then
         call pg_brdcst(7124,q(jblkqa),l3*8,0)
         call pg_brdcst(7126,q(jblkpa),l2*8,0)
         call pg_brdcst(7127,q(jblkea),l1*8,0)
      endif
c
_ENDIF
      if(numdis.eq.num8)numdis=0
      ndaf = mouta
cDEBUG
c     canonicalise orbital coefficients
c     write(*,*)'*** RHFCLM 2',num,num
c     call fixorb(q(jblkqa),num,num,num)
c     call tdown(q(i30),ilifq,q(jblkqa),ilifq,l0)
c     call llvmo(q(jblkea),q(i80),na,nocmx,l1)
c     call dscal(l1,two,q(i80),1)
c     call dmtx(q(jblkpa),q(i30),q(i80),iky,nocmx,l1,l1)
cDEBUG
      call scfsav(q(jblkqa),q(jblkpa),q(jblkea),q(i80),ndaf,l1,
     *l2,ibl3pa,ibl3ea)
_IF(drf)
cdrf
      if (oreact) then
       if ((igetden.ne.1) .and. (field .ne. ' ')
     +    .and. ocvged ) then
         call dawrit(idafh,ioda,q(jblkpa),l2,16,navh)
       endif
      endif
cdrf
_ENDIF
_IF(nbo)
c
c save a copy of the fock matrix on the dumpfile
c (temporarily use section 242)
c
      len42 = lensec(l2)
      m = 42
      call secput(isect(42),m,len42,iblk42)
      call wrt3(q(jblkh),l2,iblk42,idaf)
      lds(isect(42)) = l2
_ENDIF
      if (nprint .eq. 5 .or. nprint .eq. -5) go to 440
      call dcopy(l3,q(jblkqa),1,q(i10),1)
      call dcopy(l2,q(jblkpa),1,q(i30),1)
      call dcopy(l1,q(jblkea),1,q(i21),1)
      if (oprint(46)) then
       write (iwr,7778)
       call prsq(q(i10),lprnt,l1,l1)
      end if
      call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
      call tdown(q(i10),ilifq,q(i10),ilifq,l0)
      if(otsym) 
     +  call symass(q(i10),q(i21),q(i80),q0)
      if(oprint(25))go to 440
      write (iwr,9148)
      call prev(q(i10),q(i21),lprnt,l1,l1)
c
c-    evaluate the correlation energy using density functional
c-    (non self-consistent)
c
      if(odenfu.and..not.osc) then
c-new version of code
         call scdfoc(q(i80),nocmx)
         call scdft(q(i10),q(idf),q(iwei),edft,0)
c
c-original code - as a check only
c
         call  becke(q(i10),1)
c
      endif
c
      if(outon)go to 440
      write (iwr,9168)
      call prtril(q(i30),l1)
  440 continue 
      if(oprint(49))then
      write(iwr,9169)
      ij=jblkh
      do 452 i=1,l1
      scr(i)=q(ij)
 452  ij=ij+i+1
      write(iwr,9170)(scr(i),i=1,l1)
      endif
_IF(ga)
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.false.)
      endif
_ENDIF
      if (mpunch .ne. 0) then
c
c     ----- punch the occupied orbitals -----
c
c     -v- at q(i10)
c
      call pusql(q(jblkqa),na,l1,l1)
      endif
      if (ocvged) irest = 0
      if(outon)then
      cpu=cpulft(1)
      write(iwr,7777)cpu,charwall()
      endif
      if (ozora) then
         call gmem_free(iz_1)
         call gmem_free(iscald)
      end if
      accdi1 = accdin
      call copq_again(mouta)
c
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapa2 = gapa1
c      maxcyc = maxcyc / 2
c     endif
c
      return
c
c9831 format(/15x,'**** increase level shifter to ',f8.2,' ****'/)
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9827 format(1x,'rising diff - de = ', f20.10)
 9825 format(/' *** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffd, diffdp = ',2f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9198 format(10x,i5,f11.7,' (',2i4,') (*)',6x,f11.7,' (',2i4,')' )
 9298 format(10x,i5,f11.7,' (',2i4,')',10x,f11.7,' (',2i4,') (*)' )
 9199 format(/10x,64('+')/
     + 10x,'CONVERGENCE / TESTER ANALYSIS',12x,'From DIIS'/
     + 10x,64('+')/
     + 10x, 'Iter.', '   Tester   (domo/vmo)',
     + 10x,          '   Tester   (domo/vmo)'/)
 9197 format(10x,64('+'))
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 9008 format(/1x,'----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9028 format(/15x,'convergence data'/15x,16('=')///
     *1x,'minimise dump and scratchfile i/o'//
     +     1x,'maximum number of iterations = ',i6,/,
     +     1x,'method of convergence        = ',i6,/,
     +     1x,'convergence criterion        =1.0e-',i2,/,
_IFN1(tcx)     +     1x,'punch out option             = ',i6//1x,103('=')/
_IF1(tcx)     +     1x,'punch out option             = ',i6/
_IF1(tcx)     +     1x,'vectorisation factor in h-build',i6//1x,103('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     * 18x,'energy',10x,'energy',35x,'shift'/1x,103('='))
 9188 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,a3)
 9048 format(20x,11('*'),/20x,'fock matrix'/20x,11('*'))
 9068 format(/20x,20('*')/20x,'skeleton fock matrix'/ 20x,20('*'))
 9088 format(/20x,23('*')/20x,'symmetrized fock matrix'/20x,23(
     +     '-'))
 9148 format(//1x,100('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9168 format(/20x,14('*')/20x,'density matrix'/20x,14('*'))
 9169 format(/40x,32('=')/
     *40x,'diagonal elements of fock matrix'/40x,32('=')/)
 9170 format(/5x,7f15.7)
 9208 format(/10x,16('-')/10x,'energy converged'/,10x,16('-'))
 9288 format(/10x,30('-')/10x,'excessive number of iterations'/
     +     10x,30('-'))
 9289 format(/10x,30('-')/
     +        10x,'excessive number of iterations'/
     +        10x,'but flag SCF convergence'/
     +        10x,30('-'))
 9308 format(1x,'core assignment'/' i50, i60, i80,',
     +     ' i90, i10, i20, i21, i30, i31, i40  = '/ 10i8/
     +     1x,'last = ', i8)
 9348 format(/40x,32('*') /40x,'closed-shell rhf scf calculation'
     +       /40x,32('*'))
 9365 format(
     +10x,'kinetic energy             ',f20.10/ 
     +10x,'electronic energy          ',f20.10/ 
     +10x,'nuclear energy             ',f20.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy               ',f20.10,8x,f10.7/
     +10x,'convergence on density     ',f20.10)
 9368 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9371 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total qm energy            ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9372 format(
     +10x,'electron-crystal energy    ',f20.10/
     +10x,'nuclear-crystal energy     ',f20.10/
     +10x,'total-crystal energy       ',f20.10/
     +10x,'total electronic energy    ',f20.10/
     +10x,'total nuclear energy       ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9373 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     +,f12.2,' seconds',a10,' wall',/10x,14('-')//)
 9374 format(
     +10x,'number of quadrature points',i9/
     +10x,'integrated electron count  ',f20.10,5x,
     +    'relative error ',e10.2/
     +10x,'XC energy                  ',f20.10)
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9369 format(//20x,90('*')/
     *20x,'*',88x,'*'/20x,'*',88x,'*'/
     *20x,'*',31x,'warning  state is averaged',31x,'*'/
     =20x,'*',88x,'*'/20x,90('*')//)
9309  format(//
     *1x,'words available ',i8/
     *1x,'words requested ',i8//)
9177  format(10x,a36,f10.2)
9178  format(//
     *10x,'********************************************'/
     *10x,'Scf timing statistics              (seconds)'/
     *10x,'********************************************'/
     *10x,'fock formation                      ',f8.2/
     *10x,'fock diagonalisation                ',f8.2/
     *10x,'diis                                ',f8.2/
     *10x,'fock transformation                 ',f8.2/
     *10x,'orthogonalisation                   ',f8.2/
     *10x,'symmetrisation, damping, extrpn.    ',f8.2/
     *10x,'density matrix, etc.                ',f8.2/
     *10x,'initialisation                      ',f8.2/
     *10x,'dft                                 ',f8.2/
     *10x,'********************************************'/
     *10x,'Total                               ',f8.2/
     *10x,'********************************************'/)
7777  format(//
     *' end of closed shell scf at ',f12.2,
     *' seconds',a10,' wall'//1x,104('-')//)
7778  format(//10x,44('=')/10x,
     +'molecular orbitals -- symmetry adapted basis'/10x,44('=')//)
9366  format(20x,a17,f16.10)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'=============================='/)
      end
      subroutine diiscm(q,h0,h1,h2,h3,jblkh0,jblkqa,jblks,
     * ndaf,diff0,idomo,ivmo)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),h0(*),h1(*),h2(*),h3(*)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),iposit(20),
     +             nstore,mp,ondiis
INCLUDE(common/scfopt)
_IF(parallel)
      common/scftim/tdiag(4),tdiis,tmult(5)
_ENDIF
c
c
c  diiscm is specifically for the closed shell  cases
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written out to memory
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored in memory, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting location for the vectors .. need space for 20
c accdi1  diis comes in only when diff is less than this
c accdi2  the diis solution is only used when the residuals
c          are less than thia
c
c num    the dimension of the arrays
c nx     the size of the arrays
c
INCLUDE(common/harmon)
INCLUDE(common/prints)
c ===
c ===  6-triangle variant of diis - h0 overlays h2
c ===
_IF(parallel)
      dumtim=dclock()
_ENDIF
      derror=0.0d0
      num0 = newbas0
      nmin=3
      iblk=ndaf
      nmax=8
      len2=nx+nx
c
c --- should we begin storing diis info ???
c
      al2=nx
      idomo = 0
      ivmo = 0
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1) goto 300
      call dcopy(nx,h0,1,q(jblkh0),1)
c
c ----- sort out the vectors
c
      call tdown(h3,ilifq,q(jblkqa),ilifq,num0)
c
c ----- calculate the error vector
c
      call mult2(h3,h1,h2,num0,num0,num)
      diff0=0.0d0
      ij=0
      do 5 i=1,num0
      do 4 j=1,i
      ij=ij+1
      if(i.le.na.and.j.le.na) h1(ij)=0.0d0
      if(i.gt.na.and.j.gt.na) h1(ij)=0.0d0
      if( dabs(h1(ij)).gt.diff0) then
        diff0= dabs(h1(ij))
        idomo = j
        ivmo =  i
      endif
 4    continue
 5    h1(ij)=0.0d0
c
      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) goto 310
c
c ------ transform back to the ao basis
c
      do 10 i=1,num
      do 10 j=1,i
      dum=h3(ilifq(i)+j)
      h3(ilifq(i)+j)=h3(ilifq(j)+i)
 10   h3(ilifq(j)+i)=dum
      call mult2(h3,h2,h1,num,num,num)
      call square(h3,q(jblks),num,num)
      call mult2(h3,h1,h2,num,num,num)
      call dcopy(nx,q(jblkh0),1,h0,1)
      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.lt.1) goto 100
      do 50 i=1,ipos1
      ibl=iposit(i)+nx
      mp=mp+1
      st(mp)=tracep(h1,q(ibl),num)
 50   continue
 100  mp=mp+1
c
      st(mp)=tracep(h1,h1,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
      call dcopy(nx,h0,1,q(iposit(ipos)),1)
      call dcopy(nx,h1,1,q(iposit(ipos)+nx),1)
      ipos1=ipos+1
      if(ipos1.gt.itot) goto 115
      do 110 i=ipos1,itot
      ibl=iposit(i)+nx
      mp=mp+i-1
      st(mp)=tracep(h1,q(ibl),num)
 110  continue
 115  continue
c
c --- now solve the diis equations
c
      if(nstore.le.nmin)goto 400
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)

      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h3,h3(ndim+1),ifail)

c      write(6,*)'ct',(ct(i),i=1,ndim)

c     call fixnag(ifail,'f04atf-diiscm')
      if(ifail.ne.0.and.oprint(47))write(6,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 210
      do 114 i=1,ndim
  114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
c     del=0.0d0
c
      call vclr(h0,1,nx)
      do 150 k=1,itot
      ibl=iposit(k)
      call daxpy(nx,ct(k),q(ibl),1,h0,1)
 150  continue
      if(.not.ondiis) kcount=0
      ondiis=.true.
      go to 9999
 210  ibl=iposit(ipos)
      call dcopy(nx,q(ibl),1,h0,1)
      call dcopy(len2,q(ibl),1,q(iblk),1)
      nstore=1
      itot=1
      mp=iky(ipos)+ipos
      st(1)=st(mp)
      go to 410
 310  continue
      call dcopy(nx,q(jblkh0),1,h0,1)
 300  nstore=0
      itot=0
      mp=0
      go to 410
 400  call dcopy(nx,q(jblkh0),1,h0,1)
 410  ondiis=.false.
9999  continue
_IF(parallel)
      tdiis=tdiis+(dclock()-dumtim)
_ENDIF
      return
      end
      subroutine dampd(de,dep,deavg,damp,acurcy,diff,diffp,dmptlc)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      data dzero,pt25,pt5,two,four,fac /0.0d0,2.5d-01,0.5d0,
     + 2.0d0,4.0d0,1.6d+01/
      data pt2 /0.2d0/
      dampo = damp
      etest = acurcy*acurcy
      if ( dabs(de) .lt. etest .and.  dabs(dep) .lt. etest) go to 300
      if ( dabs(de) .lt. etest) go to 320
      if ( dabs(dep) .lt. etest) go to 340
      if ((diffp-diff) .lt. dzero) go to 100
      if (  dabs(de) .ge. acurcy .or. de .gt. dzero) go to 100
c
c     ----- converged -----
c
      damp = damp/fac
      go to 280
  100 continue
      if ( de .gt. dzero) go to 200
      if (dep .gt. dzero) go to 180
      if ( de .gt. dep) go to 140
c
c     ----- de < 0. , dep < 0. , de < dep -----
c
      if (  dabs(de) .lt. two*deavg) go to 120
      damp = fac*dmax1(damp,deavg)
      go to 280
  120 if (  dabs(de) .gt. pt5*deavg) go to 280
      damp = damp/fac
      go to 280
  140 continue
c
c     ----- de < 0. , dep < 0. , de > dep -----
c
      if (de .gt. pt25*dep) go to 160
      damp = (de/dep)**2*dmax1(damp,deavg)
      go to 280
  160 damp = damp/fac
      go to 280
  180 continue
c
c     ----- de < 0. , dep > 0. -----
c
      damp = four*dmax1(damp,deavg)
      if (-de .gt. deavg) damp = damp*fac
      if (-de+dep .ge. deavg) go to 280
      damp = damp/fac
      go to 280
  200 continue
      if (dep .gt. dzero) go to 220
c
c     ----- de > 0. , dep < 0. -----
c
      damp = four*dmax1(damp,deavg)
      if (de .gt. pt5*deavg) damp = damp*fac
      if (de-dep .ge. pt2*deavg) go to 280
      damp = damp/fac
      go to 280
  220 continue
c
c     ----- de > 0. , dep > 0. -----
c
      damp = four*dmax1(damp,deavg)
      if (de .lt. four*dep) go to 240
      damp = fac*dmax1(damp,deavg)
      go to 280
  240 if (de .gt. pt25*dep) go to 260
      damp = damp/fac
      go to 280
  260 damp = (de/dep)**2*dmax1(damp,deavg)
  280 continue
c
c     ----- if the density convergence worsened - make sure
c           that the damping can't decrease -----
c
      if ((diffp-diff) .lt. dzero) damp = dmax1(damp,dampo)
      go to 360
  300 continue
c
c     de < etest and dep < etest
      damp = damp/fac
      go to 360
c     de < etest  dep > etest
  320 continue
      damp = damp/fac
c     dep < etest  de > etest
  340 continue
      damp = dampo
      if (de .gt. dzero) damp = dmax1(two*damp,dmptlc)
  360 continue
      return
      end
      subroutine scfsav(q,p,e,pop,ndaf,l1,l2,iblkp,iblke)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
INCLUDE(common/sizes)
INCLUDE(common/runlab)
INCLUDE(common/harmon)
      dimension q(*),p(*),e(*),pop(*)
c
      data m1/1/
      l0 = newbas0
      call putq(zcom,ztitle,e,pop,l1,l1,l0,
     *m1,m1,q,ndaf,iblk)
      call wrt3(p,l2,iblkp,idaf)
      call wrt3(e,l1,iblke,idaf)
      return
      end
      subroutine wdisk(a1,a2,a3,ndaf,numscr,l2)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a1(*),a2(*),a3(*)
      call wrt3(a1,l2,ndaf,numscr)
      call wrt3s(a2,l2,numscr)
      call wrt3s(a3,l2,numscr)
      return
      end
      subroutine rdisk(a1,a2,a3,ndaf,numscr,l2)
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a1(*),a2(*),a3(*)
      call rdedx(a1,l2,ndaf,numscr)
      call reads(a2,l2,numscr)
      call reads(a3,l2,numscr)
      return
      end
      subroutine uhfop(sz,s2,q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/drfopt)
_IF(drf)
INCLUDE(../drf/comdrf/sizesrf)
INCLUDE(../drf/comdrf/gsspar)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfbem)
INCLUDE(../drf/comdrf/runpar)
c
      integer idafh, navh
      common/hdafile/idafh,navh,ioda(2,1000)
      common/enrghf/enhh,etothh,ehfhh
      logical odrf
cahv
_ENDIF
c
c     ----- unrestricted hf-scf calculation -----
c           j.a. pople and r.k. nesbet,
c           j.chem.phys. 22, 571 (1954)
c
c     irest = 3 ..... restart scf
c
c     nprint = 5   ... mo*s + convergence data printed for each
c                      iteration.
c              0   ... normal printing after convergenge or at exit
c
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/mapper)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/timez)
INCLUDE(common/scra7)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/atmol3)
INCLUDE(common/segm)
INCLUDE(common/harmon)
INCLUDE(common/fermidirac)
INCLUDE(common/runlab)
INCLUDE(common/zorac)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/datgue)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/diisd/st(210),cdiis(20),rdiis(19),derr,sdiis(20),ipos(20),
     + nsti(2),ondiis,nspce
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c           run through all code section
INCLUDE(common/nodeio)
INCLUDE(common/parcntl)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/drive_dft)
INCLUDE(common/statis)
_ENDIF
      Logical use_symmetry, force_serial
      Integer n_irreps
c
       logical logtest
       integer domoa, vmoa, domob, vmob
       common/testerp/testera(200), idomoa(200), ivmoa(200),
     +                testerb(200), idomob(200), ivmob(200),
     +               testerda(200),idomoda(200), ivmoda(200),
     +               testerdb(200),idomodb(200), ivmodb(200),
     +               logtest(200)
c
c     dimension osign(2)
      character*10 charwall
      dimension q(*)
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      data zscf/'scf'/
      data dzero,done,two     /0.0d0,1.0d0,2.0d0/
      data pt2,pt5,twopt2 /0.2d0,0.5d0,2.2d0/
      data dmptlc/1.0d-02/
      data igs /5/
      data bigd/1.0d4/
         call check_feature('uhfop')
_IF(drf)
cahv --- DRF extension ---
      odrf = .false.
cahv
_ENDIF
      odft = .false.
_IF(ccpdft)
      idum = CD_update_geom(c)
      odft = CD_active()
_ENDIF
      dft_accu = 1.0d0
c
c     variables for monitoring tester + convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmna = 0.0d0
      emrmnb = 0.0d0
      emrmn = 0.0d0
      esmeara = esmear_start
      esmearb = esmear_start
      nocca  = na
      noccb  = nb
      nocmxa = na
      nocmxb = nb
      diffdens = 0.0d0
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      oredo = .false.
      out = nprint .eq. 5
      outon = nprint .ne. -5
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(nprint.ne.-5)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)
      dfacd = acurcy*diisf
      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
_IF(parallel)
      lscf = 20+12/nav
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
_ENDIF
c
c     ----- set pointers for partitioning of core -----
c
c        -------- ---- -- ---- -- ---- -- -- -- -- -- --
c      i10       i20   i21 i30   i31 i40   i41 i50 i60 i61 i70 i71 i72
c        -------- -------- -------- -------- -- -- -- -- --
c           l2        l2        l2        l2      l1  l1  l1  ni  ni
c        -------------- --
c              l3        l1
c                            -------------- --
c                                  l3        l1
c                  -------------- --
c                        l3        l1
c
c     first evaluate total core requirements
      i10 = 0
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      i41 = i30+l3
      i50 = i40+l2
      i60 = i50+l2
      i70 = i60+l2
      i80 = i70+l2
      i81 = i80+l1
      i90 = i81+l1
      last = i90+l1
c
c     ----- now get core memory and determine actual pointers -----
c
      i10 = igmem_alloc(last)
c
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      i41 = i30+l3
      i50 = i40+l2
      i60 = i50+l2
      i70 = i60+l2
      i80 = i70+l2
      i81 = i80+l1
      i90 = i81+l1
      last = i90+l1
      if (out) write (iwr,9188) i10,i20,i21,i30,i31,i40,i41,
     +i50, i60,i70,i80,i81,i90,last
_IF(drf)
cahv --- DRF extension ---
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
cafc
c-----  set acur for use in 2-electron rf routines drfhamx
c
        acur = acurcy
        odrf = .true. 
c
c      ovlap at i110
c      dipx  at i120
c      dipy  at i130
c      dipz  at i140
c      omega(s) at i150
c      omega(op) at i160
c      iexp  at i170
c      ijbit(s) at i180
c      ijbit(op) at i190
c      fdrfa at i200
c      fdrfb at i210
c      diffdensa at i220
c      diffdensb at i230
c
        i110 = igmem_alloc(l2)
        i120 = igmem_alloc(l2)
        i130 = igmem_alloc(l2)
        i140 = igmem_alloc(l2)
        i150 = igmem_alloc(nomga)
        if (itwoeps .eq. 1) i160 = igmem_alloc(nomga)
        i200 = igmem_alloc(l2)
        i210 = igmem_alloc(l2)
        i220 = igmem_alloc(l2)
        i230 = igmem_alloc(l2)
        i170 = igmem_alloc(l2)
        i180 = igmem_alloc(l2)
        i190 = igmem_alloc(l2)
      endif
_ENDIF
      if (ozora .and. oioraz) call caserr('no iora in UHF')
c
c     ----- occupation numbers -----
c
c     -oa- at q(i80)
c     -ob- at q(i81)
c
      do 120 i = 1,l1
         popa=dzero
         popb=dzero
         if(i.le.na)popa=done
         q(i-1+i80)=popa
         if(i.le.nb)popb=done
         q(i-1+i81) = popb
  120 continue
c
c     ----- initialize variables -----
c
_IF(drf)
      if (odrf) then
        call daread(idafh,ioda,q(i110),l2,12)
        call daread(idafh,ioda,q(i120),l2,53)
        call daread(idafh,ioda,q(i130),l2,54)
        call daread(idafh,ioda,q(i140),l2,55)
c
        call daread(idafdrf,iodadrf,q(i170),l2,2)
        call daread(idafdrf,iodadrf,q(i150),nomga,50)
        call daread(idafdrf,iodadrf,q(i180),l2,56)
        if (itwoeps .eq. 1) then
          call daread(idafdrf,iodadrf,q(i160),nomga,51)
          call daread(idafdrf,iodadrf,q(i190),l2,57)
        else
          call clear(q(i190),l2)
        endif
c       iodadrf(1,72) = 0
        call clear(q(i200),l2)
        call clear(q(i210),l2)
      endif
_ENDIF
c
      if (maxcyc .lt. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
      testera(loop) = 0.0d0
      testerb(loop) = 0.0d0
      testerda(loop) = 0.0d0
      testerdb(loop) = 0.0d0
      logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1
      ehf = dzero
      ehf0 = dzero
      ek = dzero
      vir = dzero
      iter = 0
      iterv = 0
      kcount = 0
      damp = dzero
      damp0 = dzero
      if (dmpcut .le. dzero) dmpcut = dzero
      rshift = dzero
      rshfta = dzero
      rshftb = dzero
      diff = dzero
      diffp= dzero
      diffpp= dzero
      diffa = dzero
      diffb = dzero
      diffpa = dzero
      diffpb = dzero
      diffdp = bigd
      diffdpp = bigd
      diffdpa = bigd
      diffdpb = bigd
      de = dzero
      dep = dzero
      lockt = lock
c
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
      mxdiis = 3
c
      deavg = dzero
      len=lensec(l3)
      len2=lensec(l2)
      iblkh0a=ibl7la
      iblkh0b=iblkh0a+len2
      iblkqq=iblkh0b+len2
      iblkha=iblkqq+len
      iblkhb=iblkha+len2
      ndafa=iblkhb+len2
      ndafb=ndafa+3*len2
      ndafd=ndafb+3*len2
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (odamph .or. oshift) damp = done
c
      ovir = .false.
      if(oaimpac())ovir = .true.
c
      tim0=cpulft(1)
c
c     ----- nuclear energy
c
      en = enucf(nat,czan,c)
c
c     l0 = number of canonical orthonormal vectors kept.
c
      l0=newbas0
c
      n_irreps = 0
      Do i = 1, l1
         isymmo( i ) = isymao( i )
         n_irreps = Max( n_irreps, isymao( i ) )
      End Do
      use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
      if(nprint.ne.-5) then
       write (iwr,9008) en
       if (use_symmetry) write(iwr,9009)
       if(odebug(31)) then
         write (iwr,9028) maxcyc,mconv,nconv,npunch
        else
         write (iwr,9029) maxcyc,mconv,nconv,npunch
        endif
      endif
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i10),l2,ibl7st,num8)
         call declare_diis_storage(num,.true.)
         call init_diis_storage(num,q(i20),q(i10))
      endif
_ENDIF
c
c     ----- compute canonical orthonormal vectors -----
c
c     -s - at q(i10)
c     -sv- at q(i20)
c     -se- at q(i31)
c
      call rdedx(q(i10),l2,ibl7st,num8)
      if (omaster) then
      write(21)l2,(q(i10+i-1),i=1,l2)
      endif
      call qmat_symm(q,q(i10),q(i20),q(i31),q(i40),iky,l0,l1,l3,l1,out,
     +              isymmo, use_symmetry )
c
c     ----- compute initial guess density matrices -----
c
c     -da- in q(i40)
c     -db- in q(i50)
c
c
c ----- ensure vector orthogonality
c
      ipass=0
      nspin=na
      i=ibl3qa
      iocc=i80
      iblkp=ibl3pa
 801  call rdedx(q(i30),l3,i,idaf)
      call rdedx(q(i50),l2,ibl7st,num8)
c
_IF(ga)
      if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
         call porth(q(i30),q(i50), num, l0)
      else
         call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
         if (nspin.ne.0) then
            call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     +      ilifq,l0,l1,1)
         endif
      endif
_ELSE
      call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
      if (nspin.ne.0) then
         call orfog(q(i30),q(i30),q(i10),q(i20),iky,ilifq,l0,l1,1)
      endif
_ENDIF
      call wrt3(q(i30),l3,i,idaf)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (l1.ge.maxorb) call caserr('Out of range iky in dmtx(p)')
      if (ipass.eq.0) then
         nocca = na
         if (osmear) then
            call rdedx(q(i90),l1,ibl3ea,idaf)
            call fermi_smear(q(i90),q(iocc),emrmna,nocca,nocmxa,esmeara,
     &                       na,l0,iwr,oprint(60))
            nspin = nocmxa
         endif
      else
         noccb = nb
         if (osmear) then
            call rdedx(q(i90),l1,ibl3eb,idaf)
            call fermi_smear(q(i90),q(iocc),emrmnb,noccb,nocmxb,esmearb,
     &                       nb,l0,iwr,oprint(60))
            nspin = nocmxb
         endif
      endif
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
      else
         call dmtx(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
      endif
_ELSE
      call dmtx(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
_ENDIF
      if(ipass.eq.1)go to 802
      i=ibl3qb
      iocc=i81
      ipass=ipass+1
      if (.not.uhfatom) call wrt3(q(i10),l2,iblkp,idaf)
      iblkp=ibl3pb
      nspin=nb
      go to 801
 802  continue
      if (uhfatom) then
c...     get density matrices straight from atdens
         call rdedx(q(i40),l2,ibl3pa,idaf)
         call rdedx(q(i50),l2,ibl3pb,idaf)
         iter = iter - 1
      else
         call dcopy(l2,q(i10),1,q(i50),1)
         call wrt3(q(i10),l2,iblkp,idaf)
         call rdedx(q(i40),l2,ibl3pa,idaf)
      end if
c
c
_IF(drf)
c  save start densa and densb to da31
      if (odrf) then
c       call dawrit(idafdrf,iodadrf,q(i40),l2,70,navdrf)
        call dcopy(l2,q(i40),1,q(i220),1)
c       call dawrit(idafdrf,iodadrf,q(i50),l2,71,navdrf)
        call dcopy(l2,q(i50),1,q(i230),1)
      endif
c
_ENDIF
c
c     ----- find the time remaining at the start -----
c
      if (ozora) id1 = igmem_alloc(l2)
      call timrem(tlefts)
      lprnta= l1
      if(.not.oprint(20)) lprnta= min(na+5,l1)
      lprntb= l1
      if(.not.oprint(20)) lprntb= min(nb+5,l1)
c
_IF(ccpdft)
c
      idum = CD_set_2e()
c
      if (CD_2e())then
        odft = .true.
c
c  Switch to UKS mode
c
        idum = CD_uks()
c
c for simplicity, only implement full coulomb version
c here for the moment
c
        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
c
      ocoul = CD_HF_coulomb() .or. .not. odft
      o2e = .not.(odft.and..not.oexch.and..not.ocoul)
c
      if (CD_active()) then
         idum = CD_uks()
         call retrieve_spare(imemspare)
         imemfree  = igmem_max_memory() - imemspare
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q,q,iwr)
         else if (ks_bas.eq.KS_MO) then
            imemreq = CD_memreq_energy_mo(l1,na,nb,q,q,iwr)
         else if (ks_bas.eq.KS_AOMO) then
            imemreq = CD_memreq_energy(l1,na,nb,q,q,iwr)
         endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            write(iwr,600)ierror
            call caserr('Out of memory in incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
c
  180 continue
c
c     ----- start scf procedure -----
c     ----- construct a skeleton fock matrix -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c -scratch- at q(i30)
c     -da- at q(i40)
c     -db- at q(i50)
c
_IF(parallel)
c***   ***node-MPP***
      call pg_synch(4444)
_ENDIF
_IF(ccpdft)
c
c skip hstaru when not required for choice of DFT functional
c
      if(o2e) then
       idum = CD_set_2e()
       call hstaru(q(i40),q(i50),q(i10),q(i20),nopk)
      else
       call vclr(q(i10),1,l2)
       call vclr(q(i20),1,l2)
      endif
      idum = CD_reset_2e()
_ELSE
      call hstaru(q(i40),q(i50),q(i10),q(i20),nopk)
_ENDIF
_IF(parallel)
c***   ***node-MPP***
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
_IF(drf) 
c
c  reaction field contribution to 2-el part of hamiltonian
c
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
c
c  fdrf     at q(i200) and q(i210)
c  diffdens at q(i220) and q(i230)
        call drfhamu(q(i200),q(i220),q(i210),q(i230),q(i110),
     1               q(i120),q(i130),q(i140),
     2               q(i170),q(i150),q(i180),
     2               q(i160),q(i190))
c
c  add rfc to fock matrix
        do 6666 i = 1, l2
          q(i10+i-1) = q(i10+i-1) + scffact*q(i200+i-1)
          q(i20+i-1) = q(i20+i-1) + scffact*q(i210+i-1)
 6666   continue
c
c  copy densities for diffdens calculation otherwise lost
        call dcopy(l2,q(i40),1,q(i220),1)
        call dcopy(l2,q(i50),1,q(i230),1)
      endif
c
_ENDIF
      if (out) then
         write (iwr,9068)
         write (iwr,9108)
         call prtril(q(i10),l1)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c     scratch area at q(i30)
c
_IF(ccpdft)
c  if only 2e integrals were performed
      if(o2e)then
_ENDIF
      call symh(q(i10),q(i30),iky,0,0)
      call symh(q(i20),q(i30),iky,0,0)
      if (out) then
         write (iwr,9088)
         write (iwr,9108)
         call prtril(q(i10),l1)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
_IF(ccpdft)
      endif
_ENDIF
c
c     ----- read in core hamiltonian matrix
c           and calculate hf energy -----
c
c     -h0- at q(i30)
c     -fa- at q(i10)
c     -fb- at q(i20)
c     -da- at q(i40)
c     -db- at q(i50)
c
      call rdedx(q(i30),l2,ibl7f,num8)
c
      if (ozora) then
c
c... build total density (alpha + beta) for zora coulomb matrix
c
          call vadd(q(i40),1,q(i50),1,q(id1),1,l2)
          itert = iter
          if (ocvged) itert = 999999999
          if ((mod(iter,niter_z).ne.0.and.itert.ne.999999999)
     1        .or.nat_z.ne.0) then
             o_update = .false.
             call zora(q,q(id1),q(i30),'calc')
          else
             o_update = .true.
             call zora(q,q(id1),q(i30),'force')
          end if
      endif
c
      call vadd(q(i10),1,q(i30),1,q(i10),1,l2)
      call vadd(q(i20),1,q(i30),1,q(i20),1,l2)
c  save previous total
c
      if (.not.(iter.le.0.and.uhfatom)) ehf0 = ehf
_IF(ccpdft)
c
c  ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i40),q(i30),l1) + tracep(q(i50),q(i30),l1)
c
      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(o2e)then
c
c Coulomb operator is in q(i10) augmented by H_1,
c compute energy using HF expression without exchange
c
            ehf2 = tracep(q(i40),q(i10),l1) + tracep(q(i50),q(i20),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1
 
         endif
c
c Update Kohn-Sham matrix and compute fitted/integrated
c energy terms
c
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(q(isma),1,l2)
_ENDIF
         oignore = CD_ignore_accuracy()
         if (uhfatom.and.iter.lt.0) then
             ierror = CD_set_ignore_accuracy(.true.)
         endif
         if (ks_bas.eq.KS_AO) then
            idum = CD_energy_ao(c,q(i10),q(i20),q(i40),q(i50),
     +           edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_AOMO) then
            i200 = igmem_alloc(l3)
            i210 = igmem_alloc(l3)
            call rdedx(q(i200),l3,ibl3qa,idaf)
            call tdown(q(i200),ilifq,q(i200),ilifq,l0)
            call rdedx(q(i210),l3,ibl3qb,idaf)
            call tdown(q(i210),ilifq,q(i210),ilifq,l0)
            idum = CD_energy(l1,na,nb,c,q(i10),q(i20),q(i200),q(i210),
     +             q(i40),q(i50),edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +            ,q(isma)
_ENDIF
     +           )
            call gmem_free(i210)
            call gmem_free(i200)
         else if (ks_bas.eq.KS_MO) then
            i200 = igmem_alloc(l3)
            i210 = igmem_alloc(l3)
            call rdedx(q(i200),l3,ibl3qa,idaf)
            call tdown(q(i200),ilifq,q(i200),ilifq,l0)
            call rdedx(q(i210),l3,ibl3qb,idaf)
            call tdown(q(i210),ilifq,q(i210),ilifq,l0)
            idum = CD_energy_mo(l1,na,nb,c,q(i10),q(i20),
     +             q(i200),q(i210),edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +            ,q(isma)
_ENDIF
     +           )
            call gmem_free(i210)
            call gmem_free(i200)
         else
            call caserr("uhfop: illegal ks_bas")
         endif
         ierror = CD_set_ignore_accuracy(oignore)
_IF(debug_S)
         call compare_S(q(isma),q(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF
         call symm_op(q(i10))
         call symm_op(q(i20))
 
         ehf = etmp+edft
 
          if(out) then
            if(opg_root()) then
                 call CD_print_dftresults(.true.,.false.,iwr)
                 write(iwr,9407)ehf1, ehf2
            endif
          endif
         call timana(31)
         call cpuwal(begin,ebegin)
      else
c
c  Unrestricted Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
         ehfa = tracep(q(i40),q(i30),l1)+tracep(q(i40),q(i10),l1)
         ehfb = tracep(q(i50),q(i30),l1)+tracep(q(i50),q(i20),l1)
         ehf = (ehfa+ehfb)*pt5

      endif

_ELSE
      ehfa = tracep(q(i40),q(i30),l1)+tracep(q(i40),q(i10),l1)
      ehfb = tracep(q(i50),q(i30),l1)+tracep(q(i50),q(i20),l1)
      ehf = (ehfa+ehfb)*pt5
_ENDIF
      emrmn = emrmna + emrmnb
      ehf   = ehf + 0.50d0*emrmn
c
c compute kinetic energy and virial if required
c
      if(ovir)call uhfvir(q(i40),q(i50),ehf,en,ek,vir,num8,l1,l2)
c
c     ----- save -fa- and -fb- on daf -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c
      call wrt3(q(i10),l2,iblkha,num8)
      call wrt3(q(i20),l2,iblkhb,num8)
      if (omaster) then
         rewind(20)
         write(20)l2,(q(i10+i-1),i=1,l2)
         write(20)l2,(q(i20+i-1),i=1,l2)
      endif
c
      iter=iter+1
c
      etot0 = etot
      etot = ehf+en
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
c
c     ----- damp and extrapolate the alpha fock matrix -----
c
c     -f - at q(i10)     fock matrix  (n th matrix)
c     -fo- at q(i20)     old fock matrix (n-1 th matrix)
c     -fa- at q(i30)     ancient fock matrix (n-2 th matrix)
c     -fp- at q(i40)     prehistoric fock matrix (n-3 th matrix)
c
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
      if (damp .lt. dmpcut) damp = dmpcut
c
c     extrapolate alpha hamiltonian matrix
c
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafa,
     +            num8,iterv,1,1)
      call wrt3(q(i10),l2,iblkh0a,num8)
c
c     now extrapolate beta hamiltonian matrix and write to disk
c
      call rdedx(q(i10),l2,iblkhb,num8)
c
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafb,
     +            num8,iterv,2,1)
      call wrt3(q(i10),l2,iblkh0b,num8)
c
      diffpp = diffp
      diffp = max(diffpa,diffpb)
c
c     ----- invoke diis procedure
c
      if (odiis) then
_IF(ga)
        if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
          call rdedx(q(i10),l2,iblkh0a,num8)
          call rdedx(q(i20),l2,iblkh0b,num8)
          call diisu_ga(q(i10),q(i20),q(i40),q(i60),ibl3qa,ibl3qb,
     +                  diffda,diffdb,domoa,vmoa,domob,vmob)
        else
          call diisu(q(i10),q(i20),q(i40),q(i60),ndafd,diffda,diffdb,
     +               lockt,
     +               iblkh0a,iblkh0b,num8,domoa,vmoa,domob,vmob)
        endif
_ELSE
        call diisu(q(i10),q(i20),q(i40),q(i60),ndafd,diffda,diffdb,
     +             lockt,
     +             iblkh0a,iblkh0b,num8,domoa,vmoa,domob,vmob)
_ENDIF
        if (iter.le.maxcycp) then
          idomoda(iter) = domoa
          ivmoda(iter) =  vmoa
          idomodb(iter) = domob
          ivmodb(iter) =  vmob
          testerda(iter) = diffda
          testerdb(iter) = diffdb
        endif
      endif
c
c ----- take precautions if tester is increasing again
c
      if(iter.eq.0.and.uhfatom) then
         itdiis = 0
         nsti(1) = 0
         kcount = 0
         go to 223
      endif
      if(iter.eq.1) go to 223
      if (odiis) then
       if (ondiis) then
          diffd = max(diffda,diffdb)
          if (oprind) then
           write(iwr,9824)diffda,diffdpa,diffdb,diffdpb
           write(iwr,9823)diffd, dfacd
           write(iwr,*) 'de = ', de
          endif
          if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.dfacd
     1    .or.(otestd.and.de.lt.1.0d-10)) then
             if (de.lt.1.0d-6) then
              if(iter.le.maxcycp) logtest(iter) = .true.
              call wrt3(q(i20),l2,iblkhb,num8)
              diffdpa = diffda
              diffdpb = diffdb
              diffdpp = diffdp
              diffdp = max(diffdpa,diffdpb)
              go to 223
             else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
              odynamic = .true.
              accdi1 = 0.0d0
              itdiis = 0
              if(oprind) write(iwr,9828) iter
             endif
          else
             if(oprind) write(iwr,9829) iter
          endif
          nsti(1)=0
          ondiis=.false.
          call rdedx(q(i10),l2,iblkha,num8)
          call rdedx(q(i20),l2,iblkhb,num8)
          diffdp = bigd
          diffdpp = bigd
          diffdpa = diffdp
          diffdpb = diffdp
       else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c         if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c          iwild = iwild + 1
c          osign(iwild) = de.lt.0.0d0
c          if (iwild.eq.2.and..not.oreset.and.
c    +       .not.(osign(1).and.osign(2)) ) then
c           gapa1 = 5.0d0
c           gapb1 = 5.0d0
c           gapa2 = gapa1
c           gapb2 = gapb1
c           oreset = .true.
c           iwild = 0
c           write(iwr,9831) gapa1
c           maxcyc = maxcyc + maxcyc
c           maxcycp = min(200,maxcyc)
c          endif
c         endif
          if (odynamic) then
           if(diff.le.diffp.or.diffp.le.diffpp) then
            itdiis = itdiis + 1
            if (itdiis.eq.mxdiis) then
              accdi1 = dmax1(diffa,diffb)
              if (oprind) write(iwr,9826) accdi1
            endif
           else
            if (oprind) write(iwr,9827) de
            itdiis = 0
           endif
          endif
          call wrt3(q(i20),l2,iblkhb,num8)
       endif
      endif
c
223   if(.not.odiis)call rdedx(q(i10),l2,iblkh0a,num8)
c
c     first, work on alpha orbitals
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -f- at q(i10)
c     -f'-at q(i20)  transformed f matrix
c     -q- at q(i30)  orthonormalizing tranformation
c     scratch area at q(i41)
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
      if(oshift)then
         call rdedx(q(i30),l3,ibl3qa,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diffpa = diffa
      diffa=0.0d0
      ii=nocca+1
      if(ii.le.l0) then
       do 250 i=ii,l0
       ij=iky(i)+i50
       loop=idamax(min(i-1,nocmxa),q(ij),1)
       if(loop.gt.0) then
       dum = dabs(q(ij+loop-1))
        if (dum.gt.diffa) then
          diffa = dum
          if(iter.le.maxcycp) then
           testera(iter) = diffa
           idomoa(iter)  = loop
           ivmoa(iter)   = i
          endif
        endif
       endif
 250   continue
      endif
c
      if(ondiis.and..not.osmear)diffa=diffda
      dmplim=dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diffa.lt.acurcy).and.(iter.gt.1)
      if (maxcyc .eq. 0) go to 321
      if(ocvged)go to 321
c
c     shift the diagonal elements of the
c     transformed f matrix by rshift for the virtual part.
c
c
      rshift=0.0d0
      if((.not.ondiis.or.olevd).and..not.(uhfatom.and.iter.eq.0))
     *call shiftq(q(i50),nocmxa,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)
      rshfta=rshift
c
      m2=2
      diaacc = diffa*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c
c     ----- diagonalize new alpha hamiltonian matrix -----
c
c     -h- at q(i50)
c     -v- at q(i30)
c     -d- at q(i41)
c
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
c
c     ----- back-transform the alpha eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i90)
c     -q- at q(i10)
c     scratch area at q(i21)
c
      if(oshift)then
         call rdedx(q(i10),l3,ibl3qa,idaf)
      else
         call rdedx(q(i10),l3,ibl3qs,idaf)
      endif
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
      ig = mod(iter,igs)
      if ( ig .ne. 0) go to 300
c
c     ----- if vshift is true, reorthogonalize the vectors and
c           thereby the transformation matrix every igs th iteration.
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
      call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
      if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
         call porth(q(i30),q(i10), num, l0)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
         call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
c
  300 call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      nocca = na
      ig = nocmxa
      if (osmear) then 
        esmeara= max(esmear_final,
     &               min(esmeara,esmear_scale*(diffa-acurcy)))
        call fermi_smear(q(i90),q(i80),emrmna,nocca,nocmxa,esmeara,na,
     &                   l0,iwr,oprint(60))
      endif
      if(ig.lt.l0) then
_IFN1(civ)      call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civ)      call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      endif
c
      if (out) then
         write (iwr,9208)
         write (iwr,9108)
         call prev(q(i30),q(i90),lprnta,l1,l1)
      endif
c
c     ----- form alpha density matrix -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -o- at q(i80)
c
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
      else
         call dmtx(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
      endif
_ELSE
      call dmtx(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
_ENDIF
_IF(drf)
c  calculate diffdensa (old at i220, new at i20) save to da31
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        diffa = cvgden(q(i20),q(i220),l2)
c       call dawrit(idafdrf,iodadrf,q(i220),l2,70,navdrf)
      endif
_ENDIF
      if (out) then
         write (iwr,9228)
         write (iwr,9108)
         call prtril(q(i20),l1)
      endif
c
c     compute RMS convergence on alpha density
c
      if (iter.ne.1) then
       call rdedx(q(i30),l2,ibl3pa,idaf)
       dsqa = cvgdens(q(i20),q(i30),l2)
      endif
c
c     ----- save alpha mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -e- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i20),q(i90),q(i80),ndaf,l1,l2,
     * ibl3pa,ibl3ea)
c
c     ----- now work on beta set -----
c
 321  dsqb = dzero
      if (nb .eq. 0) go to 420
c
c     restore beta hamiltonian from disk
c
      if (odiis) then
        call rdedx(q(i10),l2,iblkhb,num8)
      else
        call rdedx(q(i10),l2,iblkh0b,num8)
      endif
      if(oshift)then
         call rdedx(q(i30),l3,ibl3qb,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
c
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diffpb = diffb
      diffb=0.0d0
      ii=noccb+1
      if(ii.le.l0) then
      do 253 i=ii,l0
        ij=iky(i)+i50
        loop=idamax(min(i-1,nocmxb),q(ij),1)
        if(loop.gt.0) then
          dum = dabs(q(ij+loop-1))
          if (dum.gt.diffb) then
            diffb = dum
            if (iter.le.maxcycp) then
              testerb(iter) = diffb
              idomob(iter) = loop
              ivmob(iter) = i
            endif
          endif
        endif
 253  continue
      endif
c
      if(ondiis.and..not.osmear)diffb=diffdb
      dmplim=dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diffb.lt.acurcy).and.(iter.gt.1)
      if (maxcyc .eq. 0) go to 322
      if(ocvged)go to 322
      rshift=0.0d0
      if((.not.ondiis.or.olevd).and..not.(uhfatom.and.iter.eq.0))
     *call shiftq(q(i50),nocmxb,0,l0,de,dep,iterv,2,gapb1,ibrk,gapb2)
      rshftb=rshift
      m2=2
      diaacc = diffb*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
c
c     ----- back-transform the beta eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i41)
c     -q- at q(i10)
c     scratch area at q(i21)
c
      if(oshift) then
        call rdedx(q(i10),l3,ibl3qb,idaf)
      else
        call rdedx(q(i10),l3,ibl3qs,idaf)
      endif
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
      ig = mod(iter,igs)
      if ( ig .ne. 0) go to 360
c
c     if vshift is true, reorthogonalize the beta vectors and
c     thereby the transformation matrix every igs th iteration
c
c     -q- at q(i10)
c     -s- at q(i30)
c     -v- at q(i30)
c
      call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
      if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
         call porth(q(i30),q(i10), num, l0)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
         call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
c
  360 call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      noccb = nb
      ig = nocmxb
      if (osmear) then 
        esmearb= max(esmear_final,
     &               min(esmearb,esmear_scale*(diffb-acurcy)))
        call fermi_smear(q(i90),q(i81),emrmnb,noccb,nocmxb,esmearb,nb,
     &                   l0,iwr,oprint(60))
      endif
      if(ig.lt.l0) then
_IFN1(civ)      call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civ)      call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      endif
      if (out) then
         write (iwr,9208)
         write (iwr,9128)
         call prev(q(i30),q(i90),lprntb,l1,l1)
      endif
c
c     ----- form beta density matrix -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -o- at q(i81)
c
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
      else
         call dmtx(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
      endif
_ELSE
      call dmtx(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
_ENDIF
_IF(drf)
c  calculate diffdensb (old at i230, new at i20) save to da31
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        diffb = cvgden(q(i20),q(i230),l2)
c       call dawrit(idafdrf,iodadrf,q(i230),l2,71,navdrf)
      endif
_ENDIF
      if (out) then
         write (iwr,9228)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
c
c     compute RMS convergence on beta density
c
      if (iter.ne.1) then
       call rdedx(q(i30),l2,ibl3pb,idaf)
       dsqb = cvgdens(q(i20),q(i30),l2)
      endif
c
  420 continue
c
      if (iter.ne.1) then
       dsq = (dsqa+dsqb) / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
       
       if (nb .eq. 0) call setz(q(i30),q(i20),q(i90),l1,l2,l3)
c
c     ----- save beta mo's + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -d- at q(i90)
c
      ndaf = moutb
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i20),q(i90),q(i81),ndaf,l1,l2,
     *ibl3pb,ibl3eb)
322   call rdedx(q(i40),l2,ibl3pa,idaf)
      call rdedx(q(i50),l2,ibl3pb,idaf)
      tim1=cpulft(1)
      delt=tim1-tim0
      tim0=tim1
      diffpp = diffp
      diffp = diff
      diff = dmax1(diffa,diffb)
c
c     monitor static tester with diis on - at the moment flag
c     convergence .... (removing diis is an alternative)
c     note that the thresholds for dum and dume below are
c     subject to revision.
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if (oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
c
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c     
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif
c
c     ----- print and check convergence behavior -----
c
      otrigp = maxcyc-iter.lt.10
      if((nprint.ne.-5).or.otrigp) then
       if(iter.gt.0) then
         if(odebug(31)) then
          write (iwr,9248) iter,kcount,etot,ehf,de,diff,rshfta,rshftb,
     +    damp,derr,delt,tim1
         else
          write (iwr,9249) iter,kcount,etot,ehf,de,diff,rshfta,rshftb,
     +    damp,derr
         endif
       endif
      endif
_IF(parallel)
      else
          iter = iter + 1
      endif
******
      if(ipiomode .eq. IO_NZ_S)then
         call pg_brdcst(7125,maxcyc,lscf*8,0)
c        brdcst alpha and beta density matrices
         call pg_brdcst(7123,q(i40),(l2+l2)*8,0)
         call pg_synch(4444)
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif
_ENDIF
      dmplim = dmax1(2.0d0,dmpcut)
      ocvged = (damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
c
c     now flag covergence if STATIC tester has been encountered
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
       if (ocvged)  then
       oredo = .false.
       if (ozora) then
c...         see if our coulomb matrix is up to date
          if (mod(iter-1,niter_z).ne.0.and.itert.ne.999999999.and.
     1        .not.(is_z .ne. 0.or.icoul_z .eq. 3.or.nat_z.ne.0)) then
             if (opzora) write(iwr,*) ' ** restart coulomb calc **'
c...          reset diis ...
             nsti(1) = 0
             oredo = .true.
             go to 666
          else if (oscalz) then
             descal = 0.0d0
c sf compute scaling factors
             iscal = igmem_alloc(l2)
c NOTE   the second parameter in call to zora is not used here (dummy)
             call zora(q,dummy,q(iscal),'scale') 
c
c....        scaling factors at iscal
c
c ....       alpha part first 
c ....       vectors at i10
c ....       orbital  energies at i21
c
             call rdedx(q(i10),l3,ibl3qa,idaf)
             call tdown(q(i10),ilifq,q(i10),ilifq,l0)
             call rdedx(q(i21),l1,ibl3ea,idaf)
c sf perform scaling
             call scale_uhf(q,q(i21),q(i10),q(iscal),
     1                    descal,ne,num,na)
             call wrt3(q(i21),l1,ibl3ea,idaf)
c
c ....       now beta part  
c ....       vectors at i10
c ....       orbital  energies at i21
c
             call rdedx(q(i10),l3,ibl3qb,idaf)
             call tdown(q(i10),ilifq,q(i10),ilifq,l0)
             call rdedx(q(i21),l1,ibl3eb,idaf)
c sf perform scaling
             call scale_uhf(q,q(i21),q(i10),q(iscal),
     1                    descal,ne,num,nb)
             call wrt3(q(i21),l1,ibl3eb,idaf)
             call gmem_free(iscal)
          end if
       endif
      endif
666   if (ocvged.and. .not. oredo) go to 500
      if (maxcyc .eq. 0) go to 500
      damp = dzero
c     ----- exit in case of time limit
c     ----- punch out restart data
c
      call timit(0)
      call timrem(tlefti)
      if (iter.gt.0) then
         if (tlefti .gt. skale*(tlefts-tlefti)/iter) go to 480
         irest=3
         write(iwr,9367)
         call texit(0,irest)
         go to 500
      endif
  480 if (iter .lt. maxcyc) go to 180
       if(omaxcyc) then
         write (iwr,9329)
         ocvged = .true.
       else
         write (iwr,9328)
c        write out error description to xml/punchfile
         call blkerror('excessive number of SCF iterations',0)
         etot = dzero
         irest=3
         ehf = -en
       endif
  500 lock=lockt
c
c
c     ----- print actual value of energy -----
c
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr(
     +         'Memory failure in uhfop:CD_jfit_clean2')
         endif
      endif
_ENDIF
_IF(drf)
c
cafc
c  the final densities should be written to da10
      if (field .ne. ' ') then
        ehfhh  = ehf
        enhh   = en
        etothh = etot
        call dawrit(idafh,ioda,q(i40),l2,16,navh)
        call dawrit(idafh,ioda,q(i50),l2,20,navh)
      endif
c
c     lock=lockt
_ENDIF
      if (ozora .and. oscalz ) then
         ehf = ehf - descal
         etot = etot -descal
      end if 
      if (ocvged.and.nprint.ne.-5) write (iwr,9268)
      write(iwr,9373) iter,tim1,charwall()
_IF(ccpdft)
      if (odft) then
         call CD_get_dftresults(npts,aelec,belec)
         aelec = aelec + belec
         write(iwr,9374) npts,aelec,(aelec-ne)/ne,edft
      endif
_ENDIF
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if(omodel) then
       write (iwr,9371) ehf,en,etot,diffdens
      else if (ovir)then
       write (iwr,9365) ek,ehf,en,etot,vir,diffdens
      else
       write (iwr,9368) ehf,en,etot,diffdens
      endif
c     write out energies to xml file      
      call blkscfe(ehf,en,etot)
c
_IF(drf)
       if(.not.oreact) then
_ENDIF
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       idum = min(iter,200)
       write(iwr,9199)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298)loop,testera(loop),idomoa(loop),ivmoa(loop),
     +                     testerb(loop),idomob(loop),ivmob(loop),
     +                     testerda(loop),idomoda(loop),ivmoda(loop),
     +                     testerdb(loop),idomodb(loop),ivmodb(loop)
       else
       write(iwr,9198)loop,testera(loop),idomoa(loop),ivmoa(loop),
     +                     testerb(loop),idomob(loop),ivmob(loop),
     +                     testerda(loop),idomoda(loop),ivmoda(loop),
     +                     testerdb(loop),idomodb(loop),ivmodb(loop)
       endif
       enddo
       write(iwr,9197)
       endif
_IF(drf)
       endif
_ENDIF
_IF(nbo)
c
c save a copy of the alpha- and beta-fock matrices on the dumpfile
c  = section isect(42) = for possible nbo analysis
c
      len42 = lensec(nx)
      len42 = len42 + len42
      m = 42
      call secput(isect(42),m,len42,iblk42)
      call rdedx(q(i10),l2,iblkha,num8)
      call rdedx(q(i20),l2,iblkhb,num8)
      call wrt3(q(i10),nx,iblk42,idaf)
      call wrt3s(q(i20),nx,idaf)
      lds(isect(42)) = nx + nx
_ENDIF
c
c     ----- check spin state -----
c
c     - s- at q(i10)
c     -da- at q(i40)
c     -db- at q(i50)
c     scratch area at q(i20)
c
      call spin(sz,s2,q(i40),q(i50),q(i10),q(i20),q(i30),iky,l2)
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
      if (out .or. nprint .eq. -5) go to 560
      write (iwr,9108)
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call rdedx(q(i30),l2,ibl3pa,idaf)
      call rdedx(q(i21),l1,ibl3ea,idaf)
      call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
      call tdown(q(i10),ilifq,q(i10),ilifq,l0)
      if(otsym) 
     +  call symass(q(i10),q(i21),q(i80),q)
      if(.not.oprint(25)) then
      write (iwr,9208)
      call prev(q(i10),q(i21),lprnta,l1,l1)
         if (out) then
            write (iwr,9228)
            call prtril(q(i30),l1)
         endif
      endif
      if (nb .ne. 0) then
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
         write (iwr,9128)
         call rdedx(q(i10),l3,ibl3qb,idaf)
         call rdedx(q(i30),l2,ibl3pb,idaf)
         call rdedx(q(i21),l1,ibl3eb,idaf)
         call analmo(q(i10),q(i21),q(i81),ilifq,l0,l1)
         call tdown(q(i10),ilifq,q(i10),ilifq,l0)
         if(otsym) 
     +     call symass(q(i10),q(i21),q(i81),q)
         if (.not.oprint(25))then
          write (iwr,9208)
          call prev(q(i10),q(i21),lprntb,l1,l1)
            if (out) then
               write (iwr,9228)
               call prtril(q(i30),l1)
            endif
         endif
      endif
560   continue
_IF(ga)
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.true.)
      endif
_ENDIF
      if (mpunch .ne. 0) then
c
c     ----- punch the occupied orbitals
c
c     -v- at q(i10)
c
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call pusql(q(i10),na,l1,l1)
         if (nb .ne. 0) then
            call rdedx(q(i10),l3,ibl3qb,idaf)
            call pusql(q(i10),nb,l1,l1)
         endif
      endif
c
      if (ocvged) irest = 0
      if(nprint.eq.-5)go to 7778
      cpu=cpulft(1)
      write(iwr,7777)cpu,charwall()
c
c     ----- reset core memory -----
c
 7778 continue           
_IF(drf)
      if (oreact) then
       if (odrf) then
        call gmem_free(i190)
        call gmem_free(i180)
        call gmem_free(i170)
        call gmem_free(i230)
        call gmem_free(i220)
        call gmem_free(i210)
        call gmem_free(i200)
        if (itwoeps .eq. 1) call gmem_free(i160)
        call gmem_free(i150)
        call gmem_free(i140)
        call gmem_free(i130)
        call gmem_free(i120)
        call gmem_free(i110)
       endif
      endif
_ENDIF
      if (ozora) then
          call gmem_free(id1)
      end if
      call gmem_free(i10)
      accdi1 = accdin
      call copq_again(mouta)
      call copq_again2(moutb)
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapb1 = 3.0d0
c      gapa2 = gapa1
c      gapb2 = gapb1
c      maxcyc = maxcyc / 2
c     endif
      return
c
c9831 format(/15x,'**** increase level shifter to ',f8.2,' ****'/)
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9827 format(1x,'rising diff - de = ', f20.10)
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9825 format(/1x,'*** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffda,diffdpa,diffdb,diffdpb = ',4f16.8)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9199 format(//1x,104('+')/
     + 1x,'CONVERGENCE / TESTER ANALYSIS',32x,'From DIIS'/
     + 1x,104('+')/
     + 1x, 'Iter.', '   Alpha    (somo/vmo)         Beta (somo/vmo)',
     + 5x,          '   Alpha    (somo/vmo)     Beta (somo/vmo)'/)
 9198 format(1x,i5,2(f11.7,' (',2i4,')',2x), '(*)', 2(f11.7,' (',
     +       2i4,')') )
 9298 format(1x,i5,2(f11.7,' (',2i4,')',2x), 3x, 2(f11.7,' (',2i4,')'),
     +       ' (*)' )
 9197 format(1x,104('+')//)
 7777 format(//
     *' end of uhf-scf at ',f12.2,' seconds',a10,' wall'//
     *1x,104('-')//)
 9008 format(/1x,'----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9028 format(/15x,'convergence data'/15x,16('=')//
     +     1x,'maximum number of iterations = ',i6,/,
     +     1x,'method of convergence        = ',i6,/,
     +     1x,'convergence criterion        =1.0e-',i2,/,
     +     1x,'punch out option             = ',i6/1x,124('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',5x,'virtual',1x,'damping',12x,'diis',
     * 3x,'del(t)',5x,'time'/
     * 18x,'energy',10x,'energy',38x,'shift'/
     * 72x,'  -a- ','  -b- ' /1x,124('='))
 9029 format(/15x,'convergence data'/15x,16('=')//
     +     1x,'maximum number of iterations = ',i6,/,
     +     1x,'method of convergence        = ',i6,/,
     +     1x,'convergence criterion        =1.0e-',i2,/,
     +     1x,'punch out option             = ',i6/1x,106('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',5x,'virtual',1x,'damping',12x,'diis'/
     * 18x,'energy',10x,'energy',38x,'shift'/
     * 72x,'  -a- ','  -b- ' /1x,106('='))
 9048 format(20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/ 20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23('-'))
 9108 format(//1x,104('-')//
     *50x,'----- alpha set -----')
 9128 format(//1x,104('-')//
     *50x,'----- beta set -----')
 9188 format(' core assignment '/1x,'i10, i20, i21,',
     +     ' i30, i31, i40, i50, i60, i70, i80, i81, i90 = '/ 13i8/
     +     ' last = ',i8)
 9208 format(//50x,12('=')/50x,'eigenvectors'/50x,12('='))
 9228 format(/10x,14('-')/10x,'density matrix'/10x,14('-'))
 9248 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,2f6.2,f8.3,f16.9,
     1       2f9.2)
 9249 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,2f6.2,f8.3,f16.9)
 9268 format(/10x,15('-')/10x,'energy converged'/10x,15('-'))
 9328 format(/10x,30('-')/10x,'excessive number of iterations'/
     +        10x,30('-'))
 9329 format(/10x,30('-')/
     +        10x,'excessive number of iterations'/
     +        10x,'but flag SCF convergence'/
     +        10x,30('-'))
 9348 format(//1x,104('-')//
     *40x,19('*')/40x,'uhf scf calculation'/40x,19('*') )
 9365 format(
     +10x,'kinetic energy             ',f20.10/ 
     +10x,'electronic energy          ',f20.10/ 
     +10x,'nuclear energy             ',f20.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy               ',f20.10,8x,f10.7/
     +10x,'convergence on density     ',f20.10)
 9368 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9371 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total qm energy            ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9373 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     +,f12.2,' seconds',a10,' wall'/10x,14('-')//)
 9374 format(
     +10x,'number of quadrature points',i9/
     +10x,'integrated electron count  ',f20.10,5x,
     +    'relative error ',e10.2/
     +10x,'XC energy                  ',f20.10)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'=============================='/)
      end
      subroutine spin(sz,s2,da,db,s,d,t,ia,l2)
c
c     ----- calculate expectation value of -sz- and -s**2
c           for the unrestricted hf wavefunction -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension da(*),db(*),s(*),d(*),t(*),ia(*)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
c
c     ----- read in overlap and density matrices (alpha+beta) -----
c
      data dzero,two /0.0d0,2.0d0/
      call rdedx(s,l2,ibl7s,num8)
      call rdedx(da,l2,ibl3pa,idaf)
c
c     ----- d = s*da*s -----
c
      call rdedx(db,l2,ibl3pb,idaf)
      do 280 j = 1,num
      do 200 i = 1,num
      dum = dzero
      do 180 k = 1,num
      ik = ia(max(i,k))+min(i,k)
      jk = ia(max(j,k))+min(j,k)
      dum = dum+da(ik)*s(jk)
  180 continue
  200 t(i) = dum
      do 280 i = 1,j
      dum = dzero
      do 260 k = 1,num
      ik = ia(max(i,k))+min(i,k)
      dum = dum+s(ik)*t(k)
  260 continue
      ij = ia(j)+i
c
c     ----- calculate spin quantum numbers -----
c
  280 d(ij) = dum
      sz =  dfloat(na-nb)/two
      s2 = sz*sz+ dfloat(na+nb)/two-tracep(db,d,num)
      write (iwr,9008) sz,s2
      return
 9008 format(/,10x,17('-'),/,10x,11hspin sz   = ,f6.3,/, 10x,
     +     11hs-squared = ,f6.3,/,10x,17('-'))
      end
      subroutine makfv(f,v,fv,t,ia,m1,m2,n,ndim)
c
c     ----- fv = f * v -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension fv(ndim,*),v(ndim,*),f(*),t(*),ia(*)
      data dzero /0.0d0/
      do 140 j = m1,m2
      do 120 i = 1,n
      dum = dzero
      do 100 k = 1,n
      ii = max(i,k)
      kk = min(i,k)
      ik = ia(ii)+kk
  100 dum = dum+f(ik)*v(k,j)
  120 t(i) = dum
      call dcopy(n,t,1,fv(1,j),1)
 140  continue
      return
      end
      subroutine makeij(v,fv,e,t,m,m0,n,n0,ndim)
c
c     ----- e = v * fv -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension fv(ndim,*),v(ndim,*),e(ndim,*),t(*)
      do 140 j = 1,m
      do 120 i = 1,m0
      t(i) = ddot(n,v(1,i),1,fv(1,j),1)
 120  continue
      do 140 i = 1,m0
      dum = t(i)
      if (j .le. n0 .and. i .gt. n0) dum = dum+dum
  140 e(i,j) = dum
      return
      end
      subroutine eijout(v,e,n,m,ndim,iw)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(ndim,*),e(*)
      mmax = 0
  100 mmin = mmax+1
      mmax = mmax+7
      if (mmax .gt. m) mmax = m
      write (iw,9008)
      write (iw,9028) (j,j = mmin,mmax)
      write (iw,9008)
      write (iw,9048) (e(j),j = mmin,mmax)
      write (iw,9008)
      do 120 i = 1,n
  120 write (iw,9068) i,(v(i,j),j = mmin,mmax)
      if (mmax .lt. m) go to 100
      return
 9008 format(/)
 9028 format(15x,7(6x,i3,6x))
 9048 format(15x,7f15.10)
 9068 format(10x,i5,7f15.10)
      end
      subroutine symh(f,h,ia,move,npair)
c
c     ----- symmetrize the skeleton fock matrix
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension f(*),h(*),ia(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/infoa)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/rdiisa(270),idiisa(24),rdiisb(270),idiisb(24),
     + nsss(30),
     + igvbo(maxorb),igvb1(maxorb),igsp(4),ekk(63),intci(150)
      common/scra/iso(1)
      common/blkin/pxyz(4),t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension mi(48)
      if (npair.ne.0) go to 600
      if (nt .eq. 1) go to 600
      call vclr(h,1,nx)
c
c     ----- find a block (i,j)
c
      do 520 ii = 1,nshell
      do 140 itr = 1,nt
      ish = iso(ii+iliso(itr))
      if (ish .gt. ii) go to 520
  140 mi(itr) = ish
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      do 500 jj = 1,ii
      do 200 itr = 1,nt
      jsh = iso(jj+iliso(itr))
      if (jsh .gt. ii) go to 500
      ish = mi(itr)
      if (ish .ge. jsh) go to 180
      n = ish
      ish = jsh
      jsh = n
  180 if (ish .eq. ii .and. jsh .gt. jj) go to 500
  200 continue
      ljt = ktype(jj)
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      oiandj = ii .eq. jj
c
c     ----- find the equivalent blocks -----
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
c
      jmax = maxj
      do 300 itr = 1,nt
      ntr = itr
      kk = mi(itr)
      ll = iso(jj+iliso(itr))
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 260 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
      do 260 l = minj,jmax
      kl = ia(max(lck,locl+l))+min(lck,locl+l)
      t(k,l) = f(kl)
      if (oiandj) t(l,k) = f(kl)
  260 continue
      if(lit.gt.1.or.ljt.gt.1) call rhr
      do 280 i = mini,maxi
      lci = ia(loci+i)+locj
      if (oiandj) jmax = i
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
      do 280 j = minj,jmax
      ij = lci+j
  280 h(ij) = h(ij)+t(i,j)
c
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
c
  300 continue
      do 480 itr = 2,nt
      kk = mi(itr)
      ll = iso(jj+iliso(itr))
      if (kk .ge. ll) go to 320
      k = ll
      l = kk
      go to 340
  320 k = kk
      l = ll
  340 if (k .eq. ii .and. l .eq. jj) go to 480
      ntr = itr+1
      if (ntr .gt. nt) go to 400
      do 380 it = ntr,nt
      i = mi(it)
      j = iso(jj+iliso(it))
      if (i .ge. j) go to 360
      ij = i
      i = j
      j = ij
  360 if (i .eq. k .and. j .eq. l) go to 480
  380 continue
  400 continue
      ntr = invt(itr)
      do 420 i = mini,maxi
      lci = ia(loci+i)+locj
      if (oiandj) jmax = i
      do 420 j = minj,jmax
      t(i,j) = h(lci+j)
      if (oiandj) t(j,i) = h(lci+j)
  420 continue
      if(lit.gt.1.or.ljt.gt.1) call rhr
      lock = kloc(kk)-kmin(kk)
      locl = kloc(ll)-kmin(ll)
      do 460 k = mini,maxi
      lck = lock+k
      if (oiandj) jmax = k
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
      do 460 l = minj,jmax
      kl = ia(max(lck,locl+l))+min(lck,locl+l)
  460 h(kl) = t(k,l)
  480 continue
  500 continue
  520 continue
      dum = 1.0d0/ dfloat(nt)
_IF1(civu)      call scaler(nx,dum,f,h)
_IFN1(civu)      call vsmul(h,1,dum,f,1,nx)
 600  continue
      if(move.ne.0) then
        call dcopy(nx,f,1,h,1)
      endif
      return
      end
      subroutine rhr
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/rdiisa(270),idiisa(24),rdiisb(270),idiisb(24),
     + nsss(30),
     + igvbo(maxorb),igvb1(maxorb),igsp(4),ekk(63),intci(150)
      common/blkin/pxyz(4),t(35,35),mink,maxk,lkt,minl,maxl,llt,ntr
c
      dimension v(35)
c
c     ----- right multiply  t  by  r,
c           result back in  t
c
      data dzero /0.0d0/
c
      go to (260,180,100,50,10),llt
c
c     ----- g shell
c
  10  ng=15*(ntr-1)-20
      do 30 k=mink,maxk
      do 20 l=21,35
      dum= dzero
      do 110 n=21,35
  110 dum=dum+t(k,n)*gtr(n-20,ng+l)
  20  v(l)=dum
      do 30 l=21,35
  30  t(k,l)=v(l)
      go to 260
c
c     ----- f shell
c
   50 nf=10*(ntr-1)-10
      do 60 k=mink,maxk
      do  80 l=11,20
      dum = dzero
      do  70 n=11,20
   70 dum=dum+t(k,n)*ftr(n-10,nf+l)
   80 v(l)=dum
      do 60 l=11,20
   60 t(k,l)=v(l)
      go to 260
c
c     ----- d shell
c
  100 nd = 6*(ntr-1)- 4
      do 160 k = mink,maxk
      do 140 l = 5,10
      dum = dzero
      do 120 n = 5,10
  120 dum = dum+t(k,n)*dtr(n-4,nd+l)
  140 v(l) = dum
      do 160 l = 5,10
  160 t(k,l) = v(l)
c
c     ----- p shell
c
      go to 260
  180 np = 3*(ntr-1)- 1
      do 240 k = mink,maxk
      do 220 l = 2,4
      dum = dzero
      do 200 n = 2,4
  200 dum = dum+t(k,n)*ptr(n-1,np+l)
  220 v(l) = dum
      do 240 l = 2,4
  240 t(k,l) = v(l)
c
c     ----- left multiply  t  by r
c           result back in  t
c
  260 continue
c
      go to (440,360,280,500,600),lkt
c
c     ----- g shell
c
  600 ng=15*(ntr-1)-20
      do 630 l=minl,maxl
      do 620 k=21,35
      dum=dzero
      do 610 n=21,35
  610 dum=dum+gtr(n-20,ng+k)*t(n,l)
  620 v(k)=dum
      do 630 k=21,35
  630 t(k,l)=v(k)
      go to 440
c
c     ----- f shell
c
  500 nf=10*(ntr-1)-10
      do 530 l=minl,maxl
      do 520 k=11,20
      dum = dzero
      do 510 n=11,20
  510 dum=dum+ftr(n-10,nf+k)*t(n,l)
  520 v(k)=dum
      do 530 k=11,20
  530 t(k,l)=v(k)
      go to 440
c
c     ----- d shell
c
  280 nd = 6*(ntr-1)-4
      do 340 l = minl,maxl
      do 320 k = 5,10
      dum = dzero
      do 300 n = 5,10
  300 dum = dum+dtr(n-4,nd+k)*t(n,l)
  320 v(k) = dum
      do 340 k = 5,10
  340 t(k,l) = v(k)
c
c     ----- p shell
c
      go to 440
  360 np = 3*(ntr-1)- 1
      do 420 l = minl,maxl
      do 400 k = 2,4
      dum = dzero
      do 380 n = 2,4
  380 dum = dum+ptr(n-1,np+k)*t(n,l)
  400 v(k) = dum
      do 420 k = 2,4
  420 t(k,l) = v(k)
  440 continue
      return
      end
      subroutine shiftq(h,nocc,nhocc,n,dee,dep,iterl,ncall,gapa,
     +                  ibrk,gapb)
c
c     ----- shift the diagonal elements of the fock matrix in
c           the molecular basis set by the change in
c           the energy and density matrix.
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scfopt)
INCLUDE(common/mapper)
      dimension h(*)
      data dzero/0.0d0/
      data rbase,rfact,nsame/2.0d0,0.2d0,4/
      data efact,rbade,badef/5.0d0,0.8d0,2.0d0/
      data drat,rbadd,baddf/2.0d0,0.2d0,2.0d0/
      data rdrop,rmin,rmax/2.0d0,0.2d0,1.0d1/
c
      data hoccf/0.5d0/
      temp=gapa
      if(iter.gt.ibrk)temp=gapb
      rmin=temp
      rmax=temp
      rshift=rmin
      pshift=rshift
      if (iter .le. 1) go to 140
      if (ncall .ne. 1) go to 140
      iterl = iterl+1
      oextra = mod(mconv,2) .eq. 0
c
c     ----- decrease the level after nsame iterations -----
c
      oshift = mod(mconv,8) .gt. 3
c
c     ----- use energy as a criterion to determine the level -----
c
      rshc = rbase*rfact**((iter-2)/nsame)
      rshe = dmin1(rbase, dabs(dee)*efact)
      if (dee .gt. dzero) rshe = dmax1(rshe,rbade,rshift*badef)
c
c     ----- use the density matrix difference as a criterion -----
c
      if (dep .gt. dzero .and. dee+dep .gt. dzero) rshe = 
     +     dmax1(rshe,rshift*badef)
      rshd = dzero
c
c     ----- take the highest level shift -----
c
      if (diff .gt. drat
     +     *diffsp .and. rshift .eq. pshift) rshd = 
     +     dmax1(rbadd,rshift*baddf)
      rhigh = dmax1(rshc,rshe,rshd)
c
c     ----- let the level drop after nsame iterations -----
c
      pshift = rshift
      if (mod(iterl,nsame) .eq. 0) rshift = rshift/rdrop
c
c     ----- bound the level between limits -----
c
      rshift = dmax1(rhigh,rshift)
      if (rshift .le. rmin) rshift = rmin
c
c     ----- check for close to convergence or automatic
c           switching to extrapolation -----
c
      if (rshift .gt. rmax) rshift = rmax
      if ( .not. oshift) rshift = dzero
      if (pshift .eq. dzero) go to 100
      if (oextra .and. rshift .lt. vshtol .and.  
     +         dabs(dee) .lt. dmptol) rshift = dzero
      go to 120
c
  100 if (oextra .and.  dabs(dee) .lt. exttol) rshift = dzero
c
c     ----- add the level shift to the fock matrix -----
c
  120 if (rshift .ne. pshift) iterl = 0
  140 continue
c
c     ----- shift any half filled orbitals by rshift*hoccf -----
c
      if (nhocc .eq. 0) go to 180
      is = nocc+1
      ie = nhocc
      do 160 i = is,ie
      ii = ikyp(i)
      h(ii) = h(ii)+rshift*hoccf
c
c     ----- shift any virtual orbitals by rshift -----
c
  160 continue
  180 is = nocc+nhocc+1
      ie = n
      if (is .gt. ie) go to 220
      do 200 i = is,ie
      ii = ikyp(i)
      h(ii) = h(ii)+rshift
  200 continue
  220 if (ncall .eq. 1) diffsp = diff
      return
      end
      subroutine analmo(q,eig,pop,ilifq,ncol,nbasis)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),eig(*),pop(*),ilifq(*)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/tran)
INCLUDE(common/machin)
INCLUDE(common/fpinfo)
INCLUDE(common/harmon)
      common/craypk/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      dimension nsymm(8),msymm(8)
      data m51,small/51,1.0d-3/
      data ev/27.21165d0/
c
      if (ocifp) then
        write(iwr,15)
      else
        write(iwr,10)
      endif
 10   format(/1x,56('=')/
     *'  m.o. irrep        orbital         orbital       orbital'/
     *'              energy (a.u.)   energy (e.v.)     occupancy'/
     *1x,56('='))
 15   format(/1x,31('=')/
     *'  m.o. irrep   orbital occupancy'/
     *1x,31('='))
      nav = lenwrd()
      if(otran)go to 20
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
      do 30 i=1,ncol
      ii=ilifq(i)
_IF1(v)      ibig=1
_IF1(v)      bigg=dabs(q(ii+1))
_IF1(v)      do 40 j=2,nbasis
_IF1(v)      if(dabs(q(ii+j)).gt.bigg) then
_IF1(v)      bigg=dabs(q(ii+j))
_IF1(v)      ibig=j
_IF1(v)      endif
_IF1(v) 40   continue
_IFN1(v)      ibig=idamax(nbasis,q(ii+1),1)
      isymmo(i)=isymao(ibig)
        bigg = 0.0d0
        ibig = 0
        do 50 j=1,nbasis
           if (isymmo(i).ne.isymao(j).and.dabs(q(ii+j)).gt.bigg) then
              bigg = dabs(q(ii+j))
              ibig = j
           end if
50      continue
        if (bigg.gt.small) write(iwr,55) i,ibig,bigg
55      format(' ** symmetry contamination of mo',i4,' at sabf',i4,
     *         ' of ',e12.5,'  **')
30    continue
c     final check    aos versus mos
      do i=1,8
       msymm(i)=0
       nsymm(i)=0
      enddo
c     build up ao array
      do i=1,nbasis
       j=isymao(i)
       nsymm(j)=nsymm(j)+1
      enddo
c     build up mo array
      do i=1,ncol
       j=isymmo(i)
       msymm(j)=msymm(j)+1
      enddo
c      now cross check
      do 170 i=1,8
      if(msymm(i).eq.nsymm(i))goto170
      if ((oharm.or.odepen).and.msymm(i).eq.nsym0(i)) go to 170
      write(iwr,140)i,nsymm(i),msymm(i)
 140  format(//
     *' *** warning -error in m.o. symmetry designation'//
     *4x,'irrep. ',i2,' no. of a.o.s =',i3/
     *14x,'no. of m.o.s =',i3/)
 170   continue
 20   if(.not.ocifp) then
       do i=1,ncol
       eigev = eig(i) * ev
        if(otran) then
         write(iwr,115)i,eig(i),eigev,pop(i)
        else
         write(iwr,125)i,isymmo(i),eig(i),eigev,pop(i)
        endif
       enddo
       write(iwr,135)
      else
       do  i=1,ncol
        if(otran) then
         write(iwr,110)i,pop(i)
        else
         write(iwr,120)i,isymmo(i),pop(i)
        endif
       enddo
       write(iwr,130)
      endif
c
 110  format(1x,i4,6x,f16.8,f20.7)
 120  format(1x,i4,i6,f16.8,f20.7)
 115  format(1x,i4,6x,f16.8,f16.4,f14.4)
 125  format(1x,i4,i6,f16.8,f16.4,f14.4)
c
 130  format(1x,31('='))
 135  format(1x,56('='))
c
c     revised section 190 with mo symmetries
c
      if(.not.otran)call wrt3i(mmmm,mach(13)*nav,iblk51,idaf)
c
      return
      end
      subroutine diiso(h0,h1,h2,h3,trans,ndaf,lockt,numscr)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension h0(*),h1(*),h2(*),h3(*),trans(*)
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/atmol3)
INCLUDE(common/scfwfn)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),iposit(20),
     + nstore,mp,ondiis,nspaca,stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
c
c diiso is specifically for the gvb open shell case
c the error matrix is the lower triangle of gradient
c matrix elements in the mo basis. this is then transformed back
c to the ao basis and written out to ed7
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting block for the vectors .. need space for 20
c                                           .. 40 for uhf
c accdi1  diis comes in only when diff is less than this
c
c num8   the stream for ed7
c num    the dimension of the arrays
c nx     the size of the arrays
c
INCLUDE(common/prints)
_IF(parallel)
      common/scftim/tdiag(4),tdiis,tmult(5)
      dumtim=dclock()
_ENDIF
      derror=0.0d0
      l3=num*num
      nmin=3
      iblk=ndaf
      nmax=5
      njk=2*(nham-ncores)+1
      length=lensec(nx)
      len2=lensec(l3)
c
c --- should we begin storing diis info ???
c
      ibf=msympr(1)-1
      do 5 i=1,nsos
      ibf=ibf+1
      ish=nconf(i)
      jbf=msympr(1)-1
      do 4 j=1,i
      jbf=jbf+1
      jsh=nconf(j)
      ij=iky(ibf)+jbf
      if(ibf.eq.jbf) h0(ij)=0.0d0
      if(ish.eq.jsh) h0(ij)=0.0d0
 4    h1(ij)=h0(ij)
 5    continue
      if( dabs(diff).gt.accdi1.or.diff.eq.0.0d0) goto 300
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax) nmin=nmax-1
      if(nmin.lt.1) goto 300
c
c ------ transform back to the ao basis
c
      lock=4
_IF1(f)      call mtrans(trans,1,h3,1,num,num)
_IFN1(f)      call dagger(num,num,trans,num,h3,num)
      call mult2(h3,h2,h1,num,num,num)
      call rdedx(h1,nx,ibl7s,num8)
      call square(h3,h1,num,num)
      call mult2(h3,h1,h2,num,num,num)
c
c ------ read in symmetry adapted vectors
c
      call rdedx(trans,l3,ibl3qa,idaf)
      nstore=nstore+1
      call search(iblk,numscr)
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.ge.1) then
         njkl=njk*length+len2
        do 50 i=1,ipos1
        ibl=iposit(i)+njkl
        call rdedx(h3,nx,ibl,numscr)
        mp=mp+1
 50     st(mp)=tracep(h1,h3,num)
      endif
      mp=mp+1
      st(mp)=tracep(h1,h1,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
      iposit(ipos)=iposun(numscr)
      npos=iposit(ipos)
      do 105 k=1,njk
      call rdedx(h0,nx,iojkao(k),num8)
      call wrt3(h0,nx,npos     ,numscr)
 105  npos=npos+length
      call wrt3s(trans,l3,numscr)
      call wrt3s(h1,nx,numscr)
      ipos1=ipos+1
      if(ipos1.le.itot) then
        do 110 i=ipos1,itot
        ibl=iposit(i)+njk*length+len2
        call rdedx(h3,nx,ibl,numscr)
        mp=mp+i-1
        st(mp)=tracep(h1,h3,num)
 110    continue
      endif
c
c --- now solve the diis equations
c
      if(nstore.lt.nmin)goto 400
      if(mod(nstore,nmin).ne.0) goto 400
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h2,h3,ifail)
c     call fixnag(ifail,'f04atf-diiso')
      if(ifail.ne.0.and.oprint(47))write(iwr,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 300
      do 114 i=1,ndim
 114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
      len=0
      do 160 ijk=1,njk
      call vclr(h0,1,nx)
      do 150 k=1,itot
      ibl=iposit(k)+len
      call rdedx(h1,nx,ibl,numscr)
      call daxpy(nx,ct(k),h1,1,h0,1)
 150  continue
      len=len+length
 160   call wrt3(h0,nx,iojkao(ijk),num8)
c
c ------ now interpolate the vectors
c
      call vclr(trans,1,l3)
      len=njk*length
      do 180 k=1,itot
      ibl=iposit(k)+len
      call rdedx(h3,l3,ibl,numscr)
      call daxpy(l3,ct(k),h3,1,trans,1)
 180  continue
      call wrt3(trans,l3,ibl3qa,idaf)
      if(.not.ondiis) kcount=0
      ondiis=.true.
      go to 9999
 300  nstore=0
      itot=0
      mp=0
      ondiis=.false.
      lock=lockt
      go to 9999
 400  ondiis=.false.
9999  continue
_IF(parallel)
      tdiis=tdiis+(dclock()-dumtim)
_ENDIF
      return
      end
      subroutine diisu(h0,h1,h2,h3,ndaf,diffa,diffb,lockt,
     * iblkha,iblkhb,numscr,idomoa,ivmoa,idomob,ivmob)
c
c --- subroutine sets up and solves the diis equations
c --- direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension h0(*),h1(*),h2(*),h3(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/harmon)
INCLUDE(common/prints)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),iposit(20),
     +  nstore,mp,ondiis
c
c diisu is specifically for  uhf cases
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and wriiten out to ed7
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting block for the vectors .. need space for 40
c ndafd   the block which holds the beginning of the extrapolated vecs
c accdi1  diis comes in only when diff is less than this
c accdi2  the diis solution is only used when the residuals
c          are less than thia
c
c num8   the stream for ed7
c l1     the dimension of the arrays
c l2     the size of the arrays
c
_IF(parallel)
      common/scftim/tdiag(4),tdiis,tmult(5)
      dumtim=dclock()
_ENDIF
      derror=0.0d0
      l3=num*num
      num0 = newbas0
      nmin=3
      nmax=8
      ipass=1
      length=lensec(nx)
      len2=length+length
      len4=len2+len2
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1) goto 300
      diffa=0.0d0
      idomoa = 0
      ivmoa = 0
      diffb=0.0d0
      idomob = 0
      ivmob = 0
      ns1=nstore
      iposa=mod(ns1,nmax)
      iblk=iposa*len4+ndaf
      iposit(iposa+1)=iblk
      nelec=na
      ibq=ibl3qa
      ibf=iblkha
      go to 3
 2    ibq=ibl3qb
      ibf=iblkhb
      iblk=iblk+len2
      nelec=nb
3     call rdedx(h3,l3,ibq,idaf)
c
c ----- sort out the vectors
c
      call tdown(h3,ilifq,h3,ilifq,num0)
c
c ----- calculate the error vector
c
      call rdedx(h0,nx,ibf,num8)
      call dcopy(nx,h0,1,h2,1)
      call mult2(h3,h1,h2,num0,num0,num)
      ij=0
      diff0=0.0d0
      do 5 i=1,num0
      do 4 j=1,i
      ij=ij+1
      if(i.le.nelec.and.j.le.nelec) h1(ij)=0.0d0
      if(i.gt.nelec.and.j.gt.nelec) h1(ij)=0.0d0
      if( dabs(h1(ij)).gt.diff0) then
       idumj = j
       idumi = i
       diff0= dabs(h1(ij))
      endif
 4    continue
 5    h1(ij)=0.0d0
c
c ------ transform back to the ao basis
c
      do 10 i=1,num
      do 10 j=1,i
      dum=h3(ilifq(i)+j)
      h3(ilifq(i)+j)=h3(ilifq(j)+i)
 10   h3(ilifq(j)+i)=dum
      call mult2(h3,h2,h1,num,num,num)
      call rdedx(h1,nx,ibl7s,num8)
      call square(h3,h1,num,num)
      call mult2(h3,h1,h2,num,num,num)
      call search(iblk,numscr)
      call wrt3s(h0,nx,numscr)
      call wrt3s(h1,nx,numscr)
      ipass=ipass+1
      if(ipass.le.2) then
       diffa = diff0
       idomoa = idumj
       ivmoa = idumi
       go to 2
      else
       diffb = diff0
       idomob = idumj
       ivmob = idumi
      endif
      if(diffa.gt.accdi1.or.diffb.gt.accdi1)go to 300
      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)
      ibla=iposit(ipos)+length
      call rdedx(h0,nx,ibla,numscr)
      if(ipos1.ge.1) then
        do 50 i=1,ipos1
        ibla=iposit(i)+length
        call rdedx(h2,nx,ibla,numscr)
        iblb=ibla+len2
        call rdedx(h3,nx,iblb,numscr)
        mp=mp+1
        st(mp)=tracep(h1,h3,num)+tracep(h0,h2,num)
 50   continue
      endif
      mp=mp+1
      st(mp)=tracep(h1,h1,num)+tracep(h0,h0,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
_IF1()      iposit(ipos) = iposun(numscr)
      ipos1=ipos+1
      if(ipos1.le.itot) then
        do 110 i=ipos1,itot
        ibla=iposit(i)+length
        call rdedx(h2,nx,ibla,numscr)
        iblb=ibla+len2
        call rdedx(h3,nx,iblb,numscr)
        mp=mp+i-1
        st(mp)=tracep(h1,h3,num)+tracep(h0,h2,num)
 110    continue
      endif
c
c --- now solve the diis equations
c
      if(nstore.le.nmin)goto 400
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h2,h3,ifail)
c     call fixnag(ifail,'f04atf-diisu')
      if(ifail.ne.0.and.oprint(47))write(iwr,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 210
      do 114 i=1,ndim
  114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
c     del=0.0d0
c
      call vclr(h0,1,nx)
      call vclr(h1,1,nx)
      do 150 k=1,itot
      ibla=iposit(k)
      call rdedx(h2,nx,ibla,numscr)
      iblb=iposit(k)+len2
      call rdedx(h3,nx,iblb,numscr)
      call daxpy(nx,ct(k),h3,1,h1,1)
      call daxpy(nx,ct(k),h2,1,h0,1)
 150  continue
      if(.not.ondiis) kcount=0
      ondiis=.true.
      return
 210  ibl=iposit(ipos)
      call rdedx(h0,nx,ibl,numscr)
      call reads(h1,nx,numscr)
      call reads(h2,nx,numscr)
      call reads(h3,nx,numscr)
      call wrt3(h0,nx,ndaf,numscr)
      call wrt3s(h1,nx,numscr)
      call wrt3s(h2,nx,numscr)
      call wrt3s(h3,nx,numscr)
      nstore=1
      itot=1
      mp=iky(ipos)+ipos
      st(1)=st(mp)
      ondiis=.false.
      go to 9999
 300  nstore=0
      itot=0
      mp=0
 400  ondiis=.false.
      call rdedx(h0,nx,iblkha,num8)
      call rdedx(h1,nx,iblkhb,num8)
      go to 9999
9999  continue
_IF(parallel)
      tdiis=tdiis+(dclock()-dumtim)
_ENDIF
      return
      end
      subroutine drhfcl(q,coord,atmchg)
_IF(parallel)
c
c  parallel direct SCF 
c
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c  
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes 
c           run through all code section
_ENDIF
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      character*127 text
_IF(drf)
cahv ##################################################################
INCLUDE(../drf/comdrf/sizesrf)
cahv ##################################################################
_ENDIF
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/atmol3)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/scra7)
INCLUDE(common/disc)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
INCLUDE(common/timeperiods)
INCLUDE(common/runlab)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/xfield)
INCLUDE(common/psscrf)
INCLUDE(common/harmon)
INCLUDE(common/fermidirac)
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     +             iposit(20),nsti(2),ondiis,junkj
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      dimension q(*),coord(3,*),atmchg(*)
c     dimension osign(2)
_IF(ga)
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
c ... for dummy symass section
      parameter(isymtp=99)
      character*20 label
      common/msglab/label
_ENDIF
c
       integer domo, vmo
       logical logtest
       common/testerp/tester(200), idomo(200), ivmo(200),
     +               testerd(200),idomod(200), ivmod(200),
     +               testdum(2,200),idomodum(2,200),
     +               logtest(200)
c
INCLUDE(common/drfopt)
_IF(drf)
cahv ##################################################################
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/opt)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfbem)
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
      common/enrghf/enhh,etothh,ehfhh
      logical odrf, oarf
cahv ##################################################################
_ENDIF
INCLUDE(common/zorac)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/drive_dft)
INCLUDE(common/statis)
_ENDIF
_IF(ga)
      character*1 xn,xt
_ENDIF
      character*10 charwall
      Logical use_symmetry, force_serial
      Integer n_irreps
      character *5 fnm
      character *6 snm

_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      data fnm/'scf.m'/
      data snm/'drhfcl'/
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data bigd/1.0d4/
            call check_feature('drhfcl')
_IF(drf)
cahv
      odrf = .false. 
      oarf = .false.
      enadd = dzero
cahv
_ENDIF
      odft = .false.
c
_IF(ccpdft)
      idum = CD_update_geom(c)
      odft = CD_active()
_ENDIF
      dft_accu = 1.0d0
c
c     variables for monitoring tester + convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmn = 0.0d0
      esmear = esmear_start
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      out = nprint .eq. 5
      outon = nprint .ne. -5
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)

      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)

cl
c     is delta density ever to be invoked? if not, we can
c     reduce i/o activity around section 471 ..
c
c     odelta = deltol .gt. dlog(acurcy)
c
_IF(parallel)
c
c - this should be default anyway -
c
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
*     write(6,*)'master',ipg_nodeid(),omaster
    
_ENDIF
c
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      lentri = l2
      l3 = num*num
c     l4 = num * 9
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
c
c note memory allocation is designed to make irdmat and 
c idmat contiguous - reversed
c
      ilen=0
      if(ocryst)ilen=l2
      ilen = ilen + 1 + ikyp(nshell)*3 + 6 * l2 + 2 * l1

      irdmat = igmem_alloc_inf(ilen,fnm,snm,'rdmat',IGMEM_NORMAL)

c      if(opg_root())write(6,*)'SCF start address',irdmat, 
c     &  ' len ',ilen

      iprefa=irdmat+ikyp(nshell)
      i10=iprefa+ikyp(nshell)
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      i50 = i40+l2
c ps: add extra space for brdcst buffer 
      i60 = i50+l2+1+ikyp(nshell)
      i80 = i60+l2
      i90 = i80+l1
      ix  = i90+l1
      isiz=0
      if(ocryst)isiz=l2
      last = ix + isiz

      length = last-irdmat
c
c consistency check
      if (length .ne. ilen)then
         write(iwr,*)length,last,irdmat,ilen
         call caserr('size error')
      endif

      if (nprint .eq. 5) write (iwr,9308) i10,i20,i21,i30,i31,i40,i50,
     +     i60,i80,i90,last
_IF(drf)
cahv --- DRF extension ---                                                  
      itdrf = 0
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
c
c-----  set acur for use in 2-electron rf routines drfhamx
c
        acur = acurcy
        if (iarfcal .eq. 0) then
          odrf = .true.
        else
          oarf = .true.
        endif
c
c      ovlap at i110  from da10, 12, l2
c      dipx  at i120  from da10, 53, l2
c      dipy  at i130  from da10, 54, l2
c      dipz  at i140  from da10, 55, l2
c      omega(s) at i150    da31, 50, nomega
c      omega(op) at i160   da31, 51, nomega, if itwoeps=1
c      iexp  at i170       da31, 2,  l2
c      ijbit(s) at i180    da31, 56, l2
c      ijbit(op) at i190   da31, 57, l2, if itwoeps=1 else clearit
c      fdrf      at i200   calc in drfhamc, l2
c      ddrf      at i210   l2
cafc
        i110 = igmem_alloc_inf(l2,fnm,snm,'ovlap',IGMEM_NORMAL)
        i120 = igmem_alloc_inf(l2,fnm,snm,'dipx',IGMEM_NORMAL)
        i130 = igmem_alloc_inf(l2,fnm,snm,'dipy',IGMEM_NORMAL)
        i140 = igmem_alloc_inf(l2,fnm,snm,'dipz',IGMEM_NORMAL)
        if (odrf) then
          i150 = igmem_alloc_inf(nomga,fnm,snm,'omega_s',IGMEM_NORMAL)
          if (itwoeps .eq. 1) then
            i160 = igmem_alloc_inf(nomga,fnm,snm,'omega_op',
     &                             IGMEM_NORMAL)
          else
            i160 = igmem_alloc_inf(1,fnm,snm,'omega_op',IGMEM_NORMAL)
          endif
          i180 = igmem_alloc_inf(l2,fnm,snm,'ijbit_s',IGMEM_NORMAL)
          i190 = igmem_alloc_inf(l2,fnm,snm,'ijbit_op',IGMEM_NORMAL)
        else
          i150 = igmem_alloc_inf(ndim*ndim,fnm,snm,'omega_s',
     &                           IGMEM_NORMAL)
          i160 = igmem_alloc_inf(ndim,fnm,snm,'omega_op',
     &                           IGMEM_NORMAL)
          i180 = igmem_alloc_inf(nwtr*nwtc,fnm,snm,'ijbit_s',
     &                           IGMEM_NORMAL)
          i190 = igmem_alloc_inf(ndim,fnm,snm,'ijbit_op',IGMEM_NORMAL)
        endif
        i200 = igmem_alloc_inf(l2,fnm,snm,'fdrf',IGMEM_NORMAL)
        i210 = igmem_alloc_inf(l2,fnm,snm,'ddrf',IGMEM_NORMAL)
        i170 = igmem_alloc_inf(l2,fnm,snm,'iexp',IGMEM_NORMAL)
        if (oarf) then
          i193 = igmem_alloc_inf(ndim,fnm,snm,'indmom',IGMEM_DEBUG)
          i196 = igmem_alloc_inf(ndim,fnm,snm,'sf',IGMEM_DEBUG)
        endif
      endif
cahv
_ENDIF
c ***
      dlnmxd=0.0d0
      dlntol=tolitr(3)
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i30),l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,q(i20),q(i30))
      endif
_ENDIF

      call start_time_period(TP_RDMAT)
c     specification of section isect(471)
c     only output delta-difference data if is this is to be invoked
c     i.e. only rdmat and prefac traffic involved in default
c *** section isect(471) now holds:
c *** prefac mat., reduced dmat, delta f, delta d, last f, current d

      m171t=171
      len171=  lensec(l2)
      nshtri=ikyp(nshell)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      do i=4,6
       iof171(i)=iof171(i-1)+len171
      enddo
      if(irest.ne.1.and.irest.ne.2) then
        call rdmake(q(iprefa))
        ln171=nshblk*2
        if(odelta) then
         ln171=4*len171+ln171
        endif
        call secput(isect(471),m171t,ln171,ibl171)
        if(outon) then
        write(iwr,2980) ibl171,ln171
        endif
        call wrt3(q(iprefa),nshtri,ibl171,idaf)
        if(odelta) then
         call zer171(q(i10),l2,4,ibl171+iof171(3),len171,idaf)
        endif
        call clredx
      endif
c
      call end_time_period(TP_RDMAT)
c ***
      call start_time_period(TP_TEST1)
_IF(drf)
c
cahv --- DRF extension ---                                                  
      if (odrf .or. oarf) then
c  -----  read necessary data from da10
        call daread(idafh,ioda,q(i110),l2,12)
        call daread(idafh,ioda,q(i120),l2,53)
        call daread(idafh,ioda,q(i130),l2,54)
        call daread(idafh,ioda,q(i140),l2,55)
c  -----  read necessary data from da31
        call daread(idafdrf,iodadrf,q(i170),l2,2)
        call clear(q(i200),l2)
      endif
      if (odrf) then
        call daread(idafdrf,iodadrf,q(i150),nomga,50)
        call daread(idafdrf,iodadrf,q(i180),l2,56)
        if (itwoeps .eq. 1) then
          call daread(idafdrf,iodadrf,q(i160),nomga,51)
          call daread(idafdrf,iodadrf,q(i190),l2,57)
        else
          call clear(q(i190),l2)
        endif
      else if (oarf) then
        call daread(idafdrf,iodadrf,q(i150),ndim*ndim,43)
        call daread(idafdrf,iodadrf,q(i160),ndim,44)
        call daread(idafdrf,iodadrf,q(i190),ndim,101)
      endif
cahv
_ENDIF
      if (maxcyc .le. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
       tester(loop) = 0.0d0
       testerd(loop) = 0.0d0
       logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1

      i=lensec(l2)
      iblkh0 =ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)
      if(numdis.eq.0)then
         numdis=num8
         ndafd=iblkh+i
      else
         ndafd=ibldis+2
      endif
      ndafdi=ndafd+3*i
      if(irest.eq.0.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         kcount = 0
         iterv = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      diffdp = bigd
      diffdpp = bigd
c
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (odamph .or. oshift) damp = done
      lockt=lock
      tim0=cpulft(1)
      ovir = .false.
      if(oaimpac())ovir = .true.
c
c     ----- nuclear energy
c     ----- psscrf dvsg
c
      if (opssc) then
         en = denuc(nat,atmchg,coord,itotbq,iwr)
      else if (ocryst) then
         en = crnuc(nat,czan,coord,ecrn)
      else
         en = enucf(nat,atmchg,coord)
      endif
      l0=newbas0
      call end_time_period(TP_TEST1)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
       n_irreps = 0
       Do i = 1, l1
          isymmo( i ) = isymao( i )
          n_irreps = Max( n_irreps, isymao( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
       if(outon) then
        write (iwr,9008) en
        if (use_symmetry) write(iwr,9009)
        write (iwr,9028) maxcyc,mconv,nconv,mpunch,
     1  dexp(-dlntol)
       endif
       call start_time_period(TP_ORFOG)
       call rdedx(q(i10),l2,ibl7st,num8)
       call qmat_symm(q,q(i10),q(i20),q(i31),q(i40),iky,l0,l1,l3,l1,out,
     +               isymmo, use_symmetry )
c load guess vectors..
       call rdedx(q(i30),l3,ibl3qa,idaf)
       if (.not.onoor) then
c load overlap
       call rdedx(q(i50),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i50), num, l0)
       else
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     *         ilifq,l0,l1,1)
       endif
_ELSE
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     *         ilifq,l0,l1,1)
_ENDIF
      else
        write(iwr,'(2a)') '** warning start-orbitals not orthogonal **',
     1                    '** result may not be variational         **'
c...    if need be one can pick up the occupation numbers here to
c...    calculate coulomb energies for correlated monomers (jvl,2003)
      end if
      call end_time_period(TP_ORFOG)

      call start_time_period(TP_TEST2)
      call wrt3(q(i30),l3,ibl3qa,idaf)
      call tdown(q(i30),ilifq,q(i30),ilifq,l1)
      call rdedx(q(i90),l1,ibl3ea,idaf)
      ndoc = na
      if (osmear) then
        call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                   iwr,oprint(60))
        emrmn = 2*emrmn
      else
        call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
c...
c...  construct initial density matrix
c...
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
       call wrt3(q(i50),l2,ibl3pa,idaf)
_IF(drf)
cahv --- drf extension ---
      if(odrf .and. (intdrf .eq. 0))
     1   call dcopy(l2,q(i50),1,q(i210),1)
cahv
c
_ENDIF
       call end_time_period(TP_TEST2)

_IF(parallel)
      endif
c...  other nodes receive density matrix in i50 if they were screened out
      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst idmat"
         call pg_brdcst(7129,q(i50),l2*8,0)
         label=" "
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      lprnt = l1
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      call timrem(tlefts)
c *** build rest of fock matrix from 2-e integrals
      ifock=i10
      idmat=i50
      if(irest.eq.1 .or. irest.eq.2) then
        call secget(isect(471),m171t,ibl171)
        write(iwr,2990)irest,ibl171
        call rdedx(q(iprefa),nshtri,ibl171,idaf)
        call reads(q(irdmat),nshtri,idaf)
        if(odelta) then
         call rdedx(q(i10),l2,ibl171+iof171(3),idaf)
         call reads(q(i50),l2,idaf)
        endif
        dlnmxd=-9999999.0d0
        do 1672 kkk=1,nshtri
            if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
            dlnmxd=q(irdmat+kkk-1)
1672    continue

_IFN(parallel)
        if(irest.eq.1) goto 6351
        if(irest.eq.2) goto 6352
_ENDIF
      endif
c ***
_IF(ccpdft)
      if (CD_active()) then
         call retrieve_spare(imemspare)
         idum = CD_set_2e()
         if (CD_2e() .and. (.not. CD_HF_exchange())
     &      .and. (.not. CD_HF_coulomb()) )then
            imemdhstaru = 0
         else
            call mem_dhstaru(imemdhstaru)
         endif
         idum = CD_reset_2e()
         imemfree = igmem_max_memory() - imemspare - imemdhstaru
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q,q,iwr)
         else if (ks_bas.eq.KS_MO) then
            imemreq = CD_memreq_energy_mo(l1,na,0,q,q,iwr)
         else if (ks_bas.eq.KS_AOMO) then
            imemreq = CD_memreq_energy(l1,na,0,q,q,iwr)
         endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then 
            write(iwr,600)ierror
            call caserr('Out of memory in incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
c        if (CD_jfit_incore().and.CD_HF_exchange()) then
c           write(iwr,*)'*** WARNING: direct integral memory usage not',
c    +                  ' accounted for!!!'
c           write(iwr,*)'*** WARNING: Calculation may run out of ',
c    +                  'memory!!!'
c        endif
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
 140  continue
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4)=.false.
         oswed3(8)=.false.
      endif
_ENDIF
c
c  intermediate times
c      call list_time_periods(.false.,.false.)

c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** delta f, delta d, last f, current d, reduced dmat, prefac mat.
c *** note that delta f and delta d are NOT accessed if
c *** delta-SCF has not been requested : odelta = .false.)
      call start_time_period(TP_RDMAT)
_IF(parallel)
      if(omaster) then
      if(irest.eq.1) goto 6351
*     if(irest.eq.2) goto 6352
_ENDIF
      if (odelta) then
         call vclr(q(i60),1,l2)
         call wrt3(q(i60),l2,ibl171+iof171(3),idaf)
         call rdedx(q(i60),l2,ibl171+iof171(6),idaf)
         call wrt3(q(i50),l2,ibl171+iof171(6),idaf)
         call vsub(q(i50),1,q(i60),1,q(i50),1,l2)
         call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
_IF(drf)
cahv ---  DRF extension ---
         if (odrf .and. (itdrf .gt. 0)) 
     +   call dcopy(l2,q(i50),1,q(i210),1)
cahv
_ENDIF
      else
c        call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
_IF(drf)
cahv ---  DRF extension ---
         if (odrf .and. (itdrf .gt. 0)) 
     +     diffaa = cvgden(q(i50),q(i210),l2)
cahv
_ENDIF
      endif
      if(outon) write(iwr,3030)
      if (odelta) then
       call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
      endif
c ***
c *** make sure that iter before delta calc is of full accuracy.
c ***
      if( (dlnmxd.lt.deltol) .and. (iter-1.lt.itrtol(2))) then
        dlnmxd=deltol+1.0d0
      endif
      if( (dlnmxd.gt.deltol) ) then
        if(outon) write(iwr,3000)
        if (odelta) then
           call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
           call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
           call vclr(q(i60),1,l2)
           call wrt3(q(i60),l2,ibl171+iof171(5),idaf)
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        else
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        endif
        call clredx
      endif
      call wrt3(q(irdmat),nshtri,ibl171+iof171(2),idaf)
6351  continue
c ***

_IF(parallel)
      endif
c
      if(ipiomode .eq. IO_NZ_S)then
c ps: combine to 1 brdcst
         call dcopy(nshtri,q(irdmat),1,q(idmat+l2+1),1)
         q(idmat+l2)=dlnmxd
         n4=l2+nshtri+1
         label="brdcst idmat"
         call pg_brdcst(7123,q(idmat),(l2+nshtri+1)*8,0)
         label=" "
         call dcopy(nshtri,q(idmat+l2+1),1,q(irdmat),1)
         dlnmxd= q(idmat+l2)

         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF(parallel)
      call end_time_period(TP_RDMAT)
      if(iter.lt.itrtol(1)) then
         dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
         dlntol=tolitr(2)
      else
         dlntol=tolitr(3)
      endif
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(outon) then
         write(iwr,3020) dlntol
      endif
_IF(ccpdft)
c
c skip dhstar when not required for choice of DFT functional
c
      idum = CD_set_2e()

      if(CD_2e() .and. (.not. CD_HF_exchange())
     &    .and. (.not. CD_HF_coulomb()) )then

        o2e = .false.
        if(irest.ne.1)call vclr(q(ifock),1,l2)
      else
        o2e = .true.
        call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),irest)
      endif

      idum = CD_reset_2e()

_ELSE
      call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),irest)
_ENDIF

_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
      call start_time_period(TP_TEST3)
      if(irest.eq.0) go to 6352
        write(iwr,3010)irest
        if(numdis.eq.num8) go to 460
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
        go to 460
c *** form complete 2-e skelton fmat and save it
6352  if (odelta) then
        call rdedx(q(i60),l2,ibl171+iof171(5),idaf)
        call vadd(q(i10),1,q(i60),1,q(i10),1,l2)
        call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
        call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
      else
c       call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
      endif
c ***
      if (nprint .eq. 5) then
         write (iwr,9068)
         call prtril(q(i10),l1)
      endif

c  symmetrise skeleton fock matrix 
_IF(ccpdft)
c  if only 2e integrals were performed
      if(o2e)then
_ENDIF
_IF(drf)
cahv ---  DRF extension ---
cafc  it's about time to add the drf-contributions to the fock-matrix
cafc  this is done by the routine drfhamc, which uses the diff dens
cafc  by reading it from da31 and places the results to be added at i200
cafc
      if (odrf) then
c
        itdrf = itdrf + 1
        call drfhamc(q(i200),q(i210),q(i110),q(i120),q(i130),
     1    q(i140),q(i170),q(i150),q(i180),q(i160),q(i190))
c --- add to fock matrix   (scffact * fdrf(i) )
        do 1111, i = 1, l2
          q(i10-1+i) = q(i10-1+i) + scffact*q(i200-1+i)
 1111   continue
        call dcopy(l2,q(i50),1,q(i210),1)
      else if (oarf) then
        call arfone(q(i200),q(i50),
     1  q(i110),q(i120),q(i130),q(i140),
     2  q(i150),q(i160),q(i180),q(i170),
     3  q(i190),q(i193),q(i196),enadd)
      endif
cahv
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at q(i10)
c     scratch area at q(i20)
c  symmetrise skeleton fock matrix
_ENDIF
      call symh(q(i10),q(i20),iky,0,0)
      if (nprint .eq. 5) then
         write (iwr,9088)
         call prtril(q(i10),l1)
      endif
_IF(ccpdft)
      endif
_ENDIF
c
c increment H_2 at q(i10) with H_1
c
      call rdedx(q(i20),l2,ibl7f,num8)
c...   add zora corrections
      if (ozora) call zora(q,q,q(i20),'read')
_IF(drf)
      if (oarf) then
        call dcopy(nx,q(i20),1,q(i210),1)
        call vadd(q(i20),1,q(i200),1,q(i20),1,nx)
      endif
_ENDIF
      call adonee(q(i10),q(i20),q(iprefa))
c
c  ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i50),q(i20),l1)
_IF(drf)
      if (oarf) call dcopy(nx,q(i210),1,q(i20),1)
_ENDIF
c
c ps moved this up so as not to overlap with dft timings
c disabled for integral restarts
c
      if(irest.ne.2)call end_time_period(TP_TEST3)
c
c  save previous total 
c
      ehf0 = ehf

_IF(ccpdft)
      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(o2e)then
c
c Coulomb operator is in q(i10) augmented by H_1, 
c compute energy using HF expression without exchange
c
            ehf2 = tracep(q(i50),q(i10),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1
         endif
c
c Update Kohn-Sham matrix and compute fitted/integrated 
c energy terms
c
c        dft_accu = 0.95d0*dft_accu
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(q(isma),1,l2)
_ENDIF
         if (ks_bas.eq.KS_AO) then
            idum = CD_energy_ao(c,q(i10),dum,q(i50),dum,
     +           edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_MO) then
            idum = CD_energy_mo(l1,na,0,c,q(i10),dum,q(i30),dum,
     +           edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q(isma)
_ENDIF
     +           )
         else if (ks_bas.eq.KS_AOMO) then
            idum = CD_energy(l1,na,0,c,q(i10),dum,q(i30),dum,q(i50),dum,
     +           edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q(isma)
_ENDIF
     +           )
         else
            call caserr("drhfcl: illegal ks_bas")
         endif
_IF(debug_S)
         call compare_S(q(isma),q(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF
         call symm_op(q(i10))
         ehf = etmp+edft

         if(out) then
           if(opg_root()) then
            call CD_print_dftresults(.true.,.false.,iwr)
            write(iwr,9407)ehf1, ehf2
           endif
         endif
         call timana(31)
         call cpuwal(begin,ebegin)

      else
c
c  Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
         ehf2 = tracep(q(i50),q(i10),l1)
         ehf = (ehf1+ehf2)*pt5
      endif
_ELSE
      ehf2 = tracep(q(i50),q(i10),l1)
      ehf = (ehf1+ehf2)*pt5
_ENDIF
      ehf = ehf + 0.5d0*emrmn

c for xtal field calculations, the core hamiltonian
c includes the crystal potential, but for the energy 
c all terms modelling molecule-molecule interactions 
c must be divided by two 

      if(ocryst)then
         m53=53   !!!!!!!!
         write(iwr,*)'secget',isecfx,m53
         call secget(isecfx,m53,iblk)
         call rdedx(q(ix),l2,iblk,idaf)
         ehf3 = tracep(q(i50),q(ix),l1)
         write(iwr,*)'ehf3',ehf3
         esave = ehf
         ehf = (ehf1+ehf2-ehf3)*pt5 
        ecre = esave - ehf
         write(iwr,*)'new ehf, ecre',ehf,ecre
      endif
c
c compute kinetic energy and virial if required
c
      if(ovir)call rhfvir(q(i50),ehf,en,ek,vir,num8,l1,l2)

      call wrt3(q(i10),l2,iblkh,num8)
c hack remco
c        rewind(76)
c        print *,'remco l2',l2
c        write(76)(q(i10+i-1),i=1,l2)
c einde hack remco
      if (nprint .eq. 5) then
        write (iwr,9048)
        call prtril(q(i10),l1)
        call rdedx(q(i90),l1,ibl3ea,idaf)
        iorb1 = 1
        iorb2 = nocmx
        call makfv(q(i10),q(i30),q(i30),q(i21),iky,iorb1,iorb2,l1,l1)
        call rdedx(q(i10),l3,ibl3qa,idaf)
        call makeij(q(i10),q(i30),q(i30),q(i21),l0,l0,l1,l0,l1)
        write (iwr,9228)
        call eijout(q(i30),q(i90),l0,nocmx,l1,iwr)
        write (iwr,9328) ehf1,ehf2,ehf
        call rdedx(q(i10),l2,iblkh,num8)
      endif
      iter = iter+1
      etot0 = etot
      etot = ehf+en
_IF(drf)
      etot = etot + enadd
_ENDIF
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
      if (damp .lt. dmpcut) damp = dmpcut
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafd,
     +numdis,iterv,1,1)
      diffpp=diffp
      diffp=diff

      if(odiis) then
      call start_time_period(TP_DIIS)
_IF(ga)
      if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
         call diis_ga(q(i10),q(i30),q(i50),ibl3qa,diffd,domo,vmo)
      else
         call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     *        ndafdi,numdis,diffd,lockt,domo,vmo)
      endif
_ELSE
      call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     *  ndafdi,numdis,diffd,lockt,domo,vmo)
_ENDIF
        if (iter.le.maxcycp) then
         idomod(iter) = domo
         ivmod(iter) =  vmo
         testerd(iter) = diffd
        endif
      call end_time_period(TP_DIIS)
      endif

      call start_time_period(TP_TEST4)

c
c ------ take precautions if tester is increasing again
c
      if(iter.eq.1 .or. .not.odiis) goto 223
      if(ondiis) then
        if (oprind) then
         write(iwr,9824)diffd,diffdp
         write(iwr,9823)diffd, acurcy*diisf
         write(iwr,*) 'de = ', de
        endif
        if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.acurcy*diisf
     1  .or.(otestd.and.de.lt.1.0d-10)) then
*        if (de.lt.1.0d-5) then
          if(iter.le.maxcycp) logtest(iter) = .true.
          diffdpp = diffdp
          diffdp = diffd
          go to 223
*        else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
*         odynamic = .true.
*         accdi1 = 0.0d0
*         itdiis = 0
*         if(oprind) write(iwr,9828) iter
*        endif
        else
          if(oprind) write(iwr,9829) iter
        endif
        nsti(1)=0
        ondiis=.false.
        call rdedx(q(i10),l2,iblkh,num8)
        diffdp = bigd
        diffdpp = bigd
      else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c       if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c        iwild = iwild + 1
c        osign(iwild) = de.lt.0.0d0
c        if (iwild.eq.2.and..not.oreset.and.
c    +       .not.(osign(1).and.osign(2)) ) then
c         gapa1 = 5.0d0
c         gapa2 = gapa1
c         oreset = .true.
c         iwild = 0
c         write(iwr,9831) gapa1
c         maxcyc = maxcyc + maxcyc
c         maxcycp = min(200,maxcyc)
c        endif
c       endif
c
        if (odynamic) then
         if(diff.le.diffp.or.diffp.le.diffpp) then
          itdiis = itdiis + 1
          if (itdiis.eq.3) then
            accdi1 = diff
            if (oprind) write(iwr,9826) accdi1
          endif
         else
          if (oprind) write(iwr,9827)
          itdiis = 0
         endif
        endif
      endif
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -q- at q(i30) orthonormalizing transformation
c     -h- at q(i10)
c    -h'- at q(i50) transformed h matrix
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
 223  if(oshift)then
         call rdedx(q(i30),l3,ibl3qa,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diff=0.0d0
      ii=ndoc+1
      if(ii.le.l0)then
         do 250 i=ii,l0
            ij=iky(i)+i50
            loop=idamax(min(i-1,nocmx),q(ij),1)
            if(loop.gt.0) then
             dum = dabs(q(ij+loop-1))
             if (dum.gt.diff) then
               diff = dum
               if(iter.le.maxcycp) then
                tester(iter) = diff
                idomo(iter) = loop
                ivmo(iter) = i
               endif
             endif
            endif
 250     continue
      endif
      if(ondiis.and..not.osmear) diff=diffd
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
c
cgdf  not sure if needed: convergence check here seems to be for 
c     skipping some work, when the real check comes later
c
      if (osrso) ocvged = .true.
c
_IF(parallel)
      if(ipiomode.eq.IO_NZ_S.and.odpdiag(l1)) then
c  in screened case, broadcast convergence event 
         dskip = 0.0d0 
         if (ocvged) dskip = 1.0d0  
         label="brdcst dskip"
         call pg_brdcst(7128,  dskip,  8,  0)
         label=" "
      endif
_ENDIF
      call end_time_period(TP_TEST4)
      if (ocvged) goto 321
      rshift=0.0d0
      if(.not.ondiis.or.olevd)
     *  call shiftq(q(i50),nocmx,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)

c for screened I/O check if we need broadcast on slave nodes
      if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &    call pg_brdcst(7127,  q(i50),  l1*(l1+1)/2*8,  0)

      m2=2
      diaacc = diff*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
      call start_time_period(TP_DIAG)
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
      call end_time_period(TP_DIAG)

      call start_time_period(TP_TEST5)
      if(oshift) then
        call rdedx(q(i10),l3,ibl3qa,idaf)
      else
        call rdedx(q(i10),l3,ibl3qs,idaf)
      endif

_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
      ig = mod(iter,igs)

      call end_time_period(TP_TEST5)
      if (ig .eq. 0) then
         call start_time_period(TP_ORFOG)
         call rdedx(q(i10),l2,ibl7st,num8) 
_IF(ga)
         if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
            call porth(q(i30),q(i10), num, l0)
         else
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1) 
         endif
_ELSE
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1) 
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1) 
_ENDIF
         call end_time_period(TP_ORFOG)
      endif

      call start_time_period(TP_TEST6)

      call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      ig=nocmx
      ndoc = na
      if (osmear) then 
         esmear = max(esmear_final,
     &                min(esmear,esmear_scale*(diff-acurcy)))
         call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                    iwr,oprint(60))
         emrmn = 2*emrmn
      endif
      if(ig.lt.l0) then
_IFN1(civu)        call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civu)        call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      endif
      if (.not.osmear) then
        call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
      if (nprint .eq. 5) then
         write (iwr,9148)
         call prev(q(i30),q(i90),lprnt,l1,l1)
      endif
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
*        if(opg_root())write(6,*)'p dmat'
         call dmtxp(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      else
         call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      endif
_ELSE
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
_ENDIF
c
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      if (nprint .eq. 5) then
         write (iwr,9168)
         call prtril(q(i50),l1)
      endif
c
c     compute RMS convergence on density
c
      if (iter.ne.1) then
       call rdedx(q(i10),l2,ibl3pa,idaf)
       dsq = cvgdens(q(i50),q(i10),l2)
       dsq = dsq / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
c
c     ----- save mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i50)
c     -e- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i50),q(i90),q(i80),ndaf,l1,l2,
     *ibl3pa,ibl3ea)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)

       call end_time_period(TP_TEST6)

 321  if(numdis.ne.num8) then
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
      endif
      tim1=cpulft(1)
      delt = tim1-tim0
      tim0 = tim1
c
c     monitor static tester with diis on - at the moment remove diis
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if (oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif

      otrigp = maxcyc-iter.lt.10
      if(outon.or.otrigp) then
       if(odebug(31)) then
         write(iwr,90281)
         write (iwr,9188) iter,kcount,etot,
     +   ehf,de,diff,rshift,damp,derror,delt,tim1,yavr
       else
         write(iwr,90282)
         write (iwr,9189) iter,kcount,etot,
     +   ehf,de,diff,rshift,damp,derror,yavr
       endif
      endif
      write(text,9188) iter,kcount,etot,ehf,de,diff,rshift,damp,
     &                 derror,delt,tim1,yavr
      call sptchk(text)
_IF(parallel)
_IF(diag_parallel)
c
c broadcast of back-transformed vectors 
c probably this is a waste of time
c
c      if (ipiomode .eq. IO_NZ_S .and. .not.ocvged) then
c         label="brdcst vect"
c         call pg_brdcst(7126,q(i30),l1*l1*8,0) 
c         label=" "
c      endif
_ENDIF
      else
c
c    stuff for other screened-out nodes 
c    this is only necessary if these nodes are joining in the
c    parallel diag
c    in future extra sections may be needed to enable parallel diis with
c    ipiomode = IO_NZ_S
c
         iter = iter + 1
         if(irest.gt.0) then
            write(iwr,3010)irest
            go to 460
         endif
c
      if(odpdiag(l1)) then

c  find out if energy has already converged 
          label="brdcst dskip"
          call pg_brdcst(7128,  dskip,  8,  0) 
          label=" "
          if (dskip .lt. 0.9d0) then 

c  receive diis fock from root node, and diagonalise
             label="brdcst fock"
             call pg_brdcst(7127,  q(i50),  l1*(l1+1)/2*8,  0)
             label=" "

             call start_time_period(TP_DIAG)
             m2=2
             diaacc = diff*5.0d-3
             if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
             if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c     
             call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),
     &            m2,lock,diaacc)
             call end_time_period(TP_DIAG)
          endif 
       endif 

c gdf:   gather back-transformed vectors 
c probably wasted, as slave nodes don't neet to know the vectors?
c       label="brdcst vect"
c       call pg_brdcst(7126,q(i30),l1*l1*8,0) 
c       label=" "
c 3219  continue 
c
      endif

      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst maxcyc"
         lscf = 20+12/nav
         call pg_brdcst(7125,maxcyc,lscf*8,0)
         label=" "
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif

      call pg_synch(5555)
_ENDIF

      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
cgdf
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if 
c
c     now flag covergence if STATIC tester has been encountered
c
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
c
c  update files in core or GAs
c
      if (ioupd.gt.1) call ioupdate

_IF(parallel)
      if (.not.ocvged) then
         call timit(0)
         call timrem(tlefti)
         if (tlefti .le. skale*(tlefts-tlefti)/iter) then
          write(6,*)'** time-out on node ',ipg_nodeid()
          irest=3
         endif
         call pg_igop(333,irest,1,'+')
         if (irest.gt.0) then
           irest = 3
           nindmx = 0
           write(iwr,9367)
           call texit(0,irest)
           go to 400
         endif
         if (iter .lt. maxcyc) go to 140
         if(omaxcyc) then
           write (iwr,9289)
           ocvged = .true.
         else
           write (iwr,9288)
           etot = dzero
           irest=3
           nindmx = 0
           ehf = -en
         endif
      endif
_ELSE
      if (.not.ocvged) then
         call timit(0)
         call timrem(tlefti)
         if (tlefti .le. skale*(tlefts-tlefti)/iter) then
           irest=3
           nindmx = 0
           write(iwr,9367)
           call texit(0,irest)
         else
           if (iter .lt. maxcyc) go to 140
           if(omaxcyc) then
            write(iwr,9289)
            ocvged = .true.
           else
            write (iwr,9288)
c           write out error description to xml/punchfile            
            call blkerror('excessive number of SCF iterations',0)
            etot = dzero
            irest=3
            nindmx = 0
            ehf = -en
           endif
         endif
      endif
_ENDIF
c
  400 continue
      lock = lockt
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr('Memory failure in drhfcl!')
         endif
      endif
_ENDIF
      if (ocvged.and.outon) write (iwr,9208)
      write(iwr,9373) iter,tim1,charwall()
_IF(ccpdft)
      if (odft) then
         call CD_get_dftresults(npts,aelec,belec)
         write(iwr,9374) npts,aelec,(aelec-ne)/ne,edft
      endif
_ENDIF
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if (omodel) then
       write (iwr,9371) ehf,en,etot,diffdens
      else if (ocryst) then
       write (iwr,9372) ecre,ecrn,ecre+ecrn,ehf,en,etot,
     +                  diffdens
      else if (ovir) then
         write (iwr,9365) ek,ehf,en,etot,vir,diffdens
      else
       write (iwr,9368) ehf,en,etot,diffdens
      endif
c     write out energies to xml file      
      call blkscfe(ehf,en,etot)
      if(outon.and.yavr.ne.yblnk)write(iwr,9369)
      if(numdis.eq.num8)numdis=0
c
_IF(drf)
       if(.not.oreact) then
_ENDIF
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       idum = min(iter,200)
       write(iwr,9199)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       else
       write(iwr,9198) loop, tester(loop), idomo(loop), ivmo(loop),
     +                       testerd(loop),idomod(loop),ivmod(loop)
       endif
       enddo
       write(iwr,9197)
       endif
_IF(drf)
       endif
_ENDIF
c
_IF(nbo)
c
c save a copy of the fock matrix on the dumpfile
c (temporarily use section 242)
c
      len42 = lensec(nx)
      m = 42
      call secput(isect(42),m,len42,iblk42)
      call rdedx(q(i10),l2,iblkh,num8)
      call wrt3(q(i10),nx,iblk42,idaf)
      lds(isect(42)) = nx
_ENDIF
      if (nprint .eq. 5 .or. nprint .eq. -5) go to 440
_IF(parallel)
      if(omaster) then
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
_ENDIF
       call rdedx(q(i10),l3,ibl3qa,idaf)
       call rdedx(q(i30),l2,ibl3pa,idaf)
       call rdedx(q(i21),l1,ibl3ea,idaf)
_IF(drf)
cahv ---  DRF extension ---
      if (field .ne. ' ') then
        ehfhh = ehf
        enhh = en
        etothh = etot
        call dawrit(idafh,ioda,q(i30),l2,16,navh)
      endif
cahv
_ENDIF
       call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
       call tdown(q(i10),ilifq,q(i10),ilifq,l0)
       if(otsym) call symass(q(i10),q(i21),q(i80),q)
       if(.not.oprint(25)) then
        write (iwr,9148)
        call prev(q(i10),q(i21),lprnt,l1,l1)
        if(nprint.eq.-5) then
          write (iwr,9168)
          call prtril(q(i30),l1)
        endif
       endif
_IF(parallel)
      else
c...  other roots
c...  allow for section created in symass
c...  save results of the symmetry assignment into ed3
c...  assuming this has been called
       if (otsym) then
        isymsc = isect(499)
        call secput(isymsc,isymtp,1+lensec((num-1)/nav+1)
     +                             +lensec(num),iblnum)
       endif
      endif
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
  440 continue
_IF(ga)
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.false.)
      endif
_ENDIF

_ELSE
  440 continue
_ENDIF
      if (mpunch .ne. 0) then
         call rdedx(q(i10),l3,ibl3qa,idaf)
         call pusql(q(i10),na,l1,l1)
      endif
  460 continue
      if (ocvged) irest = 0
      if(outon)then
         cpu=cpulft(1)
         write(iwr,7777)cpu,charwall()
      endif
_IF(drf)
cahv  --- DRF extension ---
      if (odrf .or. oarf) then
        if (oarf) then
          call gmem_free_inf(i196,fnm,snm,'sf')
          call gmem_free_inf(i193,fnm,snm,'indmom')
        endif
        call gmem_free_inf(i170,fnm,snm,'iexp')
        call gmem_free_inf(i210,fnm,snm,'ddrf')
        call gmem_free_inf(i200,fnm,snm,'fdrf')
        call gmem_free_inf(i190,fnm,snm,'ijbit_op')
        call gmem_free_inf(i180,fnm,snm,'ijbit_s')
        call gmem_free_inf(i160,fnm,snm,'omega_op')
        call gmem_free_inf(i150,fnm,snm,'omega_s')
        call gmem_free_inf(i140,fnm,snm,'dipz')
        call gmem_free_inf(i130,fnm,snm,'dipy')
        call gmem_free_inf(i120,fnm,snm,'dipx')
        call gmem_free_inf(i110,fnm,snm,'ovlap')
      endif 
_ENDIF
c
c  update files in core or GAs
c
      if (ioupd.eq.1) call ioupdate

      call gmem_free_inf(irdmat,fnm,snm,'rdmat')

c ***
c *** set tolitr so that if reenter in geom opt with old vecs
c *** start off with full accuracy
c ***
      tolitr(1)=tolitr(3)
      tolitr(2)=tolitr(3)

      accdi1 = accdin
      call copq_again(mouta)
c
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapa2 = gapa1
c      maxcyc = maxcyc / 2
c     endif
      return
c
c9831 format(/15x,'**** increase level shifter to ',f8.2,' ****'/)
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9827 format(1x,'rising diff - de = ', f20.10)
 9825 format(/' *** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffd, diffdp = ',2f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9198 format(10x,i5,f11.7,' (',2i4,') (*)',6x,f11.7,' (',2i4,')' )
 9298 format(10x,i5,f11.7,' (',2i4,')',10x,f11.7,' (',2i4,') (*)' )
 9199 format(/10x,64('+')/
     + 10x,'CONVERGENCE / TESTER ANALYSIS',12x,'From DIIS'/
     + 10x,64('+')/
     + 10x, 'Iter.', '   Tester   (domo/vmo)',
     + 10x,          '   Tester   (domo/vmo)'/)
 9197 format(10x,64('+'))
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 2970 format(' restart parameter in direct-scf ',i2)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))
 3010 format(' restart parameter after dhstart ',i3)
 7777 format(//
     *' end of closed shell scf at ',f12.2,
     *' seconds',a10,' wall',//1x,104('-')//)
 9369 format(//20x,90('*')/
     *20x,'*',88x,'*'/20x,'*',88x,'*'/
     *20x,'*',31x,'warning  state is averaged',31x,'*'/
     *20x,'*',88x,'*'/20x,90('*')//)
 2990 format(/1x,'restarting direct-scf in integrals :',
     *  'restart parameter =',i3/
     *  1x,'section 171 at block ',i5)
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9008 format(/' ----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9028 format(//15x,'convergence data'/15x,16('=')///
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criterion        =1.0e-',i2/
     +     ' punch out option             = ',i6/
     +     ' integral prefactor tolerance = ',e10.2/)
90281 format(/1x,126('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis',
     * 4x,'del(t)',6x,'time'/
     * 18x,'energy',10x,'energy',36x,'shift'/1x,126('='))
90282 format(/1x,103('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     * 18x,'energy',10x,'energy',36x,'shift'/1x,103('='))
 9188 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,
     1       2f10.3,a3)
 9189 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,a3)
 9048 format(20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/ 20x,23(
     +     '-'))
 9148 format(///1x,104('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9168 format(/10x,14('-')/10x,'density matrix'/10x,14('-'))
 9208 format(/10x,16('-')/10x,'energy converged'/10x,16('-'))
 9228 format(/10x,27('-')/10x,'lagrange multipliers matrix'/10x,
     +     27('-'))
 9288 format(/10x,30('-')/10x,'excessive number of iterations'/
     +     10x,30('-'))
 9289 format(/10x,30('-')/
     +        10x,'excessive number of iterations'/
     +        10x,'but flag SCF convergence'/
     +        10x,30('-'))
 9308 format(1x,'core assignment'/' i10, i20, i21,',
     +     ' i30, i31, i40, i50, i60, i80, i90  = '/ 10i8/
     +     1x,'last = ', i8)
 9328 format(/' ehf1 = ',f20.12,' ehf2 = ',f20.12,' ehf = ',f20.12)
 9348 format(//36x,39(1h*)/
     + 36x,'direct closed-shell rhf scf calculation',
     + ' - v1.1',
     +       /36x,39(1h*))
 9365 format(
     +10x,'kinetic energy             ',f20.10/
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy               ',f20.10,8x,f10.7/
     +10x,'convergence on density     ',f20.10)
 9368 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9371 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total qm energy            ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9372 format(
     +10x,'electron-crystal energy    ',f20.10/
     +10x,'nuclear-crystal energy     ',f20.10/
     +10x,'total-crystal energy       ',f20.10/
     +10x,'total electronic energy    ',f20.10/
     +10x,'total nuclear energy       ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9373 format(/10x,14('-')/
     +10x,'final energies   after',i4,' cycles at ',f12.2,' seconds',
     +a10,' wall',/10x,14('-')//)
 9374 format(
     +10x,'number of quadrature points',i9/
     +10x,'integrated electron count  ',f20.10,5x,
     +    'relative error ',e10.2/
     +10x,'XC energy                  ',f20.10)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'==============================')
      end
_IF(ga)
      subroutine drhfcl_ga(q,coord,atmchg)
c
c  parallel direct SCF 
c
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  not supported in ga-based code
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c  
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes 
c           run through all code section
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      character*127 text
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/atmol3)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/scra7)
INCLUDE(common/disc)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
INCLUDE(common/timeperiods)
INCLUDE(common/runlab)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/xfield)
INCLUDE(common/psscrf)
INCLUDE(common/harmon)
INCLUDE(common/fermidirac)
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     +             iposit(20),nsti(2),ondiis,junkj
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      dimension q(*),coord(3,*),atmchg(*)
c     dimension osign(2)
INCLUDE(common/gadiis)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
c ... for dummy symass section
      parameter(isymtp=99)
      character*20 label
      common/msglab/label
c
       integer domo, vmo
       logical logtest
       common/testerp/tester(200), idomo(200), ivmo(200),
     +               testerd(200),idomod(200), ivmod(200),
     +               testdum(2,200),idomodum(2,200),
     +               logtest(200)
c
INCLUDE(common/drfopt)
INCLUDE(common/zorac)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/statis)
      character*10 charwall
      character*1 xn,xt
      Logical use_symmetry, force_serial
      Integer n_irreps
      character*5 fnm
      character*9 snm
      data fnm,snm/'scf.m','drhfcl_ga'/
c
      data xn,xt/'n','t'/
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data bigd/1.0d4/
c
      call check_feature('drhfcl_ga')
c
*IJB - timer JUST for the scf
*
      Call start_time_period( TP_DRHFCL_GA )
*
      idum = CD_update_geom(c)
      odft = CD_active()
      dft_accu = 1.0d0
c
c     variables for monitoring tester + convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmn = 0.0d0
      esmear = esmear_start
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      out = nprint .eq. 5
      outon = nprint .ne. -5
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)

      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
c
c     is delta density ever to be invoked? if not, we can
c     reduce i/o activity around section 471 ..
c     this is controlled by the variable odelta
c     now set in routine start2
c     odelta = deltol .gt. dlog(acurcy)
c
c
c - this should be default anyway -
c
      oswed3(4) = .true.
      oswed3(8) = .true.
c
c master nodes execute whole code
c
      omaster = .true.
c
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      lentri = l2
      l3 = num*num
c     l4 = num * 9
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
c
      ilen=0
      if(ocryst) then
       call caserr('ocryst is true')
       ilen=l2
      endif
c     ilen = ilen + ikyp(nshell)*2 + 6 * l2 + 2 * l1
      ilen = ilen + ikyp(nshell)*2 + 4 * l2 + 2 * l1

      irdmat = igmem_alloc_inf(ilen,fnm,snm,'irdmat',IGMEM_NORMAL)

c      if(opg_root())write(6,*)'SCF start address',irdmat, 
c     &  ' len ',ilen

      iprefa=irdmat+ikyp(nshell)
      i10=iprefa+ikyp(nshell)
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
c     i50 = i40+l2
      i50 = i20+l2
      i60 = i50+l2
      i61 = i50+l3
      i80 = i60+l2
      i90 = i80+l1
      ix  = i90+l1
      if(ocryst) then
       last = ix + l2
      else
       last = ix
      endif

      length = last-irdmat
c
c consistency check
      if (length .ne. ilen)then
         write(iwr,*)length,last,irdmat,ilen
         call caserr('size error')
      endif

      if (nprint .eq. 5) write (iwr,9308) i10,i20,i21,i30,i31,i40,i50,
     +     i60,i80,i90,last
c ***
      dlnmxd=0.0d0
      dlntol=tolitr(3)
c     if(opg_root())write(6,*)' ipiomode,IO_NZ_S = ', 
c    +              ipiomode,IO_NZ_S
**** GA
c
c  check/create GA storage, load S
c
      call rdedx(q(i20),l2,ibl7s,num8)
      call rdedx(q(i10),l2,ibl7st,num8)
      call declare_diis_storage(num,.false.)
      call init_diis_storage(num,q(i20),q(i10))

      call start_time_period(TP_RDMAT)
c
c     specification of section isect(471)
c     only output delta-difference data if is this is to be invoked
c     i.e. only rdmat and prefac traffic involved in default
c *** section isect(471) now holds:
c *** prefac mat., reduced dmat, delta f, delta d, last f, current d
c
      m171t=171
      len171=lensec(l2)
      nshtri=ikyp(nshell)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      do i=4,6
       iof171(i)=iof171(i-1)+len171
      enddo
      if(irest.ne.1.and.irest.ne.2) then
        call rdmake(q(iprefa))
        ln171=nshblk*2
        if(odelta) then
         ln171=4*len171+ln171
        endif
        call secput(isect(471),m171t,ln171,ibl171)
        if(outon) then
         write(iwr,2980) ibl171,ln171
        endif
        call wrt3(q(iprefa),nshtri,ibl171,idaf)
        if(odelta) then
         call zer171(q(i50),l2,4,ibl171+iof171(3),len171,idaf)
        endif
        call clredx
      endif
c
      call end_time_period(TP_RDMAT)
c ***
      call start_time_period(TP_TEST1)
      if (maxcyc .le. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
       tester(loop) = 0.0d0
       testerd(loop) = 0.0d0
       logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1

      i=lensec(l2)
      iblkh0 =ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)
      if(numdis.eq.0)then
         numdis=num8
         ndafd=iblkh+i
      else
         ndafd=ibldis+2
      endif
      ndafdi=ndafd+3*i
      if(irest.eq.0.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         kcount = 0
         iterv = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      diffdp = bigd
      diffdpp = bigd
c
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (odamph .or. oshift) damp = done
      lockt=lock
      tim0=cpulft(1)
      ovir = .false.
      if(oaimpac())ovir = .true.
c
c     ----- nuclear energy
c     ----- psscrf dvsg
c
      if (opssc) then
         en = denuc(nat,atmchg,coord,itotbq,iwr)
      else if (ocryst) then
         en = crnuc(nat,czan,coord,ecrn)
      else
         en = enucf(nat,atmchg,coord)
      endif
      l0=newbas0
      call end_time_period(TP_TEST1)
       n_irreps = 0
       Do i = 1, l1
          isymmo( i ) = isymao( i )
          n_irreps = Max( n_irreps, isymao( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
       if(outon) then
        write (iwr,9008) en
        if (use_symmetry) write(iwr,9009)
        write (iwr,9028) maxcyc,mconv,nconv,mpunch,
     1  dexp(-dlntol)
       endif
       call start_time_period(TP_ORFOG)
*      call rdedx(q(i10),l2,ibl7st,num8)
       call qmat_symm(q,q(i10),q(i50),q(i61),q(i20),iky,l0,l1,l3,l1,out,
     +               isymmo, use_symmetry )
c load guess vectors..
       call rdedx(q(i10),l3,ibl3qa,idaf)
c load overlap
       call rdedx(q(i50),l2,ibl7st,num8)
**** GA
       call porth(q(i10),q(i50), num, l0)
*****
       call end_time_period(TP_ORFOG)

       call start_time_period(TP_TEST2)
       call wrt3(q(i10),l3,ibl3qa,idaf)
       call tdown(q(i10),ilifq,q(i10),ilifq,l1)
       call rdedx(q(i90),l1,ibl3ea,idaf)
       ndoc = na
       if (osmear) then
         call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                    iwr,oprint(60))
         emrmn = 2*emrmn
       else
         call llvmo(q(i90),q(i80),na,nocmx,l1)
       endif
       call dscal(l1,two,q(i80),1)
c...
c...  construct initial density matrix
c...
***** GA
       call dmtxp(q(i50),q(i10),q(i80),iky,nocmx,l1,l1)
       call wrt3(q(i50),l2,ibl3pa,idaf)
       call end_time_period(TP_TEST2)
*****
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      lprnt = l1
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      call timrem(tlefts)
c *** build rest of fock matrix from 2-e integrals
      ifock=i10
      idmat=i50
      if(irest.eq.1 .or. irest.eq.2) then
        call secget(isect(471),m171t,ibl171)
        write(iwr,2990)irest,ibl171
        call rdedx(q(iprefa),nshtri,ibl171,idaf)
        call reads(q(irdmat),nshtri,idaf)
        if(odelta) then
         call rdedx(q(i10),l2,ibl171+iof171(3),idaf)
         call reads(q(i50),l2,idaf)
        endif
        dlnmxd=-9999999.0d0
        do 1672 kkk=1,nshtri
            if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
            dlnmxd=q(irdmat+kkk-1)
1672    continue

      endif
c ***
      if (CD_active()) then
         call retrieve_spare(imemspare)
         idum = CD_set_2e()
         if (CD_2e() .and. (.not. CD_HF_exchange())
     &      .and. (.not. CD_HF_coulomb()) )then
            imemdhstaru = 0
         else
            call mem_dhstaru(imemdhstaru)
         endif
         idum = CD_reset_2e()
         imemfree = igmem_max_memory() - imemspare - imemdhstaru
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
c        if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q,q,iwr)
c        else if (ks_bas.eq.KS_MO) then
c           imemreq = CD_memreq_energy_mo(l1,na,0,q,q,iwr)
c        else if (ks_bas.eq.KS_AOMO) then
c           imemreq = CD_memreq_energy(l1,na,0,q,q,iwr)
c        endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then 
            write(iwr,600)ierror
            call caserr('Out of memory in incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
      endif
*****
 140  continue
c
c  intermediate times
c      call list_time_periods(.false.,.false.)

c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** delta f, delta d, last f, current d, reduced dmat, prefac mat.
c *** note that delta f and delta d are NOT accessed if
c *** delta-SCF has not been requested : odelta = .false.)
c *** section now contains:
c *** prefac mat., reduced dmat, delta f, delta d, last f, current d
      call start_time_period(TP_RDMAT)
      if(irest.eq.1) goto 6351
*     if(irest.eq.2) goto 6352
      if (odelta) then
         call vclr(q(i60),1,l2)
         call wrt3(q(i60),l2,ibl171+iof171(3),idaf)
         call rdedx(q(i60),l2,ibl171+iof171(6),idaf)
         call wrt3(q(i50),l2,ibl171+iof171(6),idaf)
         call vsub(q(i50),1,q(i60),1,q(i50),1,l2)
         call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
c     else
c        call wrt3(q(i50),l2,ibl171+iof171(6),idaf)
      endif
      if(outon) write(iwr,3030)
      if (odelta) then
       call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
      endif
c ***
c *** make sure that iter before delta calc is of full accuracy.
c ***
      if( (dlnmxd.lt.deltol) .and. (iter-1.lt.itrtol(2))) then
        dlnmxd=deltol+1.0d0
      endif
      if( (dlnmxd.gt.deltol) ) then
        if(outon) write(iwr,3000)
        if (odelta) then
           call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
           call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
           call vclr(q(i60),1,l2)
           call wrt3(q(i60),l2,ibl171+iof171(5),idaf)
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        else
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        endif
        call clredx
      endif
      call wrt3(q(irdmat),nshtri,ibl171+iof171(2),idaf)
6351  continue
c ***

      call end_time_period(TP_RDMAT)
      if(iter.lt.itrtol(1)) then
         dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
         dlntol=tolitr(2)
      else
         dlntol=tolitr(3)
      endif
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(outon) then
         write(iwr,3020) dlntol
      endif
c
c skip dhstar when not required for choice of DFT functional
c
        idum = CD_set_2e()

        if(CD_2e() .and. (.not. CD_HF_exchange())
     &      .and. (.not. CD_HF_coulomb()) )then

          o2e = .false.
          if(irest.ne.1)call vclr(q(ifock),1,l2)
        else
          o2e = .true.
          call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),irest)
        endif

      idum = CD_reset_2e()

       call start_time_period(TP_TEST3)
        if(irest.gt.0) then
         write(iwr,3010)irest
         if(numdis.eq.num8) go to 460
         call wrt3(en,nw1d,ibldis,numdis)
         call wrt3s(st,nw2d,numdis)
         go to 460
c *** form complete 2-e skelton fmat and save it
        endif
       if (odelta) then
         call rdedx(q(i60),l2,ibl171+iof171(5),idaf)
         call vadd(q(i10),1,q(i60),1,q(i10),1,l2)
         call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
         call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
c      else
c        call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
       endif
c ***
       if (nprint .eq. 5) then
          write (iwr,9068)
          call prtril(q(i10),l1)
       endif

c  symmetrise skeleton fock matrix 
c  if only 2e integrals were performed
       if(o2e)then
        call symh(q(i10),q(i20),iky,0,0)
        if (nprint .eq. 5) then
           write (iwr,9088)
           call prtril(q(i10),l1)
        endif
       endif
c
c increment H_2 at q(i10) with H_1
c
      call rdedx(q(i20),l2,ibl7f,num8)
      if (ozora) call zora(q,q,q(i20),'read')
      call adonee(q(i10),q(i20),q(iprefa))
c
c  ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i50),q(i20),l1)
c
c ps moved this up so as not to overlap with dft timings
c disabled for integral restarts
c
      if(irest.ne.2)call end_time_period(TP_TEST3)
c
c  save previous total 
c
      ehf0 = ehf

      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(o2e)then
c
c Coulomb operator is in q(i10) augmented by H_1, 
c compute energy using HF expression without exchange
c
            ehf2 = tracep(q(i50),q(i10),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1
         endif
c
c Update Kohn-Sham matrix and compute fitted/integrated 
c energy terms
c
c        dft_accu = 0.95d0*dft_accu
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
         idum = CD_energy_ao(c,q(i10),dum,q(i50),dum,
     +        edft,q,q,outon,dft_accu,iwr
     +        )
         call symm_op(q(i10))
         ehf = etmp+edft

         if(out) then
           if(opg_root()) then
            call CD_print_dftresults(.true.,.false.,iwr)
            write(iwr,9407)ehf1, ehf2
           endif
         endif
         call timana(31)
         call cpuwal(begin,ebegin)

      else
c
c  Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
         ehf2 = tracep(q(i50),q(i10),l1)
         ehf = (ehf1+ehf2)*pt5
      endif

      ehf = ehf + 0.5d0*emrmn

c for xtal field calculations, the core hamiltonian
c includes the crystal potential, but for the energy 
c all terms modelling molecule-molecule interactions 
c must be divided by two 

      if(ocryst)then
         m53=53   !!!!!!!!
         write(iwr,*)'secget',isecfx,m53
         call secget(isecfx,m53,iblk)
         call rdedx(q(ix),l2,iblk,idaf)
         ehf3 = tracep(q(i50),q(ix),l1)
         write(iwr,*)'ehf3',ehf3
         esave = ehf
         ehf = (ehf1+ehf2-ehf3)*pt5 
        ecre = esave - ehf
         write(iwr,*)'new ehf, ecre',ehf,ecre
      endif
c
c compute kinetic energy and virial if required
c
      if(ovir)call rhfvir(q(i50),ehf,en,ek,vir,num8,l1,l2)

      call wrt3(q(i10),l2,iblkh,num8)
c hack remco
c        rewind(76)
c         print *,'remco l2',l2
c        write(76)(q(i10+i-1),i=1,l2)
c einde hack
      if (nprint .eq. 5) then
        write (iwr,9048)
        call prtril(q(i10),l1)
        call rdedx(q(i90),l1,ibl3ea,idaf)
        iorb1 = 1
        iorb2 = nocmx
        call makfv(q(i10),q(i30),q(i30),q(i21),iky,iorb1,iorb2,l1,l1)
        call rdedx(q(i10),l3,ibl3qa,idaf)
        call makeij(q(i10),q(i30),q(i30),q(i21),l0,l0,l1,l0,l1)
        write (iwr,9228)
        call eijout(q(i30),q(i90),l0,nocmx,l1,iwr)
        write (iwr,9328) ehf1,ehf2,ehf
        call rdedx(q(i10),l2,iblkh,num8)
      endif
      iter = iter+1
      etot0 = etot
      etot = ehf+en
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
      if (damp .lt. dmpcut) damp = dmpcut
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i50),q(i60),l1,l2,ndafd,
     +numdis,iterv,1,1)
      diffpp=diffp
      diffp=diff

      if(odiis) then
      call start_time_period(TP_DIIS)
***** GA
      call diis_ga(q(i10),q(i20),q(i50),ibl3qa,diffd,domo,vmo)
*****
        if (iter.le.maxcycp) then
         idomod(iter) = domo
         ivmod(iter) =  vmo
         testerd(iter) = diffd
        endif
      call end_time_period(TP_DIIS)
      endif

      call start_time_period(TP_TEST4)

c
c ------ take precautions if tester is increasing again
c
      if(iter.eq.1 .or. .not.odiis) goto 223
      if(ondiis) then
        if (oprind) then
         write(iwr,9824)diffd,diffdp
         write(iwr,9823)diffd, acurcy*diisf
         write(iwr,*) 'de = ', de
        endif
        if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.acurcy*diisf
     1  .or.(otestd.and.de.lt.1.0d-10)) then
*        if (de.lt.1.0d-5) then
          if(iter.le.maxcycp) logtest(iter) = .true.
          diffdpp = diffdp
          diffdp = diffd
          go to 223
*        else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
*         odynamic = .true.
*         accdi1 = 0.0d0
*         itdiis = 0
*         if(oprind) write(iwr,9828) iter
*        endif
        else
          if(oprind) write(iwr,9829) iter
        endif
        nsti(1)=0
        ondiis=.false.
        call rdedx(q(i10),l2,iblkh,num8)
        diffdp = bigd
        diffdpp = bigd
      else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c       if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c        iwild = iwild + 1
c        osign(iwild) = de.lt.0.0d0
c        if (iwild.eq.2.and..not.oreset.and.
c    +       .not.(osign(1).and.osign(2)) ) then
c         gapa1 = 5.0d0
c         gapa2 = gapa1
c         oreset = .true.
c         iwild = 0
c         write(iwr,9831) gapa1
c         maxcyc = maxcyc + maxcyc
c         maxcycp = min(200,maxcyc)
c        endif
c       endif
c
        if (odynamic) then
         if(diff.le.diffp.or.diffp.le.diffpp) then
          itdiis = itdiis + 1
          if (itdiis.eq.3) then
            accdi1 = diff
            if (oprind) write(iwr,9826) accdi1
          endif
         else
          if (oprind) write(iwr,9827)
          itdiis = 0
         endif
        endif
      endif
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -q- at q(i50) orthonormalizing transformation
c     -h- at q(i10)
c    -h'- at q(i50) transformed h matrix
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
 223  if(oshift)then
         call rdedx(q(i50),l3,ibl3qa,idaf)
      else
         call rdedx(q(i50),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i50 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i50),ilifq,q(i50),ilifq,l0)
***** GA
      call load_ga_from_square(ih_vec,q(i50),l1)
      call load_ga_from_triangle(ih_scr2,q(i10),l1)
      call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
      call load_triangle_from_ga(q(i10),ih_scr2, l1)
*****
      diff=0.0d0
      ii=ndoc+1
      if(ii.le.l0)then
         do 250 i=ii,l0
            ij=iky(i)+i10
            loop=idamax(min(i-1,nocmx),q(ij),1)
            if(loop.gt.0) then
c            diff=dmax1(diff,dabs(q(ij+loop-1)))
             dum = dabs(q(ij+loop-1))
             if (dum.gt.diff) then
               diff = dum
               if(iter.le.maxcycp) then
                tester(iter) = diff
                idomo(iter) = loop
                ivmo(iter) = i
               endif
             endif
            endif
 250     continue
      endif
      if(ondiis.and..not.osmear) diff=diffd
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
c
cgdf  not sure if needed: convergence check here seems to be for 
c     skipping some work, when the real check comes later
c
      if (osrso) ocvged = .true.
c
      call end_time_period(TP_TEST4)
      if (ocvged) goto 321
      rshift=0.0d0
      if(.not.ondiis.or.olevd)
     *  call shiftq(q(i10),nocmx,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)
      m2=2
      diaacc = diff*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
      call start_time_period(TP_DIAG)
      call jacobi_symm(q(i10),iky,l0,q(i50),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
      call end_time_period(TP_DIAG)

      call start_time_period(TP_TEST5)
      if(oshift) then
        call rdedx(q(i10),l3,ibl3qa,idaf)
      else
        call rdedx(q(i10),l3,ibl3qs,idaf)
      endif

***** GA
         call load_ga_from_square(ih_scr,q(i50),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)
         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)
         call load_square_from_ga(q(i50),ih_scr2, l1)
***** GA
      ig = mod(iter,igs)

      call end_time_period(TP_TEST5)
      if (ig .eq. 0) then
         call start_time_period(TP_ORFOG)
         call rdedx(q(i10),l2,ibl7st,num8) 
***** GA
         call porth(q(i50),q(i10), num, l0)
*****
         call end_time_period(TP_ORFOG)
      endif

      call start_time_period(TP_TEST6)

      call wrt3(q(i50),l3,iblkqq,num8)
      call tdown(q(i10),ilifq,q(i50),ilifq,l0)
      ig=nocmx
      ndoc = na
      if (osmear) then
         esmear = max(esmear_final,
     &                min(esmear,esmear_scale*(diff-acurcy)))
         call fermi_smear(q(i90),q(i80),emrmn,ndoc,nocmx,esmear,na,l0,
     &                    iwr,oprint(60))
         emrmn = 2*emrmn
      endif
      if(ig.lt.l0) then
        call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
      endif
      if (.not.osmear) then
        call llvmo(q(i90),q(i80),na,nocmx,l1)
      endif
      call dscal(l1,two,q(i80),1)
      if (nprint .eq. 5) then
         write (iwr,9148)
         call prev(q(i10),q(i90),lprnt,l1,l1)
      endif
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
***** GA
      call dmtxp(q(i50),q(i10),q(i80),iky,nocmx,l1,l1)
c
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      if (nprint .eq. 5) then
         write (iwr,9168)
         call prtril(q(i50),l1)
      endif
c
c     compute RMS convergence on density
c
      if (iter.ne.1) then
       call rdedx(q(i60),l2,ibl3pa,idaf)
       dsq = cvgdens(q(i50),q(i60),l2)
       dsq = dsq / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
c
c     ----- save mo*s + density + orbital energies -----
c
c     -v- at q(i10)
c     -d- at q(i50)
c     -e- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i10),l3,iblkqq,num8)
      call scfsav(q(i10),q(i50),q(i90),q(i80),ndaf,l1,l2,
     *ibl3pa,ibl3ea)

       call end_time_period(TP_TEST6)

 321  if(numdis.ne.num8) then
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
      endif
      tim1=cpulft(1)
      delt = tim1-tim0
      tim0 = tim1
c
c     monitor static tester with diis on - at the moment remove diis
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if (oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif

      otrigp = maxcyc-iter.lt.10
      if(outon.or.otrigp) then
       if(odebug(31)) then
         write(iwr,90281)
         write (iwr,9188) iter,kcount,etot,
     +   ehf,de,diff,rshift,damp,derror,delt,tim1,yavr
       else
         write(iwr,90282)
         write (iwr,9189) iter,kcount,etot,
     +   ehf,de,diff,rshift,damp,derror,yavr
       endif
      endif
      write(text,9188) iter,kcount,etot,ehf,de,diff,rshift,damp,
     &                 derror,delt,tim1,yavr
      call sptchk(text)
c
      call pg_synch(5555)

      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
cgdf
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if 
c
c     now flag covergence if STATIC tester has been encountered
c
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
c
c  update files in core or GAs
c
      if (ioupd.gt.1) call ioupdate

      if (.not.ocvged) then
         call timit(0)
         call timrem(tlefti)
         if (tlefti .le. skale*(tlefts-tlefti)/iter) then
          write(6,*)'** time-out on node ',ipg_nodeid()
          irest=3
         endif
         call pg_igop(333,irest,1,'+')
         if (irest.gt.0) then
           irest = 3
           nindmx = 0
           write(iwr,9367)
           call texit(0,irest)
           go to 400
         endif
         if (iter .lt. maxcyc) go to 140
         if(omaxcyc) then
           write (iwr,9289)
           ocvged = .true.
         else
           write (iwr,9288)
c          write out error description to xml/punchfile            
           call blkerror('excessive number of SCF iterations',0)
           etot = dzero
           irest=3
           nindmx = 0
           ehf = -en
         endif
      endif
c
  400 continue
      lock = lockt
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr('Memory failure in drhfcl_ga!')
         endif
      endif
      if (ocvged.and.outon) write (iwr,9208)
      write(iwr,9373) iter,tim1,charwall()
      if (odft) then
         call CD_get_dftresults(npts,aelec,belec)
         write(iwr,9374) npts,aelec,(aelec-ne)/ne,edft
      endif
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if (omodel) then
       write (iwr,9371) ehf,en,etot,diffdens
      else if (ocryst) then
       write (iwr,9372) ecre,ecrn,ecre+ecrn,ehf,en,etot,
     +                  diffdens
      else if (ovir) then
         write (iwr,9365) ek,ehf,en,etot,vir,diffdens
      else
       write (iwr,9368) ehf,en,etot,diffdens
      endif
c     write out energies to xml file      
      call blkscfe(ehf,en,etot)
      if(outon.and.yavr.ne.yblnk)write(iwr,9369)
      if(numdis.eq.num8)numdis=0
c
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       idum = min(iter,200)
       write(iwr,9199)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       else
       write(iwr,9198) loop, tester(loop), idomo(loop), ivmo(loop),
     +                       testerd(loop),idomod(loop),ivmod(loop)
       endif
       enddo
       write(iwr,9197)
       endif
c
      if (nprint .eq. 5 .or. nprint .eq. -5) go to 440
       call rdedx(q(i10),l3,ibl3qa,idaf)
       call rdedx(q(i50),l2,ibl3pa,idaf)
       call rdedx(q(i21),l1,ibl3ea,idaf)
       call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
       call tdown(q(i10),ilifq,q(i10),ilifq,l0)
       if(otsym) call symass(q(i10),q(i21),q(i80),q)
       if(.not.oprint(25)) then
        write (iwr,9148)
        call prev(q(i10),q(i21),lprnt,l1,l1)
        if(nprint.eq.-5) then
          write (iwr,9168)
          call prtril(q(i50),l1)
        endif
       endif
  440 continue
***** GA
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.false.)
      endif
***** GA
      if (mpunch .ne. 0) then
         call rdedx(q(i10),l3,ibl3qa,idaf)
         call pusql(q(i10),na,l1,l1)
      endif
  460 continue
      if (ocvged) irest = 0
      if(outon)then
         cpu=cpulft(1)
         write(iwr,7777)cpu,charwall()
      endif
c
c  update files in core or GAs
c
      if (ioupd.eq.1) call ioupdate

      call gmem_free_inf(irdmat,fnm,snm,'irdmat')

c ***
c *** set tolitr so that if reenter in geom opt with old vecs
c *** start off with full accuracy
c ***
      tolitr(1)=tolitr(3)
      tolitr(2)=tolitr(3)

      accdi1 = accdin
      call copq_again(mouta)
c
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapa2 = gapa1
c      maxcyc = maxcyc / 2
c     endif
      Call end_time_period( TP_DRHFCL_GA )
      return
c
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9827 format(1x,'rising diff - de = ', f20.10)
 9825 format(/' *** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffd, diffdp = ',2f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9198 format(10x,i5,f11.7,' (',2i4,') (*)',6x,f11.7,' (',2i4,')' )
 9298 format(10x,i5,f11.7,' (',2i4,')',10x,f11.7,' (',2i4,') (*)' )
 9199 format(/10x,64('+')/
     + 10x,'CONVERGENCE / TESTER ANALYSIS',12x,'From DIIS'/
     + 10x,64('+')/
     + 10x, 'Iter.', '   Tester   (domo/vmo)',
     + 10x,          '   Tester   (domo/vmo)'/)
 9197 format(10x,64('+'))
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 2970 format(' restart parameter in direct-scf ',i2)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))
 3010 format(' restart parameter after dhstart ',i3)
 7777 format(//
     *' end of closed shell scf at ',f12.2,
     *' seconds',a10,' wall'//1x,104('-')//)
 9369 format(//20x,90('*')/
     *20x,'*',88x,'*'/20x,'*',88x,'*'/
     *20x,'*',31x,'warning  state is averaged',31x,'*'/
     *20x,'*',88x,'*'/20x,90('*')//)
 2990 format(/1x,'restarting direct-scf in integrals :',
     *  'restart parameter =',i3/
     *  1x,'section 171 at block ',i5)
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9008 format(/' ----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9028 format(//15x,'convergence data'/15x,16('=')///
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criterion        =1.0e-',i2/
     +     ' punch out option             = ',i6/
     +     ' integral prefactor tolerance = ',e10.2/)
90281 format(/1x,126('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis',
     * 4x,'del(t)',6x,'time'/
     * 18x,'energy',10x,'energy',36x,'shift'/1x,126('='))
90282 format(/1x,103('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     * 18x,'energy',10x,'energy',36x,'shift'/1x,103('='))
 9188 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,
     1       2f10.3,a3)
 9189 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,a3)
 9048 format(20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/ 20x,23(
     +     '-'))
 9148 format(///1x,104('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9168 format(/10x,14('-')/10x,'density matrix'/10x,14('-'))
 9208 format(/10x,16('-')/10x,'energy converged'/10x,16('-'))
 9228 format(/10x,27('-')/10x,'lagrange multipliers matrix'/10x,
     +     27('-'))
 9288 format(/10x,30('-')/10x,'excessive number of iterations'/
     +     10x,30('-'))
 9289 format(/10x,30('-')/
     +        10x,'excessive number of iterations'/
     +        10x,'but flag SCF convergence'/
     +        10x,30('-'))
 9308 format(1x,'core assignment'/' i10, i20, i21,',
     +     ' i30, i31, i40, i50, i60, i80, i90  = '/ 10i8/
     +     1x,'last = ', i8)
 9328 format(/' ehf1 = ',f20.12,' ehf2 = ',f20.12,' ehf = ',f20.12)
 9348 format(//36x,39(1h*)/
     + 36x,'direct closed-shell rhf scf calculation',
     + ' - v1.1',
     +       /36x,39(1h*))
 9365 format(
     +10x,'kinetic energy             ',f20.10/
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy               ',f20.10,8x,f10.7/
     +10x,'convergence on density     ',f20.10)
 9368 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9371 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total qm energy            ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9372 format(
     +10x,'electron-crystal energy    ',f20.10/
     +10x,'nuclear-crystal energy     ',f20.10/
     +10x,'total-crystal energy       ',f20.10/
     +10x,'total electronic energy    ',f20.10/
     +10x,'total nuclear energy       ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9373 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     +,f12.2,' seconds',a10,' wall',/10x,14('-')//)
 9374 format(
     +10x,'number of quadrature points',i9/
     +10x,'integrated electron count  ',f20.10,5x,
     +    'relative error ',e10.2/
     +10x,'XC energy                  ',f20.10)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'==============================')
      end
_ENDIF
_IF()
      subroutine symm_dens(d)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/atmol3)
c
c     ----- symmetrize existing density matrix
c
      if (osym_dens) then
         call symdenm(d,iky)
      endif
c
      return
      end
      subroutine symdenm(f,ia)
c
c     ----- symmetrize the density matrix
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension f(*),ia(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/infoa)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/scra/iso(1)
      common/blkin/pxyz(4),t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension mi(48),mj(48),u(35,35)
      if (nt .eq. 1) return
c
c     ----- find a block (i,j)
c
      do 520 ii = 1,nshell
       do itr = 1,nt
       ish = iso(ii+iliso(itr))
       if (ish .gt. ii) go to 520
       mi(itr) = ish
       enddo
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
c
      do 500 jj = 1,ii
      do 200 itr = 1,nt
      jsh = iso(jj+iliso(itr))
      mj(itr) = jsh
      if (jsh .gt. ii) go to 500
      ish = mi(itr)
      if (ish .ge. jsh) go to 180
      n = ish
      ish = jsh
      jsh = n
  180 if (ish .eq. ii .and. jsh .gt. jj) go to 500
  200 continue
      ljt = ktype(jj)
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
c
c     ----- apply projection operator
c
      do j =minj,maxj
       do i =mini,maxi
        u(i,j)=0.0d0
       enddo
      enddo
c
      do itr = 1,nt
       kk = mi(itr)
       mink = kmin(kk)
       maxk = kmax(kk)
       lkt = ktype(kk)
       lock = kloc(kk)-mink
       ll = mj(itr)
       minl = kmin(ll)
       maxl = kmax(ll)
       llt = ktype(ll)
       locl = kloc(ll)-minl
c
       do l = minl,maxl
        do k = mink,maxk
         lck = lock+k
         lcl = locl+l
         kl = ia(max(lck,lcl))+min(lck,lcl)
         t(k,l) = f(kl)
        enddo
       enddo
c
       call rdenr(t,mink,maxk,lkt,minl,maxl,llt,invt(itr))
c
       do l = minl,maxl
        do k = mink,maxk
         u(k,l) = u(k,l) + t(k,l)
        enddo
       enddo
c
      enddo
c
c     set correct (ii,jj) block
c
      dum = 1.0d0 / dble(nt)
      do j = minj,maxj
       do i = mini,maxi
        lci = loci + i
        lcj = locj + j
        lcij = ia(max(lci,lcj))+min(lci,lcj)
        u(i,j) = u(i,j) * dum
        f(lcij) = u(i,j)
       enddo
      enddo
c
c     set equivalent (kk,ll) blocks
c
      do itr = 1,nt
       kk = mi(itr)
       lkt=ktype(kk)
       mink=kmin(kk)
       maxk=kmax(kk)
       lock = kloc(kk)-mink
       ll = mj(itr)
       llt = ktype(ll)
       minl = kmin(ll)
       maxl = kmax(ll)
       locl = kloc(ll)-minl
c
       do j = minj, maxj
        do i = mini, maxi
         t(i,j) = u(i,j)
        enddo
       enddo
c
      call rdenr(t,mink,maxk,lkt,minl,maxl,llt,itr)
c
       do l = minl,maxl
        do k = mink,maxk
         lck = lock + k
         lcl = locl + l
         kl = ia(max(lck,lcl))+min(lck,lcl)
         f(kl) = t(k,l)
        enddo
       enddo
c
      enddo
c
  500 continue
  520 continue
c
      return
      end
      subroutine rdenr(t,mink,maxk,lkt,minl,maxl,llt,ntr)
c
      implicit REAL  (a-h,o-z)
      parameter (ndim=35)
      parameter (maxgro=48)
      common/bufb/ptr(3,3,maxgro),  dtr(6,6,maxgro),
     +            ftr(10,10,maxgro),gtr(15,15,maxgro)
      dimension t(ndim,ndim),v(ndim)
c
c     ----- right multiply  t  by  r,
c           result back in  t
c
      data zero /0.0d0/
c
c     ----- right multiply  t  by  r , result back in  t -----
c
      go to (500,400,300,200,100),llt
c
c     ----- g shell -----
c
  100 ng=15*(ntr-1)
      do k=mink,maxk
         do l=21,35
            dum=zero
            do n=21,35
               dum=dum+t(k,n)*gtr(l-20,n-20,ntr)
            enddo
            v(l)=dum
         enddo
         do l=21,35
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- f shell -----
c
  200 nf=10*(ntr-1)
      do k=mink,maxk
         do l=11,20
            dum=zero
            do n=11,20
               dum=dum+t(k,n)*ftr(l-10,n-10,ntr)
            enddo
            v(l)=dum
         enddo
         do l=11,20
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- d shell -----
c
  300 nd= 6*(ntr-1)
      do k=mink,maxk
         do l= 5,10
            dum=zero
            do n= 5,10
               dum=dum+t(k,n)*dtr(l- 4,n- 4,ntr)
            enddo
            v(l)=dum
         enddo
         do l= 5,10
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- p shell -----
c
  400 np= 3*(ntr-1)
      do k=mink,maxk
         do l= 2, 4
            dum=zero
            do n= 2, 4
               dum=dum+t(k,n)*ptr(l- 1,n- 1,ntr)
            enddo
            v(l)=dum
         enddo
         do l= 2, 4
            t(k,l)=v(l)
         enddo
      enddo
c
  500 continue
c
c     ----- left multiply  t  by r , result back in  t -----
c
      go to (1000,900,800,700,600),lkt
c
c     ----- g shell -----
c
  600 ng=15*(ntr-1)
      do l=minl,maxl
         do k=21,35
            dum=zero
            do n=21,35
               dum=dum+gtr(k-20,n-20,ntr)*t(n,l)
            enddo
            v(k)=dum
          enddo
          do k=21,35
             t(k,l)=v(k)
          enddo
      enddo
      go to 1000
c
c     ----- f shell -----
c
  700 nf=10*(ntr-1)
      do l=minl,maxl
         do k=11,20
            dum=zero
            do n=11,20
               dum=dum+ftr(k-10,n-10,ntr)*t(n,l)
            enddo
            v(k)=dum
          enddo
          do k=11,20
             t(k,l)=v(k)
          enddo
      enddo
      go to 1000
c
c     ----- d shell -----
c
  800 nd= 6*(ntr-1)
      do l=minl,maxl
         do k= 5,10
            dum=zero
            do n= 5,10
               dum=dum+dtr(k-4,n-4,ntr)*t(n,l)
            enddo
            v(k)=dum
         enddo
         do k=5,10
            t(k,l)=v(k)
         enddo
      enddo
      go to 1000
c
c     ----- p shell -----
c
  900 np= 3*(ntr-1)
      do l=minl,maxl
         do k= 2,4
            dum=zero
            do n= 2,4
               dum=dum+ptr(k-1,n-1,ntr)*t(n,l)
            enddo
            v(k)=dum
         enddo
         do k= 2, 4
            t(k,l)=v(k)
         enddo
      enddo
c
 1000 continue
      return
      end
_ENDIF
      subroutine symm_op(d)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/atmol3)
INCLUDE(common/vcore)
c
c     ----- symmetrize existing operator matrix
c
      if (osym_op) then
         need = 35 * 35
         iu = igmem_alloc(need)
         call symopm(d,Q(iu),iky)
         call gmem_free(iu)
      endif
c
      return
      end
      subroutine symopm(f,u,ia)
c
c     ----- symmetrize the operator matrix
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension f(*),ia(*),u(35,35)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/scra/iso(1)
INCLUDE(common/infoa)
      common/blkin/pxyz(4),t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension mi(48),mj(48)
      if (nt .eq. 1) return
c
c     ----- find a block (i,j)
c
      do 520 ii = 1,nshell
         do itr = 1,nt
            ish = iso(ii+iliso(itr))
            if (ish .gt. ii) go to 520
            mi(itr) = ish
         enddo
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii)-mini
c
         do 500 jj = 1,ii
            do 200 itr = 1,nt
               jsh = iso(jj+iliso(itr))
               mj(itr) = jsh
               if (jsh .gt. ii) go to 500
               ish = mi(itr)
               if (ish .ge. jsh) go to 180
               n = ish
               ish = jsh
               jsh = n
  180          if (ish .eq. ii .and. jsh .gt. jj) go to 500
  200       continue
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj)-minj
c
c     ----- apply projection operator
c
            do j =minj,maxj
               do i =mini,maxi
                  u(i,j)=0.0d0
               enddo
            enddo
c
            do itr = 1,nt
               kk = mi(itr)
               mink = kmin(kk)
               maxk = kmax(kk)
               lkt = ktype(kk)
               lock = kloc(kk)-mink
               ll = mj(itr)
               minl = kmin(ll)
               maxl = kmax(ll)
               llt = ktype(ll)
               locl = kloc(ll)-minl
c
               do l = minl,maxl
                  do k = mink,maxk
                     lck = lock+k
                     lcl = locl+l
                     kl = ia(max(lck,lcl))+min(lck,lcl)
                     t(k,l) = f(kl)
                  enddo
               enddo
c
               call ropr(t,mink,maxk,lkt,minl,maxl,llt,itr)
c
               do l = minl,maxl
                  do k = mink,maxk
                     u(k,l) = u(k,l) + t(k,l)
                  enddo
               enddo
c
            enddo
c
c     set correct (ii,jj) block
c
            dum = 1.0d0 / dble(nt)
            do j = minj,maxj
               do i = mini,maxi
                  lci = loci + i
                  lcj = locj + j
                  lcij = ia(max(lci,lcj))+min(lci,lcj)
                  u(i,j) = u(i,j) * dum
                  f(lcij) = u(i,j)
               enddo
            enddo
c
c     set equivalent (kk,ll) blocks
c
            do itr = 1,nt
               kk = mi(itr)
               lkt=ktype(kk)
               mink=kmin(kk)
               maxk=kmax(kk)
               lock = kloc(kk)-mink
               ll = mj(itr)
               llt = ktype(ll)
               minl = kmin(ll)
               maxl = kmax(ll)
               locl = kloc(ll)-minl
c
               do j = minj, maxj
                  do i = mini, maxi
                     t(i,j) = u(i,j)
                  enddo
               enddo
c
               call ropr(t,mink,maxk,lkt,minl,maxl,llt,invt(itr))
c
               do l = minl,maxl
                  do k = mink,maxk
                     lck = lock + k
                     lcl = locl + l
                     kl = ia(max(lck,lcl))+min(lck,lcl)
                     f(kl) = t(k,l)
                  enddo
               enddo
c
            enddo
c
  500    continue
  520 continue
c
      return
      end
      subroutine ropr(t,mink,maxk,lkt,minl,maxl,llt,ntr)
c
      implicit REAL  (a-h,o-z)
      parameter (ndim=35)
      parameter (maxgro=48)
      common/bufb/ptr(3,3,maxgro),  dtr(6,6,maxgro),
     +            ftr(10,10,maxgro),gtr(15,15,maxgro)
      dimension t(ndim,ndim),v(ndim)
c
c     ----- right multiply  t  by  r,
c           result back in  t
c
      data zero /0.0d0/
c
c     ----- right multiply  t  by  r , result back in  t -----
c
      go to (500,400,300,200,100),llt
c
c     ----- g shell -----
c
  100 ng=15*(ntr-1)
      do k=mink,maxk
         do l=21,35
            dum=zero
            do n=21,35
               dum=dum+t(k,n)*gtr(n-20,l-20,ntr)
            enddo
            v(l)=dum
         enddo
         do l=21,35
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- f shell -----
c
  200 nf=10*(ntr-1)
      do k=mink,maxk
         do l=11,20
            dum=zero
            do n=11,20
               dum=dum+t(k,n)*ftr(n-10,l-10,ntr)
            enddo
            v(l)=dum
         enddo
         do l=11,20
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- d shell -----
c
  300 nd= 6*(ntr-1)
      do k=mink,maxk
         do l= 5,10
            dum=zero
            do n= 5,10
               dum=dum+t(k,n)*dtr(n- 4,l- 4,ntr)
            enddo
            v(l)=dum
         enddo
         do l= 5,10
            t(k,l)=v(l)
         enddo
      enddo
      go to 500
c
c     ----- p shell -----
c
  400 np= 3*(ntr-1)
      do k=mink,maxk
         do l= 2, 4
            dum=zero
            do n= 2, 4
               dum=dum+t(k,n)*ptr(n- 1,l- 1,ntr)
            enddo
            v(l)=dum
         enddo
         do l= 2, 4
            t(k,l)=v(l)
         enddo
      enddo
c
  500 continue
c
c     ----- left multiply  t  by r , result back in  t -----
c
      go to (1000,900,800,700,600),lkt
c
c     ----- g shell -----
c
  600 ng=15*(ntr-1)
      do l=minl,maxl
         do k=21,35
            dum=zero
            do n=21,35
               dum=dum+gtr(n-20,k-20,ntr)*t(n,l)
            enddo
            v(k)=dum
          enddo
          do k=21,35
             t(k,l)=v(k)
          enddo
      enddo
      go to 1000
c
c     ----- f shell -----
c
  700 nf=10*(ntr-1)
      do l=minl,maxl
         do k=11,20
            dum=zero
            do n=11,20
               dum=dum+ftr(n-10,k-10,ntr)*t(n,l)
            enddo
            v(k)=dum
          enddo
          do k=11,20
             t(k,l)=v(k)
          enddo
      enddo
      go to 1000
c
c     ----- d shell -----
c
  800 nd= 6*(ntr-1)
      do l=minl,maxl
         do k= 5,10
            dum=zero
            do n= 5,10
               dum=dum+dtr(n-4,k-4,ntr)*t(n,l)
            enddo
            v(k)=dum
         enddo
         do k=5,10
            t(k,l)=v(k)
         enddo
      enddo
      go to 1000
c
c     ----- p shell -----
c
  900 np= 3*(ntr-1)
      do l=minl,maxl
         do k= 2,4
            dum=zero
            do n= 2,4
               dum=dum+ptr(n-1,k-1,ntr)*t(n,l)
            enddo
            v(k)=dum
         enddo
         do k= 2, 4
            t(k,l)=v(k)
         enddo
      enddo
c
 1000 continue
      return
      end
      subroutine dhstar(q,f,d,prefac,rdmat,irest)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*),f(*),q(*),prefac(*),rdmat(*)
INCLUDE(common/sizes)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/statis)
INCLUDE(common/timeperiods)
INCLUDE(common/fock_critic)
INCLUDE(common/iofile)
      common/blkin/gin(510),nint
      data pt5/0.5d0/
      inull = igmem_null()
_IF1(civ)      if(irest.ne.1)call  szero(f,num*(num+1)/2)
_IFN1(civu)      if(irest.ne.1)call  vclr(f,1,num*(num+1)/2)
      do   1 m=1,num
      nij=ikyp(m)
   1  d(nij)=d(nij)*pt5
c
c...   for debug purposes fock matrix builder may be skipped
c
      if (ofipc_nofock) go to 9999
c
      call start_time_period(TP_DHSTAR)
      call timana(5)
      call jandk(zscftp,q,f,q(inull),q(inull),d,q(inull),prefac,rdmat)
      call end_time_period(TP_DHSTAR)
      call cpuwal(begin,ebegin)
c
9999  if (ofipc_nofock) write(iwr,*) '****** Fockbuilder skipped ******'
_IF(parallel)
c
c***   ***node-MPP***
c...   gather fock-matrices and send back to  all
      call start_time_period(TP_DHSTAR_GOP)
      call pg_dgop(201,f,nx,'+')
      call end_time_period(TP_DHSTAR_GOP)
c***   ***node-MPP***
_ENDIF
      do   2 m = 1,num
      nij = ikyp(m)
    2 d(nij) = d(nij)+d(nij)
      call dscal(num*(num+1)/2,pt5,f(1),1)
      return
      end
      subroutine dscrf(q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character*127 text
c
c -----    version to perform scrf calcs  -----
c -----    common blocks for scrf         -----
c
INCLUDE(common/scrf)
INCLUDE(common/gvalue)
INCLUDE(common/cslosc)
INCLUDE(common/modj)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/atmol3)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/segm)
INCLUDE(common/timeperiods)
INCLUDE(common/runlab)
INCLUDE(common/harmon)
INCLUDE(common/zorac)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     +             iposit(20),nsti(2),ondiis,junkj
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      dimension q(*),tmol(3),tele(3),tnuc(3)
c
       integer domo, vmo
       logical logtest
       common/testerp/tester(200), idomo(200), ivmo(200),
     +               testerd(200),idomod(200), ivmod(200),
     +               testdum(2,200),idomodum(2,200),
     +               logtest(200)
c
_IF(ga)
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
_ENDIF
      character*20 label
      common/msglab/label
      character*10 charwall

c
      data zscf/'scf'/
      data dzero,two,pt5/0.0d0,2.0d0,0.5d0/
      data yblnk,yav/' ',' av'/
      data done,pt2,twopt2 /1.0d0,0.2d0,2.2d0/
      data dmptlc/1.0d-2/
      data igs/5/
      data bigd/1.0d4/
      data small/1.0d-16/
      data m3/3/
            call check_feature('dscrf')
c
c
c     variables for monitoring tester and convergence
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      iterstat = 0
c
      out = nprint .eq. 5
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)
c 
_IF(parallel)
c
c - this should be default anyway -
c
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
*     write(6,*)'master',ipg_nodeid(),omaster
_ENDIF
c
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(nprint.ne.-5)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      lentri = l2
      l3 = num*num
c     l4 = num * 9
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
c
c note memory allocation is designed to make irdmat and 
c idmat contiguous - reversed
c
c     determine memory required, ilen
      ilen = 1 + ikyp(nshell)*3 + 6 * l2 + 2 * l1
c     allocate ilen
c
      irdmat = igmem_alloc(ilen)
      iprefa=irdmat+ikyp(nshell)
      i10   = iprefa+ikyp(nshell)
      i20   = i10+l2
      i21   = i10+l3
      i30   = i20+l2
      i31   = i20+l3
      i40   = i30+l2
      i50   = i40+l2
c ps: add extra space for brdcst buffer 
      i60   = i50+l2+1+ikyp(nshell)
      i80   = i60+l2
      i90   = i80+l1
      last  = i90+l1
      if (nprint .eq. 5) write (iwr,9308) i10,i20,i21,i30,i31,i40,i50,
     +     i60,i80,i90,last
c ***
      dlnmxd=0.0d0
      dlntol=tolitr(3)
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i30),l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,q(i20),q(i30))
      endif
_ENDIF
      call start_time_period(TP_RDMAT)
      m171t=171
      len171=lensec(l2)
      nshtri=ikyp(nshell)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      do i=4,6
       iof171(i)=iof171(i-1)+len171
      enddo
      if(irest.ne.0)write(iwr,2970)irest
      if(irest.ne.1.and.irest.ne.2) then
        call rdmake(q(iprefa))
        ln171=nshblk*2
        if(odelta) then
         ln171=4*len171+ln171
        endif
        call secput(isect(471),m171t,ln171,ibl171)
        if(nprint.ne.-5) then
        write(iwr,2980) ibl171,ln171
        endif
        call wrt3(q(iprefa),nshtri,ibl171,idaf)
        if(odelta) then
         call zer171(q(i10),l2,4,ibl171+iof171(3),len171,idaf)
        endif
        call clredx
      endif
      call end_time_period(TP_RDMAT)
c ***
      call start_time_period(TP_TEST1)
      if (maxcyc .le. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, maxcycp
       tester(loop) = 0.0d0
       testerd(loop) = 0.0d0
       logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1

      i=lensec(l2)
      iblkh0 =ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)
      if(numdis.eq.0)then
         numdis=num8
         ndafd=iblkh+i
      else
         ndafd=ibldis+2
      endif
      ndafdi=ndafd+3*i
      if(irest.eq.0.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         kcount = 0
         iterv = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      diffdp = bigd
      diffdpp = bigd
c
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (odamph .or. oshift) damp = done
      lockt=lock
      tim0=cpulft(1)
      ovir = .false.
      if(oaimpac())ovir = .true.
      en = enucf(nat,czan,c)
      if(nprint.ne.-5) then
        write (iwr,9008) en
        write (iwr,9028) maxcyc,mconv,nconv,mpunch,
     1  dexp(-dlntol)
      endif
      l0=newbas0
      call end_time_period(TP_TEST1)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
       call start_time_period(TP_ORFOG)
       call rdedx(q(i10),l2,ibl7st,num8)
       call qmat(q,q(i10),q(i20),q(i31),q(i40),iky,l0,l1,l3,l1,out)
c load guess vectors..
       call rdedx(q(i30),l3,ibl3qa,idaf)
c load overlap
       call rdedx(q(i50),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i50), num, l0)
       else
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     *         ilifq,l0,l1,1)
       endif
_ELSE
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     *         ilifq,l0,l1,1)
_ENDIF
       call end_time_period(TP_ORFOG)

       call start_time_period(TP_TEST2)
       call wrt3(q(i30),l3,ibl3qa,idaf)
       call tdown(q(i30),ilifq,q(i30),ilifq,l1)
       call rdedx(q(i90),l1,ibl3ea,idaf)
       call llvmo(q(i90),q(i80),na,nocmx,l1)
       call dscal(l1,two,q(i80),1)
c...
c...  construct initial density matrix
c...
       call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
       call wrt3(q(i50),l2,ibl3pa,idaf)
       call end_time_period(TP_TEST2)
_IF(parallel)
      endif
c...  other nodes receive density matrix in i50 if they were screened out
      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst idmat"
         call pg_brdcst(7129,q(i50),l2*8,0)
         label=" "
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      lprnt = l1
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      call timrem(tlefts)
c *** build rest of fock matrix from 2-e integrals
      ifock=i10
      idmat=i50
      if(irest.eq.1 .or. irest.eq.2) then
        call secget(isect(471),m171t,ibl171)
        write(iwr,2990)irest,ibl171
        call rdedx(q(iprefa),nshtri,ibl171,idaf)
        call reads(q(irdmat),nshtri,idaf)
        if(odelta) then
         call rdedx(q(i10),l2,ibl171+iof171(3),idaf)
         call reads(q(i50),l2,idaf)
        endif
        dlnmxd=-9999999.0d0
        do 1672 kkk=1,nshtri
            if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
            dlnmxd=q(irdmat+kkk-1)
1672    continue
_IFN(parallel)
        if(irest.eq.1) goto 6351
        if(irest.eq.2) goto 6352
_ENDIF
      endif
c ***
 140  continue
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4)=.false.
         oswed3(8)=.false.
      endif
_ENDIF
c
c  intermediate times
c      call list_time_periods(.false.,.false.)

c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** delta f, delta d, last f, current d, reduced dmat, prefac mat.
c *** note that delta f and delta d are NOT accessed if
c *** delta-SCF has not been requested : odelta = .false.)
      call start_time_period(TP_RDMAT)
_IF(parallel)
      if(omaster) then
      if(irest.eq.1) goto 6351
*     if(irest.eq.2) goto 6352
_ENDIF
      if (odelta) then
         call vclr(q(i60),1,l2)
         call wrt3(q(i60),l2,ibl171+iof171(3),idaf)
         call rdedx(q(i60),l2,ibl171+iof171(6),idaf)
         call wrt3(q(i50),l2,ibl171+iof171(6),idaf)
         call vsub(q(i50),1,q(i60),1,q(i50),1,l2)
         call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
      else
c        call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
      endif
      if(nprint.ne.-5) write(iwr,3030)
      call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
c ***
c *** make sure that iter before delta calc is of full accuracy.
c ***
      if( (dlnmxd.lt.deltol) .and. (iter-1.lt.itrtol(2))) then
        dlnmxd=deltol+1.0d0
      endif
      if( (dlnmxd.gt.deltol) ) then
        if(nprint.ne.-5) write(iwr,3000)
        if (odelta) then
           call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
           call wrt3(q(i50),l2,ibl171+iof171(4),idaf)
           call vclr(q(i60),1,l2)
           call wrt3(q(i60),l2,ibl171+iof171(5),idaf)
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        else
           call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        endif
        call clredx
      endif
      call wrt3(q(irdmat),nshtri,ibl171+iof171(2),idaf)
6351  continue
c ***
_IF(parallel)
      endif
c
      if(ipiomode .eq. IO_NZ_S)then
c ps: combine to 1 brdcst
         call pg_brdcst(7122,q(irdmat),nshtri*8,0)
         call pg_brdcst(7124,dlnmxd,8,0)
         call dcopy(nshtri,q(irdmat),1,q(idmat+l2+1),1)
         q(idmat+l2)=dlnmxd
         n4=l2+nshtri+1
         label="brdcst idmat"
         call pg_brdcst(7123,q(idmat),(l2+nshtri+1)*8,0)
         label=" "
         call dcopy(nshtri,q(idmat+l2+1),1,q(irdmat),1)
         dlnmxd= q(idmat+l2)
c
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF(parallel)
      call end_time_period(TP_RDMAT)
      if(iter.lt.itrtol(1)) then
         dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
         dlntol=tolitr(2)
      else
         dlntol=tolitr(3)
      endif
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(nprint.ne.-5) then
         write(iwr,3020) dlntol
      endif
c
      call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),irest)
c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
       call start_time_period(TP_TEST3)
      if(irest.eq.0) go to 6352
        write(iwr,3010)irest
        if(numdis.eq.num8) go to 460
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
        go to 460
c *** form complete 2-e skelton fmat and save it
6352  if (odelta) then
        call rdedx(q(i60),l2,ibl171+iof171(5),idaf)
        call vadd(q(i10),1,q(i60),1,q(i10),1,l2)
        call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
        call rdedx(q(i50),l2,ibl171+iof171(6),idaf)
      else
c       call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
      endif
c ***
      if (nprint .eq. 5) then
         write (iwr,9068)
         call prtril(q(i10),l1)
      endif
c  symmetrise skeleton fock matrix 
      call symh(q(i10),q(i20),iky,0,0)
      if (nprint .eq. 5) then
         write (iwr,9088)
         call prtril(q(i10),l1)
      endif
c
c ----- start of scrf procedure dvsg 17/9/90
c
c ----- write 2e matrix to ed7
      call wrt3(q(i10),l2,iblkh,num8)
c
c ----- calculate dipole moment
c ----- taken from rfdpol.f
c
c ----- read in dipole integrals
      call rdedx(q(i10),l2,ibl7x,num8)
      call rdedx(q(i20),l2,ibl7y,num8)
      call rdedx(q(i30),l2,ibl7z,num8)
c
c ----- read in density matrix
c     call rdedx(q(i50),l2,ibl171+iof171(4),idaf)
c
c ----- electronic part
      tele(1) = -tracep(q(i50),q(i10),l1)
      tele(2) = -tracep(q(i50),q(i20),l1)
      tele(3) = -tracep(q(i50),q(i30),l1)
c
c ----- nuclear part
      tnuc(1) = ddot(nat,czan,1,c(1,1),3)
      tnuc(2) = ddot(nat,czan,1,c(2,1),3)
      tnuc(3) = ddot(nat,czan,1,c(3,1),3)
c
c ----- totals x,y,z
      do 776 iq=1,3
      tmol(iq)=tnuc(iq)+tele(iq)
776   continue
c
c ----- total dipole moment
      dtot =ddot(m3,tmol,1,tmol,1)
      if(dtot.gt.small)dtot=dsqrt(dtot)
c ----- scale dipole integrals by g's and dipole moments
      do 211 loop=1,l2
                   q(i10+loop-1)=q(i10+loop-1)*tmol(1)*gx
                   q(i20+loop-1)=q(i20+loop-1)*tmol(2)*gy
                   q(i30+loop-1)=q(i30+loop-1)*tmol(3)*gz
211   continue
c
c ----- now add x,y,z components into i30
c ----- leaves i10,i20 for fock matrices
      do 212 loop=1,l2
           q(i30+loop-1)=q(i30+loop-1)+q(i20+loop-1)+q(i10+loop-1)
212   continue
c
c ----- read back fock matrices
      call rdedx(q(i10),l2,iblkh,num8)
      call rdedx(q(i20),l2,ibl7f,num8)
c...  add zora corrections (only atomic/restore)
      if (ozora) call zora(q,q,q(i20),'read')
c
c ----- add dipole terms to 1e matrix
      do 213 loop=1,l2
                   q(i20+loop-1)=q(i20+loop-1)+q(i30+loop-1)
213   continue
c
c ----- do hf energy calculation
      call adonee(q(i10),q(i20),q(iprefa))
      ehf1 = tracep(q(i50),q(i20),l1)
      ehf2 = tracep(q(i50),q(i10),l1)
      ehf0 = ehf
      ehf = (ehf1+ehf2)*pt5
c ----- energy expression for scrf
      ehf = ehf + gx*(dtot**2)*pt5
c ----- now allow for charged systems
      ehf = ehf + ((1.0d0-dielec)/(dielec*aradix))*(ich**2)
c
c ----- resume original coding
      call wrt3(q(i10),l2,iblkh,num8)
      if (nprint .eq. 5) then
        write (iwr,9048)
        call prtril(q(i10),l1)
        call rdedx(q(i90),l1,ibl3ea,idaf)
        iorb1 = 1
        iorb2 = nocmx
        call makfv(q(i10),q(i30),q(i30),q(i21),iky,iorb1,iorb2,l1,l1)
        call rdedx(q(i10),l3,ibl3qa,idaf)
        call makeij(q(i10),q(i30),q(i30),q(i21),l0,l0,l1,l0,l1)
        write (iwr,9228)
        call eijout(q(i30),q(i90),l0,nocmx,l1,iwr)
        write (iwr,9328) ehf1,ehf2,ehf
        call rdedx(q(i10),l2,iblkh,num8)
      endif
      iter = iter+1
c
c compute kinetic energy and virial if required
c
      if(ovir)call rhfvir(q(i50),ehf,en,ek,vir,num8,l1,l2)
      etot0 = etot
      etot = ehf+en
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
      if (damp .lt. dmpcut) damp = dmpcut
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafd,
     +numdis,iterv,1,1)
      diffpp=diffp
      diffp=diff
c
      call end_time_period(TP_TEST3)
c
      if(odiis) then
      call start_time_period(TP_DIIS)
c
_IF(ga)
      if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
         call diis_ga(q(i10),q(i30),q(i50),ibl3qa,diffd,domo,vmo)
      else
         call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     *        ndafdi,numdis,diffd,lockt,domo,vmo)
      endif
_ELSE
      call diiscd(q(i10),q(i30),q(i10),q(i50),iblkh0,ibl3qa,
     *  ndafdi,numdis,diffd,lockt,domo,vmo)
_ENDIF
      if (iter.le.maxcycp) then
       idomod(iter) = domo
       ivmod(iter) =  vmo
       testerd(iter) = diffd
      endif
      call end_time_period(TP_DIIS)
      endif
c
      call start_time_period(TP_TEST4)
c
c
c ------ take precautions if tester is increasing again
c
      if(iter.eq.1 .or. .not.odiis) goto 223
      if(ondiis) then
        if (oprind) then
         write(iwr,9824)diffd,diffdp
         write(iwr,9823)diffd, acurcy*diisf
         write(iwr,*) 'de = ', de
        endif
        if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.diffd.le.acurcy*diisf
     1  .or.(otestd.and.de.lt.1.0d-10)) then
*        if (de.lt.1.0d-5) then
          if(iter.le.maxcycp) logtest(iter) = .true.
          diffdpp = diffdp
          diffdp = diffd
          go to 223
*        else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
*         odynamic = .true.
*         accdi1 = 0.0d0
*         itdiis = 0
*         if(oprind) write(iwr,9828) iter
*        endif
        else
          if(oprind) write(iwr,9829) iter
        endif
        nsti(1)=0
        ondiis=.false.
        call rdedx(q(i10),l2,iblkh,num8)
        diffdp = bigd
        diffdpp = bigd
      else
        if (odynamic) then
         if(diff.le.diffp.or.diffp.le.diffpp) then
          itdiis = itdiis + 1
          if (itdiis.eq.3) then
            accdi1 = diff
            if (oprind) write(iwr,9826) accdi1
          endif
         else
          if (oprind) write(iwr,9827)
          itdiis = 0
         endif
        endif
      endif
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -q- at q(i30) orthonormalizing transformation
c     -h- at q(i10)
c    -h'- at q(i50) transformed h matrix
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
 223  if(oshift)then
         call rdedx(q(i30),l3,ibl3qa,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diff=0.0d0
      ii=nocmx+1
      if(ii.le.l0)then
         do 250 i=ii,l0
            ij=iky(i)+i50
            loop=idamax(nocmx,q(ij),1)
            if(loop.gt.0)then
c            diff=dmax1(diff,dabs(q(ij+loop-1)))
             dum = dabs(q(ij+loop-1))
             if (dum.gt.diff) then
               diff = dum
               if(iter.le.maxcycp) then
                tester(iter) = diff
                idomo(iter) = loop
                ivmo(iter) = i
               endif
             endif
            endif
 250     continue
      endif
      if(ondiis)diff=diffd
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
c
_IF(parallel)
      if(ipiomode.eq.IO_NZ_S.and.odpdiag(l1))then
c  in screened case, broadcast convergence event 
         dskip = 0.0d0 
         if (ocvged) dskip = 1.0d0  
         label="brdcst dskip"
         call pg_brdcst(7128,  dskip,  8,  0)
         label=" "
      endif
_ENDIF

      call end_time_period(TP_TEST4)
      if (ocvged) goto 321
      rshift=0.0d0
      if(.not.ondiis)
     *  call shiftq(q(i50),nocmx,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)
c
c broadcast fock for parallel diag
      if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &    call pg_brdcst(7127,  q(i50),  l1*(l1+1)/2*8,  0)
C
      m2=2
      diaacc = diff*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
      call start_time_period(TP_DIAG)
      call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
      call end_time_period(TP_DIAG)
      call start_time_period(TP_TEST5)
      if(oshift) then
        call rdedx(q(i10),l3,ibl3qa,idaf)
      else
        call rdedx(q(i10),l3,ibl3qs,idaf)
      endif
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      ig = mod(iter,igs)
      call end_time_period(TP_TEST5)
      if (ig .eq. 0) then
         call start_time_period(TP_ORFOG)
         call rdedx(q(i10),l2,ibl7st,num8) 
_IF(ga)
         if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
            call porth(q(i30),q(i10), num, l0)
         else
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1) 
         endif
_ELSE
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1) 
            call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1) 
_ENDIF
         call end_time_period(TP_ORFOG)
      endif
c
      call start_time_period(TP_TEST6)
c
      call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      ig=nocmx+1
      if(ig.le.l0) then
_IFN1(civu)        call vsadd(q(nocmx+i90),1,-rshift,q(nocmx+i90),1,l0-nocmx)
_IF1(civu)        call vsav(l0-nocmx,-rshift,q(nocmx+i90),q(nocmx+i90))
      endif
      if (nprint .eq. 5) then
         write (iwr,9148)
         call prev(q(i30),q(i90),lprnt,l1,l1)
      endif
      call llvmo(q(i90),q(i80),na,nocmx,l1)
      call dscal(l1,two,q(i80),1)
      yavr=yblnk
      if(nocmx.gt.na)yavr=yav
      call dmtx(q(i50),q(i30),q(i80),iky,nocmx,l1,l1)
      if(.not.oprint(20)) lprnt=min(nocmx+5,l1)
      if (nprint .eq. 5) then
         write (iwr,9168)
         call prtril(q(i50),l1)
      endif
c
c     compute RMS convergence on density
c
      if (iter.ne.1) then
       call rdedx(q(i10),l2,ibl3pa,idaf)
       dsq = cvgdens(q(i50),q(i10),l2)
       dsq = dsq / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
c
c     ----- save mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i50)
c     -d- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i50),q(i90),q(i80),ndaf,l1,l2,
     *ibl3pa,ibl3ea)
c
       call end_time_period(TP_TEST6)
c
 321  if(numdis.ne.num8) then
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
      endif
      tim1=cpulft(1)
      delt = tim1-tim0
      tim0 = tim1
c
c     monitor static tester with diis on - at the moment remove diis
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if( oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif

      otrigp = maxcyc-iter.lt.10
      if((nprint.ne.-5).or.otrigp) then
         write(iwr,90281)
         write (iwr,9188) iter,kcount,etot,
     1   ehf,de,diff,rshift,damp,derror,delt,tim1,yavr
      endif
      write(text,9188) iter,kcount,etot,ehf,de,diff,rshift,damp,
     &                 derror,delt,tim1,yavr
      call sptchk(text)
_IF(parallel)
_IF(diag_parallel)
c
c broadcast of back-transformed vectors 
c probably this is a waste of time
c
c      if (ipiomode .eq. IO_NZ_S .and. .not.ocvged) then
c         label="brdcst vect"
c         call pg_brdcst(7126,q(i30),l1*l1*8,0) 
c         label=" "
c      endif
_ENDIF
      else
c
c    stuff for other screened-out nodes 
c    this is only necessary if these nodes are joining in the
c    parallel diag
c    in future extra sections may be needed to enable parallel diis with
c    ipiomode = IO_NZ_S
c
         iter = iter + 1
         if(irest.gt.0) then
            write(iwr,3010)irest
            go to 460
         endif

         if(odpdiag(l1)) then
c  
c  find out if energy has already converged 
          label="brdcst dskip"
          call pg_brdcst(7128,  dskip,  8,  0) 
          label=" "
          if (dskip .lt. 0.9d0) then 
c
c  receive diis fock from root node, and diagonalise
             label="brdcst fock"
             call pg_brdcst(7127,  q(i50),  l1*(l1+1)/2*8,  0)
             label=" "
c
             call start_time_period(TP_DIAG)
             m2=2
             diaacc = diff*5.0d-3
             if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
             if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c     
             call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),
     &            m2,lock,diaacc)
             call end_time_period(TP_DIAG)
          endif 
          endif 
c
c gdf:   gather back-transformed vectors 
c probably wasted, as slave nodes don't neet to know the vectors?
c       label="brdcst vect"
c       call pg_brdcst(7126,q(i30),l1*l1*8,0) 
c       label=" "
c 3219  continue 

      endif
c
      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst maxcyc"
         lscf = 20+12/nav
         call pg_brdcst(7125,maxcyc,lscf*8,0)
         label=" "
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif
c
      call pg_synch(5555)
c
_ENDIF
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy).and.(iter.gt.1)
      if (.not.ocvged) then
         call timit(0)
         call timrem(tlefti)
         if (tlefti .le. skale*(tlefts-tlefti)/iter) then
_IF(parallel)
         write(6,*)'** time-out on node ',ipg_nodeid()
_ENDIF
         irest=3
         nindmx = 0
         write(iwr,9367)
         call texit(0,irest)
      else
         if (iter .lt. maxcyc) go to 140
         write (iwr,9288)
c        write out error description to xml/punchfile
         call blkerror('excessive number of SCF iterations',0)
         etot = dzero
         irest=3
         nindmx = 0
         ehf = -en
      endif
      endif
c
  400 continue
      lock = lockt
      if (ocvged.and.nprint.ne.-5) write (iwr,9208)
      if (omodel) then
       write (iwr,9371) iter,tim1,charwall(),ehf,en,etot,diffdens
      else if (ovir)then
       write (iwr,9365) iter,tim1,charwall(),ek,ehf,en,etot,vir,diffdens
      else
       write (iwr,9368) iter,tim1,charwall(),ehf,en,etot,diffdens
      endif
      if(nprint.ne.-5.and.yavr.ne.yblnk)write(iwr,9369)
      if(numdis.eq.num8)numdis=0
c
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       idum = min(iter,200)
       write(iwr,9199)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298) loop, tester(loop), idomo(loop), ivmo(loop),
     +                     testerd(loop),idomod(loop),ivmod(loop)
       else
       write(iwr,9198) loop, tester(loop), idomo(loop), ivmo(loop),
     +                       testerd(loop),idomod(loop),ivmod(loop)
       endif
       enddo
       write(iwr,9197)
       endif
c
      if (nprint .eq. 5 .or. nprint .eq. -5) go to 440
_IF(parallel)
      if(omaster) then
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
_ENDIF
       call rdedx(q(i10),l3,ibl3qa,idaf)
       call rdedx(q(i30),l2,ibl3pa,idaf)
       call rdedx(q(i21),l1,ibl3ea,idaf)
       call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
       call tdown(q(i10),ilifq,q(i10),ilifq,l0)
       call prev(q(i10),q(i21),lprnt,l1,l1)
       if(nprint.eq.-5) then
         write (iwr,9168)
         call prtril(q(i30),l1)
       endif
_IF(parallel)
      else
c...  other roots
c...  allow for section created in symass
c...  save results of the symmetry assignment into ed3
c...  assuming this has been called
      endif
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
  440 continue
_IF(ga)
      call destroy_diis_storage(.false.)
_ENDIF
_ELSE
  440 continue
_ENDIF
      if (mpunch .ne. 0) then
         call rdedx(q(i10),l3,ibl3qa,idaf)
         call pusql(q(i10),na,l1,l1)
      endif
  460 continue
      if (ocvged) irest = 0
      if(nprint.ne.-5)then
         cpu=cpulft(1)
         write(iwr,7777)cpu,charwall()
      endif
      call gmem_free(irdmat)
c ***
c *** set tolitr so that if reenter in geom opt with old vecs
c *** start off with full accuracy
c ***
      tolitr(1)=tolitr(3)
      tolitr(2)=tolitr(3)
c
      accdi1 = accdin
      call copq_again(mouta)
      return
c
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9827 format(1x,'rising diff - de = ', f20.10)
 9825 format(/' *** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffd, diffdp = ',2f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9198 format(10x,i5,f11.7,' (',2i4,') (*)',6x,f11.7,' (',2i4,')' )
 9298 format(10x,i5,f11.7,' (',2i4,')',10x,f11.7,' (',2i4,') (*)' )
 9199 format(/10x,64('+')/
     + 10x,'CONVERGENCE / TESTER ANALYSIS',12x,'From DIIS'/
     + 10x,64('+')/
     + 10x, 'Iter.', '   Tester   (domo/vmo)',
     + 10x,          '   Tester   (domo/vmo)'/)
 9197 format(10x,64('+'))
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 2970 format(' restart parameter in direct-scf ',i2)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))
 3010 format(' restart parameter after dhstart ',i3)
 7777 format(//
     *' end of direct closed shell scrf scf at ',f12.2,
     *' seconds',2x,a10,' wall',//1x,104('-')//)
 9369 format(//20x,90('*')/
     *20x,'*',88x,'*'/20x,'*',88x,'*'/
     *20x,'*',31x,'warning  state is averaged',31x,'*'/
     *20x,'*',88x,'*'/20x,90('*')//)
 2990 format(/1x,'restarting direct-scf in integrals :',
     *  'restart parameter =',i3/
     *  1x,'section 171 at block ',i5)
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'scf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 9008 format(/' ----- nuclear energy ----- = ',f20.12)
 9028 format(//15x,'convergence data'/15x,16('=')//
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criterion        =1.0e-',i2/
     +     ' punch out option             = ',i6/
     +     ' integral prefactor tolerance = ',e10.2/)
90281 format(/1x,123('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',3x,'virtual',1x,'damping',11x,'diis',
     * 4x,'del(t)',6x,'time'/
     * 18x,'energy',10x,'energy',36x,'shift'/1x,123('='))
 9188 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,f10.3,f8.3,f15.9,
     1       2f10.3,a3)
 9048 format(/20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23(
     +     '-'))
 9148 format(//1x,104('-')//
     *  50x,12('-')/50x,'eigenvectors'/50x,12('-'))
 9168 format(/20x,14('-')/20x,'density matrix'/20x,14('-'))
 9208 format(/20x,16('-')/20x,'energy converged'/20x,16('-'))
 9228 format(/20x,27('-')/20x,'lagrange multipliers matrix'/20x,
     +     27('-'))
 9288 format(/20x,30('-')/20x,'excessive number of iterations'/
     +     20x,30('-'))
 9308 format(1x,'core assignment'/' i10, i20, i21,',
     +     ' i30, i31, i40, i50, i60, i80, i90  = '/10i8/
     +     ' last = ', i8)
 9328 format(/' ehf1 = ',f20.12,' ehf2 = ',f20.12,' ehf = ',f20.12)
 9348 format(/36x,36('*')/
     +        36x,'direct closed-shell scrf calculation', ' - v1.1'/
     +        36x,36('*'))
 9365 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     +,f12.2,' seconds',2x,a10,' wall',/
     +10x,14('-')/
     +10x,'kinetic energy         ',f18.10/
     +10x,'electronic energy      ',f18.10/
     +10x,'nuclear energy         ',f18.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy           ',f18.10,8x,f10.7/
     +10x,'convergence on density ',f18.10)
 9368 format(/20x,14('-')/
     +20x,'final energies   after',i4,' cycles at ',f12.2,' seconds ',
     +a10,' wall',/,20x,14('-')//
     +20x,'electronic energy      ',f18.10/
     +20x,'nuclear energy         ',f18.10/
     +20x,'total energy           ',f18.10/
     +20x,'convergence on density ',f18.10)
 9371 format(/20x,14('-')/
     +20x,'final energies   after',i4,' cycles at ',f12.2,' seconds',
     +2x,a10,' wall',/20x,14('-')//
     +20x,'electronic energy      ',f18.10/
     +20x,'nuclear energy         ',f18.10/
     +20x,'total qm energy        ',f18.10/
     +20x,'convergence on density ',f18.10)
      end
      subroutine rdmake(prefac)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension prefac(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/prints)
INCLUDE(common/iofile)
c ***
c *** form reduced integral prefactor matrix
c *** over shells
c *** form is max over contractions in a pair of shells
c *** ( 0.5*pi**2.5 )**2 * ca * cb * exp(-rab**2 * ab/(a+b))
c ***    /   (a+b) * 2**0.5 * (a+b)**0.25
c ***
      pi125=3.141592654d0 ** 1.25d0
      nnshl=iky(nshell)+nshell
_IFN1(cfuivx)      call dfill(nnshl,-1.d20,prefac,1)
_IF1(fx)      call vfill(-1.d20,prefac,1,nnshl)
_IF1(cu)      call setsto(nnshl,-1.d20,prefac)
_IF1(iv)      call setstr(nnshl,-1.d20,prefac)
      ij=1
      do 10 i=1,nshell
          ipt=kstart(i)
          iat=katom(i)
          do 20 j=1,i
              jpt=kstart(j)
              jat=katom(j)
              rijsq=(c(1,iat)-c(1,jat))**2+(c(2,iat)-c(2,jat))**2+
     1              (c(3,iat)-c(3,jat))**2
              iipt=ipt
              do 30 ii=1,kng(i)
                 atmp = dmax1(dabs(cs(iipt)),dabs(cp(iipt)))
                 atmp = dmax1(atmp,dabs(cd(iipt)),dabs(cf(iipt)))
                 ci = dmax1(atmp,dabs(cg(iipt)))
                 ei=ex(iipt)
                 jjpt=jpt
                 do 40 jj=1,kng(j)
                    atmp = dmax1(dabs(cs(jjpt)),dabs(cp(jjpt)))
                    atmp = dmax1(atmp,dabs(cd(jjpt)),dabs(cf(jjpt)))
                    cj= dmax1(atmp,dabs(cg(jjpt)))
                    ej=ex(jjpt)
                    atmp = dlog(ci*cj*pi125) - 1.25d0*dlog(ei+ej) -
     2                   rijsq*ei*ej/(ei+ej)
                    prefac(ij)=dmax1(prefac(ij),atmp)
                    jjpt=jjpt+1
 40              continue
                 iipt=iipt+1
30            continue
              ij=ij+1
20        continue
10    continue
      if(oprint(57)) then
        write(iwr,3000)
 3000   format(/30x,'integral prefactor matrix over shells'/
     *  30x,37('='))
         call writel(prefac,nshell)
      endif
      return
      end
      subroutine mkrdmt(zscftp,rdmat,dmat,l2,nprint)
c ***
c *** make reduced density matrix over shells
c ***
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      character *(*) zscftp
c
INCLUDE(common/sizes)
INCLUDE(common/prints)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/nshel)
INCLUDE(common/cslosc)
c 
      dimension dmat(l2,*),rdmat(*)
c
      ij=1
      ijkl=ikyp(nshell)
_IFN1(cuiv)      call vfill(1.0d-20,rdmat,1,ijkl)
_IF1(cu)      call setsto(ijkl,1.0d-20,rdmat)
_IF1(iv)      call setstr(ijkl,1.0d-20,rdmat)
      do 10 i=1,nshell
      mini=kmin(i)
      maxi=kmax(i)
          iloc=kloc(i)-mini
          do 20 j=1,i
          minj=kmin(j)
          maxj=kmax(j)
              jloc=kloc(j)-minj
              dmax = rdmat(ij)
              do 30 ii=mini,maxi
                  ibfn=iloc+ii
                  do 40 jj=minj,maxj
                    jbfn=jloc+jj
                    ijbfn=iky(max(ibfn,jbfn))+min(ibfn,jbfn)
                    if (zscftp.eq.'uhf') then
                        tmpa = dmat(ijbfn,1)
                        tmpb = dmat(ijbfn,2)
                        if (dabs(tmpa+tmpb).gt.dmax)
     +                      dmax = dabs(tmpa+tmpb)
                        if (dabs(tmpa-tmpb).gt.dmax)
     +                      dmax = dabs(tmpa-tmpb)
                    else if (zscftp.eq.'gvb'.or.zscftp.eq.'grhf') then
                      do 80 ishl = 1,nsheld
                      if (dabs(dmat(ijbfn,ishl)).gt.dmax)
     +                    dmax = dabs(dmat(ijbfn,ishl))
80                    continue
                    else
_IF1()*                      rdmat(ij)=dmax1(rdmat(ij), dabs(dmat(ijbfn,1)))
                      if (dabs(dmat(ijbfn,1)).gt.dmax) 
     +                    dmax = dabs(dmat(ijbfn,1))
                    endif
40                continue
30            continue
              rdmat(ij) = dmax
20            ij=ij+1
10    continue
      do 60 i=1,ijkl
   60 rdmat(i)=dlog(rdmat(i))
_IF1(v)      dlnmxd=-20.0d0
_IF1(v)      do 1 kk=1,ijkl
_IF1(v)   1      dlnmxd=dmax1(rdmat(kk),dlnmxd)
_IFN1(v)      loop=idmax(ijkl,rdmat,1)
_IFN1(v)      dlnmxd=rdmat(loop)
      if(oprint(57)) then
        write(iwr,1000)
1000    format(/30x,'reduced density matrix over shells'/
     *  30x,34('='))
        call writel(rdmat,nshell)
      endif
      if(nprint.ne.-5.and.zscftp.ne.'guess') write(iwr,70) dlnmxd
 70   format(5x,'dlnmxd =  ',f15.10)
      return
      end
      subroutine adonee(fock,tpv,prefac)
c
c *** add one electron integrals to the fock matrix.
c *** have to also select on these to conserve charge.
c ***
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension fock(*),tpv(*),prefac(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/cslosc)
INCLUDE(common/nshel)
c
      ij=1
      do 10 i=1,nshell
      mini=kmin(i)
      maxi=kmax(i)
          iloc=kloc(i)-mini
          do 20 j=1,i
          minj=kmin(j)
c         maxj=kmax(j)
              jloc=kloc(j)-minj
              if(prefac(ij)+dlntol.le.0.0d0) goto 20
              do 30 ii=mini,maxi
                  ibfn=iloc+ii
                  jjhi=kmax(j)
                  if(i.eq.j) jjhi=ii
                  do 40 jj=minj,jjhi
                      jbfn=jloc+jj
                      ijbfn=iky(max(ibfn,jbfn))+min(ibfn,jbfn)
                      fock(ijbfn)=fock(ijbfn)+tpv(ijbfn)
40                continue
30            continue
20            ij=ij+1
10    continue
      return
      end
      subroutine zer171(q,l2,ntri,ibl171,len171,idaf)
      implicit REAL  (a-h,o-z)
      dimension q(*)
      call vclr(q,1,l2)
      ib=ibl171
      do 667 kq1=1,ntri
          call wrt3(q,l2,ib,idaf)
667       ib=ib + len171
      return
      end
_IF(drf)
      function cvgden(d0,d1,l)
      implicit REAL  (a-h,o-z),integer  (i-n)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c
      dimension d0(*),d1(*)
      data dzero /0.0d+00/
c
      l2 = (num*(num+1))/2
      dmax=dzero
      do 10 i=1,l
c
c  -----  difference density in d1 (to be saved in drf extensions)
c
        d1(i) = d0(i) - d1(i)
        if (i .le. l2) then
          dum = abs(d1(i))
          if(dum .gt. dmax) dmax = dum
        endif
   10 continue
      cvgden = dmax
      return
      end
_ENDIF
      function cvgdens(d0,d1,l2)
      implicit REAL  (a-h,o-z),integer  (i-n)
c
      dimension d0(l2),d1(l2)
c
      dsq = 0.0d0
      do i = 1,l2
       dum = dabs(d0(i) - d1(i))
       dsq = dum * dum + dsq
      enddo
      cvgdens = dsq
      return
      end
      subroutine duhfop(sz,s2,q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     ----- unrestricted direct-scf calculation -----
c
INCLUDE(common/sizes)
INCLUDE(common/fsymas)
INCLUDE(common/modj)
INCLUDE(common/mapper)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/timez)
INCLUDE(common/scra7)
INCLUDE(common/infoa)
INCLUDE(common/scfopt)
INCLUDE(common/restar)
INCLUDE(common/atmol3)
INCLUDE(common/segm)
INCLUDE(common/timeperiods)
INCLUDE(common/runlab)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/drfopt)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/datgue)
INCLUDE(common/harmon)
INCLUDE(common/fermidirac)
_IF(drf)
cahv --- DRF extension ---
INCLUDE(../drf/comdrf/sizesrf)
INCLUDE(../drf/comdrf/drfdaf)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/drfbem)
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
      common/enrghf/enhh,etothh,ehfhh
      logical odrf
cahv
_ENDIF
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c  i/o vs message passing logic controlled by ipiomode
c
c  IO
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
c ... for dummy symass section
      parameter(isymtp=99)
      character*20 label
      common/msglab/label
_ENDIF
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      Logical use_symmetry, force_serial
      Integer n_irreps
c
       logical logtest
       integer domoa, vmoa, domob, vmob
       common/testerp/testera(200), idomoa(200), ivmoa(200),
     +                testerb(200), idomob(200), ivmob(200),
     +               testerda(200),idomoda(200), ivmoda(200),
     +               testerdb(200),idomodb(200), ivmodb(200),
     +               logtest(200)
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/st(210),cdiis(20),rdiis(19),derr,sdiis(20),
     + ipos(20),nsti(2),ondiis,nspce
c
      dimension q(*)
c     dimension osign(2)
INCLUDE(common/zorac)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/drive_dft)
INCLUDE(common/statis)
_ENDIF
      character*10 charwall
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
      data zscf/'scf'/
      data dzero,done,two     /0.0d0,1.0d0,2.0d0/
      data pt2,pt5,twopt2 /0.2d0,0.5d0,2.2d0/
      data dmptlc/1.0d-02/
      data igs /5/
      data bigd/1.0d4/
            call check_feature('duhfop')
_IF(drf)
cahv
      odrf = .false.
cahv
_ENDIF
      odft = .false.
_IF(ccpdft)
      idum = CD_update_geom(c)
      odft = CD_active()
_ENDIF
      dft_accu = 1.0d0
      out = nprint .eq. 5
      outon = nprint .ne. -5
c
c     variables for monitoring tester +  convergence
c     and for tracking inadequate level shifting
c
      ofalsec = .false.
      oprind = optester
      etot = 0.0d0
      etot0 = 0.0d0
      emrmna = 0.0d0
      emrmnb = 0.0d0
      emrmn = 0.0d0
      esmeara = esmear_start
      esmearb = esmear_start
      nocca  = na
      noccb  = nb
      nocmxa = na
      nocmxb = nb
      diffdens = 0.0d0
      iterstat = 0
c
c     iwild = 0
c     oreset = .false.
c
      nav = lenwrd()
      diisf = 10.0d0**(idiisf)
      call secget(isect(490),51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
c
c     is delta density ever to be invoked? if not, we can
c     reduce i/o activity around section 471 ..
c
      odelta = deltol .gt. dlog(acurcy)
c
_IF(parallel)
c
c - this should be default anyway -
c
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
*     write(6,*)'master',ipg_nodeid(),omaster

_ENDIF
      skale=2.0d0
      if(zruntp.eq.zscf)skale=1.1d0
      if(nprint.ne.-5)write (iwr,9348)
      l1 = num
      l2 = (num*(num+1))/2
      l22= l2+l2
      l3 = num*num
c
      nshtri=ikyp(nshell)
c
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
c
c   ----- set pointers for partitioning of core -----
c
c    i10    i20    i30    i40    i50    i60    i70    i80  i81  i90
c    -----  -----  -----  -----  -----  -----  -----  --   --   --
c    l2     l2     l2     l2     l2     l2     l2     l1   l1   l1
c
c
c     first evaluate total core requirements
      irdmat = 0
      iprefa = irdmat+nshtri
      i10 = iprefa+nshtri
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      i41 = i30+l3
      i50 = i40+l2
      i60 = i50+l2
      i70 = i60+l2
      i80 = i70+l2
      i81 = i80+l1
      i90 = i81+l1
      last = i90+l1
c
c     ----- get core memory and determine actual pointers -----
c
      irdmat = igmem_alloc(last)
c
      iprefa = irdmat+nshtri
      i10 = iprefa+nshtri
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      i31 = i20+l3
      i40 = i30+l2
      i41 = i30+l3
      i50 = i40+l2
      i60 = i50+l2
      i70 = i60+l2
      i80 = i70+l2
      i81 = i80+l1
      i90 = i81+l1
      last = i90+l1
      if (out) write (iwr,9188) i10,i20,i21,i30,i31,i40,i41,
     +i50, i60,i70,i80,i81,i90,last
_IF(drf)
cahv --- DRF extension ---
      itdrf = 0
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
c
c-----  set acur for use in 2-electron rf routines drfhamx
c
        acur = acurcy
        odrf = .true.
c
c      ovlap at i110  from da10, 12, l2
c      dipx  at i120  from da10, 53, l2
c      dipy  at i130  from da10, 54, l2
c      dipz  at i140  from da10, 55, l2
c      omega(s) at i150    da31, 50, nomega
c      omega(op) at i160   da31, 51, nomega, if itwoeps=1
c      iexp  at i170       da31, 2,  l2
c      ijbit(s) at i180    da31, 56, l2
c      ijbit(op) at i190   da31, 57, l2, if itwoeps=1 else clearit
c      fadrf     at i200   calc in drfhamu, l2
c      fbdrf     at i210   calc in drfhamu, l2
c      dadrf     at i220   calc in drfhamu, l2
c      dbdrf     at i230   calc in drfhamu, l2
c
        i110 = igmem_alloc(l2)
        i120 = igmem_alloc(l2)
        i130 = igmem_alloc(l2)
        i140 = igmem_alloc(l2)
        i150 = igmem_alloc(nomga)
        i200 = igmem_alloc(l2)
        i210 = igmem_alloc(l2)
        i220 = igmem_alloc(l2)
        i230 = igmem_alloc(l2)
        if (itwoeps .eq. 1)  i160 = igmem_alloc(nomga)
        i170 = igmem_alloc(l2)
        i180 = igmem_alloc(l2)
        i190 = igmem_alloc(l2)
      endif
cahv
_ENDIF
c
c     ----- occupation numbers -----
c
c     -oa- at q(i80)
c     -ob- at q(i81)
c
      do 120 i = 1,l1
      popa=dzero
      popb=dzero
      if(i.le.na)popa=done
      q(i-1+i80)=popa
      if(i.le.nb)popb=done
  120 q(i-1+i81) = popb
c 
      dlnmxd=0.0d0
      dlntol=tolitr(3)
      call start_time_period(TP_RDMAT)
      m171t=171
      len171=lensec(l2)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      do i=4,6
       iof171(i)=iof171(i-1)+len171+len171
      enddo
      if(irest.ne.0)write(iwr,2970)irest
      if(irest.ne.1.and.irest.ne.2) then
        call rdmake(q(iprefa))
        ln171=nshblk*2
        if(odelta) then
         ln171=8*len171+ln171
        endif
        call secput(isect(471),m171t,ln171,ibl171)
        if(nprint.ne.-5) then
           write(iwr,2980) ibl171,ln171
        endif
        call wrt3(q(iprefa),nshtri,ibl171,idaf)
        if(odelta)then
         call zer171(q(i10),l2,8,ibl171+iof171(3),len171,idaf)
        endif
        call clredx
      endif
      call end_time_period(TP_RDMAT)
_IF(drf)
cahv --- DRF extension ---
      if (odrf) then
c  -----  read necessary data from da10
        call daread(idafh,ioda,q(i110),l2,12)
        call daread(idafh,ioda,q(i120),l2,53)
        call daread(idafh,ioda,q(i130),l2,54)
        call daread(idafh,ioda,q(i140),l2,55)
c  -----  read necessary data from da31
        call daread(idafdrf,iodadrf,q(i170),l2,2)
        call daread(idafdrf,iodadrf,q(i150),nomga,50)
        call daread(idafdrf,iodadrf,q(i180),l2,56)
        if (itwoeps .eq. 1) then
          call daread(idafdrf,iodadrf,q(i160),nomga,51)
          call daread(idafdrf,iodadrf,q(i190),l2,57)
        else
          call clear(q(i190),l2)
        endif
        call clear(q(i200),l2)
        call clear(q(i210),l2)
      endif
cahv
_ENDIF
c
c     ----- initialize variables -----
c
      call start_time_period(TP_TEST1)
      if (maxcyc .lt. 0) maxcyc = 30
      maxcycp = min(200,maxcyc)
      do loop = 1, 200
      testera(loop) = 0.0d0
      testerb(loop) = 0.0d0
      testerda(loop) = 0.0d0
      testerdb(loop) = 0.0d0
      logtest(loop) = .false.
      enddo
      mpunch = 1
      if (npunch .ne. 1) mpunch = 0
      if (nprint .eq. 7) mpunch = 1
      if (numdis. eq. 0) numdis = num8
      if(irest.ne.3.or.numdis.eq.num8) then
         ehf = dzero
         ehf0 = dzero
         ek = dzero
         vir = dzero
         iter = 0
         iterv = 0
         kcount = 0
         damp = dzero
         damp0 = dzero
         rshift = dzero
         rshfta = dzero
         rshftb = dzero
         diff = dzero
         diffp = dzero
         diffpp = dzero
         de = dzero
         dep = dzero
         lockt = lock
         deavg = dzero
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      diffa = dzero
      diffb = dzero
      diffpa = dzero
      diffpb = dzero
      diffdpa = bigd
      diffdpb = bigd
      diffdp = bigd
      diffdpp = bigd
c
      accdin = accdi1
c
      if (odynamic) then
c
c     now make accdi1 (diis onset) dynamic
c
      accdi1 = 0.0d0
c
      endif
      itdiis = 0
c
      if (dmpcut .le. dzero) dmpcut = dzero
c
      len=lensec(l3)
      len2=lensec(l2)
      iblkh0a=ibl7la
      iblkh0b=iblkh0a+len2
      iblkqq=iblkh0b+len2
      iblkha=iblkqq+len
      iblkhb=iblkha+len2
      if(numdis.eq.num8)then
         ndafa=iblkhb+len2
      else
         ndafa=ibldis+2
      endif
      ndafb=ndafa+3*len2
      ndafd=ndafb+3*len2
      ocvged = .false.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (odamph .or. oshift) damp = done
      ovir = .false.
      if(oaimpac())ovir = .true.
c
      tim0=cpulft(1)
c
c     ----- nuclear energy
c
      en = enucf(nat,czan,c)
c
c     ----- compute canonical orthonormal vectors -----
c
c     -s - at q(i10)
c     -sv- at q(i20)
c     -se- at q(i31)
c
c     l0 = number of canonical orthonormal vectors kept.
c
      l0=newbas0
      call end_time_period(TP_TEST1)
c
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
       call start_time_period(TP_ORFOG)
c
       n_irreps = 0
       Do i = 1, l1
          isymmo( i ) = isymao( i )
          n_irreps = Max( n_irreps, isymao( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
       if(nprint.ne.-5) then
        write (iwr,9008) en
        if(use_symmetry) write(iwr,9009)
        write (iwr,9030) maxcyc,mconv,nconv,npunch
       endif
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(q(i20),l2,ibl7s,num8)
         call rdedx(q(i10),l2,ibl7st,num8)
         call declare_diis_storage(num,.true.)
         call init_diis_storage(num,q(i20),q(i10))
      endif
_ENDIF
      call rdedx(q(i10),l2,ibl7st,num8)
      if (omaster) then
      write(21)l2,(q(i10+i-1),i=1,l2)
      endif
      call qmat_symm(q,q(i10),q(i20),q(i31),q(i40),iky,l0,l1,l3,l1,out,
     +              isymmo, use_symmetry )
c
c     ----- compute initial guess density matrices -----
c
c     -da- in q(i40)
c     -db- in q(i50)
c
c ----- ensure vector orthogonality
c
       ipass=0
       nspin=na
       i=ibl3qa
       iocc=i80
       iblkp=ibl3pa
 801   call rdedx(q(i30),l3,i,idaf)
       call rdedx(q(i50),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i50), num, l0)
       else
          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
          if (nspin.ne.0) then
             call orfog(q(i30),q(i30),q(i10),q(i20),iky,
     +       ilifq,l0,l1,1)
          endif
       endif
_ELSE
       call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
       if (nspin.ne.0) then
          call orfog(q(i30),q(i30),q(i10),q(i20),iky,ilifq,l0,l1,1)
       endif
_ENDIF
      call wrt3(q(i30),l3,i,idaf)
      call tdown(q(i30),ilifq,q(i30),ilifq,l1)
      if (l1.ge.maxorb) call caserr('Out of range iky in dmtx(p)')
      if (ipass.eq.0) then
         nocca = na
         if (osmear) then
            call rdedx(q(i90),l1,ibl3ea,idaf)
            call fermi_smear(q(i90),q(iocc),emrmna,nocca,nocmxa,esmeara,
     &                       na,l0,iwr,oprint(60))
            nspin = nocmxa
         endif
      else
         noccb = nb
         if (osmear) then
            call rdedx(q(i90),l1,ibl3eb,idaf)
            call fermi_smear(q(i90),q(iocc),emrmnb,noccb,nocmxb,esmearb,
     &                       nb,l0,iwr,oprint(60))
            nspin = nocmxb
         endif
      endif
_IF(ga)
       if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
          call dmtxp(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
       else
          call dmtx(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
       endif
_ELSE
       call dmtx(q(i10),q(i30),q(iocc),iky,nspin,l1,l1)
_ENDIF
       if (ipass.eq.0) then
          i=ibl3qb
          iocc=i81
          ipass=ipass+1
          if (.not.uhfatom) call wrt3(q(i10),l2,iblkp,idaf)
          iblkp=ibl3pb
          nspin=nb
          go to 801
       endif
c
       if (uhfatom) then
c....     get density bmateices straight from atdens
          call rdedx(q(i40),l2,ibl3pa,idaf)
          call rdedx(q(i50),l2,ibl3pb,idaf)
          iter = iter - 1
       else
          call dcopy(l2,q(i10),1,q(i50),1)
          call wrt3(q(i10),l2,iblkp,idaf)
          call rdedx(q(i40),l2,ibl3pa,idaf)
       end if
_IF(drf)
cahv --- DRF extension ---
      if (odrf) then
        call dcopy(l2,q(i40),1,q(i220),1)
        call dcopy(l2,q(i50),1,q(i230),1)
      endif
cahv
_ENDIF
       call end_time_period(TP_ORFOG)
_IF(parallel)
      endif
c...  other nodes receive density matrices in i40 and 
c...  i50 if they were screened out
      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst idmat"
         call pg_brdcst(7129,q(i40),l22*8,0)
         label=" "
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
c
c     ----- start scf procedure -----
c
c     ----- construct a skeleton fock matrix -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c -scratch- at q(i30)
c     -da- at q(i40)
c     -db- at q(i50)
c
c     ----- find the time remaining at the start -----
c
      call timrem(tlefts)
      lprnta= l1
      if(.not.oprint(20)) lprnta= min(na+5,l1)
      lprntb= l1
      if(.not.oprint(20)) lprntb= min(nb+5,l1)
c
c *** build rest of fock matrix from 2-e integrals
      ifock=i10
      idmat=i40
      if(irest.eq.1 .or. irest.eq.2) then
        call secget(isect(471),m171t,ibl171)
        write(iwr,2990)irest,ibl171
        call rdedx(q(iprefa),nshtri,ibl171,idaf)
        call reads(q(irdmat),nshtri,idaf)
        if (odelta) then
          call rdedx(q(i10),l2,ibl171+iof171(3),idaf)
          call reads(q(i20),l2,idaf)
          call reads(q(i40),l2,idaf)
          call reads(q(i50),l2,idaf)
        endif
        dlnmxd=-9999999.0d0
        do 1672 kkk=1,nshtri
            if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
            dlnmxd=q(irdmat+kkk-1)
1672    continue
_IFN(parallel)
        if(irest.eq.1) goto 6351
        if(irest.eq.2) goto 6352
_ENDIF
      endif
c ***
_IF(ccpdft)
c
c     All memory allocations have been done apart from the ones
c     in SUBROUTINE DHSTARU. However, we have to reserve the memory
c     for the coulomb fit stuff here!
c
      if (CD_active()) then
         call retrieve_spare(imemspare)
         idum = CD_set_2e()
         if (CD_2e() .and. (.not. CD_HF_exchange())
     &      .and. (.not. CD_HF_coulomb()) )then
            imemdhstaru = 0
         else 
            call mem_dhstaru(imemdhstaru)
         endif
         idum = CD_reset_2e()
         idum = CD_uks()
         imemfree = igmem_max_memory() - imemspare - imemdhstaru
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q,q,iwr)
         else if (ks_bas.eq.KS_MO) then
            imemreq = CD_memreq_energy_mo(l1,na,nb,q,q,iwr)
         else if (ks_bas.eq.KS_AOMO) then
            imemreq = CD_memreq_energy(l1,na,nb,q,q,iwr)
         endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            write(iwr,600)ierror
            call caserr('Not enough memory for incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +             '3-center integrals in core')
c        if (CD_jfit_incore().and.CD_HF_exchange()) then
c           write(iwr,*)'*** WARNING: direct integral memory usage not',
c    +                  ' accounted for!!!'
c           write(iwr,*)'*** WARNING: Calculation may run out of ',
c    +                  'memory!!!'
c        endif
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
180   continue
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S) then
         oswed3(4)=.false.
         oswed3(8)=.false.
      endif
      call pg_synch(4444)
_ENDIF
c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** prefac matrix, red. dmat mat.
c *** delta fa,fb , delta da,db, last fa,fb, current da,db, 
c *** note that delta f and delta d are NOT accessed if
c *** delta-SCF has not been requested : odelta = .false.)
c ***
      call start_time_period(TP_RDMAT)
_IF(parallel)
      if(omaster) then
      if(irest.eq.1) goto 6351
*     if(irest.eq.2) goto 6352
_ENDIF
      if (odelta) then
        call vclr(q(i60),1,l2)
c       delta fa,fb
        call wrt3(q(i60),l2,ibl171+iof171(3),idaf)
        call wrt3s(q(i60),l2,idaf)
c       delta da,db
        call rdedx(q(i60),l2,ibl171+iof171(6),idaf)
        call reads(q(i70),l2,idaf)
        call wrt3(q(i40),l2,ibl171+iof171(6),idaf)
        call wrt3s(q(i50),l2,idaf)
        call vsub(q(i40),1,q(i60),1,q(i40),1,l22)
        call wrt3(q(i40),l2,ibl171+iof171(4),idaf)
        call wrt3s(q(i50),l2,idaf)
_IF(drf)
cahv --- DRF extension ---
      if (odrf .and. (itdrf .ne. 0)) then
        call dcopy(l2,q(i40),1,q(i220),1)
        call dcopy(l2,q(i50),1,q(i230),1)
      endif
cahv
_ENDIF
      else
c       call wrt3(q(i40),l2,ibl171+iof171(4),idaf)
c       call wrt3s(q(i50),l2,idaf)
_IF(drf)
cahv --- DRF extension ---
      if (odrf  .and. (itdrf .ne. 0)) then
        diffaa = cvgden(q(i40),q(i220),l2)
        diffbb = cvgden(q(i50),q(i230),l2)
      endif
cahv
_ENDIF
      endif
c
      if(nprint.ne.-5) write(iwr,3030)
      if (odelta) then
        call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
      endif
c ***
c *** make sure that iter before delta calc is of full accuracy.
c ***
      if( (dlnmxd.lt.deltol) .and. (iter-1.lt.itrtol(2))) then
        dlnmxd=deltol+1.0d0
      endif
      if( (dlnmxd.gt.deltol) ) then
        if(nprint.ne.-5) write(iwr,3000)
        if (odelta) then
          call rdedx(q(i40),l2,ibl171+iof171(6),idaf)
          call reads(q(i50),l2,idaf)
          call wrt3(q(i40),l2,ibl171+iof171(4),idaf)
          call wrt3s(q(i50),l2,idaf)
          call vclr(q(i60),1,l2)
          call wrt3(q(i60),l2,ibl171+iof171(5),idaf)
          call wrt3s(q(i60),l2,idaf)
          call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        else
          call mkrdmt(zscftp,q(irdmat),q(idmat),l2,nprint)
        endif
        call clredx
      endif
      call wrt3(q(irdmat),nshtri,ibl171+iof171(2),idaf)
6351  continue
c ***
_IF(parallel)
      endif
c
      if(ipiomode .eq. IO_NZ_S)then
         call pg_brdcst(7122,q(irdmat),nshtri*8,0)
         call pg_brdcst(7123,q(idmat),l22*8,0)
         call pg_brdcst(7124,dlnmxd,8,0)
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF(parallel)
      call end_time_period(TP_RDMAT)
      if(iter.lt.itrtol(1)) then
        dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
        dlntol=tolitr(2)
      else
        dlntol=tolitr(3)
      endif
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(nprint.ne.-5) then
      write(iwr,3020) dlntol
      endif
c
_IF(ccpdft)
c
c  Switch to UKS mode
c
        idum = CD_uks()
c
c skip dhstaru when not required for choice of DFT functional
c
        idum = CD_set_2e()

        if(CD_2e() .and. (.not. CD_HF_exchange())
     &      .and. (.not. CD_HF_coulomb()) )then

          o2e = .false.
          if(irest.ne.1)then
             call vclr(q(ifock),1,l2)
             call vclr(q(i20),1,l2)
          endif
       else
          o2e = .true.
          call dhstaru(q,q(idmat),q(i50),q(ifock),q(i20),
     &                 q(iprefa),q(irdmat),irest)
      endif

      idum = CD_reset_2e()

_ELSE
      call dhstaru(q,q(idmat),q(i50),q(ifock),q(i20),
     &             q(iprefa),q(irdmat),irest)
_ENDIF

c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
       call start_time_period(TP_TEST3)
      if(irest.eq.0) go to 6352
        write(iwr,3010)irest
        if(numdis.eq.num8) go to 580
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
        go to 580
c *** form complete 2-e skelton fmat and save it
6352  if (odelta) then
         call rdedx(q(i60),l2,ibl171+iof171(5),idaf)
         call reads(q(i70),l2,idaf)
         call vadd(q(i10),1,q(i60),1,q(i10),1,l22)
         call wrt3(q(i10),l2,ibl171+iof171(5),idaf)
         call wrt3s(q(i20),l2,idaf)
         call rdedx(q(i40),l2,ibl171+iof171(6),idaf)
         call reads(q(i50),l2,idaf)
      else
c        call rdedx(q(i40),l2,ibl171+iof171(6),idaf)
c        call reads(q(i50),l2,idaf)
      endif
c ***
      if (out) then
         write (iwr,9068)
         write (iwr,9108)
         call prtril(q(i10),l1)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
_IF(drf)
cahv --- DRF extension ---
c
c  reaction field contribution to 2-el part of hamiltonian
c
      if (odrf) then
        itdrf = itdrf + 1
c
c  fdrf     at q(i200) and q(i210)
c  diffdens at q(i220) and q(i230)
        call drfhamu(q(i200),q(i220),q(i210),q(i230),q(i110),
     1               q(i120),q(i130),q(i140),q(i170),q(i150),q(i180),
     2               q(i160),q(i190))
c
c  add rfc to fock matrix
        do 6666 i = 1, l2
          q(i10+i-1) = q(i10+i-1) + scffact*q(i200+i-1)
          q(i20+i-1) = q(i20+i-1) + scffact*q(i210+i-1)
 6666   continue
c
c  copy densities for diffdens calculation otherwise lost
        call dcopy(l2,q(i40),1,q(i220),1)
        call dcopy(l2,q(i50),1,q(i230),1)
      endif
cahv
_ENDIF
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c     scratch area at q(i30)
c
_IF(ccpdft)
c  if only 2e integrals were performed
      if(o2e)then
_ENDIF
      call symh(q(i10),q(i30),iky,0,0)
      call symh(q(i20),q(i30),iky,0,0)
      if (out) then
         write (iwr,9088)
         write (iwr,9108)
         call prtril(q(i10),l1)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
_IF(ccpdft)
      endif
_ENDIF
c
c     ----- read in core hamiltonian matrix
c           and calculate hf energy -----
c
c     -h0- at q(i30)
c     -fa- at q(i10)
c     -fb- at q(i20)
c     -da- at q(i40)
c     -db- at q(i50)
c
      call rdedx(q(i30),l2,ibl7f,num8)
c...   add zora corrections
       if (ozora) call zora(q,q,q(i30),'read')
      call adonee(q(i10),q(i30),q(iprefa))
      call adonee(q(i20),q(i30),q(iprefa))
c
c  ehf1 = TrP(P*H_1)
c
      ehf1 = tracep(q(i40),q(i30),l1) + tracep(q(i50),q(i30),l1)
c
c ps moved this up so as not to overlap with dft timings
c disabled for integral restarts
c
      if(irest.ne.2)call end_time_period(TP_TEST3)
c
c  save previous total 
c
      if (.not.(iter.le.0.and.uhfatom)) ehf0 = ehf

_IF(ccpdft)
      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(o2e)then
c
c Coulomb operator is in q(i10) augmented by H_1, 
c compute energy using HF expression without exchange
c
            ehf2 = tracep(q(i40),q(i10),l1) + tracep(q(i50),q(i20),l1)
            etmp = (ehf1+ehf2)*pt5
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            ehf2 = 0.0d0
            etmp = ehf1

         endif
c     
c Update Kohn-Sham matrix and compute fitted/integrated 
c energy terms
c
         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))
         call timana(5)
         call cpuwal(begin,ebegin)
_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(q(isma),1,l2)
_ENDIF
         oignore = CD_ignore_accuracy()
         if (uhfatom.and.iter.lt.0) then
             ierror = CD_set_ignore_accuracy(.true.)
         endif
         if (ks_bas.eq.KS_AO) then
            idum = CD_energy_ao(c,q(i10),q(i20),q(i40),q(i50),
     +           edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +           ,q(isma)
_ENDIF
     +           )
         else
            call caserr("duhfop: illegal ks_bas")
         endif
         ierror = CD_set_ignore_accuracy(oignore)
_IF(debug_S)
         call compare_S(q(isma),q(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF
         call symm_op(q(i10))
         call symm_op(q(i20))

         ehf = etmp+edft

          if(out) then
            if(opg_root()) then
                 call CD_print_dftresults(.true.,.false.,iwr)
                 write(iwr,9407)ehf1, ehf2
            endif
          endif
         call timana(31)
         call cpuwal(begin,ebegin)

      else
c
c  Unrestricted Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
c         ehf2 = tracep(q(i50),q(i10),l1)
c         ehf = (ehf1+ehf2)*pt5

         ehfa = tracep(q(i40),q(i30),l1)+tracep(q(i40),q(i10),l1)
         ehfb = tracep(q(i50),q(i30),l1)+tracep(q(i50),q(i20),l1)
         ehf = (ehfa+ehfb)/two

      endif
_ELSE

c      ehf2 = tracep(q(i50),q(i10),l1)
c      ehf = (ehf1+ehf2)*pt5

      ehfa = tracep(q(i40),q(i30),l1)+tracep(q(i40),q(i10),l1)
      ehfb = tracep(q(i50),q(i30),l1)+tracep(q(i50),q(i20),l1)
      ehf = (ehfa+ehfb)/two

_ENDIF
      emrmn = emrmna + emrmnb
      ehf   = ehf + 0.50d0*emrmn
c
c compute kinetic energy and virial if required
c
      if(ovir)call uhfvir(q(i40),q(i50),ehf,en,ek,vir,num8,l1,l2)
c
c     ----- save -fa- and -fb- on daf -----
c
c     -fa- at q(i10)
c     -fb- at q(i20)
c
      call wrt3(q(i10),l2,iblkha,num8)
      call wrt3(q(i20),l2,iblkhb,num8)
      if (omaster) then
         rewind(20)
         write(20)l2,(q(i10+i-1),i=1,l2)
         write(20)l2,(q(i20+i-1),i=1,l2)
      endif
c
      iter=iter+1
c
      etot0 = etot
      etot = ehf+en
      dep = de
      de = ehf-ehf0
      if (iter .eq. 1) deavg = dzero
      if (iter .eq. 2) deavg =  dabs(de)
      if (iter .ge. 3) deavg = (  dabs(de)+  dabs(dep)+pt2*deavg)/twopt2
c
c     ----- damp and extrapolate the alpha fock matrix -----
c
c     -f - at q(i10)     fock matrix  (n th matrix)
c     -fo- at q(i20)     old fock matrix (n-1 th matrix)
c     -fa- at q(i30)     ancient fock matrix (n-2 th matrix)
c     -fp- at q(i40)     prehistoric fock matrix (n-3 th matrix)
c
      if (iter .gt. 2) call dampd(de,dep,deavg,damp,acurcy,diff,diffp,
     +     dmptlc)
      if (damp .lt. dmpcut) damp = dmpcut
c
c     extrapolate alpha hamiltonian matrix
c
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafa,
     +  numdis,iterv,1,1)
      call wrt3(q(i10),l2,iblkh0a,num8)
c
c     now extrapolate beta hamiltonian matrix and write to disk
c
      call rdedx(q(i10),l2,iblkhb,num8)
c
      call extrpd(de,damp,damp0,q(i10),q(i20),q(i30),q(i40),l1,l2,ndafb,
     +numdis,iterv,2,1)
      call wrt3(q(i10),l2,iblkh0b,num8)
c
      diffpp = diffp
      diffp = max(diffpa,diffpb)
c
c     ----- invoke diis procedure
c
      if (odiis) then
        call start_time_period(TP_DIIS)
_IF(ga)
        if((l1.ge.idpdiis).and.(ipiomode.ne.IO_NZ_S))then
          call rdedx(q(i10),l2,iblkh0a,num8)
          call rdedx(q(i20),l2,iblkh0b,num8)
          call diisu_ga(q(i10),q(i20),q(i40),q(i60),ibl3qa,ibl3qb,
     +                  diffda,diffdb,domoa,vmoa,domob,vmob)
        else
          call diisu(q(i10),q(i20),q(i40),q(i60),ndafd,diffda,diffdb,
     +               lockt,
     +               iblkh0a,iblkh0b,numdis,domoa,vmoa,domob,vmob)
        endif
_ELSE
        call diisu(q(i10),q(i20),q(i40),q(i60),ndafd,diffda,diffdb,
     +             lockt,
     +             iblkh0a,iblkh0b,numdis,domoa,vmoa,domob,vmob)
_ENDIF
        if (iter.le.maxcycp) then
          idomoda(iter) = domoa
          ivmoda(iter) =  vmoa
          idomodb(iter) = domob
          ivmodb(iter) =  vmob
          testerda(iter) = diffda
          testerdb(iter) = diffdb
        endif
        call end_time_period(TP_DIIS)
      endif
c
      call start_time_period(TP_TEST4)
c
c ----- take precautions if tester is increasing again
c
      if(iter.eq.0.and.uhfatom) then
         itdiis = 0
         nsti(1) = 0
         kcount = 0
         go to 223
      endif
      if(iter.eq.1) go to 223
      if (odiis) then
       if (ondiis) then
           diffd = max(diffda,diffdb)
           if (oprind) then
            write(iwr,9824)diffda,diffdpa,diffdb,diffdpb
            write(iwr,9823)diffd, acurcy*diisf
            write(iwr,*) 'de = ', de
           endif
           if(diffd.lt.diffdp.or.diffdp.lt.diffdpp.or.
     1        diffd.le.acurcy*diisf.or.(otestd.and.de.lt.1.0d-10)) then
             if (de.lt.1.0d-5) then
              if(iter.le.maxcycp) logtest(iter) = .true.
              call wrt3(q(i20),l2,iblkhb,num8)
              diffdpa = diffda
              diffdpb = diffdb
              diffdpp = diffdp
              diffdp = max(diffdpa,diffdpb)
              go to 223
             else
c trouble .. energy is increasing with diis on
c switch to level shifting with dynamic diis onset
              odynamic = .true.
              accdi1 = 0.0d0
              itdiis = 0
              if(oprind) write(iwr,9828) iter
             endif
           else
            if(oprind) write(iwr,9829) iter
           endif
          nsti(1)=0
          ondiis=.false.
          call rdedx(q(i10),l2,iblkha,num8)
          call rdedx(q(i20),l2,iblkhb,num8)
          diffdp = bigd
          diffdpp = bigd
          diffdpa = diffdp
          diffdpb = diffdp
       else
c
c     first look for wild values of energy and tester
c     and re-level shift (only once ..). Problem here is
c     that this has to be recognised at the outset, or it
c     is not rescuable. Increase maxcyc to compensate.
c     After initial testing, I've deactivated this (for the
c     time being) - it triggers on too many occasions ...
c
c         if (dabs(de).gt.done.and.diff.gt.pt5.and.iter.gt.4) then
c          iwild = iwild + 1
c          osign(iwild) = de.lt.0.0d0
c          if (iwild.eq.2.and..not.oreset.and.
c    +         .not.(osign(1).and.osign(2)) ) then
c           gapa1 = 5.0d0
c           gapb1 = 5.0d0
c           gapa2 = gapa1
c           gapb2 = gapb1
c           oreset = .true.
c           iwild = 0
c           write(iwr,9831) gapa1
c           maxcyc = maxcyc + maxcyc
c           maxcycp = min(200,maxcyc)
c          endif
c         endif
          if (odynamic) then
           if(diff.le.diffp.or.diffp.le.diffpp) then
            itdiis = itdiis + 1
            if (itdiis.eq.3) then
              accdi1 = dmax1(diffa,diffb)
              if (oprind) write(iwr,9826) accdi1
            endif
           else
            if(oprind) write(iwr,9827) de
            itdiis = 0
           endif
           endif
          call wrt3(q(i20),l2,iblkhb,num8)
       endif
      endif
c
223   if (.not.odiis) call rdedx(q(i10),l2,iblkh0a, num8)
c
c     first, work on alpha orbitals
c
c     ----- read in fock tranformation matrix -----
c           transform hamiltonian matrix
c
c     -f- at q(i10)
c     -f'-at q(i20)  transformed f matrix
c     -q- at q(i30)  orthonormalizing tranformation
c     scratch area at q(i41)
c
c     ----- if vshift is true, use the previous set of
c           molecular orbitals to transform the hamiltonian matrix
c
      if(oshift)then
         call rdedx(q(i30),l3,ibl3qa,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diffpa = diffa
      diffa=0.0d0
      ii=nocca+1
      if(ii.le.l0) then
        do 250 i=ii,l0
          ij=iky(i)+i50
          loop=idamax(min(i-1,nocmxa),q(ij),1)
          if(loop.gt.0) then
            dum = dabs(q(ij+loop-1))
            if (dum.gt.diffa) then
              diffa = dum
              if(iter.le.maxcycp) then
                testera(iter) = diffa
                idomoa(iter)  = loop
                ivmoa(iter)   = i
              endif
            endif
          endif
 250    continue
      endif
c
      if(ondiis.and..not.osmear)diffa=diffda
      dmplim=dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diffa.lt.acurcy).and.(iter.gt.1)
      call end_time_period(TP_TEST4)
      if (maxcyc .eq. 0) go to 321
      if(ocvged)go to 321
c
c       shift the diagonal elements of the
c       transformed f matrix by rshift for the virtual part.
c
      rshift=0.0d0
      if((.not.ondiis.or.olevd).and..not.(uhfatom.and.iter.eq.0))
     *call shiftq(q(i50),nocmxa,0,l0,de,dep,iterv,1,gapa1,ibrk,gapa2)
      rshfta=rshift
c
      call start_time_period(TP_DIAG)
      m2=2
      diaacc = diffa*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c
c     ----- diagonalize new hamiltonian matrix -----
c
c     -h- at q(i50)
c     -v- at q(i30)
c     -d- at q(i41)
c
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
      call end_time_period(TP_DIAG)
c
      call start_time_period(TP_TEST5)
c
c     ----- back-transform the eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i90)
c     -q- at q(i10)
c     scratch area at q(i21)
c
      if (oshift) then
         call rdedx(q(i10),l3,ibl3qa,idaf)
      else
         call rdedx(q(i10),l3,ibl3qs,idaf)
      endif
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
      call end_time_period(TP_TEST5)
      ig = mod(iter,igs)
      if ( oshift .and. ig .eq. 0) then
        call start_time_period(TP_ORFOG)
c
c     ----- if vshift is true, reorthogonalize the vectors and
c           thereby the transformation matrix every igs th iteration.
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
        call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
        if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
           call porth(q(i30),q(i10), num, l0)
        else
           call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
           call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
        endif
_ELSE
        call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
        call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
c
        call end_time_period(TP_ORFOG)
      endif
      call start_time_period(TP_TEST6)
      call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      nocca = na
      ig = nocmxa
      if (osmear) then 
        esmeara= max(esmear_final,
     &               min(esmeara,esmear_scale*(diffa-acurcy)))
        call fermi_smear(q(i90),q(i80),emrmna,nocca,nocmxa,esmeara,na,
     &                   l0,iwr,oprint(60))
      endif
      if(ig.lt.l0) then
_IFN1(civ)      call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civ)      call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      endif
      if (out) then
         write (iwr,9208)
         write (iwr,9108)
         call prev(q(i30),q(i90),lprnta,l1,l1)
      endif
c
c     ----- form alpha density matrix -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -o- at q(i80)
c
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
      else
         call dmtx(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
      endif
_ELSE
      call dmtx(q(i20),q(i30),q(i80),iky,nocmxa,l1,l1)
_ENDIF
      if (out) then
         write (iwr,9228)
         write (iwr,9108)
         call prtril(q(i20),l1)
      endif
c
c     compute RMS convergence on alpha density
c
      if (iter.ne.1) then
       call rdedx(q(i30),l2,ibl3pa,idaf)
       dsqa = cvgdens(q(i20),q(i30),l2)
      endif
c
c     ----- save alpha mo*s + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -d- at q(i90)
c
      ndaf = mouta
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i20),q(i90),q(i80),ndaf,l1,l2,
     * ibl3pa,ibl3ea)
      call end_time_period(TP_TEST6)
c
c     ----- now work on beta set -----
c
 321  dsqb = dzero
      if (nb .eq. 0) go to 420
      call start_time_period(TP_TEST4)
      if (odiis) then
        call rdedx(q(i10),l2,iblkhb,num8)
      else
        call rdedx(q(i10),l2,iblkh0b,num8)
      endif
c
      if(oshift)then
         call rdedx(q(i30),l3,ibl3qb,idaf)
      else
         call rdedx(q(i30),l3,ibl3qs,idaf)
      endif
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymao,
     +                      n_irreps, isymmo, ierr )
         if( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
         ierr = 999
      end if
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i50),ih_scr2, l1)
      else
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
      endif
_ELSE
      call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
      diffpb = diffb
      diffb=0.0d0
      ii=noccb+1
      if(ii.le.l0) then
        do 253 i=ii,l0
          ij=iky(i)+i50
          loop=idamax(min(i-1,nocmxb),q(ij),1)
          if(loop.gt.0) then
          dum = dabs(q(ij+loop-1))
           if (dum.gt.diffb) then
             diffb = dum
             if (iter.le.maxcycp) then
              testerb(iter) = diffb
              idomob(iter) = loop
              ivmob(iter) = i
             endif
           endif
          endif
 253    continue
      endif
c
      if(ondiis.and..not.osmear)diffb=diffdb
      dmplim=dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diffb.lt.acurcy).and.(iter.gt.1)
      call end_time_period(TP_TEST4)
      if (maxcyc .eq. 0) go to 322
      if(ocvged)go to 322
      rshift=0.0d0
      if((.not.ondiis.or.olevd).and..not.(uhfatom.and.iter.eq.0))
     *call shiftq(q(i50),nocmxb,0,l0,de,dep,iterv,2,gapb1,ibrk,gapb2)
      rshftb=rshift
      call start_time_period(TP_DIAG)
      m2=2
      diaacc = diffb*5.0d-3
      if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
      if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
      call jacobi_symm(q(i50),iky,l0,q(i30),ilifq,l0,q(i90),
     +     m2,lock,diaacc, isymmo, use_symmetry)
*     call jacobi(q(i50),iky,l0,q(i30),ilifq,l1,q(i90),m2,lock,diaacc)
      call end_time_period(TP_DIAG)
c
c     ----- back-transform the eigenvectors -----
c
c     -v- at q(i30)
c     -d- at q(i41)
c     -q- at q(i10)
c     scratch area at q(i21)
c
      call start_time_period(TP_TEST5)
      if(oshift) then
         call rdedx(q(i10),l3,ibl3qb,idaf)
      else
         call rdedx(q(i10),l3,ibl3qs,idaf)
      endif
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i50),l0,l1,l1)
_ENDIF
      call end_time_period(TP_TEST5)
      ig = mod(iter,igs)
      if ( oshift .and. ig .eq. 0) then
c
c     ----- if vshift is true, reorthogonalize the vectors and
c           thereby the transformation matrix every igs th iteration
c
c     -q- at q(i10)
c     -s- at q(i30)
c     -v- at q(i30)
c
        call start_time_period(TP_ORFOG)
        call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
        if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
           call porth(q(i30),q(i10), num, l0)
        else
           call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
           call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
        endif
_ELSE
        call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
        call orfog(q(i30),q(i30),q(i50),q(i60),iky,ilifq,l0,l1,1)
_ENDIF
        call end_time_period(TP_ORFOG)
      endif
      call start_time_period(TP_TEST6)
      call wrt3(q(i30),l3,iblkqq,num8)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      noccb = nb
      ig = nocmxb
      if (osmear) then 
        esmearb= max(esmear_final,
     &               min(esmearb,esmear_scale*(diffb-acurcy)))
        call fermi_smear(q(i90),q(i81),emrmnb,noccb,nocmxb,esmearb,nb,
     &                   l0,iwr,oprint(60))
      endif
      if(ig.lt.l0) then
_IFN1(civ)      call vsadd(q(ig+i90),1,-rshift,q(ig+i90),1,l0-ig)
_IF1(civ)      call vsav(l0-ig,-rshift,q(ig+i90),q(ig+i90))
      endif
      if (out) then
         write (iwr,9208)
         write (iwr,9128)
         call prev(q(i30),q(i90),lprntb,l1,l1)
      endif
c
c     ----- form beta density matrix -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -o- at q(i81)
c
_IF(ga)
      if(l1.ge.idpmult2 .and. ipiomode .ne. IO_NZ_S)then
         call dmtxp(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
      else
         call dmtx(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
      endif
_ELSE
      call dmtx(q(i20),q(i30),q(i81),iky,nocmxb,l1,l1)
_ENDIF
      if (out) then
         write (iwr,9228)
         write (iwr,9128)
         call prtril(q(i20),l1)
      endif
c
c     compute RMS convergence on beta density
c
      if (iter.ne.1) then
       call rdedx(q(i30),l2,ibl3pb,idaf)
       dsqb = cvgdens(q(i20),q(i30),l2)
      endif
c
      call end_time_period(TP_TEST6)

  420 continue
c
c
      if (iter.ne.1) then
       dsq = (dsqa+dsqb) / dfloat(l2)
       diffdens = dsqrt(dsq)
       if (oprind) write(iwr,9876) diffdens
      endif
c
      if (nb .eq. 0) call setz(q(i30),q(i20),q(i90),l1,l2,l3)
c
c     ----- save mo's + density + orbital energies -----
c
c     -v- at q(i30)
c     -d- at q(i20)
c     -d- at q(i90)
c
      ndaf = moutb
c
      call rdedx(q(i30),l3,iblkqq,num8)
      call scfsav(q(i30),q(i20),q(i90),q(i81),ndaf,l1,l2,
     *ibl3pb,ibl3eb)

322   continue
      if(numdis.ne.num8) then
        call wrt3(en,nw1d,ibldis,numdis)
        call wrt3s(st,nw2d,numdis)
      endif
c
      call rdedx(q(i40),l2,ibl3pa,idaf)
      call rdedx(q(i50),l2,ibl3pb,idaf)
      tim1=cpulft(1)
      delt=tim1-tim0
      tim0=tim1
      diffpp = diffp
      diffp = diff
      diff = dmax1(diffa,diffb)
c
c     monitor static tester with diis on - at the moment flag
c     convergence .... (removing diis is an alternative)
c     note that the thresholds for dum and dume below are
c     subject to revision.
c
      dum = dabs(diffp-diff)
      dume = dabs(etot-etot0)
      if (oprind) write(iwr,7878) dum, dume
      if(dum.lt.1.0d-6.and.dume.lt.1.0d-8) then
       iterstat = iterstat + 1
      else
       iterstat = 0
      endif
c
      if(iterstat.gt.8.and.ondiis) then
c
c     trouble .. both energy and tester are stationary but are above
c     the convergence threshold ... time to quit ...
c    
        write(iwr,9825)
        nsti(1)=0
        ondiis=.false.
        iterstat = 0
        ofalsec = .true.
      endif
c
c     ----- print and check convergence behavior -----
c
      otrigp = maxcyc-iter.lt.10
      if((nprint.ne.-5).or.otrigp) then
       if(iter.gt.0) then
         if(odebug(31)) then
          write (iwr,9028)
          write (iwr,9248) iter,kcount,etot,ehf,de,diff,rshfta,rshftb,
     +    damp,derr,delt,tim1
         else
          write (iwr,9029)
          write (iwr,9249) iter,kcount,etot,ehf,de,diff,rshfta,rshftb,
     +    damp,derr
         endif
       endif
      endif
_IF(parallel)
      else
c    stuff for other screened-out nodes
c    this will only be necessary if these nodes are joining in the
c    parallel diag
c    in future extra sections will be needed to enable parallel 
c    diis with ipiomode = IO_NZ_S
       iter = iter + 1
       if(irest.gt.0) then
        write(iwr,3010)irest
        go to 580
       endif
      endif
      if(ipiomode .eq. IO_NZ_S)then
         label="brdcst maxcyc"
         lscf = 20+12/nav
         call pg_brdcst(7125,maxcyc,lscf*8,0)
         label=" "
         oswed3(4)=.true.
         oswed3(8)=.true.
      endif
      call pg_synch(5555)
_ENDIF
      dmplim = dmax1(2.0d0,dmpcut)
      ocvged = (damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
c
c     now flag covergence if STATIC tester has been encountered
      if (ofalsec) then
       ocvged = .true.
       write(iwr,9191) diff
       oprind = .true.
      endif
      if (maxcyc .eq. 0) go to 500
      if (ocvged) go to 500
      damp = dzero
c
      call timit(0)
      call timrem(tlefti)
c
      if (iter.eq.0) go to 180
      if (tlefti .le. skale*(tlefts-tlefti)/iter) then
c
c     ----- exit as time limit exceeded
c
         irest=3
         write(iwr,9367)
         call texit(0,irest)
      else
         if (iter .lt. maxcyc) go to 180
          if(omaxcyc) then
            write (iwr,9329)
            ocvged = .true.
          else
            write (iwr,9328)
c        write out error description to xml/punchfile
         call blkerror('excessive number of SCF iterations',0)
            etot = dzero
            irest=3
            ehf = -en
          endif
      endif
c
c     ----- print actual value of energy -----
c
  500 lock=lockt
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr('Memory failure in duhfop')
         endif
      endif
_ENDIF
      if (ocvged.and.nprint.ne.-5) write (iwr,9268)
      write(iwr,9373) iter,tim1,charwall()
_IF(ccpdft)
      if (odft) then
         call CD_get_dftresults(npts,aelec,belec)
         aelec = aelec + belec
         write(iwr,9374) npts,aelec,(aelec-ne)/ne,edft
      endif
_ENDIF
c
c     At this point the electronic energy has been extrapolated to
c     zero Kelvin. For any gradients to be consistent with the
c     energy we need to compute the finite temperature energies.
c
      ehf  = ehf + 0.5d0*emrmn
      if (etot.lt.dzero) etot = etot + 0.5d0*emrmn
c
      if(omodel) then
       write (iwr,9371) ehf,en,etot,diffdens
      else if (ovir)then
       write (iwr,9365) ek,ehf,en,etot,vir,
     +                  diffdens
      else
       write (iwr,9368) ehf,en,etot,diffdens
      endif
_IF(drf)
       if(.not.oreact) then
_ENDIF
       if(nprint.ne.-5.or.oprind.or.otrigp) then
       idum = min(iter,200)
       write(iwr,9199)
       do loop =1,idum
       if (logtest(loop)) then
       write(iwr,9298)loop,testera(loop), idomoa(loop), ivmoa(loop),
     +                     testerb(loop), idomob(loop), ivmob(loop),
     +                     testerda(loop), idomoda(loop), ivmoda(loop),
     +                     testerdb(loop), idomodb(loop), ivmodb(loop)
       else
       write(iwr,9198)loop,testera(loop), idomoa(loop), ivmoa(loop),
     +                     testerb(loop), idomob(loop), ivmob(loop),
     +                     testerda(loop), idomoda(loop), ivmoda(loop),
     +                     testerdb(loop), idomodb(loop), ivmodb(loop)
       endif
       enddo
       write(iwr,9197)
       endif
_IF(drf)
       endif
_ENDIF
_IF(parallel)
      if(omaster) then
        if(ipiomode .eq. IO_NZ_S)then
           oswed3(4) = .false.
           oswed3(8) = .false.
        endif
_ENDIF
_IF(nbo)
c
c save a copy of the alpha- and beta-fock matrices on the dumpfile
c  = section isect(42) = for possible nbo analysis
c
      len42 = lensec(nx)
      len42 = len42 + len42
      m = 42
      call secput(isect(42),m,len42,iblk42)
      call rdedx(q(i10),l2,iblkha,num8)
      call rdedx(q(i20),l2,iblkhb,num8)
      call wrt3(q(i10),nx,iblk42,idaf)
      call wrt3s(q(i20),nx,idaf)
      lds(isect(42)) = nx + nx
_ENDIF
c
c     ----- check spin state -----
c
c     - s- at q(i10)
c     -da- at q(i40)
c     -db- at q(i50)
c     scratch area at q(i20)
c
      call spin(sz,s2,q(i40),q(i50),q(i10),q(i20),q(i30),iky,l2)
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
_IF(parallel)
        if(ipiomode .eq. IO_NZ_S)then
           oswed3(4) = .true.
           oswed3(8) = .true.
        endif
      endif
_ENDIF
      if (numdis.eq.num8) numdis = 0
      if (out .or. nprint .eq. -5) go to 560
_IF(parallel)
      if(omaster) then
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
_ENDIF
      write (iwr,9108)
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call rdedx(q(i30),l2,ibl3pa,idaf)
      call rdedx(q(i21),l1,ibl3ea,idaf)
_IF(drf)
cahv --- DRF extension ---
      if (field .ne. ' ') then
        ehfhh = ehf
        enhh = en
        etothh = etot
        call dawrit(idafh,ioda,q(i30),l2,16,navh)
      endif
cahv
_ENDIF
      call analmo(q(i10),q(i21),q(i80),ilifq,l0,l1)
      call tdown(q(i10),ilifq,q(i10),ilifq,l0)
      if(otsym) 
     +  call symass(q(i10),q(i21),q(i80),q)
      if(.not.oprint(25))then
         write (iwr,9208)
         call prev(q(i10),q(i21),lprnta,l1,l1)
         if(out) then
            write (iwr,9228)
            call prtril(q(i30),l1)
         endif
      endif
      if (nb .gt. 0) then
c
c     -v- at q(i10)
c     -d- at q(i30)
c     -d- at q(i21)
c
        write (iwr,9128)
        call rdedx(q(i10),l3,ibl3qb,idaf)
        call rdedx(q(i30),l2,ibl3pb,idaf)
        call rdedx(q(i21),l1,ibl3eb,idaf)
_IF(drf)
cahv --- DRF extension ---
      if (field .ne. ' ') then
        call dawrit(idafh,ioda,q(i30),l2,20,navh)
      endif
cahv
_ENDIF
        call analmo(q(i10),q(i21),q(i81),ilifq,l0,l1)
        call tdown(q(i10),ilifq,q(i10),ilifq,l0)
        if(otsym) 
     +    call symass(q(i10),q(i21),q(i81),q)
        if(.not.oprint(25))then
           write (iwr,9208)
           call prev(q(i10),q(i21),lprntb,l1,l1)
           if (out) then
              write (iwr,9228)
              call prtril(q(i30),l1)
           endif
        endif
      endif
_IF(parallel)
      else
c...  other roots
c...  allow for section created in symass
c...  save results of the symmetry assignment into ed3
c...  assuming this has been called
       if (otsym) then
        isymsc = isect(499)
        call secput(isymsc,isymtp,1+lensec((num-1)/nav+1)
     +                             +lensec(num),iblnum)
       endif
      endif
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
  560 continue
_IF(ga)
      if ( zruntp .ne. 'saddle'   .and. zruntp.ne.'optxyz'
     + .and. zruntp .ne. 'optimize' ) then
          call destroy_diis_storage(.true.)
      endif
_ENDIF
      if (mpunch .gt. 0) then
c
c     ----- punch the occupied orbitals
c
c     -v- at q(i10)
c
         call rdedx(q(i10),l3,ibl3qa,idaf)
         call pusql(q(i10),na,l1,l1)
         if (nb .gt. 0) then
            call rdedx(q(i10),l3,ibl3qb,idaf)
            call pusql(q(i10),nb,l1,l1)
         endif
      endif
  580 continue
      if (ocvged) irest = 0
      if(nprint.ne.-5) then
        cpu=cpulft(1)
        write(iwr,7777)cpu,charwall()
      endif
c
c     ----- reset core memory -----
c
_IF(drf)
      if ((field(5:) .eq. 'scf') .and. (intdrf .eq. 0)) then
        call gmem_free(i190)
        call gmem_free(i180)
        call gmem_free(i170)
        if (itwoeps .eq. 1) call gmem_free(i160)
        call gmem_free(i230)
        call gmem_free(i220)
        call gmem_free(i210)
        call gmem_free(i200)
        call gmem_free(i150)
        call gmem_free(i140)
        call gmem_free(i130)
        call gmem_free(i120)
        call gmem_free(i110)
      endif
_ENDIF
      call gmem_free(irdmat)
c ***
c *** set tolitr so that if reenter in geom opt with old vecs
c *** start off with full accuracy
c ***
      tolitr(1)=tolitr(3)
      tolitr(2)=tolitr(3)
c
      accdi1 = accdin
      call copq_again(mouta)
      call copq_again2(moutb)
c     if (oreset) then
c
c     reset maxcyc and level shifters to reasonable values
c     if set  dynamically (for geom. optimisation etc).
c
c      oreset = .false.
c      gapa1 = 3.0d0
c      gapb1 = 3.0d0
c      gapa2 = gapa1
c      gapb2 = gapb1
c      maxcyc = maxcyc / 2
c     endif
      return
c
c9831 format(/15x,'**** increase level shifter to ',f8.2,' ****'/)
 9876 format(15x,' **** convergence on density = ', f20.10)
 9828 format(1x,'**** diis off **** at cycle',i4,' due to energy rise')
 9829 format(1x,'**** diis off **** at cycle',i4,' due to tester rise')
 9827 format(1x,'rising diff - de = ', f20.10)
 9826 format(/1x,'diis initiated at tester = ', f15.9)
 9825 format(/1x,'*** STATIC tester .. flag convergence  ****')
 9824 format(1x,'diffda,diffdpa,diffdb,diffdpb = ',4f15.9)
 9823 format(1x,'diffd, acurcy*diisf = ', 2f15.9)
 7878 format(1x,'tester difference = ', f15.9,
     +       1x,'energy difference = ', f15.9)
 9191 format(/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,'+ convergence monitored at tester = ', f15.10,'  +'/
     + 10x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
 9199 format(//1x,104('+')/
     + 1x,'CONVERGENCE / TESTER ANALYSIS',32x,'From DIIS'/
     + 1x,104('+')/
     + 1x, 'Iter.', '   Alpha    (somo/vmo)         Beta (somo/vmo)',
     + 5x,          '   Alpha    (somo/vmo)     Beta (somo/vmo)'/)
 9198 format(1x,i5,2(f11.7,' (',2i4,')',2x), '(*)', 2(f11.7,' (',
     +       2i4,')') )
 9298 format(1x,i5,2(f11.7,' (',2i4,')',2x), 3x, 2(f11.7,' (',2i4,')'),
     +       ' (*)' )
 9197 format(1x,104('+')//)
 2990 format(/1x,'restarting direct-uhf in integrals :',
     *  'restart parameter =',i3/
     *  1x,'section 171 at block ',i5)
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))
 3010 format(' restart parameter after dhstart ',i3)
 9008 format(/' ----- nuclear energy ----- = ',f20.12)
 9009 format(/1x,'use symmetry adapted jacobi diagonalisation')
 9367 format(/10x,26('*')/
     *10x,'*** warning ***'/
     *10x,'uhf has not converged '/
     *10x,'this job must be restarted'/
     *10x,26('*')/)
 7777 format(//
     *' end of direct-uhf at ',f12.2,' seconds',a10,' wall',
     *1x,104('-')//)
 9030 format(/15x,'convergence data'/15x,16('=')//
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criterion        =1.0e-',i2/
     +     ' punch out option             = ',i6)
 9028 format(/1x,124('=')/
     * 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',5x,'virtual',1x,'damping',12x,'diis',
     * 3x,'del(t)',5x,'time'/
     * 18x,'energy',10x,'energy',38x,'shift'/
     * 72x,'  -a- ','  -b- ' /1x,124('='))
 9029 format(/1x,106('=')/
     + 3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     + 9x,'tester',5x,'virtual',1x,'damping',12x,'diis'/
     + 18x,'energy',10x,'energy',38x,'shift'/
     + 72x,'  -a- ','  -b- ' /1x,106('='))
 9048 format(20x,11('-')/20x,'fock matrix'/20x,11('-'))
 9068 format(/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23('-'))
 9108 format(//1x,104('-')//
     *50x,'----- alpha set -----')
 9128 format(//1x,104('-')//
     *50x,'----- beta set -----')
 9188 format(1x,'core assignment'/' i10, i20, i21,',
     +     ' i30, i31, i40, i50, i60, i70, i80, i81, i90 = '/13i8/
     +     1x,'last = ',i8)
 9208 format(//50x,12('=')/50x,'eigenvectors'/50x,12('='))
 9228 format(/20x,14('-')/20x,'density matrix'/20x,14('-'))
 9248 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,2f6.2,f8.3,f16.9,
     1       f9.2)
 9249 format(i3,1x,i3,1x,g20.14,g16.10,g12.5,g15.5,2f6.2,f8.3,f16.9)
 9268 format(/10x,15('-')/10x,'energy converged'/10x,15('-'))
 9328 format(/10x,30('-')/10x,'excessive number of iterations'/
     +        10x,30('-'))
 9329 format(/10x,30('-')/
     +        10x,'excessive number of iterations'/
     +        10x,'but flag SCF convergence'/
     +        10x,30('-'))
 9348 format(//1x,104('-')//
     *40x,22('*')/40x,'direct uhf calculation'/40x,22('*') )
 9365 format(
     +10x,'kinetic energy             ',f20.10/
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10,5x,'virial ratio(v/2e)'/
     +10x,'total energy               ',f20.10,8x,f10.7/
     +10x,'convergence on density     ',f20.10)
 9368 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total energy               ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9371 format(
     +10x,'electronic energy          ',f20.10/
     +10x,'nuclear energy             ',f20.10/
     +10x,'total qm energy            ',f20.10/
     +10x,'convergence on density     ',f20.10)
 9373 format(/10x,14('-')/10x,'final energies   after',i4,' cycles at '
     +,f12.2,' seconds',2x,a10,' wall',/10x,14('-')//)
 9374 format(
     +10x,'number of quadrature points',i9/
     +10x,'integrated electron count  ',f20.10,5x,
     +    'relative error ',e10.2/
     +10x,'XC energy                  ',f20.10)
 2970 format(' restart parameter in direct-uhf ',i2)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
9407  format(1x,'EHF1:               ',f16.8/
     +       1x,'EHF2:               ',f16.8/
     +       1x,'==============================')
      end
      subroutine mem_dhstaru(imem)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     Returns the estimated amount of memory needed by dhstaru
c     in words
c
c     In:
c         common/restar/nopk
c         common/restar/intg76
c         common/restar/nindmx
c         common/zorac/oint_zora
c         common/sortp/oschw
c         common/nshel/nshell
c         common/mapper/ikyp
c         common/sortp/osortp
c         common/atmblk/num2ep
c         common/ijlab/opfbas
c         common/ijlab/opgbas
c         common/infoa/num
c         common/symtry/nw196
c
c     Out:
c         imem : the amount of memory used by dhstaru
c
c     Functions:
c         lenwrd()
c
INCLUDE(common/sizes)
INCLUDE(common/restar)
INCLUDE(common/zorac)
INCLUDE(common/sortp)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
INCLUDE(common/ijlab)
INCLUDE(common/infoa)
INCLUDE(common/symtry)
      external lenwrd
_IF(ma)
      parameter (MT_DBL = 1013)
      external MA_sizeof_overhead
_ENDIF
      integer iohdmem
_IF(ma)
      iohdmem = MA_sizeof_overhead(MT_DBL)
_ELSEIF(dynamic_memory)
      iohdmem = 0
_ELSE
      iohdmem = 4
_ENDIF
      imem = 0
      icmem = 0
      if (nopk.ne.1) then
         imem_gout = 151875+iohdmem
      else
         imem_gout =  50625+iohdmem
      endif
      if (oint_zora) then
         imem_zora = nwiso_z+iohdmem
      else
         imem_zora = nw196(5)+iohdmem
      endif
      icmem = icmem + imem_gout + imem_zora
      imem  = max(imem,icmem)
      if (oschw) then
         imem_schw = ikyp(nshell) + iohdmem
         icmem     = icmem + imem_schw
         imem      = max(imem,icmem)
         icmem     = icmem - imem_schw
      endif
      if (nopk.ne.1) then
         if (intg76.eq.0.and.nindmx.ne.0) then
            if (oschw) then
c              call pkints
               imem_schw = ikyp(nshell) + iohdmem
               if (osortp) then
                  n__gbf = num*6
                  if (opfbas) n__gbf = num*10
                  if (opgbas) n__gbf = num*15
                  imem_sortp = n__gbf*num2ep+(n__gbf*num2ep+1)/lenwrd()
     +                       + (n__gbf+1)/lenwrd() + iohdmem
               else
                  imem_sortp = 0
               endif
               icmem = icmem + imem_schw + imem_sortp
               imem  = max(imem,icmem)
               icmem = icmem - imem_schw - imem_sortp
            else
c              call pkinta
               if (osortp) then
                  n__gbf = num*6
                  if (opfbas) n__gbf = num*10
                  if (opgbas) n__gbf = num*15
                  imem_sortp = n__gbf*num2ep+(n__gbf*num2ep+1)/lenwrd()
     +                       + (n__gbf+1)/lenwrd() + iohdmem
               else
                  imem_sortp = 0
               endif
               icmem = icmem + imem_sortp
               imem  = max(imem,icmem)
               icmem = icmem - imem_sortp
            endif
         else
c           call pkin70
            if (oschw) then
               imem_schw = ikyp(nshell) + iohdmem
            else
               imem_schw = 0
            endif
            if (osortp) then
               imem_sortp = 4*num*num2ep+(4*num*num2ep+1)/lenwrd()
     +                    + (4*num+1)/lenwrd() + iohdmem
            else 
               imem_sortp = 0
            endif
            icmem = icmem + imem_schw
            icmem = icmem + imem_sortp
            imem  = max(imem,icmem)
            icmem = icmem - imem_sortp
            icmem = icmem - imem_schw
         endif
      else
         if (oschw) then
            imem_schw = ikyp(nshell) + iohdmem
         else
            imem_schw = 0
         endif
         icmem = icmem + imem_schw
         imem  = max(imem,icmem)
         icmem = icmem - imem_schw
      endif
      icmem = icmem - imem_zora - imem_gout
      if (icmem.ne.0) then
         write(*,*)'*** MEM_DHSTARU icmem = ',icmem
         call caserr('*** mem_dhstaru: Fatal memory error')
      endif
      return
      end
      subroutine dhstaru(q,da,db,fa,fb,prefac,rdmat,irest)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension da(*),fa(*),db(*),fb(*),q(*)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/maxlen)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/statis)
INCLUDE(common/timeperiods)

      data pt5/0.5d0/
c     data two/2.0d0/
c
      inull = igmem_null()
      if(irest.ne.1) then
        call vclr(fa,1,nx)
        call vclr(fb,1,nx)
      endif
c
c     construct fock matrices using (Pa+Pb) and (Pa-Pb)
c
      do 100 m = 1,nx
      duma=da(m)
      dumb=db(m)
      da(m)=duma+dumb
 100  db(m)=duma-dumb
c
c     scale diagonal of density MEs
c
      do 1   m=1,num
      nij=ikyp(m)
      da(nij)=da(nij)*pt5
   1  db(nij)=db(nij)*pt5
c
      call start_time_period(TP_DHSTAR)
      call timana(5)
      call jandk(zscftp,q,fa,fb,q(inull),da,db,prefac,rdmat)
      call end_time_period(TP_DHSTAR)
      call cpuwal(begin,ebegin)
c
      do 700 i=1,nx
      duma=fa(i)
      dumb=fb(i)
      fa(i)=duma-dumb
 700  fb(i)=duma+dumb
c
      do   6 m = 1,num
      nij = ikyp(m)
      da(nij) = da(nij)+da(nij)
    6 db(nij) = db(nij)+db(nij)
      do 420 m = 1,nx
      duma = da(m)
      dumb = db(m)
      da(m) = (duma+dumb)*pt5
 420  db(m) = (duma-dumb)*pt5
c
      call dscal(nx,pt5,fa,1)
      call dscal(nx,pt5,fb,1)
_IF(parallel)
c***   ***node-MPP***
c...    gather fock-matrices and distribute (watch double size)
      call start_time_period(TP_DHSTAR_GOP)
      call pg_dgop(700,fa,nx,'+')
      call pg_dgop(701,fb,nx,'+')
      call end_time_period(TP_DHSTAR_GOP)
c***   ***node-MPP***
_ENDIF
c
      return
      end
      subroutine drhfgvb(q,ogrhf)
c
c     ----- direct grhf and gvb program driver -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/tran)
INCLUDE(common/scra7)
INCLUDE(common/timez)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/scfwfn)
INCLUDE(common/cslosc)
INCLUDE(common/field)
INCLUDE(common/zorac)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      common/blkcore/correv(512),array(10)
      dimension q(*),o1e(6)
c
      call check_feature('drhfgvb')
c
c     ----- set up length values -----
c
c     nav = lenwrd()
      l1 = num
      l2 = nx
      l3 = l1*l1
      lentri = nx
      nshtri = nshell*(nshell+1)/2
c
c     ----- set up grhf-gvb data from /scfwfn/ -----
c
      call gvbset(nl2max,l1)
c
c     ----- allocate dynamic storage space -----
c
c     -rdmat-  at  q(irdmat)
c     -prefa-  at  q(iprefa)
c     -e-      at  q(i10)
c     -trans-  at  q(i20)
c     -q   -   at  q(i50)
c     -dens-   at  q(i60)
c     -xx,eig- at  q(i40)
c     -coul-   at  q(i70)
c     -exch-   at  q(i80)
c
c     first determine total memory available, and use
c     this in deriving data storage, no of passes of
c     fock build etc. Free allocated storage, and reasssign
c     with amount actually needed
c
      nmaxly  = igmem_max_memory()
c
      mc = 10*l2+2*nshtri
      left = nmaxly - mc 
      if (left .lt. 1) left = 0
      nwscmu = left + 6*l2
      nwscmu = (nwscmu/l2)*l2
      nwscmu = 2*((nwscmu/l2)/3) + ((nwscmu/l2)/3)
      if (nwscmu .ge. nl2max) then
         nwscmu = nl2max
         nl2max = 0
      endif
      nl2max = nl2max - nwscmu
      nl2max = 3*((nl2max+2)/3)
      nwscmu = nwscmu*l2
      nwscmu = max(nwscmu,6*l2)
      nwscm4 = ((nwscmu/l2)/3)*l2
      nwscm4 = max(nwscm4,l1+l3)
c
c     determine total requirements
c
      irdmat = 0
      iprefa=irdmat+nshtri
      i10  = iprefa+nshtri
      i20  = i10+l1
      i30  = i10+l2
      i40  = i30+l2
      i50  = i40+l1
      i60  = i50+l3
      idmat= i60
      i70  = i60+nwscm4
      ifock= i70
      i80  = i70+nwscm4
      iexch= i80
      last = i80+nwscm4
c
      if (ofield) then
         i51  = i40+l2
         i61  = i51+l2
         i71  = i61+l2
         lastf  = i71+l2
         last   = max(last,lastf)
      endif
c
c     ----- get core memory -----
c
      irdmat = igmem_alloc(last)
      iprefa = irdmat+nshtri
      i10  = iprefa+nshtri
      i20  = i10+l1
      i30  = i10+l2
      i40  = i30+l2
      i50  = i40+l1
      i60  = i50+l3
      idmat= i60
      i70  = i60+nwscm4
      ifock= i70
      i80  = i70+nwscm4
      iexch= i80
      last = i80+nwscm4
c
      if (ofield) then
         i51  = i40+l2
         i61  = i51+l2
         i71  = i61+l2
         lastf  = i71+l2
         last   = max(last,lastf)
      endif
      if (nprint .eq. 5) write (iwr,9008)irdmat,iprefa, i10,
     +    i20,i30,i40,i50 ,i60,i70,i80, last
c
c ----- first restore s,t,f from section 192 and store on ed7
c ----- are field values to be introduced
c
      if(ofield) then
c
       do loop =1,6
        o1e(loop) = .true.
       enddo
       call getmat(q(i10),q(i30 ),q(i40),q(i51),q(i61),q(i71),
     * array,num,o1e,ionsec)
       do 220 loop=1,l2
       q(i40+loop-1) = q(i40+loop-1) - fieldx * q(i51+loop-1)
     *                               - fieldy * q(i61+loop-1)
     *                               - fieldz * q(i71+loop-1)
 220   continue
      else
c
       do loop =1,3
        o1e(loop) = .true.
        o1e(loop+3) = .false.
       enddo
       call getmat(q(i10),q(i30 ),q(i40),q(i40),q(i40),q(i40),
     * array,num,o1e,ionsec)
      endif
c
      call wrt3(q(i10),l2,ibl7s,num8)
      call wrt3(q(i30),l2,ibl7t,num8)
      call wrt3(q(i40),l2,ibl7f,num8)
c
c ----- transform s-matrix
c
      call tranp(q(i10),q(i30 ))
      call wrt3(q(i30),l2,ibl7st,num8)
      call dgvbitr(q,q(irdmat),q(iprefa),q(i20),q(i50),q(idmat),
     +  q(ifock),q(iexch),q(i10),q(i40),nwscmu,l1,l2,l3)
c
c     ----- reset memory -----
c
      call gmem_free(irdmat)
c
      if(nprint.eq.-5)return
      cpu=cpulft(1)
      write(iwr,8888)cpu
c
c ***
c *** set tolitr so that if reenter in geom opt with old vecs
c *** start off with full accuracy
c ***
      tolitr(1)=tolitr(3)
      tolitr(2)=tolitr(3)
      return
 8888 format(/' end of direct grhf-gvb scf at ',f12.2,' seconds ',
     +  a10,' wall',//1x,104('-')/)
 9008 format(' core assignment'/
     +' irdmat, iprefa,    i10,   i20,   i30,   i40,   i50,'
     +     , '   i60,   i70,   i80'      /i8,9i7 /' last = ',
     +     i7)
      end
      subroutine dgvbitr(core,rdmat,prefac,vect,temp,
     +  dens,coul,exch,e,pp,nwscmu,l1,l2,l3)
c
c     ----- direct grhf-gvb driver -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension vect(l1,*),temp(l1,*),dens(*),coul(*),exch(*)
      dimension e(*),rdmat(*),prefac(*),core(*)
      dimension pp(*)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scfwfn)
INCLUDE(common/dm)
INCLUDE(common/timez)
INCLUDE(common/restar)
INCLUDE(common/scfopt)
INCLUDE(common/infoa)
INCLUDE(common/harmon)
INCLUDE(common/atmol3)
INCLUDE(common/machin)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/runlab)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/timeperiods)
cgdf
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
c
      common/diisd/st(210),ccc(20),rrr(19),derror,scale(20),iposit(20),
     +        nstore,mp,ondiis,nspaca,
     +        stb(210),cccb(20),rrrb(19),derrb,scaleb(20),iposb(24),
     +        ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     +        nsytyp,nsos,
     +        tti(25),e1i(25),asym,
     +        cilow(12),jpair,kone,ktwo,kcorb(2,12),
     +        iojk(49),ioham(25),iojkao(49)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
_IF(ga)
      character*1 xn,xt
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
c
c  i/o vs message passing logic controlled by ipiomode
c
c  IO_NZ_S  one copy of ed3 and ed7 is maintained by node 0
c           in this case, to minimise comms, some work is performed
c           only on node 0. parallel diag will work (with extra
c           brdcsts) but  orfog,diis,mult2 wont
c
c  IO_NZ    one copy of ed3 and ed7 is maintained by node 0
c           all nodes execute all sections
c
c  IO_A     all nodes use private copies of ed3, ed7, so all nodes
c           run through all code section
c
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
INCLUDE(common/fsymas)
c ... for dummy symass section
      parameter(isymtp=99)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      data zscf/'scf'/
      data done,two,ten/1.0d0,2.0d0,1.0d+01/
      data pt2,twopt2/0.2d0,2.2d0/
      data m23,m20/23,20/
      data dmptlc/1.0d-02/
_IF(ga)
      data xn,xt/'n','t'/
_ENDIF
c
      out = nprint .eq. 5
      outon = nprint .ne. -5
      nav = lenwrd()
c
_IF(parallel)
c
c - this should be default anyway -
c
      if(ipiomode .eq. IO_NZ .or. ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
c
c master nodes execute whole code
c
      omaster = .true.
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
*     write(6,*)'master',ipg_nodeid(),omaster
_ENDIF
      call start_time_period(TP_TEST1)
      skale=2.0d0
      lprnt = l1
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
      if(.not.oprint(20))lprnt = min(norb+5,l1)
      if(zruntp.eq.zscf)skale=1.1d0
      if(outon)write (iwr,9008)
_IF(ccpdft)
      if(CD_active())then
       call caserr('DFT not available for GRHF-GVB, use UHF')
      endif
_ENDIF
      tim=0.0d0
      en=enucf(nat,czan,c)
_IF(ga)
c
c  check/create GA storage, load S
c
      if(ipiomode .ne. IO_NZ_S)then
         call rdedx(vect,l2,ibl7s,num8)
         call rdedx(dens,l2,ibl7st,num8)
         call declare_diis_storage(num,.false.)
         call init_diis_storage(num,vect,dens)
      endif
_ELSE
      call rdedx(dens,l2,ibl7st,num8)
_ENDIF
      call rdedx(vect,l3,ibl3qa,idaf)
      call vclr(e,1,norb)
      if(out) then
        call prev(vect,e,lprnt,l1,l1)
        call prtri(dens,l1)
      endif
c
c     ----- initialize individual timers for the five sections
c           of grhf-gvb-----
c
      timed = 0.0d0
      timef = 0.0d0
      timit = 0.0d0
      timeg = 0.0d0
      timem = 0.0d0
      timeo = 0.0d0
c 
      dlnmxd=0.0d0
      dlntol=tolitr(3)
      m171t=171
      len171=lensec(l2)
      nshtri=ikyp(nshell)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      iof171(4)=iof171(1)+len171*(nham+nham)
      iof171(5)=iof171(2)+len171*nham
      iof171(6)=iof171(3)+len171*(nham+nham)
c
      if(irest.ne.0) write(iwr,2970)irest
      if(irest.ne.1.and.irest.ne.2) then
        call rdmake(prefac)
        ln171=nshblk*2
        if(odelta) then
         ln171=6*nham*len171+ln171
        endif
        call secput(isect(471),m171t,ln171,ibl171)
        if(outon) then
         write(iwr,2980) ibl171,ln171
        endif
         call wrt3(prefac,nshtri,ibl171,idaf)
         if(odelta) then
          call zer171(coul,l2,6*nham,ibl171+iof171(3),len171,idaf)
         endif
         call clredx
      endif
c
c     ---- initialize some constants for convergence control -----
c
      iprint=nprint
      if (nconv .le. 0) nconv = 5
      acurcy = ten**(-nconv)
c
c     is delta density ever to be invoked? if not, we can
c     reduce i/o activity around section 471 ..
c
c     odelta = deltol .gt. dlog(acurcy)
      if(numdis.eq.0) numdis = num8
      if(irest.ne.3.or.numdis.eq.num8) then
         damp = 0.0d0
         damp0 = 0.0d0
c        sqcdf = 0.0d0
         iter = 0
         iterv = 0
         kcount = 0
         rshift = 0.0d0
         diff = 0.0d0
         diffo = diff
         de = 0.0d0
         dep = 0.0d0
         deavg = 0.0d0
         ehf = 0.0d0
         lockt = lock
         ehf0 = 0.0d0
         ek = 0.0d0
         vir = 0.0d0
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0
c
c     ----- initialize disk for extrapolation -----
c
      ndaf = ibl7la
      lnoc = l1*norb
      if (lnoc .lt. l2) lnoc = l2
      lene=lensec(lnoc)
      lene=lene*3
      lenq=lensec(l3)
      ibl7qa=ndaf+lene
      iblko=ibl7qa+lenq
      lenp=lensec(l2)
      ibl7la=iblko+lenp
c
c ----- allocate space on ed7 for j,k,ham
c
      lenh=lensec(l2)
      iojk(1)=ibl7la
      ij=ibl7la+lenh
      nham1=nham-ncores
      if(nham1.gt.0) then
       do 3000 i=1,nham1
       ii=i+i
       iojk(ii)=ij
       iojk(ii+1)=ij+lenh
 3000  ij=ij+lenh+lenh
      endif
      iojkao(1)=ij
      ij=ij+lenh
      if(nham1.gt.0) then
       do 3003 i=1,nham1
       ii=i+i
       iojkao(ii)=ij
       iojkao(ii+1)=ij+lenh
 3003  ij=ij+lenh+lenh
      endif
      do 3002 i=1,nham
      ioham(i)=ij
 3002 ij=ij+lenh
      if(numdis.eq.num8)then
         ndafd=ij
      else
         ndafd=ibldis+2
      endif
      ocvged = .false.
c     orotdm = .false.
c     if (mconv .gt. 7) orotdm = .true.
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
      if (.not.odiis) then
       if (odamph .or. oshift) damp = dmax1(done,dmpcut)
      endif
      call wdisk(vect,vect,vect,ndaf,num8,lnoc)
c
c     ----- compute a set of orthonormal orbitals from the
c           overlap matrix -----
c
      l0=newbas0
      call qmat(core,dens,coul,e,exch,iky,l0,l1,l3,l1,out)
      lprnt = min(lprnt,l0)
      l0 = l1
      ltran = l0
c
c     ----- get initial time remaining for time check purposes -----
c
      call timrem(tlefts)
      call end_time_period(TP_TEST1)
c
c     ----- if iter .eq. 0 , skip the convergence routines -----
c
      timit=cpulft(1)
 100  continue
      call start_time_period(TP_TEST2)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
      if (iter .gt. 0) then
        if(iter.eq.maxcyc-1) nprint=15
c
c     ----- damp and extrapolate the eigenvectors if necessary -----
c
        if (iter .eq. 1) deavg = 0.0d0
        if (iter .eq. 2) deavg =  dabs(de)
        if (iter .ge. 3) deavg = ( dabs(de)+ dabs(dep)+pt2*deavg)/
     +                             twopt2
        dmptst = acurcy
        if ((diffo-diff) .lt. 0.0d0) dmptst = diffo-diff
        if (iter .gt. 2) call dampd(de,dep,deavg,damp,dmptst,diff,diffp,
     +       dmptlc)
        if (damp .lt. dmpcut) damp = dmpcut
        call extrpd(de,damp,damp0,vect,dens,coul,exch,
     +       l1,lnoc,ndaf,numdis,iterv,1,2)
        ltran=l0
      else
        ltran=norb
      endif
c
c     ----- reorthonormalize if damping, extrapolation or
c           level shifting was done -----
c
      if(iter.eq.0)go to 142
      if(mod(iter,5).eq.0) goto 142
      if(.not.odiis.and. kcount.eq.0) goto 142
      if(damp.ne.0.0d0) goto 142
      goto 160
 142  call ortho1(dens,dens,vect,pp,iky,newbas0,l0,l1,l2,l3,l1)
      ltran=l0
      call rdedx(dens,l3,ibl3qs,idaf)
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,vect,l1)
         call load_ga_from_square(ih_vec,dens,l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(vect,ih_scr2, l1)

      else
          call tfsqc(vect,dens,coul,l0,l1,l1)
      endif
_ELSE
      call tfsqc(vect,dens,coul,l0,l1,l1)
_ENDIF
 160  if (out) call prev(vect,e,lprnt,l1,l1)
c
c     ----- update the contents of the transformation vector
c           file -----
c
      call wrt3(vect,l3,ibl3qa,idaf)
      time0=cpulft(1)
      call tdown(vect,ilifq,vect,ilifq,ltran)
c
c     ----- set up the data for the virial calculation -----
c
      call virset(vect,dens,dens(l2+1),l1,l2,num8,core)
c
c     ----- set up data for direct jkform (reduced dmat etc)
c     ----- assumes density matrices in dens
c     ----- form the j and k fock matrices -----
c
_IF(parallel)
      else
c
c...  other roots ** idle **
c...  
      endif
c... 
      if(ipiomode .eq. IO_NZ_S)then
        call pg_brdcst(7124,ltran,8,0)
        call pg_brdcst(7123,vect,l3*8,0)
        call pg_synch(4444)
        oswed3(4)=.true.
        oswed3(8)=.true.
      endif
      call pg_synch(5555)
c...
_ENDIF
      call end_time_period(TP_TEST2)
      if(irest.eq.1 .or. irest.eq.2) then
        call secget(isect(471),m171t,ibl171)
        write(iwr,2990)irest,ibl171
      endif
_IF(parallel)
c***   ***node-MPP***
      call djkform(core,vect,temp,rdmat,prefac,dens,coul,exch,
     +     nwscmu,npass,ltran,irest,iter,nprint,omaster)
c***   ***node-MPP***
_ELSE
      call djkform(core,vect,temp,rdmat,prefac,dens,coul,exch,
     +     nwscmu,npass,ltran,irest,iter,nprint)
_ENDIF
      timdum=cpulft(1)
      timef = timef +timdum - timit
      timit = timdum
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
      if(irest.ne.0) then
        write(iwr,3010)irest
        if(numdis.ne.num8) then
           call wrt3(en,nw1d,ibldis,numdis)
           call wrt3s(st,nw2d,numdis)
        endif
        return
      endif
c
c     ----- determine the optimal ci coefficients -----
c
      call optci(coul,exch,pp,iky,l1,l2,nprint)
      timdum=cpulft(1)
      timeg = timeg + timdum - timit
      timit = timdum
c
c ------ determine gvb orbital gradients
c
      diffpp=diffp
      diffp=diff
      diff=0.0d0
      call gvbgrd(temp,coul,coul(l2+1),iky,l2)
c
c ------ monitor convergence here
c
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
c
cgdf  not sure if needed: convergence check here seems to be for 
c     skipping some work, when the real check comes later
c
      if (osrso) ocvged = .true.
c
c ---- additional code for diis2 method
c
      if(.not.odiis) goto 180
c
c ------ switch off diis if tester rises again
c
      if(diff.gt.diffp.and.diffp.gt.diffpp) then
*      if(diff.ge.acurcy*1000.0d0.or.iter.eq.0) then
        nstore=0
        goto 180
*      endif
      endif
      if (ocvged) go to 270
      call start_time_period(TP_DIIS)
      call diiso(temp,dens,coul,exch,vect,ndafd,lockt,numdis)
      call end_time_period(TP_DIIS)
c
c ------ read in symmetry adapted vectors
c
      call rdedx(vect,l3,ibl3qa,idaf)
      if (ondiis) then
       call start_time_period(TP_ORFOG)
        njkmat=2*(nham-ncores)+1
        call rdedx(dens,l2,ibl7st,num8)
_IF(ga)
        if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
           call porth(vect,dens, num, newbas0)
        else
           call mult2(vect,coul,dens,l0,l0,l1)
           call orfog(vect,vect,coul,coul(l2+1),iky,ilifq,newbas0,l1,1)
        endif
_ELSE
        call mult2(vect,coul,dens,l0,l0,l1)
        call orfog(vect,vect,coul,coul(l2+1),iky,ilifq,newbas0,l1,1)
_ENDIF
        call wrt3(vect,l3,ibl3qa,idaf)
        call tdown(vect,ilifq,vect,ilifq,newbas0)
        call end_time_period(TP_ORFOG)
c
c ------ transform diis generated jk matrices to mo basis
c
        call start_time_period(TP_TEST3)
        do 163 jk=1,njkmat
        call rdedx(coul,l2,iojkao(jk),num8)
_IF(ga)
        if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
           call load_ga_from_square(ih_vec,vect,l1)
           call load_ga_from_triangle(ih_scr2,coul,l1)
           call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
           call load_triangle_from_ga(dens,ih_scr2, l1)
        else
           call mult2(vect,dens,coul,l1,l1,l1)
        endif
_ELSE
        call mult2(vect,dens,coul,l1,l1,l1)
_ENDIF
 163    call wrt3(dens,l2,iojk(jk)  ,num8)
c
c     ----- determine the optimal ci coefficients -----
c
        call optci(coul,exch,pp,iky,l1,l2,nprint)
        call end_time_period(TP_TEST3)
      else
        call tdown(vect,ilifq,vect,ilifq,l1)
      endif
 180  timdum=cpulft(1)
      timed = timed + timdum - timit
      timit = timdum
c
c     ----- perform orbital mixing among occupied orbitals -----
c
c     ----- do the two by two rotations among the occupied orbitals ----
c
      call start_time_period(TP_TEST4)
      call rdedx(vect,l3,ibl3qa,idaf)
      call rotb(dens,vect,coul(l2+1),coul,iky,l1,l2,gapa1,ibrk,gapa2,
     +          nprint)
      timdum=cpulft(1)
      timem = timem + timdum - timit
      timit = timdum
c
c      ----- compute the fock matrices in the sab basis
c
      call gvbham(dens,dens(l2+1),coul,nwscmu,l2)
      call end_time_period(TP_TEST4)
c
c     ----- perform the ocbse step - mix in virtuals -----
c
      call search(ioham(1),num8)
      ncall=1
      do 260 iham = 1,nham
      call reads(coul,l2,num8)
      if (out) call prtri(coul,l1)
      call start_time_period(TP_DIAG)
      call ocbse(iham,vect,dens,coul,temp,pp,e,de,dep,iterv,ncall,
     1           l1,l2)
      call end_time_period(TP_DIAG)
      ncall=0
  260 continue
      call wrt3(vect,l3,ibl3qa,idaf)
      call tdown(vect,ilifq,vect,ilifq,l1)
 270  continue
      timdum=cpulft(1)
      timeo = timeo + timdum - timit
      timit = timdum
c
c     ----- calculate the reduced density matrices -----
c
      call redden(vect,dens,dens(l2+1),nconf,l1,l2)
c
c     ----- combine the alpha and beta parts into one density
c           matrix. then get the old density matrix and save the
c           new density matrix -----
c
      call vadd(dens,1,dens(l2+1),1,dens,1,l2)
      if (iter .gt. 0) call rdedx(dens(l2+1),l2,iblko,num8)
      call wrt3(dens,l2,iblko,num8)
      call rdedx(vect,l3,ibl3qa,idaf)
      if (out) call prev(vect,e,norb,l1,l1)
      time1=cpulft(1)
      shift2=rshift
      shift1=gapa1
      if(iter.gt.ibrk) shift1=gapa2
c     if (iter .gt. maxcyc) go to 340
      ek = ek+ek
      iter = iter + 1
      dep = de
      de = ehf - ehf0
      ehf0 = ehf
      etot = ehf+en
      delt = time1 - time0
      if (iter.eq.1.and.outon) then
        write(iwr,9048) maxcyc,mconv,nconv,npunch,en
      endif
      diffo = diff
      vir = (etot-ek)/(two*etot)
      if(outon.or.(maxcyc-iter.lt.10)) then
       if (odebug(31)) then
        write (iwr,9028)
        write (iwr,9068) iter,kcount,etot,ehf,de,diff,shift1,shift2,
     +   damp,derror,delt,time1
       else
        write (iwr,9029)
        write (iwr,9069) iter,kcount,etot,ehf,de,diff,shift1,shift2,
     +   damp,derror
       endif
      endif
_IF(parallel)
      else
          iter = iter + 1
      endif
******
      if(ipiomode .eq. IO_NZ_S)then
        lscf = 20+12/nav
        call pg_brdcst(7125,maxcyc,lscf*8,0)
        call pg_brdcst(7126,vect,l3*8,0)
        call pg_synch(4444)
        oswed3(4)=.true.
        oswed3(8)=.true.
      endif
      call pg_synch(5555)
_ENDIF
c
c ------ monitor convergence here
c
      dmplim = dmax1(dmpcut,2.0d0)
      ocvged=(damp.lt.dmplim).and.(diff.lt.acurcy).and.(iter.gt.1)
cgdf
      if (osrso) then
        write(iwr,9512)
 9512 format(
     + /,1x,'****************************************************'
     +,/,1x,'* a single scf iteration has been computed'
     +,/,1x,'* now choosing more powerful second-order method'
     +,/,1x,'* using the Newton-Raphson solver in MASSCF'
     +,/,1x,'****************************************************'
     +,/)
        ocvged = .true.
      end if
c
      if (.not.ocvged) then
c
c     ----- check for time remaining -----
c
        call timrem(tlefti)
        if (tlefti .gt. skale*(tlefts-tlefti)/iter) then
c       sufficient time to usefully continue
          if (iter .lt. maxcyc) go to 100
          write (iwr,9108)
c         write out error description to xml/punchfile
          call blkerror('excessive number of SCF iterations',0)
          etot = 0.0d0
          irest=3
          ehf0 =ehf
          ehf = -en
        else
c     insufficient time to usefully continue
          write (iwr,9088)
          irest=3
          call texit(0,irest)
          etot = 0.0d0
          ehf = -en
        endif
      else
c
c     ----- energy converged -----
c
      if(outon)write (iwr,9128)
c
      endif
c
c     ----- calculate the five time steps per iteration -----
c
      timed = timed/ dfloat(iter)
      timef = timef/ dfloat(iter)
      timeg = timeg/ dfloat(iter)
      timem = timem/ dfloat(iter)
      timeo = timeo/ dfloat(iter)
      if(outon)write(iwr,9148)npass,timef,timeg,timem,
     +                              timeo,timed
c
c     ----- save the final /scfwfn/ -----
c
      m4=lensec(mach(5))
      i =m4+m4
      call secput(isect(504),m23,i,iblk23)
      iblk23 = iblk23 + m4
      call wrt3(cicoef,mach(5),iblk23,idaf)
c
c     ----- get the true eigenvalues -----
c
      do 400 i = 1,nham
      do 380 j = 1,norb
      if (nconf(j) .ne. i) go to 380
      e(j) = two*e(j)*f(i)
  380 continue
  400 continue
c
      l2norb = norb*norb
      len2=lensec(l2norb)
      call secput(iseclg,m20,len2,iblk20)
_IF(parallel)
      if(omaster) then
       if(ipiomode .eq. IO_NZ_S)then
          oswed3(4) = .false.
          oswed3(8) = .false.
       endif
_ENDIF
c
c     ----- save the final eigenvalues and orbitals -----
c
       call start_time_period(TP_TEST5)
       call gvbsav(vect,e,mouta,l1,ibl3ea)
       call rdedx(vect,l3,ibl3qa,idaf)
       call lagrng(l1,norb,l2,coul,dens,vect,iky)
c
c     ----- save the lagrangian multpliers in section iseclg -----
c
       call wrt3(dens,l2norb,iblk20,idaf)
       if(ocvged)irest=0
       nprint=iprint
       if(numdis.eq.num8)numdis=0
       call rdedx(vect,l3,ibl3qa,idaf)
c
       call gvbout(vect,temp,e,l1,nprint)
c
c     canonicalise gvb mos
c
       call canon(core,vect,dens,coul,e,l0,l1,l2,l3,nprint)
       call end_time_period(TP_TEST5)
_IF(parallel)
      else
c...  other roots
c...  allow for sections created in canon (moutb) and symass
c...  save results of the symmetry assignment into ed3
c...  assuming this has been called
      len1=lensec(mach(8))
      len2=lensec(mach(9))
      j=len2+lenq+len1+1
      call secput(moutb,3,j,iblnum)
       if (otsym) then
        isymsc = isect(499)
        call secput(isymsc,isymtp,1+lensec((num-1)/nav+1)
     +                             +lensec(num),iblnum)
       endif
      endif
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
c
      call copq_again(mouta)
      return
c
 9008 format(//40x,31('*')/ 40x,
     *'direct grhf-gvb scf calculation'/40x,31('*'))
 9028 format(/1x,127('=')/
     *3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',4x,'level shift',2x,'damping',11x,'diis',
     * 3x,'del(t)',5x,'time'/
     *18x,'energy',10x,'energy'/1x,127('='))
 9029 format(/1x,109('=')/
     *3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     * 9x,'tester',4x,'level shift',2x,'damping',11x,'diis'/
     *18x,'energy',10x,'energy'/1x,109('='))
 9068 format(1x,i3,1x,i3,3f16.8,f15.8,3f8.3,f15.8,2f9.2)
 9069 format(1x,i3,1x,i3,3f16.8,f15.8,3f8.3,f15.8)
 9048 format( /15x,'convergence data'/15x,16('=')//
     +     ' maximum number of iterations = ',i6/
     +     ' method of convergence        = ',i6/
     +     ' convergence criteria         =1.0e-',i2,/
     +     ' punch out option             = ',i6,//
     +    ' ----- nuclear energy ----- = ',f20.12/)
 9088 format(//10x,26('*')//
     *10x,'*** warning ***'/
     *10x,'scf has not converged yet'/
     *10x,'this job must be restarted'/
     *10x,26('*')//)
 9108 format(/10x,30('-')/10x,'excessive number of iterations'/
     +     10x,30('-'))
 9128 format(/10x,16('-')/10x,'energy converged'/10x,16('-'),/)
 9148 format(/10x,28('-')/
     +        10x,'scf statistics per iteration'/
     +        10x,28('-')/
     +        10x, 'number of integral passes  ',i5/
     +        10x,'j+k  formation time    ',f10.3/
     +        10x,'geminal opt time       ',f10.3/
     +        10x,'mixorb  opt time       ',f10.3/
     +        10x,'ocbse   opt time       ',f10.3/
     +        10x,'orbital grad. and diis ',f10.3)
 2970 format(' restart parameter in direct-gvb ',i2)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
 2990 format(/1x,'restarting direct-gvb in integrals :',
     *  'restart parameter =',i3/
     *  1x,'section 171 at block ',i5)
 3010 format(' restart parameter after djkform',i3)
      end
      subroutine djkform(core,trans,temp,rdmat,prefac,
     + dens,coul,exch,nwscmu,npass,ltran,irest,iter,
_IF(parallel)
     + nprint,omaster)
_ELSE
     + nprint)
_ENDIF
c
c  djkform drives the two electron integral part of the fock matrix
c  formation for direct-RHF/GVB calculations.
c  multiple passes over the integrals are done only if needed. 
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/scfwfn)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/cslosc)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/zorac)
INCLUDE(common/timeperiods)
      common/diisd/sta(270),iposa(24),stb(270),iposb(24),
     + ns(15),melsym(15),nconf(maxorb),msympr(maxorb),norb,nham,
     + nsytyp,nsos,
     + tti(25),e1i(25),asym,
     + cilow(12),jpair,kone,ktwo,kcorb(2,12),
     + iojk(49),ioham(25),iojkao(49)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
_IF(parallel)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
      character*20 label
      common/msglab/label
_ENDIF
c
      dimension trans(num,*),core(*),rdmat(*),prefac(*)
      dimension dens(*),coul(*),exch(*),temp(*)
c
      nshtri = ikyp(nshell)
c
      nscmu = nwscmu/nx
      kham = 1
      npass = 0
      iscmf = 1 + ncores
      nleft = nham
  100 continue
      npass = npass + 1
      nclf = 0
      nscmf = nscmu/3
      if (npass .le . 1) then
         nclf = ncores
         nscmf = (nscmu-ncores*3)/3
         nleft = nleft - ncores
      endif
      if (nleft .le. nscmf) then
         nscmf = nleft
      endif
c
      nleft = nleft - nscmf
      nsheld = nclf + nscmf
c
c     load up the density matrices
c
      ioff = 1
      if (nclf .gt. 0) then
c
c      closed shell part
c
       call dengvb(trans,dens,1,iky,num)
       ioff=ioff+nx
      endif
c
c     open shell part
c
      if (nscmf .gt. 0) then
       ilow = iscmf
       ihi = ilow + nscmf - 1
       do 200 i = ilow,ihi
       call dengvb(trans,dens(ioff),i,iky,num)
 200   ioff = ioff + nx
      endif
c
c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** delta coul,exch , delta dens, last coul,exch, current dens
c *** red. dmat, prefac mat.
c ***
      if(irest.eq.1 .or. irest.eq.2) then
        call rdedx(prefac,nshtri,ibl171,idaf)
        call reads(rdmat,nshtri,idaf)
        if (odelta) then
         call rdnsh(coul,nx,nsheld,ibl171+iof171(3),idaf)
         call rdnshs(exch,nx,nsheld,idaf)
         call rdnshs(dens,nx,nsheld,idaf)
        endif
        dlnmxd=-9999999.0d0
        do 1672 kkk=1,nshtri
            if(rdmat(kkk).le.dlnmxd) goto 1672
            dlnmxd=rdmat(kkk)
1672    continue
_IFN(parallel)
        if(irest.eq.1) goto 6351
        if(irest.eq.2) goto 6352
_ENDIF
      endif
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4)=.false.
         oswed3(8)=.false.
      endif
_ENDIF
      ndzero=nsheld*nx
      call start_time_period(TP_RDMAT)
_IF(parallel)
      if(omaster) then
      if(irest.eq.1) goto 6351
_ENDIF
c     delta coul,exch (nsheld)
      if (odelta) then
        call vclr(coul,1,ndzero)
        call vclr(exch,1,ndzero)
        call wrtnsh(exch,nx,nsheld,ibl171+iof171(3),idaf)
        call wrtnshs(exch,nx,nsheld,idaf)
        call rdnsh(exch,nx,nsheld,ibl171+iof171(6),idaf)
        call wrtnsh(dens,nx,nsheld,ibl171+iof171(6),idaf)
        call vsub(dens,1,exch,1,dens,1,ndzero)
        call wrtnsh(dens,nx,nsheld,ibl171+iof171(4),idaf)
      else
c       call wrtnsh(dens,nx,nsheld,ibl171+iof171(6),idaf)
      endif
      if(nprint.ne.-5) write(iwr,3030)
      if (odelta) then
        call mkrdmt(zscftp,rdmat,dens,nx,nprint)
      endif
c ***
c *** make sure that iter before delta calc is of full accuracy.
c ***
      if( (dlnmxd.lt.deltol) .and. (iter-1.lt.itrtol(2))) then
        dlnmxd=deltol+1.0d0
      endif
      if( (dlnmxd.gt.deltol) ) then
        if(nprint.ne.-5) write(iwr,3000)
        if (odelta) then
           call rdnsh(dens,nx,nsheld,ibl171+iof171(6),idaf)
           call wrtnsh(dens,nx,nsheld,ibl171+iof171(4),idaf)
           call vclr(exch,1,ndzero)
           call wrtnsh(exch,nx,nsheld,ibl171+iof171(5),idaf)
           call wrtnshs(exch,nx,nsheld,idaf)
           call mkrdmt(zscftp,rdmat,dens,nx,nprint)
        else
           call mkrdmt(zscftp,rdmat,dens,nx,nprint)
        endif
        call clredx
      endif
      call wrt3(rdmat,nshtri,ibl171+iof171(2),idaf)
6351  continue
_IF(parallel)
      endif
c
      if(ipiomode .eq. IO_NZ_S)then
c ps: combine to 1 brdcst
         call pg_brdcst(7122,rdmat,nshtri*8,0)
         call pg_brdcst(7124,dlnmxd,8,0)
         label="brdcst idmat"
         call pg_brdcst(7123,dens,ndzero*8,0)
         label=" "
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF(parallel)
      call end_time_period(TP_RDMAT)
c ***
      if(iter.lt.itrtol(1)) then
        dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
        dlntol=tolitr(2)
      else
        dlntol=tolitr(3)
      endif
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(nprint.ne.-5) then
      write(iwr,3020) dlntol
      endif
c
      call dhstarg(core,coul,exch,dens,prefac,rdmat,nclf,nscmf,
     +     irest,nprint)
c
c *** form complete j+K matrices and save them to section 471
c
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
      if(omaster) then
_ENDIF
6352  continue
      if(irest.ne.0) return
      call start_time_period(TP_RDMAT)
      if (odelta) then
         call rdnsh(dens,nx,nsheld,ibl171+iof171(5),idaf)
         call vadd(coul,1,dens,1,coul,1,ndzero)
         call wrtnsh(coul,nx,nsheld,ibl171+iof171(5),idaf)
         iblexc = iposun(idaf)
         call rdnshs(dens,nx,nsheld,idaf)
         call vadd(exch,1,dens,1,exch,1,ndzero)
         call wrtnsh(exch,nx,nsheld,iblexc,idaf)
      else
c        call wrtnsh(coul,nx,nsheld,ibl171+iof171(5),idaf)
c        call wrtnshs(exch,nx,nsheld,idaf)
      endif
      call end_time_period(TP_RDMAT)
c
c     symmetrize the fock matrix and transform to molecular orbital
c     operators
c
c     ----- closed shell part -----
c
      call start_time_period(TP_TEST6)
      ioff = 1
      if (npass .le. 1 ) then
       call rdedx(dens,nx,ibl7f,num8)
c
c...  add zora corrections if required
c
      if (ozora) call zora(core,core,dens,'read')
c
        if (nco .le. 0) then
           call dcopy(nx,dens,1,temp,1)
           call wrt3(temp,nx,iojkao(kham),num8)
           call mult2(trans,dens,temp,ltran,ltran,num)
           call wrt3(dens,nx,iojk(kham),num8)
        else
           if (nprint.eq.5) then
            write(iwr,*)' 2J-K closed shell'
            call prtri(coul,num)
           endif
           call symh(coul,temp,iky,1,npair)
           call adonee(temp,dens,prefac)
           call wrt3(temp,nx,iojkao(kham),num8)
           call mult2(trans,dens,temp,ltran,ltran,num)
           call wrt3(dens,nx,iojk(kham),num8)
           ioff=ioff + nx
        endif
      endif
      if (nscmf .ne. 0) then
c
c     ----- open shell part -----
c
         do 280 ifo = 1,nscmf
c
         if (nprint.eq.5) then
          write(iwr,*)' coulomb matrix, ifo = ',ifo
          call prtri(coul(ioff),num)
         endif
         call symh(coul( ioff),temp,iky,1,npair)
         call wrt3(temp,nx,iojkao(kham+1),num8)
         call mult2(trans,dens,temp,ltran,ltran,num)
         call wrt3(dens,nx,iojk(kham+1),num8)
c
         if (nprint.eq.5) then
          write(iwr,*)' exchange matrix, ifo = ',ifo
          call prtri(exch(ioff),num)
         endif
         call symh(exch(ioff),temp,iky,1,npair)
         call wrt3(temp,nx,iojkao(kham+2),num8)
         call mult2(trans,dens,temp,ltran,ltran,num)
         call wrt3(dens,nx,iojk(kham+2),num8)
c
         ioff = ioff + nx
  280    kham = kham+2
c
c     ----- get the original trans back -----
c
      endif
      call end_time_period(TP_TEST6)
_IF(parallel)
      endif
_ENDIF
      if (nleft .lt. 1) go to 400
      iscmf = iscmf + nscmf
      go to 100
c
 400  return
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))
      end
      subroutine dhstarg(core,coul,exch,dens,prefac,rdmat,
     +           nclf,nscmf,irest,nprint)
c
c     dhstarg- forms multiple fock matrices 
c     doing only one integral pass 
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/maxlen)
INCLUDE(common/statis)
INCLUDE(common/timeperiods)
INCLUDE(common/runlab)

      dimension coul(*),exch(*),dens(*)
      dimension core(*)
c
      data pt5,two/0.5d0,2.0d0/
c
      nsheld=nclf+nscmf
      inull=igmem_null()
c
c     ----- divide the diagonal density elements by 2
c
      ioff=0
      do 300 i = 1,nsheld
      if (nprint.eq.5) then
       write(iwr,*)' density matrix, shell = ', i
       call prtri(dens(ioff+1),num)
      endif
      do 310 j = 1,num
      nij = ikyp(j) + ioff
 310  dens(nij)=dens(nij)*pt5
 300  ioff=ioff+nx
c
       if (irest.ne.1) then
        call vclr(coul,1,nx*nsheld)
        call vclr(exch,1,nx*nsheld)
       endif
c
c      now compute integrals 
c
      call start_time_period(TP_DHSTAR)
      call timana(5)
      call jandk(zscftp,core,coul,core(inull),exch,dens,core(inull),
     +           prefac,rdmat)
      call end_time_period(TP_DHSTAR)
      call cpuwal(begin,ebegin)
c
      if(nclf.gt.0) then
       if (nprint.eq.5) then
        write(iwr,*)' closed shell J matrix'
        call prtri(coul,num)
        write(iwr,*)' closed shell K matrix'
        call prtri(exch,num)
       endif
c
c     (2J-K) for closed shell
c
       call dscal(nx,two,coul,1)
       call vsub(coul,1,exch,1,coul,1,nx)
c
      endif
_IF(parallel)
c***   ***node-MPP***
      call start_time_period(TP_DHSTAR_GOP)
c...  gop results together
      ioff = 1
      if (nclf.ne.0) then
         call pg_dgop(300,coul(ioff),nx,'+')
         ioff = ioff + nx
      end if
      do 741 ifo=1,nscmf
         call pg_dgop(400+ioff,coul(ioff),nx,'+')
         call pg_dgop(500+ioff,exch(ioff),nx,'+')
741   ioff = ioff + nx
      call end_time_period(TP_DHSTAR_GOP)
c***   ***node-MPP***
_ENDIF
      return
      end
      subroutine wrtnsh(q,l2,nham,iblock,idaf)
      implicit REAL  (a-h,o-z)
      dimension q(*)
      call wrt3(q,l2,iblock,idaf)
      ii = l2+1
      ntri = nham-1
      do 10 itri=1,ntri
      call wrt3s(q(ii),l2,idaf)
10    ii=ii + l2
      return
      end
      subroutine wrtnshs(q,l2,nham,idaf)
      implicit REAL  (a-h,o-z)
      dimension q(*)
      call wrt3s(q,l2,idaf)
      ii = l2+1
      ntri = nham-1
      do 10 itri=1,ntri
      call wrt3s(q(ii),l2,idaf)
10    ii=ii + l2
      return
      end
      subroutine rdnsh(q,l2,nham,iblock,idaf)
      implicit REAL  (a-h,o-z)
      dimension q(*)
      call rdedx(q,l2,iblock,idaf)
      ii = l2+1
      ntri = nham-1
      do 10 itri=1,ntri
      call reads(q(ii),l2,idaf)
10    ii=ii + l2
      return
      end
      subroutine rdnshs(q,l2,nham,idaf)
      implicit REAL  (a-h,o-z)
      dimension q(*)
      call reads(q,l2,idaf)
      ii = l2+1
      ntri = nham-1
      do 10 itri=1,ntri
      call reads(q(ii),l2,idaf)
10    ii=ii + l2
      return
      end
_IF(debug_S)
      subroutine compare_S(sa,sb,l2,num)
      implicit none
      integer l2, num
      REAL sa(l2), sb(l2)
c
c...  Commons
c
INCLUDE(common/scra7)
INCLUDE(common/iofile)
c
c...  Local variables
c
      REAL    eomax,vomax,edmax,vdmax
      integer iomax,jomax,idmax,jdmax
      integer n,i,j
c
      logical opg_root
      external opg_root
c
      call pg_dgop(3362,sa,l2,'+')
      call rdedx(sb,l2,ibl7s,num8)
      if (opg_root()) then
         n = 0
         iomax = 0
         jomax = 0
         eomax = 0.0d0
         vomax = 0.0d0
         idmax = 0
         jdmax = 0
         edmax = 0.0d0
         vdmax = 0.0d0
         do i = 1, num
            do j = 1, i-1
               n = n+1
               if (dabs(eomax).lt.dabs(sa(n)-sb(n))) then
                  iomax = i
                  jomax = j
                  eomax = sa(n)-sb(n)
                  vomax = sb(n)
               endif
            enddo
            n = n+1
            if (dabs(edmax).lt.dabs(sa(n)-sb(n))) then
               idmax = i
               jdmax = j
               edmax = sa(n)-sb(n)
               vdmax = sb(n)
            endif
         enddo
         write(iwr,'("Max error in S-matrix off-diag:",2i5,2f20.15)')
     +        iomax,jomax,eomax,vomax
         write(iwr,'("Max error in S-matrix diag    :",2i5,2f20.15)')
     +        idmax,jdmax,edmax,vdmax
      endif 
      return
      end
_ENDIF
_IF(ga)
      subroutine denscf_ga(q,zscf)
c
c     performs one cycle scf with the density matrix from
c     routine denat
c     limited at the moment to direct scf closed shells
c     modified to use 4, rather than 6 triangles with GAs (jun 2003)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*),nsymm(8),msymm(8)
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
INCLUDE(common/machin)
INCLUDE(common/tran)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
INCLUDE(common/gjs)
INCLUDE(common/atmol3)
INCLUDE(common/cslosc)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/sta(210),cca(20),ra(20),scalea(20),iposa(20),
     + nstora,mpa,odiisa,
     * nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     * mpb,odiisb,nspacb,nsss(30),igvbo(maxorb),igvb1(maxorb),igsp(4),
     * ekk(63),intci(150)
      common/scra/iso(mxshel,48)
INCLUDE(common/mapper)
INCLUDE(common/datgue)
INCLUDE(common/nshel)
INCLUDE(common/prints)
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/runlab)
INCLUDE(common/scfopt)
INCLUDE(common/scra7)
INCLUDE(common/segm)
INCLUDE(common/symtry)
INCLUDE(common/scfwfn)
INCLUDE(common/timeperiods)
INCLUDE(common/harmon)
INCLUDE(common/zorac)
INCLUDE(common/gadiis)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/statis)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     + iacsct(508)

      Integer n_irreps
      Logical use_symmetry, force_serial

      character*5 fnm
      character*9 snm
      data fnm,snm/'scf.m','denscf_ga'/

      character*1 xn
      data xn/'n'/

      data m51,small/51,1.0d-3/
      data m3,m167 /3,167/
c
      call check_feature('denscf_ga')
c
c..
c..    core partitioning (like hcore + extra-s)
c..    core for hstar seprate (i10/iexch/idmat)
c
      Call start_time_period( TP_DENSCF )
      iter = 0
      nav = lenwrd()
      out = nprint.eq.2
      l1 = num
      l2 = num*(num+1)/2
      lentri = l2
      l3 = num*num
      l0 = newbas0
c
      oswed3(4) = .true.
      oswed3(8) = .true.

      if (oatdo.and.isecat.eq.-1) then
c...   for atorbs most is not needed
c...   i10 orbs i30 backtranformed orbs / i80 occupations / i70 eps
       i10 = igmem_alloc_inf(2*l3+2*l1,fnm,snm,'i10',IGMEM_DEBUG)
       irdmat = i10
       i30 = i10 + l3
       i70 = i30 + l3
       i80 = i70 + l1
       ndum = nprint
       nprint = -5
       call getq(q(i10),q(i70),q(i80),l1,l1,3,ieig,ipop,mouta,'natorb') 
       nprint = ndum
c...    set noc1
       noc1 = 0
       do i=1,l1
         if (q(i80+i-1).gt.0.0d0) noc1 = noc1 + 1
       end do
       go to 10
      end if
c
c...  Begin constructing the core partitioning:
c...  Compute the amount of memory needed to store the local 
c...  datastructures.
c
      nss = 1
      ilen = ikyp(nshell)*2 +  4 * l2 + 2 * l1
      irdmat = igmem_alloc_inf(ilen,fnm,snm,'irdmat',IGMEM_NORMAL)
      iprefa=irdmat+ikyp(nshell)
      i10 = iprefa+ikyp(nshell)
      ifock = i10
      i20 = i10 + l2
      i30 = i20 + l2
      idmat = i30
      i40 = i30 + l2
****  6 to 4 trangles triggered here
      i50 = i30 
      i60 = i50 + l2
****
      i70 = i60 + l2
      i80 = i70 + l1
      last = i80 + l1
c
c...  Now <last> is the amount of memory needed.
      length = last-irdmat
c
c consistency check
      if (length .ne. ilen)then
         write(iwr,*)length,last,irdmat,ilen
         call caserr('size error')
      endif

      ierror = CD_update_geom(c)
      if (CD_active().and.odenscfdft) then
         call retrieve_spare(imemspare)
         imemfree = igmem_max_memory() - imemspare
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         orks = CD_is_rks()
         ierror = CD_rks()
c        if (ks_bas.eq.KS_AO) then
            imemreq = CD_memreq_energy_ao(q,q,iwr)
c        else if (ks_bas.eq.KS_MO) then
c           imemreq = CD_memreq_energy_mo(l1,na,0,q,q,iwr)
c        else if (ks_bas.eq.KS_AOMO) then
c           imemreq = CD_memreq_energy(l1,na,0,q,q,iwr)
c        endif
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            call caserr('Out of memory in incore coulomb fit')
         endif
         if (.not.orks) ierror=CD_uks()
      endif
c
      omaster = .true.

         call start_time_period(TP_RDMAT)
         m171t=171
         len171=  lensec(l2)
         nshtri=ikyp(nshell)
         nshblk=lensec(nshtri)
         lent = len171
         if(odelta) then
          ntri = 4
         else
          ntri = 0
         endif
         iof171(1)=0
         iof171(2)= iof171(1) + nshblk
         iof171(3)= iof171(2) + nshblk
         do i=4,6
          iof171(i)=iof171(i-1)+lent
         enddo
         if(irest.eq.1.or.irest.eq.2) then
            call secget(isect(471),m171t,ibl171)
            write(iwr,2990)irest,ibl171
 2990       format(/1x,'restarting denscf in integrals :',
     *      'restart parameter =',i3/
     *      1x,'section 171 at block ',i5)
            call rdedx(q(iprefa),nshtri,ibl171,idaf)
            call reads(q(irdmat),nshtri,idaf)
            call rdedx(q(ifock),l2,ibl171+iof171(3),idaf)
            call rdedx(q(idmat),l2,ibl171+iof171(6),idaf)
            dlnmxd=-9999999.0d0
            do 1672 kkk=1,nshtri
               if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
               dlnmxd=q(irdmat+kkk-1)
1672        continue
         endif
         ln171=ntri*len171+nshblk*2
         call secput(isect(471),m171t,ln171,ibl171)
         if(irest.eq.1) go to 6353
         call rdedx(q(idmat),l2,ibl3pa,idaf)
         call rdmake(q(iprefa))
         if(out) then
            write(iwr,2980) ibl171,ln171
 2980       format(/1x,
     &      'output section 171 to block ',i5/1x,
     &      'section length              ',i5/)
         endif
         call wrt3(q(iprefa),nshtri,ibl171,idaf)
         call mkrdmt('guess',q(irdmat),q(idmat),l2,nprint)
         call wrt3s(q(irdmat),nshtri,idaf)
         if(odelta) then
          call zer171(q(i10),l2,ntri,ibl171+iof171(3),len171,idaf)
          call wrt3(q(idmat),l2,ibl171+iof171(6),idaf)
         endif
6353     continue
c
6351  continue
      if (odscf) then
c...      perform closed-shell direct scf
         if( out ) write(iwr,3030)
3030     format(/5x,25('*'))
         dlntol = tolitr(1) - dmin1(dmax1(dlnmxd,delfac),0.0d0)
         if ( out ) write (iwr,6010) dlntol
c ...
c         we must reset "zscftp" in /runlab/ to ensure
c         the closed shell fock builder is called from intega etc
c
          zsave = zscftp
          zscftp = 'rhf'
          call end_time_period(TP_RDMAT)
      endif
c
c    ----- construct a skeleton fock matrix -----
c    h       at q(i10)
c    density at q(idmat)
c    exch    at q(iexch) (scratch)
c
c Modify fock builder options
c
      if (odenscfdft) then
         idum = CD_set_2e()
      endif
          call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),irest)
          zscftp = zsave
          if(irest.ne.0) go to 6352
c
c
c restore fock builder options
c
      if (odenscfdft) then
         idum = CD_reset_2e()
      endif
c
      if (out) then
         write (iwr,6060)
         call prtril(q(i10),l1)
      end if
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at q(i10)
c     scratch area at q(i20)
c
      call symh(q(i10),q(i20),iky,0,0)
c
      if (out) then
         write (iwr,6070)
         call prtril(q(i10),l1)
      end if
c
c     ---- read in core hamiltonian matrix, one electron ints in q(i20)
c     - h - at q(i10)
c
      call rdedx(q(i20),l2,ibl7f,num8)
c
      if (ozora) call zora(q,q(idmat),q(i20),'read')
c
      call vadd(q(i10),1,q(i20),1,q(i10),1,nx)
c
      if(CD_active().and.odenscfdft)then
c
c...     Store current spintyp setting
c
         orks = CD_is_rks()
c
c...     Make sure we are running closed shell DFT here
c
         ierror = CD_rks()
         oignore = CD_ignore_accuracy()
         ierror = CD_set_ignore_accuracy(.true.)
         idum = CD_energy_ao(c,q(i10),dum,q(idmat),dum,
     +        edft,q,q,outon,dft_accu,iwr
     +        )
         ierror = CD_set_ignore_accuracy(oignore)
         call symm_op(q(i10))
c
c...     Restore original spintyp setting
c
         if (.not.orks) ierror=CD_uks()
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr(
     +         'Memory failure in rhfcld:CD_jfit_clean2')
         endif
      endif
c
      if (out) then
         write (iwr,6050)
         call prtril(q(i10),l1)
      end if
***
***   GA
       call load_ga_from_triangle(ih_scr2,q(i10),l1)
***
c
c      make orthormalizing vectors at q(i30)
c..
c..    first get canonical orthonormal vectors like hcore
c..    read s-matrix in i10, vectors end up at i30 and on idaf (ibl3qs)
c..    this is symmetry adapted to keep scf happy
c
c     ----- read in fock transformation matrix -----
c           transform hamiltonian matrix
c
c       at q(i10) overlapmatrix
c       at q(i30) orthonormalised vectors, output
c       at q(i70) and q(i20) scratch
c
      call start_time_period(TP_ORFOG)
      call rdedx(q(i10),l2,ibl7st,num8)

c serialise the diag here until modified further
c         itemp = idpdiag
c         idpdiag=99999999
c
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
c
      n_irreps = 0
      Do i = 1, l1
         isymmos( i ) = isymaos( i )
         n_irreps = Max( n_irreps, isymaos( i ) )
      End Do
      use_symmetry = n_irreps .GT. 1 .and. symm_diag
*
      call qmat_symm(q,q(i10),q(i30),q(i70),q(i20),iky,l0,l1,l3,l1,out,
     +               isymmos, use_symmetry )
c
      call wrt3(q(i30),l3,ibl7la,num8)
c
c     - q - at q(i30) orthonormalizing transformation vectors
c     - h'- at q(i30) transformed h matrix
c     (both h and h' require space of square)
c
      if (symm_diag) then
         call characterize_mo( l1, l0, q( i30 ), isymaos, 
     +                         n_irreps, isymmos, ierr )
         If( ierr .NE. 0 ) Then
            If( opg_root() ) Then
               Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                       'symmetry of mo''s for use in the diag.'
               Write( 6, * ) 'Ignoring symmetry in diag'
            end if
         end if
      else
        ierr = 999
      endif
      use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. symm_diag
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
**** GA
         call load_ga_from_square(ih_vec,q(i30),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i30),ih_scr2, l1)
      call end_time_period(TP_ORFOG)
c
c     ----- diagonalize new hamiltonian matrix -----
c
c      at q(i30) transformed fock matrix
c      at q(i10) eigenvectors
c      at q(i70) eigenvalues
c
      m2 = 2
c
c     WARNING: current setting of 1.0d-8 could cause potential problems 
c     in some high symmetry cases .. consider resetting diaacc to as low
c     as 1.0d-11 (the lower bound at SCF convergence)
c
      diaacc = 1.0d-8
c
c  broadcast fock to all nodes 
c
      call start_time_period(TP_DIAG)
      call jacobi_symm(q(i30),iky,l0,q(i10),ilifq,l0,q(i70),
     +     m2,2,diaacc, isymmos, use_symmetry)
      call end_time_period(TP_DIAG)
c
c     ----- back-transform the eigenvectors -----
c
c     eigenvectors at q(i10)
c     the famous orthonormalising vectors from qmat at q(i30)
c
      call rdedx(q(i30),l3,ibl7la,num8)
*** GA
 
         call load_ga_from_square(ih_scr,q(i10),l1)
         call load_ga_from_square(ih_vec,q(i30),l1)
 
         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)
 
         call load_square_from_ga(q(i10),ih_scr2, l1)
 
c
c...  re-entry oatdo
c     need i10 (orbitals) i70 (e's) i80 (occupations)
c
10    continue
c
c...  select mo-s
c
      noc = na
      nsav = mouta
      mswap = nswapa
c
c...   swap orbitals if requested
      if (mswap.gt.0) call swap(q(i10),q(i70),num,mswap,nswap)
      if (oatdo.and.isecat.eq.-1) then
c..       swap occupations along and do not replace them for atorbs
         do i=1,mswap
            temp = q(i80-i+nswap(i))
            q(i80-i+nswap(i)) = q(i80-i+nswap(i+1))
            q(i80-i+nswap(i+1)) = temp
         end do
      else  
c
c..     set up occupations in i80 like hcore
c
         call llvmo(q(i70),q(i80),noc,noc1,l0)
         pop = 2.0d0
         call dscal(l0,pop,q(i80),1)
      end if
c
c...   save ** mo-s and eigenvalues / orbital occupations **
c...          (i10)      (i70)             (i80)
c
      call putq(zcom,ztitle,q(i70),q(i80),l1,l1,l0,1,1,q(i10),nsav,
     +          ibl3qa)
      call copq_again(nsav)
c
      if (oprint(45)) write(iwr,6080)
      if(otran)go to 20
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
      do 30 i=1,l0
      ii=i10+ilifq(i)
      ibig=idamax(l1,q(ii),1)
      isymmos(i)=isymaos(ibig)
        bigg = 0.0d0
        ibig = 0
        do 50 j=1,l1
         if (isymmos(i).ne.isymaos(j).and.dabs(q(ii+j-1)).gt.bigg) then
            bigg = dabs(q(ii+j-1))
            ibig = j
         end if
50      continue
        if (bigg.gt.small) write(iwr,6090) i,ibig,bigg
30    continue
c     final check    aos versus mos
      do 70 i=1,8
      msymm(i)=0
 70   nsymm(i)=0
c     build up ao array
      do 80 i=1,l1
      j=isymaos(i)
 80   nsymm(j)=nsymm(j)+1
c     build up mo array
      do 90 i=1,l1
      j=isymmos(i)
 90   msymm(j)=msymm(j)+1
c      now cross check
      do 170 i=1,8
      if(msymm(i).eq.nsymm(i))goto170
      if ((oharm.or.odepen).and.msymm(i).eq.nsym0(i)) go to 170
      write(iwr,6100)i,nsymm(i),msymm(i)
 170   continue
 20   if (oprint(45)) then
       do 100 i=1,l0
       if(otran)write(iwr,6110)i,q(i70+i-1),q(i80+i-1)
       if(.not.otran)write(iwr,6120)i,isymmos(i),q(i70+i-1),q(i80+i-1)
 100   continue
       write(iwr,6130)
      endif
c
c     revised section 190 with mo symmetries
c
      if(.not.otran)call wrt3i(mmmm,mach(13)*nav,iblk51,idaf)
      call setsto(1360,0,mmmm)
c
c = end of analysis =
c
c..     transform back to original basis and print if requested
c
      call tdown(q(i30),ilifq,q(i10),ilifq,l0)
c
      if (oprint(45)) then
         write (iwr,6020)
         call prev(q(i30),q(i70),l0,l1,l1)
      end if
c
c...  now do density matrices
c...  i10 resulting d-matrix
c...  i30 orbitals
c...  i80 occupations
c...  i70 eigenvalues
c
      call dmtxp(q(i30),q(i10),q(i80),iky,noc1,l1,l1)
      call wrt3(q(i10),l2,ibl3pa,idaf)
      call wrt3(q(i70),l1,ibl3ea,idaf)
c all nodes here
      irest = 0
c
c...   reset core
c
6352  continue
      call gmem_free_inf(irdmat,fnm,snm,'irdmat')
c
      Call end_time_period( TP_DENSCF )
      return
 6010 format (5x,'dlntol =  ',f15.10/5x,25('*'))
 6020 format (//30x,28('-'),/,30x,'   initial guess orbitals   ',/,30x,
     +        28('-')//)
 6050 format (20x,11('-')/20x,'fock matrix'/,20x,11('-'))
 6060 format (/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 6070 format (/20x,23('-')/20x,'symmetrized fock matrix'/20x,
     +        23('-'))
 6080 format(/1x,47('=')/
     *' m.o.  irrep  orbital energy   orbital occupancy'/
     *1x,47('=')/)
 6090 format(' ** symmetry contamination of mo',i4,' at sabf',i4,
     *         ' of ',e12.5,'  **')
 6100 format(/
     *' *** warning -error in m.o. symmetry designation'//
     *4x,'irrep. ',i2,' no. of a.o.s =',i3/
     *14x,'no. of m.o.s =',i3/)
 6110 format(1x,i3,7x,f16.8,f20.7)
 6120 format(1x,i3,i7,f16.8,f20.7)
 6130 format(/1x,47('=')/)
      end
_ENDIF
      subroutine ver_scf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/scf.m,v $
     +     "/
      data revision /"$Revision: 6059 $"/
      data date /"$Date: 2009-09-09 12:17:28 +0200 (Wed, 09 Sep 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
