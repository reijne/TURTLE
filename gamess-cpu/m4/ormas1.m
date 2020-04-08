c
c  ormas stands for occupation restricted multiple active space.
c  i originally called this full class ci, which explains why
c  some things are labelled fc---.
c
c*module ormas1  *deck fcinput
c     -----------------------------------------
      subroutine fcinput(nprint,gpname,gpname2,x,molab)
c     -----------------------------------------
c
      implicit double precision(a-h,o-z)
c
      logical some,pures,goparr,dskwrk,maswrk,abel
c
      parameter (mxatm=2000, mxrt=100, mxao=8192, mxsh=5000)
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
INCLUDE(common/fccwfn)
      integer ncore_f,nact_f,nels_f
      common /fccwfn2 / ncore_f,nact_f,nels_f
c
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      common /ijpair/ ia(mxao)
      common /infoa_f / nat,ich,mul,num,nqmt,ne,ma,mb,
     *                zan(mxatm),c(3,mxatm),ian(mxatm)
      common /iofile_f/ ir,iw,ip,is,ijkt,idaf,nav,ioda(950)
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
      common /symmol/ group,complex,igroup,naxis,ilabmo,abel
      common /symtry_f/ mapshl(mxsh,48),mapctr(mxatm,48),
     *                t(432),invt(48),nt
      common /wfnopt_f/ scftyp,cityp,dftype,cctyp,mplevl,mpctyp
c
c     /fmcom/ unused here as x is passed in under the G-UK scheme
c     common /fmcom / x(1)
c     gamess-uk memory management, replaces /fmcom/
c
      dimension x(*), molab(*)
c
      parameter (zero=0.0d+00, one=1.0d+00, two=2.0d+00)
c
      dimension fant(8),lfant(8),gant(27),lgant(8)
c
      data det,cidet/8hdet     ,8hcidet   /
      data hess/8hhessian /
      data rnone/8hnone    /
      data fant/8hc1      ,8hci      ,8hcs      ,8hc2      ,
     *          8hd2      ,8hc2v     ,8hc2h     ,8hd2h     /
      data lfant/1,1,1,1,2,2,2,3/
      data gant/8ha       ,8hag      ,8hau      ,8ha'      ,
     *          8ha"      ,8ha       ,8hb       ,8ha       ,
     *          8hb1      ,8hb2      ,8hb3      ,8ha1      ,
     *          8ha2      ,8hb1      ,8hb2      ,8hag      ,
     *          8hbg      ,8hbu      ,8hau      ,8hag      ,
     *          8hb1g     ,8hb2g     ,8hb3g     ,8hau      ,
     *          8hb1u     ,8hb2u     ,8hb3u     /
      data lgant/0,1,3,5,7,11,15,19/
c
      some = maswrk  .and.  nprint.ne.-5  .and.  nprint.ne.-23
c
      if(some) write(iw,9000)
c
c     set up input to specify the ci space
c     select point group and electron/orbital counts
c
      grpdet = fant(1)
      if(igroup.eq.1)                grpdet = fant(1)
      if(igroup.eq.3)                grpdet = fant(2)
      if(igroup.eq.2)                grpdet = fant(3)
      if(igroup.eq.4.and.naxis.eq.2) grpdet = fant(4)
      if(igroup.eq.8.and.naxis.eq.2) grpdet = fant(5)
      if(igroup.eq.7.and.naxis.eq.2) grpdet = fant(6)
      if(igroup.eq.6.and.naxis.eq.2) grpdet = fant(7)
      if(igroup.eq.9.and.naxis.eq.2) grpdet = fant(8)
      if (nt.eq.1) grpdet=fant(1)
c
      ncore  = ncore_f
      nact   = nact_f
      nels   = nels_f
c
c     ensure orbital labels are unit if symmetry has been 
c     switched off, this may occur for srso calc. but we 
c     wish to preserve symmetry data for other parts of calc.
c
      if (group.eq.fant(1)) then
        do i = 1, num
          molab(i) = 1
        end do
      end if
c
c     by convention, set the z-component of spin to the total spin
c
      s = 0.5d0*dble(mul-1)
      sz = s
c
c     set up input to control the diagonalization
c
      prttol = 0.05d+00
c
c     set up input to control the first order density computation
c
      do 5 i=1,mxrt
         nflgdm(i) = 0
    5 continue
      nflgdm(1)=1
c
c     set up input to control the second order density computation
cgdf  in all likelihood we will not want to change pures but keep
c     here for future option, same goes for fdirct,qcorr
c
      pures = .true.
c
      igpdet = -1
c
c     set up defaults for fdirct, qcorr.
c
      if (gpname.eq.cidet) fdirct=.true.
      if (gpname.eq.det)   fdirct=.false.
c
      qcorr =.true.
      if (gpname.eq.det) qcorr=.false.
c
c     the input for c2h is supposed to be identical to the guga
c     order, namely 1,2,3,4=ag,bu,bg,au, but the ci code wants
c     the order of  1,2,3,4=ag,bg,bu,au.  see also gajasw routine.
c
      if (grpdet.eq.fant(7)) then
         modi = kstsym
         if(kstsym.eq.2) modi=3
         if(kstsym.eq.3) modi=2
         kstsym=modi
      end if
      do i=1,8
         if (grpdet.eq.fant(i)) then
            igpdet=lfant(i)
            stsym = gant(lgant(i)+kstsym)
         endif
      end do
      if (igpdet.eq.-1) then
         if(maswrk) write(iw,*) '$det point group is unrecognized!'
         call abrt
         stop
      endif
      if (grpdet.eq.fant(1).and.kstsym.gt.1) then
         if(maswrk) write(iw,*) '$det state symm is not correct irrep'
         call abrt
         stop
      endif
      if (kstsym.gt.(2**igpdet)) then
         if(maswrk) write(iw,*)
     *       '$det state symmetry is too large for this group'
         call abrt
         stop
      endif
c
      l0 = nqmt
      l1 = num
      l2 = (l1*l1+l1)/2
      l3 =  l1*l1
c
      call valfm(loadfm)
      lmolab = loadfm + 1
      lmoirp = lmolab + l1
      lvec   = lmoirp + l1
      ls     = lvec   + l3
      lq     = ls     + l2
      lwrk   = lq     + l3
      lmodeg = lwrk   + l1
      last   = lmodeg + l1
      needd   = last - loadfm - 1
      call getfm(needd)
c
c     this is the reason it is simplest to have 2 versions 
c     of the symmetry labels (see masscf.m)
c
      call gajasw(molab,num,grpdet)
c
c     1.  set ncorsv,ncor,nact,norb,na,nb for determinant specification
c         and copy into internally used variable names.
c         ncor will be set to zero to drop cores, so ncorsv saves this.
c
      ncorsv = ncore
      ncor   = ncore
      norb   = ncore + nact
      nhigh = int(sz+sz+0.0001d+00)
      nb = (nels-nhigh)/2
      na = nb+nhigh
      ma = na+ncorsv
      mb = nb+ncorsv
      neltot = 2*ncor+na+nb
      nerr=0
      if(neltot.ne.ne)  nerr=1
      if(nels.ne.na+nb) nerr=1
      if(na.lt.nb)      nerr=1
      if(na.le.0)       nerr=1
      if(nb.lt.0)       nerr=1
      if(nprint.ne.-23  .and.  nerr.gt.0) then
         if(maswrk) write(iw,9030) ncore,nels,sz,ich,mul
         call abrt
         stop
      end if
c
      if(nprint.ne.-23  .and.  nstate.gt.mxrt) then
         if(maswrk) write(iw,9035) nstate,mxrt
         call abrt
         stop
      end if
c
c     3. setup for 1st order density computation
c
      if(iroot.gt.nstate) then
         if(maswrk) write(iw,9036) iroot,nstate
         call abrt
         stop
      end if
      if(iroot.gt.mxrt) then
         if(maswrk) write(iw,9037) iroot,mxrt
         call abrt
         stop
      end if
      if(nflgdm(iroot).eq.0) nflgdm(iroot)=1
c
c     4. setup for state-averaging 1st and 2nd order densities.
c     -iwts- indexes the non-zero elements of -wstate-
c
      ipures=0
      if(pures) ipures=1
      mxstat=0
      wsum = zero
      do 15 i=1,mxrt
         if(wstate(i).gt.zero) then
            if(i.le.nstate) then
               wsum = wsum + wstate(i)
               mxstat = mxstat+1
               iwts(mxstat) = i
            else
               if(maswrk) write(iw,9040) nstate
               call abrt
               stop
            end if
         end if
         if(wstate(i).lt.zero) then
            if(maswrk) write(iw,9050)
            call abrt
            stop
         end if
   15 continue
      scale = one/wsum
      call dscal(mxrt,scale,wstate,1)
c
c     5. check if determination of davidson correction is meaningful
c
      if (qcorr) then
         if (nspace.gt.3.or.nspace.lt.2) go to 2000
         if (maxi(nspace).ne.2) go to 2000
         if (mini(nspace).ne.0) go to 2000
         if (maxi(1)-mini(1).ne.2) go to 2000
         go to 3000
 2000    qcorr=.false.
         write(iw,9055)
 3000    continue
      endif
c
      if(some) then
         write(iw,9100) grpdet,stsym,ncor,nact,
     *                  na+ncor,na,nb+ncor,nb,norb
c
         write(iw,9110) nstate,kst,maxp,maxw1,niter,crit
         if(gpname.eq.cidet) then
            write(iw,9120) iroot
            write(iw,9130) (nflgdm(ii),ii=1,nstate)
         end if
         if(gpname.eq.det) then
            write(iw,9140) pures
            write(iw,9150) (iwts(ii),wstate(iwts(ii)),ii=1,mxstat)
            if(mxstat.gt.1) write(iw,9159) iroot
         end if
      end if
c
c     if (some) then
c        write(iw,9157) fdirct
c        write(iw,9158) qcorr
c     endif
c
cgdf  mosypr needs the symmetry label -symbols- (e.g. a1,ag, etc)
c
c     if(some) then
c        write(iw,9160) ncor,nact
c        call mosypr(x(lmolab),ncor,nact)
c     endif
c
c    determine and write out the ormas information
c
      if(some) then
         write(iw,9170) nspace,(msta(i),i=1,nspace)
         write(iw,9180) (mnum(i),i=1,nspace)
         write(iw,9190) (mini(i),i=1,nspace)
         write(iw,9195) (maxi(i),i=1,nspace)
      endif
c
c     send ormas information for checking and determination of
c     other important information.
c
      call fccheck(iw,some,
     *             nspace,mnum,mini,maxi,iami,iama,ibmi,ibma,
     *             na,nb)
c
      call retfm(needd)
c
      return
c
 9000 format(/5x,60(1h-)/
     *   5x,'  direct determinant ormas-ci input sorter'/
     *   5x,'  program written by joe ivanic and mike schmidt'/
     *   5x,'  ormas = occupation restricted multiple active space'/
     *   5x,60(1h-))
 9030 format(/1x,'**** error, input ncore=',i4,' nels=',i4,
     *          ' sz=',f6.3/
     *       1x,' is inconsistent with input ich=',i4,
     *          ' mult=',i4)
 9035 format(/1x,'***** error, requested number of ci roots=',i5/
     *        1x,'exceeds the dimension limit for number of states',i5)
 9036 format(/1x,'**** error, your state selected for properties=',i5/
     *        1x,'exceeds the number of roots you requested=',i5)
 9037 format(/1x,'**** error, your state selected for properties=',i5/
     *        1x,'exceeds the dimension limit for number of states',i5)
 9040 format(/1x,'**** error, weights assigned to states higher',
     *          ' than nstate=',i5)
 9050 format(/1x,'**** error, negative value for -wstate- ???')
 9055 format(/1x,'**** davidson + q correction meaningless for ',
     *           'current ormas paritioning.')
c
 9100 format(/1x,'the point group                  =',3x,a8/
     *       1x,'the state symmetry               =',3x,a8/
     *       1x,'number of core orbitals          =',i5/
     *       1x,'number of active orbitals        =',i5/
     *       1x,'number of alpha electrons        =',i5,
     *          ' (',i4,' active)'/
     *       1x,'number of beta electrons         =',i5,
     *          ' (',i4,' active)'/
     *       1x,'number of occupied orbitals      =',i5)
 9110 format(/1x,'number of ci states requested    =',i5/
     *       1x,'number of ci starting vectors    =',i5/
     *       1x,'max. no. of ci expansion vectors =',i5/
     *       1x,'size of initial ci guess matrix  =',i5/
     *       /1x,'max. no. of ci iters/state       =',i5/
     *       1x,'ci diagonalization criterion     =',1p,e9.2)
c
 9120 format(1x,'ci properties will be found for root number',i4)
 9130 format(1x,'1e- density matrix options are',20i2)
 9140 format(1x,'pure spin state averaged 1e- and 2e- density matrix',
     *          ' option=.',l1,'.')
 9150 format(2(1x,'state=',i4,' dm2 weight=',f10.5,4x,:))
c9157 format(/1x,'fully direct option              =',l5)
c9158 format(/1x,'calc. of davidson q correction   =',l5)
 9159 format(1x,'state-averaged mcscf using iroot=',i3,
     *          ' as its state-specific e')
c
c9160 format(/1x,'symmetries for the',i4,' core,',i4,' active are')
c
 9170 format(/1x,'the number of spaces             =',i4/
     *        1x,'each space starts at orbital     =',50i4)
 9180 format(1x,'no of orbitals in each space     =',50i4)
 9190 format(1x,'min no of elecs in each space    =',50i4)
 9195 format(1x,'max no of elecs in each space    =',50i4)
c
      end
c
c*module ormas1  *deck fccheck
c     -------------------------------------------------------------
      subroutine fccheck(iw,some,
     *             nspace,mnum,mini,maxi,iami,iama,ibmi,ibma,
     *             na,nb)
c     -------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical some
      dimension mnum(nspace),mini(nspace),maxi(nspace)
      dimension iami(nspace),iama(nspace)
      dimension ibmi(nspace),ibma(nspace)
c
      ntot = na + nb
c
c   preliminary checks:  abort if specifications do not pass.
c
c  1/ check to see number of orbitals in each space > 0
c
      do ii=1,nspace
         if (mnum(ii).le.0) then
            if (some) write(iw,9000) (mnum(i),i=1,nspace)
            call abrt
            stop
         endif
      enddo
c
c  2/ check to see that minimum number of electrons in each space i not
c     greater than 2*mnum(i)
c
      do ii=1,nspace
         if (mini(ii).gt.2*mnum(ii)) then
            if (some) write(iw,9010) ii,mini(ii),mnum(ii)
            call abrt
            stop
         endif
      enddo
c
c  3a/ check to see that sum of minimum number of electrons in each
c  space is less than or equal to ntot.
c
c  3b/ also check to see that sum of maximum number of electrons in
c  each space is greater than or equal to ntot.
c
      isum = 0
      isum2 = 0
      do ii=1,nspace
         isum = isum + mini(ii)
         isum2 = isum2 + maxi(ii)
      enddo
      if (isum.gt.ntot) then
         if (some) write(iw,9020) isum,ntot
         call abrt
         stop
      endif
      if (isum2.lt.ntot) then
         if (some) write(iw,9030) isum2,ntot
         call abrt
         stop
      endif
c
c  4a/ check to see that maximum number of electrons in each space i
c     is not greater than capacity of that space.
c
      do ii=1,nspace
         if (maxi(ii).gt.2*mnum(ii)) then
            if (some) write(iw,9040) ii,maxi(ii),mnum(ii)
            call abrt
            stop
         endif
      enddo
c
c  4b/ check to see that maximum number of electrons in each space i
c      is not greater than ntot
c
      do ii=1,nspace
         if (maxi(ii).gt.ntot) then
            if (some) write(iw,9050) ii,maxi(ii),ntot
            call abrt
            stop
         endif
      enddo
c
c  5/ check to see that all maxi(i).ge.mini(i).
c
      do ii=1,nspace
         if (maxi(ii).lt.mini(ii)) then
            if (some) write(iw,9060) ii,maxi(ii),mini(ii)
            call abrt
            stop
         endif
      enddo
c
c  preliminary checks passed.
c
      if (some) write(iw,9070)
c
c  now for reduncancy checks.
c
c  6a/ check to see that maximums are not impossible, if
c     so adjust accordingly.
c
      istat = 0
c
      nr = ntot
      do ii=1,nspace
         nr = nr - mini(ii)
      enddo
c
      do ii=1,nspace
         itemp = maxi(ii) - mini(ii)
         imax = min(itemp,nr) + mini(ii)
         if (imax.ne.maxi(ii)) then
            maxi(ii) = imax
            istat=1
         endif
      enddo
c
c  6b/ check to see that minimums are not redundant, if
c      so adjust accordingly.
c
      do ii=1,nspace
         it1 = ntot
         im1 = mini(ii)
         do 50 kk=1,nspace
            if (kk.eq.ii) go to 50
            it1 = it1 - maxi(kk)
   50    continue
         im2 = max(it1,im1)
         if (im2.ne.mini(ii)) then
            mini(ii) = im2
            istat=1
         endif
      enddo
c
      if (istat.eq.1) then
         if (some) write(iw,9080) (maxi(i),i=1,nspace)
         if (some) write(iw,9085) (mini(i),i=1,nspace)
      endif
c
c  now to determine the minimum and maximum alpha and beta parts.
c
      do ii=1,nspace
         it1 = min(mnum(ii),maxi(ii))
         iama(ii) = min(it1,na)
         ibma(ii) = min(it1,nb)
      enddo
c
      do ii=1,nspace
         ita1 = max(0,mini(ii)-ibma(ii))
         itb1 = max(0,mini(ii)-iama(ii))
         ita2 = na
         itb2 = nb
         do 100 kk=1,nspace
            if (kk.eq.ii) go to 100
            ita2 = ita2 - iama(kk)
            itb2 = itb2 - ibma(kk)
  100    continue
         iami(ii) = max(ita1,ita2)
         ibmi(ii) = max(itb1,itb2)
      enddo
c
      if (some) then
      write(iw,9090) (iami(i),i=1,nspace)
      write(iw,9100) (iama(i),i=1,nspace)
      write(iw,9110) (ibmi(i),i=1,nspace)
      write(iw,9120) (ibma(i),i=1,nspace)
      endif
c
      return
c
 9000 format(/1x,'no. orbitals in each space must be bigger than zero'/
     *        1x,'these were input as :',50i4)
 9010 format(/1x,'minimum no of electrons in space ',i4/
     *        1x,'is larger than capacity of this space'/
     *        1x,'input minimum        = ',i4/
     *        1x,'input no of orbitals = ',i4)
 9020 format(/1x,'total number of minimum electrons is larger than',
     *        1x,'number of active electrons'/
     *        1x,'total number of minimum electrons = ',i4/
     *        1x,'total number of active  electrons = ',i4)
 9030 format(/1x,'total number of maximum electrons is less than',
     *        1x,'number of active electrons'/
     *        1x,'total number of maximum electrons = ',i4/
     *        1x,'total number of active  electrons = ',i4)
 9040 format(/1x,'maximum no of electons in space ',i4/
     *        1x,'is larger than capacity of this space'/
     *        1x,'input maximum        = ',i4/
     *        1x,'input no of orbitals = ',i4)
 9050 format(/1x,'maximum no of electrons in space ',i4/
     *        1x,'is larger than number of active electrons'/
     *        1x,'input maximum                    = ',i4/
     *        1x,'total number of active electrons = ',i4)
 9060 format(/1x,'error in specification for space ',i4/
     *        1x,'maximum number of electrons is less than',
     *        1x,'minimum number.'/
     *        1x,'maxi = ',i4/
     *        1x,'mini = ',i4)
 9070 format(/1x,'space specifications have passed preliminary checks')
 9080 format(/1x,'reduncancies found, maximum/minimum numbers of',
     *        1x,'electrons have been adjusted.'//
     *        1x,'max no of elecs in each space    = ',50i4)
 9085 format(1x, 'min no of elecs in each space    = ',50i4)
 9090 format(/1x,'all valid determinants will be expressed',
     *        1x,'as pairs of alpha and beta strings.'//
     *        1x,'min no of alpha elecs = ',50i4)
 9100 format(1x,'max no of alpha elecs = ',50i4)
 9110 format(/1x,'min no of beta  elecs = ',50i4)
 9120 format(1x,'max no of beta  elecs = ',50i4)
c
      end
c
c*module ormas1  *deck defcci_init
      subroutine defcci_init(nprint,x,molab)
c
      implicit double precision(a-h,o-z)
c
      logical clabel
      logical some,goparr,dskwrk,maswrk,fdirct,qcorr
c
      parameter (mxrt=100, mxatm=2000)
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      common /enrgys/ enucr,eelct,etot,stot,ssquar,ecore,escf,eerd,
     *                e1,e2,ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
c     common /fmcom / x(1)
      common /infoa_f / nat,ich,mul,num,nqmt,ne,ma,mb,
     *                zan(mxatm),c(3,mxatm),ian(mxatm)
      common /iofile_f/ ir,iw,ip,is,ijkt,idaf,nav,ioda(950)
      common /machin_f/ nwdvar,maxfm,maxsm,limfm,limsm
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
      common /wfnopt_f/ scftyp,cityp,dftype,cctyp,mplevl,mpctyp
c
c     addresses needed by defcci calls in here
c
      common /defaddr/ econst,lsint2,lwk,m1,m2,m4,ksp,kel,nsym,lbox1
     *,lbox2,lbox3,lbox4,lbox5,kcoeff,kab,kq,kb,kef,kf,kec,kgr,iposa,
     * ipera,iind1,igroa,iwrk1,iwrk2,isd,iso,index,immc,ihmcon,lgmul,
     * lktab,lcon,lcon1,lcon2,lcon3,landet,lbndet,nast,nbst,lsyma,
     * lsymb,lgcom,lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,itga,itgb,
     * iast,ibst,nb1ex,ibst2,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy,
     * istat,need1,need2,need3,need4,ntot,ntco
c
c     gamess-uk memory management, replaces /fmcom/
c
      dimension x(*), molab(*)
c
      parameter (zero=0.0d+00, one=1.0d+00)
c
      data amcscf/8hmcscf   /
      data check/8hcheck   /
c
      istat=0
c
c     ----- driver for full class ci computation -----
c
      some = maswrk  .and.  nprint.ne.-5
c
c     core contribution to the energy is obtained from -ecore-,
c     and from modifications to the transformed 1e- integrals.
c     this effectively removes core orbitals from the computation.
c
      ntot = nact + ncorsv
      ntco = ncorsv
      norb = nact
      ncor = 0
      nsym = 2**igpdet
c
c     modify msta to get rid of core.
c
      do ii=1,nspace+1
         msta(ii) = msta(ii) - ncorsv
      enddo
c
c     compute the total number of determinants in this full class ci.
c     decide necessary double/integer working storage -ids- and -iis-
c
c     first need to store nspace sets of binomial coefficients for each
c     sub-space.  there is reason behind this, trust me.  just ask me.
c     lbst(i) will store where in x binomial arrays for space i start.
c     idim(i) will store the maximum number of alpha or beta electrons.
c
c
      call valfm(loadfm)
      lbst(1)  = loadfm + 1
      do ii=2,nspace
         idim(ii-1) = max(iama(ii-1),ibma(ii-1))
         lbst(ii) = lbst(ii-1)+((mnum(ii-1)+1)*(idim(ii-1)+1))
      enddo
      idim(nspace) = max(iama(nspace),ibma(nspace))
      last = lbst(nspace)+((mnum(nspace)+1)*(idim(nspace)+1))
      need1 = last - loadfm - 1
      if (some) write(iw,9000) need1
      call getfm(need1)
c
      do ii=1,nspace
         call binom8(x(lbst(ii)),mnum(ii),idim(ii))
      enddo 
c
c     now to make tables for the ormas problem.
c     make memory available for table info and
c     tables themselves.
c
      call valfm(loadfm)
      lbox1 = loadfm + 1
      lbox2 = lbox1 + nspace/nwdvar + 1
      lbox3 = lbox2 + nspace/nwdvar + 1
      lbox4 = lbox3 + nspace/nwdvar + 1
      lbox5 = lbox4 + nspace/nwdvar + 1
      libo = lbox5 + nspace/nwdvar + 1
      last = libo + ntot/nwdvar + 1
      need2 = last - loadfm - 1
      call getfm(need2)
c
c     work out how many alpha and beta groups there are.
c
      call totalco(x(lbox1),nspace,na,iama,iami,x(lbox2),itga)
      call totalco(x(lbox1),nspace,nb,ibma,ibmi,x(lbox2),itgb)
c
      if (some) write(iw,9010) itga,itgb
c
c     determine the total number of alpha and beta strings.
c
      call resetco(x(lbox1),nspace,na,iama,iami,x(lbox2))
      iast = 0
      do ii=1,itga
      call totdetg(x(lbox1),nspace,lbst,x(lbst(1)),need1,
     *    mnum,idim,itot)
      iast = iast + itot
      call pushco(x(lbox1),nspace,na,iama,iami,x(lbox2),iend)
      enddo
c
      call resetco(x(lbox1),nspace,nb,ibma,ibmi,x(lbox2))
      ibst = 0
      do ii=1,itgb
      call totdetg(x(lbox1),nspace,lbst,x(lbst(1)),need1,
     *    mnum,idim,itot)
      ibst = ibst + itot
      call pushco(x(lbox1),nspace,nb,ibma,ibmi,x(lbox2),iend)
      enddo
c
      if (some) write(iw,9020) iast,ibst
c
c     make storage to hold symmetry information of strings.
c     the labels below appear the same in routine maktabs
c     and this routine also explains what they are used for.
c
      call valfm(loadfm)
      lwrk = loadfm + 1
      lcoa = lwrk + 43/nwdvar + 1
      lcob = lcoa + (nsym*itga)/nwdvar + 1
      lgmul = lcob  + (nsym*itgb)/nwdvar + 1
      lktab = lgmul  + (nsym*nsym)/nwdvar + 1
      lcon  = lktab + nsym/nwdvar + 1
      lcon1 = lcon + na/nwdvar + 1
      lcon2 = lcon1 + na/nwdvar + 1
      lcon3 = lcon2 + na/nwdvar + 1
      landet = lcon3 + na/nwdvar + 1
      lbndet = landet + (itga*nspace)/nwdvar + 1
      nast = lbndet + (itgb*nspace)/nwdvar + 1
      nbst = nast + (itga+1)/nwdvar + 1
      lsyma = nbst + (itgb+1)/nwdvar + 1
      lsymb = lsyma + iast/nwdvar + 1
      lgcom = lsymb  + ibst/nwdvar + 1
      lspa = lgcom + (itga*itgb)/nwdvar + 1
      lspb = lspa + iast/nwdvar + 1
      ldisb = lspb + ibst/nwdvar + 1
      lsas = ldisb + (nsym*itgb*itga)/nwdvar + 1
      lsbs = lsas + ((nsym+1)*itga)/nwdvar + 1
      lsac = lsbs + ((nsym+1)*itgb)/nwdvar + 1
      lsbc = lsac + iast/nwdvar + 1
      last = lsbc + ibst/nwdvar + 1
      need3 = last - loadfm - 1
      if (some) write(iw,9030) need2 + need3
      call getfm(need3)
c
c     make integer tables for ormas. get number of dets, etc.
c
      if (some) call tsecnd(e0)
c
      call maktabs(iw,some,x(lbox1),x(lbox2),x(lbox3),nspace,na,nb,
     *            lbst,x(lbst(1)),
     *    need1,iama,iami,ibma,ibmi,mnum,idim,msta,molab(ncorsv+1),
     *            igpdet,kstsym,nsym,
     *            nact,x(lwrk),x(lktab),x(lgmul),
     *            x(lcon),
     *            x(lcoa),x(lcob),x(landet),x(lbndet),
     *            x(nast),x(nbst),x(lsyma),x(lsymb),x(lgcom),
     *            mini,maxi,x(lspa),x(lspb),x(ldisb),
     *            x(lsas),x(lsbs),x(lsac),x(lsbc),
     *            itga,itgb,iast,ibst,nci,na1ex,nb1ex,fdirct)
c
      if(some) then
         call tsecnd(e1)
         elap = e1 - e0
         e0 = e1
         write(iw,9107) elap
      endif
c
      if(some) then
         write(iw,9110) stsym,grpdet,sz,nci
      end if
c
c    end of integer tables.  now to determine the extra memory
c    requirements, double precision and integer.
c    used to do this by calling memci or something similar, i am
c    just going to do it here.
c
      m1 = nact
      m2 = (m1*m1+m1)/2
      m4 = (m2*m2+m2)/2
      l1 = num
      l2 = (l1*l1+l1)/2
      l3 = l1*l1
c
c     integral storage requirements first.
c     ormas double precision next.
c
      call valfm(loadfm)
      lsint2 = loadfm + 1
      lwk    = lsint2 + m4
      kcoeff = lwk    + l2
      kab    = kcoeff + maxp*nci
      kq     = kab    + maxp*nci
      kb     = kq     + nci
      kef    = kb     + 8*maxw1
      kf     = kef    + maxw1*maxw1
      kec    = kf     + (maxw1*(maxw1+1))/2
      kgr    = kec    + maxp*maxp
      kel    = kgr    + maxp
      ksp    = kel    + maxw1
      last   = ksp    + maxw1
c
c     ormas integer storage is last.
c
      iposa  = last
      ipera  = iposa  + (na*(nact-na)*nsym)/nwdvar + 1
      iind1  = ipera  + (na*(nact-na)*nsym)/nwdvar + 1
      igroa  = iind1  + (na*(nact-na)*nsym)/nwdvar + 1
      iwrk1  = igroa  + (na*(nact-na)*nsym)/nwdvar + 1
      iwrk2  = iwrk1  + (2*maxw1)/nwdvar + 1
      isd    = iwrk2  + maxw1/nwdvar + 1
      iso    = isd    + (na+nb)/nwdvar + 1
      index  = iso    + na/nwdvar + 1
      immc   = index  + ((nact*(nact+1))/2+1)/nwdvar + 1
      ihmcon = immc + nsym
c
c     leave space for the single excited storage.
c     note, if fdirct = .true. then nb1ex = itga.
c
      jb1gr = ihmcon + nstate
      jb1pe = jb1gr + nb1ex/nwdvar + 1
      jb1in = jb1pe + nb1ex/nwdvar + 1
      jb1po = jb1in + nb1ex/nwdvar + 1
      jb1sy = jb1po + nb1ex/nwdvar + 1
      jb1st = jb1sy + (nb*(nact-nb))/nwdvar + 1
c
      if (.not.fdirct) then
         last   = jb1st + ((ibst+1)*(nsym+1))/nwdvar + 1
         ibst2 = (ibst+1)*(nsym+1)
      else
         last = jb1st + 1
         nb1ex = itga
         ibst2 =  1
      endif
c
      need4  = last - loadfm - 1
      needci = need1 + need2 + need3 + need4
      if (some) then
          write(iw,9120) kcoeff-loadfm-1,
     *          need4-(kcoeff-loadfm-1),maxp,needci
          ifdm = ((nb1ex/nwdvar) + 1)*4 +
     *    (((ibst+1)*(nsym+1))/nwdvar)+1 + (nb*(nact-nb))/nwdvar + 1
          if (.not.fdirct) write(iw,9125) needci-ifdm+6
      endif
c
      call getfm(need4)
c
      if(exetyp.eq.check) then
         do ist=1,min(nstate,mxrt)
            spins(ist) = s
            estate(ist) = zero
         enddo
         lcivec = kcoeff
         call vclr(x(kcoeff),1,nstate*nci)
c        go to 450
      end if
c
      return
c
 9000 format(/5x,50(1h-)/
     *       5x,'            direct determinant ormas-ci '/
     *       5x,'           program written by joe ivanic'/
     *       5x,50(1h-)//
     *       1x,'storage of binomial coefficients requires',
     *       i12,' words')
 9010 format(/1x,'total number of alpha groups  = ',i12/
     *        1x,'total number of beta  groups  = ',i12)
 9020 format(/1x,'total number of alpha strings = ',i12/
     *        1x,'total number of beta  strings = ',i12)
 9030 format(/1x,'storage of tables requires              ',
     *       i12,' words')
 9107 format(/1x,'time for setting up table set 1 :',f13.1)
 9110 format(/1x,'the number of determinants having space symmetry ',a3/
     *        1x,'in point group ',a4,' with sz=',f5.1,' is',i15)
 9120 format(/1x,'integral storage requires        ',i12,' words'/
     *        1x,'extra ormas storage requires     ',i12,' words'/
     * 1x,'(extra ormas storage includes that for mxxpan =',i3,')'/
     *        1x,'total ormas calculation requires ',i12,' words')
 9125 format(/1x,'fully direct would require       ',i12,' words')
 9140 format(/1x,'..... done with ormas-ci computation .....')
 9150 format(1x,'ci computation did not converge, job cannot continue')
 9230 format(/1x,'1-particle density matrix in mo basis')
 9240 format(1x,'..... done with one particle density matrix .....')
 9250 format(1x,'determining reference weight for ',
     *          'davidson + q correction')
c
      end
c
c*module ormas1  *deck defcci
      subroutine defcci(nprint,clabel,ndm1,ndm2,
     *                  npri5,npri6,x,sint1,molab,dm1,dm2)
c
      implicit double precision(a-h,o-z)
c
      logical clabel
      logical some,goparr,dskwrk,maswrk,fdirct,qcorr
c
      parameter (mxrt=100, mxatm=2000)
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      common /enrgys/ enucr,eelct,etot,stot,ssquar,ecore,escf,eerd,
     *                e1,e2,ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
      common /infoa_f / nat,ich,mul,num,nqmt,ne,ma,mb,
     *                zan(mxatm),c(3,mxatm),ian(mxatm)
      common /iofile_f/ ir,iw,ip,is,ijkt,idaf,nav,ioda(950)
      common /machin_f/ nwdvar,maxfm,maxsm,limfm,limsm
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
      common /wfnopt_f/ scftyp,cityp,dftype,cctyp,mplevl,mpctyp
c
c     addresses needed by defcci calls in here
c
      common /defaddr/ econst,lsint2,lwk,m1,m2,m4,ksp,kel,nsym,lbox1
     *,lbox2,lbox3,lbox4,lbox5,kcoeff,kab,kq,kb,kef,kf,kec,kgr,iposa,
     * ipera,iind1,igroa,iwrk1,iwrk2,isd,iso,index,immc,ihmcon,lgmul,
     * lktab,lcon,lcon1,lcon2,lcon3,landet,lbndet,nast,nbst,lsyma,
     * lsymb,lgcom,lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,itga,itgb,
     * iast,ibst,nb1ex,ibst2,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy,
     * istat,need1,need2,need3,need4,ntot,ntco
c
c     gamess-uk memory management, replaces /fmcom/
c     common /fmcom / x(1)
c
c     x(lsint1) replaced by sint1 passed from ormas to avoid i/o
c     x(libo) replaced by molab passed from fciput to avoid i/o
c
      dimension x(*), sint1(*), molab(*), dm1(*), dm2(*)
c
      parameter (zero=0.0d+00, one=1.0d+00)
c
      data amcscf/8hmcscf   /
      data check/8hcheck   /
c
c     flag for 1st masscf iter
c
      logical first_fnr
      common /flag_fnr/ first_fnr
c
c     problems up to 10,10 best executed sequentially
cgdf  ncibig could be input ?
c
      parameter (ncibig = 63504)
c
c     core contribution to the energy is obtained from -ecore-,
c     and from modifications to the transformed 1e- integrals.
c     this effectively removes core orbitals from the computation.
c
      econst = ecore + enucr
c
c     print flag on for first iteration
c
      some = maswrk  .and.  nprint.ne.-5 .and. first_fnr
c
c     -- obtain 2 e- transformed integrals over active orbitals --
c     calling argument -clabel- governs whether transformed integrals
c     include the core orbitals or not.  it is assumed
c     that no core integrals are in -ijkt-, so this variable tells if
c     the active orbitals start from 1,2,3... or ncorsv+1,ncorsv+2,...
c
      ncore = 0
      if(clabel) ncore=ncorsv
      call rdci12(x(lsint2),ncore,m1,m2,m4,x(lwk))
c
c     switch off parallel for small cases...
c
      if (nci.lt.ncibig) goparr = .false.
c
c     ----- compute full class ci wavefunction -----
c
      call dafcci(iw,some,econst,istat,
     *              sint1 ,x(lsint2),m2,m4,nact,nci,na,nb,
     *           nstate,kst,maxp,maxw1,niter,crit,x(lbst(1)),need1,
     *           x(ksp),x(kel),nsym,molab,
     *           x(lbox1),x(lbox2),x(lbox3),x(lbox4),x(lbox5),
     *           x(kcoeff),x(kab),x(kq),x(kb),x(kef),x(kf),
     *           x(kec),x(kgr),
     *           x(iposa),x(ipera),x(iind1),x(igroa),
     *           x(iwrk1),x(iwrk2),
     *           x(isd),x(iso),
     *           x(index),x(immc),x(ihmcon),
     *           x(lgmul),x(lktab),
     *           x(lcon),x(lcon1),x(lcon2),x(lcon3),
     *           x(landet),x(lbndet),x(nast),x(nbst),
     *           x(lsyma),x(lsymb),x(lgcom),x(lspa),x(lspb),
     *           x(ldisb),x(lsas),x(lsbs),x(lsac),x(lsbc),
     *           itga,itgb,iast,ibst,
     *           nb1ex,ibst2,
     *           x(jb1gr),x(jb1pe),x(jb1in),x(jb1po),x(jb1st),
     *           x(jb1sy))
c
c     ...otherwise masscf assumed to be parallel
c
      goparr = .true.
c
      do i=1,min(nstate,mxrt)
         estate(i) = x(kel-1+i)+econst
         spins(i)  = x(ksp-1+i)
      enddo
c
c     save energy quantities
c
  450 continue
      etot = estate(iroot)
      eelct = etot - enucr
      stot = spins(iroot)
      ssquar = stot*(stot+one)
      statn = nstate
c
c     if asked for a davidson + q correction, then need to determine
c     weight of the reference in the total mr-cisd wavefunction
c
      if(qcorr) then
         write(iw,9250)
         call refwe(x(kcoeff),x(nast),itga,itgb,x(lsyma),
     *              iast,x(lgcom),x(lsbs),nsym,x(lktab))
      endif
c
c     copy ci vectors to ab as the ones in kcoeff are destroyed
c     during printing.
c
      call cicopy(x(kcoeff),x(kab),nci*nstate)
c
c     print results of the ci calculation
c
      call masprt(iw,some,
     *       x(kab),x(nast),itga,itgb,x(lsyma),
     *       iast,ibst,x(lgcom),x(lsbs),nsym,x(lsbc),
     *       x(lbox1),x(lbox2),x(lbox3),x(lcon1),x(lcon2),
     *       x(lktab))
c
      if(some) write(iw,9140)
      if(some) call timit(1)
c
c     determine 1st order density matrix if requested
c
 1234 if (ndm1.gt.0) then
         call masdm1(iw,npri5,
     *              dm1   ,m2,nact,nci,na,nb,iroot,
     *           x(lbst(1)),need1,x(kcoeff+(nci*(iroot-1))),
     *           x(index),nsym,molab,
     *           x(lbox1),x(lbox2),x(lbox4),x(lbox5),
     *           x(lktab),x(lcon),x(lcon1),x(lcon2),x(lcon3),
     *           x(landet),x(lbndet),x(nast),x(nbst),
     *           x(lsyma),x(lsymb),x(lgcom),x(lspa),x(lspb),
     *           x(ldisb),x(lsas),x(lsbs),x(lsac),
     *           itga,itgb,iast,ibst)
         if(some  .and.  nflgdm(ist).eq.2) then
            write(iw,9230)
            call prtri(x(lsint1),norb)
         end if
         if (some) write(iw,9240)
         if (some) call timit(1)
      endif
c
c     determine if state averaged 1e- and 2e- density matrices requested
c
      if (ndm2.gt.0) then
         if (.not.fdirct) then
            idim1 = nsym+1
            idim2 = ibst+1
         else
            idim1 = 1
            idim2 = 1
         endif
c
         call masdm2(iw,npri6,iwts,wstate,spins,ipures,
     *           s,nstate,grpdet,ncorsv,
     *              dm1   ,   dm2   ,m2,m4,nact,nci,na,nb,
     *           x(kcoeff),x(kab),
     *           x(lbst(1)),need1,
     *           x(index),nsym, molab ,
     *           x(lbox1),x(lbox2),x(lbox3),x(lbox4),x(lbox5),
     *           x(lgmul),x(lktab),
     *           x(lcon),x(lcon1),x(lcon2),x(lcon3),
     *           x(landet),x(lbndet),x(nast),x(nbst),
     *           x(lsyma),x(lsymb),x(lgcom),x(lspa),x(lspb),
     *           x(ldisb),x(lsas),x(lsbs),x(lsac),x(lsbc),
     *           itga,itgb,iast,ibst,
     *           x(iposa),x(ipera),x(iind1),x(igroa),x(immc),
     *           nb1ex,
     *           x(jb1gr),x(jb1pe),x(jb1in),x(jb1po),x(jb1st),
     *           idim1,idim2,x)
      endif
c
      if(exetyp.ne.check  .and.  istat.ne.0 .and. scftyp.ne.amcscf) then
         if(maswrk) write(iw,9150)
         call abrt
         stop
      end if
c
      return
c
 9000 format(/5x,50(1h-)/
     *       5x,'            direct determinant ormas-ci '/
     *       5x,'           program written by joe ivanic'/
     *       5x,50(1h-)//
     *       1x,'storage of binomial coefficients requires',
     *       i12,' words')
 9010 format(/1x,'total number of alpha groups  = ',i12/
     *        1x,'total number of beta  groups  = ',i12)
 9020 format(/1x,'total number of alpha strings = ',i12/
     *        1x,'total number of beta  strings = ',i12)
 9030 format(/1x,'storage of tables requires              ',
     *       i12,' words')
 9107 format(/1x,'time for setting up table set 1 :',f13.1)
 9110 format(/1x,'the number of determinants having space symmetry ',a3/
     *        1x,'in point group ',a4,' with sz=',f5.1,' is',i15)
 9120 format(/1x,'integral storage requires        ',i12,' words'/
     *        1x,'extra ormas storage requires     ',i12,' words'/
     * 1x,'(extra ormas storage includes that for mxxpan =',i3,')'/
     *        1x,'total ormas calculation requires ',i12,' words')
 9125 format(/1x,'fully direct would require       ',i12,' words')
 9140 format(/1x,'..... done with ormas-ci computation .....')
 9150 format(1x,'ci computation did not converge, job cannot continue')
 9230 format(/1x,'1-particle density matrix in mo basis')
 9240 format(1x,'..... done with one particle density matrix .....')
 9250 format(1x,'determining reference weight for ',
     *          'davidson + q correction')
c
      end
c
c*module ormas1  *deck defcci_end
c
      subroutine defcci_end(iw,nprint,x)
c
      implicit double precision(a-h,o-z)
c
      logical some,goparr,dskwrk,maswrk 
c
      parameter (mxrt=100)
c
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      common /defaddr/ econst,lsint2,lwk,m1,m2,m4,ksp,kel,nsym,lbox1
     *,lbox2,lbox3,lbox4,lbox5,kcoeff,kab,kq,kb,kef,kf,kec,kgr,iposa,
     * ipera,iind1,igroa,iwrk1,iwrk2,isd,iso,index,immc,ihmcon,lgmul,
     * lktab,lcon,lcon1,lcon2,lcon3,landet,lbndet,nast,nbst,lsyma,
     * lsymb,lgcom,lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,itga,itgb,
     * iast,ibst,nb1ex,ibst2,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy,
     * istat,need1,need2,need3,need4,ntot,ntco
c
      dimension x(*)
c
      some = maswrk  .and.  nprint.ne.-5  .and.  nprint.ne.-23
c
c     print the final ci vectors
c
      call masprt(iw,some,
     *       x(kab),x(nast),itga,itgb,x(lsyma),
     *       iast,ibst,x(lgcom),x(lsbs),nsym,x(lsbc),
     *       x(lbox1),x(lbox2),x(lbox3),x(lcon1),x(lcon2),
     *       x(lktab))
c
c     deallocations and resets from end of defcci
c
      call retfm(need4)
      call retfm(need3)
      call retfm(need2)
      call retfm(need1)
c
      do ii=1,nspace+1
         msta(ii) = msta(ii) + ncorsv
      enddo
c 
      end
c
c*module ormas1   *deck reset
c     ---------------------------------------------
      subroutine reset(ibox,nbox,nr,mb)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mb(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mb(i) = maximum capacity of box i.
c
c     subroutine makes first group in the full list,
c     assuming that minimum(i) = 0
c
      nt = nr
      do i=1,nbox
         ibox(i) = min(mb(i),nt)
         nt = nt - ibox(i)
      enddo
c
      return
      end
c
c*module ormas1   *deck push
c     ---------------------------------------------
      subroutine push(ibox,nbox,nr,mb,iend)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mb(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mb(i) = maximum capacity of box i.
c
c     subroutine returns next group in the full list.
c     if (iend.eq.0) then last group was input in ibox.
c     if (iend.eq.1) then push was successful.
c
      iend = 0
      nu = nr
      do 100 k=nbox-1,1,-1
         if (ibox(k).gt.0) then
         do 200 m=k+1,nbox
            if (ibox(m).lt.mb(m)) then
               ibox(k) = ibox(k)-1
               do 300 ii=1,m-1
                  nu = nu - ibox(ii)
  300          continue
               call reset(ibox(m),nbox-m+1,nu,mb(m))
               iend = 1
               return
            endif
  200    continue
         endif
  100 continue
c
      return
      end
c
c*module ormas1   *deck total
c     ---------------------------------------------
      subroutine total(ibox,nbox,nr,mb,itot)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mb(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.  scratch.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mb(i) = maximum capacity of box i.
c
c     returned -----
c
c     itot = number of groups.
c
      iy = 0
      do ii=1,nbox
         iy = iy + mb(ii)
      enddo
      if (iy.lt.nr) then
         itot = 0
         return
      endif
c
      itot = 0
      call reset(ibox,nbox,nr,mb)
c
  100 continue
      itot = itot + 1
      call push(ibox,nbox,nr,mb,iend)
      if (iend.eq.1) go to 100
c
      return
      end
c
c*module ormas1   *deck getbox
c     ---------------------------------------------
      subroutine getbox(ibox,icon,nr,msta,nbox)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),icon(nr),msta(nbox+1)
c
      ista=1
      do ii=1,nbox
         ibox(ii)=0
      enddo
      do 10 ii=1,nr
         ia=icon(ii)
         do jj=ista,nbox
            if (msta(jj+1).gt.ia) then
               ibox(jj)=ibox(jj)+1
               ista=jj
               go to 10
            endif
         enddo
   10 continue
c
      return
      end
c
c*module ormas1   *deck posit_f
c     ---------------------------------------------
      subroutine posit_f(ibox,nbox,nr,mb,jbox,ipos)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mb(nbox),jbox(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.  scratch.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mb(i) = maximum capacity of box i.
c     jbox = group to have it position determined
c
c     returned -----
c
c     ipos = position of jbox in full set.
c
      ipos=0
      call reset(ibox,nbox,nr,mb)
  100 continue
      ipos = ipos + 1
      if (icompa(ibox,jbox,nbox).eq.0) return
      call push(ibox,nbox,nr,mb,iend)
      if (iend.ne.1) then
         ipos=0
         return
      endif
      go to 100
c
      end
c
c*module ormas1   *deck icompa
c     ---------------------------------------------
      integer function icompa(ia,ib,n)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ia(n),ib(n)
c
c     if ia = ib, icompa=0, =1 otherwise.
c
      icompa=1
      do i=1,n
         if (ia(i).ne.ib(i)) return
      enddo
      icompa=0
c
      return
      end
c
c*module ormas1   *deck resetco
c     ---------------------------------------------
      subroutine resetco(ibox,nbox,nr,mab,mib,isc)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mab(nbox),mib(nbox),isc(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mab(i) = maximum capacity of box i.
c     mib(i) = minimum capacity of box i.
c     isc(i) = scratch array
c
c     subroutine makes first group in the full list,
c     assuming that mib(i).ne.0
c
      np = nr
      do i=1,nbox
         isc(i) = mab(i) - mib(i)
         np = np - mib(i)
      enddo
      call reset(ibox,nbox,np,isc)
      do i=1,nbox
         ibox(i) = ibox(i) + mib(i)
      enddo
c
      return
      end
c
c*module ormas1   *deck pushco
c     ------------------------------------------------
      subroutine pushco(ibox,nbox,nr,mab,mib,isc,iend)
c     ------------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mab(nbox),mib(nbox),isc(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mab(i) = maximum capacity of box i.
c     mib(i) = minimum capacity of box i.
c     isc(i) = scratch array
c
c     subroutine returns next group in the full list.
c     if (iend.eq.0) then last group was input in ibox.
c     if (iend.eq.1) then push was successful.
c
      np = nr
      do i=1,nbox
         isc(i) = mab(i) - mib(i)
         ibox(i) = ibox(i) - mib(i)
         np = np - mib(i)
      enddo
      call push(ibox,nbox,np,isc,iend)
      do i=1,nbox
         ibox(i) = ibox(i) + mib(i)
      enddo
c
      return
      end
c
c*module ormas1   *deck totalco
c     -------------------------------------------------
      subroutine totalco(ibox,nbox,nr,mab,mib,isc,itot)
c     -------------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mab(nbox),mib(nbox),isc(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.  scratch.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mab(i) = maximum capacity of box i.
c     mib(i) = minimum capacity of box i.
c     isc(i) = scratch array
c
c     returned -----
c
c     itot = number of groups.
c
      np = nr
      do i=1,nbox
         isc(i) = mab(i) - mib(i)
         np = np - mib(i)
      enddo
      call total(ibox,nbox,np,isc,itot)
c
      return
      end
c
c*module ormas1   *deck positco
c     ------------------------------------------------------
      subroutine positco(ibox,nbox,nr,mab,mib,isc,jbox,ipos)
c     ------------------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),mab(nbox),jbox(nbox)
      dimension mib(nbox),isc(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.  scratch.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     mab(i) = maximum capacity of box i.
c     mib(i) = minimum capacity of box i.
c     isc(i) = scratch array.
c     jbox = group to have it position determined
c
c     returned -----
c
c     ipos = position of jbox in full set.
c
      np = nr
      do i=1,nbox
         isc(i) = mab(i) - mib(i)
         jbox(i) = jbox(i) - mib(i)
         np = np - mib(i)
      enddo
      call posit_f(ibox,nbox,np,isc,jbox,ipos)
      do i=1,nbox
         jbox(i) = jbox(i) + mib(i)
      enddo
c
      return
      end
c
c*module ormas1  *deck binom8
c     --------------------------
      subroutine binom8(ifa,n,m)
c     --------------------------
      implicit double precision(a-h,o-z)
      integer ifa(0:n,0:m)
c
c     returns all binomial numbers (i,j) for i=0,n and j=0,min(i,m)
c     in ifa. the binomial number (i,j) is stored in ifa(i,j).
c
      do 11 ii=0,n
         do 12 jj=0,m
            ifa(ii,jj) = 0
   12    continue
   11 continue
c
      do 13 ii=0,n
         ifa(ii,0)  = 1
   13 continue
c
      do 14 ii=0,m
         ifa(ii,ii) = 1
   14 continue
c
      do 113 iy = 2, n
         iz = min(iy-1,m)
         do 114 ix = 1, iz
            ifa(iy,ix) = ifa(iy-1,ix-1) + ifa(iy-1,ix)
  114    continue
  113 continue
c
c      do ii=0,n
c         write(6,'(100i12)') (ifa(ii,jj),jj=0,m)
c      enddo
c
      return
      end
c
c*module ormas1  *deck totdetg
c     --------------------------------------------------------
      subroutine totdetg(ibox,nbox,lbst,x,nx,mnum,idim,itot)
c     --------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ibox(nbox),lbst(nbox)
      dimension x(nx)
      dimension mnum(nbox),idim(nbox)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c     lbst(i) = where binomial coefficient for space i
c     start in x.
c     nx is the size of x.  this is very specific.
c     mnum(i) = number of orbitals in space i
c     idim(i) = row dimension of binomial coefficients for space i.
c
c     returned -----
c
c     itot = total number of determinants in group ibox.
c
      itot = 1
      ista1 = lbst(1)
      do ii=1,nbox
         ista2 = lbst(ii) - ista1 + 1
         itot = itot *
     *   ispadet(x(ista2),mnum(ii),idim(ii),ibox(ii))
      enddo
c
      return
      end
c
c*module ormas1  *deck ispadet
c     --------------------------------------------
      integer function ispadet(ifa,norb,ndim,nele)
c     --------------------------------------------
      implicit integer(a-z)
      dimension ifa(0:norb,0:ndim)
c
      ispadet = ifa(norb,nele)
      return
      end
c
c*module ormas1  *deck moveup
c     ---------------------------------------------
      subroutine moveup(con,nele,norb)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension con(nele)
c
      if (nele.lt.1) return
c
      if (con(nele).ne.norb) then
         con(nele) = con(nele)+1
         return
      else
         do 50 i=nele-1,1,-1
            if (con(i+1)-con(i).ne.1) then
               con(i) = con(i) + 1
               do 40 j=i+1,nele
                  con(j) = con(j-1) + 1
   40          continue
               return
            endif
   50    continue
c
      endif
      return
      end
c
c*module ormas1   *deck resetde
c     ---------------------------------------------
      subroutine resetde(ibox,nbox,nr,msta,icon)
c     ---------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),msta(nbox)
      dimension icon(nr)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     msta(i) = where orbital space i starts.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c
c  ---------------
c  returned
c
c     icon = returned reset string.
c
      iste = 0
      do ii=1,nbox
         isto = msta(ii)-1
         inum = ibox(ii)
         do jj=1,inum
            icon(jj+iste) = jj+isto
         enddo
         iste = iste + inum
      enddo
c
      return
      end
c
c*module ormas1   *deck moveup2
c     ----------------------------------------------------
      subroutine moveup2(ibox,nbox,nr,msta,icon)
c     ----------------------------------------------------
      implicit integer(a-z)
      dimension ibox(nbox),msta(nbox+1)
      dimension icon(nr)
c
c     ibox(i),i=1,nbox stores number of particles
c                          in each box i.
c     msta(i) = where orbital space i starts.
c     nbox = total number of boxes.
c     nr   = total number of particles.
c
c   ---------------------
c   returned
c
c     icon = advanced total string.
c
      nc = nr
      do ii=nbox,1,-1
         iele=ibox(ii)
         if (iele.ne.0) then
            nb = nc-iele+1
            imst=msta(ii+1)
            iend = imst-iele
            if (icon(nb).ne.iend) then
            call moveup(icon(nb),iele,imst-1)
            if (ii.ne.nbox)
     *   call resetde(ibox(ii+1),nbox-ii,nr-nc,msta(ii+1),icon(nc+1))
            return
            endif
         endif
         nc = nc - iele
      enddo
c
      return
      end
c
c*module ormas1  *deck idpos1
c     -----------------------------------------------
      integer function idpos1(nact,noel,con,ifa,ndim)
c     -----------------------------------------------
      implicit integer(a-z)
      dimension con(noel)
      dimension ifa(0:nact,0:ndim)
c
c     nact   = no. of orbitals.
c     noel   = no. of electrons.
c     con(i) = orbital occupied by electron i.
c     ifa    = binomial coefficients.
c     ndim   = row dimension of the binomial coefficients.
c
c     returns position of string con in full string list.
c
      ipos1 = 0
      idpos1 = 1
      do 33 i=1,noel
         do 55 j=ipos1+1,con(i)-1
            idpos1 = idpos1 + ifa(nact-j,noel-i)
   55    continue
         ipos1 = con(i)
   33 continue
c
      return
      end
c
c*module ormas1  *deck idpost
c     --------------------------------------------------------
      subroutine idpost(icon,na,ibox,nbox,msta,idim,x,nx,lbst,
     *                  lnum,jcon,ipos)
c     --------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension icon(na),ibox(nbox),msta(nbox)
      dimension idim(nbox),lbst(nbox)
      dimension x(nx)
      dimension lnum(nbox),jcon(na)
c
c     icon(i) = orbital occupied by electron i.
c     na      = no. of electrons.
c     ibox(i),i=1,nbox stores number of electrons
c                          in each box i.
c     nbox = total number of boxes.
c     msta(i) = where orbital space i starts.
c     idim(i) = row dimension of binomial coefficients for space i.
c     lbst(i) = where binomial coefficient for space i
c     start in x.
c     nx is the size of x.  this is very specific.
c     lnum(i) = total number of strings for space i.
c     jcon = scratch array.
c
c     returns position of string icon in group string list.
c
      ipos = 1
      ista1 = lbst(1)
      ieles = 0
      do ii=1,nbox
         ne = ibox(ii)
         ista2 = lbst(ii) - ista1 + 1
         nact = msta(ii+1) - msta(ii)
         do jj=1,ne
            jcon(jj) = icon(jj+ieles)-msta(ii)+1
         enddo
         ipos1 = idpos1(nact,ne,jcon,x(ista2),idim(ii)) - 1
         do jj=ii+1,nbox
            ipos1 = ipos1 * lnum(jj)
         enddo
         ipos = ipos + ipos1
         ieles = ieles + ne
      enddo
c
      return
      end
c
c*module ormas1  *deck maktabs
c     ------------------------------------------------------------------
      subroutine maktabs(iw,some,lbox1,lbox2,lbox3,nspace,na,nb,
     *            lbst,x,
     *            nx,iama,iami,ibma,ibmi,mnum,idim,msta,ibo,
     *            idsym,isym1,nsym,
     *            nact,lwrk,ktab,lgmul,
     *            lcon,lcoa,lcob,
     *            landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *            mini,maxi,lspa,lspb,ldisb,
     *            lsas,lsbs,lsac,lsbc,
     *            itga,itgb,iast,ibst,nci,na1ex,nb1ex,fdirct)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      logical some
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace),lbst(nspace)
      dimension x(nx)
      dimension iama(nspace),iami(nspace),ibma(nspace),ibmi(nspace)
      dimension mnum(nspace),idim(nspace),msta(nspace+1),ibo(nact)
      dimension lwrk(43),ktab(nsym),lgmul(nsym,nsym)
      dimension lcon(na)
      dimension lcoa(nsym,itga),lcob(nsym,itgb)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst)
      dimension lgcom(itgb,itga)
      dimension mini(nspace),maxi(nspace)
      dimension lspa(iast),lspb(ibst)
      dimension ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast),lsbc(ibst)
c
      logical fdirct
c
c 1/
c   make symmetry tables
c
      if (idsym.gt.0) then
      call gtab(idsym,isym1,ktab,lwrk(1),lwrk(4),lwrk(7),lwrk(10))
      call gmul(idsym,lgmul,lwrk(1),lwrk(4),lwrk(7),lwrk(10))
      else
      call gtab(1,1,ktab,lwrk(1),lwrk(4),lwrk(7),lwrk(10))
      call gmul(1,lgmul,lwrk(1),lwrk(4),lwrk(7),lwrk(10))
      endif
c
      do ii=1,itga
         do jj=1,nsym
            lcoa(jj,ii) = 0
         enddo
      enddo
c
      do ii=1,itgb
         do jj=1,nsym
            lcob(jj,ii) = 0
         enddo
      enddo
c
c  2/
c    make landet, lbndet.  landet(i,j) says how many alpha strings
c    there are for group j, space i.  analogous for lbndet(i,j).
c
c    loop through all alpha groups:
c
      ista1 = lbst(1)
      call resetco(lbox1,nspace,na,iama,iami,lbox2)
      do ii=1,itga
         do jj=1,nspace
            ista2 = lbst(jj) - ista1 + 1
            landet(jj,ii)=ispadet(x(ista2),mnum(jj),idim(jj),lbox1(jj))
         enddo
      call pushco(lbox1,nspace,na,iama,iami,lbox2,iend)
      enddo
c
c    loop through all beta groups, do same as above essentially.
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox2)
      do ii=1,itgb
         do jj=1,nspace
            ista2 = lbst(jj) - ista1 + 1
            lbndet(jj,ii)=ispadet(x(ista2),mnum(jj),idim(jj),lbox1(jj))
         enddo
      call pushco(lbox1,nspace,nb,ibma,ibmi,lbox2,iend)
      enddo
c
c  3/
c    make nast, nbst.  nast(i) says where alpha strings of group i
c    start in full string list - 1.  analogous for nbst.
c
      nast(1) = 0
      nbst(1) = 0
      do ii=1,itga
         itot = 1
         do jj=1,nspace
            itot = itot * landet(jj,ii)
         enddo
         nast(ii+1) = nast(ii) + itot
      enddo
c
      do ii=1,itgb
         itot = 1
         do jj=1,nspace
            itot = itot * lbndet(jj,ii)
         enddo
         nbst(ii+1) = nbst(ii) + itot
      enddo
c
c  4/
c       make lgcom.  if lgcom(i,j).ne.0 then beta group i
c       and alpha group j are compatible.
c
      icomp = 0
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      do jj=1,itga
         call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
         do ii=1,itgb
            lgcom(ii,jj) = 0
            do kk=1,nspace
               ioc = lbox1(kk) + lbox2(kk)
               if (ioc.gt.maxi(kk).or.ioc.lt.mini(kk)) go to 100
            enddo
c
            lgcom(ii,jj) = 1
            icomp = icomp + 1
  100       call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
         enddo
         call pushco(lbox1,nspace,na,iama,iami,lbox3,iend)
      enddo
c
      if (some) write(iw,9000) itga*itgb,icomp
c
c  5/
c    a) make lsymb.  lsymb(i) is symmetry of beta string i.
c       make lcob.  lcob(i,j) is number of determinants
c       of symmetry i in group j.
c       make lspb.  lspb(i) is beta string i's symmetry position in
c       it's own group.
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox2)
c
c    loop through groups.
c
      icount = 0
      do ii=1,itgb
c
c    loop through all strings in group ii.
c
         call resetde(lbox1,nspace,nb,msta,lcon)
c
         itot = 1
         do jj=1,nspace
            itot = itot * lbndet(jj,ii)
         enddo
c
         do kk=1,itot
            icount = icount + 1
            call getsym1(iw,lcon,nact,nb,ibo,idsym,isym,
     *      lwrk(1),lwrk(4),lwrk(7),lwrk(10))
            lsymb(icount) = isym
            lcob(isym,ii) = lcob(isym,ii) + 1
            lspb(icount) = lcob(isym,ii)
c
            call moveup2(lbox1,nspace,nb,msta,lcon)
         enddo
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox2,iend)
      enddo
c
c    b) - make lsyma.  lsyma(i) is symmetry of alpha string i.
c       - make lcoa.  lcoa(i,j) is number of determinants
c         of symmetry i in group j.
c       - make lspa.  lspa(i) is where alpha string i starts in
c         the full list of determinants.  each alpha string is coupled
c         with relevant beta strings (in beta group and symmetry order).
c       - make ldisb.  ldisb(isym,i,j) says where beta group i, and
c         symmetry isym starts when coupled to alpha group j.
c
      call resetco(lbox1,nspace,na,iama,iami,lbox2)
c
c    loop through groups.
c
      icount = 0
      nci = 0
      do ii=1,itga
c
c    loop through all strings in group ii.
c
         call resetde(lbox1,nspace,na,msta,lcon)
         itot = 1
         do jj=1,nspace
            itot = itot * landet(jj,ii)
         enddo
c
         do kk=1,itot
            icount = icount + 1
            call getsym1(iw,lcon,nact,na,ibo,idsym,isym,
     *      lwrk(1),lwrk(4),lwrk(7),lwrk(10))
            lsyma(icount) = isym
            lcoa(isym,ii) = lcoa(isym,ii) + 1
c
            lspa(icount) = nci
            jsym = ktab(isym)
            do 200 ll=1,itgb 
               if (lgcom(ll,ii).eq.0) go to 200
               nci = nci + lcob(jsym,ll)
  200       continue
c
            call moveup2(lbox1,nspace,na,msta,lcon)
         enddo
c
         do kk=1,nsym
            lwrk(kk) = 0
         enddo
         do 300 jj=1,itgb
            do kk=1,nsym
               ldisb(kk,jj,ii) = 0
            enddo
            if (lgcom(jj,ii).eq.0) go to 300
            do kk=1,nsym
               ldisb(kk,jj,ii) = lwrk(kk)
               lwrk(kk) = lwrk(kk) + lcob(kk,jj)
            enddo
  300    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox2,iend)
      enddo
c
c  6/
c    a) make lsas,lsbs,and lsac, lsbc.
c       lsas(i,j) says where alpha strings
c       of group j and symmetry i, start in lsac.
c       analogous for lsbs and lsbc.
c
      ipla = 1
      do ii=1,itga
         do jj=1,nsym
            lsas(jj,ii) = ipla
            ipla = ipla + lcoa(jj,ii)
         enddo
         lsas(nsym+1,ii) = ipla
      enddo
c
      iplb = 1
      do ii=1,itgb
         do jj=1,nsym
            lsbs(jj,ii) = iplb
            iplb = iplb + lcob(jj,ii)
         enddo
         lsbs(nsym+1,ii) = iplb
      enddo
c
      do ii=1,itga
         itot = 1
         do jj=1,nsym
            lwrk(jj) = 0
         enddo
         do jj=1,nspace
            itot = itot * landet(jj,ii)
         enddo
         idisa = nast(ii)
         do kk=1,itot
            jsym = lsyma(idisa+kk)
            lwrk(jsym) = lwrk(jsym)+1
            lsac(lsas(jsym,ii)+lwrk(jsym)-1) = kk
         enddo
      enddo
c
      do ii=1,itgb
         itot = 1
         do jj=1,nsym
            lwrk(jj) = 0
         enddo
         do jj=1,nspace
            itot = itot * lbndet(jj,ii)
         enddo
         idisb = nbst(ii)
         do kk=1,itot
            jsym = lsymb(idisb+kk)
            lwrk(jsym) = lwrk(jsym)+1
            lsbc(lsbs(jsym,ii)+lwrk(jsym)-1) = kk
         enddo
      enddo
c
cc *******************************************
cc  test of string ordering system.
cc  do not delete this code as it is excellent
cc  for possible debugs.  joe
cc *******************************************
cc
c      write(6,*) '***',nci
c      call resetco(lbox1,nspace,na,iama,iami,lbox3)
cc
cc    loop through groups.
cc
c      icount = 0
c      do ii=1,itga
cc
cc    loop through all strings in group ii.
cc
c         call resetde(lbox1,nspace,na,msta,lcon)
cc
c         do kk=nast(ii)+1,nast(ii+1)
c            jasym = lsyma(kk)
c            ksym = ktab(jasym)
cc
cc    loop through all applicable beta groups
cc
c         call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
cc
c         do iib = 1,itgb
c            if (lgcom(iib,ii).eq.0) go to 600
cc
c            call resetde(lbox2,nspace,nb,msta,lcon1)
c            ista = 1
c            do kkb=lsbs(ksym,iib),lsbs(ksym+1,iib)-1
c               ibst = lsbc(kkb)
c               do iiz=ista,ibst-1
c                  call moveup2(lbox2,nspace,nb,msta,lcon1)
c               enddo
c               ista = ibst
c               icount = icount + 1
cc
c         call positco(lwrk,nspace,na,iama,iami,lbox3,lbox1,iga)
c         call positco(lwrk,nspace,nb,ibma,ibmi,lbox3,lbox2,igb)
c
c            call idpost(lcon,na,lbox1,nspace,msta,idim,x,nx,lbst,
c     *      landet(1,iga),lcon3,iposa)
c            call idpost(lcon1,nb,lbox2,nspace,msta,idim,x,nx,lbst,
c     *      lbndet(1,igb),lcon3,iposb)
c            iposa = iposa + nast(iga)
c            iposb = iposb + nbst(igb)
c            ibsym = lsymb(iposb)
cc
cc here is the position formula.  complicated yes but i don't know a
cc better way yet.
cc
c            itotp = lspa(iposa)+ldisb(ibsym,igb,iga)+lspb(iposb)
cc
c            write(6,*) itotp,icount
c            if (icount.ne.itotp) then
c                stop
c            endif
cc
c            enddo
cc
c  600       call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
c         enddo
cc
c            call moveup2(lbox1,nspace,na,msta,lcon)
c         enddo
cc
c         call pushco(lbox1,nspace,na,iama,iami,lbox3,iend)
c      enddo
cc
cc ********************************
cc end of test code.
cc ********************************
c
      na1ex = 0
      nb1ex = itga
      if (fdirct) return
      nb1ex = 0
c
c ****************************************
c   determine total number of single beta excitations
c   where b' > b.
c ****************************************
c
c    loop through all beta groups
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox3)
c
      do 1000 iib = 1,itgb
c
         call resetde(lbox1,nspace,nb,msta,lcon)
c
         do 900 kkb=nbst(iib)+1,nbst(iib+1)
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over all single excitations, checking to see
c  if it is valid.
c
c  loop spaces to excite electrons from.
c
            iebs = nb+1
            do 890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) go to 890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
               do 885 ib1=iebe,iebs,-1
                  io1 = lcon(ib1)
c
                  igbe = iebe - lbox1(ispb1)
c
c  loop over possible spaces to excite into.
c
               do 880 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) go to 870
                  if (lbox2(ispb2).gt.ibma(ispb2)) go to 870
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = lcon(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = lcon(igba)-1
                  endif
c
c  loop over gaps
c
                  do 860 igap=igba,igbe+1
c
                     do 850 jj=ista,iend
c
                     nb1ex = nb1ex + 1
c
c  ****** all the work has to be done in here. *******
c
  850             continue
c
                  ista = lcon(igap)+1
                  iend = lcon(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
  860             continue
c
  870             lbox2(ispb2) = lbox2(ispb2) - 1
  880          continue
c
  885          continue
c
               lbox2(ispb1) = lbox2(ispb1) + 1
  890       continue
c
            call moveup2(lbox1,nspace,nb,msta,lcon)
  900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox3,iend)
c
 1000 continue
c
      if (some) write(iw,9010) nb1ex
      return
c
 9000 format(/1x,'total number of alpha-beta group ',
     *           'combinations = ',i12
     *       /1x,'of these the allowed number of ',
     *           'combinations   = ',i12)
c
 9010 format(/1x,'total number of (b -> b) where  (b > b)  = ',i12)
c
      end
c
c
c*module ormas1  *deck getqfcc
c     -----------------------------------------------------
      subroutine getqfcc(si1,si2,nact,nci,na,nb,
     *           iacon1,ibcon1,index,q,lbox1,lbox2,lbox3,
     *           nast,nbst,lsyma,lsymb,nsym,lspa,lspb,
     *           lgcom,lsas,lsbs,ktab,lsac,lsbc,ldisb,
     *           itga,itgb,iast,ibst)
c     -----------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      dimension si1(*),si2(*)
      dimension iacon1(na),ibcon1(na)
      dimension index((nact*(nact+1))/2+1),q(nci)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension nast(nsym+1),nbst(nsym+1),lsyma(iast),lsymb(ibst)
      dimension lspa(iast),lspb(ibst)
      dimension lgcom(itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension ktab(nsym)
      dimension lsac(iast),lsbc(ibst),ldisb(nsym,itgb,itga)
c
      do ii=1,nci
         q(ii)=0.0d+00
      enddo
c
c   big loop over alpha
c
      icount=0
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      do 1000 iga=1,itga
c
c   loop over all strings isa, in group iga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
         do 990 isa=nast(iga)+1,nast(iga+1)
            jasym=lsyma(isa)
            ksym=ktab(jasym)
c
c  make alpha part of diagonal energy.
c
            c=0.0d+00
            do ii=1,na
               i1=iacon1(ii)
               ind1=index(i1+1)
               c=c+si1(ind1)
               do jj=1,ii-1
                  i2=iacon1(jj)
                  ind2=index(i2+1)
                  indm=index(i1)+i2
                  j1=index(ind1)+ind2
                  j2=index(indm+1)
                  c=c+si2(j1)-si2(j2)
               enddo
            enddo
c
c  loop over relevant beta
c
         call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
         do 800 igb=1,itgb
            if (lgcom(igb,iga).ne.1) go to 750
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            istab=1
            do 700 isb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               icount=icount+1
c
               iendb=lsbc(isb)
               do kk=istab,iendb-1
                  call moveup2(lbox2,nspace,nb,msta,ibcon1)
               enddo
               istab=iendb
c
               d=0.0d+00
               do jj=1,nb
                  i2=ibcon1(jj)
                  ind2=index(i2+1)
                  do kk=1,na
                     i1=iacon1(kk)
                     ind1=index(i1+1)
                  ima = max(ind1,ind2)
                  imi = min(ind1,ind2)
                     j2=index(ima) + imi
                     d=d+si2(j2)
                  enddo
               enddo
               t = c+d
               q(icount)=q(icount)+t
c
  700       continue
  750       call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
  800    continue
c
            if (nb.eq.0) q(isa)=q(isa)+c
            call moveup2(lbox1,nspace,na,msta,iacon1)
  990    continue
         call pushco(lbox1,nspace,na,iama,iami,lbox3,iend)
 1000 continue
c
c  now big loop over beta
c
      call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
      do 2000 igb=1,itgb
c
c  loop over all strings isb, in group igb
c
         call resetde(lbox2,nspace,nb,msta,ibcon1)
         do 1990 isb=nbst(igb)+1,nbst(igb+1)
            jbsym=lsymb(isb)
            ksym=ktab(jbsym)
            ipb1=lspb(isb)
c
            c=0.0d+00
            do ii=1,nb
               i1=ibcon1(ii)
               ind1=index(i1+1)
               c=c+si1(ind1)
               do jj=1,ii-1
                  i2=ibcon1(jj)
                  ind2=index(i2+1)
                  indm=index(i1)+i2
                  j1=index(ind1)+ind2
                  j2=index(indm+1)
                  c=c+si2(j1)-si2(j2)
               enddo
            enddo
c
            do 1800 iga=1,itga
               if (lgcom(igb,iga).eq.0) go to 1800
               ipb2=ldisb(jbsym,igb,iga)+ipb1
               ipa2=nast(iga)
               do isa=lsas(ksym,iga),lsas(ksym+1,iga)-1
                  ipa=lsac(isa)+ipa2
                  ipo=lspa(ipa)
                  icit=ipo+ipb2
                  q(icit)=q(icit)+c
               enddo
 1800       continue
c
            call moveup2(lbox2,nspace,nb,msta,ibcon1)
 1990    continue
         call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 2000 continue
c
      return
      end
c
c*module ormas1   *deck fccisrt
c     ----------------------------------------------------
      subroutine fccisrt(ipica,ipicb,ci2,npic)
c     ----------------------------------------------------
      implicit integer(a-z)
      dimension ipica(npic),ipicb(npic)
      dimension ci2(npic)
c
c    code to sort fcci determinants.
c
c    ipica,ipicb contain a list of determinants, alpha and
c    beta strings desired for a full class ci.
c    ci2 contains the actual position of the determinants.
c
c    ipica,ipicb will contain the list of determinants ordered
c    first according to alpha and then beta.
c    ci2 will contain the actual position of the determinants
c    which will now be in ascending order.
c
      n = npic
      l=n/2+1
      ir=n
c
   10 continue
         if (l.gt.1) then
            l=l-1
            rra=ci2(l)
            irrb=ipica(l)
            irrc=ipicb(l)
         else
            rra=ci2(ir)
            irrb=ipica(ir)
            irrc=ipicb(ir)
            ci2(ir)=ci2(1)
            ipica(ir)=ipica(1)
            ipicb(ir)=ipicb(1)
c
            ir=ir-1
            if (ir.eq.1) then
               ci2(1)=rra
               ipica(1)=irrb
               ipicb(1)=irrc
c
               go to 122
            endif
         endif
         i=l
         j=l+l
   20    if (j.le.ir) then
            if (j.lt.ir) then
                   if (ci2(j).lt.ci2(j+1)) j=j+1
            endif
            if (rra.lt.ci2(j)) then
               ci2(i)=ci2(j)
               ipica(i)=ipica(j)
               ipicb(i)=ipicb(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         go to 20
         endif
         ci2(i)=rra
         ipica(i)=irrb
         ipicb(i)=irrc
      go to 10
c
  122 continue
c
      return
      end
c
c*module ormas1  *deck fccihe
c     --------------------------------------------------------
      subroutine fccihe(sint1,sint2,norb,naele,nbele,
     *       iacon1,ibcon1,iacon2,ibcon2,ij,ji,index,elem)
c     --------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension iacon1(naele),ibcon1(nbele)
      dimension iacon2(naele),ibcon2(nbele)
      dimension sint1(*),sint2(*)
      dimension index(norb*(norb+1)/2+1)
      integer diff1(2),diff2(2),ipos1(2),ipos2(2)
c
c     returns the matrix element < k | h | l > where k, l are
c     determinants.  alpha and beta occupations are stored in
c     iacon1, ibcon1 for k and iacon2, ibcon2 for l.
c
      elem = 0.0d+00
c
c    if determinants are same
c
      if (ij.eq.ji) then
         do 200 i=1,naele
            ia = iacon1(i)
            i1 = index(ia+1)
            elem = elem + sint1(i1)
            do 197 j=i+1,naele
               ia1 = iacon1(j)
               i2 = index(ia1+1)
               it = index(i2)+i1
               elem = elem + sint2(it)
               ic = index(ia1)+ia
               it = index(ic+1)
               elem = elem - sint2(it)
c
  197       continue
            do 198 j=1,nbele
               ib1 = ibcon1(j)
               i2 = index(ib1+1)
               ima = max(i1,i2)
               imi = min(i1,i2)
               it = index(ima) + imi
               elem = elem + sint2(it)
  198       continue
  200    continue
c
         do 210 i=1,nbele
            ib = ibcon1(i)
            i1 = index(ib+1)
            elem = elem + sint1(i1)
            do 204 j=i+1,nbele
               ib1 = ibcon1(j)
               i2 = index(ib1+1)
               it = index(i2)+i1
               elem = elem + sint2(it)
               ic = index(ib1)+ib
               it = index(ic+1)
               elem = elem - sint2(it)
  204       continue
  210    continue
      return
      endif
c
      idea=0
      ideb=0
c
c     different orbitals in first deterinant
c
      do 20 i=1,naele
         do 15 j=1,naele
            if (iacon1(i).eq.iacon2(j)) go to 20
   15    continue
         idea = idea + 1
         if (idea.gt.2) return
         diff1(idea) = iacon1(i)
         ipos1(idea) = i
   20 continue
c
      do 30 i=1,nbele
         do 25 j=1,nbele
            if (ibcon1(i).eq.ibcon2(j)) go to 30
   25    continue
         ideb = ideb + 1
         if (idea+ideb.gt.2) return
         diff1(idea+ideb) = ibcon1(i)
         ipos1(idea+ideb) = i
   30 continue
c
c
c    to find the different orbitals in second determinant
c
      ist = 1
      do 63 ii=1,idea
            do 50 i=ist,naele
               do 45 j=1,naele
                  if (iacon2(i).eq.iacon1(j)) go to 50
   45          continue
               go to 60
   50       continue
c
   60       diff2(ii) = iacon2(i)
            ipos2(ii) = i
            ist = i+1
   63 continue
c
      ist = 1
      do 163 ii=1,ideb
            do 150 i=ist,nbele
               do 145 j=1,nbele
                  if (ibcon2(i).eq.ibcon1(j)) go to 150
  145          continue
               go to 160
  150       continue
  160       diff2(ii+idea) = ibcon2(i)
            ipos2(ii+idea) = i
            ist = i+1
  163 continue
c
c    if determinants differ by one orbital
c
       if (idea+ideb.eq.1) then
c
c   one particle density contribution
c
          id1 = diff1(1)
          id2 = diff2(1)
          ida = max(id1,id2)
          idi = min(id1,id2)
          ind1 = index(ida) + idi
          iperm = (-1)**(ipos1(1)-ipos2(1))
          elem = elem + iperm*sint1(ind1)
c
c    two particle density contribution
c
c    if different orbitals are alpha spin orbs
c
         if (idea.eq.1) then
           do 673 k=1,naele
              nk = iacon1(k)
              if (nk.eq.id1) go to 673
              ind2 = index(nk+1)
              ima = max(ind1,ind2)
              imi = min(ind1,ind2)
              indx = index(ima) + imi
              elem = elem + sint2(indx)*iperm
              ima = max(id1,nk)
              imi = min(id1,nk)
              i1 = index(ima) + imi
              ima = max(id2,nk)
              imi = min(id2,nk)
              i2 = index(ima) + imi
              ima = max(i1,i2)
              imi = min(i1,i2)
              inx = index(ima) + imi
              elem = elem - iperm*sint2(inx)
  673     continue
c
           do 678 k=1,nbele
              nk = ibcon1(k)
              ind2 = index(nk+1)
              ima = max(ind1,ind2)
              imi = min(ind1,ind2)
              indx = index(ima) + imi
              elem = elem + iperm*sint2(indx)
  678     continue
c
        else
c
c     different orbitals are beta spin orbs
c
           do 732 k=1,naele
              nk = iacon1(k)
              ind2 = index(nk+1)
              ima = max(ind1,ind2)
              imi = min(ind1,ind2)
              indx = index(ima) + imi
              elem = elem + iperm*sint2(indx)
  732      continue
c
           do 752 k=1,nbele
              nk = ibcon1(k)
              if (nk.eq.id1) go to 752
              ind2 = index(nk+1)
              ima = max(ind1,ind2)
              imi = min(ind1,ind2)
              indx = index(ima) + imi
              elem = elem + iperm*sint2(indx)
              ima = max(id1,nk)
              imi = min(id1,nk)
              i1 = index(ima) + imi
              ima = max(id2,nk)
              imi = min(id2,nk)
              i2 = index(ima) + imi
              ima = max(i1,i2)
              imi = min(i1,i2)
              inx = index(ima) + imi
              elem = elem - iperm*sint2(inx)
  752      continue
c
          endif
          return
c
      else
c
c     two orbitals are different
c     contribution only to 2-particle density matrix.
c     differing orbitals in diff1(1),diff1(2) for con1 and
c     diff2(1),diff2(2) for con2.  position stored in
c     ipos1(1),ipos1(2) and ipos2(1),ipos2(2).
c
         iperm = (-1)**(ipos1(1)-ipos2(1)+ipos1(2)-ipos2(2))
         i11 = diff1(1)
         i12 = diff2(1)
         i21 = diff1(2)
         i22 = diff2(2)
         ima = max(i11,i12)
         imi = min(i11,i12)
         i1 = index(ima) + imi
         ima = max(i21,i22)
         imi = min(i21,i22)
         i2 = index(ima) + imi
         ima = max(i1,i2)
         imi = min(i1,i2)
         inx = index(ima) + imi
         elem = elem + iperm*sint2(inx)
c
c     if all differing orbitals are or same spin then
c     have extra matrix elements.
c
         if (idea.eq.2.or.ideb.eq.2) then
            ima = max(i11,i22)
            imi = min(i11,i22)
            i1 = index(ima) + imi
            ima = max(i12,i21)
            imi = min(i12,i21)
            i2 = index(ima) + imi
            ima = max(i1,i2)
            imi = min(i1,i2)
            inx = index(ima) + imi
            elem = elem - iperm*sint2(inx)
         endif
      endif
c
      return
      end
c
c*module ormas1  *deck initfcc
c     -----------------------------------------------------
      subroutine initfcc(iw,some,b,nci,na,nb,nact,iacon1,ibcon1,
     *           iacon2,ibcon2,isd,ido,ci,iwrk1,maxwx,kst,
     *           index,f,el,ef,si1,si2,iwrk2,imark,ab,
     *           lgcom,nast,nbst,lsyma,lsymb,
     *           lsbs,lsbc,lspa,lspb,ldisb,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           x,nx,landet,lbndet,
     *           itga,itgb,iast,ibst,ktab,nsym)
c     -----------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      logical some
      dimension si1(*),si2(*),iwrk2(maxwx),ab(nci,kst)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension iacon1(na),ibcon1(na),iacon2(na),ibcon2(na)
      dimension ci(nci,kst),iwrk1(maxwx,2)
      dimension index((nact*(nact+1))/2+1)
      dimension isd(na+nb),ido(na)
      dimension f((maxwx*(maxwx+1))/2),ef(maxwx,maxwx),el(maxwx)
      dimension lgcom(itgb,itga),nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst)
      dimension lsbs(nsym+1,itgb)
      dimension lsbc(ibst)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension ktab(nsym)
c
      maxw1=maxwx
      ibg=1
      imark=0
      if (nci.le.maxw1) then
         maxw1=nci
         imark=1
         if(some) write(iw,*)
     *        'small ci matrix, just using incore diagonalization...'
      else
         if(nb.eq.0) then
            write(iw,9020) nci,maxw1,nci
            call abrt
            stop
         end if
      endif
      if (kst.gt.nci) then
         if(some) write(iw,9010) kst,nci
         call abrt
         stop
      endif
c
c   sort of a loop structure here, keep coming
c   back to 999 until all initial determinant space <=maxw1
c   are found.
c
  999 continue
      if (ibg.gt.maxw1) go to 9999
      pmin = 100.0d+00
      ind = 0
      do ii=1,nci
         if (ci(ii,1).lt.pmin) then
            ind = ii
            pmin = ci(ind,1)
         endif
      enddo
c
c for nb=0
c
      if (nb.eq.0) then
         iwrk1(ibg,1) = ind
         iwrk1(ibg,2) = 1
         iwrk2(ibg) = ind
         ci(ind,1) = 101.0d+00
         ibg = ibg + 1
         go to 999
      endif
c
c  store the alpha and beta indices.
c
      ici=0
      do 100 iga=1,itga
         do 80 isa=nast(iga)+1,nast(iga+1)
            jasym=lsyma(isa)
            ksym=ktab(jasym)
            do 60 igb=1,itgb
               if (lgcom(igb,iga).ne.1) go to 60
               do 40 isb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ici=ici+1
                  if (ici.ne.ind) go to 40
                  nend=lsbc(isb)+nbst(igb)
                  go to 200
   40          continue
   60       continue
   80    continue
  100 continue
c
  200 continue
c
      isb = nend
c
c  make the alpha and beta parts.
c
      idea = isa - nast(iga)
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      do ii=1,iga-1
         call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
      enddo
      call resetde(lbox1,nspace,na,msta,iacon1)
      do ii=1,idea-1
         call moveup2(lbox1,nspace,na,msta,iacon1)
      enddo
c
      ideb = isb - nbst(igb)
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox3)
      do ii=1,igb-1
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox3,ise)
      enddo
      call resetde(lbox1,nspace,nb,msta,ibcon1)
      do ii=1,ideb-1
         call moveup2(lbox1,nspace,nb,msta,ibcon1)
      enddo
c
c  find doubly occupied orbitals, put in beginning of isd.
c
      nss = 0
      nsd = 0
      do ii=1,na
         ia = iacon1(ii)
         do jj=1,nb
            if (ia.eq.ibcon1(jj)) then
               nsd = nsd + 1
               isd(nsd) = ia
            endif
         enddo
      enddo
c
c   check to see if all beta orbitals are paired.
c
      if (nsd.eq.nb) then
         iwrk1(ibg,1) = isa
         iwrk1(ibg,2) = isb
         iwrk2(ibg) = ind
         ci(ind,1) = 101.0d+00
         ibg = ibg + 1
         if (ibg.le.maxw1) go to 999
         go to 9999
      endif
c
c      find singly occupied orbs now, put in end of isd, beta first
c      then alpha.
c
      do 20 ii=1,nb
         ib = ibcon1(ii)
         do 24 jj=1,nsd
            if (ib.eq.isd(jj)) go to 20
   24    continue
         nss = nss + 1
         isd(nss+nsd) = ib
   20 continue
c
      do 30 ii=1,na
         ia = iacon1(ii)
         do 34 jj=1,nsd
            if (ia.eq.isd(jj)) go to 30
   34    continue
         nss = nss + 1
         isd(nss+nsd) = ia
   30 continue
c
c       reorder the things.
c
      do ii=1,nss-1
         do jj=ii+1,nss
            if (isd(jj+nsd).lt.isd(ii+nsd)) then
               kk=isd(ii+nsd)
               isd(ii+nsd) = isd(jj+nsd)
               isd(jj+nsd) = kk
            endif
         enddo
      enddo
c
c  determine how many determinants there are with the
c  same space function.
c
      nspa=na-nsd
      ntop=1
      nbot=1
      do ii=nspa+1,nss
         ntop=ntop*ii
      enddo
      do ii=1,nss-nspa
         nbot=nbot*ii
      enddo
      node = ntop/nbot
c
      if (node+ibg-1.gt.maxw1) go to 9999
      do ii=1,nspa
         ido(ii) = ii
      enddo
c
c     now to store positions of all possible determinants with
c     same space function.  alpha first.
c
      do 3000 ijk=1,node
         do ii=1,nsd
            iacon1(ii) = isd(ii)
         enddo
         do ii=1,nspa
            iacon1(ii+nsd) = isd(nsd+ido(ii))
         enddo
c
c   must reorder here.
c
         do ii=1,na-1
            do jj=ii+1,na
               if (iacon1(jj).lt.iacon1(ii)) then
                  kk=iacon1(ii)
                  iacon1(ii) = iacon1(jj)
                  iacon1(jj) = kk
               endif
            enddo
         enddo
c
c  determine the position of the alpha string
c
         call getbox(lbox1,iacon1,na,msta,nspace)
         call positco(lbox2,nspace,na,iama,iami,lbox3,lbox1,iga)
         call idpost(iacon1,na,lbox1,nspace,msta,idim,x,nx,lbst,
     *        landet(1,iga),ibcon1,iposa)
c
         iwrk1(ijk+ibg-1,1) = iposa+nast(iga)
c
         call advanc(ido,nspa,nss)
 3000 continue
c
c  now for the beta
c
      nspb=nb-nsd
      do ii=1,nspb
         ido(ii)=ii
      enddo
c
      do 4000 ijk=1,node
         do ii=1,nsd
            ibcon1(ii) = isd(ii)
         enddo
         do ii=1,nspb
            ibcon1(ii+nsd) = isd(nsd+ido(ii))
         enddo
c
c  must reorder here
c
         do ii=1,nb-1
            do jj=ii+1,nb
               if (ibcon1(jj).lt.ibcon1(ii)) then
                  kk=ibcon1(ii)
                  ibcon1(ii) = ibcon1(jj)
                  ibcon1(jj) = kk
               endif
            enddo
         enddo
c
c  determine the position of the beta string.
c
         call getbox(lbox1,ibcon1,nb,msta,nspace)
         call positco(lbox2,nspace,nb,ibma,ibmi,lbox3,lbox1,igb)
         call idpost(ibcon1,nb,lbox1,nspace,msta,idim,x,nx,lbst,
     *        lbndet(1,igb),iacon1,iposb)
c
         iwrk1(ibg+node-ijk,2) = iposb+nbst(igb)
c
         call advanc(ido,nspb,nss)
 4000 continue
c
c  zero all diagonal elements just found.
c
      do ii=1,node
         ina=iwrk1(ii+ibg-1,1)
         do 110 iga=1,itga
            if (nast(iga+1).ge.ina) go to 120
  110    continue
  120    continue
         inb=iwrk1(ii+ibg-1,2)
         do 130 igb=1,itgb
            if (nbst(igb+1).ge.inb) go to 140
  130    continue
  140    continue
         isymb = lsymb(inb)
         ind = lspa(ina)+ldisb(isymb,igb,iga)+lspb(inb)
         iwrk2(ii+ibg-1) = ind
         ci(ind,1) = 101.0d+00
      enddo
c
      ibg = ibg + node
      go to 999
c
 9999 continue
      nsize = ibg - 1
c
c  now to reorder according to position of determinants.
c
      if (nsize.gt.1) call fccisrt(iwrk1,iwrk1(1,2),iwrk2,nsize)
c
c  now to form the hamiltonian.
c
      ixi=0
      igas=1
      inas=1
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      call resetde(lbox1,nspace,na,msta,iacon1)
      igbs=1
      inbs=1
      call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
      call resetde(lbox2,nspace,nb,msta,ibcon1)
c
      do 6000 ijk=1,nsize
c a
         inae=iwrk1(ijk,1)
         do ii=igas,itga
            if (nast(ii+1).ge.inae) go to 150
         enddo
  150    continue
         if (ii.gt.igas) then
            do jj=igas,ii-1
               call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
            enddo
            call resetde(lbox1,nspace,na,msta,iacon1)
            igas=ii
            inas=1
         endif
         inae = inae - nast(igas)
         do kk=inas,inae-1
            call moveup2(lbox1,nspace,na,msta,iacon1)
         enddo
         inas = inae
c b
         inbe=iwrk1(ijk,2)
         do ii=1,itgb
            if (nbst(ii+1).ge.inbe) go to 160
         enddo
  160    continue
         if (ii.lt.igbs) then
            call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
            do jj=1,ii-1
               call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            igbs=ii
            inbs=1
         elseif (ii.gt.igbs) then
            do jj=igbs,ii-1
               call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            igbs=ii
            inbs=1
         endif
         inbe = inbe - nbst(igbs)
         if (inbs.gt.inbe) then
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            inbs=1
         endif
         do kk=inbs,inbe-1
            call moveup2(lbox2,nspace,nb,msta,ibcon1)
         enddo
         inbs = inbe
c
      igas2=1
      inas2=1
      call resetco(lbox4,nspace,na,iama,iami,lbox3)
      call resetde(lbox4,nspace,na,msta,iacon2)
      igbs2=1
      inbs2=1
      call resetco(lbox5,nspace,nb,ibma,ibmi,lbox3)
      call resetde(lbox5,nspace,nb,msta,ibcon2)
c
      do 5000 kji=1,ijk
         ixi=ixi+1
c a
         inae2=iwrk1(kji,1)
         do ii=igas2,itga
            if (nast(ii+1).ge.inae2) go to 250
         enddo
  250    continue
         if (ii.gt.igas2) then
            do jj=igas2,ii-1
               call pushco(lbox4,nspace,na,iama,iami,lbox3,ise)
            enddo
            call resetde(lbox4,nspace,na,msta,iacon2)
            igas2=ii
            inas2=1
         endif
         inae2 = inae2 - nast(igas2)
         do kk=inas2,inae2-1
            call moveup2(lbox4,nspace,na,msta,iacon2)
         enddo
         inas2 = inae2
c b
         inbe2=iwrk1(kji,2)
         do ii=1,itgb
            if (nbst(ii+1).ge.inbe2) go to 260
         enddo
  260    continue
         if (ii.lt.igbs2) then
            call resetco(lbox5,nspace,nb,ibma,ibmi,lbox3)
            do jj=1,ii-1
               call pushco(lbox5,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            igbs2=ii
            inbs2=1
         elseif (ii.gt.igbs2) then
            do jj=igbs2,ii-1
               call pushco(lbox5,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            igbs2=ii
            inbs2=1
         endif
         inbe2 = inbe2 - nbst(igbs2)
         if (inbs2.gt.inbe2) then
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            inbs2=1
         endif
         do kk=inbs2,inbe2-1
            call moveup2(lbox5,nspace,nb,msta,ibcon2)
         enddo
         inbs2=inbe2
c
         call fccihe(si1,si2,nact,na,nb,iacon1,ibcon1,
     *             iacon2,ibcon2,ijk,kji,index,elem)
c
         f(ixi) = elem
c
 5000 continue
 6000 continue
c
      call evvrsp(-1,nsize,nsize,(nsize*(nsize+1))/2,maxwx,
     *              f,b,iwrk1,el,ef,0,ierr)
      if (ierr.ne.0) then
         if(some) write(iw,*) 'error in small diagonalization'
         if(some) write(iw,*) 'ierr = ',ierr
         return
      endif
c
      do 347 ii=1,nci
         do 450 jj=1,kst
         ci(ii,jj) = 0.0d+00
         ab(ii,jj) = 0.0d+00
  450    continue
  347 continue
c
      do 799 ijk=1,kst
         do 899 ii=1,nsize
            ki = iwrk2(ii)
            ci(ki,ijk) = ef(ii,ijk)
  899    continue
  799 continue
c
      return
c
 9010 format(/1x,'***** error *****'/
     *       1x,'input nstate=',i4,' exceeds hamiltonian dimension',i5)
 9020 format(/1x,'***** error *****'/
     *   1x,'this job has no beta electrons, and more determinants=',i8/
     *   1x,'than the initial hamiltonian matrix guess size=',i8,'.'/
     *   1x,'please increase -nhgss- in $det to ',i8,' and rerun.')
c
      end
c
c*module ormas1  *deck fccspin
c     ----------------------------------------------------------
      subroutine fccspin(na,nb,iacon1,iacon2,ibcon1,ibcon2,
     *           isd,ido,ci,ab,nv,nci,spin,iwrk1,iwrk2,maxw1,
     *           lgcom,nast,nbst,lsyma,lsymb,
     *           lsbs,lsbc,lspa,lspb,ldisb,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           x,nx,landet,lbndet,
     *           itga,itgb,iast,ibst,ktab,nsym)
c     ----------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      dimension ci(nci,nv),ab(nci,nv)
      dimension iacon1(na),ibcon1(nb)
      dimension iacon2(na),ibcon2(nb)
      dimension isd(na+nb),ido(na)
      dimension iwrk1(maxw1,2),iwrk2(maxw1)
      dimension spin(nv)
      dimension x(nx),landet(nspace,itga),lbndet(nspace,itgb)
      dimension lgcom(itgb,itga),nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst)
      dimension lsbs(nsym+1,itgb)
      dimension lsbc(ibst)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension ktab(nsym)
c
c     returns spin of ci vector i in spin(i) for full space
c     determinantal wavefunction.
c     nbig is the number of highest ci coefficients to take.
c
      crit = 1.0d-10
      nbig=min(5,nci)
c
      do ii=1,nv
         spin(ii) = 0.0d+00
      enddo
c
      do ii=1,nv
         do jj=1,nci
            ab(jj,ii) = ci(jj,ii)
         enddo
      enddo
c
c  loop over states
c
      do 9000 jjj=1,nv
         znorm = 0.0d+0
c
c  loop over nbig
c
         do 8000 iii=1,nbig
            ndone=0
c
c  find the largest ci coefficient.
c
            ind = 1
            pmin = 0.0d+00
            do ii=1,nci
               if (abs(ab(ii,jjj)).gt.pmin) then
                  pmin = abs(ab(ii,jjj))
                  ind = ii
               endif
            enddo
            ab(ind,jjj) = 0.0d+00
            if (pmin.lt.crit) go to 8050
c
c  store the alpha and beta indices.
c
      ici=0
      do 100 iga=1,itga
         do 80 isa=nast(iga)+1,nast(iga+1)
            jasym=lsyma(isa)
            ksym=ktab(jasym)
            do 62 igb=1,itgb
               if (lgcom(igb,iga).ne.1) go to 62
               do 40 isb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ici=ici+1
                  if (ici.ne.ind) go to 40
                  nend=lsbc(isb)+nbst(igb)
                  go to 200
   40          continue
   62       continue
   80    continue
  100 continue
c
  200 continue
c
      isb = nend
c
c  make the alpha and beta parts.
c
      idea = isa - nast(iga)
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      do ii=1,iga-1
         call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
      enddo
      call resetde(lbox1,nspace,na,msta,iacon1)
      do ii=1,idea-1
         call moveup2(lbox1,nspace,na,msta,iacon1)
      enddo
c
      ideb = isb - nbst(igb)
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox3)
      do ii=1,igb-1
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox3,ise)
      enddo
      call resetde(lbox1,nspace,nb,msta,ibcon1)
      do ii=1,ideb-1
         call moveup2(lbox1,nspace,nb,msta,ibcon1)
      enddo
c
c
c  find doubly occupied orbitals, put in beginning of isd.
c
      nss = 0
      nsd = 0
      do ii=1,na
         ia = iacon1(ii)
         do jj=1,nb
            if (ia.eq.ibcon1(jj)) then
               nsd = nsd + 1
               isd(nsd) = ia
            endif
         enddo
      enddo
c
c  nsd is how many doubly occupied orbitals there are.
c  isd contains a list of the doubly occuppied orbitals
c
c      find singly occupied orbs now, put in end of isd, beta first
c      then alpha.
c
            do 320 ii=1,nb
                  ib = ibcon1(ii)
               do 324 jj=1,nsd
                  if (ib.eq.isd(jj)) go to 320
  324          continue
               nss = nss + 1
               isd(nss+nsd) = ib
  320       continue
c
             do 330 ii=1,na
                 ia = iacon1(ii)
                 do 334 jj=1,nsd
                  if (ia.eq.isd(jj)) go to 330
  334          continue
                nss = nss + 1
                 isd(nss+nsd) = ia
  330       continue
c
c       reorder the things.
c
            do 340 ii=1,nss-1
                 do 342 jj=ii+1,nss
               if (isd(jj+nsd).lt.isd(ii+nsd)) then
                  kk=isd(ii+nsd)
                  isd(ii+nsd) = isd(jj+nsd)
                  isd(jj+nsd) = kk
               endif
  342          continue
  340       continue
c
c  determine how many determinants there are with the
c  same space function.
c
      nspa=na-nsd
      ntop=1
      nbot=1
      do ii=nspa+1,nss
         ntop=ntop*ii
      enddo
      do ii=1,nss-nspa
         nbot=nbot*ii
      enddo
      node = ntop/nbot
c
      do ii=1,nspa
         ido(ii) = ii
      enddo
c
c     now to store positions of all possible determinants with
c     same space function.  alpha first.
c
      do 3000 ijk=1,node
         do ii=1,nsd
            iacon1(ii) = isd(ii)
         enddo
         do ii=1,nspa
            iacon1(ii+nsd) = isd(nsd+ido(ii))
         enddo
c
c   must reorder here.
c
         do ii=1,na-1
            do jj=ii+1,na
               if (iacon1(jj).lt.iacon1(ii)) then
                  kk=iacon1(ii)
                  iacon1(ii) = iacon1(jj)
                  iacon1(jj) = kk
               endif
            enddo
         enddo
c
c  determine the position of the alpha string
c
         call getbox(lbox1,iacon1,na,msta,nspace)
         call positco(lbox2,nspace,na,iama,iami,lbox3,lbox1,iga)
         call idpost(iacon1,na,lbox1,nspace,msta,idim,x,nx,lbst,
     *        landet(1,iga),ibcon1,iposa)
c
         iwrk1(ijk,1) = iposa+nast(iga)
c
         call advanc(ido,nspa,nss)
 3000 continue
c
c  now for the beta
c
      nspb=nb-nsd
      do ii=1,nspb
         ido(ii)=ii
      enddo
c
      do 4000 ijk=1,node
         do ii=1,nsd
            ibcon1(ii) = isd(ii)
         enddo
         do ii=1,nspb
            ibcon1(ii+nsd) = isd(nsd+ido(ii))
         enddo
c
c  must reorder here
c
         do ii=1,nb-1
            do jj=ii+1,nb
               if (ibcon1(jj).lt.ibcon1(ii)) then
                  kk=ibcon1(ii)
                  ibcon1(ii) = ibcon1(jj)
                  ibcon1(jj) = kk
               endif
            enddo
         enddo
c
c  determine the position of the beta string.
c
         call getbox(lbox1,ibcon1,nb,msta,nspace)
         call positco(lbox2,nspace,nb,ibma,ibmi,lbox3,lbox1,igb)
         call idpost(ibcon1,nb,lbox1,nspace,msta,idim,x,nx,lbst,
     *        lbndet(1,igb),iacon1,iposb)
c
         iwrk1(node-ijk+1,2) = iposb+nbst(igb)
c
         call advanc(ido,nspb,nss)
 4000 continue
c
c  zero all diagonal elements just found.
c
      do ii=1,node
         ina=iwrk1(ii,1)
         do 110 iga=1,itga
            if (nast(iga+1).ge.ina) go to 120
  110    continue
  120    continue
         inb=iwrk1(ii,2)
         do 130 igb=1,itgb
            if (nbst(igb+1).ge.inb) go to 140
  130    continue
  140    continue
         isymb = lsymb(inb)
         ind = lspa(ina)+ldisb(isymb,igb,iga)+lspb(inb)
         iwrk2(ii) = ind
         ab(ind,jjj) = 0.0d+00
      enddo
c
c  now to reorder according to position of determinants.
c
      if (node.gt.1) call fccisrt(iwrk1,iwrk1(1,2),iwrk2,node)
c
c  now to include these determinants contribution to the spin.
c
      igas=1
      inas=1
      call resetco(lbox1,nspace,na,iama,iami,lbox3)
      call resetde(lbox1,nspace,na,msta,iacon1)
      igbs=1
      inbs=1
      call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
      call resetde(lbox2,nspace,nb,msta,ibcon1)
c
      do 6000 ijk=1,node
         znorm = znorm + ci(iwrk2(ijk),jjj)**2.0d+00
c a
         inae=iwrk1(ijk,1)
         do ii=igas,itga
            if (nast(ii+1).ge.inae) go to 150
         enddo
  150    continue
         if (ii.gt.igas) then
            do jj=igas,ii-1
               call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
            enddo
            call resetde(lbox1,nspace,na,msta,iacon1)
            igas=ii
            inas=1
         endif
         inae = inae - nast(igas)
         do kk=inas,inae-1
            call moveup2(lbox1,nspace,na,msta,iacon1)
         enddo
         inas = inae
c b
         inbe=iwrk1(ijk,2)
         do ii=1,itgb
            if (nbst(ii+1).ge.inbe) go to 160
         enddo
  160    continue
         if (ii.lt.igbs) then
            call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
            do jj=1,ii-1
               call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            igbs=ii
            inbs=1
         elseif (ii.gt.igbs) then
            do jj=igbs,ii-1
               call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            igbs=ii
            inbs=1
         endif
         inbe = inbe - nbst(igbs)
         if (inbs.gt.inbe) then
            call resetde(lbox2,nspace,nb,msta,ibcon1)
            inbs=1
         endif
         do kk=inbs,inbe-1
            call moveup2(lbox2,nspace,nb,msta,ibcon1)
         enddo
         inbs = inbe
c
      igas2=1
      inas2=1
      call resetco(lbox4,nspace,na,iama,iami,lbox3)
      call resetde(lbox4,nspace,na,msta,iacon2)
      igbs2=1
      inbs2=1
      call resetco(lbox5,nspace,nb,ibma,ibmi,lbox3)
      call resetde(lbox5,nspace,nb,msta,ibcon2)
c
c  loop over left side of the operator,ie determinants i=1,node
c
      do 5000 kji=1,node
         if (kji.eq.ijk) go to 5000
c a
         inae2=iwrk1(kji,1)
         do ii=igas2,itga
            if (nast(ii+1).ge.inae2) go to 250
         enddo
  250    continue
         if (ii.gt.igas2) then
            do jj=igas2,ii-1
               call pushco(lbox4,nspace,na,iama,iami,lbox3,ise)
            enddo
            call resetde(lbox4,nspace,na,msta,iacon2)
            igas2=ii
            inas2=1
         endif
         inae2 = inae2 - nast(igas2)
         do kk=inas2,inae2-1
            call moveup2(lbox4,nspace,na,msta,iacon2)
         enddo
         inas2 = inae2
c b
         inbe2=iwrk1(kji,2)
         do ii=1,itgb
            if (nbst(ii+1).ge.inbe2) go to 260
         enddo
  260    continue
         if (ii.lt.igbs2) then
            call resetco(lbox5,nspace,nb,ibma,ibmi,lbox3)
            do jj=1,ii-1
               call pushco(lbox5,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            igbs2=ii
            inbs2=1
         elseif (ii.gt.igbs2) then
            do jj=igbs2,ii-1
               call pushco(lbox5,nspace,nb,ibma,ibmi,lbox3,ise)
            enddo
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            igbs2=ii
            inbs2=1
         endif
         inbe2 = inbe2 - nbst(igbs2)
         if (inbs2.gt.inbe2) then
            call resetde(lbox5,nspace,nb,msta,ibcon2)
            inbs2=1
         endif
         do kk=inbs2,inbe2-1
            call moveup2(lbox5,nspace,nb,msta,ibcon2)
         enddo
         inbs2=inbe2
c
c  now have the two determinants.  have to find out whether
c  they differ by only two orbitals, one in alpha and one in
c  beta.  if differ by more than no contribution to spin.
c
         id1=0
         do 167 ix=1,na
            ia1=iacon1(ix)
            do 169 iy=1,na
               ia2=iacon2(iy)
               if (ia1.eq.ia2) go to 167
  169       continue
            id1=id1+1
            if (id1.gt.1) go to 5000
            ja1 = ix
  167    continue
         if (id1.ne.1) go to 5000
c
         id1=0
         do 171 ix=1,na
            ia1=iacon2(ix)
            do 172 iy=1,na
               ia2=iacon1(iy)
               if (ia1.eq.ia2) go to 171
  172       continue
            id1=id1+1
            if(id1.gt.1) go to 5000
            ja2 = ix
  171    continue
         if (id1.ne.1) go to 5000
c
         id1=0
         do 175 ix=1,nb
            ib1=ibcon1(ix)
            do 177 iy=1,nb
               ib2=ibcon2(iy)
               if (ib1.eq.ib2) go to 175
  177       continue
            id1=id1+1
            if (id1.gt.1) go to 5000
            jb1 = ix
  175    continue
         if (id1.ne.1) go to 5000
c
         id1=0
         do 178 ix=1,nb
            ib1=ibcon2(ix)
            do 179 iy=1,nb
               ib2=ibcon1(iy)
               if (ib1.eq.ib2) go to 178
  179       continue
            id1=id1+1
            if (id1.gt.1) go to 5000
            jb2=ix
  178    continue
         if (id1.ne.1) go to 5000
c
         if (iacon1(ja1).ne.ibcon2(jb2)) go to 5000
         if (iacon2(ja2).ne.ibcon1(jb1)) go to 5000
c
         ipt = abs(ja1-ja2) + abs(jb1-jb2) + 1
         iper = (-1)**ipt
c
      spin(jjj)=spin(jjj)+
     *      ci(iwrk2(ijk),jjj)*ci(iwrk2(kji),jjj)*iper
c
 5000 continue
c
      spin(jjj)=spin(jjj)-
     *  (ci(iwrk2(ijk),jjj)*ci(iwrk2(ijk),jjj))*nsd
c
 6000 continue
c
      ndone = ndone + node
      if (ndone.ge.nci) go to 8050
 8000 continue
c
 8050 continue
      spin(jjj)=spin(jjj)/znorm
c
 9000 continue
c
      do ii=1,nv
         spin(ii)=spin(ii)+((na-nb)/2.0d+00)**2.0d+00 + (na+nb)/2.0d+00
      enddo
      do ii=1,nv
         srt = sqrt(4.0d+00*abs(spin(ii)) + 1.0d+00)
         spin(ii) = (srt-1.0d+00)/2.0d+00
      enddo
c
      return
      end
c
c*module ormas1  *deck fccsup
c     ------------------------------------------------------------------
      subroutine fccsup(iw,na,nb,nact,iacon1,ibcon1,ibcon2,
     *           index,nbst,lspb,
     *           lbox1,lbox2,lbox3,lbox4,
     *           x,nx,lbndet,lsymb,nsym,itgb,ibst,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      dimension iacon1(na),ibcon1(nb),ibcon2(nb)
      dimension index((nact*(nact+1))/2+1)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace)
      dimension nbst(itgb+1),lspb(ibst)
      dimension x(nx)
      dimension lbndet(nspace,itgb),lsymb(ibst)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(nsym+1,ibst+1),jb1sy(nb*(nact-nb))
c
c  make and store all b -> b' data, where b' > b.
c
      nb1ch = 0
      inb = 0
c
c  loop through all beta groups
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox3)
c
      do 1000 iib = 1,itgb
c
         call resetde(lbox1,nspace,nb,msta,ibcon1)
c
         do 900 kkb=nbst(iib)+1,nbst(iib+1)
            inb = inb + 1
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
            jb1st(1,inb) = nb1ch+1
            kbst = nb1ch+1
c
c  loop over all single excitations, checking to see if it is valid.
c
c  loop over spaces to excite electron from.
c
            iebs = nb+1
            do 890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) go to 890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
               do 885 ib1 = iebe,iebs,-1
                  io1 = ibcon1(ib1)
                  igbe = iebe - lbox1(ispb1)
c
c  loop over possible spaces to excite into.
c
               do 880 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) go to 870
                  if (lbox2(ispb2).gt.ibma(ispb2)) go to 870
c
c  get group number
c
          call positco(lbox4,nspace,nb,ibma,ibmi,lbox3,lbox2,igb)
          nias = nbst(igb)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 860 igap=igba,igbe+1
c
                     do 850 jj=ista,iend
c
                        nb1ch = nb1ch + 1
                        jb1gr(nb1ch) = igb
                        jb1in(nb1ch) = index(jj)+io1
c
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb),iacon1,iposb)
c
                        jb1pe(nb1ch) = (-1)**iper
                        jb1po(nb1ch) = lspb(iposb+nias)
                        jb1sy(nb1ch-kbst+1) = lsymb(iposb+nias)
c
  850                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
  860             continue
c
  870             lbox2(ispb2) = lbox2(ispb2) - 1
  880          continue
c
  885          continue
c
               lbox2(ispb1) = lbox2(ispb1) + 1
  890       continue
c
c  order these excitations for beta inb according to symmetry.
c
            knum = nb1ch - kbst + 1
            call fccsrt2(jb1gr(kbst),jb1pe(kbst),jb1in(kbst),
     *                   jb1po(kbst),jb1sy,knum)
c
c  make the symmetry starting points
c ----
            ist=1
            do 100 ii=1,nsym
               do 200 jj=ist,knum
                  if (jb1sy(jj).ne.ii) go to 180
  200          continue
  180          jb1st(ii+1,inb) = jj+kbst-1
               ist = jj
  100       continue
c ----
c
c  reorder within each symmetry by position.
c
            ijst=jb1st(1,inb)
            do 300 ii=1,nsym
               jsta=jb1st(ii,inb)
               jend=jb1st(ii+1,inb)
               jnum=jend-jsta
               if (jnum.le.1) go to 300
               call fccsrt2(jb1gr(jsta),jb1pe(jsta),jb1in(jsta),
     *                      jb1sy(jsta-ijst+1),jb1po(jsta),jnum)
  300       continue
c
            call moveup2(lbox1,nspace,nb,msta,ibcon1)
  900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox3,iend)
c
 1000 continue
c
      if (nb1ch.ne.nb1ex) then
         write(iw,9000)
         call abrt
         stop
      endif
      jb1st(1,ibst+1) = nb1ch+1
c
      return
c
 9000 format(/1x,'error in calculation of single beta excites !!! ',
     *           'panic !!!!!!'/1x,'call joe and tell him off.')
c
      end
c
c*module ormas1   *deck fccsrt2
c     ------------------------------------------------------
      subroutine fccsrt2(jb1gr,jb1pe,jb1in,jb1po,jb1sy,n)
c     ------------------------------------------------------
      implicit integer(a-z)
      dimension jb1gr(n),jb1pe(n),jb1in(n),jb1po(n),jb1sy(n)
c
c    sorting garbage.
c
      if (n.lt.2) return
      l=n/2+1
      ir=n
c
   10 continue
         if (l.gt.1) then
            l=l-1
            irra=jb1sy(l)
            irrb=jb1gr(l)
            irrc=jb1pe(l)
            irrd=jb1in(l)
            irre=jb1po(l)
         else
            irra=jb1sy(ir)
            irrb=jb1gr(ir)
            irrc=jb1pe(ir)
            irrd=jb1in(ir)
            irre=jb1po(ir)
            jb1sy(ir)=jb1sy(1)
            jb1gr(ir)=jb1gr(1)
            jb1pe(ir)=jb1pe(1)
            jb1in(ir)=jb1in(1)
            jb1po(ir)=jb1po(1)
c
            ir=ir-1
            if (ir.eq.1) then
               jb1sy(1)=irra
               jb1gr(1)=irrb
               jb1pe(1)=irrc
               jb1in(1)=irrd
               jb1po(1)=irre
c
               go to 122
            endif
         endif
         i=l
         j=l+l
   20    if (j.le.ir) then
            if (j.lt.ir) then
                   if (jb1sy(j).lt.jb1sy(j+1)) j=j+1
            endif
            if (irra.lt.jb1sy(j)) then
               jb1sy(i)=jb1sy(j)
               jb1gr(i)=jb1gr(j)
               jb1pe(i)=jb1pe(j)
               jb1in(i)=jb1in(j)
               jb1po(i)=jb1po(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         go to 20
         endif
         jb1sy(i)=irra
         jb1gr(i)=irrb
         jb1pe(i)=irrc
         jb1in(i)=irrd
         jb1po(i)=irre
      go to 10
c
  122 continue
c
      return
      end
c
c*module ormas1   *deck fccsrt3
c     ------------------------------------------------------
      subroutine fccsrt3(jb1gr,jb1pe,jb1in,jb1sy,n)
c     ------------------------------------------------------
      implicit integer(a-z)
      dimension jb1gr(n),jb1pe(n),jb1in(n),jb1sy(n)
c
c    sorting garbage.
c
      if (n.lt.2) return
      l=n/2+1
      ir=n
c
   10 continue
         if (l.gt.1) then
            l=l-1
            irra=jb1sy(l)
            irrb=jb1gr(l)
            irrc=jb1pe(l)
            irrd=jb1in(l)
         else
            irra=jb1sy(ir)
            irrb=jb1gr(ir)
            irrc=jb1pe(ir)
            irrd=jb1in(ir)
            jb1sy(ir)=jb1sy(1)
            jb1gr(ir)=jb1gr(1)
            jb1pe(ir)=jb1pe(1)
            jb1in(ir)=jb1in(1)
c
            ir=ir-1
            if (ir.eq.1) then
               jb1sy(1)=irra
               jb1gr(1)=irrb
               jb1pe(1)=irrc
               jb1in(1)=irrd
c
               go to 122
            endif
         endif
         i=l
         j=l+l
   20    if (j.le.ir) then
            if (j.lt.ir) then
                   if (jb1sy(j).lt.jb1sy(j+1)) j=j+1
            endif
            if (irra.lt.jb1sy(j)) then
               jb1sy(i)=jb1sy(j)
               jb1gr(i)=jb1gr(j)
               jb1pe(i)=jb1pe(j)
               jb1in(i)=jb1in(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         go to 20
         endif
         jb1sy(i)=irra
         jb1gr(i)=irrb
         jb1pe(i)=irrc
         jb1in(i)=irrd
      go to 10
c
  122 continue
c
      return
      end
c
c*module ormas1  *deck dafcci
c     ------------------------------------------------------------------
      subroutine dafcci(iw,some,econst,istat,
     *           si1,si2,m2,m4,nact,nci,na,nb,k,kst,maxp,maxw1,
     *           niter,crit,
     *           x,nx,
     *           spin,el,nsym,iob,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           ci,ab,q,b,ef,f,ec,gr,
     *           iposa,ipera,iind1,igroa,iwrk1,iwrk2,
     *           isd,ido,index,immc,ihmcon,
     *           lgmul,ktab,iacon1,iacon2,ibcon1,ibcon2,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           nb1ex,ibst2,
     *           jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical some,goparr,dskwrk,maswrk,dsksav,fdirct,qcorr
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension si1(m2),si2(m4),x(nx),spin(kst),el(maxw1)
      dimension iob(nact),lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension ci(nci,maxp),ab(nci,maxp),q(nci),b(8*maxw1)
      dimension ef(maxw1,maxw1),f((maxw1*(maxw1+1))/2)
      dimension ec(maxp,maxp),gr(maxp)
      dimension iposa(na*(nact-na)*nsym)
      dimension ipera(na*(nact-na)*nsym)
      dimension iind1(na*(nact-na)*nsym)
      dimension igroa(na*(nact-na)*nsym)
      dimension iwrk1(2*maxw1),iwrk2(maxw1),isd(na+nb),ido(na)
      dimension index((nact*(nact+1))/2+1)
      dimension immc(nsym),ihmcon(k)
      dimension lgmul(nsym,nsym),ktab(nsym)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast),lsbc(ibst)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1st(ibst2)
      dimension jb1po(nb1ex),jb1sy(nb*(nact-nb))
c
c     flag for 1st full-NR iter
c
      logical first_fnr
      common /flag_fnr/ first_fnr
c 
      if (.not.fdirct) then
         idim1 = nsym+1
         idim2 = ibst+1
      else
         idim1 = 1
         idim2 = 1
      endif
c
      if (some) call tsecnd(e0)
      istat=0
c
      do 7 i=1,(nact*(nact+1))/2 + 1
         index(i) = (i*(i-1))/2
    7 continue
c
c     keep ci vectors from previous fullnr iteration 
c     can be the 'guess' vector for this calc.
c     provided they are not zeroed here and were 
c     not deallocated at the end of the last run.
c
c
      do 20 ii=1,kst
         do 30 jj=1,nci
            ab(jj,ii) = 0.0d+00
   30    continue
   20 continue
c
      if (first_fnr) then
c
        do ii=1,kst
          do jj=1,nci
            ci(jj,ii) = 0.0d+00
          end do
        end do
c
      end if
c
c     initital setup, work out diagonal elements.
c
      call getqfcc(si1,si2,nact,nci,na,nb,iacon1,ibcon1,
     *     index,q,lbox1,lbox2,lbox3,
     *     nast,nbst,lsyma,lsymb,nsym,lspa,lspb,
     *     lgcom,lsas,lsbs,ktab,lsac,lsbc,ldisb,
     *     itga,itgb,iast,ibst)
c
      critsv = crit
      nvecs = 0
c
c     determine the initial guess vectors.
c
c     skip to 3333 if not first fullnr iteration since we have 
c     guess vectors from the previous fullnr calc.
c
      if (.not.first_fnr) go to 3333
c
      do ii=1,nci
         ci(ii,1) = q(ii)
      enddo
c
      if(some) call tsecnd(e0)
c
      call initfcc(iw,some,b,nci,na,nb,nact,iacon1,ibcon1,
     *     iacon2,ibcon2,isd,ido,ci,iwrk1,maxw1,kst,
     *     index,f,el,ef,si1,si2,iwrk2,imark,ab,
     *     lgcom,nast,nbst,lsyma,lsymb,
     *     lsbs,lsbc,lspa,lspb,ldisb,
     *     lbox1,lbox2,lbox3,lbox4,lbox5,
     *     x,nx,landet,lbndet,
     *     itga,itgb,iast,ibst,
     *     ktab,nsym)
c
      if(some) then
         call tsecnd(e1)
         elap = e1 - e0
         e0 = e1
         if (nvecs.eq.0) write(iw,9010) elap
      end if
c
c     check if we have finished the ci by doing the first
c     diagonalization.
c
      if (imark.eq.1) then
         call  fccspin(na,nb,iacon1,iacon2,ibcon1,ibcon2,
     *      isd,ido,ci,ab,k,nci,spin,iwrk1,iwrk2,maxw1,
     *           lgcom,nast,nbst,lsyma,lsymb,
     *           lsbs,lsbc,lspa,lspb,ldisb,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           x,nx,landet,lbndet,
     *           itga,itgb,iast,ibst,ktab,nsym)
c
         if (.not.fdirct) then
         call fccsup(iw,na,nb,nact,iacon1,ibcon1,ibcon2,
     *           index,nbst,lspb,
     *           lbox1,lbox2,lbox3,lbox4,
     *           x,nx,lbndet,lsymb,nsym,
     *           itgb,ibst,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy)
         endif
         return
      endif
c
 3333 continue
c
c
c     if not fully direct, call setup for storage of b -> b' data.
c     this is needed for density matrix determination if not here
c     in the ci
c
      if (.not.fdirct) then
         call fccsup(iw,na,nb,nact,iacon1,ibcon1,ibcon2,
     *           index,nbst,lspb,
     *           lbox1,lbox2,lbox3,lbox4,
     *           x,nx,lbndet,lsymb,nsym,
     *           itgb,ibst,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,jb1sy)
c
         if(some) then
            call tsecnd(e1)
            elap = e1 - e0
            e0 = e1
            write(iw,9015) elap
         endif
      endif
c
      if (na.eq.nb) then
         call  fccspin(na,nb,iacon1,iacon2,ibcon1,ibcon2,
     *      isd,ido,ci,ab,kst,nci,spin,iwrk1,iwrk2,maxw1,
     *           lgcom,nast,nbst,lsyma,lsymb,
     *           lsbs,lsbc,lspa,lspb,ldisb,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           x,nx,landet,lbndet,
     *           itga,itgb,iast,ibst,ktab,nsym)
c
         do i=1,kst
            iwrk1(i)=int(spin(i) + 0.3d+00)
            iwrk2(i)=i
         enddo
c
         if(some) then
            call tsecnd(e1)
            elap = e1 - e0
            e0 = e1
            write(iw,9013) elap
         endif
      endif
c
c     if ms=0, impose restriction on the initial ci coefficients.
c
      if (na.eq.nb) then
      ip = kst
      ntcon = k
      do 5511 iga=1,itga
         do 5522 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
c
            do 5533 igb=1,iga
            if (lgcom(igb,iga).ne.1) goto 5533
            icc1 = jpza1 + ldisb(ksym,igb,iga)
            icc2 = ldisb(jasym,iga,igb) + lspb(kka)
            nibs=nbst(igb)
            do 5544 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               ibpos = lsbc(kkb) + nibs
               if (ibpos.gt.kka) goto 5522
               ici1 = icc1+lspb(ibpos)
               if (ibpos.eq.kka) then
                  do 5532 kj=1,ntcon
                  ips = (iwrk1(kj)/2)
                  if ((ips+ips).ne.iwrk1(kj)) ci(ici1,kj+ip) = 0.0d+00
 5532             continue
                  goto 5522
               endif
               ici2 = lspa(ibpos) + icc2
               do 5531 kj=1,ntcon
               is = (-1)**iwrk1(kj)
               ci(ici2,kj+ip) = is*ci(ici1,kj+ip)
 5531          continue
 5544       continue
 5533       continue
 5522    continue
 5511 continue
      end if
c
c  the routine names have meaning.  fchc is common to all of them.
c  the 5th character can be a 0 or x.  0 means ms=0, x means otherwise.
c  the 6th character can be a 1 or y.  1 means a single hc vector is
c  determined, y means more than one is determined.
c  the 7th and final character can be a d or s.  d means fully direct
c  and s means semi-direct.  d is yet to be implemented, i will
c  see if there is a need for it first.
c
      if (na.eq.nb) then
      if (kst.gt.1) then
      call fchc0ys(si1,si2,index,nact,na,nb,ci,ab,kst,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,
     *           iwrk1,iwrk2)
      else
      call fchc01s(si1,si2,index,nact,na,nb,ci,ab,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,idim1,idim2,
     *           iwrk1,iwrk2)
      endif
      else
      if (kst.gt.1) then
      call fchcxys(si1,si2,index,nact,na,nb,ci,ab,kst,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
      else
      call fchcx1s(si1,si2,index,nact,na,nb,ci,ab,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
      endif
      endif
c
      if(some) then
         call tsecnd(e1)
         elap = e1 - e0
         e0 = e1
         write(iw,9020) elap
      endif
c
      do 13 ii=1,kst
         el(ii) = 0.0d+00
         do 15 kk=1,nci
            el(ii) = el(ii) + ci(kk,ii)*ab(kk,ii)
   15    continue
   13 continue
      do 555 ii=1,maxp
         do 677 jj=1,ii-1
            ec(ii,jj) = 0.0d+00
            ef(ii,jj) = 0.0d+00
            ec(jj,ii) = 0.0d+00
            ef(jj,ii) = 0.0d+00
  677    continue
         ec(ii,ii) = 1.0d+00
         ef(ii,ii) = 1.0d+00
  555 continue
c
c     now to get into the loop, set some loop constants here
c     ip is the current number of ci vectors being dealt with.
c     il is the current root being optimized.
c
      ipxt = -1
      ip = kst
      il = 1
      ntcon = k
      do 4599 kl = 1,k
      ihmcon(kl) = kl
      iwrk2(kl) = kl
 4599 continue
c
c     loop over number of roots, iterations for each root.
c
      if(some) write(iw,9040)
  333 continue
c
      do 1315 iter=0,niter
c
      ipxt = ipxt + 1
      if(ipxt.eq.3) crit = critsv
c
c     check to see if ip = maxp, if
c     so then transform the first kst vectors in ci and ab
c     and start over with ip = k
c
      if (ip+ntcon.gt.maxp) then
         call tran_f(ci,nci,maxw1,ef,ip,ec,kst)
         call tran_f(ab,nci,maxw1,ef,ip,ec,kst)
         ip = kst
         do 1396 ii=1,maxp
            ec(ii,ii) = 1.0d+00
            ef(ii,ii) = 1.0d+00
            do 1398 jj=1,ii-1
               ec(ii,jj) = 0.0d+00
               ef(ii,jj) = 0.0d+00
               ec(jj,ii) = 0.0d+00
               ef(jj,ii) = 0.0d+00
 1398       continue
 1396    continue
      endif
c
c     make gradient vectors, put in ci(ip+1) -> ci(ip+ntcon)
c
      do 4588 kk=1,ntcon
      il = ihmcon(kk)
      gr(kk) = 0.0d+00
      do 80 ii=1,nci
         ci(ii,ip+kk) = 0.0d+00
         do 70 jj=1,ip
      ci(ii,ip+kk) = ci(ii,ip+kk) +
     *     ef(jj,il)*(ab(ii,jj)-el(il)*ci(ii,jj))
   70    continue
         gr(kk) = gr(kk) + (ci(ii,ip+kk)*ci(ii,ip+kk))
   80 continue
      gr(kk) = sqrt(gr(kk))
      if(some) then
         write(iw,9050) iter,el(il)+econst,gr(kk)
         call flshbf(iw)
      end if
 4588 continue
c
      if (some.and.ntcon.gt.1) write(iw,*)
      if (iter.eq.niter) then
         call tran_f(ci,nci,maxw1,ef,ip,ec,kst)
         goto 9890
      endif
c
c     check for convergence of any state, if converged, transform
c     all ip vectors in ci and ab, modify ihmcon and ntcon.
c     start with  ip = kst again.
c
      numc = 0
      do 4255 ii=1,ntcon
         iwrk2(ii) = ihmcon(ii)
 4255 continue
      do 4522 kk=1,ntcon
      il = ihmcon(kk)
      if (gr(kk).le.crit) then
         if(some) write(iw,9060) il,el(il)+econst,ipxt
         do 3233 ii=kk-numc,ntcon-numc
         iwrk2(ii) = iwrk2(ii+1)
 3233    continue
         numc = numc + 1
      endif
 4522 continue
c
      if (numc.gt.0) then
         call tran_f(ci,nci,maxw1,ef,ip,ec,kst)
         call tran_f(ab,nci,maxw1,ef,ip,ec,kst)
         ntcon = ntcon - numc
         do ii=1,ntcon
            ihmcon(ii) = iwrk2(ii)
         enddo
         do 74 ii=1,maxp
            do 75 jj=1,maxp
               ef(ii,jj) = 0.0d+00
               ec(ii,jj) = 0.0d+00
   75       continue
            ef(ii,ii) = 1.0d+00
            ec(ii,ii) = 1.0d+00
   74    continue
         ip = kst
         if (ntcon.ne.0) goto 333
         if(some) write(iw,*) 'all states converged.'
c
         call  fccspin(na,nb,iacon1,iacon2,ibcon1,ibcon2,
     *      isd,ido,ci,ab,k,nci,spin,iwrk1,iwrk2,maxw1,
     *           lgcom,nast,nbst,lsyma,lsymb,
     *           lsbs,lsbc,lspa,lspb,ldisb,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           x,nx,landet,lbndet,
     *           itga,itgb,iast,ibst,ktab,nsym)
c
         return
      endif
c
      do 68 jj=ip+1,ip+ntcon
         il = ihmcon(jj-ip)
         do 63 ii=1,nci
            ci(ii,jj) = ci(ii,jj)/(el(il) - q(ii))
   63    continue
   68 continue
c
c     if ms=0, impose restriction on the ci coefficients.
c
      if (na.eq.nb) then
      do 1111 iga=1,itga
         do 1122 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
c
            do 1133 igb=1,iga
            if (lgcom(igb,iga).ne.1) goto 1133
            icc1 = jpza1 + ldisb(ksym,igb,iga)
            icc2 = ldisb(jasym,iga,igb) + lspb(kka)
            nibs=nbst(igb)
            do 1144 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               ibpos = lsbc(kkb) + nibs
               if (ibpos.gt.kka) goto 1122
               ici1 = icc1+lspb(ibpos)
               if (ibpos.eq.kka) then
                  do 3232 kj=1,ntcon
                  nv=ihmcon(kj)
                  ips = (iwrk1(nv)/2)
                  if ((ips+ips).ne.iwrk1(nv)) ci(ici1,kj+ip) = 0.0d+00
 3232             continue
                  goto 1122
               endif
               ici2 = lspa(ibpos) + icc2
               do 3331 kj=1,ntcon
               nv = ihmcon(kj)
               is = (-1)**iwrk1(nv)
               ci(ici2,kj+ip) = is*ci(ici1,kj+ip)
 3331          continue
 1144       continue
 1133       continue
 1122    continue
 1111 continue
      endif
c
c     make the new vectors (ip+1 -> ip+ntcon).
c
c     assume the new vectors are bi, have to orthogonalize
c     these vectors to all others and then renormalize.
c
      do 97 kk=ip+1,ip+ntcon
         spin(kk) = spin(ihmcon(kk-ip))
         do 86 ii=1,kk-1
            rov = 0.0d+00
            do 81 jj=1,nci
               rov = rov + ci(jj,kk)*ci(jj,ii)
   81       continue
            do 90 jj=1,nci
              ci(jj,kk) = ci(jj,kk) - rov*ci(jj,ii)
   90       continue
   86    continue
c
         rnor = 0.0d+00
         do 40 ii=1,nci
            rnor = rnor + ci(ii,kk)*ci(ii,kk)
   40    continue
         rnor = sqrt(rnor)
         do 42 ii=1,nci
            ci(ii,kk) = ci(ii,kk)/rnor
   42    continue
   97 continue
c
      ip = ip + 1
c
c     now to return the new part of ab
c
      if (na.eq.nb) then
      if (kst.gt.1) then
      call fchc0ys(si1,si2,index,nact,na,nb,ci(1,ip),ab(1,ip),
     *           ntcon,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,
     *           iwrk1,iwrk2)
      else
      call fchc01s(si1,si2,index,nact,na,nb,ci(1,ip),ab(1,ip),
     *           q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,idim1,idim2,
     *           iwrk1,iwrk2)
      endif
c
      else
c
      if (kst.gt.1) then
      call fchcxys(si1,si2,index,nact,na,nb,ci(1,ip),ab(1,ip),
     *           ntcon,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
      else
      call fchcx1s(si1,si2,index,nact,na,nb,ci(1,ip),ab(1,ip),
     *           q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
      endif
      endif
c
      ip = ip + ntcon - 1
c
c     make the new matrix elements between the ci vectors.
c
      ix = 0
      do 103 ii=1,ip
          do 102 jj=1,ii
          ix = ix + 1
            f(ix) = 0.0d+00
            do 115 kk=1,nci
               f(ix) = f(ix) + ci(kk,ii)*ab(kk,jj)
  115       continue
  102    continue
  103 continue
c
c     diagonalize small matrix
c
      call evvrsp(-1,ip,ip,(ip*(ip+1))/2,maxw1
     *              ,f,b,iwrk2,el,ef,0,ierr)
      if (ierr.ne.0) then
         write(iw,*) 'error in small diagonalization'
         write(iw,*) ierr
         return
      endif
c
c     check to see if any states have skipped in
c
      do 700 ijk=1,kst
         do 705 ii=1,ntcon
            if (ijk.eq.ihmcon(ii)) goto 700
  705    continue
         idxc = 0
         pov = 0.0d+00
         do 713 jj=1,kst
            uit = 0.0d+00
            do 715 kk=1,ip
               uit = uit + ef(kk,jj)*ec(kk,ijk)
  715       continue
            pov = abs(uit)
            idxc = jj
  713    continue
         do 720 ii=1,ntcon
            if (jj.eq.ihmcon(ii)) then
               ihmcon(ii) = ijk
               goto 700
            endif
  720    continue
  700 continue
c
c     check to see where the spins occur now
c
      if (na.eq.nb) then
c
      do 800 ii=1,ip
         pov = 0.0d+00
         do 813 jj=1,ip
            uit = 0.0d+00
            do 823 kk=1,ip
               uit = uit + ef(kk,ii)*ec(kk,jj)
  823       continue
            if (abs(uit).gt.pov) then
               pov = abs(uit)
               idxc = jj
            endif
  813    continue
         if (idxc.ne.ii) then
            gr(ii) = spin(idxc)
         else
            gr(ii) = spin(ii)
         endif
800   continue
c
      do 786 kk=1,ip
         spin(kk) = gr(kk)
  786 continue
c
      endif
c
cji change below for state tracking.
c
c     check to see where the converged states are now.
c
      ncon = 0
      do 805 ii=1,k
         do 850 jj=1,ntcon
            if (ii.eq.ihmcon(jj)) goto 805
  850    continue
         do 852 jj=1,ip
            uit = 0.0d+00
            do 853 kk=1,ip
               uit = uit + ef(kk,jj)*ec(kk,ii)
  853       continue
            if (abs(uit).gt.0.6d+00) then
               ncon = ncon + 1
               iwrk1(ncon) = jj
               goto 805
            endif
  852    continue
c
  805 continue
c
      ntcon1 = 0
      do 860 ii=1,k
         do 864 jj=1,ncon
            if (ii.eq.iwrk1(jj)) goto 860
  864    continue
         ntcon1 = ntcon1 + 1
         ihmcon(ntcon1) = ii
  860 continue
c
      if (ntcon1.ne.ntcon) then
         write(6,*)
         write(6,*) 'converged states have skipped out'
         write(6,*) 'number of unconverged states now = ',ntcon1
         ntcon = ntcon1
      endif
c
c   end of check
c
cji end of changes
c
      do 543 ii=1,kst
         iwrk1(ii) = int(spin(ii) + 0.3d+00)
  543 continue
c
      do 55 ii=1,ip
         do 66  jj=1,ip
            ec(ii,jj) = ef(ii,jj)
   66    continue
   55 continue
c
 1315 continue
c
 9890 continue
c
      if (maswrk) write(iw,9080) il-1
c
      call  fccspin(na,nb,iacon1,iacon2,ibcon1,ibcon2,
     *   isd,ido,ci,ab,kst,nci,spin,iwrk1,iwrk2,maxw1,
     *        lgcom,nast,nbst,lsyma,lsymb,
     *        lsbs,lsbc,lspa,lspb,ldisb,
     *        lbox1,lbox2,lbox3,lbox4,lbox5,
     *        x,nx,landet,lbndet,
     *        itga,itgb,iast,ibst,ktab,nsym)
      istat=1
c
      return
c
 9005 format(/1x,'error, number of vectors stored=',i5,
     *           ' not equal to nci=',i5/
     *        1x,'this may be due to a garbage -civectr- file',
     *           ' left over from an earlier run.')
 9007 format(/1x,'initial ormas vectors read from disk')
 9010 format(/1x,'initial ormas vector guess time :',f13.1)
 9013 format(1x,'ormas spin calculation time     :',f13.1)
 9015 format(1x,'setting up calculation time     :',f13.1)
 9020 format(1x,'initial ormas ci iteration time :',f13.1)
 9040 format(/1x,'iteration',6x,'energy',11x,'gradient')
 9050 format(1x,i5,f20.10,f15.8)
 9060 format(/1x,'converged state',i5,' energy=',f20.10,' in',
     *           i5,' iters'/)
 9080 format(1x,'determinant ormas ci converged only',i4,' roots.')
      end
