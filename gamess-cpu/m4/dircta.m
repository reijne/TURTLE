c
c  $Author: jmht $
c  $Date: 2010-08-13 21:08:12 +0200 (Fri, 13 Aug 2010) $
c  $Locker:  $
c  $Revision: 6183 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dircta.m,v $
c  $State: Exp $
c
c ******************************************************
c ******************************************************
c             =   dircta    =
c ******************************************************
c ******************************************************
**==cizero.f
      subroutine cizero(first,lword,ocifp)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer excis,z95,zz11
_ELSEIF(convex,i8drct)
      integer *8 excis,z95,zz11
_ELSE
      REAL  excis,z95,zz11
_ENDIF
c
c   this routine replaces block data in direct
c
      logical ocifp,first
      character *8 textc1,texta,text1,text2
INCLUDE(common/sizes)
INCLUDE(common/natorb)
INCLUDE(common/symchk)
      common/pager/iparit,khbjco,nwolf,iclogo,icval,izval,
     *mxx,khbzsx,khbztx,khbcsx,khbctx,mumpzx,mumpcx,
     *myy,khbzsy,khbzty,khbcsy,khbcty,mumpzy,mumpcy
INCLUDE(common/hold)
INCLUDE(common/trial)
INCLUDE(common/dir_erg)
INCLUDE(common/timanb)
INCLUDE(common/cic)
      common/xynumb/xnumb(3),ynumb(3)
INCLUDE(common/helpr)
      common/symcon/nornor(26),ndimax
INCLUDE(common/corctl)
INCLUDE(common/ccntl)
INCLUDE(common/orb)
INCLUDE(common/cntrl)
INCLUDE(common/diactl)
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/table)
INCLUDE(common/disktl)
INCLUDE(common/presrt)
INCLUDE(common/mapp)
INCLUDE(common/ccepa)
INCLUDE(common/comrjb)
INCLUDE(common/modeci)
      common/apollo/pad1,pad2,pad3
_IF(cray,t3e,i8)
      integer pad1,pad2,pad3
_ELSEIF(convex,i8drct)
      integer *8 pad1,pad2,pad3
_ELSE
      REAL pad1,pad2,pad3
_ENDIF
c
      dimension text1(14),text2(14),msym(36)
      data msym/
     *1,
     *2,1,
     *3,4,1,
     *4,3,2,1,
     *5,6,7,8,1,
     *6,5,8,7,2,1,
     *7,8,5,6,3,4,1,
     *8,7,6,5,4,3,2,1/
      data text1/'(ij/kl)','(ij/kl)','(ij/kl)','(ij/ka)','(ij/ka)',
     *'(ij/ab)','(ij/ab)','(ia/jb)','(ia/jb)','(ia/jb)','(ia/ib)',
     *'(ia/bc)','(ab/cd)','other'/
      data text2/'n  /n','n-1/n-1','n-2/n-2','n-1/n','n-2/n-1',
     *'n-1/n-1','n-2/n-2','n-1/n-1','n-2/n-2','n-2/n','n-2/n',
     *'n-2/n-1','n-2/n-2','other'/
      data textc1,texta/'c1','a'/
_IF(cray)
      data excis,z95,zz11/137b,137b,13b/
_ELSEIF(t3e,convex,i8,i8drct)
      data excis,z95,zz11/z'5f',z'5f',z'b'/
_ELSEIF(ibm)
      data excis,z95,zz11/z5f,z5f,zb/
_ELSE
_ENDIF
c
_IFN(ibm,convex,cray,t3e,i8,i8drct)
c
_IF(bits8)
      call pad8(1,pad1)
      call pad8(2,pad2)
      call pad8(3,pad3)
      call pad8(95,excis)
      call pad8(95,z95)
      call pad8(11,zz11)
_ELSE
      pad1 = pad(1)
      pad2 = pad(2)
      pad3 = pad(3)
      excis = pad(95)
      z95 = pad(95)
      zz11 = pad(11)
_ENDIF
_ELSEIF(ibm)
      pad1 = pad(1)
      pad2 = pad(2)
      pad3 = pad(3)
_ELSE
      pad1 = 1
      pad2 = 2
      pad3 = 3
_ENDIF
c
_IF(cyber205)
      if (nd200.ne.200) call caserr('nd200 =>  200 for assembler IHINT')
_ENDIF
c
c//// modeci //// 
      in_store = .false.
      out_store = .false.
      odiag_1 = .false.
      odiag_2 = .false.
      odiag_3 = .false.
c////  cepa  ////
      ifcepa = 0
      if (first) then
        icepa = -999
        itcepa = 0
        critcep = 0.0d0
c////  combuen  ////  mrdci-option
        ifbuen = .false.
        cradd(1) = 0.01d0
        cradd(2) = 0.0016d0
        cradd(3) = 0.01d0
        iex1 = z95
        iex2 = zz11
        ncv = 0
        nbuenk = 0
        mbuenk = 2
        weighb = .95d0
c////  natorb  ////
        latorb = .true.
        ispind = .true.
        ispacg = 0
        isping = 0
_IF1()        iprwh=1
_IF1()        iprcon=1
        iprwh = 199999999
        iprcon = 199999999
      end if
c////  pager  ////
      iparit = -1
c////  hold  ////
      if (first) then
        screan = .true.
        mxstri = 10
        nesvac = 0
        nesuno = 0
        nesdue = 0
        nesref = 0
_IF(cray,t3e,i8)
       call setsto(80,0,mestri)
_ELSE
       call setstz(80,0,mestri)
_ENDIF
      end if
c ///// corctl /////
      ngot = lword
c//// symchk ////
      lsymd = .false.
      if(.not.ocifp) lsadap = .false.
      if (first) then
        crtsym = 1.0d-5
        nrep = 0
        excit = excis
c..     since usually adapt info is taken no default symdet
        lsymd = .false.
        lsadap = .true.
c////  trial  ////
        nvec = 1
        ivec(1) = 1
        tvec(1) = 1.0d0
        mxvec = 20
        iselvc = -1
        istvc = 1
        iprvc = 0
        nselvc = -1
      end if
c////  diactl  ////
      if (first) mrstor = 0
      msort = .false.
      lvar = .false.
      instor = .true.
      moded = .true.
      modev = .true.
      modest = .true.
      lcall = 50
c//// cntrl ////
      nbig = 0
      mbig = 0
      nsingl = 0
      ntripl = 0
      ntotal = 0
      jbig = 0
      if (first) then
       nref0 = 0
      endif
c//// symcon ////
      ndimax = 4
c//// erg ////
      icyc = 0
      if (first) then
        do 40 i = 1 , 3
          vdstpr(i) = 1000.0d0
 40     continue
        malt = .false.
        maxcyc = 50
        moddia = -1
        thresh = 3.0d-5
        eshift = 0.0d0
        gprin = 3.0d-2
        scream = 0.71d0
        mprin = 1000
        qester = 1.0d-8
      end if
      tester = 0.0d0
      e11 = 0.0d0
c//// auxh ////
      do 50 i = 1 , 8
        nspips(i) = 0
        nspipt(i) = 0
 50   continue
      do 60 i = 1 , 3
        maxsin(i) = 0
        maxpa(i) = 0
 60   continue
      maxpa(4) = 0
      do 70 i = 1 , 8
        nconst(i) = 0
        ncond(i) = 0
 70   continue
c//// presrt ////
      do 80 i = 1 , 78
        link(i) = -1
        nwbuf(i) = 0
 80   continue
      iblkp = 0
      indxi = 0
      indxj = 0
      indxk = 0
      indxl = 0
c////  ccntl ////
      nmin1 = 0
      nmin2 = 0
      if (first) then
        nref = 0
        mspin = 0
        nelec = 0
c////   orb  ////
        nint = 0
        next = -1
        ninnex = 0
        nw = 5
c///   helpr   ///
        call setsto(nd200,1,iro)
      end if
      m1e = 0
      m2e = 0
      m1ia = 0
      m2coul = 0
      m2exch = 0
      m2psup = 0
c///  sym    ////
      nubtr = 0
      nubsq = 0
      nubsi = 0
      nvmax = 0
      norbsi(1) = 0
      norbtr(1) = 0
      norbsq(1) = 0
      call setsto(64,0,ibassi)
c.... ****  cray-1 & cyber 205  ****
      n127 = 127
      nincr = 64
      maxsq = maxorb - 1
c.... **** cyber 174 ****
c     data n127,nincr,maxsq/31,32,127/
c////  table ////
      if (first) then
        jfrep = 1
        ifrep = 0
        call isquar(mult,msym,8,8)
        rname(1) = texta
        sname = textc1
      end if
c//// xynumb ////
      xnumb(1) = 0.0d0
      xnumb(2) = 1.0d0
      xnumb(3) = 0.0d0
      ynumb(1) = 0.0d0
      ynumb(2) = 0.0d0
      ynumb(3) = 1.0d0
c//// disktl ////
_IF(parallel)
c MPP
      numcip = numci
      numcic = 16
c     numcic = numcip
c MPP
_ENDIF
      call setsto(199,-1,index)
c... index 1 --- coupling coefficients
c... index 2 --- (ij/ka) integrals
c... index 3 --- fock(ia) for doublets and n-2 states
c... index 4 --- (ia/ib) + (ij/ab) integrals
c... index 5 --- (ia/jb)   p and q   supermatrices
c... index 6 --- fock(ab) for n-2 states
c... index 7 --- fock(ab) for doublets
c... index 8 --- (ia/jb) integrals
c... index 9 --- (ia/bc) p and q supermatrix
c... index 10 to 17 (ab/cd) p supermatrices
c... index 18 to 25 (ab/cd) q supermatrices
c... index 26 to 145 c,z vectors (max. 60 of them)
c... index 146,147 highest possible disk address for best c,z
c... index 148,149 best c,z vectors
c... index 197 coupling coeffs for space and spin nos
c... index 198 diagonal of hamiltonian(doublets)
c... index 199 diagonal of hamiltonian
c////  map  ////
      if (first.or.ocifp) modei = .true.
      if (ocifp) then
       do loop = 1 , nd200
         mapie(loop) = loop
         mapei(loop) = loop
       enddo
      endif
c////    timanb     ////
      do 90 loop = 1 , 14
        timsec(loop) = 0.0d0
        text11(loop) = text1(loop)
        text22(loop) = text2(loop)
 90   continue
      return
      end
**==nsior.f
      function nsior(conf)
      implicit REAL  (a-h,o-z),integer (i-n)
c     count the number of singly occupied orbitals in
c     the configuration described by 'nw' words
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt
_ENDIF
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSE
      REAL conf
_ENDIF
      dimension conf(*)
INCLUDE(common/cdcryz)
INCLUDE(common/orb)
c.... count number of electrons (nel) and number of doubly occupied
c.... orbitals  (ndoc)  ==> number of singles : nsin = nel-2*ndoc
      nel = 0
      ndoc = 0
      do 20 i = 1 , nwm1
        nel = nel + popcnt(conf(i))
_IF(bits8)
        call dand8(conf(i),n32m,tmp8)
        ndoc = ndoc + popcnt(tmp8)
_ELSEIF(cray,t3e,i8,i8drct)
        ndoc = ndoc + popcnt(iand(conf(i),n32m))
_ELSE
        ndoc = ndoc + popcnt(dand(conf(i),n32m))
_ENDIF
 20   continue
      nsior = nel - ndoc - ndoc
      return
      end
**==direcg.f
      subroutine direcg(cr,icr,energy)
_IF(parallel)
c MPP
c MPP    parallel version
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer (i-n)
c
c...  direct CI module
c...  V.R. Saunders, J.H. van Lenthe
c...  Daresbury (1980 ..)
c
c     **note instore and outstore mode answers differ very slightly **
c
_IF(cray,t3e,i8)
      integer icr
_ELSEIF(convex,i8drct)
      integer *8 icr
_ELSE
      REAL icr
_ENDIF
      logical oprn
      character *8 text,txcepa
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/cntrl)
INCLUDE(common/auxh)
INCLUDE(common/diactl)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
      common/stoctl/ntotw,khbz,khbc,khbb
INCLUDE(common/disktl)
INCLUDE(common/timez)
INCLUDE(common/restar)
INCLUDE(common/dir_erg)
INCLUDE(common/corctl)
INCLUDE(common/comrjb)
_IF(parallel)
INCLUDE(common/blksiz2)
INCLUDE(common/atmblk)
_ENDIF
INCLUDE(common/modeci)
INCLUDE(common/prints)
INCLUDE(common/mr_mp2)
INCLUDE(common/mr_comci)
INCLUDE(common/ccepa)
INCLUDE(common/mr_cepctl)
INCLUDE(common/mr_pairs)
INCLUDE(common/disc)
      common/remc_ed19/ninini,ipipip
      dimension cr(*),icr(5,*)
      dimension h(3),h2(3),hci(3),xn1(2)
c...   h/h2 are the effective elements (incl. shift); hci : straight ci
      logical ocepa_first
      data h,h2/6*0.0d0/,xn1/1.0d0,1.0d0/
      data text/'       '/,txcepa/'  cepa  '/
      data ntdiis,diiscr/987654321,1.0d0/
c
      call mpmask(cr,nref,nwm1,nw)
      z11 = 0.0d0
      ocepa_first = .true.
      ncycep = 0 
c
_IFN(ga)
      if (icepa.ge.0) 
     1 call set_docfind(cr,cr(nref*nw+1),nref,mrdoc,'setmask')  
_ENDIF
c
c....  c o n f i g u r a t i o n    g e n e r a t i o n
c
      call itint(cr,nw)
      if (next.ne.0) call ext(cr,nw)
      call compre(cr,nw)
      kres = nwp1*nw
      nuse = kres
c
c....  c o n f i.g u r a t i o n   s e l e c t i o n
c
      call selct(cr,cr(kres+1),nw)
c
c.... get everything in proper position
c
      ntotw = nlast2*5
_IF(i8drct)
      do i=1,nlast2
       do j=1,5
        icr(j,i) = icr(j,i+nwp1)
       enddo   
      enddo   
_ELSEIF(linux)
      call dcopy(ntotw,cr(kres+1),1,cr,1)
_ELSE
      call fmove(cr(kres+1),cr,ntotw)
_ENDIF
c
c.... sort configurations for internally contracted stuff
c
      if (mp.ne.0.or.icci.ne.0) call mpsrcc(cr,cr(nlast2*5+1))
c
c....  p r i n t  c o n f i g u r a t i o n s
c
      call prsymb(icr)
_IFN(ga)
      if (icepa.ge.0) call cepinf(cr,5)
_ENDIF
c
c.... make space to store the reference vector for the internally
c.... contracted stuff
c
      krefmp = 0
      if (mp.ne.0.or.icci.ne.0) then
         krefmp = ntotw
         ntotw  = ntotw + nrfcsf
      endif
c
c---------------  e n d   o f   c o n g e n  ---------------------------
c.... h - build
      call concal
c
c... pre sort
c
      if (mp.eq.2) then
c...     In the MP2 case there are considerable savings 
c...     so call some special routines
         call mpgijk(cr,icr)
         call mphbil(cr,icr)
      else
         call gijkl(cr,icr)
         call hbill(cr,icr)
      endif
      if (next.ne.0) then
c... post sort
        call gijka(cr,icr)
        call gijab(cr,icr)
        call giabc(cr)
_IF(parallel)
c
c     p-sort file settings
c
      call setsrtp(nsz,o255i,
     +  nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     +  nsz510,nsz341,nsz342,nsz680)
      nsstat = -1
      nsmax = 0
c
_ENDIF
        call gabcd(cr,icr)
      end if
      nxblk = iposun(numci)
c
c...  save the configurations on numci
c
      nblcf  = nxblk
      iblkcc = nxblk
      call wrt3(cr,5*nlast2,nxblk,numci)
      nxblk = iposun(numci)
c
c...  The config information is different for the uncontracted and
c...  the contracted stuff. So we store both tables on file and load
c...  them in the appropriate locations at the right time.
c
      if (mp.ne.0.or.icci.ne.0) then
c...     set block to store config at for uncontracted ci
         call sciblk(nxblk)
         idiff = nxblk - nblcf
         nxblk = nxblk + idiff
c...     set block to store config at for contracted mp and ci
         call smpblk(nxblk)
         nxblk  = nxblk + idiff
         call mpstar(cr,cr,cr)
         nxblk=iposun(numci)
      end if

c---------------  e n d   o f   h - b u i l d  -------------------------
c
c... diagonalize
c
      top = cpulft(1)
      tlefts = timlim - top
      call diabas(icr,iwr)
      khbccc = khbc
      if (mp.ne.2.and.maxcyc.gt.0) call diagon(cr,icr,iwr)
      if (icci.ne.0) call sttria
      if (icepa.ge.0) call inicepa(cr,icr)
_IFN(ga)
      call genc1(cr,icr,iwr)
_ELSE
      oadd(numcz)=.false.
      call switch_oswed(numcz,.true.)
      call genc1(cr,icr,iwr)
      oadd(numcz)=.true.
_ENDIF
_IF(parallel)
      call switch_oswed(numcz,.false.)
      call genz(cr,icr,index(26))
      call pg_synch(1234)
      call switch_oswed(numcz,.true.)
_ELSE
      call genz(cr,icr,index(26),h11)
_ENDIF
      if (instor) call dumpcz(icr,cr(khbz+1),index(27))
      modev = nnst.eq.0
      moded = nmin1.eq.0
      modest = nmin2.eq.0
_IF(parallel)
      if (lcall.ne.2) then
       call davips(cr,icr,cr(khbccc+1),cr(ntotw+1),iwr)
       call timprt(iwr)
      else
       call caserr('2X2 disabled by J.H. van Lenthe')
      endif
_ENDIF
_IFN(parallel)
c
c...  At this point we have enough to do the MP stuff
c
      if (mp.ne.0) then
         if (mp.le.3) then
            call mpmr(cr,cr,h11)
         else
            call mpmr2(cr,cr)
            call getcz(cr,cr(khbc+1),index(26))
            call genz(cr,icr,index(26),h11)
            if(instor)call dumpcz(icr,cr(khbz+1),index(27))
         endif
c...     if no ci to be started we are through
         if (maxcyc.eq.0) go to 4444
c...     else we have vector in index(26) and z in index(27)
c...     and we can start ci
      else if (icci.ne.0) then
         write(iwr,2000)
         call davidc(cr,icr,cr(kcpcc+1))
         go to 4444
      end if
      write(iwr,2000)
2000  format(/3x,77('=')/
     *3x,'cycle',2x,'time(secs)',10x,'energy(au)',14x,'tester',12x,
     *'variance'/3x,77('=')/)
c
      if (lcall.eq.3) then
c...    lcall =3 switched on 3*3 / parameters are fixed here
c...    start right away (unless things look **really** bad)
        ntdiis = 1
        diiscr = 1.0d0
c...    temp disable for cepa (may be fixed with mrdcepa
        if (icepa.gt.0) then
           ntdiis = 987654321
        else 
           go to 20
        end if
      end if
      if (lcall.ne.2) then
        if (instor) then
          call david3(cr,icr,cr(khbccc+1),cr(ntotw+1),iwr)
          go to 50
        else
          naw = (ngot-ntotw)/ntotal
          if (naw.ge.2) then
            if (odiag_1) then
              call david1(cr,icr,cr(ntotw+1),iwr)
            else
              call david2(cr,icr,cr(khbccc+1),cr(ntotw+1),iwr)
            endif
            go to 50
          else if (naw.ge.1) then
            call david1(cr,icr,cr(ntotw+1),iwr)
            go to 50
          end if
        end if
      end if
c... start iterative loop
 20   if (icepa.ge.0 .and. tester.le.critcep .and. icyc.ge.itcepa)
     +    ifcepa = 1
      h11o = h11
c
      if (ifcepa.gt.0) then
        text = txcepa
        call cntrci(icr,cr(kcepc+1),cr(kcepz+1),index(26),index(27),
     1  cr(kcpcc+1),1,2,h11ci)
        hci(1) = h11ci
        call dumpcc(cr(kcpcc+1),1,1)
        call cepair(cr,icr,index(26))
        call cepac(cr(khbc+1),icr)
c...    shifted h11
        call hh2cep(cr(kcpcc+1),ncpcc,icr,h11,hci,z11,h2ci,xn1,1)
c...    cepa micro-iteration ?? (first time special c1.c2 not prepared)
        ncycep = ncycep + 1
        if (ocepa_first.and.mcycep.gt.0) then
           if (instor) call getcz(cr,cr(kcepc+1),index(28))
           if (instor) call getcz(cr,cr(kcepz+1),index(29))
           go to 117
        else if (ncycep.le.mcycep.and.
     *           dabs(h11-h11o).gt.tester*crycep) then
           go to 113
        end if
c
        ncycep = 0
c...      restore c/z vector for instore
        if (instor) call getcz(cr,cr(khbc+1),index(26))
        if (instor) call getcz(cr,cr(khbz+1),index(27))
      end if
cmp
cmp   cg/diis acceleration of 2*2 davidson (genc3) dumped the c/z vector
cmp
cmp   run genc2 once just to get a tester and check if diisci may start
cmp
      if (icyc.eq.1.and.icyc.ge.ntdiis.and.ngot-ntotw.ge.2*ntotal) then
         call genc2(cr,icr)
         if (instor) call getcz(icr,cr(khbc+1),index(26))
         if (instor) call getcz(icr,cr(khbz+1),index(27))
      end if
c
      if (ifcepa.ne.0.and.ntdiis.ne.123456789) then
        ntdiis = 123456789
cjvl         write(iwr,604)
cjvl 604     format(/' DIIS not yet compatible with CEPA ')
      endif
c
      if (tester.lt.diiscr.and.icyc.ge.ntdiis) then
         if (ngot-ntotw.lt.2*ntotal) then
            write(iwr,103)
103         format(/' *** DIIS not started due to lack of core ***'/)
            ntdiis = 98765432
         else
            call diisci(cr(ntotw+1),cr(ntotw+ntotal+1),index(30),
     *                  ntotal,index(26),index(27),
c    *                  h11,cr,ifcepa,nref,numci,iwr)
     *                  h11,cr,nref,numci,iwr)
            if (instor) call getcz(icr,cr(khbc+1),index(26))
            if (instor) call getcz(icr,cr(khbz+1),index(27))
         end if
      end if
cmp
      call genc2(cr,icr)
      s22s = s22
c
      top = cpulft(1)
      if (.not.oprint(28))
     + write (iwr,6010) icyc , top , h11 , tester , z11 , text
      if (tester.lt.thresh .and. (icepa.lt.0 .or. ifcepa.eq.1)) 
     + goto 777
      if (icyc.ge.maxcyc) goto 440
          tlefti = timlim - top
         if (.not.(tlefti.ge.2.0d0*(tlefts-tlefti)/dfloat(icyc))) 
     + goto 440
c
c...       be sure not to do jacobi-davidson and any kind of cepa,
c...       because it fails (you change the operator behind 
c...       JD's back)
c
              if (icepa.lt.0) then
c
c...
c...  possible jacobi davidson preconditioning on vacuum space
c...  here khbill is the first-1 word available and c starts there
c...  (and is ntotal long)
c
c...   read original pert vector first
c...   if in store mode all is different
c...   Since jacdav decides itself if it should be called, we have to
c...   calculated some (vacuum) info, which may not be needed (sigh)
c
             njac = 2*ntotal
             khbcj = khbc
             khbzj = khbz
             if (.not.instor) then
                njac = 2*nval
                khbcj = ntotw
                khbzj = khbcj + nval
c...            read c-vector (vacuum part) )in store in instor mode)
                call rdedx(cr(khbcj+1),nval,index(28),numci)
             end if
c
c...         calculate the vacuum s22,s12,h12 ...
c...         (the z-vector is available) 
c
             s22v = ddot(nval,cr(khbcj+1),1,cr(khbcj+1),1)
             call rdedx(cr(khbzj+1),nval,index(26),numci)
             s12v = ddot(nval,cr(khbcj+1),1,cr(khbzj+1),1)
             call rdedx(cr(khbzj+1),nval,index(27),numci)
             h12v = ddot(nval,cr(khbcj+1),1,cr(khbzj+1),1)
c
             call gdvdmr(index(1),index(199),index(26),index(27),
     +                   ntotw+njac,nval,maxpa(4),
     +                   cr(khbcj+1),cr(khbzj+1),h11,cr,icr,
     +                   numci,ncycg)
c
c...     adjust s22, h12, s12 if jacobi-davidson is done
c
             if (ncycg.ne.0) then
                s22 = s22 - s22v + 
     +                ddot(nval,cr(khbcj+1),1,cr(khbcj+1),1)
                call rdedx(cr(khbzj+1),nval,index(26),numci)
                s12 = s12 - s12v + 
     +                ddot(nval,cr(khbcj+1),1,cr(khbzj+1),1)
                call rdedx(cr(khbzj+1),nval,index(27),numci)
                h12 = h12 - h12v + 
     +                ddot(nval,cr(khbcj+1),1,cr(khbzj+1),1)
c...            write the c-vector (vacuum part) 
                if (.not.instor) 
     1             call wrt3(cr(khbcj+1),nval,index(28),numci)
              end if
c
                     end if
c
_IF(parallel)
            call genz(cr,icr,index(28))
_ELSE
            call genz(cr,icr,index(28),h22)
_ENDIF
            hci(3) = h22
c
            if (icepa.lt.0) goto 300
            if (instor) call dumpcz(cr,cr(kcepc+1),index(28))
            if (instor) call dumpcz(cr,cr(kcepz+1),index(29))
            if (ifcepa.le.0) go to 300
117         continue
c...    <c2.c2>  /  <c2.z2>
            call cntrci(icr,cr(kcepc+1),cr(kcepz+1),
     *       index(28),index(29),cr(kcpcc+1),1,2,temp)
            call dumpcc(cr(kcpcc+1),2,2)
c
113         continue
            if (lprcep.gt.0.and.ncycep.gt.0) 
     *      write(6,114) ncycep,h11,h11o,dabs(h11-h11o)
114   format(/' cepa micro iteration ',i3,' new energy ',e20.12,
     *'   prev. energy ',e20.12,' diff ',e10.3)
c...    <c2.c1>      /  s12 
            if (icpsto) call getcz(cr,cr(kcepz+1),index(26))
            call cntrci(icr,cr(kcepc+1),cr(kcepz+1),
     *       index(28),index(26),cr(kcpcc+1),0,1,s12)
c...    <z2.c1> (for vmin) /  h12
            if (icpsto) call getcz(cr,cr(kcepc+1),index(29))
            call cntrci(icr,cr(kcepc+1),cr(kcepz+1),
     *       index(29),index(26),cr(kcpcc+1),0,3,hci(2))
c...           <c2.z1> / <z1.z1>  (summed => z11)
               if (icpsto) call getcz(cr,cr(kcepc+1),index(27))
               if (icpsto) call getcz(cr,cr(kcepz+1),index(28))
               call cntrci(icr,cr(kcepz+1),cr(kcepc+1),
     *          index(27),index(28),cr(kcpcc+1),4,2,temp)
               call setstr(ncpcc,1.0d0,cr(kcpcc+4*ncpcc+1))
c...           redetermine z11
               h2(1) = ddot(ncpcc,cr(kcpcc+4*ncpcc+1),1,
     *          cr(kcpcc+3*ncpcc+1),1)
            if (moddia.gt.0) then 
c...           compute z12 and z22 for variance minimizer
               if (instor) call getcz(icr,cr(khbz+1),index(29))
               call genz3(cr,icr)
               h2(2) = z12
               h2(3) = z22
            end if
c...
            call dumpcc(cr(kcpcc+1),2,1)
c...   compute shifted 2*2 h-matrix  (and z-matrix for vmin)
            call hh2cep(cr(kcpcc+1),ncpcc,icr,h,hci,h2,h2,xn1,2)
            h12 = h(2)
            s22 = s22s
            h22 = h(3)
            z11 = h2(1) - h11*h11
            z12 = h2(2)
            z22 = h2(3)
            ocepa_first = .false.
c...
300         s22i = 1.0d0/s22
            s22isq = dsqrt(s22i)
            s12 = s12*s22isq
            s12s12 = s12 + s12
            h12 = h12*s22isq
            s12sq = s12*s12
            s121 = 1.0d0/(1.0d0-s12sq)
            s121sq = dsqrt(s121)
            h22 = (h22*s22i+s12sq*h11-s12s12*h12)*s121
            h12 = (h12-s12*h11)*s121sq
            h12h12 = h12 + h12
            if (moddia.le.0) then
              bbb = b22(h12,h22)
            else
c... minimize variance
ccepa..............................................................
              if (ifcepa.le.0) call genz3(cr,icr)
ccepa..............................................................
              z12 = z12*s22isq
              x11 = z11 + h11*h11
              z22 = (z22*s22i+s12sq*x11-s12s12*z12)*s121 - h22*h22
              z12 = (z12-s12*x11)*s121sq
              b0 = z12 - h12h12*h11
              b1sq = h11 - h22
              top = h12h12*h12h12 - b1sq*b1sq
              b4 = z22 - z11
              b1 = b4 - top
              b2 = 3.0d0*h12h12*b1sq
              b3 = b4 + top
              b4 = h12h12*h22 - z12
              bbb = 0.0d0
              bsq = b2 + b2
              b33 = 3.0d0*b3
              b44 = 4.0d0*b4
c... find appropriate root of quartic equation
              do 30 loop = 1 , 100
                top = (b1+bbb*(bsq+bbb*(b33+bbb*b44)))
                b1sq = (b0+bbb*(b1+bbb*(b2+bbb*(b3+bbb*b4))))
     +                 /(0.01d0+dabs(top))
                bbb = bbb - b1sq/(1.0d0+10.0d0*dabs(b1sq))
                if (dabs(b1sq).lt.1.0d-11) go to 40
 30           continue
            end if
 40         bsq = bbb*bbb
            b1 = 1.0d0/(1.0d0+bsq)
            b1sq = dsqrt(b1)
            h11 = (h11+bsq*h22+bbb*h12h12)*b1
            aaa = (1.0d0-bbb*s12*s121sq)*b1sq
            bbb = bbb*b1sq*s121sq*s22isq
ccepa
      if (instor.and.ifcepa.gt.0.and.icepa.ge.0) 
     *  call getcz(cr,cr(khbc+1),index(28))
      if (instor.and.ifcepa.gt.0.and.icepa.ge.0) 
     *  call getcz(cr,cr(khbz+1),index(29))
ccepa
            call genc3(cr,icr)
            if (malt) eshift = -eshift
            go to 20
777   if (.not.oprint(28))
     +  write (iwr,6020) icyc
      irest = 0
      goto 50
440     write (iwr,6030)
        irest = 9
        tim = timlim + 0.1d0
 50   call timprt(iwr)
_ENDIF
c---------------- end of diagonalizer ----------------------------------
4444  call revise
      mrstor = 0
      call vecpri(cr,icr,iblkcc,iwr)
_IFN(parallel)
      if (ifcepa.gt.0) call prcepa(iwr)
_ENDIF
      energy = h11
      oprn = nprint.ne.-5
      call nosas(cr,icr,oprn)
c
c...  print (mr)-size-consistency corrections if no cepa
c
      if (ifcepa.le.0.and.(mp.eq.0.or.maxcyc.gt.0)) then
        call sccor(cr(ntotw+1),icr,ngot-ntotw)
_IF(ga)
        call pg_synch(129)
_ENDIF
      end if
c
      call revind
      call clredx
      return
 6010 format (1x,i7,f12.3,3e20.11,a8)
 6020 format (///20x,26('*')/20x,'convergence at cycle',i6/20x,26('*'))
 6030 format (//20x,14('*')/20x,'no convergence'/20x,14('*'))
      end
**==compre.f
      subroutine compre(iconf,nword)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf,z1
_ELSEIF(convex,i8drct)
      integer *8 iconf,z1
_ELSE
      REAL iconf,z1
      logical eq
_ENDIF
      dimension iconf(nword,*)
INCLUDE(common/hold)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
_IF(cray,t3e,convex,i8,i8drct)
      data z1/1/
_ELSE
_IF(bits8)
      call pad8(1,z1)
_ELSE
      z1 = pad(1)
_ENDIF
_ENDIF
c...   keep track of nref
      nreft = 0
      nv = 0
      nstrt1 = nnst + 1
      nlast1 = nnst + nmin1
      nstrt2 = nlast1 + 1
      nlast2 = nlast1 + nmin2
      do 20 loop = 1 , nnst
_IF(cray,t3e,convex,i8,i8drct)
        if (iand(iconf(nw,loop),z1).eq.z1)then
_ELSEIF(bits8)
        call dand8(iconf(nw,loop),z1,tmp8)
        if(eq(tmp8,z1)) then
_ELSE
        if (eq(dand(iconf(nw,loop),z1),z1))then
_ENDIF
          if (irrcon(iconf(1,loop)).eq.1) then
            if (ihold(iconf(1,loop),mestri(1,1,1),minstr(1,1),
     +          maxstr(1,1),nesvac).ne.0) then
              nv = nv + 1
              if (loop.le.nref) nreft = nreft + 1
_IF(i8drct)
                do i = 1,nwm1
                 iconf(i,nv) = iconf(i,loop)
                enddo
_ELSEIF(linux)
              call dcopy(nwm1,iconf(1,loop),1,iconf(1,nv),1)
_ELSE
              call fmove(iconf(1,loop),iconf(1,nv),nwm1)
_ENDIF
            end if
          end if
        end if
 20   continue
      nnst = nv
      nref = nreft
      if (nmin1.ne.0) then
        do 30 loop = nstrt1 , nlast1
_IF(cray,t3e,convex,i8,i8drct)
          if(iand(iconf(nw,loop),z1).eq.z1) then
_ELSEIF(bits8)
          call dand8(iconf(nw,loop),z1,tmp8)
          if(eq(tmp8,z1)) then
_ELSE
          if(eq(dand(iconf(nw,loop),z1),z1))then
_ENDIF
            if (ihold(iconf(1,loop),mestri(1,1,2),minstr(1,2),
     +          maxstr(1,2),nesuno).ne.0) then
              nv = nv + 1
_IF(cray,t3e,convex,i8,i8drct)
              iconf(nw,nv) = irrcon(iconf(1,loop))
_ELSEIF(bits8)
              call pad8(irrcon(iconf(1,loop)),iconf(nw,nv))
_ELSE
              iconf(nw,nv) = pad(irrcon(iconf(1,loop)))
_ENDIF
_IF(i8drct)
              do i=1,nwm1
               iconf(i,nv) = iconf(i,loop)
              enddo
_ELSEIF(linux)
              call icopy(nwm1*2,iconf(1,loop),1,iconf(1,nv),1)
_ELSE
              call fmove(iconf(1,loop),iconf(1,nv),nwm1)
_ENDIF
            end if
          end if
 30     continue
        nmin1 = nv - nnst
      end if
      nwp1 = nv
      if (nmin2.ne.0) then
        do 40 loop = nstrt2 , nlast2
          nwp1 = nwp1 + 1
_IF(cray,t3e,convex,i8,i8drct)
          iconf(nw,nwp1) = irrcon(iconf(1,loop))
_ELSEIF(bits8)
          call pad8(irrcon(iconf(1,loop)),iconf(nw,nwp1))
_ELSE
          iconf(nw,nwp1) = pad(irrcon(iconf(1,loop)))
_ENDIF
_IF(i8drct)
          do i=1,nwm1
           iconf(i,nwp1) = iconf(i,loop)
          enddo
_ELSEIF(linux)
          call dcopy(nwm1,iconf(1,loop),1,iconf(1,nwp1),1)
_ELSE
          call fmove(iconf(1,loop),iconf(1,nwp1),nwm1)
_ENDIF
 40     continue
        nmin2 = nwp1 - nv
      end if
      return
      end
**==ihold.f
      function ihold(iconf,mestri,minstr,maxstr,nes)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer mestri,iconf,mas
_ELSEIF(convex,i8drct)
      integer *8 mestri,iconf,mas
_ELSE
      REAL mestri,iconf,mas
_ENDIF
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt
_ENDIF
      dimension iconf(*),mestri(2,*),minstr(*),maxstr(*)
c... to apply occupation number restrictions
INCLUDE(common/orb)
      ihold = 1
      if (nes.ne.0) then
        do 30 loop = 1 , nes
          icount = 0
          do 20 koop = 1 , nwm1
_IF(bits8)
            call dand8(iconf(koop),mestri(koop,loop),mas)
_ELSEIF(cray,t3e,convex,i8,i8drct)
            mas = iand(iconf(koop),mestri(koop,loop))
_ELSE
            mas = dand(iconf(koop),mestri(koop,loop))
_ENDIF
            icount = icount + popcnt(mas)
 20       continue
          if (icount.lt.minstr(loop) .or. icount.gt.maxstr(loop))
     +        go to 40
 30     continue
      end if
      return
 40   ihold = 0
      return
      end
**==cpack.f
      subroutine cpack(iconf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf,i013
_ELSEIF(convex)
      integer *8 iconf,i013
_ELSEIF(i8drct)
      integer *8 iconf,i013
      integer *8 shift
      external shift
_ELSE
      REAL iconf,i013
_ENDIF
c     form the occupation-scheme of a configuration
c     ( 2 bits/orbital ; 1 electron/bit)
c     conf : resulting occ. pattern ((nint-1)/n32+1 words)
c     nint : number of internal orbitals
c     iorb : occupations for variably occupied orbitals
c           which are read here
      dimension iconf(*),i013(3)
INCLUDE(common/work)
INCLUDE(common/orb)
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data i013/0,1,3/
_ELSE
_IF(bits8)
      call pad8(0,i013(1))
      call pad8(1,i013(2))
      call pad8(3,i013(3))
_ELSE
      i013(1) = pad(0)
      i013(2) = pad(1)
      i013(3) = pad(3)
_ENDIF
_ENDIF
      ns = 64
      iw = 1
_IF(cray,t3e,i8)
      call setsto(nwm1+nwm1,0,iconf)
_ELSE
      call setstz(nwm1+nwm1,0,iconf)
_ENDIF
      do 20 i = 1 , nint
        if (ns.le.0) then
c.... new word
          ns = 64
          iw = 2
        end if
        call rinpi(iorb)
        if (iorb.gt.2 .or. iorb.lt.0) call errout(istrt(jrec))
        ns = ns - 2
_IF1() PAUL Try ishftc version here on cray
_IF(cray,t3e,convex,i8)
        iconf(iw) = ior(iconf(iw),ishftc(i013(iorb+1),ns,64))
_ELSEIF(i8drct)
        iconf(iw) = ior(iconf(iw),shift(i013(iorb+1),ns))
_ELSEIF(bits8)
        call shift8(i013(iorb+1),ns,tmp8)
        call dor8(iconf(iw),tmp8,iconf(iw))
_ELSE
        iconf(iw) = dor(iconf(iw),SHIFT(i013(iorb+1),ns))
_ENDIF
 20   continue
      return
      end
**==upack.f
_IFN(bits8)
      subroutine upack(iconf,nint,iorb)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer z3,iconf,joker
_ELSEIF(convex)
      integer *8 z3,iconf,joker
_ELSEIF(i8drct)
      integer *8 z3,iconf,joker,shift
      external shift
_ELSE
      REAL iconf,joker
_ENDIF
c     unpack occupation scheme in 'conf' into 'iorb('nint')'
c     conf has 2 bits/orbital ; 1 bit/electron
c     so room for 32 orbitals / word
      dimension i0102(4),iconf(*),iorb(*)

_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
      data i0102/0,1,0,2/
_IF(cray,t3e,convex,i8,i8drct)
      data z3/3/
_ELSE
      z3 = pad(3)
_ENDIF
      mm = min(32,nint)
      joker = iconf(1)
      ii = 1
 20   do 30 i = ii , mm
_IF(convex,i8)
       joker = ishftc(joker,2,64)
_ELSEIF(i8drct)
       joker = shift(joker,2)
_ELSE
       joker = SHIFT(joker,2)
_ENDIF
_IF(cray,t3e,convex,i8,i8drct)
        iorb(i) = i0102(iand(joker,z3)+1)
_ELSE
        iorb(i) = i0102(ipad(dand(joker,z3))+1)
_ENDIF
 30   continue
      if (mm.eq.nint) return
      mm = nint
      joker = iconf(2)
      ii = 33
      go to 20
      end
_ELSE
      subroutine upack(iconf,nint,iorb)
      implicit REAL  (a-h,o-z),integer (i-n)
      REAL  iconf,joker
c     unpack occupation scheme in 'conf' into 'iorb('nint')'
c     conf has 2 bits/orbital ; 1 bit/electron
c     so room for 32 orbitals / word
      dimension i0102(4),iconf(*),iorb(*)
c
      data m1/1/
      data i0102/0,1,0,2/
c
      data z3/3/
      call pad8(3,z3)
      mm = min(32,nint)
c***** the statement below fails with iconf = FFF3C00000000000
c***** under linux/g77 V5.0.24
c**** this gets reset to FFFBC00000000000 by the statement:
c     joker = iconf(1)
c     dcopy does not fail (at least the lastest ASCI/blas version)
      call dcopy(m1,iconf(1),1,joker,1)
      ii = 1
 20   do 30 i = ii , mm
        call shift8(joker,2,joker)
        call dand8(joker,z3,tmp8)
        itmp = ipad(tmp8)+1
        iorb(i) = i0102(itmp)
 30   continue
      if (mm.eq.nint) return
      mm = nint
      call dcopy(m1,iconf(2),1,joker,1)
      ii = 33
      go to 20
      end
_ENDIF(bits8)
**==irrcon.f
      function irrcon(conf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,z1
_ELSEIF(convex)
      integer *8 conf,z1
_ELSEIF(i8drct))
      integer *8 conf,z1,shift
      external shift
_ELSE
      REAL conf,z1
      logical eq
_ENDIF
      dimension conf(*)
c     determine irreducible representation for the configuration
c     given as occupation pattern in conf (nint orbitals)
c     iro => representations of the orbitals
INCLUDE(common/sizes)
INCLUDE(common/cdcryi)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/table)
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z1/1/
_ELSE
_IF(bits8)
      call pad8(1,z1)
_ELSE
      z1 = pad(1)
_ENDIF
_ENDIF
      irrcon = jfrep
      iw = 1
      is = 1
      do 30 i = 1 , nint
        do 20 j = 1 , 2
_IF1() PAUL Try ishftc version here on cray
_IF(cray,t3e)
          if (iand(shift(conf(iw),is),z1).eq.z1)
_ELSEIF(convex,i8)
          if (iand(ishftc(conf(iw),is,64),z1).eq.z1)
_ELSEIF(i8drct)
          if (iand(shift(conf(iw),is),z1).eq.z1)
_ELSEIF(bits8)
          call shift8(conf(iw),is,tmp8)
          call dand8(tmp8,z1,tmp88)
          if (eq(tmp88,z1))
_ELSE
          if (eq(dand(SHIFT(conf(iw),is),z1),z1))
_ENDIF
     +        irrcon = mult(irrcon,iro(i))
          is = is + 1
 20     continue
        if (is.ge.n64) then
          iw = iw + 1
          is = 1
        end if
 30   continue
      return
      end
**==prcon.f
      subroutine prcon(conf,iseq,itype)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer z7,conf,joker
_ELSEIF(convex,i8drct)
      integer *8 z7,conf,joker
      integer *8 shiftr
      external shiftr
_ELSE
      REAL z7,conf,joker
_ENDIF
c     itype gives the type of the configuration :
c     itype = 1  : reference configuration + excitation pattern
c             2 : (n) or (n-1) state with no. of spin couplings
c             3 : (n-2) state with no. of singlets/triplets
c     itype = 0 : just the occupation pattern with sequence number
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/prints)
INCLUDE(common/natorb)
INCLUDE(common/mapp)
      dimension conf(*)
      dimension joct(22)
      common/moco/iorb(64)
INCLUDE(common/orb)
_IF(cray,t3e,convex,i8,i8drct)
      data z7/7/
_ELSEIF(bits8)
      call pad8(7,z7)
_ELSE
      z7 = pad(7)
_ENDIF
c.... print occupations in batches of 32 / first batch is special
c.... choose proper action
      if (itype.ge.2) then
c.... (n) and (n-1) states
        if (iseq.ne.iprwh) go to 30
        iprwh = iprwh + iprcon
        nn = min(32,nint)
        call upack(conf,nn,iorb)
        if (itype.gt.2) then
c.... (n-2) states
*****
c         write(6,7879) conf(3)
c7879     format(' *** conf(3) = ', z16)
****
          call upack2(conf(3),jj,ii)
          if (.not.oprint(28)) then
           if (modei) then
             write (iwr,6040) iseq , ii , jj , (iorb(j),j=1,nn)
           else
             write (iwr,6040) iseq , ii , jj , (iorb(mapei(j)),j=1,nn)
           end if
          end if
        else if (modei) then
          if (.not.oprint(28)) then
_IF(cray,t3e,convex,i8,i8drct)
           write (iwr,6030) iseq,conf(3), (iorb(j),j=1,nn)
_ELSE
           iiii=ipad(conf(3))
           write (iwr,6030) iseq, iiii, (iorb(j),j=1,nn)
_ENDIF
          endif
        else
          if (.not.oprint(28)) then
_IF(cray,t3e,convex,i8,i8drct)
           write (iwr,6030) iseq,conf(3), (iorb(mapei(j)),j=1,nn)
_ELSE
           iiii=ipad(conf(3))
           write (iwr,6030) iseq, iiii, (iorb(mapei(j)),j=1,nn)
_ENDIF
          endif
        end if
      else if (itype.eq.1) then
        nn = min(32,nint)
        call upack(conf,nn,iorb)
c.... reference conf
_IF(i8drct)
        joker = conf(nw)
_ELSEIF(linux)
        call dcopy(1,conf(nw),1,joker,1)
_ELSE
        joker = conf(nw)
_ENDIF
        do 20 loop = 1 , 22
_IF(cray,t3e,convex,i8,i8drct)
          joct(loop)=iand(joker,z7)
_ELSEIF(bits8)
          call dand8(joker,z7,tmp8)
          joct(loop)=ipad(tmp8)
_ELSE
          joct(loop)=ipad(dand(joker,z7))
_ENDIF
_IF(convex,i8)
          joker=ishft(joker,-3)
_ELSEIF(bits8)
          call shiftr8(joker,3,joker)
_ELSE
          joker=shiftr(joker,3)
_ENDIF
 20     continue
        if(.not.oprint(28)) then
         if (modei) then
           write (iwr,6020) iseq , joct , (iorb(j),j=1,nn)
         else
           write (iwr,6020) iseq , joct , (iorb(mapei(j)),j=1,nn)
         end if
        endif
      else if (itype.eq.0) then
        nn = min(32,nint)
        call upack(conf,nn,iorb)
        if (modei) then
          if (iseq.gt.0) write (iwr,6050) iseq , (iorb(j),j=1,nn)
          if (iseq.le.0) write (iwr,6010) (iorb(j),j=1,nn)
        else
          if (iseq.gt.0) write (iwr,6050) iseq , (iorb(mapei(j)),j=1,nn)
          if (iseq.le.0) write (iwr,6010) (iorb(mapei(j)),j=1,nn)
        end if
      end if
      if (nint.gt.32) then
        nn = nint - 32
        call upack(conf(2),nn,iorb)
        if (.not.oprint(28)) then
         if (modei) then
           write (iwr,6010) (iorb(j),j=1,nn)
         else
           write (iwr,6010) (iorb(mapei(j)),j=1,nn)
         end if
        end if
      end if
 30   return
c.... formats
 6010 format (36x,32i2)
 6020 format (2x,i6,3x,22i1,3x,32i2)
 6030 format (2x,i6,5x,'(# ',i9,7x,')',3x,32i2)
 6040 format (2x,i6,5x,'(#s',i5,' / #t',i5,' )',3x,32i2)
 6050 format (25x,i6,5x,32i2)
      end
**==iann.f
      function iann(ref,i,conf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer z1,z2,z3,conf,ref
_ELSEIF(convex)
      integer *8 z1,z2,z3,conf,ref
      integer *8 imas
      integer *8 shiftr
_ELSEIF(i8drct)
      integer *8 z1,z2,z3,conf,ref
      integer *8 imas
      integer *8 shiftr,shift
      external shift,shiftr
_ELSE
      REAL z1,z2,z3,imas,conf,ref
      logical eq
_ENDIF
c     annihilate an electron in orbital i in configuration ref
c     to produce configuration conf, both using nw words
c     conf contains the occupation-scheme (2 bits per orbital)
c     left-justified , to be read from left to right ,
c     while the number of 1-bits gives the occupation (11b = 2)
c     so each word has room for  n32 orbitals
c     upon return : iann = 1 => success  / iann = 0 => failure
      dimension ref(*),conf(*)
INCLUDE(common/cdcryi)
      common/dopey/nshif
INCLUDE(common/orb)
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z1,z2,z3/1,2,3/
_ELSE
_IF(bits8)
      call pad8(1,z1)
      call pad8(2,z2)
      call pad8(3,z3)
_ELSE
      z1 = pad(1)
      z2 = pad(2)
      z3 = pad(3)
_ENDIF
_ENDIF
c.... see in which word and which position  i is
      kc = (i-1)/n32
      ks = (i-kc*n32)*2
      kc = kc + 1
c.... see if orbital was just filled
_IF1() PAUL Try ishft version here on cray
_IF(cray,t3e,i8drct)
      if(iand(shiftr(ref(nwm1+kc),n64-ks),z1).ne.z1) then
_ELSEIF(convex,i8)
      if(iand(ishft(ref(nwm1+kc),ks-n64),z1).ne.z1) then
_ELSEIF(bits8)
      call shiftr8(ref(nwm1+kc),n64-ks,tmp8)
      call dand8(tmp8,z1,tmp88)
      if(.not.eq(tmp88,z1)) then
_ELSE
      if(.not.(eq(dand(shiftr(ref(nwm1+kc),n64-ks),z1),z1))) then
_ENDIF
c.... check occupation of i
_IF1() PAUL Try ishftc version here on cray
_IF(cray,t3e)
        iocc = iand(shift(ref(kc),ks),z3)
_ELSEIF(convex,i8)
        iocc = iand(ishftc(ref(kc),ks,64),z3)
_ELSEIF(i8drct)
        iocc = iand(shift(ref(kc),ks),z3)
_ELSEIF(bits8)
        call shift8(ref(kc),ks,tmp8)
        call dand8(tmp8,z3,tmp88)
        iocc = ipad(tmp88)
_ELSE
        iocc = ipad(dand(SHIFT(ref(kc),ks),z3))
_ENDIF
        if (iocc.lt.1) go to 20
        if (iocc.eq.1) then
c.... 1 electron in i
          imas = z1
        else
c.... 2 electrons in i
          imas = z2
        end if
        iann = 1
_IF(i8drct)
        do loop =1,nw
         conf(loop) = ref(loop)
        enddo
_ELSEIF(linux)
        call dcopy(nw,ref,1,conf,1)
_ELSE
        call fmove(ref,conf,nw)
_ENDIF
_IF1() PAUL missing cray .. try first clause (previously 
_IF1() PAUL would have used third
_IF(cray,t3e,convex,i8)
        conf(nwm1+kc) = ior(conf(nwm1+kc),ishftc(z1,n64-ks+1,64))
        conf(nw)=ishft(ref(nw),-nshif)
        conf(kc) = ieor(conf(kc),ishftc(imas,n64-ks,64))
_ELSEIF(bits8)
        call shift8(z1,n64-ks+1,tmp8)
        call dor8(conf(nwm1+kc),tmp8,conf(nwm1+kc))
        call shiftr8(ref(nw),nshif,conf(nw))
        call shift8(imas,n64-ks,tmp8)
        call dxor8(conf(kc),tmp8,conf(kc))
_ELSEIF(i8drct)
        conf(nwm1+kc) = ior(conf(nwm1+kc),shift(z1,n64-ks+1))
        conf(nw)=shiftr(ref(nw),nshif)
        conf(kc) = ieor(conf(kc),shift(imas,n64-ks))
_ELSE
        conf(nwm1+kc) = dor(conf(nwm1+kc),SHIFT(z1,n64-ks+1))
        conf(nw)=shiftr(ref(nw),nshif)
        conf(kc) = dxor(conf(kc),SHIFT(imas,n64-ks))
_ENDIF
        return
      end if
c.... no electrons in i => failure
 20   iann = 0
      return
      end
**==icre.f
      function icre(ref,i,conf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer z1,z2,z3,conf,ref
_ELSEIF(convex)
      integer *8 z1,z2,z3,conf,ref
      integer *8 imas
      integer *8 shiftr
_ELSEIF(i8drct)
      integer *8 z1,z2,z3,conf,ref
      integer *8 imas
      integer *8 shiftr, shift
      external shift,shiftr
_ELSE
      REAL z1,z2,z3,imas,conf,ref
      logical eq
_ENDIF
INCLUDE(common/cdcryi)
      dimension ref(*),conf(*)
INCLUDE(common/orb)
c     create an electron in orbital i in configuration conf
c     see iann
c.... determine which word and which position
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z1,z2,z3/1,2,3/
_ELSE
_IF(bits8)
      call pad8(1,z1)
      call pad8(2,z2)
      call pad8(3,z3)
_ELSE
      z1 = pad(1)
      z2 = pad(2)
      z3 = pad(3)
_ENDIF
_ENDIF
      kc = (i-1)/n32
      ks = (i-kc*n32)*2
      kc = kc + 1
c.... check if orbital was just emptied
_IF1() PAUL Try ishft version here on cray
_IF(cray,t3e,i8drct)
      if(iand(shiftr(ref(nwm1+kc),n64-ks+1),z1).ne.z1) then
_ELSEIF(convex,i8)
      if(iand(ishft(ref(nwm1+kc),ks-n64-1),z1).ne.z1) then
_ELSEIF(bits8)
      call shiftr8(ref(nwm1+kc),n64-ks+1,tmp8)
      call dand8(tmp8,z1,tmp88)
      if(.not.eq(tmp88,z1)) then
_ELSE
      if(.not.eq(dand(shiftr(ref(nwm1+kc),n64-ks+1),z1),z1)) then
_ENDIF
c.... check occupation of i
_IF(cray,t3e)
        iocc = iand(shift(ref(kc),ks),z3)
_ELSEIF(convex,i8)
        iocc = iand(ishftc(ref(kc),ks,64),z3)
_ELSEIF(i8drct)
        iocc = iand(shift(ref(kc),ks),z3)
_ELSEIF(bits8)
        call shift8(ref(kc),ks,tmp8)
        call dand8(tmp8,z3,tmp88)
        iocc = ipad(tmp88)
_ELSE
        iocc = ipad(dand(SHIFT(ref(kc),ks),z3))
_ENDIF
        if (iocc.lt.1) then
c.... no electrons in orbital i
          imas = z1
        else if (iocc.eq.1) then
c.... 1 electron in orbital i
          imas = z2
        else
          go to 20
        end if
        icre = 1
_IF(i8drct)
        do loop=1,nw
         conf(loop) = ref(loop)
        enddo
_ELSEIF(linux)
        call dcopy(nw,ref,1,conf,1)
_ELSE
        call fmove(ref,conf,nw)
_ENDIF
_IF1() PAUL try first clause on cray
_IF(cray,t3e,convex,i8)
        conf(nwm1+kc) = ior(conf(nwm1+kc),ishftc(z1,n64-ks,64))
        conf(kc) = ieor(conf(kc),ishftc(imas,n64-ks,64))
_ELSEIF(i8drct)
        conf(nwm1+kc) = ior(conf(nwm1+kc),shift(z1,n64-ks))
        conf(kc) = ieor(conf(kc),shift(imas,n64-ks))
_ELSEIF(bits8)
        call shift8(z1,n64-ks,tmp8)
        call dor8(conf(nwm1+kc),tmp8,conf(nwm1+kc))
        call shift8(imas,n64-ks,tmp8)
        call dxor8(conf(kc),tmp8,conf(kc))
_ELSE
        conf(nwm1+kc) = dor(conf(nwm1+kc),SHIFT(z1,n64-ks))
        conf(kc) = dxor(conf(kc),SHIFT(imas,n64-ks))
_ENDIF
        return
      end if
c.... 2 electrons in i => failure
 20   icre = 0
      return
      end
**==ichh.f
      function ichh(ref,nref,conf,nword)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer ref,conf
_ELSEIF(convex,i8drct)
      integer *8 ref,conf
_ELSE
      REAL ref,conf
_ENDIF
c     compare the (nword-1) occ.-scheme words of conf to those
c     of the reference set to see if conf is really new
c     if so : ** return ichh = 1 **
c     if conf equals a reference state => combine their excitation
c     patterns and   ** return ichh = 0 **
      dimension ref(nword,*),conf(*)
INCLUDE(common/orb)
      if (nwm1.ne.1) then
        k = 0
        jr = 0
 20     ir = locate(ref(1,k+1),nref-k,nword,conf(1))
        if (ir.ne.0) then
          if (ir.eq.jr) then
c.... found identical configuration  / adjust excit. allowance
            ir = ir + k
            go to 30
          else
            k = k + ir - 1
            ir = locate(ref(2,k+1),nref-k,nword,conf(2))
            if (ir.ne.0) then
              if (ir.eq.1) then
                ir = ir + k
                go to 30
              else
                k = k + ir - 1
                jr = 1
                go to 20
              end if
            end if
          end if
        end if
      else
        ir = locate(ref(1,1),nref,nword,conf(1))
        if (ir.ne.0) go to 30
      end if
c.... no match found
      ichh = 1
      return
_IF(cray,t3e,convex,i8,i8drct)
 30   ref(nword,ir) = ior(ref(nword,ir),conf(nword))
_ELSEIF(bits8)
 30   call dor8(ref(nword,ir),conf(nword),ref(nword,ir))
_ELSE
 30   ref(nword,ir) = dor(ref(nword,ir),conf(nword))
_ENDIF
c... adjust annhilate/create mask
      do 40 k = 1 , nwm1
        jr = k + nwm1
_IF(cray,t3e,convex,i8,i8drct)
        ref(jr,ir)=iand(ref(jr,ir),conf(jr))
_ELSEIF(bits8)
        call dand8(ref(jr,ir),conf(jr),ref(jr,ir))
_ELSE
        ref(jr,ir)=dand(ref(jr,ir),conf(jr))
_ENDIF
 40   continue
      ichh = 0
      return
      end
**==itint.f
      subroutine itint(conf,nw)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,mask,masin,z0
_ELSEIF(convex,i8drct)
      integer *8 conf,mask,masin,z0
_ELSE
      REAL conf,mask,masin,z0
      logical eq
_ENDIF
      logical nruse
c     *** perform internal excitations ***
c     starting from 'nref' 'conf'igurations we produce
c     a total of nnst, each consisting of 'nw' words
c     the first (nw-1) words contain the occupation scheme
c     (left-justified).
c     conf(nw,i) indicates what excitations are allowed :
c      counting from the right each three bits give the allowed
c      external excitations for each successive internal one, i.e.
c      each bit indicates (1=yes,0=no) if a internal-external
c      excitation-pair is appreciated ,so the following pattern
c      (2,2),(2,1),(2,0),(1,2),(1,1),(1,0),(0,2),(0,1),(0,0)
c        1     0     0     0     0     0     0     0     1
c      means that only the ref (0,0) and a quadruple exc.(2,2) is
c      allowed ;  default probably is :  1011111 (all s+d)
INCLUDE(common/iofile)
INCLUDE(common/prints)
INCLUDE(common/cdcryi)
      dimension conf(nw,*)
      common/dopey/nshif
      common/maskc/mask(64)
INCLUDE(common/ccntl)
_IF(cray,t3e,convex,i8,i8drct)
      data z0/0/
_ELSEIF(bits8)
      call pad8(0,z0)
_ELSE
      z0 = pad(0)
_ENDIF
      top = cpulft(1)
      if (.not.oprint(28)) write (iwr,6010) top
c.... masin  : mask for internal excitation allowance
      masin = mask(n64-3)
      nshif = 3
      nnst = nref
 20   nruse = .false.
      irend = nnst
c.... loop over the present set of n-conf.

      do 30 ir = 1 , irend
c.... check if we want an internal excitation from this one
_IF(cray,t3e,convex,i8,i8drct)
        if (iand(conf(nw,ir),masin).ne.z0) then
_ELSEIF(bits8)
        call dand8(conf(nw,ir),masin,tmp8)
        if (.not.(eq(tmp8,z0))) then
_ELSE
        if (.not.(eq(dand(conf(nw,ir),masin),z0))) then
_ENDIF
c.... try to produce some ('nn') new configurations
          nruse = .true.
          call exint(conf(1,ir),conf(1,nnst+1),nn,nw,conf(1,1),nnst)
          nnst = nnst + nn
        end if
 30   continue
      if (nruse) go to 20
      return
 6010 format (/1x,80('=')//' **** configuration generator called at',
     +        f9.2,' secs')
      end
**==exint.f
      subroutine exint(ref,conf,nc,nw,old,nold)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,ref,tcr,old,scr
      integer z1,z7
_ELSEIF(convex,i8drct,i8drct)
      integer *8 conf,ref,tcr,old,scr
      integer *8 z1,z7
_ELSE
      REAL conf,ref,tcr,old,scr
      REAL z1,z7
      logical eq
_ENDIF
c     perform single (internal) replacements of orbitals
c     producing from the reference configuration 'ref'
c     'nc' new  'conf's , each using nw words
c     .. note : the nw'st word contains the exc. pattern ..
c     every new formed conf. is checked against the nold old
c     ones to see if it is really new
INCLUDE(common/sizes)
INCLUDE(common/corctl)
      common/orb/nint
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/hold)
      dimension ref(*),conf(nw,*),scr(5),tcr(5),old(nw,*)
_IF(cray,t3e,convex,i8,i8drct)
      data z1,z7/1,7/
_ELSEIF(bits8)
      call pad8(1,z1)
      call pad8(7,z7)
_ELSE
      z1 = pad(1)
      z7 = pad(7)
_ENDIF
      nc = 0
_IF(i8drct)
      do i=1,nw
       tcr(i) = ref(i)
      enddo
_ELSEIF(linux)
      call dcopy(nw,ref,1,tcr,1)
_ELSE
      call fmove(ref,tcr,nw)
_ENDIF

_IF(cray,t3e,convex,i8,i8drct)
      ref(nw) = iand(ref(nw),z7)
_ELSEIF(bits8)
      call dand8(ref(nw),z7,ref(nw))
_ELSE
      ref(nw) = dand(ref(nw),z7)
_ENDIF
      irrx = irrcon(tcr)
      do 30 j = 1 , nint
c.... produce a (n+1) state
        if (icre(tcr,j,scr).ne.0) then
          irry = mult(iro(j),irrx)
          do 20 i = 1 , nint
            if (i.ne.j) then
c.... annihilate an electron => n-state
              call corfai('congen')
              if (iann(scr,i,conf(1,nc+1)).ne.0) then
_IF(cray,t3e,convex,i8,i8drct)
                if(conf(nw,nc+1).eq.z1) then
_ELSE
                if (eq(conf(nw,nc+1),z1)) then
_ENDIF
c... generated conf is a pure vacuum state---no further excitations
c... are allowed ------ therefore see if it is of appropriate symmetry
c... or obeys any occupation number restrictions in force
                  if (mult(iro(i),irry).ne.1) go to 20
                  if (ihold(conf(1,nc+1),mestri(1,1,1),minstr(1,1),
     +                maxstr(1,1),nesvac).eq.0) go to 20
                end if
c.... check if we had this one before
                if (ichh(old,nold,conf(1,nc+1),nw).ne.0) then
                  nc = nc + 1
                  nuse = nuse + nw
                end if
              end if
            end if
 20       continue
        end if
 30   continue
      return
      end
**==ext.f
      subroutine ext(conf,nw)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,z0,z2,z6
_ELSEIF(convex,i8drct)
      integer *8 conf,z0,z2,z6
_ELSE
      REAL conf,z0,z2,z6
      logical eq
_ENDIF
c     produce (n-1) and (n-2) states from the 'nnst' reference
c     configurations, produced by int (see there).
      dimension conf(nw,*)
      common/dopey/nshif,nlev,nes
INCLUDE(common/ccntl)
INCLUDE(common/hold)
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z2,z6/0,2,6/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(2,z2)
      call pad8(6,z6)
_ELSE
      z0 = pad(0)
      z2 = pad(2)
      z6 = pad(6)
_ENDIF
      nes = nesuno
      nlev = 2
      nshif = 1
c.... (n-1) states
c.... loop over the 'nnst' (n)-states
      do 20 ir = 1 , nnst
_IF(cray,t3e,convex,i8,i8drct)
        if (iand(conf(nw,ir),z6).ne.z0) then
_ELSEIF(bits8)
        call dand8(conf(nw,ir),z6,tmp8)
        if (.not.(eq(tmp8,z0))) then
_ELSE
        if (.not.(eq(dand(conf(nw,ir),z6),z0))) then
_ENDIF
          call exext(conf(1,ir),conf(1,nnst+nmin1+1),nn,nw,
     +               conf(1,nnst+1),nmin1)
          nmin1 = nmin1 + nn
        end if
 20   continue
c... (n-2) states
      if (nmin1.le.0) return
      kkkk = nnst + nmin1
      nes = nesdue
      nlev = 3
c.... loop over nmin1 (n-1)-states
      do 30 ir = 1 , nmin1
_IF(cray,t3e,convex,i8,i8drct)
        if (iand(conf(nw,nnst+ir),z2).ne.0) then
_ELSEIF(bits8)
        call dand8(conf(nw,nnst+ir),z2,tmp8)
        if (.not.(eq(tmp8,z0))) then
_ELSE
        if (.not.(eq(dand(conf(nw,nnst+ir),z2),z0))) then
_ENDIF
          call exext(conf(1,nnst+ir),conf(1,kkkk+nmin2+1),nn,nw,
     +               conf(1,kkkk+1),nmin2)
          nmin2 = nmin2 + nn
        end if
 30   continue
c.... end of excitation process
      return
      end
**==exext.f
      subroutine exext(ref,conf,nc,nw,old,nold)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,ref,old,z1
_ELSEIF(convex,i8drct)
      integer *8 conf,ref,old,z1
_ELSE
      REAL conf,ref,old,z1
      logical eq
_ENDIF
c     produce a set of 'n-1' states 'conf' from the 'n' state
c     'ref' by annihilating an electron in the internal space
c     the number of 'conf's produced is 'nc', each using 'nw'
c     words / the nw'st word contains the excitation scheme, see
c     routines  int/exint.
c     compare the new found conf's to the 'nold' 'old' ones
      dimension ref(*),conf(nw,*),old(nw,*)
       common/orb/nint
INCLUDE(common/corctl)
       common/dopey/nshif,nlev,nes
INCLUDE(common/hold)
_IF(cray,t3e,convex,i8,i8drct)
      data z1/1/
_ELSE
_IF(bits8)
      call pad8(1,z1)
_ELSE
      z1 = pad(1)
_ENDIF
_ENDIF
      nc = 0
      do 20 i = 1 , nint
        call corfai('congen')
c.... annihilate an electron => n-1 state
        if (iann(ref,i,conf(1,nc+1)).ne.0) then
_IF(cray,t3e,convex,i8,i8drct)
          if(conf(nw,nc+1).eq.z1) then
_ELSE
          if (eq(conf(nw,nc+1),z1)) then
_ENDIF
c... no further excitations allowed---screen for occ number restrictions
            if (ihold(conf(1,nc+1),mestri(1,1,nlev),minstr(1,nlev),
     +          maxstr(1,nlev),nes).eq.0) go to 20
          end if
c.... check if we had this one before
          if (ichh(old,nold,conf(1,nc+1),nw).ne.0) then
            nuse = nuse + nw
            nc = nc + 1
          end if
        end if
 20   continue
      return
      end
**==rdsym.f
      subroutine rdsym
c      initializes the multiplication tables
c     for group, whose symbol it recognizes
c**** only d2h and subgroups thereof ****
c.... symr  : group-name, dimension
c.... name  : names of irr. rep.
c.... mst   : group-multiplication table  in triangular form
c
      implicit REAL  (a-h,o-z),integer (i-n)
      character *8 isymr,ites
      character *8 name
      dimension isymr(8),isyms(8),name(27)
INCLUDE(common/sizes)
INCLUDE(common/symci)
INCLUDE(common/symchk)
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/cic)
INCLUDE(common/orb)
c.... nst is the number of entries in the symr-table
      data nst/8/
c-----------------------------------------------------------------
c c1
      data isymr(1),isyms(1)   /'c1'  ,  1 /
      data  name(1)           /'a'   /
c cs,ci,c2
      data isymr(2),isyms(2)   /'cs'  ,2/
      data (name(i),i=2,3  )  /'a'''  ,'a''''' /
      data isymr(3),isyms(3)   /'ci'  ,2/
      data (name(i),i=4,5  )  /'ag'  ,'au'  /
      data isymr(4),isyms(4)   /'c2'  ,2/
      data (name(i),i=6,7  )  /'a'   ,'b'   /
c d2,c2v,c2h
      data isymr(5),isyms(5)   /'d2'  ,4/
      data (name(i),i=8,11 )  /'a1' ,'b1'  ,'b2'  ,'b3'  /
      data isymr(6),isyms(6)   /'c2v' ,4/
      data (name(i),i=12,15)  /'a1','a2'  ,'b1'  ,'b2'  /
      data isymr(7),isyms(7)   /'c2h' ,4/
      data (name(i),i=16,19)  /'ag'  ,'bg'  ,'au'  ,'bu'  /
c d2h
      data isymr(8),isyms(8)   /'d2h' ,8/
      data (name(i),i=20,27)  /'ag'  ,'b1g' ,'b2g' ,'b3g' ,
     #                         'au'  ,'b1u' ,'b2u' ,'b3u' /
c-----------------------------------------------------------------
c
c
      if (ninnex.le.0) call caserr('invalid number of active orbitals')
      call rchar(sname)
      mbas = 0
      do 20 i = 1 , nst
        if (sname.eq.isymr(i)) then
           nirr = isyms(i)
           go to 30
        end if
        mbas = mbas + isyms(i)
 20   continue
c
c...  in case we have a symmetry one or symmetry 1 directive
c
      do i=1,8
         if (sname.eq.anumt(i).or.sname.eq.bnumt(i)) then
            ifrep = i
            return
         end if
      end do
c
      call caserr('invalid group name in symmetry directive')
c....
c.... symbol recognised
c....
 30   lsymd = .false.
      lsadap = .false.
      do 40 i = 1 , nirr
        rname(i) = name(mbas+i)
 40   continue
      call inpa(ites)
      ifrep = locatc(rname,nirr,ites)
      nintt = 0
      nextt = 0
      do 50 nsym = 1 , nirr
        call input
        call inpi(kkk)
        nintt = nintt + kkk
        call inpi(lll)
        nextt = nextt + lll
        if (nintt.gt.nint .or. nextt.gt.next .or. kkk.lt.0 .or. 
     +      lll.lt.0) go to 60
        call setsto(kkk,nsym,iro(nintt-kkk+1))
        call setsto(lll,nsym,iro(nint+nextt-lll+1))
        norbi(nsym) = kkk
        norbe(nsym) = lll
 50   continue
      if (nintt.eq.nint .and. nextt.eq.next) go to 70
 60   call caserr(
     + 'invalid number of orbitals specified in symmetry directive')
 70   return
      end
**==cidata.f
      subroutine cidata(conf,odci)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer prorb,conf,z3
_ELSEIF(convex)
      integer *8 prorb,conf,z3
_ELSEIF(i8drct)
      integer *8 prorb,conf,z3
      integer *8 shift
      external shift
_ELSE
      REAL prorb,conf,z3
_ENDIF
      logical odci, lext 
      character *4 itext
      character *8 csftyp,to,print
      character *8 frodo,bilbo
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
INCLUDE(common/cdcryi)
      common/dopey/nshif
c     **** input routine for scep-ci  symbolic part ****
      dimension csftyp(4)
      dimension conf(*)
INCLUDE(common/diactl)
      common/moco/iprorb(64),prorb(64)
INCLUDE(common/hold)
      dimension nestri(4)
      equivalence (nestri(1),nesvac)
c
INCLUDE(common/direc)
INCLUDE(common/machin)
INCLUDE(common/infoa)
c     common/data1/vlist(400),newpro(206),louta(2,maxorb),norbt,
c    * norbta(maxorb),norbtp,norbc
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
      logical omulla,omullo,omullg,omullp,omullu,oactrn
c
INCLUDE(common/scfwfn)
INCLUDE(common/restar)
INCLUDE(common/mapp)
INCLUDE(common/natorb)
      common/expanc/lext(nd200)
INCLUDE(common/dir_erg)
INCLUDE(common/symci)
INCLUDE(common/helpr)
INCLUDE(common/work)
INCLUDE(common/table)
INCLUDE(common/cic)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
INCLUDE(common/corctl)
INCLUDE(common/symchk)
INCLUDE(common/iofile)
INCLUDE(common/trial)
INCLUDE(common/modeci)
INCLUDE(common/cjdavid)
INCLUDE(common/ccepa)
_IF(ga)
      logical opcivect
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
      common/remc_ed19/nini,ipip,ioff,ibasp,dfrac_core,ncoretot
_ENDIF
      data ii/0/
      data csftyp/'vac','n-1','n-2','ref'/
      data to,print/'to','print'/
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z3/3/
_ELSE
_IF(bits8)
      call pad8(3,z3)
_ELSE
      z3 = pad(3)
_ENDIF
_ENDIF
      top = cpulft(1)
_IF(ga)
      opcivect = .false.
      dfrac_core = 0.75d0
_ENDIF
      ifjajd = -1
      do loop = 1 , nd200
        iky(loop) = ii
        ii = ii + loop
        mapie(loop) = loop
        mapei(loop) = loop
      enddo
      write (iwr,6010) top
      mspin = mul
      if (jump.eq.1 .or. odci) then
        nint = ne - nopen - npair - npair
        nint = nint/2 - norbc
        if (nseto.ne.0) then
          do i = 1 , nseto
            nint = nint + no(i)
          enddo
        end if
        nint = nint + npair + npair
        next = norbt - nint 
        nelec = ne - norbc - norbc
c...     adjust next (cf norbt adjust in master.m 
c...     which unfortunately comes after this
        if (.not.(oactrn)) then
           next = next - norbc
        else
c...       adjust active (take out cores)
           do i=1,next
              do j=1,norbc
                 if (norbta(i).eq.norbca(j)) then
                    next = next - 1
                 end if
              end do
           end do
         end if
        go to 40
      end if
      call inpi(nelec)
      call inpi(nint)
      call inpi(next)
 40   continue
c     add default settings for natorb - spinfree NOS to
c     section 11 + print, spin nos to section 12.
c     if mspin is set ; correct ...
       ispacg = 11
       if (mspin.ne.1) then
        isping = 12
        ispind = .false.
       else
        isping = 0
       endif
       latorb = .false.
      write (iwr,6020) nelec , nint , next
      ninnex = nint + next
      if (nint.le.0 .or. nint.gt.62)
     +     call caserr('invalid number of internal orbitals')
      if (next.lt.0 .or. ninnex.gt.nd200)
     +     call caserr('invalid number of active orbitals')
      if (nelec.le.0 .or. nelec.gt.(nint+nint))
     +     call caserr('invalid number of electrons')
      nintwo = nint + nint
      nint21 = nintwo + 1
      nintp1 = nint + 1
      norbi(1) = nint
      norbe(1) = next
c.... determine number of words/conf
c...   now fixed in gamess to 5
      nwm1 = (nint-1)/n32 + 1
c****   extra words for keeping track of excitation
      nwp1 = nwm1 + 1
c     nw=nwm1+nwp1
      nw = 5
      nuse = nw
      if (odci) go to 460
 50   call input
c.... see what password is on it
 60   call pass(ii,itext)
      if (ii.le.0) then
c...   pass not recognised / leave to gamess
        jrec = jrec - 1
        ii = locatc(ydd(101),limit2,itext)
        if (ii.ne.0) go to 460
        call caserr(
     +    'unrecognised directive or invalid directive sequence')
      end if
c.... go to proper place
 70   go to (100,110,120,140,150,160,180,210,220,230,250,260,270,50,320,
     +       340,350,360,390,400,410,420,430,440,440,450,60,90,90,490,
     +        80,81,470,480,485,486) , ii
c
c.... =forcediag=
c
 81   call inpa(bilbo)
      if(bilbo.eq.' ')go to 50
      if(bilbo.eq.'instore') in_store = .true.
      if(bilbo.eq.'outstore') out_store = .true.
      if(bilbo.eq.'2x2') then
c     force emin
       moddia = -1
       lcall = 2
       go to 81
      endif
      if(bilbo.eq.'david1') odiag_1 = .true.
      if(bilbo.eq.'david2') odiag_2 = .true.
      if(bilbo.eq.'david3') odiag_3 = .true.
      go to 81
c.... =empty= 
 80   call caserr('was cepaz')
      go to 50
c.... =printvar= or =varprint=
 90   lvar = .true.
      go to 50
c.... symmetry
 100  call rdsym
      go to 50
c.... mrdci / read buenker info
 110  call buin
      go to 50
c.... conf   .. read reference configurations
 120  if (ninnex.le.0) call caserr(
     +    'invalid number of orbitals .. faulty directive ordering')
c
      call inpa4(itext)
      if (itext.eq.'file'.or.itext.eq.'card'.or.itext.eq.'asci') then
         call confile(conf,excit,nelec,nw,nref,nuse)
         go to 50
      end if  
c
      call input
 130  nrnw = nuse - nw
      call corfai('CI input')
      call cpack(conf(nrnw+1))
c.... add excitation-pattern
_IF(i8drct)
      conf(nuse) = excit
_ELSEIF(linux)
      call dcopy(1,excit,1,conf(nuse),1)
_ELSE
      conf(nuse) = excit
_ENDIF
c       call prcon(conf(nrnw+1),1,1)
      if (ichh(conf,nref,conf(nrnw+1),nw).ne.0) then
        nref = nref + 1
        nuse = nuse + nw
      end if
      call pass(ii,itext)
      if (ii.lt.0) then
        call errout(istrt(jrec))
        go to 400
      else if (ii.eq.0) then
        jrec = jrec - 1
        ii = locatc(ydd(101),limit2,itext)
        if (ii.eq.0) go to 130
        go to 460
      else
        go to 70
      end if
c.... excit   .. allowed excitation pattern
 140  call rexc(excit)
      go to 60
c.... spin   .. read spin-multiplicity
 150  call rspin(mspin)
      if (mspin.eq.1) ispind = .true.
      go to 50
c    emin
 160  moddia = -1
 170  call inpi(lcall)
      if (lcall.lt.2 .or. lcall.gt.50) lcall = 50
      go to 50
c... reorder
 180  if (ninnex.le.0) call caserr(
     +          'invalid number of active orbitals at reorder')
      call setsto(nd200,nd200,mapei)
_IF(cray)
      call setsto(nd200,.false.,lext)
_ELSE
      call setstl(nd200,.false.,lext)
_ENDIF
      call input
      mcount = 0
 190  call pass(ii,itext)
      if (ii.ne.0) then
        if (mcount.ne.ninnex)
     +       call caserr('invalid number or orbitals in reorder data')
        write (iwr,6030)
        write (iwr,6040) (mapie(loop),loop,loop=1,ninnex)
        modei = .false.
        go to 70
      else
        jrec = jrec - 1
        call inpi(iii)
        call rchar(bilbo)
        if (bilbo.eq.to) then
          call rinpi(jjj)
        else
          jjj = iii
          jrec = jrec - 1
        end if
        if (iii.le.0 .or. jjj.gt.nd200 .or. iii.gt.jjj .or. 
     +      (mcount+jjj-iii).ge.ninnex)
     +       call caserr(' invalid orbital specified in reorder data')
        do 200 loop = iii , jjj
          if (lext(loop)) 
     +     call caserr('repeated orbital in reorder data')
          mcount = mcount + 1
          mapie(mcount) = loop
          mapei(loop) = mcount
          lext(loop) = .true.
 200    continue
        go to 190
      end if
c...   shift
 210  call inpf(eshift)
      go to 50
c...  =lock=
 220  moddia = 0
      go to 170
c...  =trial=
 230  call inpa(frodo)
      if (frodo.eq.'diag') then
c.... trial diag 'select'/'first' 'ref'/'vac'/nselvc 'state' istate 'print'
 235     call inpa(frodo)
         if (frodo.eq.'first') then
            iselvc = 1
         else if (frodo.eq.'select') then
            iselvc = 2
         else if (frodo.eq.'state') then
            call inpi(istvc)
            istvc = max(istvc,1)
         else if (frodo.eq.'print') then
            iprvc = iprvc + 10
         else if (frodo.eq.'ref') then
            iselvc = -1
         else if (frodo.eq.'vac') then
            iselvc = -2
         else if (frodo.eq.' ') then
            go to 237
         else
            jrec = jrec - 1
            call inpi(nselvc)
            nselvc = max(nselvc,1)
         end if
         go to 235
c...  trial diag print
237      frodo = 'selected'
         if (iselvc.eq.1) frodo = '  first '
         write(iwr,6051) frodo
 6051    format (//4x,'trial ci vector by diagonalising ',a8,' states')
         if (iselvc.eq.-1) then
             write(iwr,6052)
 6052        format(4x,'taken from the reference space')
         else if (iselvc.eq.-2) then
             write(iwr,6053)
 6053        format(4x,'taken from the vacuum space')
         else 
             write(iwr,6054) nselvc
 6054        format(4x,i5,' states will be selected ')
         end if
         frodo = 'no print'
         if (iprvc.gt.0) frodo = '  print '
         if (iprvc.gt.10) frodo = 'ppprint '
         write(iwr,6055) istvc,frodo
 6055    format(4x,'state choosen is',i4,' and ',a8,' requested')
         go to 50
c...  end trial diag
      end if
      jrec = jrec -1
      iselvc = 0
      write (iwr,6050)
      nvec = 0
      temp = 0.0d0
 240  call input
      call pass(ii,itext)
      if (ii.ne.0) then
        temp = ddot(nvec,tvec,1,tvec,1)
        if (temp.lt.(1.0d-20)) call caserr('invalid trial data')
        temp = dsqrt(1.0d0/temp)
        call dscal(nvec,temp,tvec,1)
        go to 70
      else
        nvec = nvec + 1
        if (nvec.gt.mxvec) call caserr(
     +            'invalid number of terms in trial directive')
        jrec = 0
        call inpi(ivec(nvec))
        call inpf(tvec(nvec))
        write (iwr,6060) ivec(nvec) , tvec(nvec)
        go to 240
      end if
c...  =maxcyc=
 250  call inpi(maxcyc)
      go to 50
c...  =vmin=
 260  moddia = 1
      go to 170
c... =restrict=
 270  if (ninnex.le.0) call caserr(
     +         'invalid number of active orbitals at restrict')
      write (iwr,6070)
 280  call input
      call pass(ii,itext)
      if (ii.ne.0) go to 70
      jrec = 0
      call inpa(frodo)
      ics = locatc(csftyp,4,frodo)
      if (ics.eq.0) then
        jrec = 0
        call inpi(ics)
        if (ics.lt.0 .or. ics.gt.3) call caserr(
     +      'illegal csftype parameter in RESTRICT directive')
        ics = ics + 1
      end if
      nest = nestri(ics) + 1
      if (nest.gt.mxstri) call caserr(
     +        'maximum of 10 restrictions in restrict exceeded')
      nestri(ics) = nest
      call inpi(minstr(nest,ics))
      call inpi(maxstr(nest,ics))
      narsil = 0
 290  call rinpi(l)
      if (l.ne.0) then
        m = l
        call rchar(frodo)
        if (frodo.ne.to) then
          jrec = jrec - 1
        else
          call rinpi(m)
        end if
        if (m.lt.l .or. l.lt.1 .or. m.gt.nint)
     +       call caserr('illegal mo label in RESTRICT directive')
        do 300 loop = l , m
          narsil = narsil + 1
          iprorb(narsil) = loop
          kc = (loop-1)/n32
          ks = (loop-kc*n32)*2
_IF1() PAUL try first clause on cray
_IF(cray,t3e,convex,i8)
          mestri(kc+1,nest,ics) = 
     *   ior(mestri(kc+1,nest,ics),ishftc(z3,n64-ks,64))
_ELSEIF(i8drct)
          mestri(kc+1,nest,ics) = 
     *   ior(mestri(kc+1,nest,ics),shift(z3,n64-ks))
_ELSEIF(bits8)
          call shift8(z3,n64-ks,tmp8)
          call dor8(mestri(kc+1,nest,ics),tmp8,
     +              mestri(kc+1,nest,ics))
_ELSE
          mestri(kc+1,nest,ics) = 
     *   dor(mestri(kc+1,nest,ics),SHIFT(z3,n64-ks))
_ENDIF
 300    continue
        go to 290
      else
        do 310 loop = 1 , narsil
_IF(cray,t3e,convex,i8,i8drct)
          prorb(loop)=iprorb(loop)
_ELSE
_IF(bits8)
          call pad8(iprorb(loop),prorb(loop))
_ELSE
          prorb(loop)=pad(iprorb(loop))
_ENDIF
_ENDIF
 310    continue
        if (narsil.gt.31) then
          write (iwr,6080) csftyp(ics) , minstr(nest,ics) , 
     +                     maxstr(nest,ics) , (iprorb(loop),loop=1,31)
          write (iwr,6090) (iprorb(loop),loop=32,narsil)
        else
          write (iwr,6080) csftyp(ics) , minstr(nest,ics) , 
     +                     maxstr(nest,ics) , 
     +                     (iprorb(loop),loop=1,narsil)
        end if
        go to 280
      end if
c...   =vprint=
 320  call inpi(mprin)
      call inpf(gprin)
      write (iwr,6100) mprin , gprin
      do 330 i = 1 , 3
        call inpa(frodo)
        j = locatc(csftyp,3,frodo)
        if (j.ne.0) then
          vdstpr(j) = 5.0d-9
          write (iwr,6110) frodo
        end if
 330  continue
      go to 50
c...   =thresh=
 340  call inpf(temp)
      call inpi(j)
      thresh = temp*(0.1d0**j)
      go to 50
c...  =screen=
 350  screan = .false.
      write (iwr,6120)
      go to 50
c... =refgen=
 360  if (nref.eq.0) call caserr(
     +            'REFGEN appears before first conf directive')
      nvv = nref
      nshif = 0
      write (iwr,6130)
 370  call pass(ii,itext)
      if (ii.ne.0) then
        nref = nvv
        go to 70
      else
        jrec = jrec - 1
        ii = locatc(ydd(101),limit2,itext)
        if (ii.gt.0) then
          nref = nvv
          go to 460
        end if
        call inpi(i)
        call rinpi(j)
        write (iwr,6140) i , j
        if (i.le.0 .or. j.le.0 .or. i.gt.nint .or. j.gt.nint .or. 
     +      i.eq.j) call caserr('invalid a/c label in refgen directive')
        joop = 1
        do 380 loop = 1 , nref
          if (icre(conf(joop),j,prorb).ne.0) then
            call corfai('CI input')
            nrnw = nuse - nw
            if (iann(prorb,i,conf(nrnw+1)).ne.0) then
_IF(cray,t3e,i8)
             call setsto(nwm1,0,conf(nrnw+nwp1))
_ELSE
             call setstz(nwm1,0,conf(nrnw+nwp1))
_ENDIF
_IF(i8drct)
              conf(nuse) = excit
_ELSEIF(linux)
              call dcopy(1,excit,1,conf(nuse),1)
_ELSE
              conf(nuse) = excit
_ENDIF
              if (ichh(conf,nvv,conf(nrnw+1),nw).ne.0) then
                nvv = nvv + 1
                nuse = nuse + nw
              end if
            end if
          end if
          joop = joop + nw
 380    continue
        go to 370
      end if
c...  =alternate= directive
 390  malt = .true.
      write (iwr,6150)
      go to 50
c...   =blank=   directive
 400  call inpf(temp)
      call inpi(j)
      qester = temp*(0.1d0**j)
      go to 50
c... =natorb=
 410  call inpi(ispacg)
      call inpi(isping)
      call inpa(bilbo)
      latorb = bilbo.ne.print
      go to 50
c..... =cepa=
 420  call cepin
c...   switch off jacobi davidson for cepa
       ifjajd = 0
      go to 50
c...    =prconf=
 430  call inpi(iprcon)
      call inpa(bilbo)
      if (bilbo.eq.'debug') odebug(40) = .true.
      if (iprcon.eq.0) iprcon = 1
      if (iprcon.lt.0) iprcon = 199999999
      write (iwr,6160) iprcon
      iprwh = iprcon
      go to 50
c...    =symdet=  or =detsym=     or  =symdet adapt=
 440  lsymd = .true.
      lsadap = .false.
      if (jrec.lt.jump) then
        call inpa(bilbo)
        if (bilbo.eq.'adapt') then
          lsymd = .false.
          lsadap = .true.
        else
          jrec = jrec - 1
          call inpf(crtsym)
        end if
      end if
      call inpa(bilbo)
      ifrep = locatc(anumt,8,bilbo)
      if (ifrep.eq.0) ifrep = locatc(bnumt,8,bilbo) 
      go to 50
c...  =spin-density=
 450  ispind = .false.
      write (iwr,6170)
      go to 50
c...  =jacdav= 'off' 'on' 'shift' shift 'thresh' thresh 
c...      ...  'maxcyc' maxcyc 'print' 'select' maxsjd
c... (maxsdj may be ref or vacuum and is rounded to whole configurations) 
 470  call inpa(bilbo)
 471  if (ifjajd.eq.-1) then
c..      shift -1000 signifies dynamic shifting
         eshijd = -1000.0d0
         iprjd = 0
         maxcjd = 100
         ifjajd = 1
         threjd = thresh*0.5d0
         crijd = 0.0d0
         maxsjd = -10
         if (bilbo.eq.'pp1436vd')  then
c...        signal a default situation (without input)
            ifjajd = 10
            go to 1000
         end if
c...       if somebody guesses this, he will crash
      end if
      if (bilbo.eq.'off') ifjajd = 0
      if (bilbo.eq.'on') ifjajd = 1
      if (bilbo.eq.'print') then
          iprjd = iprjd + 10
      else if (bilbo.eq.'shift') then
          call inpa(frodo)
          if (frodo.eq.'dynamic') then
             eshijd = -1000.0d0
          else
             jrec = jrec - 1
             call inpf(eshijd)
          end if
      else if (bilbo.eq.'thresh') then
          call inpf(threjd)
      else if (bilbo.eq.'crit') then
          call inpf(crijd)
      else if (bilbo.eq.'select') then
          call inpa(frodo)
          if (frodo(1:3).eq.'ref') then
             maxsjd = -1
          else if (frodo(1:3).eq.'vac') then
             maxsjd = -2
          else
             jrec = jrec - 1
             call inpi(maxsjd)
          end if
c...   crijd not used in direct CI
      else if (bilbo.eq.'maxcyc') then
          call inpi(maxcjd)
      else if (bilbo.eq.' ') then
          go to 50
      end if
      go to 470
c
c...  =casgen= 
c...  real processing in direct-call
c
480   call casgcc('input',conf,idum,idum,iwr)
      go to 50
c
c...  =pciv==
c...  parallel option civector parallel in iabc (default integrals)
c
485   continue
_IF(ga)
      opcivect = .true.
_ENDIF
      go to 50
c
c...  =glob==
c...  parallel option memory of total core for global array
c
486   call inpf(dfrac_core)
      go to 50
c
c...   =mp=
c...   multi reference  mp2/3
490    call mpin(iwr)
       go to 50
c
c======================================================================
c
 460  if (irest.gt.8) mrstor = 1
c
      bilbo = 'pp1436vd'
      if (ifjajd.eq.-1) go to 471
c....  471 starts a set of post-direct-input initialisations, see there
c....  to add another one // They finish by jumping to 1000
1000  continue
c
      return
 6010 format (//1x,80('=')//' **** direct-ci input processor called at'
     +        ,f9.2,' secs')
 6020 format (/' # electrons=',i6,//' # internal orbitals=',i6,
     +        //' # external orbitals=',i6)
 6030 format (/5x,'reordered orbital sequence'/5x,26('-')
     +        //4('     function  e/i label'))
 6040 format (4(i13,i11))
 6050 format (//4x,'trial ci vector'/4x,15('*')/4x,'csf coefficient')
 6060 format (i7,f12.5)
 6070 format (//' occupation number restrictions for csfs'/1x,39('-')
     +        //'  csf   # of electrons'/
     +        ' type  minimum maximum  internal orbitals')
 6080 format (2x,a3,i9,i8,2x,31i3)
 6090 format (24x,31i3)
 6100 format (/' ci vector print control parameters=',i8,e16.4)
 6110 format (/' all ',a4,'ci coefficients will be printed')
 6120 format (/' reference conf. screened for space/spin symmetry')
 6130 format (//' reference generator'/1x,19('-')//'  annil create')
 6140 format (2i7)
 6150 format (/' level shift sign alternation requested')
 6160 format (/,' every ',i8,'th configuration is printed ')
 6170 format (/,' spin-density option requested ')
      end
      subroutine confile(conf,excit,nelec,nw,nref,nuse)
      implicit REAL  (a-h,o-z),integer (i-n)
c
c...  generate reference config set from cards file (e.g. cards casscf)
c...  default civecs.ascii
c...  e.g.        1   222220000                    .9644755897
c...  the various spin contributios have to be added
c...   *note* ndoc/nvoc are strictly local
c
_INCLUDE(common/iofile)
      parameter (maxc=2001)
      character*15 ccctxt
      common/junk/ccc(maxc)
      common/junkc/ ccctxt(maxc)
c...   maxc corresponds to dimensioning of ccctxt in casb used for punching
c...   (only coefficients gt cipuch=1.0d-10 are punched)
      dimension conf(*)
      character*4 ytext
      character*44 fname
c
      fname = 'civecs.ascii'
      crit = 1.0d-3
      ndoc = 0
      nfrz = 0
c
1     call inpa4(ytext)
      if (ytext.eq.'file') then
         call inpan(fname)
      else if (ytext.eq.'coef') then
         call inpf(crit)
         crit = crit*crit
      else if (ytext.eq.'weig') then
         call inpf(crit)
      else if (ytext.eq.'doc') then
         call inpi(ndoc)
      else if (ytext.eq.'frz') then
         call inpi(nfrz)
      else if (ytext.eq.' ') then
         go to 3
      else 
         call caserr('unrecognised key in confile')
      end if
c    
3     open(1,file=fname,form='formatted')
      rewind 1
204   format(1x,i6,3x,a15,5x,f20.10)
      nc = 0
5      nc = nc + 1
       if (nc.gt.maxc) call caserr('too many configs on ascii file')
       read(1,204,end=10) ii,ccctxt(nc),ccc(nc)
       go to 5
10    nc = nc - 1
      close(1)
c
c...  check
c
      nel = 0
      do i=1,15
         if (ccctxt(1)(i:i).eq.' ') then
            nvoc = i-1
            go to 20
         end if
         if (ccctxt(1)(i:i).eq.'2') nel = nel + 2
         if (ccctxt(1)(i:i).eq.'1') nel = nel + 1
      end do
      nvoc = 15
20    ndoc = ndoc - nfrz
      nel= nel + ndoc*2 
      if (nel.ne.nelec) then
         write(iwr,21) nel,ndoc,nelec
21       format(' # electrons in ascii including',i3,' doc: ',i3,
     #          ' # electrons in gamess "',i4)
         call caserr('wrong reference set')
      end if 
c...   gather different spin-path contributions
      do ic=1,nc
         ccc(ic) = ccc(ic)*ccc(ic)
         do jc=ic+1,nc
            if (ccctxt(jc).eq.ccctxt(ic)) then
               ccc(ic) = ccc(ic) + ccc(jc)*ccc(jc)
               ccc(jc) = 0.0d0
            end if
         end do
      end do
c...   select reference configurations
      nref0 = nref
      cccc = 0.0d0
      do ic=1,nc
         if (ccc(ic).gt.crit) then
            nref = nref + 1
            cccc = cccc + ccc(ic)
            nrnw = nuse - nw
            call corfai('confile')
            do i=1,nw
               conf(nrnw+i) = 0.0d0
            end do
            do i=1,ndoc
              if (icrer(conf(nrnw+1),i,2).eq.0) call caserr('confile1')
            end do
            do i=1,nvoc
             if (i+ndoc.gt.0) then
              if(ccctxt(ic)(i:i).eq.'1') ii=icrer(conf(nrnw+1),i+ndoc,1)
              if(ccctxt(ic)(i:i).eq.'2') ii=icrer(conf(nrnw+1),i+ndoc,2)
              if (ii.eq.0) call caserr('create confile2')
             end if
            end do
            call dcopy(1,excit,1,conf(nuse),1)
            nuse = nuse + nw
         end if
      end do
c
      write(iwr,600) fname,nc,ndoc,nvoc,crit,nref-nref0,cccc
600   format(/' **** reference configurations from ascii file ****',/,
     1       1x,'file : ',a44,' with ',i4,' states',/
     2       1x,'doc ',i4,' voc ',i3,' minimum weight ',1pe10.3,/,
     3       1x,'selected ',i4,' configurations total weight ',0pf10.6)
c
      end
**==pass.f
      subroutine pass(ip,iword)
      implicit REAL  (a-h,o-z),integer (i-n)
      character *8 cnt,word
      character *4 iword,ytrunc
c     get data field of a card and see what it is
c     ip =  0  : not recognised
c     ip >  0  : it was item ip of the list
c.... if the card is empty, take next until something's found
      parameter (ncnt=36)
      dimension cnt(ncnt)
      data cnt/'symmetry','mrdci','conf','excit','spin','emin',
     +         'reorder', 'shift','lock','trial','maxcyc','vmin',
     +         'restrict','end','vprint','thresh','screen','refgen',
     +         'alternat','blank','natorb','cepa','prconf','symdet',
     +         'detsym','density','diagmode','printvar','varprint',
     +         'mp','     ','forced' , 'jacdav','casgen','pciv',
     +         'glob'/
      call rchar(word)
      iword = ytrunc(word)
      ip = locatc(cnt,ncnt,word)
      return
      end
**==rexc.f
      subroutine rexc(excit)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer excit,sdp
_ELSEIF(convex)
      integer *8 excit,sdp,i1
_ELSEIF(i8drct)
      integer *8 excit,sdp,i1
      integer *8 shift
      external shift
_ELSE
      REAL excit,sdp
_ENDIF
      character *8 aw,octal,oct
c     read excitation pattern (or restore to default)
INCLUDE(common/work)
      dimension oct(8)
      data octal/'octal'/,oct/'0','1','2','3','4','5','6','7'/
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
c...       sdp /137b/
_IF(cray,t3e,convex,i8,i8drct)
      data sdp/95/
_ELSEIF(bits8)
      call pad8(95,sdp)
_ELSE
      sdp = pad(95)
_ENDIF
      kss = 1
      ks = 0
_IF(cray,t3e,convex,i8,i8drct)
      excit = 0
_ELSEIF(bits8)
      call pad8(0,excit)
_ELSE
      excit = pad(0)
_ENDIF
      limit = 2
      call rchar(aw)
      if (aw.eq.octal) then
        kss = 3
        limit = 8
        call rchar(aw)
      end if
 20   i = locatc(oct,limit,aw)
      if (i.ne.0) then
c.... fill in excitation allowance
_IF1() PAUL try second clause on cray
_IF(cray,t3e)
      excit = or(shift(i-1,ks),excit)
_ELSEIF(convex,i8)
      i1=i-1
      excit = ior(ishftc(i1,ks,64),excit)
_ELSEIF(i8drct)
      i1=i-1
      excit = ior(shift(i1,ks),excit)
_ELSEIF(bits8)
      i1 = i-1
      call pad8(i1,pad1)
      call shift8(pad1,ks,tmp8)
      call dor8(tmp8,excit,excit)
_ELSE
      excit = dor(SHIFT(pad(i-1),ks),excit)
_ENDIF
        ks = ks + kss
        call rchar(aw)
        go to 20
      else
c.... end of read (unknown symbol encountered)
        jrec = jrec - 1
        if (ks.eq.0) excit = sdp
        return
      end if
      end
**==rspin.f
      subroutine rspin(mspin)
c     read the spin-multiplicity 'mspin'
      implicit REAL  (a-h,o-z),integer (i-n)
      character *8 sp,word
INCLUDE(common/work)
      dimension sp(9)
      data sp/'singlet','doublet','triplet','quartet','quintupl',
     #        'sextet','septaple','octaplet','nonaplet'/
      call rchar(word)
c.... see if we recognise the word
      mspin = locatc(sp,9,word)
      if (mspin.eq.0) then
c.... not recognised : perhaps a number
        jrec = jrec - 1
        call inpi(mspin)
      end if
      return
      end
**==direct.f
      subroutine direct(core,energy)
      implicit REAL  (a-h,o-z),integer (i-n)
      logical again, otmpr
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/statis)
INCLUDE(common/disktl)
INCLUDE(common/corctl)
INCLUDE(common/ccntl)
INCLUDE(common/restri)
INCLUDE(common/file12)
INCLUDE(common/restar)
INCLUDE(common/iofile)
INCLUDE(common/timez)
INCLUDE(common/symchk)
      common/craypk/ nact(5),nacta(2,maxorb),iseca
INCLUDE(common/machin)
INCLUDE(common/comrjb)
INCLUDE(common/fpinfo)
INCLUDE(common/prints)
INCLUDE(common/dir_erg)
INCLUDE(common/mr_mp2)
_IF(ga)
      common/remc_ed19/nini,ipip,ioff,ibasp,dfrac_core,ncoretot
_ENDIF
      character*10 charwall
c
      dimension core(*)
      data nword/5/
c
c     set print flags for optimization
c
      otmpr = oprint(28)
      if (nprint.eq.-5) oprint(28) = .true.
      write (iwr,6010)
      call cpuwal(begin,ebegin)
c
      nsect2 = m6file
      do 20 i = 1 , nsect2
        notap2(i) = m6tape(i)
        iblk2(i) = m6blk(i)
        lblk2(i) = m6last(i)
 20   continue
_IFN(cray,t3e)
      call smask
_ENDIF
      call secget(isect(470),1005,iiii)
      nav = lenwrd()
      call readi(nact,mach(14)*nav,iiii,idaf)
      isec3 = iseca
      nw = nref*nword
c
c     allocate core
c     allocate all available memory
c
c use only 75% for direct ...
_IFN(ga)
      ibase = igmem_alloc_all(lword)
c     i10 = ibase + nw
_ELSE
      lword=igmem_max_memory()
      ncoretot=lword
      lword=lword*dfrac_core
      ibase=igmem_alloc(lword)
      ibasp=ibase
_ENDIF
c
      ngot = lword
c
c...   restart here if buenker is working
c
      iblkcr = mblkci
 30   nw = nref*nword
c
c     restore refs
c
      if (nw.ne.0) then
_IF(i8drct)
        call readi8(core(ibase),nw,iblkcr,numci)
_ELSE
        call rdedx(core(ibase),nw,iblkcr,numci)
_ENDIF
        if (ocifp) iblkcr = mblkci+ lensec(nw)
      endif
c
c..    in symsad some core is used for temp store of vectors
c..    this is not checked . since symsad may switch on symdet
c..    the order symsad .. sym12 should **not** be changed
c
      if (lsadap) call symsap(core(ibase+nw))
      if (lsymd) call sym12
      nreft = nref
      call endinp(core(ibase),excit)
      nw = nref*nword
_IF(i8drct)
      call wrt3i8(core(ibase),nw,iblkcr,numci)
_ELSE
      call wrt3(core(ibase),nw,iblkcr,numci)
_ENDIF
      nxblk = iposun(numci)
      nref0 = nref
      if (ncv.gt.0) nxblk = iblkcv + lensec(ncv)
      if (lsymd .or. lsadap) call inreor(core(ibase),nword,nref,iwr)
c
c...  possible casgen post processing ...
c
      call casgcc('execute',core(ibase),nreft,again,iwr)
c
      call direcg(core(ibase),core(ibase),energy)
      if (tim.le.timlim) then
c
_IF(parallel)
      again = .false.
c MPP
_ELSE
        call bucntl(core(ibase),core(ibase),again,iwr)
_ENDIF
        if (again.or.ocifp) call cizero(.false.,ngot,ocifp)
        if (ocifp) nref = nref0
        if (again) go to 30
      end if
c
      cpu = cpulft(1)
      if (mp.ne.0.and.maxcyc.le.0) then
         write (iwr,6021) cpu ,charwall()
      else
         write (iwr,6020) cpu ,charwall()
      end if
      call timana(12)
c     reset print flag
      oprint(28) = otmpr
c
c     reset core allocation
c
_IF(ga)
      ipip=ipip+ibasp-1
      call gmem_free(ipip)
_ENDIF
      call gmem_free(ibase)
c
      return
 6010 format (/1x,80('=')//40x,21('*')/40x,'direct-ci calculation'/40x,
     +        21('*'))
 6020 format (/' end of direct-ci calculation at ',f8.2,' seconds'
     +        ,a10,' wall'//1x,80('-')/)
 6021 format (/' end of (mr)-mp calculation at ',f8.2,' seconds'
     +        ,a10,' wall'//1x,80('-')/)
      end
**==casgen begin
**==casgen end
**==endinp.f
      subroutine endinp(conf,excit)
      implicit REAL  (a-h,o-z),integer (i-n)
       logical nerr
_IF(cray,t3e,i8)
      integer conf,excit
_ELSEIF(convex,i8drct)
      integer *8 conf,excit
_ELSE
      REAL conf,excit
_ENDIF
INCLUDE(common/sizes)
c     finish handling of input-data
c     fill in missing information (if possible) and check consistency
c     of the input, that was provided
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt
_ENDIF
      dimension conf(*)
INCLUDE(common/natorb)
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/cic)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
INCLUDE(common/corctl)
INCLUDE(common/file12)
INCLUDE(common/iofile)
INCLUDE(common/prints)
INCLUDE(common/discc)
      k = 1
      nerr = .false.
      if (ninnex.le.0) call caserr(
     +                  'invalid number of orbitals in ENDINP')
c... if spin not given set to minimum for nelec electrons
      if (mspin.le.0) mspin = nelec - ((nelec/2)*2) + 1
      if (mod(mspin-1,2).ne.mod(nelec,2)) call caserr(
     +    'invalid spin multiplicity detected in ci preprocessor')
c.... check if any conf. have been input, if not => init. obvious one
      if (nref.le.0) then
        nref = 1
        nuse = nuse + nw
_IF(cray,t3e,i8)
       call setsto(nwm1+nwm1,0,conf)
_ELSE
       call setstz(nwm1+nwm1,0,conf)
_ENDIF
        nsinel = mspin + 1
 20     nsinel = nsinel - 2
        nd = (nelec-nsinel)/2
        if ((nd+nsinel).gt.nint) go to 20
        if (nd.ne.0) then
          do 30 i = 1 , nd
            iflop = icre(conf,i,conf)
            iflop = icre(conf,i,conf)
 30       continue
        end if
        if (nsinel.ne.0) then
          do 40 i = 1 , nsinel
            iflop = icre(conf,nd+i,conf)
 40       continue
        end if
_IF(cray,t3e,i8)
       call setsto(nwm1,0,conf(nwm1+1))
_ELSE
       call setstz(nwm1,0,conf(nwm1+1))
_ENDIF
_IF(i8drct)
        conf(nw) = excit
_ELSEIF(linux)
        call dcopy(1,excit,1,conf(nw),1)
_ELSE
        conf(nw) = excit
_ENDIF
      end if
c... if symmetry not given set to that of first conf
      if (ifrep.lt.1) ifrep = irrcon(conf)
      jfrep = ifrep
c... eliminate unwanted reference states
      call weed(conf,nw)
      if(.not.oprint(28)) write (iwr,6010) nref
      if (nref.eq.0) then
       write(iwr,6010) nref
       call caserr(
     +            'invalid number of reference configurations')
      endif
c... check number of electrons in reference configurations
      do 60 i = 1 , nref
        nn = 0
        call prcon(conf(k),i,1)
        do 50 iw = 1 , nwm1
          nn = nn + popcnt(conf(k))
          k = k + 1
 50     continue
        if (nn.ne.nelec) then
c.... error   .. number of electrons in conf. wrong
          nerr = .true.
          write (iwr,6020) nn
        end if
        k = k + nw - nwm1
 60   continue
      if (nerr) call caserr('invalid number of electrons detected')
      if(.not.oprint(28)) then
       write (iwr,6030) sname , nirr , (i,rname(i),i=1,nirr)
       write (iwr,6040)
       write (iwr,6050) (i,rname(iro(i)),i=1,nint)
       if (next.ne.0) then
         write (iwr,6060)
         write (iwr,6050) (i,rname(iro(i)),i=nintp1,ninnex)
       end if
       if (ifrep.gt.nirr) then 
          write (iwr,6070) anumt(ifrep) , mspin
          call caserr('symmetry requested > # representations')
       end if
       write (iwr,6070) rname(ifrep) , mspin
       write (iwr,6080) (yed(notap2(loop)),loop=1,nsect2)
       write (iwr,6090) (iblk2(loop),loop=1,nsect2)
       write (iwr,6100) (lblk2(loop),loop=1,nsect2)
      endif
      mspin1 = mspin.eq.1
      return
 6010 format (///15x,'list of',i6,' reference configurations',/15x,
     +        38('-')/)
 6020 format (' # of electrons wrong in above configuration=',i6)
 6030 format (//' point group= ',a8,//,' dimension',i3,
     +        ' representations ',(5x,8(i3,1x,a6)))
 6040 format (/' symmetry-partitioning for internal orbitals',/1x,
     +        43('-')/)
 6050 format ((1x,8(2x,i3,2x,a5)))
 6060 format (/' symmetry-partitioning for external orbitals',/1x,
     +        43('-')/)
 6070 format (/' *********************************',/,
     +         ' *** state symmetry    = ',a5,' ***',/,
     +         ' *** spin multiplicity = ',i3,2x,' ***',/,
     +         ' *********************************')
 6080 format (/' 2-electron integral files'/' *************************'
     +        /12x,'lfns ',20(3x,a3))
 6090 format (' starting blocks ',20i6)
 6100 format (' terminal blocks ',20i6)
      end
**==selct.f
      subroutine selct(conf,res,nword)
c     perform selection on spin and symmetry of the
c     'conf'igurations, to produce the final list 'res'
c     'ncons' returns the number of unique configurations/symmetry
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf,res,z0,z1,pack2
_ELSEIF(convex,i8drct)
      integer *8 conf,res,z0,z1
      integer *8 pack2
_ELSEIF(bits8)
      REAL conf,res,z0,z1
_ELSE
      REAL conf,res,z0,z1,pack2
_ENDIF
INCLUDE(common/auxh)
INCLUDE(common/table)
INCLUDE(common/orb)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/corctl)
INCLUDE(common/ccepa)
      dimension conf(nword,*),res(*)
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z1/0,1/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(1,z1)
_ELSE
      z0 = pad(0)
      z1 = pad(1)
_ENDIF
_IF(i8)
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
_ENDIF
      nr = 0
c      keep track of nref  and nrfcsf for mp2/cepa
      nreft = 0
      nrfcsf = 0
      jres = 0
c.... 'final' symmetry partitioning
      call finsym
c.... select vacuum states
      if (nnst.ne.0) then
        do 20 i = 1 , nnst
          np1 = nspip(conf(1,i),1)
          if (np1.ne.0) then
c... found one
            nuse = nuse + 5
            call corfai('congen')
            res(jres+4) = z1
_IF(cray,t3e,convex,i8,i8drct)
            res(jres+3)=np1
_ELSEIF(bits8)
            call pad8(np1,res(jres+3))
_ELSE
            res(jres+3)=pad(np1)
_ENDIF
            res(jres+2) = z0
_IF(i8drct)
                do loop = 1,nwm1
                 res(jres+loop) = conf(loop,i)
                enddo
_ELSEIF(linux)
            call dcopy(nwm1,conf(1,i),1,res(jres+1),1)
_ELSE
            call fmove(conf(1,i),res(jres+1),nwm1)
_ENDIF
            jres = jres + 5
            ntotal = ntotal + np1
            nr = nr + 1
            if (i.le.nref) then
               nreft = nreft + 1
               nrfcsf = nrfcsf + np1
            end if
          end if
 20     continue
      end if
      nval = ntotal
      kkmin = nnst + 1
      kkmax = nnst + nmin1
      nnst = nr
      nref = nreft
      nstrt1 = nr + 1
c.... select doublets ( (n-1)-states )
c.... loop over symmetries
      if (nmin1.ne.0) then
        nmin1 = 0
        do 40 isym = 1 , nirr
          nr = 0
          nq1 = norbe(isym)
          do 30 i = kkmin , kkmax
_IF(cray,t3e,convex,i8,i8drct)
            if(conf(nw,i).eq.isym) then
_ELSE
            if(ipad(conf(nw,i)).eq.isym) then
_ENDIF
              np1 = nspip(conf(1,i),2)
              jumper = np1*nq1
              if (jumper.ne.0) then
c... found one
                nuse = nuse + 5
                call corfai('congen')
_IF(cray,t3e,convex,i8,i8drct)
                res(jres+4)=isym
                res(jres+3)=np1
                res(jres+2)=0
_ELSEIF(bits8)
                call pad8(isym,res(jres+4))
                call pad8(np1,res(jres+3))
                res(jres+2)=z0
_ELSE
                res(jres+4)=pad(isym)
                res(jres+3)=pad(np1)
                res(jres+2)=z0
_ENDIF
_IF(i8drct)
                do loop = 1,nwm1
                 res(jres+loop) = conf(loop,i)
                enddo
_ELSEIF(linux)
                call dcopy(nwm1,conf(1,i),1,res(jres+1),1)
_ELSE
                call fmove(conf(1,i),res(jres+1),nwm1)
_ENDIF
                jres = jres + 5
                ntotal = ntotal + jumper
                nr = nr + 1
              end if
            end if
 30       continue
          nmin1 = nmin1 + nr
          ncond(isym) = nr
 40     continue
      end if
      ndoub = ntotal - nval
      nvd = ntotal
      nlast1 = nnst + nmin1
      nstrt2 = nlast1 + 1
      if (nmin2.ne.0) then
        kkmin = kkmax + 1
        kkmax = kkmax + nmin2
        nmin2 = 0
c.... select singlets/triplets
c.... loop over symmetries
        do 60 isym = 1 , nirr
          nr = 0
          nq1 = norbsi(isym)
          nq3 = norbtr(isym)
          do 50 i = kkmin , kkmax
_IF(cray,t3e,convex,i8,i8drct)
            if (conf(nw,i).eq.isym) then
_ELSE
            if (ipad(conf(nw,i)).eq.isym) then
_ENDIF
              np1 = nspip(conf(1,i),1)
              np3 = 0
              if (nq3.ne.0) np3 = nspip(conf(1,i),3)
              npq1 = np1*nq1
              npq3 = np3*nq3
              jumper = npq1 + npq3
              if (jumper.ne.0) then
c... found one
                nuse = nuse + 5
                call corfai('congen')
                nspips(isym) = nspips(isym) + np1
                nspipt(isym) = nspipt(isym) + np3
                nsingl = nsingl + npq1
                ntripl = ntripl + npq3
_IF(bits8)
                call pack2(0,isym,res(jres+4))
                call pack2(np3,np1,res(jres+3))
                call pack2(0,0,res(jres+2))
_ELSEIF(i8drct)
                res(jres+4) = pack2(0,isym)
                res(jres+3) = pack2(np3,np1)
                res(jres+2) = pack2(0,0)
_ELSE
                res(jres+4) = pack2(0,isym)
                res(jres+3) = pack2(np3,np1)
                res(jres+2) = pack2(0,0)
_ENDIF
_IF(i8drct)
                do loop = 1,nwm1
                 res(jres+loop) = conf(loop,i)
                enddo
_ELSEIF(linux)
                call dcopy(nwm1,conf(1,i),1,res(jres+1),1)
_ELSE
                call fmove(conf(1,i),res(jres+1),nwm1)
_ENDIF
                jres = jres + 5
                ntotal = ntotal + jumper
                nr = nr + 1
              end if
            end if
 50       continue
          nmin2 = nmin2 + nr
          nconst(isym) = nr
 60     continue
      end if
      nlast2 = nlast1 + nmin2
      if (nlast2.le.0) call caserr(
     + 'no csfs in final ci list .. possible restrict problem')
      nmin12 = nmin1*nmin2
      nminv2 = nnst*nmin2
      nminv1 = nnst*nmin1
      nlst21 = nlast2 + 1
      nlastx = nlast1 + nconst(1)
      nvds = nvd + nsingl
      nsitri = nsingl + ntripl
      return
      end
**==nspip.f
      function nspip(conf,ispb)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSE
      REAL conf
_ENDIF
c     this function returns the number of spin-path's
c     for the function described by the orbital-occupation
c     in 'conf' ('nw' words).
c     2 bits/orbital : number of 1's is number of electrons
c     only spin-path's starting at 'ispb (=2s+1) ' and ending at 'mspin'
      dimension conf(*)
c.... auxilary array to help in recursion (note : ia(2) is the bottom)
c....  (i.e maximal about 80 singly occupied orbitals)
      common/moco/ia(64)
INCLUDE(common/ccntl)
c.... determine number of singly occupied orbitals
      nsin = nsior(conf)
      ndbt = iabs(ispb-mspin)
      if (mod(nsin,2).eq.mod(ndbt,2)) then
c.... very foolish : can't go from ispb to mspin with nsin singles
        if (nsin.lt.ndbt) then
        else if (nsin.eq.ndbt) then
c.... nsin = spin-difference : just one pass
          nspip = 1
          return
        else
c.... nsin > spin-difference => go along branching diagram
          call setsto(64,0,ia)
          ia(ispb+1) = 1
c.... loop over singles
          ibot = mod(ispb,2)
          do 30 is = 1 , nsin
            ibot1 = ibot + 1
            idb = max(ibot1,ispb-is,mspin-(nsin-is)) + 1
            ibot = mod(ibot+1,2)
            idt = min(ispb+is,mspin+(nsin-is)) + 1
            do 20 id = idb , idt , 2
              ia(id) = ia(id+1) + ia(id-1)
 20         continue
 30       continue
c.... collect result
          nspip = ia(mspin+1)
          return
        end if
      end if
c.... nsin < spin-difference : no spin-possibility
      nspip = 0
      return
      end
**==prsymb.f
      subroutine prsymb(iconf)
c     print final result for configuration-generator
c     conf(nw,..) : occupation-scheme of the configurations
c                   the nw'st word contains the number of triplet
c                   and singlet couplings (shift n32/ shift 0)
c     ncons       : contains the symmetry-partitioning for 'conf'
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf,mask,iplus,jplus
_ELSEIF(convex,i8drct)
      integer *8 shiftr
      integer *8 iconf,mask,iplus,jplus
      external shiftr
_ELSE
      REAL iconf,mask,iplus,jplus
_ENDIF
      character *1 leftb,righb
      character *8 ror1,ror2
INCLUDE(common/cdcryi)
      common/symcon/norb,norb1,norb16,nword,nword1,nel,mspi
INCLUDE(common/auxh)
      common/moco/ror1(32),ror2(32)
INCLUDE(common/table)
INCLUDE(common/cic)
      common/maskc/mask(64)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/natorb)
INCLUDE(common/iofile)
INCLUDE(common/cntrl)
INCLUDE(common/prints)
INCLUDE(common/sizes)
      dimension iconf(5,*)
      data leftb,righb/'(',')'/
      ibig = max(nval,ndoub)
c... print vacuum states
      if(.not.oprint(28)) write (iwr,6010) ntotal , nnst , nval
      if (nnst.gt.0) then
        do 20 i = 1 , nnst
          call prcon(iconf(1,i),i,2)
          maxsin(1) = max(maxsin(1),nsior(iconf(1,i)))
_IF(cray,t3e,convex,i8,i8drct)
          loop8=iconf(3,i)
_ELSE
          loop8 = ipad(iconf(3,i))
_ENDIF
          maxpa(4) = max(maxpa(4),loop8)
 20     continue
      end if
c.... print doublet ( (n-1) ) states
      if(.not.oprint(28)) write (iwr,6020) nmin1 , ndoub
      joop = nnst
      kk = nstrt1
      nel = nelec
      mspi = mspin
      norb = nintp1 + 1
      nword = nintp1/n32
      nmask = norb + nintp1 - nword*n64
      nword = nword + 1
_IF(convex,i8)
      iplus=ishft(mask(1),-nmask)
_ELSEIF(bits8)
      call shiftr8(mask(1),nmask,iplus)
_ELSE
      iplus=shiftr(mask(1),nmask)
_ENDIF
      nww = nint/n32
      nmask = nintp1 + nint - nww*n64
      nww = nww + 1
_IF(convex,i8)
      jplus=ishft(mask(1),-nmask)
_ELSEIF(bits8)
      call shiftr8(mask(1),nmask,jplus)
_ELSE
      jplus=shiftr(mask(1),nmask)
_ENDIF
      do 40 isym = 1 , nirr
        nj = ncond(isym)
        if (nj.ne.0) then
          if (iprcon.ne.199999999) write (iwr,6030) isym , rname(isym) , 
     +        nj , norbe(isym)
          do 30 j = 1 , nj
            joop = joop + 1
            call prcon(iconf(1,kk),joop,2)
            maxsin(2) = max(maxsin(2),nsior(iconf(1,kk)))
_IF(cray,t3e,convex,i8,i8drct)
            loop8=iconf(3,kk)
_ELSE
            loop8=ipad(iconf(3,kk))
_ENDIF
            maxpa(3) = max(maxpa(3),loop8)
_IF(cray,t3e,convex,i8,i8drct)
            iconf(nword,kk) = ior(iconf(nword,kk),iplus)
_ELSEIF(bits8)
            call dor8(iconf(nword,kk),iplus,iconf(nword,kk))
_ELSE
            iconf(nword,kk) = dor(iconf(nword,kk),iplus)
_ENDIF
            kk = kk + 1
 30       continue
        end if
 40   continue
c.... print (n-2)-states  (singlet and triplet)
      if(.not.oprint(28)) write (iwr,6040) nmin2 , nsingl , ntripl
      do 70 isym = 1 , nirr
        nj = nconst(isym)
        if (nj.ne.0) then
         if (iprcon.ne.199999999) write (iwr,6050) isym , rname(isym) , 
     +        nj
          norsi = norbsi(isym)
          nortr = norbtr(isym)
          if (iprcon.ne.199999999) write (iwr,6060) norsi , nortr
          k = 0
          do 50 i = 1 , nirr
            multi = mult(i,isym)
            if (i.le.multi) then
              k = k + 1
              ror1(k) = rname(i)
              ror2(k) = rname(multi)
            end if
 50       continue
          if (iprcon.ne.199999999) then
           if(.not.oprint(28)) then
            write (iwr,6070) (leftb,ror1(i),ror2(i),righb,i=1,k)
            write (iwr,6080)
           endif
          end if
          maxpa1 = 0
          maxpa2 = 0
          do 60 j = 1 , nj
            joop = joop + 1
            call prcon(iconf(1,kk),joop,3)
            call upack2(iconf(3,kk),lamt,lams)
            maxpa1 = max(maxpa1,lams)
            maxpa2 = max(maxpa2,lamt)
            maxsin(3) = max(maxsin(3),nsior(iconf(1,kk)))
_IF(cray,t3e,convex,i8,i8drct)
            iconf(nword,kk) = ior(iconf(nword,kk),iplus)
            iconf(nww,kk) = ior(iconf(nww,kk),jplus)
_ELSEIF(bits8)
            call dor8(iconf(nword,kk),iplus,iconf(nword,kk))
            call dor8(iconf(nww,kk),jplus,iconf(nww,kk))
_ELSE
            iconf(nword,kk) = dor(iconf(nword,kk),iplus)
            iconf(nww,kk) = dor(iconf(nww,kk),jplus)
_ENDIF
            kk = kk + 1
 60       continue
          maxpa(1) = max(maxpa(1),maxpa1)
          maxpa(2) = max(maxpa(2),maxpa2)
          mubig = maxpa1*norsi
          nubig = maxpa2*nortr
          nbig = max(nbig,mubig,nubig)
          mbig = max(mbig,mubig)
          jbig = max(jbig,mubig+nubig)
        end if
 70   continue
      kbig = max(ibig,jbig)
      lbig = max(ibig,nbig)
      maxpat = max(maxpa(1),maxpa(2),maxpa(3),maxpa(4))
      mxvvpa = maxpa(4)*maxpa(4)
      mxdvpa = maxpa(3)*maxpa(4)
      mxddpa = maxpa(3)*maxpa(3)
      max12 = maxpa(1) + maxpa(2)
      maxpsq = max12*max12
      mxv2pa = max12*maxpa(4)
      mx12pa = max12*maxpa(3)
      mxsspa = maxpa(1)*maxpa(1)
      mxttpa = maxpa(2)*maxpa(2)
      mxstpa = maxpa(1)*maxpa(2)
c
      if (maxpat**2.gt.mxcan1) then
         write(iwr,6090) maxpat,mxcan1,maxpat**2
6090     format(' the number of spincases/config is excessive ',
     1          '(',i9,')',/,' it is wise to raise mxcan1 (in sizes)',
     2          ', which is',i9,', to at least the square of this ',
     3          '(',i9,')',/' and mxcan2 to twice that')
         call caserr('mxcan problem expected in direct-CI')
      end if
      return
 6010 format (/' *** total # states=',i12//' *** vacuum states (# conf',
     +        i8,'   # states',i8,') ***',/,1x)
 6020 format (//' *** doublet states (# conf',i8,'   # states',i8,
     +        ') ***',/,1x)
 6030 format (/' sym',i3,2x,a7,' # conf',i8,'   external dimension',i8)
 6040 format (//' *** (n-2) - states (# conf',i8,'   # singlet states',
     +        i10,'   # triplet states',i10,') ***',/,1x)
 6050 format (/' sym',i3,2x,a7,' # conf',i8)
 6060 format (' singlet states external dimension',
     +        i8/' triplet states external dimension',i8)
 6070 format (1x,'external symmetries ',(t23,5(a2,a6,'-',1x,a6,a1,4x)))
 6080 format (1x)
      end
**==finsym.f
      subroutine finsym
      implicit REAL  (a-h,o-z),integer (i-n)
c     - initialise  part of common /sym/
c... generates constants to help in integral retreival
      common/stoctl/khbill,khbm1e,khbm2e,ijklno
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/orb)
INCLUDE(common/symci)
c
      ibass = 0
      ij = 0
      nik = 0
      ijkl = 0
      muj = 0
c...    above former data-statement
      jbass = nint
      do 40 ioop = 1 , nirr
        niky(ioop) = nik - ioop
        nik = nik + ioop*ioop
        ni = norbi(ioop)
        nbase(ioop) = m1e - muj
        m1e = m1e + iky(ni+1)
        mbasi(ioop) = muj
        mnbasi(ioop) = -muj - norbi(ioop) - 1
        nv = norbe(ioop)
        call ibasgn(ni,1,1,norbas(ibass+1))
        ibass = ibass + ni
        call ibasgn(nv,1,1,norbas(jbass+1))
        jbass = jbass + nv
        nvmax = max(nv,nvmax)
        ibassi(ioop,1) = norbsi(1)
        ibastr(ioop,1) = norbtr(1)
        ibassq(ioop,1) = norbsq(1)
        ikynv = iky(nv+1)
        norbsi(1) = norbsi(1) + ikynv
        norbtr(1) = norbtr(1) + ikynv - nv
        ikynv = nv*nv
        norbsq(1) = norbsq(1) + ikynv
        call ibasgn(ni,m1ia,nv,iap(muj+1))
        m1ia = nv*ni + m1ia
        muj = muj + ni
        do 30 joop = 1 , ioop
          ij = ij + 1
          if (joop.ne.ioop) then
            nij = ni*norbi(joop)
          else
            nij = iky(ni+1)
          end if
          npairi(ij) = nij
          ijsym = mult(ioop,joop)
          do 20 koop = 1 , ioop
            ijkl = ijkl + 1
            mbase(ijkl) = m2e
            loop = mult(koop,ijsym)
            if (loop.le.koop) then
              kl = iky(koop) + loop
              if (kl.ne.ij) then
                nijkl = nij*npairi(kl)
              else
                nijkl = (nij*(nij+1))/2
              end if
              m2e = m2e + nijkl
            end if
 20       continue
 30     continue
 40   continue
c... m1e=no. of fully internal 1-e ints
c... m2e no. of fully internal 2-e ints
      mubsi = iky(nvmax+1)
      mubsq = nvmax*nvmax
      ijklno = ijkl
      if (nirr.ne.1) then
        do 60 ioop = 2 , nirr
          muj = 0
          ij = 0
          do 50 joop = 1 , nirr
            koop = mult(ioop,joop)
            nv = norbe(joop)
            nij = nv*norbe(koop)
            if (koop.gt.joop) then
              ibassi(joop,ioop) = muj
              mbassi(joop,ioop) = muj - nv
              ibastr(joop,ioop) = muj
              mubsi = max(nij,mubsi)
              muj = muj + nij
            end if
            ibassq(joop,ioop) = ij
            ij = ij + nij
 50       continue
          norbsq(ioop) = ij
          norbsi(ioop) = muj
          norbtr(ioop) = muj
 60     continue
      end if
      return
      end
**==weed.f
      subroutine weed(iconf,nw)
c... weeds out unwanted reference cofigurations
c... by occupation number restrictions (if any)
c... and by space/spin symmetry  if  =screen= directive  presented
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iiii,iconf
_ELSE
      REAL iconf
_ENDIF
      dimension iconf(nw,*)
INCLUDE(common/hold)
INCLUDE(common/ccntl)
      nvv = 0
      do 30 loop = 1 , nref
        if (ihold(iconf(1,loop),mestri(1,1,4),minstr(1,4),maxstr(1,4),
     +      nesref).eq.0) go to 30
_IF(cray,t3e,convex,i8,i8drct)
        iiii = iconf(nw,loop)
_ELSE
        iiii = ipad(iconf(nw,loop))
_ENDIF
        if (iiii.lt.1) go to 30
        if (iiii.eq.1) then
          if (ihold(iconf(1,loop),mestri(1,1,1),minstr(1,1),maxstr(1,1),
     +        nesvac).eq.0) go to 30
        else if (screan) then
          go to 20
        end if
        if (irrcon(iconf(1,loop)).ne.1) go to 30
        if (nspip(iconf(1,loop),1).eq.0) go to 30
 20     nvv = nvv + 1
_IF(i8drct)
        do i=1,nw
         iconf(i,nvv) = iconf(i,loop)
        enddo
_ELSEIF(linux)
        call dcopy(nw,iconf(1,loop),1,iconf(1,nvv),1)
_ELSE
        call fmove(iconf(1,loop),iconf(1,nvv),nw)
_ENDIF
 30   continue
      nref = nvv
      return
      end
**==buin.f
      subroutine buin
c
c...  input unscrambler for mrdci (buenker) optie
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer joct,joker,z0,z7
_ELSEIF(convex,i8drct)
      integer *8 joct,joker,z0,z7
      integer *8 shiftr
      external shiftr
_ELSE
      REAL joct,joker,z0,z7
_ENDIF
      character *8 ites,itt
INCLUDE(common/comrjb)
INCLUDE(common/work)
INCLUDE(common/iofile)
      dimension itt(6),joct(22),ioct(22)
      data itt/'crit','max','weight','exc1','exc2','end'/
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z7/0,7/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(7,z7)
_ELSE
      z0 = pad(0)
      z7 = pad(7)
_ENDIF
c
      ifbuen = .true.
c
c...  see if anything more on the card
c
 20   if (jrec.ge.jump) go to 100
      call inpa(ites)
      ll = locatc(itt,6,ites)
      if (ll.le.0) go to 20
c
 30   go to (40,60,70,80,90,100) , ll
c...  crit / max. 3 fp-numbers
 40   do 50 i = 1 , 3
        if (jrec.ge.jump) go to 100
        call inpa(ites)
        ll = locatc(itt,3,ites)
        if (ll.gt.0) go to 30
        jrec = jrec - 1
        call inpf(cradd(i))
 50   continue
      go to 20
c...  max
 60   call inpi(mbuenk)
      go to 20
c...  weight
 70   call inpf(weighb)
      go to 20
c...   exc1
 80   call rexc(iex1)
      go to 20
c...   exc2
 90   call rexc(iex2)
      go to 20
c
c...  print
c
 100  write (iwr,6010) mbuenk , weighb
      do 130 i = 1 , 2
c...  unpack excitation pattern
_IF(cray,t3e,i8)
       call setsto(22,0,joct)
_ELSE
       call setstz(22,0,joct)
_ENDIF
        if (i.eq.1) then
           joker = iex2
        else
           joker = iex1
        endif
        do 110 loop = 1 , 22

_IF(cray,t3e,convex,i8,i8drct)
          joct(loop) = iand(joker,z7)
_ELSEIF(bits8)
          call dand8(joker,z7,joct(loop))
_ELSE
          joct(loop) = dand(joker,z7)
_ENDIF

_IF1() PAUL shiftr is outmoded on cray,t3e,  try ishft version
_IF(convex,i8)
          joker = ishft(joker,-3)
_ELSEIF(bits8)
          call shiftr8(joker,3,joker)
_ELSE
          joker = shiftr(joker,3)
_ENDIF
 110    continue
c...     print criterium + excit-allowance
        nn = locatz(joct,22,z0) - 1
        do 120 loop = 1 , nn
_IF(cray,t3e,convex,i8,i8drct)
          ioct(loop) = joct(loop)
_ELSE
          ioct(loop) = ipad(joct(loop))
_ENDIF
 120    continue
        ccccc = cradd(3-i)
        write (iwr,6020) ccccc , (ioct(mfg),mfg=1,nn)
 130  continue
      write (iwr,6030) cradd(3)
      write (iwr,6040)
c
      return
 6010 format (/1x,77('='),//,'  mrdci-option requested for ',i3,
     +        ' iterations  or weight gt ',f8.5)
 6020 format ('  weight >',e12.5,' excitation-pattern',22i2)
 6030 format ('  reference configurations > ',e12.5,' keep their old',
     +        ' excitation pattern')
 6040 format (/,1x,77('='))
      end
**==cepin.f
      subroutine cepin
      implicit REAL  (a-h,o-z),integer (i-n)
      character *4 ites
c
c... read control info for cepa
c
INCLUDE(common/iofile)
INCLUDE(common/ccepa)
ccepa
ccepa  .. ifcepa > 0 => cepa is on
ccepa  .. icepa  > -1 => cepa  may be switched on
ccepa  .. icepa = 0,1,2   : closed-chell cepa
ccepa  ..       =   10    : mrdcepa 
ccepa  ..       =   11    : mrcepa(1) doc-epv
ccepa  ..       =   12    : mrcepa(1) all-epv
ccepa  ..       =   20    : mrcepa1
ccepa  ..       =  100    : mrcepa00
ccepa  ..       =  101    : AQCC   
ccepa  ..       =  102    : ACPF   
crz    .. ibrst = 1 : include brillouin states
ccepa  .. icorr_mr = 0 : Ecorr  = sum pair energies (=~projection)
ccepa  .. icorr_mr = 1 : use projection to calc Ecorr
ccepa  .. icorr_mr = 2 : use variational ref. function to calc Ecorr
ccepa
c
c...  set defaults
      lcepaul = 0
      lsinsh = 1
      itcepa = 3
      critcep = 1.0d-2
      lprcep = 0
crwah
      ndoc = 0
      mcycep = 3
      crycep = 0.01d0
      ibrst = 1
      icorr_mr = 0
crwah
c
c...  read CEPA variant
c
      call inpa4(ites)
c...     cepa 0 1 2
      if (ites.eq.'0') then
         icepa = 0 
         lsinsh = 1
      else if (ites.eq.'1') then
         icepa = 1
         lsinsh = 1
      else if (ites.eq.'2') then
         icepa = 2 
         lsinsh = 0
      else if (ites.eq.'mrd') then
c...     mrdcepa (only VI)
         icepa = 10
         write(iwr,6002)
      else if (ites.eq.'mr(1') then
c...     mrdcepa (1)
         icepa = 12
         write(iwr,6003)
      else if (ites.eq.'mr-1') then
c...     mrdcepa-1
         icepa = 20
         write(iwr,6007)
      else if (ites.eq.'mr0') then
c...     mrdcepa (00)
         icepa = 100
         write(iwr,6004)
      else if (ites.eq.'aqcc') then
c...     aqcc 
         icepa = 101
         icorr_mr = 2
         write(iwr,6005)
      else if (ites.eq.'acpf') then
c...     acpf 
         icepa = 102
         icorr_mr = 2
         write(iwr,6006)
      else
         call caserr('first the cepa variant should be specified')
      end if
      if (icepa.lt.10)write(iwr,6001) icepa
 6001 format (/'  ***************************************************'/
     +         '  *   closed shell cepa',i1,
     +         ' option requested           *')
 6002 format(//,
     1       '  ***************************************************',/,
     2       '  *           MRDCEPA    option requested           *',/,
     3       '  *    P.J.A. Ruttink et al. JCP94,7212 (1991)      *')
 6003 format(//,
     1       '  ***************************************************',/,
     2       '  *      MRCEPA(1)(EPV+VI)   option requested       *',/,
     3       '  *    P.J.A. Ruttink et al. to be pubished         *')
 6004 format(//,
     1       '  ***************************************************',/,
     2       '  *           MRDCEPA 00 option requested           *')
 6005 format(//,
     1       '  ***************************************************',/,
     2       '  *           MR-AQCC    option requested           *',/,
     3       '  *    P.G.Szalay,R.J.Bartlett, CPL 214,481(1993)   *')
 6006 format(//,
     1       '  ***************************************************',/,
     2       '  *           MR-ACPF    option requested           *',/,
     3       '  *    R.Gdanitz,R.Ahlrichs, CPL 143,413(1988)      *')
 6007 format(//,
     1       '  ***************************************************',/,
     2       '  *          MRCEPA-1      option requested         *',/,
     3       '  *    P.J.A. Ruttink et al. to be pubished         *')
c
c...  read additional options
c
 20   call inpa4(ites)
      if (ites.eq.'sin') then
         lsinsh = 1
      else if (ites.eq.'nosi') then
         lsinsh = 0
      else if (ites.eq.'paul') then
         lcepaul = 1
       if (icepa.eq.2) write (iwr,6012)
 6012  format('  *    paul s recommended version of cepa2          *')
      else if (ites.eq.'it') then
         call inpi(itcepa)
      else if (ites.eq.'crit') then
         call inpf(critcep)
      else if (ites.eq.'prin') then
         lprcep = 1
      else if (ites.eq.'micr') then
         call inpi(mcycep)
         call inpf(crycep)
         if (crycep.eq.0.0d0) crycep = 1.0d-15
      else if (ites.eq.'brst') then
         ibrst = 1
         write(iwr,*) ' brst is the only reasonable option'
      else if (ites.eq.'pair') then
         icorr_mr = 0 
      else if (ites.eq.'proj') then
         icorr_mr = 1 
      else if (ites.eq.'vari') then
         icorr_mr = 2 
      else if (ites.eq.'doc') then
         if (icepa.ne.11.and.icepa.ne.12) 
     1   call caserr('doc subdirective only for mrcepa(1)')
         icepa = 11
      else if (ites.eq.'all') then
         if (icepa.ne.11.and.icepa.ne.12) 
     1   call caserr('all subdirective only for mrcepa(1)')
         icepa = 12
      else if (ites.eq.' ') then 
         go to 30
      else
         call caserr('unrecognised CEPA subdirective')
      end if
      go to 20
c...  print info
30    if (icepa.ge.10) then
c...    always single shift  / no consistent way to turn it off yet
         if (lsinsh.eq.0) 
     1   write(iwr,6013)
6013       format(//,' *** WARNING   *********************************',
     1           /,'                                                ',
     2           /,' The DOUBLET configurations will not be shifted ',
     3           /,' in the MRCEPA. This is slightly different from ',
     4           /,' not shifting the singles.                      ',
     5           /,'                                                ',
     6           /,' We hope YOU KNOW what you are doing, WE DONT !!',
     7           /,'                                                ',
     8           /,' *** GOOD LUCK *********************************',
     9             //)
      endif
crwah
      if (icepa.eq.11) write(iwr,6019)
 6019 format('  *    VI/EPV for inactive (DOC) only               *',/,
     1       '  * Cf. L.Fusti-Molnar,P.G.Szalay JPC100,6288(1996) *')
      if (icepa.eq.12) write(iwr,6020)
 6020 format('  *    VI/EPV for ALL                               *')
      if (lsinsh.gt.0) then
         write (iwr,6010)
 6010   format ('  *    include shifting of singles                  *')
      else 
         write (iwr,6011)
 6011  format ('  *    exclude shifting of singles                  *')
      end if
      if (lprcep.ge.1) write (iwr,6014)
 6014 format ('  *    print requested  ',6x,'                      *')
      if (ibrst.gt.0) write (iwr,6015)
 6015 format('  *    include brillouin states  -- renate          *')
      write(iwr,6016) mcycep,crycep
 6016 format('  *    # micro iterations            ',i6,'         *',/,
     1       '  *      until relative difference < ',f8.5,'       *')
      if (icorr_mr.eq.0) write(iwr,6029)
 6029  format ('  *    Ecorr from the sum of the PAIR energies      *')
      if (icorr_mr.eq.1) write(iwr,6017)
 6017  format ('  *    Ecorr from the projection of the CI-vector   *')
      if (icorr_mr.eq.2) write(iwr,6018)
 6018  format ('  *    Ecorr from the variation of reference space  *')
      write (iwr,6030) itcepa , critcep
 6030 format ('  *    start cepa after it ',i4,
     +        '                     *'/'  *            and tester <',
     +        f9.6,'                *'/
     +        '  ***************************************************')
c
      return
c
      end
**==casgen
      subroutine casgcc(call,conf,nreft,again,iwr)
c
c...  calling routine for casgen directive
c
      implicit integer (a-z)
      character*(*) call
_IF(cray,t3e,i8)
      integer excitc,conf(*)
_ELSEIF(convex,i8drct)
      integer *8 excitc,conf(*)
_ELSE
      REAL excitc,conf(*)
_ENDIF
      logical again
INCLUDE(common/symchk)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
c
      dimension ndsc(8,3),ndc(3),norbc(3),excitc(3),maxexc(3),
     1          ifsymc(3),ifspic(3)
      save ndsc,ndc,norbc,excitc,maxexc,ifsymc,ifspic,ncasset,icaspr
      data ncasset,icaspr/0,0/
c
      if (call.eq.'input') then
         ii = ncasset + 1
         call cassin(ncasset,ndsc(1,ii),ndc(ii),norbc(ii),
     +            excitc(ii),maxexc(ii),ifsymc(ii),ifspic(ii),excit,
     +            icaspr,iwr)
      else if (call.eq.'execute') then
         if (ncasset.le.0) return
c...   ** if no ref conf was given nref is 0 (1 is set by endinp though)
         nref = nreft
         again = .false.
         do ii=1,ncasset
            call casset(conf,nw,ndsc(1,ii),ndc(ii),norbc(ii),
     +                  excitc(ii),maxexc(ii),ifsymc(ii),ifspic(ii),
     +                  iwr)
         enddo
         if(icaspr.gt.0)then
           write(iwr,78)nref
78         format(///////15x,'list of',i6,' reference configurations',/,
     *                   15x,'                generated by casgen',
     *                  /15x,38('-')/)
           do ii = 1, nref
             call prcon(conf((ii-1)*nw+1),ii,1)
           enddo
         endif
      else
         call caserr(' error in casgcc ')
      end if
c
      return
      end

      subroutine cassin(ncasset,nds,nd,norbc,excit,maxex,ifsym,
     *                  ifspin,cexcit,icaspr,iwr)
c
c     input handler for cas- and fullci function generator
c
c     ncasset  .. number of times casgen is called
c     nds(8)   .. # doubly occupied/symmetry
c     nd       .. total # doubly occupied
c     norbc    .. # orbitals for casgen
c     excit    .. excitation pattern / 0 => default of 7 for all configs
c     maxex    .. maximum excitation level to current referene set
c     ifsym    .. 1 => use symmetry / 0 dont
c
c     cexcit is the current excitation allowance (input parameter)
c
c     icaspr   .. if(icaspr.gt.0) print the generated reference
c                                 configurations
c     iwr      .. print unit
c
c****       *** implementation ***
c     inpro : dimension ndsc(8,3),ndc(3),norbc(3),excitc(3),maxexc(3)
c            1          ifsymc(3),ifspic(3)
c             ncasset = 0
c     add 'casgen' to the list in pass / label 37 in inpro:
c     37  ii = ncasset + 1
c         call cassin(ncasset,ndsc(1,ii),ndc(ii),norbc(ii),
c        +            excitc(ii),maxexc(ii),ifsymc(ii),ifspic(ii),excit,
c        +            iwr)
c         go to 93
c     na call endinp(excit)  in inpro
c     c
c           if (ncasset.le.0) go to 77
c           do 76 ii=1,ncasset
c              call casset(conf,nword,ndsc(1,ii),ndc(ii),norbc(ii),
c          +               excitc(ii),maxexc(ii),ifsymc(ii),ifspic(ii),
c          +               iwr)
c     76    continue
c     77    continue
c
c=======================================================================
c     the directive reads
c     casgen  'keyword' params ... etc.  (1 card)
c     input parameters are :
c     'doc'   ndoc   : specifeis # orbitals that remain doubly occupied
c     'sdoc'  nd1 nd2  ... : as above but per symmetry
c     ** doc and sdoc are mutually exclusive **
c     ** doc orbitals are taken to be first
c     'norb'  norbc  : # orbitals taken into account (first in order)
c     'maxex' maxex  : maximum excitation level with respect to the
c                      current reference set
c     'excit'        : use current excitation allowance for the newly
c                      generated configurations only
c     ** by default the excitation mask  7b (no internal excitations)
c     'nosym'        : do not select the configs on symmetry (default)
c     'sym'          : select config on symnmetry
c     'spin'         : select configs on spin-symmetry
c     'print'        : print the generated reference configurations
c=======================================================================
c
c*      implicit REAL (a-h,o-z)
      implicit integer (a-z)
_IF(cray,t3e,i8)
      integer excit,cexcit,z7
      data z7/7/
_ELSEIF(convex,i8drct)
      integer *8 excit,cexcit,z7
      data z7/7/
_ELSEIF(bits8)
      REAL excit,cexcit
_ELSE
      REAL excit,cexcit
      REAL pad
_ENDIF
c
INCLUDE(common/orb)
INCLUDE(common/work)
c
      character*8 dir,ites
      dimension dir(9),    nds(8)
      data dir/'doc','sdoc','norb','maxex','excit','nosym','sym','spin',
     1         'print'/
c
c...  initialise
c
      ncasset = ncasset + 1
      if (ncasset.gt.3) call caserr(' too many sets in casgen')
c...  maximal 3 casgen directives
      call setsto(8,0,nds)
      jsdoc = 0
      nd = 0
      norbc = nint
      ifsym = 0
      ifspin = 0
_IF(cray,t3e,convex,i8,i8drct)
      excit = z7
_ELSEIF(bits8)
      call pad8(7,excit)
_ELSE
      excit = pad(7)
_ENDIF
      maxex = 999999
c
1     if (jrec.ge.jump) return
c
      call inpa(ites)
      i = locatc(dir,9,ites)
      if (i.eq.0) call caserr('unrecognised directive in casgen')
2     go to (10,20,30,40,50,60,70,80,90),i
c...  doc
10    call inpi(nd)
      if (jsdoc.ne.0) call caserr('not sdoc and doc please')
      go to 1
c...  sdoc
20    call inpa(ites)
      i = locatc(dir,9,ites)
      if (i.gt.0) go to 2
      jrec = jrec - 1
      jsdoc = jsdoc + 1
      call inpi(nds(jsdoc))
      if (nd.ne.0) call caserr('not sdoc and doc please')
      if (jsdoc.eq.8.or.jrec.ge.jump) go to 1
      go to 20
c...  norb
30    call inpi(norbc)
      go to 1
c...  maxex
40    call inpi(maxex)
      go to 1
c...  excit
50    excit = cexcit
      go to 1
c...  nosym
60    ifsym = 0
      write(iwr,61)
61    format(' **casgen : nosym is default **')
      go to 1
c...  sym
70    ifsym = 1
      go to 1
c...  spin
80    ifspin = 1
      go to 1
c...  print
90    icaspr = icaspr + 1
      go to 1
c
      end
      subroutine casset(noc,nword,nds,ndc,norb,excit,maxex,ifsym,ifspin,
     +                  iwr)
c
c     routine to produce (full) ci csf-set (d2h symmetry adapted)
c     written by p.j.a.ruttink (utrecht,1985)
c     atmol adaption j.h. v. lenthe
c
c     noc : array to contain csf occupations etc. (atmol format)
c     nword : # words used to store each config in atmol4 ci
c     maxex : maximum excitation level with respect to reference set
c     nds   : # of doubly occupied orbitals per symmetry
c     ndc     : # doubly occupied orbitals / instead of nds
c     ifsym :   => symmetry selected set is produced
c     ifspin:   => spin-impossible configs are eliminated
c     norb  : # (internal) orbitals involved (max. 62)
c     excit : excitation level for configs
c     iwr   : print unit
c
c     *** popcnt = sum1s ***  ==> kbit
c
      implicit REAL  (a-h,o-z),integer (i-n)
c
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt
_ENDIF
      integer occ, stms, sym, et, stsym
_IF(cray,t3e,i8)
      integer z0,z7,joker,excit,noc
      integer occ0,occ1,occ2,occ3,occ4,occ5,occ6,occ7,occ8
_ELSEIF(convex,i8drct)
      integer *8 z0,z7,joker,excit,noc
      integer *8 occ0,occ1,occ2,occ3,occ4,occ5,occ6,occ7,occ8
      integer *8 shiftr
      external shiftr
_ELSE
      REAL z0,z7,joker,excit,noc
_ENDIF
c
INCLUDE(common/sizes)
INCLUDE(common/orb)
INCLUDE(common/ccntl)
INCLUDE(common/table)
c
INCLUDE(common/helpr)
INCLUDE(common/corctl)
INCLUDE(common/mapp)
c
      dimension occ(62),ms(8),msd(8),et(8),stms(8),joct(22)
      dimension noc(nword,*),nds(8)
c
      dimension occ0(2),occ1(2),occ2(2),occ3(2),occ4(2),
     1          occ5(2),occ6(2),occ7(2),occ8(2)
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z7/0,7/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(7,z7)
_ELSE
      z0 = pad(0)
      z7 = pad(7)
_ENDIF
c
c...  if 1st config is made by endinp give it the same excitation pattern
c...  as the rest
      if (nref.eq.0) then
         nref = 1
_IF(i8drct)
         noc(nword,1) = excit
_ELSEIF(linux)
         call dcopy(1,excit,1,noc(nword,1),1)
_ELSE
         noc(nword,1) = excit
_ENDIF
      endif
c
c...  sort out nds / ndc
c
      if (ndc.eq.0) go to 4
      do 3 i=1,ndc
        do 3 j = 1, nint
          if(mapie(j).eq.i) then
            nds(iro(j)) = nds(iro(j)) + 1
          endif
3     continue
c
c # of orbitals per irrep
c
4     icon = nref
_IF(cray,t3e,convex,i8,i8drct)
      occ0(1) = 0
      occ0(2) = 0
_ELSE
      occ0(1) = 0.0d0
      occ0(2) = 0.0d0
_ENDIF
      do i=1,8
         msd(i) = nds(i)
      end do
      call setsto(norb,0,occ)
      call setsto(8,0,ms)
      call setsto(8,0,et)
      stsym = ifrep
      ndoc   = nds(1) + nds(2) + nds(3) + nds(4) +
     1         nds(5) + nds(6) + nds(7) + nds(8)
c
      do 5 iorb = 1,norb
      sym = iro(iorb)
      if (msd(sym).eq.0) goto 6
      msd(sym) = msd(sym) - 1
      icrerr = icrer(occ0,iorb,2)
      if (icrerr.eq.0) call caserr('not able to create in casset')
      occ(iorb) = 2
      goto 5
6     ms(sym) = ms(sym) + 1
5     continue
c
      nelv = nelec
      stms(1) = nds(1) + 1
      do 7 i = 1,8
      msd(i) = 2*ms(i)
      if (i.gt.1) stms(i) = stms(i-1) + ms(i-1) + nds(i)
      nelv = nelv - 2*nds(i)
7     continue
      if (nelv.le.0) call caserr('wrong # electrons in casset')
c
c initialize the partitioning of the electrons over the irreps
c
      nr = nelv
c
      do 10 i = 1,8
      et(i) = nr
      if (nr.le.msd(i)) goto 20
      et(i) = msd(i)
      nr = nr - et(i)
10    continue
20    continue
c
c symmetry for this partitioning
c
75    sym = 1
      do 30 i = 1,8
c*ap     use of integer *and*
      if (IAND32(et(i),1).eq.1) sym = mult(sym,i)
30    continue
c
c check the symmetry
c
      if (sym.eq.stsym.or.ifsym.eq.0) goto 100
c
c next weight partitioning
c
40    nr  = 0
      nri = 0
c
      do 50 i = 1,7
      if (et(i).eq.0) goto 55
      nri = et(i)
      nr = nr + nri
55    if (nri.eq.0.or.et(i+1).eq.msd(i+1)) goto 50
      et(i+1) = et(i+1) + 1
      nr = nr - 1
      if (i.gt.0) call setsto(i,0,et)
      do 60 j = 1,i
      et(j) = nr
      if (nr.le.msd(j)) goto 75
      et(j) = msd(j)
      nr = nr - et(j)
60    continue
      goto 75
50    continue
c
c ready
c
      goto 300
c
c occupation number schemes for this partitioning
c
100   continue
      call stocc(occ(stms(1)),ms(1),et(1),l1)
110   call pock(occ(stms(1)),ms(1),stms(1),occ0,occ1,l1)
      call stocc(occ(stms(2)),ms(2),et(2),l2)
120   call pock(occ(stms(2)),ms(2),stms(2),occ1,occ2,l2)
      call stocc(occ(stms(3)),ms(3),et(3),l3)
130   call pock(occ(stms(3)),ms(3),stms(3),occ2,occ3,l3)
      call stocc(occ(stms(4)),ms(4),et(4),l4)
140   call pock(occ(stms(4)),ms(4),stms(4),occ3,occ4,l4)
      call stocc(occ(stms(5)),ms(5),et(5),l5)
150   call pock(occ(stms(5)),ms(5),stms(5),occ4,occ5,l5)
      call stocc(occ(stms(6)),ms(6),et(6),l6)
160   call pock(occ(stms(6)),ms(6),stms(6),occ5,occ6,l6)
      call stocc(occ(stms(7)),ms(7),et(7),l7)
170   call pock(occ(stms(7)),ms(7),stms(7),occ6,occ7,l7)
      call stocc(occ(stms(8)),ms(8),et(8),l8)
180   call pock(occ(stms(8)),ms(8),stms(8),occ7,occ8,l8)
c
c check the number of open shells
c
      nsin = nsior(occ8)
      if (ifspin.eq.1.and.nsin.lt.mspin-1) goto 190
c
c check the excitation level
c
      minnex = nelec
      do 200 iref = 1,nref
_IF1() PAUL try first version on cray
_IF(cray,t3e,convex,i8,i8drct)
      nex = popcnt(ieor(occ8(1),noc(1,iref)))
      if (nwm1.eq.2) nex = nex + popcnt(ieor(occ8(2),noc(2,iref)))
_ELSEIF(bits8)
      call dxor8(occ8(1),noc(1,iref),tmp8)
      nex = popcnt(tmp8)
      if (nwm1.eq.2) then
        call dxor8(occ8(2),noc(2,iref),tmp8)
        nex = nex + popcnt(tmp8)
      endif
_ELSE
      nex = popcnt(dxor(occ8(1),noc(1,iref)))
      if (nwm1.eq.2) nex = nex + popcnt(dxor(occ8(2),noc(2,iref)))
_ENDIF
      minnex = min(minnex,nex/2)
200   continue
      if (minnex.eq.0.or.minnex.gt.maxex) goto 190
c
c add a configuration
c
      icon = icon + 1
      nuse = nuse + nword
      call corfai('congen')
      do i=1,nword
         noc(i,icon) = z0
      end do
_IF(i8drct)
      do i=1,nwm1
       noc(i,icon) = occ8(i)
      enddo
_ELSEIF(linux)
      call dcopy(nwm1,occ8,1,noc(1,icon),1)
_ELSE
      call fmove(occ8,noc(1,icon),nwm1)
_ENDIF
      
c***  set excitation pattern
_IF(i8drct)
      noc(nword,icon) = excit
_ELSEIF(linux)
      call dcopy(1,excit,1,noc(nword,icon),1)
_ELSE
      noc(nword,icon) = excit
_ENDIF
c
c next weight
c
190   call shifw(occ(stms(8)),ms(8),l8)
      if (l8.eq.0) goto 180
      call shifw(occ(stms(7)),ms(7),l7)
      if (l7.eq.0) goto 170
      call shifw(occ(stms(6)),ms(6),l6)
      if (l6.eq.0) goto 160
      call shifw(occ(stms(5)),ms(5),l5)
      if (l5.eq.0) goto 150
      call shifw(occ(stms(4)),ms(4),l4)
      if (l4.eq.0) goto 140
      call shifw(occ(stms(3)),ms(3),l3)
      if (l3.eq.0) goto 130
      call shifw(occ(stms(2)),ms(2),l2)
      if (l2.eq.0) goto 120
      call shifw(occ(stms(1)),ms(1),l1)
      if (l1.eq.0) goto 110
c
      goto 40
c
300   nref = icon
c...  unpack excitation pattern
_IF(i8drct)
      joker = excit
_ELSEIF(linux)
      call dcopy(1,excit,1,joker,1)
_ELSE
      joker = excit
_ENDIF
      do 4321 i=1,22

_IF1() PAUL  Try ishft version on cray,t3e
_IF(cray,t3e,i8drct)
         joct(22-i+1) = iand(joker,z7)
         joker = shiftr(joker,3)
_ELSEIF(convex,i8)
         joct(22-i+1) = iand(joker,z7)
         joker = ishft(joker,-3)
_ELSEIF(bits8)
         call dand8(joker,z7,tmp8)
         joct(22-i+1) = ipad(tmp8)
         call shiftr8(joker,3,joker)
_ELSE
         joct(22-i+1) = ipad(dand(joker,pad(7)))
         joker = shiftr(joker,3)
_ENDIF
4321  continue
c
      write(iwr,601) nref,maxex,joct,nds,norb
601   format(//' ************************************************',
     1        /' *           **cas** ci generator               *',
     2        /' *  total # configs ',i6,'                      *',
     3        /' *  max. excitation level ',i6,'                *',
     4        /' *  excitation pattern  ',22i1,              '  *',
     5        /' *  ndoc  ',8i3,                 '              *',
     6        /' *  # orbitals considered ',i6,'                *')
      if (ifspin.eq.1) write(iwr,604)
604   format(  ' *  **spin adapted**                            *')
      if (ifspin.eq.0) write(iwr,605)
605   format(  ' *  **not necessarily spin adapted**            *')
      if (ifsym.eq.1) write(iwr,602)
602   format(  ' *  **symmetry adapted**                        *',
     1        /' ************************************************')
      if (ifsym.eq.0) write(iwr,603)
603   format(  ' *  no symmetry used                            *',
     1        /' ************************************************')
      return
      end
      subroutine stocc(occ,nm,nel,last)
      implicit REAL (a-h,o-z)
      integer occ(nm)
c
c construct starting occupations
c for nel electrons in nm orbitals
c
      if (nel.gt.2*nm) call caserr('too many electrons in stocc')
      last = 1
      if (nm.eq.0) return
      last = 0
      nd = nel/2
      if (nd.eq.0) goto 20
      call setsto(nd,2,occ)
      if (nd.eq.nm) return
20    occ(nd+1) = nel - nd*2
      if (nm.eq.nd+1) return
      call setsto(nm-nd-1,0,occ(nd+2))
      return
      end
      subroutine pock(occ,nm,stms,occin,occout,last)
c*ap*      implicit REAL (a-h,o-z)
      implicit integer (a-z)
      REAL  occin,occout
      dimension occ(nm),occin(2),occout(2)
c
c add the electrons in occ
c
      occout(1) = occin(1)
      occout(2) = occin(2)
      if (last.ne.0) return
c
      stm0 = stms - 1
      do 10 iorb = 1,nm
         ni = occ(iorb)
         if (ni.eq.0) goto 10
         icrerr = icrer(occout,stm0+iorb,ni)
         if (icrerr.eq.0) call caserr('not able to create in pock')
10    continue
      return
      end
**==icre.f
      function icrer(conf,i,nel)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer z(2),conf,imas
_ELSEIF(convex)
      integer *8 z(2),conf
      integer *8 imas
_ELSEIF(i8drct)
      integer *8 z(2),conf
      integer *8 imas
      integer *8 shift
      external shift
_ELSE
      REAL  imas,conf
_ENDIF
INCLUDE(common/cdcryi)
      dimension conf(*)
INCLUDE(common/orb)
c     create nel electrons in orbital i in configuration conf
c     see iann
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z/ 1, 3/
      imas = z(nel)
_ELSE
_IF(bits8)
      call pad8(nel,imas)
      if(nel.eq.2) call pad8(3,imas)
_ELSE
      imas = pad(nel)
      if (nel.eq.2) imas = pad(3)
_ENDIF
_ENDIF
c.... determine which word and which position
      kc = (i-1)/n32
      ks = (i-kc*n32)*2
      kc = kc + 1
c.... check occupation of i
_IF1() PAUL try second clause on cray
_IF(cray,t3e)
      iocc = iand(shift(conf(kc),ks),z(2))
_ELSEIF(convex,i8)
      iocc = iand(ishftc(conf(kc),ks,64),z(2))
_ELSEIF(i8drct)
      iocc = iand(shift(conf(kc),ks),z(2))
_ELSEIF(bits8)
      call shift8(conf(kc),ks,tmp8)
      call pad8(3,pad3)
      call dand8(tmp8,pad3,tmp88)
      iocc = ipad(tmp88)
_ELSEIF(sgi,ipsc)
      iocc = ipad(dand(shiftc(conf(kc),ks),pad(3)))
_ELSE
      iocc = ipad(dand(shift(conf(kc),ks),pad(3)))
_ENDIF
      if (iocc.gt.0) then
c.... electrons in i => failure
         icrer = 0
         return
      end if
      icrer = 1
_IF1() PAUL try first clause on cray
_IF(cray,t3e,convex,i8)
      conf(kc) = ior(conf(kc),ishftc(imas,n64-ks,64))
_ELSEIF(i8drct)
      conf(kc) = ior(conf(kc),shift(imas,n64-ks))
_ELSEIF(bits8)
      call shift8(imas,n64-ks,tmp8)
      call dor8(conf(kc),tmp8,conf(kc))
_ELSEIF(sgi,ipsc)
      conf(kc) = dor(conf(kc),shiftc(imas,n64-ks))
_ELSE
      conf(kc) = dor(conf(kc),shift(imas,n64-ks))
_ENDIF

      return
      end
      subroutine shifw(occ,nm,last)
      implicit REAL (a-h,o-z)
      integer occ(nm)
c
c construct the occupations with the next weight
c
      if (last.ne.0) return
c
      last = 1
      if (nm.eq.1) return
c
      iorbl = nm - 1
      nr = 0
      do 10 iorb = 1,iorbl
      nr = nr + occ(iorb)
      if (occ(iorb).eq.0.or.occ(iorb+1).eq.2) goto 10
      occ(iorb+1) = occ(iorb+1) + 1
      call stocc(occ,iorb,nr-1,last)
      goto 20
c
10    continue
c
20    return
      end
**==symb22.f
      subroutine symb22(mexlhs,mexrhs)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension jlif6(16),jlif7(16),jlif8(16),mtrm6(16),mtrm7(16),
     *mtrm8(16)
      dimension mexlhs(*),mexrhs(*)
      common/scra  /mexlh(mxcan2),mexrh(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),
     *mrange(mxcan1),x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif6(16),ilif7(16),ilif8(16),ntrm6(16),ntrm7(16),ntrm8(16)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,
     *min8,max8,ndimax
      equivalence (jlif6(1),ilif2(1)),(jlif7(1),ilif3(1)),
     *(jlif8(1),jlif2(1)),(mtrm6(1),ntrm2(1)),
     *(mtrm7(1),ntrm3(1)),(mtrm8(1),mtrm2(1))
      data oneneg/-1.0d0/
      kounte = 0
c****
c**** double orthogonality
c**** compute spin coupling coefficient for
c**** <i j/k l>=direct and <i l/k j>=exchange forms
c**** where i and k more occupied on rhs
c**** and   j and l more occupied on lhs
c**** and i lt j and i le k and j le l
c****
      if ((npopr(i)-npopl(i)).ne.2) then
        if ((npopl(j)-npopr(j)).ne.2) then
          call toprh
          if (k.lt.j) then
c...
c...   ***case 1***
c...
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,
     +                 k)
c...
c...
c... generate rh double at level k
c... see     *****table 5*****
c...
c...
            ibase = i16(k)
            call min678(k)
            kont = 0
            if (npopl(k).ne.0) then
c...
c... occupations lhs=1  rhs=2
c...
              if (min6.lt.min2) min6 = min2
              if (max6.gt.max2) max6 = max2
              if (min8.lt.min3) min8 = min3
              if (max8.gt.max3) max8 = max3
              if (max6.ge.min6) then
c... generate type 6 points
                do 30 ib = min6 , max6 , 2
                  ilif6(ib) = kont
                  ntrm = ntrm2(ib)
                  ntrm6(ib) = ntrm
                  ilif = ilif2(ib)
                  do 20 loop = 1 , ntrm
                    jexlhs(kont+loop) = lexlhs(ilif+loop)
                    jexrhs(kont+loop) = lexrhs(ilif+loop)
                    jrange(kont+loop) = mrange(ilif+loop)
                    xj(kont+loop) = x(ilif+loop)*vneg(ib+1)
 20               continue
                  kont = kont + ntrm
 30             continue
              end if
              if (max7.ge.min7) then
c... generate type 7 points
                do 60 ib = min7 , max7 , 2
                  ilif7(ib) = kont
                  if (ib.ge.min3 .and. ib.le.max3) then
c... coupling lhs=2 rhs=3
                    ntrm = ntrm3(ib)
                    ilif = ilif3(ib)
                    do 40 loop = 1 , ntrm
                      jexlhs(kont+loop) = lexlhs(ilif+loop)
                      jexrhs(kont+loop) = lexrhs(ilif+loop)
                      jrange(kont+loop) = mrange(ilif+loop)
                      xj(kont+loop) = x(ilif+loop)*p2(ib-1)
                      yj(kont+loop) = x(ilif+loop)*tneg
 40                 continue
                    kont = kont + ntrm
                  end if
                  if (ib.ge.min2 .and. ib.le.max2) then
c... coupling lhs=1 rhs=3
                    ntrm = ntrm2(ib)
                    ilif = ilif2(ib)
                    do 50 loop = 1 , ntrm
                      jexlhs(kont+loop) = lexlhs(ilif+loop)
     +                  + n2l(ib+ibase)
                      jexrhs(kont+loop) = lexrhs(ilif+loop)
                      jrange(kont+loop) = mrange(ilif+loop)
                      xj(kont+loop) = x(ilif+loop)*q2(ib)
                      yj(kont+loop) = x(ilif+loop)*tneg
 50                 continue
                    kont = kont + ntrm
                  end if
                  ntrm7(ib) = kont - ilif7(ib)
 60             continue
              end if
              if (max8.ge.min8) then
c... generate type 8 points
                do 80 ib = min8 , max8 , 2
                  ilif8(ib) = kont
                  ntrm = ntrm3(ib)
                  ntrm8(ib) = ntrm
                  ilif = ilif3(ib)
                  do 70 loop = 1 , ntrm
                    jexlhs(kont+loop) = lexlhs(ilif+loop)
     +                                  + n2l(ib+ibase-2)
                    jexrhs(kont+loop) = lexrhs(ilif+loop)
                    jrange(kont+loop) = mrange(ilif+loop)
                    xj(kont+loop) = x(ilif+loop)*wneg(ib-1)
 70               continue
                  kont = kont + ntrm
 80             continue
              end if
            else
c...
c... occupations lhs=0 rhs=1
c...
              if (min6.lt.min2) min6 = min2 - 1
              if (max6.ge.max2) max6 = max2 - 1
              if (min8.le.min3) min8 = min3 + 1
              if (max8.gt.max3) max8 = max3 + 1
              if (max6.ge.min6) then
c... generate type 6 points
                do 100 ib = min6 , max6 , 2
                  ilif6(ib) = kont
                  ntrm = ntrm2(ib+1)
                  ntrm6(ib) = ntrm
                  ilif = ilif2(ib+1)
                  do 90 loop = 1 , ntrm
                    jexlhs(kont+loop) = lexlhs(ilif+loop)
                    jexrhs(kont+loop) = lexrhs(ilif+loop)
     +                                  + n2r(ib+ibase)
                    jrange(kont+loop) = mrange(ilif+loop)
                    xj(kont+loop) = x(ilif+loop)
 90               continue
                  kont = kont + ntrm
 100            continue
              end if
              if (max7.ge.min7) then
c... generate type 7 points
                min3m1 = min3 - 1
                max2p1 = max2 + 1
                do 130 ib = min7 , max7 , 2
                  v2ib = v2(ib)
                  ilif7(ib) = kont
                  if (ib.lt.max3 .and. ib.ge.min3m1) then
c... coupling lhs=0 rhs=1
                    ntrm = ntrm3(ib+1)
                    ilif = ilif3(ib+1)
                    do 110 loop = 1 , ntrm
                      jexlhs(kont+loop) = lexlhs(ilif+loop)
                      jexrhs(kont+loop) = lexrhs(ilif+loop)
     +                  + n2r(ib+ibase)
                      jrange(kont+loop) = mrange(ilif+loop)
                      xj(kont+loop) = x(ilif+loop)*w2neg(ib)
                      yj(kont+loop) = x(ilif+loop)*v2ib
 110                continue
                    kont = kont + ntrm
                  end if
                  if (ib.gt.min2 .and. ib.le.max2p1) then
c... coupling lhs=0  rhs=2
                    ntrm = ntrm2(ib-1)
                    ilif = ilif2(ib-1)
                    do 120 loop = 1 , ntrm
                      jexlhs(kont+loop) = lexlhs(ilif+loop)
                      jexrhs(kont+loop) = lexrhs(ilif+loop)
                      jrange(kont+loop) = mrange(ilif+loop)
                      yj(kont+loop) = x(ilif+loop)*w2(ib)
                      xj(kont+loop) = x(ilif+loop)*v2ib
 120                continue
                    kont = kont + ntrm
                  end if
                  ntrm7(ib) = kont - ilif7(ib)
 130            continue
              end if
              if (max8.ge.min8) then
c... generate type 8 points
                do 150 ib = min8 , max8 , 2
                  ilif8(ib) = kont
                  ntrm = ntrm3(ib-1)
                  ntrm8(ib) = ntrm
                  ilif = ilif3(ib-1)
                  do 140 loop = 1 , ntrm
                    jexlhs(kont+loop) = lexlhs(ilif+loop)
                    jexrhs(kont+loop) = lexrhs(ilif+loop)
                    jrange(kont+loop) = mrange(ilif+loop)
                    xj(kont+loop) = x(ilif+loop)
 140              continue
                  kont = kont + ntrm
 150            continue
              end if
            end if
c...
c...
c... recur through rh double segments from k+1 to j-1
c...   see   ****table 6****
c...
c...
            kp1 = k + 1
            jm1 = j - 1
            if (jm1.ge.kp1) then
              do 330 kk = kp1 , jm1
                if (npopl(kk).eq.1) then
c... skip over empty or doubly occupied levels
                  max6p1 = max6 + 1
                  max7p1 = max7 + 1
                  max8p1 = max8 + 1
                  min6p1 = min6 + 1
                  min7p1 = min7 + 1
                  min8p1 = min8 + 1
                  max6m1 = max6 - 1
                  max7m1 = max7 - 1
                  max8m1 = max8 - 1
                  min6m1 = min6 - 1
                  min7m1 = min7 - 1
                  min8m1 = min8 - 1
                  ibase = i16(kk)
                  call min678(kk)
                  kont = 0
                  if (max6.ge.min6) then
c...
c... generate type 6 points
c...
                    do 190 ib = min6 , max6 , 2
                      jlif6(ib) = kont
                      n2rib = n2r(ib+ibase)
                      if (ib.le.max6m1 .and. ib.ge.min6m1) then
c... coupling lhs=1  rhs=1
                        ntrm = ntrm6(ib+1)
                        ilif = ilif6(ib+1)
                        do 160 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
     +                      + n2l(ib+ibase+2)
                          kexrhs(kont+loop) = jexrhs(ilif+loop) + n2rib
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)
 160                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.ge.min6p1 .and. ib.le.max6p1) then
c... coupling lhs=2  rhs=2
                        ntrm = ntrm6(ib-1)
                        ilif = ilif6(ib-1)
                        do 170 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*f(ib)
 170                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.le.max7m1 .and. ib.ge.min7m1) then
c... coupling lhs=2 rhs=1
                        ntrm = ntrm7(ib+1)
                        ilif = ilif7(ib+1)
                        do 180 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
                          kexrhs(kont+loop) = jexrhs(ilif+loop) + n2rib
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*factng(ib)
 180                    continue
                        kont = kont + ntrm
                      end if
                      mtrm6(ib) = kont - jlif6(ib)
 190                continue
                  end if
                  if (max7.ge.min7) then
c...
c... generate type 7 points
c...
                    do 240 ib = min7 , max7 , 2
                      jlif7(ib) = kont
                      n2lib = n2l(ib+ibase)
                      n2rib = n2r(ib+ibase)
                      if (ib.le.max7m1 .and. ib.ge.min7m1) then
c... coupling  lhs=1  rhs=1
                        ntrm = ntrm7(ib+1)
                        ilif = ilif7(ib+1)
                        do 200 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop) + n2lib
                          kexrhs(kont+loop) = jexrhs(ilif+loop) + n2rib
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*f(ib)
                          x(kont+loop) = yj(ilif+loop)
 200                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.ge.min7p1 .and. ib.le.max7p1) then
c... coupling lhs=2  rhs=2
                        ntrm = ntrm7(ib-1)
                        ilif = ilif7(ib-1)
                        do 210 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*f(ib-1)
                          x(kont+loop) = yj(ilif+loop)
 210                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.ge.min6p1 .and. ib.le.max6p1) then
c... coupling lhs=1 rhs=2
                        ntrm = ntrm6(ib-1)
                        ilif = ilif6(ib-1)
                        do 220 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop) + n2lib
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*factn(ib)
                          x(kont+loop) = 0.0d0
 220                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.le.max8m1 .and. ib.ge.min8m1) then
c... coupling lhs=2 rhs=1
                        ntrm = ntrm8(ib+1)
                        ilif = ilif8(ib+1)
                        do 230 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
                          kexrhs(kont+loop) = jexrhs(ilif+loop) + n2rib
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*factng(ib-1)
                          x(kont+loop) = 0.0d0
 230                    continue
                        kont = kont + ntrm
                      end if
                      mtrm7(ib) = kont - jlif7(ib)
 240                continue
                  end if
                  if (max8.ge.min8) then
c...
c... generate type 8 points
c...
                    do 280 ib = min8 , max8 , 2
                      jlif8(ib) = kont
                      n2lib = n2l(ib+ibase-2)
                      if (ib.le.max8m1 .and. ib.ge.min8m1) then
c... coupling lhs=1  rhs=1
                        ntrm = ntrm8(ib+1)
                        ilif = ilif8(ib+1)
                        do 250 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop) + n2lib
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
     +                      + n2r(ib+ibase)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*f(ib-1)
 250                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.ge.min8p1 .and. ib.le.max8p1) then
c... coupling lhs=2  rhs=2
                        ntrm = ntrm8(ib-1)
                        ilif = ilif8(ib-1)
                        do 260 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop)
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)
 260                    continue
                        kont = kont + ntrm
                      end if
                      if (ib.ge.min7p1 .and. ib.le.max7p1) then
c... coupling  lhs=1  rhs=2
                        ntrm = ntrm7(ib-1)
                        ilif = ilif7(ib-1)
                        do 270 loop = 1 , ntrm
                          kexlhs(kont+loop) = jexlhs(ilif+loop) + n2lib
                          kexrhs(kont+loop) = jexrhs(ilif+loop)
                          krange(kont+loop) = jrange(ilif+loop)
                          xk(kont+loop) = xj(ilif+loop)*factn(ib-1)
 270                    continue
                        kont = kont + ntrm
                      end if
                      mtrm8(ib) = kont - jlif8(ib)
 280                continue
c... shuffle generated data back to correct locations
                    do 290 ib = min8 , max8 , 2
                      ntrm8(ib) = mtrm8(ib)
                      ilif8(ib) = jlif8(ib)
 290                continue
                  end if
                  if (max7.ge.min7) then
                    do 300 ib = min7 , max7 , 2
                      ntrm7(ib) = mtrm7(ib)
                      ilif7(ib) = jlif7(ib)
 300                continue
                  end if
                  if (max6.ge.min6) then
                    do 310 ib = min6 , max6 , 2
                      ntrm6(ib) = mtrm6(ib)
                      ilif6(ib) = jlif6(ib)
 310                continue
                  end if
                  do 320 loop = 1 , kont
                    jexlhs(loop) = kexlhs(loop)
                    jexrhs(loop) = kexrhs(loop)
                    jrange(loop) = krange(loop)
                    xj(loop) = xk(loop)
                    yj(loop) = x(loop)
 320              continue
                end if
 330          continue
            end if
c...
c...
c... transition from rh double to rh single at level j
c...   see ****table 7****
c...
c...
            ibase = i16(j)
            call min23(j)
            kont = 0
            if (npopr(j).ne.0) then
c...
c... occupations lhs=2 rhs=1
c...
              if (max2.ge.min2) then
c... generate type 2 points
                min7m1 = min7 - 1
                max6p1 = max6 + 1
                do 360 ib = min2 , max2 , 2
                  ilif2(ib) = kont
                  if (ib.lt.max7 .and. ib.ge.min7m1) then
c... coupling  lhs=3  rhs=1
                    ntrm = ntrm7(ib+1)
                    ilif = ilif7(ib+1)
                    do 340 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
     +                  + n2r(ib+ibase)
                      mrange(kont+loop) = jrange(ilif+loop)
                      t1 = xj(ilif+loop)*p2(ib)
                      t2 = yj(ilif+loop)*tneg
                      x(kont+loop) = t1 + t2
                      xk(kont+loop) = t2 - t1
 340                continue
                    kont = kont + ntrm
                  end if
                  if (ib.gt.min6 .and. ib.le.max6p1) then
c... coupling lhs=3  rhs=2
                    ntrm = ntrm6(ib-1)
                    ilif = ilif6(ib-1)
                    do 350 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      xk(kont+loop) = xj(ilif+loop)*w(ib)
                      x(kont+loop) = -xk(kont+loop)
 350                continue
                    kont = kont + ntrm
                  end if
                  ntrm2(ib) = kont - ilif2(ib)
 360            continue
              end if
              if (max3.ge.min3) then
c... generate type 3 points
                min8m1 = min8 - 1
                max7p1 = max7 + 1
                do 390 ib = min3 , max3 , 2
                  ilif3(ib) = kont
                  if (ib.lt.max8 .and. ib.ge.min8m1) then
c... coupling  lhs=3  rhs=1
                    ntrm = ntrm8(ib+1)
                    ilif = ilif8(ib+1)
                    do 370 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
     +                  + n2r(ib+ibase)
                      mrange(kont+loop) = jrange(ilif+loop)
                      xk(kont+loop) = xj(ilif+loop)*v(ib)
                      x(kont+loop) = -xk(kont+loop)
 370                continue
                    kont = kont + ntrm
                  end if
                  if (ib.gt.min7 .and. ib.le.max7p1) then
c... coupling lhs=3  rhs=2
                    ntrm = ntrm7(ib-1)
                    ilif = ilif7(ib-1)
                    do 380 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      t1 = xj(ilif+loop)*q2(ib-1)
                      t2 = yj(ilif+loop)*tneg
                      x(kont+loop) = t1 + t2
                      xk(kont+loop) = t2 - t1
 380                continue
                    kont = kont + ntrm
                  end if
                  ntrm3(ib) = kont - ilif3(ib)
 390            continue
              end if
            else
c...
c... occupations lhs=1 rhs=0
c...
              if (max2.ge.min2) then
c... generate type 2 points
                do 420 ib = min2 , max2 , 2
                  ilif2(ib) = kont
                  if (ib.ge.min6 .and. ib.le.max6) then
c... coupling lhs=1 rhs=0
                    ntrm = ntrm6(ib)
                    ilif = ilif6(ib)
                    do 400 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
     +                  + n2l(ib+ibase+1)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      x(kont+loop) = xj(ilif+loop)
                      xk(kont+loop) = -xj(ilif+loop)
 400                continue
                    kont = kont + ntrm
                  end if
                  if (ib.ge.min7 .and. ib.le.max7) then
c... coupling lhs=2  rhs=0
                    ntrm = ntrm7(ib)
                    ilif = ilif7(ib)
                    do 410 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      t2 = yj(ilif+loop)*v2(ib)
                      t1 = xj(ilif+loop)*w2(ib)
                      x(kont+loop) = t2 - t1
                      xk(kont+loop) = t2 + t1
 410                continue
                    kont = kont + ntrm
                  end if
                  ntrm2(ib) = kont - ilif2(ib)
 420            continue
              end if
              if (max3.ge.min3) then
c... generate type 3 points
                do 450 ib = min3 , max3 , 2
                  ilif3(ib) = kont
                  if (ib.ge.min7 .and. ib.le.max7) then
c... coupling lhs=1  rhs=0
                    ntrm = ntrm7(ib)
                    ilif = ilif7(ib)
                    do 430 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
     +                  + n2l(ib+ibase-1)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      t2 = yj(ilif+loop)*w2(ib)
                      t1 = xj(ilif+loop)*v2(ib)
                      x(kont+loop) = t1 + t2
                      xk(kont+loop) = t2 - t1
 430                continue
                    kont = kont + ntrm
                  end if
                  if (ib.ge.min8 .and. ib.le.max8) then
c... coupling lhs=2 rhs=0
                    ntrm = ntrm8(ib)
                    ilif = ilif8(ib)
                    do 440 loop = 1 , ntrm
                      lexlhs(kont+loop) = jexlhs(ilif+loop)
                      lexrhs(kont+loop) = jexrhs(ilif+loop)
                      mrange(kont+loop) = jrange(ilif+loop)
                      x(kont+loop) = xj(ilif+loop)
                      xk(kont+loop) = -xj(ilif+loop)
 440                continue
                    kont = kont + ntrm
                  end if
                  ntrm3(ib) = kont - ilif3(ib)
 450            continue
              end if
            end if
c...
c... complete exchange part
c...
            do 460 loop = 1 , kont
              kexlhs(loop) = lexlhs(loop)
              kexrhs(loop) = lexrhs(loop)
              krange(loop) = mrange(loop)
 460        continue
            if (max2.ge.min2) then
              do 470 ib = min2 , max2 , 2
                mtrm2(ib) = ntrm2(ib)
                jlif2(ib) = ilif2(ib)
 470          continue
            end if
            if (max3.ge.min3) then
              do 480 ib = min3 , max3 , 2
                mtrm3(ib) = ntrm3(ib)
                jlif3(ib) = ilif3(ib)
 480          continue
            end if
            call midrh(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,mtrm3,jlif3,
     +                 j,l)
            call botrh(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,mtrm3,jlif3,
     +                 mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,l,kounte)
            call point1(mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,l,norb1,
     +                  kounte)
c...
c... complete direct term
c...
            call min23(j)
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,j,
     +                 l)
          else
c...
c...  ***case 2 or 3***
c...
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,
     +                 j)
            call rhlhd
            call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +                 kexlhs,kexrhs,krange,xk,mtrm2,jlif2,j,junke)
            if (k.lt.l) then
c...
c...   ***case 2***
c... first complete exchange term
c...
              call midlhd(j,k)
              call lhdrh(k)
              call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +                   k,l)
              call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +                   mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,l,kounte)
              call point1(mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,l,norb1,
     +                    kounte)
c...
c... now for direct term
c...
              call minmax(j)
              call point1(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,j,k,junke)
c...
c...
c... top segment of 2nd rh single range at level k
c... see     *****table 1*****
c...
              kont = 0
              ibase = i16(k)
              call min23(k)
              if (npopl(k).ne.0) then
c...
c... occupations  lhs=1 rhs=2
c...
                if (min2.lt.mmin) min2 = mmin
                if (min3.lt.mmin) min3 = mmin
                if (max2.gt.mmax) max2 = mmax
                if (max3.gt.mmax) max3 = mmax
                if (max2.ge.min2) then
c... generate type 2 points
                  do 500 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    ntrm = mtrm2(ib)
                    ntrm2(ib) = ntrm
                    ilif = jlif2(ib)
                    do 490 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)*v(ib)
 490                continue
                    kont = kont + ntrm
 500              continue
                end if
                if (max3.ge.min3) then
c... generate type 3 points
                  do 520 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    ntrm = mtrm2(ib)
                    ntrm3(ib) = ntrm
                    ilif = jlif2(ib)
                    do 510 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
     +                  + n2l(ib+ibase-1)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)*w(ib)
 510                continue
                    kont = kont + ntrm
 520              continue
                end if
              else
c...
c... occupations  lhs=0  rhs=1
c...
                if (min2.lt.mmin) min2 = mmin - 1
                if (min3.le.mmin) min3 = mmin + 1
                if (max2.ge.mmax) max2 = mmax - 1
                if (max3.gt.mmax) max3 = mmax + 1
                if (max2.ge.min2) then
c... generate type 2 points
                  do 540 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    ntrm = mtrm2(ib+1)
                    ntrm2(ib) = ntrm
                    ilif = jlif2(ib+1)
                    do 530 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
     +                  + n2r(ib+ibase)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)
 530                continue
                    kont = kont + ntrm
 540              continue
                end if
                if (max3.ge.min3) then
c... generate type 3 points
                  do 560 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    ntrm = mtrm2(ib-1)
                    ntrm3(ib) = ntrm
                    ilif = jlif2(ib-1)
                    do 550 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)
 550                continue
                    kont = kont + ntrm
 560              continue
                end if
              end if
              call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +                   k,l)
            else
c...
c...    ***case 3***
c... first complete exchange term
c...
c...
              call midlhd(j,l)
c...
c... transition from lh double to lh single at level l
c...   see   ****table 15****
c...
              ibase = i16(l)
              call min23(l)
              kont = 0
              if (npopr(l).ne.0) then
c...
c... occupations  lhs=2  rhs=1
c...
                max6p1 = max6 + 1
                max7p1 = max7 + 1
                min7m1 = min7 - 1
                min8m1 = min8 - 1
                if (max2.ge.min2) then
c... generate type 4 points
                  do 590 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    if (ib.gt.min6 .and. ib.le.max6p1) then
                      ntrm = ntrm6(ib-1)
                      ilif = ilif6(ib-1)
                      do 570 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*w(ib)
 570                  continue
                      kont = kont + ntrm
                    end if
                    if (ib.lt.max7 .and. ib.ge.min7m1) then
                      ntrm = ntrm7(ib+1)
                      ilif = ilif7(ib+1)
                      do 580 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
     +                    + n2r(ib+ibase)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*v2neg(ib)
     +                                 + yj(ilif+loop)*w2neg(ib+1)
 580                  continue
                      kont = kont + ntrm
                    end if
                    ntrm2(ib) = kont - ilif2(ib)
 590              continue
                end if
                if (max3.ge.min3) then
c... generate type 5 points
                  do 620 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    if (ib.gt.min7 .and. ib.le.max7p1) then
                      ntrm = ntrm7(ib-1)
                      ilif = ilif7(ib-1)
                      do 600 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*w2(ib)
     +                                 + yj(ilif+loop)*v2neg(ib-1)
 600                  continue
                      kont = kont + ntrm
                    end if
                    if (ib.lt.max8 .and. ib.ge.min8m1) then
                      ntrm = ntrm8(ib+1)
                      ilif = ilif8(ib+1)
                      do 610 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
     +                    + n2r(ib+ibase)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*v(ib)
 610                  continue
                      kont = kont + ntrm
                    end if
                    ntrm3(ib) = kont - ilif3(ib)
 620              continue
                end if
              else
c...
c... occupations  lhs=1  rhs=0
c...
                if (max2.ge.min2) then
c... generate type 4 points
                  do 650 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    if (ib.ge.min6 .and. ib.le.max6) then
                      ntrm = ntrm6(ib)
                      ilif = ilif6(ib)
                      do 630 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
     +                    + n2l(ib+ibase+1)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*aneg(ib+1)
 630                  continue
                      kont = kont + ntrm
                    end if
                    if (ib.ge.min7 .and. ib.le.max7) then
                      ntrm = ntrm7(ib)
                      ilif = ilif7(ib)
                      do 640 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*factj2(ib+1)
     +                                 + yj(ilif+loop)*tneg
 640                  continue
                      kont = kont + ntrm
                    end if
                    ntrm2(ib) = kont - ilif2(ib)
 650              continue
                end if
                if (max3.ge.min3) then
c... generate type 5 points
                  do 680 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    if (ib.ge.min7 .and. ib.le.max7) then
                      ntrm = ntrm7(ib)
                      ilif = ilif7(ib)
                      do 660 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
     +                    + n2l(ib+ibase-1)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*factk2(ib-1)
     +                                 + yj(ilif+loop)*tneg
 660                  continue
                      kont = kont + ntrm
                    end if
                    if (ib.ge.min8 .and. ib.le.max8) then
                      ntrm = ntrm8(ib)
                      ilif = ilif8(ib)
                      do 670 loop = 1 , ntrm
                        lexlhs(kont+loop) = jexlhs(ilif+loop)
                        lexrhs(kont+loop) = jexrhs(ilif+loop)
                        mrange(kont+loop) = jrange(ilif+loop)
                        x(kont+loop) = xj(ilif+loop)*aneg(ib-1)
 670                  continue
                      kont = kont + ntrm
                    end if
                    ntrm3(ib) = kont - ilif3(ib)
 680              continue
                end if
              end if
              call lhgen(mexlhs,mexrhs,nrange,xp,l,kounte)
c...
c... now for direct term
c...
              call minmax(j)
              call point1(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,j,l,junke)
c...
c... top segment of lh single range at level l
c...   see   ****table 12****
c...
              call min23(l)
              kont = 0
              if (npopr(l).ne.0) then
c...
c... occupations  lhs=2  rhs=1
c...
                if (min2.lt.mmin) min2 = mmin - 1
                if (min3.le.mmin) min3 = mmin + 1
                if (max2.ge.mmax) max2 = mmax - 1
                if (max3.gt.mmax) max3 = mmax + 1
                if (max2.ge.min2) then
c... generate type 4 points
                  do 700 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    ntrm = mtrm2(ib+1)
                    ntrm2(ib) = ntrm
                    ilif = jlif2(ib+1)
                    do 690 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
     +                  + n2r(ib+ibase)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)*w(ib+1)
 690                continue
                    kont = kont + ntrm
 700              continue
                end if
                if (max3.ge.min3) then
c... generate type 5 points
                  do 720 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    ntrm = mtrm2(ib-1)
                    ntrm3(ib) = ntrm
                    ilif = jlif2(ib-1)
                    do 710 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)*v(ib-1)
 710                continue
                    kont = kont + ntrm
 720              continue
                end if
              else
c...
c... occupations are lhs=1  rhs=0
c...
                if (min2.lt.mmin) min2 = mmin
                if (min3.lt.mmin) min3 = mmin
                if (max2.gt.mmax) max2 = mmax
                if (max3.gt.mmax) max3 = mmax
                if (max2.ge.min2) then
c... generate type 4 points
                  do 740 ib = min2 , max2 , 2
                    ilif2(ib) = kont
                    ntrm = mtrm2(ib)
                    ntrm2(ib) = ntrm
                    ilif = jlif2(ib)
                    do 730 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)
 730                continue
                    kont = kont + ntrm
 740              continue
                end if
                if (max3.ge.min3) then
c... generate type 5 points
                  do 760 ib = min3 , max3 , 2
                    ilif3(ib) = kont
                    ntrm = mtrm2(ib)
                    ntrm3(ib) = ntrm
                    ilif = jlif2(ib)
                    do 750 loop = 1 , ntrm
                      lexlhs(kont+loop) = kexlhs(ilif+loop)
     +                  + n2l(ib+ibase-1)
                      lexrhs(kont+loop) = kexrhs(ilif+loop)
                      mrange(kont+loop) = krange(ilif+loop)
                      x(kont+loop) = xk(ilif+loop)
 750                continue
                    kont = kont + ntrm
 760              continue
                end if
              end if
              call lhgen(mexlhs(kounte+1),mexrhs(kounte+1),
     +                   nrange(kounte+1),xp(kounte+1),l,kount)
              return
            end if
          end if
          call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +               mexlhs(kounte+1),mexrhs(kounte+1),nrange(kounte+1),
     +               xp(kounte+1),mtrm2,jlif2,l,kount)
          call point1(mexlhs(kounte+1),mexrhs(kounte+1),nrange(kounte+1)
     +                ,xp(kounte+1),mtrm2,jlif2,l,norb1,kount)
          return
        else
c...
c... unorthodox top segments of lh single range at level i
c... see       *****table 16*****
c...
c...
          ibase = i16(i)
          call min23(i)
          kont = 0
          if (npopl(i).ne.0) then
c...
c... occupations lhs=1  rhs=2
c...
            if (max2.ge.min2) then
c... generate type 4 points
              do 770 loop = min2 , max2 , 2
                kont = kont + 1
                lexlhs(kont) = 0
                lexrhs(kont) = 0
                x(kont) = oneneg
                ntrm2(loop) = 1
                ilif2(loop) = numb(kont)
                mrange(kont) = n1r(loop+ibase)
 770          continue
            end if
            if (max3.ge.min3) then
c... generate type 5 points
              do 780 loop = min3 , max3 , 2
                kont = kont + 1
                ilif3(loop) = numb(kont)
                ntrm3(loop) = 1
                lexrhs(kont) = 0
                x(kont) = oneneg
                lexlhs(kont) = n2l(loop+ibase-1)
                mrange(kont) = n1l(loop+ibase-1)
 780          continue
            end if
          else
c...
c... occupations lhs=0 rhs=1
c...
            if (max2.ge.min2) then
c... generate type 4 points
              do 790 loop = min2 , max2 , 2
                kont = kont + 1
                ilif2(loop) = numb(kont)
                ntrm2(loop) = 1
                lexlhs(kont) = 0
                x(kont) = w(loop+1)
                lexrhs(kont) = n2r(loop+ibase)
                mrange(kont) = n1r(loop+ibase)
 790          continue
            end if
            if (max3.ge.min3) then
c... generate type 5 points
              do 800 loop = min3 , max3 , 2
                kont = kont + 1
                ilif3(loop) = numb(kont)
                ntrm3(loop) = 1
                lexlhs(kont) = 0
                lexrhs(kont) = 0
                mrange(kont) = n2r(loop+ibase)
                x(kont) = v(loop-1)
 800          continue
            end if
          end if
c...
c... complete lh range
c...
          call lhgen(mexlhs,mexrhs,nrange,xp,i,kount)
          return
        end if
c...
c...
c... unorthodox entry to rh double range
c...
c...
      else if ((npopl(j)-npopr(j)).ne.2) then
c...
c...
c... unorthodox top segment of rh single range
c... at level j   see    *****table 8*****
c...
c...
        ibase = i16(j)
        call min23(j)
        kont = 0
        if (npopr(j).ne.0) then
c...
c... occupations are lhs=2  rhs=1
c...
          if (max2.ge.min2) then
c... generate type 2 points
            do 810 loop = min2 , max2 , 2
              kont = kont + 1
              lexlhs(kont) = 0
              x(kont) = oneneg
              ntrm2(loop) = 1
              ilif2(loop) = numb(kont)
              mrange(kont) = n1r(loop+ibase)
              lexrhs(kont) = n2r(loop+ibase)
 810        continue
          end if
          if (max3.ge.min3) then
c... generate type 3 points
            do 820 loop = min3 , max3 , 2
              kont = kont + 1
              ilif3(loop) = numb(kont)
              ntrm3(loop) = 1
              lexrhs(kont) = 0
              lexlhs(kont) = 0
              x(kont) = oneneg
              mrange(kont) = n2r(loop+ibase)
 820        continue
          end if
        else
c...
c... occupations are  lhs=1  rhs=0
c...
          if (max2.ge.min2) then
c... generate type 2 points
            do 830 loop = min2 , max2 , 2
              kont = kont + 1
              lexlhs(kont) = 0
              lexrhs(kont) = 0
              ntrm2(loop) = 1
              ilif2(loop) = numb(kont)
              x(kont) = v(loop)
              mrange(kont) = n1r(loop+ibase)
 830        continue
          end if
          if (max3.ge.min3) then
c... generate type 3 points
            do 840 loop = min3 , max3 , 2
              kont = kont + 1
              ilif3(loop) = numb(kont)
              ntrm3(loop) = 1
              lexrhs(kont) = 0
              x(kont) = w(loop)
              lexlhs(kont) = n2l(ibase+loop-1)
              mrange(kont) = n1l(ibase+loop-1)
 840        continue
          end if
        end if
        call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,j,l)
        call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +             mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,l,kount)
        call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,l,norb1,kount)
        return
      else
c...
c...
c... unorthodox entry to and exit from rh double range
c...
c...
        kount = 1
        nrange(1) = n1l(norb16+1)
        mexlhs(1) = 0
        mexrhs(1) = 0
        xp(1) = 1.0d0
        return
      end if
      end
**==lhdrh.f
      subroutine lhdrh(k)
c...   see   ****table 11****
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/junk/llhs(mxcan1),lrhs(mxcan1),lr(mxcan1),x(mxcan1),
     *i2(16),n2(16),i3(16),n3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jlhs(mxcan1),jrhs(mxcan1),jr(mxcan1),xj(mxcan1),yj(mxcan1),
     *i6(16),i7(16),i8(16),n6(16),n7(16),n8(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,
     *min8,max8
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      kont = 0
      ibase = i16(k)
      call min23(k)
      if (npopl(k).ne.0) then
c...
c... occupations lhs=1  rhs=2
c...
        if (max2.ge.min2) then
c... generate type 2 points
          do 40 ib = min2 , max2 , 2
            i2(ib) = kont
            if (ib.ge.min6 .and. ib.le.max6) then
              ntrm = n6(ib)
              ilif = i6(ib)
              do 20 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop) + n2l(ib+ibase+1)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*v(ib+1)
 20           continue
              kont = kont + ntrm
            end if
            if (ib.ge.min7 .and. ib.le.max7) then
              ntrm = n7(ib)
              ilif = i7(ib)
              do 30 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*w2(ib+1) + yj(ilif+loop)
     +                         *v2neg(ib)
 30           continue
              kont = kont + ntrm
            end if
            n2(ib) = kont - i2(ib)
 40       continue
        end if
        if (max3.ge.min3) then
c... generate type 3 points
          do 70 ib = min3 , max3 , 2
            i3(ib) = kont
            if (ib.ge.min7 .and. ib.le.max7) then
              ntrm = n7(ib)
              ilif = i7(ib)
              do 50 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop) + n2l(ib+ibase-1)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*v2neg(ib-1) + yj(ilif+loop)
     +                         *w2neg(ib)
 50           continue
              kont = kont + ntrm
            end if
            if (ib.ge.min8 .and. ib.le.max8) then
              ntrm = n8(ib)
              ilif = i8(ib)
              do 60 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*w(ib-1)
 60           continue
              kont = kont + ntrm
            end if
            n3(ib) = kont - i3(ib)
 70       continue
        end if
      else
c...
c... occupations are lhs=0 rhs=1
c...
        max6p1 = max6 + 1
        max7p1 = max7 + 1
        min7m1 = min7 - 1
        min8m1 = min8 - 1
        if (max2.ge.min2) then
c... generate type 2 points
          do 100 ib = min2 , max2 , 2
            i2(ib) = kont
            if (ib.gt.min6 .and. ib.le.max6p1) then
              ntrm = n6(ib-1)
              ilif = i6(ib-1)
              do 80 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*aneg(ib)
 80           continue
              kont = kont + ntrm
            end if
            if (ib.lt.max7 .and. ib.ge.min7m1) then
              ntrm = n7(ib+1)
              ilif = i7(ib+1)
              do 90 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop) + n2r(ib+ibase)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*factk2(ib) + yj(ilif+loop)
     +                         *tneg
 90           continue
              kont = kont + ntrm
            end if
            n2(ib) = kont - i2(ib)
 100      continue
        end if
        if (max3.ge.min3) then
c... generate type 3 points
          do 130 ib = min3 , max3 , 2
            i3(ib) = kont
            if (ib.gt.min7 .and. ib.le.max7p1) then
              ntrm = n7(ib-1)
              ilif = i7(ib-1)
              do 110 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*factj2(ib) + yj(ilif+loop)
     +                         *tneg
 110          continue
              kont = kont + ntrm
            end if
            if (ib.lt.max8 .and. ib.ge.min8m1) then
              ntrm = n8(ib+1)
              ilif = i8(ib+1)
              do 120 loop = 1 , ntrm
                llhs(kont+loop) = jlhs(ilif+loop)
                lrhs(kont+loop) = jrhs(ilif+loop) + n2r(ib+ibase)
                lr(kont+loop) = jr(ilif+loop)
                x(kont+loop) = xj(ilif+loop)*aneg(ib)
 120          continue
              kont = kont + ntrm
            end if
            n3(ib) = kont - i3(ib)
 130      continue
          return
        end if
      end if
      return
      end
**==rhlhd.f
      subroutine rhlhd
c...   see   ****table 9****
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),
     *jlif2(16),mtrm2(16),jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif6(16),ilif7(16),ilif8(16),ntrm6(16),ntrm7(16),ntrm8(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kkonte(2),
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,
     *min8,max8
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      ibase = i16(j)
      call min678(j)
      kounte = 0
      if (npopr(j).ne.0) then
c...
c... occupations  lhs=2  rhs=1
c...
        if (min6.lt.min2) min6 = min2 - 1
        if (max6.ge.max2) max6 = max2 - 1
        if (min8.le.min3) min8 = min3 + 1
        if (max8.gt.max3) max8 = max3 + 1
        if (max6.ge.min6) then
c... generate type 9 positions
          do 30 ib = min6 , max6 , 2
            ilif6(ib) = kounte
            ntrm = ntrm2(ib+1)
            ntrm6(ib) = ntrm
            ilif = ilif2(ib+1)
            do 20 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop) + n2r(ib+ibase)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*wneg(ib+2)
 20         continue
            kounte = kounte + ntrm
 30       continue
        end if
c... generate type 10 positions
        max2p1 = max2 + 1
        min3m1 = min3 - 1
c... generate type 10 points
        do 60 ib = min7 , max7 , 2
          ilif7(ib) = kounte
          if (ib.gt.min2 .and. ib.le.max2p1) then
c... coupling lhs=3 rhs=2
            ntrm = ntrm2(ib-1)
            ilif = ilif2(ib-1)
            do 40 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*factm2(ib)
              yj(kounte+loop) = x(ilif+loop)*w2(ib)
 40         continue
            kounte = kounte + ntrm
          end if
          if (ib.ge.min3m1 .and. ib.lt.max3) then
c... coupling  lhs=3  rhs=1
            ntrm = ntrm3(ib+1)
            ilif = ilif3(ib+1)
            do 50 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop) + n2r(ib+ibase)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*factl2(ib)
              yj(kounte+loop) = x(ilif+loop)*v2(ib)
 50         continue
            kounte = kounte + ntrm
          end if
          ntrm7(ib) = kounte - ilif7(ib)
 60     continue
        if (max8.ge.min8) then
c... generate type 11 points
          do 80 ib = min8 , max8 , 2
            ilif8(ib) = kounte
            ntrm = ntrm3(ib-1)
            ntrm8(ib) = ntrm
            ilif = ilif3(ib-1)
            do 70 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*vneg(ib-2)
 70         continue
            kounte = kounte + ntrm
 80       continue
        end if
      else
c...
c... occupations  lhs=1  rhs=0
c...
        if (min6.lt.min2) min6 = min2
        if (max6.gt.max2) max6 = max2
        if (min8.lt.min3) min8 = min3
        if (max8.gt.max3) max8 = max3
        if (max6.ge.min6) then
c... generate type 9 points
          do 100 ib = min6 , max6 , 2
            ilif6(ib) = kounte
            ntrm = ntrm2(ib)
            ntrm6(ib) = ntrm
            ilif = ilif2(ib)
            do 90 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)
 90         continue
            kounte = kounte + ntrm
 100      continue
        end if
c... generate type 10 points
        do 130 ib = min7 , max7 , 2
          ilif7(ib) = kounte
          if (ib.ge.min2 .and. ib.le.max2) then
            ntrm = ntrm2(ib)
            ilif = ilif2(ib)
            do 110 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop) + n2l(ib+ibase)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*factj2(ib)
              yj(kounte+loop) = x(ilif+loop)*t
 110        continue
            kounte = kounte + ntrm
          end if
          if (ib.ge.min3 .and. ib.le.max3) then
            ntrm = ntrm3(ib)
            ilif = ilif3(ib)
            do 120 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)*factk2(ib)
              yj(kounte+loop) = x(ilif+loop)*t
 120        continue
            kounte = kounte + ntrm
          end if
          ntrm7(ib) = kounte - ilif7(ib)
 130    continue
        if (max8.ge.min8) then
c... generate type 11 positions
          do 150 ib = min8 , max8 , 2
            ilif8(ib) = kounte
            ntrm = ntrm3(ib)
            ntrm8(ib) = ntrm
            ilif = ilif3(ib)
            do 140 loop = 1 , ntrm
              jexlhs(kounte+loop) = lexlhs(ilif+loop) + n2l(ib+ibase-2)
              jexrhs(kounte+loop) = lexrhs(ilif+loop)
              jrange(kounte+loop) = mrange(ilif+loop)
              xj(kounte+loop) = x(ilif+loop)
 140        continue
            kounte = kounte + ntrm
 150      continue
          return
        end if
      end if
      return
      end
**==midlhd.f
      subroutine midlhd(j,k)
c...   see   ****table 10****
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *ilhs(mxcan1),irhs(mxcan1),ir(mxcan1),xi(mxcan1),yi(mxcan1),
     *il9(16),il10(16),il11(16),mt9(16),mt10(16),mt11(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/lsort/jlhs(mxcan1),jrhs(mxcan1),jr(mxcan1),
     *xj(mxcan1),
     *yj(mxcan1),nt9(16),jl9(16)
      common/moco/nt10(16),jl10(16),nt11(16),jl11(16)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),
     *mmin,mmax,min2,max2,min3,max3,min9,max9,min10,max10,
     *min11,max11
INCLUDE(common/icon)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      jp1 = j + 1
      km1 = k - 1
      if (km1.ge.jp1) then
c...
c...
c... recur through middle segments of lh double range
c...
c...
        do 190 kk = jp1 , km1
          if (npopr(kk).eq.1) then
c... skip over empty or doubly occupied levels
            max9m1 = max9 - 1
            min9m1 = min9 - 1
            mx10m1 = max10 - 1
            mn10m1 = min10 - 1
            mx11m1 = max11 - 1
            mn11m1 = min11 - 1
            max9p1 = max9 + 1
            min9p1 = min9 + 1
            mx10p1 = max10 + 1
            mn10p1 = min10 + 1
            mx11p1 = max11 + 1
            mn11p1 = min11 + 1
            ibase = i16(kk)
            call min678(kk)
            kount = 0
            if (max9.ge.min9) then
c...
c... generate type 9 points
c...
              do 50 ib = min9 , max9 , 2
                jl9(ib) = kount
                n2rib = n2r(ib+ibase)
                if (ib.le.max9m1 .and. ib.ge.min9m1) then
c... coupling  lhs=1  rhs=1
                  ntrm = mt9(ib+1)
                  ilif = il9(ib+1)
                  do 20 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop) + n2l(ib+ibase+2)
                    jrhs(kount+loop) = irhs(ilif+loop) + n2rib
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*a(ib+2)
 20               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min9p1 .and. ib.le.max9p1) then
c... coupling lhs=2 rhs=2
                  ntrm = mt9(ib-1)
                  ilif = il9(ib-1)
                  do 30 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop)
                    jrhs(kount+loop) = irhs(ilif+loop)
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*a(ib)
 30               continue
                  kount = kount + ntrm
                end if
                if (ib.le.mx10m1 .and. ib.ge.mn10m1) then
c... coupling lhs=2 rhs=1
                  ntrm = mt10(ib+1)
                  ilif = il10(ib+1)
                  do 40 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop)
                    jrhs(kount+loop) = irhs(ilif+loop) + n2rib
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*g2(ib)
 40               continue
                  kount = kount + ntrm
                end if
                nt9(ib) = kount - jl9(ib)
 50           continue
            end if
            if (max10.ge.min10) then
c...
c... generate type 10 points
c...
              do 100 ib = min10 , max10 , 2
                jl10(ib) = kount
                n2lib = n2l(ib+ibase)
                n2rib = n2r(ib+ibase)
                cib = c(ib)
                eib = e(ib)
                if (ib.le.mx10m1 .and. ib.ge.mn10m1) then
c... coupling lhs=1 rhs=1
                  ntrm = mt10(ib+1)
                  ilif = il10(ib+1)
                  do 60 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop) + n2lib
                    jrhs(kount+loop) = irhs(ilif+loop) + n2rib
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*cib
                    yj(kount+loop) = yi(ilif+loop)
 60               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.mn10p1 .and. ib.le.mx10p1) then
c... coupling  lhs=2  rhs=2
                  ntrm = mt10(ib-1)
                  ilif = il10(ib-1)
                  do 70 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop)
                    jrhs(kount+loop) = irhs(ilif+loop)
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*cib
                    yj(kount+loop) = yi(ilif+loop)
 70               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min9p1 .and. ib.le.max9p1) then
c... coupling lhs=1  rhs=2
                  ntrm = mt9(ib-1)
                  ilif = il9(ib-1)
                  do 80 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop) + n2lib
                    jrhs(kount+loop) = irhs(ilif+loop)
                    jr(kount+loop) = ir(ilif+loop)
                    yj(kount+loop) = 0.0d0
                    xj(kount+loop) = xi(ilif+loop)*eib
 80               continue
                  kount = kount + ntrm
                end if
                if (ib.le.mx11m1 .and. ib.ge.mn11m1) then
c... coupling  lhs=2  rhs=1
                  ntrm = mt11(ib+1)
                  ilif = il11(ib+1)
                  do 90 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop)
                    jrhs(kount+loop) = irhs(ilif+loop) + n2rib
                    jr(kount+loop) = ir(ilif+loop)
                    yj(kount+loop) = 0.0d0
                    xj(kount+loop) = xi(ilif+loop)*eib
 90               continue
                  kount = kount + ntrm
                end if
                nt10(ib) = kount - jl10(ib)
 100          continue
            end if
            if (max11.ge.min11) then
c...
c... generate type 11 points
c...
              do 140 ib = min11 , max11 , 2
                jl11(ib) = kount
                n2lib = n2l(ib+ibase-2)
                if (ib.le.mx11m1 .and. ib.ge.mn11m1) then
c... coupling lhs=1  rhs=1
                  ntrm = mt11(ib+1)
                  ilif = il11(ib+1)
                  do 110 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop) + n2lib
                    jrhs(kount+loop) = irhs(ilif+loop) + n2r(ib+ibase)
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*a(ib)
 110              continue
                  kount = kount + ntrm
                end if
                if (ib.ge.mn11p1 .and. ib.le.mx11p1) then
c... coupling lhs=2 rhs=2
                  ntrm = mt11(ib-1)
                  ilif = il11(ib-1)
                  do 120 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop)
                    jrhs(kount+loop) = irhs(ilif+loop)
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*a(ib-2)
 120              continue
                  kount = kount + ntrm
                end if
                if (ib.ge.mn10p1 .and. ib.le.mx10p1) then
c... coupling lhs=1 rhs=2
                  ntrm = mt10(ib-1)
                  ilif = il10(ib-1)
                  do 130 loop = 1 , ntrm
                    jlhs(kount+loop) = ilhs(ilif+loop) + n2lib
                    jrhs(kount+loop) = irhs(ilif+loop)
                    jr(kount+loop) = ir(ilif+loop)
                    xj(kount+loop) = xi(ilif+loop)*g2(ib-2)
 130              continue
                  kount = kount + ntrm
                end if
                nt11(ib) = kount - jl11(ib)
 140          continue
c... shuffle generated data back to correct locations
              do 150 ib = min11 , max11 , 2
                mt11(ib) = nt11(ib)
                il11(ib) = jl11(ib)
 150          continue
            end if
            if (max10.ge.min10) then
              do 160 ib = min10 , max10 , 2
                mt10(ib) = nt10(ib)
                il10(ib) = jl10(ib)
 160          continue
            end if
            if (max9.ge.min9) then
              do 170 ib = min9 , max9 , 2
                mt9(ib) = nt9(ib)
                il9(ib) = jl9(ib)
 170          continue
            end if
            do 180 loop = 1 , kount
              ilhs(loop) = jlhs(loop)
              irhs(loop) = jrhs(loop)
              ir(loop) = jr(loop)
              yi(loop) = yj(loop)
              xi(loop) = xj(loop)
 180        continue
          end if
 190    continue
      end if
      return
      end
**==lhgen.f
      subroutine lhgen(mexlhs,mexrhs,nrange,xp,ll,kont)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension mexlhs(*),mexrhs(*),nrange(*),xp(*)
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),ilif2(16),ntrm2(16),ilif3(16),ntrm3(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kkonte(2),
     *mmin,mmax,min2,max2,min3,max3
INCLUDE(common/icon)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      common/moco/jlif2(16),mtrm2(16),jlif3(16),mtrm3(16)
c... recur through lh single range
c...   see   ****table 13****
      lp1 = ll + 1
      km1 = k - 1
      if (km1.ge.lp1) then
        do 130 kk = lp1 , km1
          if (npopr(kk).eq.1) then
c... skip over empty or doubly occupied levels
            max3m1 = max3 - 1
            min3m1 = min3 - 1
            max2m1 = max2 - 1
            min2m1 = min2 - 1
            max3p1 = max3 + 1
            min3p1 = min3 + 1
            max2p1 = max2 + 1
            min2p1 = min2 + 1
            ibase = i16(kk)
            call min23(kk)
            kount = 0
            if (max2.ge.min2) then
c...
c... generate type 4 points
c...
              do 50 ib = min2 , max2 , 2
                jlif2(ib) = kount
                nib = n2r(ib+ibase)
                if (ib.le.max2m1 .and. ib.ge.min2m1) then
c... coupling  lhs=1  rhs=1
                  ntrm = ntrm2(ib+1)
                  ilif = ilif2(ib+1)
                  do 20 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
     +                + n2l(ib+ibase+1)
                    mexrhs(kount+loop) = lexrhs(ilif+loop) + nib
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*aneg(ib+1)
 20               continue
                  kount = kount + ntrm
                end if
                if (ib.le.max3m1 .and. ib.ge.min3m1) then
c... coupling is  lhs=2  rhs=1
                  ntrm = ntrm3(ib+1)
                  ilif = ilif3(ib+1)
                  do 30 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop) + nib
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*bneg(ib+1)
 30               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min2p1 .and. ib.le.max2p1) then
c... coupling  lhs=2  rhs=2
                  ntrm = ntrm2(ib-1)
                  ilif = ilif2(ib-1)
                  do 40 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = -x(ilif+loop)
 40               continue
                  kount = kount + ntrm
                end if
                mtrm2(ib) = kount - jlif2(ib)
 50           continue
            end if
            if (max3.ge.min3) then
c...
c... generate type 5 points
c...
              do 90 ib = min3 , max3 , 2
                jlif3(ib) = kount
                nib = n2l(ib+ibase-1)
                if (ib.le.max3m1 .and. ib.ge.min3m1) then
c... coupling  lhs=1  rhs=1
                  ntrm = ntrm3(ib+1)
                  ilif = ilif3(ib+1)
                  do 60 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop) + nib
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
     +                + n2r(ib+ibase)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = -x(ilif+loop)
 60               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min2p1 .and. ib.le.max2p1) then
c... coupling  lhs=1  rhs=2
                  ntrm = ntrm2(ib-1)
                  ilif = ilif2(ib-1)
                  do 70 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop) + nib
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*b(ib-1)
 70               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min3p1 .and. ib.le.max3p1) then
c... coupling lhs=2  rhs=2
                  ntrm = ntrm3(ib-1)
                  ilif = ilif3(ib-1)
                  do 80 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*aneg(ib-1)
 80               continue
                  kount = kount + ntrm
                end if
                mtrm3(ib) = kount - jlif3(ib)
 90           continue
c... shuffle data to correct location
              do 100 ib = min3 , max3 , 2
                ilif3(ib) = jlif3(ib)
                ntrm3(ib) = mtrm3(ib)
 100          continue
            end if
            if (max2.ge.min2) then
              do 110 ib = min2 , max2 , 2
                ilif2(ib) = jlif2(ib)
                ntrm2(ib) = mtrm2(ib)
 110          continue
            end if
            do 120 loop = 1 , kount
              lexlhs(loop) = mexlhs(loop)
              lexrhs(loop) = mexrhs(loop)
              mrange(loop) = nrange(loop)
              x(loop) = xp(loop)
 120        continue
          end if
 130    continue
      end if
c...
c... bottom segment of lh single range at level k
c...   see   ****table 14****
c...
      ibase = i16(k)
      call minmax(k)
      kount = 0
      if (npopl(k).ne.0) then
c...
c... occupations lhs=1 rhs=2
c...
        do 160 ib = mmin , mmax , 2
          jlif2(ib) = kount
          if (ib.ge.min2 .and. ib.le.max2) then
c... coupling  lhs=1  rhs=3
            ntrm = ntrm2(ib)
            ilif = ilif2(ib)
            do 140 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop) + n2l(ib+ibase)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)*v(ib)
 140        continue
            kount = kount + ntrm
          end if
          if (ib.ge.min3 .and. ib.le.max3) then
c... coupling lhs=2 rhs=3
            ntrm = ntrm3(ib)
            ilif = ilif3(ib)
            do 150 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)*w(ib)
 150        continue
            kount = kount + ntrm
          end if
          mtrm2(ib) = kount - jlif2(ib)
 160    continue
      else
c...
c... occupations are  lhs=0 rhs=1
c...
        max2p1 = max2 + 1
        min3m1 = min3 - 1
        do 190 ib = mmin , mmax , 2
          jlif2(ib) = kount
          if (ib.gt.min2 .and. ib.le.max2p1) then
c... coupling  lhs=0  rhs=2
            ntrm = ntrm2(ib-1)
            ilif = ilif2(ib-1)
            do 170 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)
 170        continue
            kount = kount + ntrm
          end if
          if (ib.ge.min3m1 .and. ib.lt.max3) then
c... coupling lhs=0  rhs=1
            ntrm = ntrm3(ib+1)
            ilif = ilif3(ib+1)
            do 180 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop) + n2r(ibase+ib)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)
 180        continue
            kount = kount + ntrm
          end if
          mtrm2(ib) = kount - jlif2(ib)
 190    continue
      end if
      kont = kount
      call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,k,norb1,kont)
      return
      end
**==min678.f
      subroutine min678(j)
      implicit REAL  (a-h,o-z),integer (i-n)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,
     *min8,max8
      maxl = maxlhs(j)
      maxr = maxrhs(j)
      minl = minlhs(j)
      minr = minrhs(j)
      max6 = min(maxr,maxl-2)
      max7 = min(maxr,maxl)
      max8 = min(maxr,maxl+2)
      min6 = max(minr,minl-2)
      min7 = max(minr,minl)
      min8 = max(minr,minl+2)
      return
      end
**==symb21.f
      subroutine symb21(mexlhs,mexrhs)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension mexlhs(*),mexrhs(*)
      common/scra  /mexlh(mxcan2),mexrh(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif6(16),ilif7(16),ilif8(16),ntrm6(16),ntrm7(16),ntrm8(16)
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,
     *min8,max8
      kount = 0
      kounte = 0
c****
c**** single orthogonality
c**** compute spin coupling coefficient for
c**** <i j/k k>=direct form and <i k/j k>=exchange form
c**** where i lt j
c****
      if (k.eq.i) then
        if (npopr(i).eq.2) go to 90
      else if (k.eq.j) then
        if (npopl(j).eq.2) go to 90
      else if (npopl(k).lt.1) then
      else if (npopl(k).eq.1) then
c...
c... level k singly occupied
c...
        if (k.lt.i) then
c...
c... *****case 1*****
c...
c... exchange term
          call toplhd(k)
          call midlhd(k,i)
          call lhdrh(i)
          call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,j)
          call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +               mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,j,kounte)
          call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,j,norb1,
     +                kounte)
c... direct term
          go to 90
        else
          call toprh
          if (k.lt.j) then
c...
c... *****case 2*****
c...
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,
     +                 k)
c...
c...  complete exchange term
c...
c... rh single to rh single transition
c... logic in ****table 2****
c... but coupling coeffs.  unity except for coupling modes
c... t2-t2 via 22 and t3-t3 via 11 which are zero
c...
            max3m1 = max3 - 1
            min3m1 = min3 - 1
            max2m1 = max2 - 1
            min2m1 = min2 - 1
            max3p1 = max3 + 1
            min3p1 = min3 + 1
            max2p1 = max2 + 1
            min2p1 = min2 + 1
            ibase = i16(k)
            call min23(k)
            kont = 0
            if (max2.gt.max3m1) max2 = max3m1
            if (min2.lt.min2m1) min2 = min2m1
            if (max3.gt.max3p1) max3 = max3p1
            if (min3.lt.min2p1) min3 = min2p1
            if (max2.ge.min2) then
c... generate type 2 points
              do 40 ib = min2 , max2 , 2
                jlif2(ib) = kont
                nib = n2r(ib+ibase)
                if (ib.le.max2m1) then
c... coupling lhs=1 rhs=1
                  ntrm = ntrm2(ib+1)
                  ilif = ilif2(ib+1)
                  do 20 loop = 1 , ntrm
                    kexlhs(kont+loop) = lexlhs(ilif+loop)
     +                                  + n2l(ib+ibase+1)
                    kexrhs(kont+loop) = lexrhs(ilif+loop) + nib
                    krange(kont+loop) = mrange(ilif+loop)
                    xk(kont+loop) = x(ilif+loop)
 20               continue
                  kont = kont + ntrm
                end if
                if (ib.ge.min3m1) then
c... coupling lhs=2 rhs=1
                  ntrm = ntrm3(ib+1)
                  ilif = ilif3(ib+1)
                  do 30 loop = 1 , ntrm
                    kexlhs(kont+loop) = lexlhs(ilif+loop)
                    kexrhs(kont+loop) = lexrhs(ilif+loop) + nib
                    krange(kont+loop) = mrange(ilif+loop)
                    xk(kont+loop) = x(ilif+loop)
 30               continue
                  kont = kont + ntrm
                end if
                mtrm2(ib) = kont - jlif2(ib)
 40           continue
            end if
            if (max3.ge.min3) then
c... generate type 3 points
              do 70 ib = min3 , max3 , 2
                jlif3(ib) = kont
                if (ib.le.max2p1) then
c... coupling lhs=1  rhs=2
                  ntrm = ntrm2(ib-1)
                  ilif = ilif2(ib-1)
                  do 50 loop = 1 , ntrm
                    kexlhs(kont+loop) = lexlhs(ilif+loop)
     +                                  + n2l(ib+ibase-1)
                    kexrhs(kont+loop) = lexrhs(ilif+loop)
                    krange(kont+loop) = mrange(ilif+loop)
                    xk(kont+loop) = x(ilif+loop)
 50               continue
                  kont = kont + ntrm
                end if
                if (ib.ge.min3p1) then
c... coupling lhs=2   rhs=2
                  ntrm = ntrm3(ib-1)
                  ilif = ilif3(ib-1)
                  do 60 loop = 1 , ntrm
                    kexlhs(kont+loop) = lexlhs(ilif+loop)
                    kexrhs(kont+loop) = lexrhs(ilif+loop)
                    krange(kont+loop) = mrange(ilif+loop)
                    xk(kont+loop) = x(ilif+loop)
 60               continue
                  kont = kont + ntrm
                end if
                mtrm3(ib) = kont - jlif3(ib)
 70           continue
            end if
            call midrh(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,mtrm3,jlif3,
     +                 k,j)
            call botrh(kexlhs,kexrhs,krange,xk,mtrm2,jlif2,mtrm3,jlif3,
     +                 mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,j,kounte)
            call point1(mexlhs,mexrhs,nrange,xp,ntrm6,ilif6,j,norb1,
     +                  kounte)
c... complete direct term
            km1 = k - 1
            call min23(km1)
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +                 km1,j)
          else
c...
c... *****case 3*****
c...
            call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,
     +                 j)
c... complete exchange term
            call rhlhd
            call midlhd(j,k)
            call botlhd(mexlhs,mexrhs)
c... complete direct term
            call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,k,norb1,
     +                  kounte)
          end if
          go to 100
        end if
      else
c...
c... level k doubly occupied
c...
        call toprh
        call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,j)
        call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +             mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,j,kount)
        call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,j,norb1,kount)
        kounte = kount
c... now scale for direct and exchange cases
        do 80 loop = 1 , kount
          mexlhs(kount+loop) = mexlhs(loop)
          mexrhs(kount+loop) = mexrhs(loop)
          nrange(kount+loop) = nrange(loop)
          xp(kount+loop) = xp(loop)*2.0d0
          xp(loop) = -xp(loop)
 80     continue
      end if
      return
 90   call toprh
      call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,j)
 100  call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +           mexlhs(kounte+1),mexrhs(kounte+1),nrange(kounte+1),
     +           xp(kounte+1),mtrm2,jlif2,j,kount)
      call point1(mexlhs(kounte+1),mexrhs(kounte+1),nrange(kounte+1),
     +            xp(kounte+1),mtrm2,jlif2,j,norb1,kount)
      return
      end
**==symb20.f
      subroutine symb20
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      logical ieqk
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif6(16),ilif7(16),ilif8(16),ntrm6(16),ntrm7(16),ntrm8(16)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,
     *mmin,mmax,min2,max2,min3,max3,min6,max6,min7,max7,min8,max8
      data oneneg,twoneg/-1.0d0,-2.0d0/
      ieqk = i.eq.k
      nstate = n1l(norb16+1)
      npopi = npopl(i) - 1
      npopk = npopl(k) - 1
      kount = 0
      kounte = 0
c****
c**** no orthogonalities
c**** compute spin coupling coefficient for
c**** <i i/k k>=direct form and <i k/i k>=exchange form
c**** where i le k
c****
      if (npopi.lt.0) go to 40
      if (npopi.eq.0) then
        if (npopk.lt.0) go to 40
        if (npopk.eq.0) then
c...
c... both levels singly occupied
c...
          if (ieqk) go to 40
          call toplhd(i)
          call midlhd(i,k)
          call botlhd(mexlhs,mexrhs)
c... complete exchange term
          call point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,k,norb1,
     +                kounte)
c... complete direct term
          xp(kounte+1) = 1.0d0
          go to 30
        else
          go to 50
        end if
      else
c...
c... level i doubly occupied
c...
        if (npopk.lt.0) go to 40
        if (npopk.eq.0) go to 50
c...
c... orbital i and k doubly occupied
c...
        if (ieqk) then
          xp(kounte+1) = 1.0d0
          go to 30
        else
          xp(1) = twoneg
          xp(2) = 4.0d0
        end if
      end if
 20   kounte = 1
      mexlhs(1) = 0
      mexrhs(1) = 0
      nrange(1) = nstate
 30   nrange(kounte+1) = nstate
      mexlhs(kounte+1) = 0
      mexrhs(kounte+1) = 0
      kount = 1
 40   return
c...
c... one level singly occupied one level doubly occupied
c...
 50   xp(1) = oneneg
      xp(2) = 2.0d0
      go to 20
      end
**==toplhd.f
      subroutine toplhd(i)
c...   see   ****table 17****
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif9(16),ilif10(16),ilif11(16),ntrm9(16),ntrm10(16),ntrm11(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
INCLUDE(common/icon)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),
     *mmin,mmax,min2,max2,min3,max3,min9,max9,min10,max10,
     *min11,max11
      kont = 0
      ibase = i16(i)
      call min678(i)
      if (max9.ge.min9) then
c... generate type 9 points
        do 20 ib = min9 , max9 , 2
          kont = kont + 1
          jexlhs(kont) = 0
          xj(kont) = 1.0d0
          jrange(kont) = n1r(ib+ibase)
          jexrhs(kont) = n2r(ib+ibase)
          ilif9(ib) = numb(kont)
          ntrm9(ib) = 1
 20     continue
      end if
      if (max11.ge.min11) then
c... generate type 11 points
        do 30 ib = min11 , max11 , 2
          kont = kont + 1
          jexrhs(kont) = 0
          xj(kont) = 1.0d0
          jrange(kont) = n1l(ib+ibase-2)
          jexlhs(kont) = n2l(ib+ibase-2)
          ilif11(ib) = numb(kont)
          ntrm11(ib) = 1
 30     continue
      end if
      if (max10.ge.min10) then
c... generate type 10 points
        do 40 ib = min10 , max10 , 2
          ilif10(ib) = kont
          m = n2l(ib+ibase)
          if (m.ne.0) then
c... coupling lhs=2 rhs=2
            kont = kont + 1
            jexlhs(kont) = 0
            jexrhs(kont) = 0
            xj(kont) = factk2(ib)
            yj(kont) = t
            jrange(kont) = m
          end if
          n = n1l(ib+ibase)
          if (n.ne.0) then
c... coupling lhs=1 rhs=1
            kont = kont + 1
            jexlhs(kont) = m
            jexrhs(kont) = m
            xj(kont) = factj2(ib)
            yj(kont) = t
            jrange(kont) = n
          end if
          ntrm10(ib) = kont - ilif10(ib)
 40     continue
      end if
      return
      end
**==botlhd.f
      subroutine botlhd(mexlhs,mexrhs)
c...   see   ****table 18****
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension mexlhs(*),mexrhs(*)
      common/scra  /mexlh(mxcan2),mexrh(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),jlif2(16),mtrm2(16),
     *jlif3(16),mtrm3(16),
     *kexlhs(mxcan1),kexrhs(mxcan1),krange(mxcan1),xk(mxcan1),
     *jexlhs(mxcan1),jexrhs(mxcan1),jrange(mxcan1),xj(mxcan1),
     *yj(mxcan1),
     *ilif9(16),ilif10(16),ilif11(16),ntrm9(16),ntrm10(16),ntrm11(16)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,
     *mmin,mmax,min2,max2,min3,max3,min9,max9,min10,max10,
     *min11,max11
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      max9p1 = max9 + 1
      mx10p1 = max10 + 1
      mn10m1 = min10 - 1
      mn11m1 = min11 - 1
      ibase = i16(k)
      call minmax(k)
      do 60 ib = mmin , mmax , 2
        jlif2(ib) = kounte
        anegib = aneg(ib)
        n2lib = n2l(ib+ibase)
        n2rib = n2r(ib+ibase)
        if (ib.gt.min9 .and. ib.le.max9p1) then
c... coupling lhs=1 rhs=2
          ntrm = ntrm9(ib-1)
          ilif = ilif9(ib-1)
          do 20 loop = 1 , ntrm
            mexlhs(kounte+loop) = jexlhs(ilif+loop) + n2lib
            mexrhs(kounte+loop) = jexrhs(ilif+loop)
            nrange(kounte+loop) = jrange(ilif+loop)
            xp(kounte+loop) = xj(ilif+loop)*anegib
 20       continue
          kounte = kounte + ntrm
        end if
        if (ib.lt.max10 .and. ib.ge.mn10m1) then
c... coupling lhs=1 rhs=1
          ntrm = ntrm10(ib+1)
          ilif = ilif10(ib+1)
          do 30 loop = 1 , ntrm
            mexlhs(kounte+loop) = jexlhs(ilif+loop) + n2lib
            mexrhs(kounte+loop) = jexrhs(ilif+loop) + n2rib
            nrange(kounte+loop) = jrange(ilif+loop)
            xp(kounte+loop) = xj(ilif+loop)*factk2(ib) + yj(ilif+loop)
     +                        *tneg
 30       continue
          kounte = kounte + ntrm
        end if
        if (ib.gt.min10 .and. ib.le.mx10p1) then
c... coupling lhs=2 rhs=2
          ntrm = ntrm10(ib-1)
          ilif = ilif10(ib-1)
          do 40 loop = 1 , ntrm
            mexlhs(kounte+loop) = jexlhs(ilif+loop)
            mexrhs(kounte+loop) = jexrhs(ilif+loop)
            nrange(kounte+loop) = jrange(ilif+loop)
            xp(kounte+loop) = xj(ilif+loop)*factj2(ib) + yj(ilif+loop)
     +                        *tneg
 40       continue
          kounte = kounte + ntrm
        end if
        if (ib.lt.max11 .and. ib.ge.mn11m1) then
c... coupling lhs=2 rhs=1
          ntrm = ntrm11(ib+1)
          ilif = ilif11(ib+1)
          do 50 loop = 1 , ntrm
            mexlhs(kounte+loop) = jexlhs(ilif+loop)
            mexrhs(kounte+loop) = jexrhs(ilif+loop) + n2rib
            nrange(kounte+loop) = jrange(ilif+loop)
            xp(kounte+loop) = xj(ilif+loop)*anegib
 50       continue
          kounte = kounte + ntrm
        end if
        mtrm2(ib) = kounte - jlif2(ib)
 60   continue
      return
      end
**==symb2.f
      subroutine symb2(ii,ll)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,minmax(12)
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      if (ndiff.lt.2) then
c...
c... no orthogonality
c...
        i = ii
        k = ll
        call symb20
        return
      else if (ndiff.eq.2) then
c...
c... single orthogonality
c...
        k = ii
        if (mode.lt.0) then
          call symb21(mexrhs,mexlhs)
          return
        else
          call symb21(mexlhs,mexrhs)
          return
        end if
c...
c... double orthogonality
c...
      else if (mode.lt.0) then
        call symb22(mexrhs,mexlhs)
        return
      else
        call symb22(mexlhs,mexrhs)
        return
      end if
      end
**==symb1.f
      subroutine symb1
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16),mtrm(16),jlif(16)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,mmin,mmax,min2,max2,min3,max3
c****
c**** single orthogonality
c**** compute spin coupling coefficient for <i / j>
c**** where i more occupied on rhs
c**** and   j more occupied on lhs
c**** and i lt j
c****
      call toprh
      call midrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,j)
      if (mode.lt.0) then
        call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +             mexrhs,mexlhs,nrange,xp,mtrm,jlif,j,kount)
        call point1(mexrhs,mexlhs,nrange,xp,mtrm,jlif,j,norb1,kount)
        return
      else
        call botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     +             mexlhs,mexrhs,nrange,xp,mtrm,jlif,j,kount)
        call point1(mexlhs,mexrhs,nrange,xp,mtrm,jlif,j,norb1,kount)
        return
      end if
      end
**==sym1tr.f
      subroutine sym1tr(c1,ifl,nshif1,nshif2)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension c1(*)
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kount,kounte,mmin,mmax,min2,max2,min3,max3
_IF(cray,t3e,i8)
      integer con1,con2
_ELSEIF(convex,i8drct)
      integer *8 con1,con2
_ELSE
      REAL con1,con2
_ENDIF
      common/expanc/con1(5),con2(5),
     *npop1(64),min1(64),max1(64),n1(1024),n2(1024)
      common/lsort/nfrnt2,nstat2,nfrnt1,nstat1
c...
c... translate coupling coefficients to rectangular sparse matrix
c...
      if (nshif2.eq.0) then
c... singlet/doublet/valence in 2nd configuration
        call upack2(con2(3),nfrnt2,nstat2)
        nfrnt2 = 0
      else
c... triplet in 2nd configuration
        call upack2(con2(3),nstat2,nfrnt2)
      end if
      if (nshif1.eq.0) then
c... singlet/doublet/valence in 1st configuration
        call upack2(con1(3),nfrnt1,nstat1)
        nfrnt1 = 0
      else
c... triplet in first configuration
        call upack2(con1(3),nstat1,nfrnt1)
      end if
      call transl(c1,mexrhs,mexlhs,nrange,xp,kount,ifl)
      return
      end
**==toprh.f
      subroutine toprh
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/junk/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ilif2(16),ntrm2(16),ilif3(16),ntrm3(16)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,kkonte(2),mmin,mmax,min2,max2,min3,
     *max3
c...
c...
c... top segment of rh single range
c... level i more occupied on rhs
c... see     *****table 1*****
c...
c...
      kount = 0
      ibase = i16(i)
      call min23(i)
      if (npopl(i).ne.0) then
c...
c... occupatons are lhs=1 rhs=2
c...
        if (max2.ge.min2) then
c... generate type 2 points
          do 20 loop = min2 , max2 , 2
            kount = kount + 1
            lexlhs(kount) = 0
            lexrhs(kount) = 0
            ntrm2(loop) = 1
            ilif2(loop) = numb(kount)
            x(kount) = v(loop)
            mrange(kount) = n1r(loop+ibase)
 20       continue
        end if
        if (max3.ge.min3) then
c... generate type 3 points
          do 30 loop = min3 , max3 , 2
            kount = kount + 1
            ilif3(loop) = numb(kount)
            ntrm3(loop) = 1
            lexrhs(kount) = 0
            x(kount) = w(loop)
            lexlhs(kount) = n2l(loop+ibase-1)
            mrange(kount) = n1l(loop+ibase-1)
 30       continue
        end if
      else
c...
c... occupations are lhs=0 rhs=1
c...
        if (max2.ge.min2) then
c... generate type 2 points
          do 40 loop = min2 , max2 , 2
            kount = kount + 1
            lexlhs(kount) = 0
            x(kount) = 1.0d0
            ntrm2(loop) = 1
            ilif2(loop) = numb(kount)
            mrange(kount) = n1r(loop+ibase)
            lexrhs(kount) = n2r(loop+ibase)
 40       continue
        end if
        if (max3.ge.min3) then
c... generate type 3 points
          do 50 loop = min3 , max3 , 2
            kount = kount + 1
            ilif3(loop) = numb(kount)
            ntrm3(loop) = 1
            lexrhs(kount) = 0
            lexlhs(kount) = 0
            x(kount) = 1.0d0
            mrange(kount) = n2r(loop+ibase)
 50       continue
          return
        end if
      end if
      return
      end
**==midrh.f
      subroutine midrh(
     *lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,i,j)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension ntrm2(*),ilif2(*),ntrm3(*),ilif3(*)
      dimension lexlhs(*),lexrhs(*),mrange(*),x(*)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/lsort/mexlhs(mxcan1),mexrhs(mxcan1),nrange(mxcan1),
     *xp(mxcan1),
     *jlif2(16),mtrm2(16),jlif3(16),mtrm3(16)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),mmin,mmax,min2,max2,min3,max3
c...
c...
c... loop over middle segments of rh single range
c... see     *****table 2*****
c...
c...
      ip1 = i + 1
      jm1 = j - 1
      if (jm1.ge.ip1) then
        do 130 kk = ip1 , jm1
c... skip over empty or doubly occupied levels
          if (npopl(kk).eq.1) then
            max3m1 = max3 - 1
            min3m1 = min3 - 1
            max2m1 = max2 - 1
            min2m1 = min2 - 1
            max3p1 = max3 + 1
            min3p1 = min3 + 1
            max2p1 = max2 + 1
            min2p1 = min2 + 1
            ibase = i16(kk)
            call min23(kk)
            kount = 0
            if (max2.ge.min2) then
c...
c... generate type 2 points
c...
              do 50 ib = min2 , max2 , 2
                jlif2(ib) = kount
                nib = n2r(ib+ibase)
                if (ib.le.max2m1 .and. ib.ge.min2m1) then
c... coupling  lhs=1 rhs=1
                  ntrm = ntrm2(ib+1)
                  ilif = ilif2(ib+1)
                  do 20 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
     +                + n2l(ib+ibase+1)
                    mexrhs(kount+loop) = lexrhs(ilif+loop) + nib
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = -x(ilif+loop)
 20               continue
                  kount = kount + ntrm
                end if
                if (ib.le.max3m1 .and. ib.ge.min3m1) then
c... coupling  lhs=2  rhs=1
                  ntrm = ntrm3(ib+1)
                  ilif = ilif3(ib+1)
                  do 30 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop) + nib
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*b(ib)
 30               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min2p1 .and. ib.le.max2p1) then
c... coupling  lhs=2  rhs=2
                  ntrm = ntrm2(ib-1)
                  ilif = ilif2(ib-1)
                  do 40 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*aneg(ib)
 40               continue
                  kount = kount + ntrm
                end if
                mtrm2(ib) = kount - jlif2(ib)
 50           continue
            end if
            if (max3.ge.min3) then
c...
c... generation of type 3 points
c...
              do 90 ib = min3 , max3 , 2
                jlif3(ib) = kount
                nib = n2l(ib+ibase-1)
                if (ib.le.max3m1 .and. ib.ge.min3m1) then
c... coupling  lhs=1  rhs=1
                  ntrm = ntrm3(ib+1)
                  ilif = ilif3(ib+1)
                  do 60 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop) + nib
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
     +                + n2r(ib+ibase)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*aneg(ib)
 60               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min2p1 .and. ib.le.max2p1) then
c... coupling lhs=1 rhs=2
                  ntrm = ntrm2(ib-1)
                  ilif = ilif2(ib-1)
                  do 70 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop) + nib
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = x(ilif+loop)*bneg(ib)
 70               continue
                  kount = kount + ntrm
                end if
                if (ib.ge.min3p1 .and. ib.le.max3p1) then
c... coupling lhs=2 rhs=2
                  ntrm = ntrm3(ib-1)
                  ilif = ilif3(ib-1)
                  do 80 loop = 1 , ntrm
                    mexlhs(kount+loop) = lexlhs(ilif+loop)
                    mexrhs(kount+loop) = lexrhs(ilif+loop)
                    nrange(kount+loop) = mrange(ilif+loop)
                    xp(kount+loop) = -x(ilif+loop)
 80               continue
                  kount = kount + ntrm
                end if
                mtrm3(ib) = kount - jlif3(ib)
 90           continue
c... shuffle generated data back to correct location
              do 100 ib = min3 , max3 , 2
                ilif3(ib) = jlif3(ib)
                ntrm3(ib) = mtrm3(ib)
 100          continue
            end if
            if (max2.ge.min2) then
              do 110 ib = min2 , max2 , 2
                ilif2(ib) = jlif2(ib)
                ntrm2(ib) = mtrm2(ib)
 110          continue
            end if
            do 120 loop = 1 , kount
              lexlhs(loop) = mexlhs(loop)
              lexrhs(loop) = mexrhs(loop)
              mrange(loop) = nrange(loop)
              x(loop) = xp(loop)
 120        continue
          end if
 130    continue
      end if
      return
      end
**==botrh.f
      subroutine botrh(lexlhs,lexrhs,mrange,x,ntrm2,ilif2,ntrm3,ilif3,
     *mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,j,kont)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension lexlhs(*),lexrhs(*),mrange(*),x(*),ntrm2(*),ilif2(*),
     *ntrm3(*),ilif3(*)
      dimension mexlhs(*),mexrhs(*),nrange(*),xp(*),mtrm2(*),jlif2(*)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
INCLUDE(common/icon)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),mmin,mmax,min2,max2,min3,max3
c...
c...
c...  bottom segment of rh single range
c... rhs function less occupied at  level  j
c...   see ****table 3****
c...
c...
      ibase = i16(j)
      call minmax(j)
      kount = 0
      if (npopr(j).ne.0) then
c...
c... occupations are  lhs=2 rhs=1
c...
        min3m1 = min3 - 1
        max2p1 = max2 + 1
        do 40 ib = mmin , mmax , 2
          jlif2(ib) = kount
          if (ib.lt.max3 .and. ib.ge.min3m1) then
c... coupling lhs=3 rhs=1
            ntrm = ntrm3(ib+1)
            ilif = ilif3(ib+1)
            do 20 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop) + n2r(ib+ibase)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)*v(ib)
 20         continue
            kount = kount + ntrm
          end if
          if (ib.gt.min2 .and. ib.le.max2p1) then
c... coupling  lhs=3  rhs=2
            ntrm = ntrm2(ib-1)
            ilif = ilif2(ib-1)
            do 30 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)*w(ib)
 30         continue
            kount = kount + ntrm
          end if
          mtrm2(ib) = kount - jlif2(ib)
 40     continue
      else
c...
c... occupations are  lhs=1  rhs=0
c...
        do 70 ib = mmin , mmax , 2
          jlif2(ib) = kount
          if (ib.ge.min2 .and. ib.le.max2) then
c... coupling  lhs=1  rhs=0
            ntrm = ntrm2(ib)
            ilif = ilif2(ib)
            do 50 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop) + n2l(ib+ibase)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)
 50         continue
            kount = kount + ntrm
          end if
          if (ib.ge.min3 .and. ib.le.max3) then
c... coupling   lhs=2  rhs=0
            ntrm = ntrm3(ib)
            ilif = ilif3(ib)
            do 60 loop = 1 , ntrm
              mexlhs(kount+loop) = lexlhs(ilif+loop)
              mexrhs(kount+loop) = lexrhs(ilif+loop)
              nrange(kount+loop) = mrange(ilif+loop)
              xp(kount+loop) = x(ilif+loop)
 60         continue
            kount = kount + ntrm
          end if
          mtrm2(ib) = kount - jlif2(ib)
 70     continue
      end if
      kont = kount
      return
      end
**==point1.f
      subroutine point1(mexlhs,mexrhs,nrange,xp,mtrm2,jlif2,
     *j,k,kont)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension mexlhs(*),mexrhs(*),nrange(*),xp(*),mtrm2(*),jlif2(*)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/lsort/lexlhs(mxcan1),lexrhs(mxcan1),mrange(mxcan1),
     *x(mxcan1),
     *ntrm2(16),ilif2(16)
INCLUDE(common/icon)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),konte(2),mmin,mmax
      kount = kont
      jp1 = j + 1
      km1 = k - 1
      if (jp1.le.km1) then
c...
c...
c... recursion through type 1 points
c...
c...
        do 70 kk = jp1 , km1
          if (npopl(kk).eq.1) then
            minm1 = mmin - 1
            maxp1 = mmax + 1
            mint = mmin
            maxt = mmax
            call minmax(kk)
            ibase = i16(kk)
            kount = 0
            do 40 ib = mmin , mmax , 2
              ilif2(ib) = kount
              if (ib.lt.maxt .and. ib.ge.minm1) then
c... coupling lhs=1  rhs=1
                ntrm = mtrm2(ib+1)
                ilif = jlif2(ib+1)
                do 20 loop = 1 , ntrm
                  lexlhs(kount+loop) = mexlhs(ilif+loop) + n2l(ib+ibase)
                  lexrhs(kount+loop) = mexrhs(ilif+loop) + n2r(ib+ibase)
                  mrange(kount+loop) = nrange(ilif+loop)
                  x(kount+loop) = xp(ilif+loop)
 20             continue
                kount = kount + ntrm
              end if
              if (ib.gt.mint .and. ib.le.maxp1) then
c... coupling  lhs=2  rhs=2
                ntrm = mtrm2(ib-1)
                ilif = jlif2(ib-1)
                do 30 loop = 1 , ntrm
                  lexlhs(kount+loop) = mexlhs(ilif+loop)
                  lexrhs(kount+loop) = mexrhs(ilif+loop)
                  mrange(kount+loop) = nrange(ilif+loop)
                  x(kount+loop) = xp(ilif+loop)
 30             continue
                kount = kount + ntrm
              end if
              ntrm2(ib) = kount - ilif2(ib)
 40         continue
c... shuffle generated data to correct position
            do 50 ib = mmin , mmax , 2
              mtrm2(ib) = ntrm2(ib)
              jlif2(ib) = ilif2(ib)
 50         continue
            do 60 loop = 1 , kount
              mexlhs(loop) = lexlhs(loop)
              mexrhs(loop) = lexrhs(loop)
              nrange(loop) = mrange(loop)
              xp(loop) = x(loop)
 60         continue
          end if
 70     continue
        kont = kount
      end if
      return
      end
**==min23.f
      subroutine min23(k)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),mmin,mmax,min2,max2,min3,max3
      maxl = maxlhs(k)
      maxr = maxrhs(k)
      minl = minlhs(k)
      minr = minrhs(k)
      max2 = min(maxr,maxl-1)
      max3 = min(maxr,maxl+1)
      min2 = max(minr,minl-1)
      min3 = max(minr,minl+1)
      return
      end
**==minmax.f
      subroutine minmax(k)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kkonte(2),mmin,mmax
      mmin = max(minlhs(k),minrhs(k))
      mmax = min(maxlhs(k),maxrhs(k))
      return
      end
**==concal.f
      subroutine concal
      implicit REAL  (a-h,o-z),integer  (i-n)
       common/lsort/vnum(18),snum(18)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),konte(2),minmax(12)
INCLUDE(common/auxh)
INCLUDE(common/icon)
      common/fcon/t,tneg,v(16),w(16),vneg(16),
     *wneg(16),a(16),b(16),aneg(16),bneg(16),c(16),e(16),f(16),
     *g2(16),factj2(16),factk2(16),factl2(16),factm2(16),
     *factn(16),factng(16),p2(16),q2(16),
     *v2(16),w2(16),w2neg(16),v2neg(16)
      common/spicon/gspin(171),hspin(171),id2(33)
      iposs = 1
      if ((max(maxsin(1),maxsin(2)+1,maxsin(3)+2)).gt.16)
     +     call caserr('CSF with more than 16 singly occupied mos')
      t = dsqrt(0.5d0)
      tneg = -t
      do 30 ib = 1 , 17
        id2(ib+ib-1) = ib
        bilbo = ib
        vnum(ib) = bilbo
        bilbo = 1.0d0/bilbo
        ic = ib
        id = ib + 1
        ixxx = 0
        do 20 mb = 1 , id
          iposs = iposs + 1
          frodo = ic
          gspin(iposs) = dsqrt(frodo*bilbo)
          frodo = ixxx
          hspin(iposs) = dsqrt(frodo*bilbo)
          ixxx = ixxx + 1
          ic = ic - 1
 20     continue
 30   continue
      do 40 ib = 2 , 18
      snum(ib) = dsqrt(vnum(ib-1))
 40   continue
      snum(1) = 0.0d0
      hot = 1.0d0/t
      do 50 ib = 1 , 16
        b(ib) = 1.0d0/vnum(ib)
        bneg(ib) = -b(ib)
        a(ib) = snum(ib)*snum(ib+2)*b(ib)
        aneg(ib) = -a(ib)
        c(ib) = aneg(ib)*aneg(ib)
        e(ib) = aneg(ib)*b(ib)*hot
        top = 1.0d0/snum(ib+1)
        v(ib) = snum(ib+2)*top
        vneg(ib) = -v(ib)
        v2(ib) = v(ib)*t
        v2neg(ib) = -v2(ib)
        wneg(ib) = snum(ib)*top
        w(ib) = -wneg(ib)
        w2(ib) = w(ib)*t
        w2neg(ib) = -w2(ib)
        factk2(ib) = (1.0d0+b(ib))*tneg
        factj2(ib) = (1.0d0-b(ib))*t
        factl2(ib) = vneg(ib)*factj2(ib)
        factm2(ib) = wneg(ib)*factk2(ib)
 50   continue
      do 60 ib = 1 , 15
        f(ib) = aneg(ib)*aneg(ib+1)
        factn(ib) = v(ib)*b(ib+1)*hot
        factng(ib) = -factn(ib)
        p2(ib) = v(ib)*v(ib+1)*tneg
        q2(ib) = w(ib)*w(ib+1)*t
 60   continue
      do 70 ib = 1 , 14
        g2(ib) = (b(ib)+b(ib+2))*t
 70   continue
      do 80 loop = 1 , 32
        numb(loop) = loop - 1
 80   continue
      m = 0
      do 90 loop = 1 , 64
        i16(loop) = m
        m = m + 16
 90   continue
      norb1 = norb + 1
      norb16 = i16(norb)
      return
      end
**==expan1.f
      subroutine expan1(conf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer conf,con1,con2
_ELSEIF(convex,i8drct)
      integer *8 conf,con1,con2
_ELSE
      REAL conf,con1,con2
_ENDIF
      dimension conf(*)
      common/expanc/con1(5),con2(5),npop1(64),min1(64),max1(64),
     *n1(1024),n2(1024)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,konte(2),minmax(12)
_IF(i8drct)
      do lp = 1 , 3
        con1(lp) = conf(lp)
      enddo
_ELSEIF(linux)
      call dcopy(3,conf,1,con1,1)
_ELSE
      do lp = 1 , 3
        con1(lp) = conf(lp)
      enddo
_ENDIF
      call upack(con1,norb,npop1)
      call genlim(con1,npop1,min1,max1,n1,n2)
      return
      end
**==expan2.f
      subroutine expan2(conf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer conf,con1,con2
_ELSEIF(convex,i8drct)
      integer *8 conf,con1,con2
_ELSE
      REAL conf,con1,con2
_ENDIF
      dimension conf(*)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),
     *npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,i,j,l,k,konte(2),minmax(12),ndimax
INCLUDE(common/icon)
      common/expanc/con1(5),con2(5),npop1(64),min1(64),max1(64),
     *n1(1024),n2(1024)
      ndiff = ndif(con1,conf,i)
      if (ndiff.le.ndimax) then
_IF(i8drct)
      do ii = 1 , 3
        con2(ii) = conf(ii)
      enddo
_ELSEIF(linux)
      call dcopy(3,conf,1,con2,1)
_ELSE
        do ii = 1 , 3
          con2(ii) = conf(ii)
        enddo
_ENDIF
        if (ndiff.ne.0) then
          call upack(con2,norb,npopr)
          mode = npopr(i) - npop1(i)
          if (mode.lt.0) then
            do 40 ii = 1 , norb
              npopl(ii) = npopr(ii)
              npopr(ii) = npop1(ii)
              m = min1(ii)
              n = max1(ii)
              minrhs(ii) = m
              maxrhs(ii) = n
              ibase = i16(ii)
              do 30 ib = m , n , 2
                n1r(ib+ibase) = n1(ib+ibase)
                n2r(ib+ibase) = n2(ib+ibase)
 30           continue
 40         continue
            call genlim(con2,npopl,minlhs,maxlhs,n1l,n2l)
            return
          else
            do 60 ii = 1 , norb
              npopl(ii) = npop1(ii)
              m = min1(ii)
              n = max1(ii)
              minlhs(ii) = m
              maxlhs(ii) = n
              ibase = i16(ii)
              do 50 ib = m , n , 2
                n1l(ib+ibase) = n1(ib+ibase)
                n2l(ib+ibase) = n2(ib+ibase)
 50           continue
 60         continue
            call genlim(con2,npopr,minrhs,maxrhs,n1r,n2r)
          end if
        else
          do 80 ii = 1 , norb
            npopl(ii) = npop1(ii)
            npopr(ii) = npop1(ii)
            m = min1(ii)
            n = max1(ii)
            minlhs(ii) = m
            minrhs(ii) = m
            maxlhs(ii) = n
            maxrhs(ii) = n
            ibase = i16(ii)
            do 70 ib = m , n , 2
              n1l(ib+ibase) = n1(ib+ibase)
              n1r(ib+ibase) = n1(ib+ibase)
              n2l(ib+ibase) = n2(ib+ibase)
              n2r(ib+ibase) = n2(ib+ibase)
 70         continue
 80       continue
          return
        end if
      end if
      return
      end
_EXTRACT(genlim,mips4)
**==genlim.f
      subroutine genlim(conf,npop,min,max,n1,n2)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSE
      REAL conf
_ENDIF
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt
_ENDIF
      dimension conf(*),npop(*),min(*),max(*),n1(*),n2(*)
INCLUDE(common/icon)
INCLUDE(common/cdcryz)
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),konte(2),minmax(12)
      common/junk/mlink(16)
c... compute no. of singly occupied orbitals
_IF(cray,t3e,convex,i8,i8drct)
      ndoc = popcnt(iand(conf(1),n32m)) + popcnt(iand(conf(2),n32m))
_ELSEIF(bits8)
      call dand8(conf(1),n32m,tmp8)
      call dand8(conf(2),n32m,tmp88)
      ndoc = popcnt(tmp8) + popcnt(tmp88)
_ELSE
      ndoc = popcnt(dand(conf(1),n32m)) + popcnt(dand(conf(2),n32m))
_ENDIF
      nsing = nel - ndoc - ndoc
c... set up data to begin recursion
      mn = mspin
      mx = mspin
      mlink(mn) = 1
c... do the recursions
      do 60 ii = 1 , norb
        ibase = i16(ii)
        if (npop(ii).ne.1) then
c... level empty or doubly occupied
          do 20 ib = mn , mx , 2
            n2(ib+ibase) = 0
            n1(ib+ibase) = mlink(ib)
 20       continue
        else
c... level singly occupied
          mnt = mn
          mxt = mx
          mn = mnt - 1
          mn1 = mnt + 1
          if (mn.lt.1) mn = mn1
          mx = mxt + 1
          mx1 = mxt - 1
          if (mx.gt.nsing) mx = mx1
          nsing = nsing - 1
          n1(mx+ibase) = 0
          n2(mn+ibase) = 0
          if (mx1.ge.mn) then
            do 30 ib = mn , mx1 , 2
              n1(ib+ibase) = mlink(ib+1)
 30         continue
          end if
          if (mx.ge.mn1) then
            do 40 ib = mn1 , mx , 2
              n2(ib+ibase) = mlink(ib-1)
 40         continue
          end if
          do 50 ib = mn , mx , 2
            mlink(ib) = n1(ib+ibase) + n2(ib+ibase)
 50       continue
        end if
        min(ii) = mn
        max(ii) = mx
 60   continue
      return
      end
_ENDEXTRACT
**==symbtr.f
      subroutine symbtr(cj,ck,ifl,nshif1,nshif2)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension cj(*),ck(*)
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
      common/lsort/nfrnt2,nstat2,nfrnt1,nstat1
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,
     *ndiff,ijkl(4),kount,kounte,minmax(12)
_IF(cray,t3e,i8)
      integer con1,con2
_ELSEIF(convex,i8drct)
      integer *8 con1,con2
_ELSE
      REAL con1,con2
_ENDIF
      common/expanc/con1(5),con2(5),
     *npop1(64),min1(64),max1(64),n1(1024),n2(1024)
      if (nshif2.eq.0) then
c... singlet/doublet/valence in 2nd configuration
        call upack2(con2(3),nfrnt2,nstat2)
        nfrnt2 = 0
      else
c... triplet case in 2nd configuration
        call upack2(con2(3),nstat2,nfrnt2)
      end if
      if (nshif1.eq.0) then
c... singlet/doublet/valence in 1st configuration
        call upack2(con1(3),nfrnt1,nstat1)
        nfrnt1 = 0
      else
c... triplet case in 1st configuration
        call upack2(con1(3),nstat1,nfrnt1)
      end if
c... translate exchange terms
      call transl(ck,mexrhs,mexlhs,nrange,xp,kounte,numbe)
c... translate direct terms
      call transl(cj,mexrhs(kounte+1),mexlhs(kounte+1),nrange(kounte+1),
     +            xp(kounte+1),kount,numb)
      ifl = numb + numbe
      return
      end
**==transl.f
      subroutine transl(c,mex2,mex1,nrange,xp,kount,numb)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension c(*)
      dimension mex2(*),mex1(*),nrange(*),xp(*)
c...
c... to translate coupling coefficients from
c... packed lexical form to rectangular array form
c...
      common/lsort/nfrnt2,nstat2,nfrnt1,nstat1
      numb = 0
      ndim = nstat2*nstat1
      if (ndim.gt.0) then
        call vclr(c,1,ndim)
        if (kount.gt.0) then
          jump = nstat2 + 1
          do 30 kk = 1 , kount
            top = xp(kk)
            if (dabs(top).ge.(1.0d-11)) then
              nr = nrange(kk)
              mexf = mex2(kk) - nfrnt2
              mexs = mex1(kk)
              if (mexf.lt.0) then
                nr = nr + mexf
                mexs = mexs - mexf
                mexf = 0
              end if
              mexs = mexs - nfrnt1
              if (mexs.lt.0) then
                nr = nr + mexs
                mexf = mexf - mexs
                mexs = 0
              end if
              mm = nstat2 - nr - mexf
              if (mm.lt.0) then
                nr = nr + mm
              end if
              mm = nstat1 - nr - mexs
              if (mm.lt.0) then
                nr = nr + mm
              end if
              if (nr.gt.0) then
                numb = numb + nr
                mm = mexs*nstat2 + mexf
c... states derived from conf 2 define most
c... rapidly varying index
                do 20 loop = 1 , nr
                  c(mm+1) = top
                  mm = mm + jump
 20             continue
_IF1()          call vfill(top,c(mm+1),jump,nr)
              end if
            end if
 30       continue
        end if
      end if
      return
      end
**==symbt.f
      subroutine symbt(c,kounte,kount)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension c(*)
      common/scra  /mexl(mxcan2),mexr(mxcan2),nr(mxcan2),x(mxcan2)
      common/dopey/ndim,ncol,nrow,isingl(62)
c... to translate all coupling coeffs to rectangular array
      call vclr(c,1,ndim)
      if (kount.ne.0) then
        nrow1 = nrow + 1
        nst = kounte + 1
        nlst = kounte + kount
        do 30 i = nst , nlst
          mr = mexr(i)
          top = x(i)
          if (mr.lt.nrow .and. dabs(top).ge.(1.0d-11)) then
            ibase = mexl(i)*nrow + mr
            mr = nr(i)
            do 20 loop = 1 , mr
              c(ibase+1) = top
              ibase = ibase + nrow1
 20         continue
_IF1()      call vfill(top,c(ibase+1),nrow1,mr)
          end if
 30     continue
      end if
      return
      end
**==ndif.f
      function ndif(con1,con2,iorr)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer c,mask,con1,con2,z3
_ELSEIF(convex)
      integer *8 shiftr
      integer *8 z3,c,mask,con1,con2
_ELSEIF(i8drct)
      integer *8 shiftr,shift
      integer *8 z3,c,mask,con1,con2
      external shift,shiftr
_ELSE
      REAL z3,c,mask,con1,con2
_ENDIF
      integer popcnt
_IF(osf,sgi,hpux11,linux)
      external popcnt, leadz
_ENDIF
c     determine number ('ndif') and positions ('ior') of orbital
c     differences between configurations 'con1' and 'con2'
c     the indices are ordened to give charge-cloud notation, where
c     the exhange part can be obtained by interchanging no's 2 and 3
      dimension con1(*),con2(*),iorr(*)
INCLUDE(common/cdcryi)
      common/symcon/norb,norb1,norb16,nword,jspppp(22),ndimax
      common/maskc/mask(64)
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z3/3/
_ELSEIF(bits8)
      call pad8(3,z3)
_ELSE
      z3 = pad(3)
_ENDIF
      ndif = 0
      k = 0
      ityp1 = 0
      do 30 i = 1 , nword
_IF(cray,t3e,convex,i8,i8drct)
        c = ieor(con1(i),con2(i))
_ELSEIF(bits8)
        call dxor8(con1(i),con2(i),c)
_ELSE
        c = dxor(con1(i),con2(i))
_ENDIF
        ndif = ndif + popcnt(c)
        if (ndif.gt.ndimax) return
 20     if (popcnt(c).ne.0) then
          k = k + 1
          ks = leadz(c)
_IF1() PAUL try ishftc version on cray, t3e
_IF(cray,t3e)
          mmmm =         iand(shift(con1(i),ks+2),3)  -
     +                   iand(shift(con2(i),ks+2),3)
_ELSEIF(convex,i8)
          mmmm =         iand(ishftc(con1(i),ks+2,64),z3)  -
     +                   iand(ishftc(con2(i),ks+2,64),z3)
_ELSEIF(i8drct)
          mmmm =         iand(shift(con1(i),ks+2),z3)  -
     +                   iand(shift(con2(i),ks+2),z3)
_ELSEIF(bits8)
          call shift8(con1(i),ks+2,tmp8)
          call dand8(tmp8,z3,tmpf)
          mmmm = ipad(tmpf)
          call shift8(con2(i),ks+2,tmp8)
          call dand8(tmp8,z3,tmpf)
          mmmm = mmmm - ipad(tmpf)
_ELSE
          mmmm =  ipad (dand(SHIFT(con1(i),ks+2),z3))  -
     +            ipad (dand(SHIFT(con2(i),ks+2),z3))
_ENDIF
          ityp = isign(1,mmmm)
          kk = k
          if (ityp.eq.ityp1.and.ndimax.le.4) then
            kk = 4
            k = k - 1
          end if
          if (kk.eq.1) ityp1 = ityp
          iorr(kk) = ks/2 + (i-1)*n32 + 1
_IF1() PAUL try first clause on cray
_IF(cray,t3e,convex,i8)
          c = ieor(ishft(mask(1),-ks),c)
_ELSEIF(i8drct)
          c = ieor(shiftr(mask(1),ks),c)
_ELSEIF(bits8)
          call shiftr8(mask(1),ks,tmp8)
          call dxor8(tmp8,c,c)
_ELSE
          c = dxor(shiftr(mask(1),ks),c)
_ENDIF
          go to 20
        end if
 30   continue
      return
      end
**==spinn.f
      subroutine spinn
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer ipack,jpack,ipa,jpa,ipaa,z1
_ELSEIF(convex,i8drct)
      integer *8 ipack,jpack,ipa,jpa,ipaa,z1
      integer *8 shiftr
      integer *8 modl,modr
      external shiftr
_ELSE
      REAL ipack,jpack,ipa,jpa,ipaa,z1
_ENDIF
      logical igtj,ieqj,npri2,nplj2
INCLUDE(common/sizes)
      common/moco/ibj,ibj1
      common/spicon/gspin(171),hspin(171),id2(33)
INCLUDE(common/helpr)
      common/junk/ipack(mxcan1),jpack(mxcan1)
      common/junk3/xfac(19),yfac(19),xxfac(19),yyfac(19)
     *,ilif(17),ntrm(17)
      common/scra  /mexlhs(mxcan2),mexrhs(mxcan2),
     *nrange(mxcan2),xp(mxcan2)
INCLUDE(common/icon)
      common/three/
     *n1l(1024),n2l(1024),minlhs(64),maxlhs(64),npopl(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),npopr(64),mode
      common/symcon/norb,norb1,norb16,nw,nwm1,nel,mspin,ndiff,
     *i,j,k,l,kount,kounte,mmin,mmax
_IF(cray,t3e,convex,i8,i8drct)
      data z1/1/
_ELSEIF(bits8)
      call pad8(1,z1)
_ELSE
      z1=pad(1)
_ENDIF
c...
c... to compute coupling coefficients for
c...   2 * (1-particle space operator) * sz
c... to allow generation of spin natural orbitals
c...
      ieqj = i.eq.j
      npri2 = npopr(i).eq.2
      ibas = i16(i)
      nplj2 = npopl(j).eq.2
      jbas = i16(j)
      nsingr = 0
      kount = 0
      do 20 kk = 1 , j
        if (npopr(kk).eq.1) nsingr = nsingr + 1
 20   continue
      nsingl = nsingr
      if (npri2) then
        nsingl = nsingl + 1
      else
        nsingl = nsingl - 1
      end if
      if (nplj2) then
        nsingl = nsingl - 1
      else
        nsingl = nsingl + 1
      end if
      call minmax(j)
      jm1 = j - 1
      im1 = i - 1
      ip1 = i + 1
      igtj = ip1.gt.jm1
      do 640 ibj = mmin , mmax , 2
        ibj1 = 1 - ibj
        ilif(ibj) = kount
        call spigen(npopl,nsingl,ipack,istate)
        call spigen(npopr,nsingr,jpack,jstate)
        do 630 ioop = 1 , istate
          ipaa = ipack(ioop)
          do 620 joop = 1 , jstate
            jpa = jpack(joop)
            ipa = ipaa
            lexl = 0
            lexr = 0
            xfac(1) = 1.0d0
            mspinl = mspin
            mspinr = mspin
            mspinx = mspin
            mrange = 1
            if (im1.ne.0) then
c...
c... recur through type 1 spin segments
c...
              do 130 kk = 1 , im1
                if (npopr(kk).ne.1) go to 130
_IF(cray,t3e,convex,i8,i8drct)
                modl = iand(ipa,z1)
                modr = iand(jpa,z1)
_ELSEIF(bits8)
                call dand8(ipa,z1,tmp8)
                modl = ipad(tmp8)
                call dand8(jpa,z1,tmp8)
                modr = ipad(tmp8)
_ELSE
                modl = ipad(dand(ipa,z1))
                modr = ipad(dand(jpa,z1))
_ENDIF

_IF1() PAUL try ishft on cray,t3e
_IF(convex,i8)
                ipa = ishft(ipa,-1)
                jpa = ishft(jpa,-1)
_ELSEIF(bits8)
                call shiftr8(ipa,1,tmp8)
                ipa = tmp8
                call shiftr8(jpa,1,tmp8)
                jpa = tmp8
_ELSE
                ipa = shiftr(ipa,1)
                jpa = shiftr(jpa,1)
_ENDIF
                ibase = i16(kk)
                yfac(1) = 0.0d0
                mbias = 0
                icase = modl + modr + modr + 1
                go to (30,50,70,90) , icase
c... case 1 = coupling 11
 30             ilbas = iky(mspinl) + id2(mspinl-mspinx+1)
                irbas = iky(mspinr) + id2(mspinr-mspinx+1)
                do loop = 1 , mrange
                  bilbo = xfac(loop)
                  yfac(loop) = hspin(ilbas)*hspin(irbas)
     +                         *bilbo + yfac(loop)
                  yfac(loop+1) = gspin(ilbas)*gspin(irbas)*bilbo
                  ilbas = ilbas + 1
                  irbas = irbas + 1
                enddo
                mbias = 1
                if (mrange.eq.mspinx) mrange = mrange - 1
                mspinx = mspinx - 1
                mspinr = mspinr - 1
                lexr = lexr + n2r(ibase+mspinr)
                mspinl = mspinl - 1
                lexl = lexl + n2l(ibase+mspinl)
                go to 110
c... case 2 = coupling 21
 50             ilbas = iky(mspinl+2) + id2(mspinl-mspinx+3)
                irbas = iky(mspinr) + id2(mspinr-mspinx+1)
                do loop = 1 , mrange
                  bilbo = xfac(loop)
                  yfac(loop) = gspin(ilbas)*hspin(irbas)
     +                         *bilbo + yfac(loop)
                  yfac(loop+1) = -hspin(ilbas)*gspin(irbas)*bilbo
                  ilbas = ilbas + 1
                  irbas = irbas + 1
                enddo
                mrange = mrange + 1
                if (mspinr.le.mspinx) then
                  mbias = 1
                  mrange = mrange - 1
                  if (mrange.eq.mspinx) mrange = mrange - 1
                end if
                mspinl = mspinl + 1
                mspinr = mspinr - 1
                lexr = lexr + n2r(ibase+mspinr)
                mspinx = min(mspinl,mspinr)
                go to 110
c... case 3 = coupling 12
 70             ilbas = iky(mspinl) + id2(mspinl-mspinx+1)
                irbas = iky(mspinr+2) + id2(mspinr-mspinx+3)
                do loop = 1 , mrange
                  bilbo = xfac(loop)
                  yfac(loop) = hspin(ilbas)*gspin(irbas)
     +                         *bilbo + yfac(loop)
                  yfac(loop+1) = -gspin(ilbas)*hspin(irbas)*bilbo
                  ilbas = ilbas + 1
                  irbas = irbas + 1
                enddo
                mrange = mrange + 1
                if (mspinl.le.mspinx) then
                  mbias = 1
                  mrange = mrange - 1
                  if (mrange.eq.mspinx) mrange = mrange - 1
                end if
                mspinr = mspinr + 1
                mspinl = mspinl - 1
                lexl = lexl + n2l(ibase+mspinl)
                mspinx = min(mspinl,mspinr)
                go to 110
c... case 4 = coupling 22
 90             ilbas = iky(mspinl+2) + id2(mspinl-mspinx+3)
                irbas = iky(mspinr+2) + id2(mspinr-mspinx+3)
                do loop = 1 , mrange
                  bilbo = xfac(loop)
                  yfac(loop) = gspin(ilbas)*gspin(irbas)
     +                         *bilbo + yfac(loop)
                  yfac(loop+1) = hspin(ilbas)*hspin(irbas)*bilbo
                  ilbas = ilbas + 1
                  irbas = irbas + 1
                enddo
                mrange = mrange + 1
                mspinx = mspinx + 1
                mspinr = mspinr + 1
                mspinl = mspinl + 1
 110            do loop = 1 , mrange
                  xfac(loop) = yfac(loop+mbias)
                enddo
 130          continue
            end if
            frodo = 0.0d0
            if (ieqj) then
c...
c...
c... occupation patterns in bra and ket equal
c... so integrate over type 2 level i=j
c...
c...
_IF(cray,t3e,convex,i8,i8drct)
              modl = iand(ipa,z1)
              modr = iand(jpa,z1)
_ELSEIF(bits8)
              call dand8(ipa,z1,tmp8)
              modl = ipad(tmp8)
              call dand8(jpa,z1,tmp8)
              modr = ipad(tmp8)
_ELSE
              modl = ipad(dand(ipa,z1))
              modr = ipad(dand(jpa,z1))
_ENDIF
              icase = modl + modr + modr + 1
              go to (520,540,570,590) , icase
            else
c...
c... top rh spin segment at level i
c...
              mrangf = mrange
              mbias = 0
              if (npri2) then
c... level i doubly occupied in the ket
_IF(cray,t3e,convex,i8,i8drct)
                modl = iand(ipa,z1)
_ELSEIF(bits8)
                call dand8(ipa,z1,tmp8)
                modl = ipad(tmp8)
_ELSE
                modl = ipad(dand(ipa,z1))
_ENDIF
_IF(convex,i8)
                ipa  =  ishft(ipa,-1)
_ELSEIF(bits8)
                call shiftr8(ipa,1,ipa)
_ELSE
                ipa = shiftr(ipa,1)
_ENDIF
                if (modl.ne.0) then
c... coupling   23
                  ilbas = iky(mspinl+2) + id2(mspinl-mspinx+3)
                  yfac(1) = 0.0d0
                  irbas = 0
                  if (mspinl.lt.mspinr) irbas = 1
                  mrangf = mrangf + irbas
                  do loop = 1 , mrange
                    bilbo = xfac(loop)
                    yfac(loop+irbas) = bilbo*hspin(ilbas)
                    xfac(loop) = bilbo*gspin(ilbas)
                    ilbas = ilbas + 1
                  enddo
                  mspinl = mspinl + 1
                else
c... coupling   13
                  ilbas = iky(mspinl) + id2(mspinl-mspinx+1)
                  mspinl = mspinl - 1
                  lexl = lexl + n2l(ibas+mspinl)
                  if (mspinl.lt.mspinx) then
                    if (mrangf.eq.mspinx) mrangf = mrangf - 1
                    mrange = mrange - 1
                    mbias = 1
                  end if
                  do loop = 1 , mrangf
                    yfac(loop) = -xfac(loop)*gspin(ilbas+loop-1)
                  enddo
                  if (mrange.ne.0) then
                    do loop = 1 , mrange
                      xfac(loop) = xfac(loop+mbias)
     +                             *hspin(ilbas+loop+mbias-1)
                    enddo
                  end if
                end if
              else
c... level i singly occupied in the ket
_IF(cray,t3e,convex,i8,i8drct)
                modr = iand(jpa,z1)
_ELSEIF(bits8)
                call dand8(jpa,z1,tmp8)
                modr = ipad(tmp8)
_ELSE
                modr = ipad(dand(jpa,z1))
_ENDIF
_IF(convex,i8)
                jpa = ishft(jpa,-1)
_ELSEIF(bits8)
                call shiftr8(jpa,1,jpa)
_ELSE
                jpa = shiftr(jpa,1)
_ENDIF
                if (modr.ne.0) then
c... coupling   02
                  irbas = iky(mspinr+2) + id2(mspinr-mspinx+3)
                  do loop = 1 , mrange
                    bilbo = xfac(loop)
                    yfac(loop) = bilbo*gspin(irbas)
                    xxfac(loop) = -bilbo*hspin(irbas)
                    irbas = irbas + 1
                  enddo
                  xfac(1) = 0.0d0
                  irbas = 0
                  if (mspinl.gt.mspinr) irbas = 1
                  do loop = 1 , mrange
                    loop1 = loop + irbas
                    xfac(loop1) = xxfac(loop)
                  enddo
                  mrange = mrange + irbas
                  mspinr = mspinr + 1
                  go to 220
                else
c... coupling   01
                  irbas = iky(mspinr) + id2(mspinr-mspinx+1)
                  mspinr = mspinr - 1
                  lexr = lexr + n2r(ibas+mspinr)
                  if (mspinr.lt.mspinx) then
                    if (mrange.eq.mspinx) mrange = mrange - 1
                    mrangf = mrangf - 1
                    mbias = 1
                    if (mrangf.eq.0) go to 200
                  end if
                  do loop = 1 , mrangf
                    yfac(loop) = xfac(loop+mbias)
     +                           *hspin(irbas+loop+mbias-1)
                  enddo
                end if
 200            do loop = 1 , mrange
                  xfac(loop) = xfac(loop)*gspin(irbas+loop-1)
                enddo
              end if
 220          if (.not.(igtj)) then
c...
c... recur through mid rh spin segments
c...
                do 430 kk = ip1 , jm1
                  if (npopl(kk).ne.1) go to 430
_IF(cray,t3e,convex,i8,i8drct)
                  modl = iand(ipa,z1)
                  modr = iand(jpa,z1)
_ELSEIF(bits8)
                  call dand8(ipa,z1,tmp8)
                  modl = ipad(tmp8)
                  call dand8(jpa,z1,tmp8)
                  modr = ipad(tmp8)
_ELSE
                  modl = ipad(dand(ipa,z1))
                  modr = ipad(dand(jpa,z1))
_ENDIF
_IF(convex,i8)
                  ipa = ishft(ipa,-1)
                  jpa = ishft(jpa,-1)
_ELSEIF(bits8)
                  call shiftr8(ipa,1,ipa)
                  call shiftr8(jpa,1,jpa)
_ELSE
                  ipa = shiftr(ipa,1)
                  jpa = shiftr(jpa,1)
_ENDIF
                  ibase = i16(kk)
                  mbias = 0
                  nbias = 0
                  xxfac(1) = 0.0d0
                  yyfac(1) = 0.0d0
                  mrp1 = mspinr + 1
                  mrm1 = mspinr - 1
                  mbiga = min(mspinl,mrp1)
                  mbigb = min(mspinl,mrm1)
                  icase = modl + modr + modr + 1
                  go to (230,270,310,350) , icase
c... case 1 = coupling 11
 230              mspinx = min(mspinl,mspinr)
                  if (mrange.ne.0) then
                    ilbas = iky(mspinl) + id2(mspinl-mbiga+1)
                    irbas = iky(mspinr) + id2(mspinr-mbiga+2)
                    do loop = 1 , mrange
                      bilbo = xfac(loop)
                      xxfac(loop) = xxfac(loop) - hspin(ilbas)
     +                              *hspin(irbas)*bilbo
                      xxfac(loop+1) = -gspin(ilbas)*gspin(irbas)*bilbo
                      ilbas = ilbas + 1
                      irbas = irbas + 1
                    enddo
                    mbias = 1
                    if (mrange.eq.mspinx) mrange = mrange - 1
                    if (mrangf.eq.0) go to 260
                  end if
                  ilbas = iky(mspinl) + id2(mspinl-mbigb+1)
                  irbas = iky(mspinr) + id2(mspinr-mbigb)
                  do loop = 1 , mrangf
                    bilbo = yfac(loop)
                    yyfac(loop) = yyfac(loop) - hspin(ilbas)
     +                            *hspin(irbas)*bilbo
                    yyfac(loop+1) = -gspin(ilbas)*gspin(irbas)*bilbo
                    ilbas = ilbas + 1
                    irbas = irbas + 1
                  enddo
                  nbias = 1
                  if (mrangf.eq.mspinx) mrangf = mrangf - 1
 260              mspinl = mspinl - 1
                  lexl = lexl + n2l(ibase+mspinl)
                  go to 300
c... case 2 = coupling 21
 270              if (mrange.ne.0) then
                    ilbas = iky(mspinl+2) + id2(mspinl-mbiga+3)
                    irbas = iky(mspinr) + id2(mspinr-mbiga+2)
                    do loop = 1 , mrange
                      bilbo = xfac(loop)
                      xxfac(loop) = xxfac(loop) - gspin(ilbas)
     +                              *hspin(irbas)*bilbo
                      xxfac(loop+1) = hspin(ilbas)*gspin(irbas)*bilbo
                      ilbas = ilbas + 1
                      irbas = irbas + 1
                    enddo
                    mrange = mrange + 1
                    if (mbiga.ge.mrp1) then
                      mbias = 1
                      mrange = mrange - 1
                      if (mrange.eq.mspinr) mrange = mrange - 1
                    end if
                    if (mrangf.eq.0) then
                      mspinl = mspinl + 1
                      go to 300
                    end if
                  end if
                  ilbas = iky(mspinl+2) + id2(mspinl-mbigb+3)
                  irbas = iky(mspinr) + id2(mspinr-mbigb)
                  do loop = 1 , mrangf
                    bilbo = yfac(loop)
                    yyfac(loop) = yyfac(loop) - gspin(ilbas)
     +                            *hspin(irbas)*bilbo
                    yyfac(loop+1) = hspin(ilbas)*gspin(irbas)*bilbo
                    ilbas = ilbas + 1
                    irbas = irbas + 1
                  enddo
                  mrangf = mrangf + 1
                  if (mbigb.ge.mrm1) then
                    nbias = 1
                    mrangf = mrangf - 1
                    if (mrangf.eq.mspinr) mrangf = mrangf - 1
                  end if
                  mspinl = mspinl + 1
 300              mspinr = mrm1
                  lexr = lexr + n2r(ibase+mspinr)
                  go to 400
c... case 3 = coupling 12
 310              if (mrange.ne.0) then
                    ilbas = iky(mspinl) + id2(mspinl-mbiga+1)
                    irbas = iky(mspinr+2) + id2(mspinr-mbiga+4)
                    do loop = 1 , mrange
                      xxfac(loop) = xxfac(loop) - hspin(ilbas)
     +                              *gspin(irbas)*xfac(loop)
                      xxfac(loop+1) = gspin(ilbas)*hspin(irbas)
     +                                *xfac(loop)
                      ilbas = ilbas + 1
                      irbas = irbas + 1
                    enddo
                    mrange = mrange + 1
                    if (mbiga.ge.mspinl) then
                      mbias = 1
                      mrange = mrange - 1
                      if (mrange.eq.mspinl) mrange = mrange - 1
                    end if
                    if (mrangf.eq.0) go to 340
                  end if
                  ilbas = iky(mspinl) + id2(mspinl-mbigb+1)
                  irbas = iky(mspinr+2) + id2(mspinr-mbigb+2)
                  do loop = 1 , mrange
                    yyfac(loop) = yyfac(loop) - hspin(ilbas)
     +                            *gspin(irbas)*yfac(loop)
                    yyfac(loop+1) = gspin(ilbas)*hspin(irbas)*yfac(loop)
                    ilbas = ilbas + 1
                    irbas = irbas + 1
                  enddo
                  mrangf = mrangf + 1
                  if (mbigb.ge.mspinl) then
                    nbias = 1
                    mrangf = mrangf - 1
                    if (mrangf.eq.mspinl) mrangf = mrangf - 1
                  end if
 340              mspinl = mspinl - 1
                  lexl = lexl + n2l(ibase+mspinl)
                  go to 390
c... case 4 = coupling 22
 350              if (mrange.ne.0) then
                    ilbas = iky(mspinl+2) + id2(mspinl-mbiga+3)
                    irbas = iky(mspinr+2) + id2(mspinr-mbiga+4)
                    do loop = 1 , mrange
                      xxfac(loop) = xxfac(loop) - gspin(ilbas)
     +                              *gspin(irbas)*xfac(loop)
                      xxfac(loop+1) = -hspin(ilbas)*hspin(irbas)
     +                                *xfac(loop)
                      ilbas = ilbas + 1
                      irbas = irbas + 1
                    enddo
                    mrange = mrange + 1
                    if (mrangf.eq.0) go to 380
                  end if
                  ilbas = iky(mspinl+2) + id2(mspinl-mbigb+3)
                  irbas = iky(mspinr+2) + id2(mspinr-mbigb+2)
                  do loop = 1 , mrangf
                    yyfac(loop) = yyfac(loop) - gspin(ilbas)
     +                            *gspin(irbas)*yfac(loop)
                    yyfac(loop+1) = -hspin(ilbas)*hspin(irbas)
     +                              *yfac(loop)
                    ilbas = ilbas + 1
                    irbas = irbas + 1
                  enddo
                  mrangf = mrangf + 1
 380              mspinl = mspinl + 1
 390              mspinr = mrp1
c... move data to correct locations
 400              if (mrange.ne.0) then
                    do loop = 1 , mrange
                      xfac(loop) = xxfac(loop+mbias)
                    enddo
                    if (mrangf.eq.0) go to 430
                  end if
                  do loop = 1 , mrangf
                    yfac(loop) = yyfac(loop+nbias)
                  enddo
 430            continue
              end if
c...
c... bottom segment of rh spin range at level j
c...
              mbiga = min(mspinl,mspinr+1)
              mbigb = min(mspinl,mspinr-1)
              if (nplj2) then
c... level j doubly occupied in the bra
_IF(cray,t3e,convex,i8,i8drct)
                jpa = iand(jpa,z1)
                if (jpa.ne.0) then
_ELSEIF(bits8)
                call dand8(jpa,z1,jpa)
                if (ipad(jpa).ne.0) then
_ELSE
                jpa = dand(jpa,z1)
                if (ipad(jpa).ne.0) then
_ENDIF
c... coupling 32
                  ilbas = iky(mspinr+2) + id2(mspinr-mbiga+4)
                  do loop = 1 , mrange
                    frodo = frodo - gspin(ilbas+loop-1)*xfac(loop)
                  enddo
                  if (mrangf.ne.0) then
                    ilbas = iky(mspinr+2) + id2(mspinr-mbigb+2)
                    do loop = 1 , mrangf
                      frodo = frodo + hspin(ilbas+loop-1)*yfac(loop)
                    enddo
                  end if
                else
c... coupling 31
                  if (mrange.ne.0) then
                    ilbas = iky(mspinr) + id2(mspinr-mbiga+2)
                    do loop = 1 , mrange
                      frodo = frodo - hspin(ilbas+loop-1)*xfac(loop)
                    enddo
                    if (mrangf.eq.0) then
                      lexr = lexr + n2r(jbas+mspinr-1)
                      go to 610
                    end if
                  end if
                  ilbas = iky(mspinr) + id2(mspinr-mbigb)
                  do loop = 1 , mrangf
                    frodo = frodo - gspin(ilbas+loop-1)*yfac(loop)
                  enddo
                  lexr = lexr + n2r(jbas+mspinr-1)
                end if
              else
c... level j singly occupied in the bra
_IF(cray,t3e,convex,i8,i8drct)
                ipa = iand(ipa,z1)
                if (ipa.ne.0) then
_ELSEIF(bits8)
                call dand8(ipa,z1,ipa)
                if (ipad(ipa).ne.0) then
_ELSE
                ipa = dand(ipa,z1)
                if (ipad(ipa).ne.0) then
_ENDIF
c... coupling 20
                  if (mrange.ne.0) then
                    ilbas = iky(mspinl+2) + id2(mspinl-mbiga+3)
                    do loop = 1 , mrange
                      frodo = frodo - hspin(ilbas+loop-1)*xfac(loop)
                    enddo
                    if (mrangf.eq.0) go to 610
                  end if
                  ilbas = iky(mspinl+2) + id2(mspinl-mbigb+3)
                  do loop = 1 , mrangf
                    frodo = frodo - gspin(ilbas+loop-1)*yfac(loop)
                  enddo
                else
c... coupling 10
                  if (mrange.ne.0) then
                    ilbas = iky(mspinl) + id2(mspinl-mbiga+1)
                    do loop = 1 , mrange
                      frodo = gspin(ilbas+loop-1)*xfac(loop) + frodo
                    enddo
                    if (mrangf.eq.0) then
                      lexl = lexl + n2l(jbas+mspinl-1)
                      go to 610
                    end if
                  end if
                  ilbas = iky(mspinl) + id2(mspinl-mbigb+1)
                  do loop = 1 , mrangf
                    frodo = frodo - hspin(ilbas+loop-1)*yfac(loop)
                  enddo
                  lexl = lexl + n2l(jbas+mspinl-1)
                end if
              end if
              go to 610
            end if
c... case 1 = coupling 11
 520        ilbas = iky(mspinl) + id2(mspinl-mspinx+1)
            do loop = 1 , mrange
              frodo = (gspin(ilbas)**2-hspin(ilbas)**2)*xfac(loop)
     +                + frodo
              ilbas = ilbas + 1
            enddo
            lexl = lexl + n2l(jbas+mspinl-1)
            go to 560
c... case 2 = coupling 21
 540        ilbas = iky(mspinl+2) + id2(mspinl-mspinx+3)
            irbas = iky(mspinr) + id2(mspinr-mspinx+1)
            do loop = 1 , mrange
              frodo = frodo - (hspin(ilbas)*gspin(irbas)+gspin(ilbas)
     +                *hspin(irbas))*xfac(loop)
              ilbas = ilbas + 1
              irbas = irbas + 1
            enddo
 560        lexr = lexr + n2r(jbas+mspinr-1)
            go to 610
c... case 3 = coupling 12
 570        ilbas = iky(mspinl) + id2(mspinl-mspinx+1)
            irbas = iky(mspinr+2) + id2(mspinr-mspinx+3)
            do loop = 1 , mrange
              frodo = frodo - (hspin(ilbas)*gspin(irbas)+gspin(ilbas)
     +                *hspin(irbas))*xfac(loop)
              ilbas = ilbas + 1
              irbas = irbas + 1
            enddo
            lexl = lexl + n2l(jbas+mspinl-1)
            go to 610
c... case 4 = coupling 22
 590        ilbas = iky(mspinl+2) + id2(mspinl-mspinx+3)
            do loop = 1 , mrange
              frodo = (hspin(ilbas)**2-gspin(ilbas)**2)*xfac(loop)
     +                + frodo
              ilbas = ilbas + 1
            enddo
 610        if (dabs(frodo).ge.(1.0d-11)) then
              kount = kount + 1
              xp(kount) = frodo
              mexlhs(kount) = lexl
              mexrhs(kount) = lexr
            end if
 620      continue
 630    continue
        ntrm(ibj) = kount - ilif(ibj)
 640  continue
      if (kount.ne.0) then
c...
c... complete recursion over spin free segments to bottom of graph
c...
_IF(cray,t3e)
       nav = lenwrd()
       kount2 = (kount-1)/nav + 1
_ENDIF
        call setsto(kount,1,nrange)
        call point1(mexlhs,mexrhs,nrange,xp,ntrm,ilif,j,norb1,kount)
        if (mode.lt.0) then
_IF(cray,t3e)
         call fmove(mexlhs,ipack,kount2)
         call fmove(mexrhs,mexlhs,kount2)
         call fmove(ipack,mexrhs,kount2)
_ELSE
         call iswap(kount,mexlhs,1,mexrhs,1)
_ENDIF
        end if
      end if
      return
      end
**==spigen.f
      subroutine spigen(npopl,nsingl,ipack,nstate)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer ipack,jpack,ipa,jpa,z1,z0
_ELSEIF(convex)
      integer *8 ipack,jpack,ipa,jpa,z1,z0
_ELSEIF(i8drct)
      integer *8 ipack,jpack,ipa,jpa,z1,z0
      integer *8 shift
      external shift
_ELSE
      REAL ipack,jpack,ipa,jpa,z1,z0
_ENDIF
INCLUDE(common/sizes)
      dimension npopl(*),ipack(*)
      common/moco/ibj,ibj1
      common/symcon/norb(6),mspin,ndiff,i,j
      common/lsort/ispin(mxcan1),jspin(mxcan1),jpack(mxcan1)
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z1/0,1/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(1,z1)
_ELSE
      z0 = pad(0)
      z1 = pad(1)
_ENDIF
c...
c... to generate all spin possibilities connecting top of graph
c... to spin ibj at level j in packed binary format where
c... 0=spin coupling 1  .. 1=spin coupling 2
c... empty and doubly occupied levels ignored in this process
c...
      nsing = nsingl
      nstate = 1
      jpa = z1
      ipack(1) = z0
      ispin(1) = mspin - ibj
      do 40 kk = 1 , j
        if (npopl(kk).eq.1) then
c... level singly occupied so interesting
          nsing = nsing - 1
          mstate = 0
          do 20 istate = 1 , nstate
            isp = ispin(istate)
            ipa = ipack(istate)
c... try mode 2 coupling
            jsp = isp + 1
            if (nsing.ge.iabs(jsp)) then
              mstate = mstate + 1
              jspin(mstate) = jsp
_IF(cray,t3e,convex,i8,i8drct)
              jpack(mstate) = ior(ipa,jpa)
_ELSEIF(bits8)
              call dor8(ipa,jpa,jpack(mstate))
_ELSE
              jpack(mstate) = dor(ipa,jpa)
_ENDIF
            end if
c... try mode 1 coupling
            jsp = isp - 1
            if (jsp.ge.ibj1 .and. nsing.ge.iabs(jsp)) then
              mstate = mstate + 1
              jspin(mstate) = jsp
              jpack(mstate) = ipa
            end if
 20       continue
_IF(convex,i8)
          jpa = ishftc(jpa,1,64)
_ELSEIF(bits8)
          call shift8(jpa,1,jpa)
_ELSEIF(sgi,ipsc)
          jpa = shiftc(jpa,1)
_ELSE
          jpa = shift(jpa,1)
_ENDIF
c... move data to correct locations
          nstate = mstate
          do istate = 1 , mstate
            ipack(istate) = jpack(istate)
            ispin(istate) = jspin(istate)
          enddo
        end if
 40   continue
      return
      end
**==hbill.f
      subroutine hbill(cr,icr)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer icr
_ELSEIF(convex,i8drct)
      integer *8 icr
_ELSE
      REAL icr
_ENDIF
INCLUDE(common/sizes)
c... model=1 end of coupling coefficients
c... model=2 ijkl vacuum/vacuum
c... model=3 ijkl doublet/doublet
c... model=4 ijkl n-2/n-2
c... model=5 ijka+fock(ia) doublet/vacuum
c... model=6 ijka+fock(ia) n-2/doublet
c... model=7 ijab+iaib+fock(ab) doublet/doublet
c... model=8 ijab+iaib+fock(ab) n-2/n-2
c... model=9 iajb doublet/doublet
c... model=10 iajb n-2/n-2
c... model=11 iajb p+q supermatrix  n-2/vacuum
c... model=12 iaib n-2/vacuum
c... model=13 iabc p+q supermatrix  n-2/doublet
      dimension cr(*),icr(*)
INCLUDE(common/corctl)
INCLUDE(common/ccntl)
      common/symcon/norbbb(7),ndiff,iiii(18),ndimax
INCLUDE(common/auxh)
INCLUDE(common/table)
      common/stoctl/khbill,khbm1e,khbm2e
      common/spew/rmode(6),word(4)
INCLUDE(common/iofile)
INCLUDE(common/disktl)
INCLUDE(common/disc)
INCLUDE(common/prints)
INCLUDE(common/ccepa)
      common /bufb / ifi,ifo,iiiblk,jjjblk
      equivalence (rmode(1),model), (rmode(2),ic), (rmode(3),ir),
     *            (rmode(4),inti) , (rmode(5),intj), (rmode(6),intk)
      equivalence (word(1),m1)
c... spin coupling coefficient evaluation control routine
      m1 = 1
      top = cpulft(1)
      if(.not.oprint(28)) write (iwr,6010) top
      call brtsb(nxblk)
      index(1) = nxblk
c...
c... (ij/kl) and (ij)
c...
c.... vacuum-vacuum
      kcj = khbm2e + 1
      nbuf = maxpat*maxpat
      kck = kcj + nbuf
      kbuf = kck + nbuf
      lbuf = kbuf + nbuf
      nuse = lbuf + nbuf - 1
      call corfai('H-build')
_IFN(ga)
ccepa     compute vacuum-vacuum h-matrix-elements for mrdcepa
      if (icepa.lt.10) go to 300
         indcep = nxblk
         model = 22
         call cepijkl(cr(kcj),cr(kck),cr(kbuf),cr(lbuf),icr,cr)
         model = 1
         call wrtsb(rmode(1),1)
         call ertsb
c...  reset wrtsb system
         nxblk = ipos(numci)
         call brtsb(nxblk)
         index(1) = nxblk
300   continue
ccepa
_ENDIF
      model = 2
      call ijkl(cr,icr,cr(kcj),cr(kck),cr(kbuf),cr(lbuf),1,nnst)
      call whersb(1)
      nn1 = nnst
      model = 3
      do 20 isym = 1 , nirr
        ns = nn1 + 1
        nn1 = nn1 + ncond(isym)
c.... doublet-doublet
        call ijkl(cr,icr,cr(kcj),cr(kck),cr(kbuf),cr(lbuf),ns,nn1)
 20   continue
      call whersb(2)
      model = 4
      do 30 isym = 1 , nirr
        ns = nn1 + 1
        nn1 = nn1 + nconst(isym)
c... singlet-singlet and triplet-triplet
        call ijkl(cr,icr,cr(kcj),cr(kck),cr(kbuf),cr(lbuf),ns,nn1)
 30   continue
      call whersb(3)
c...
c... (ij/ka) and (ia)
c...
      call basmak
      kcj = khbill + 1
      kck = kcj + nbuf
      ifo = kck + nbuf - 1
      ifi = ifo
      iiiblk = 1
c... doublet - vacuum
      if (nmin1.ne.0) then
        call ijka(icr,cr(kcj),nstrt1,nlast1,1,nnst)
c... singlet/triplet - doublet
        call ijka(icr,cr(kcj),nstrt2,nlast2,nstrt1,nlast1)
        call whersb(5)
c...
c... (ij/ab) (ia/jb) and (ab)
c...
c... doublet - doublet
        model = 7
        call ijab(cr,icr,cr(kcj),cr(kck),nstrt1,nlast1,2)
        call whersb(6)
      end if
c.. n-2  -  n-2
      if (nmin2.ne.0) then
        model = 8
        call ijab(cr,icr,cr(kcj),cr(kck),nstrt2,nlast2,1)
        call whersb(7)
      end if
      call ertsc(cr)
      call whersb(9)
      if (nmin2.ne.0) then
c... singlet/triplet - vacuum
        ifi = ifo
        iiiblk = 1
        call iajb(cr,icr,cr(kcj))
        call ertsc(cr)
c...
c... (ia/bc)
c...
c... singlet/triplet - doublet
        call iabc(icr,cr(kcj))
      end if
      call wrtsb(word,1)
      call ertsb
c...
c... coupling coefficients for natural orbital generator
c...
c... model=1 end of coefficients
c...
c... spinfree natural orbitals
c... model=2 vacuum/vacuum  (ij)
c... model=3 vacuum/vacuum diagonal cases  (ii)
c... model=4 doublet/doublet (ij)
c... model=5 doublet/doublet diagonal cases (ii)
c... model=6 n-2/n-2 (ij)
c... model=7 n-2/n-2 diagonal cases (ii)
c... model=8 doublet/vacuum (ia)
c... model=9 n-2/doublet (ia)
c...
c... spin natural orbitals
c... model=10 vacuum/vacuum (ii) + (ij)
c... model=12 doublet/doublet (ii) + (ij)
c... model=14 n-2/n-2 (ii) + (ij)
c... model=15 doublet/vacuum (ia)
c... model=16 n-2/doublet (ia)
c... doublet/doublet and n-2/n-2 (ab) not labelled by model
      ndimax = 2
      index(197) = iposun(numci)
      call brtsb(index(197))
c... vacuum/vacuum (ij)
      call natij(icr,cr(kcj),1,nnst,2)
c... doublet/doublet (ij)
      nn1 = nnst
      do 40 isym = 1 , nirr
        ns = nn1 + 1
        nn1 = nn1 + ncond(isym)
        call natij(icr,cr(kcj),ns,nn1,4)
 40   continue
c... n-2/n-2 (ij)
      do 50 isym = 1 , nirr
        ns = nn1 + 1
        nn1 = nn1 + nconst(isym)
        call natij(icr,cr(kcj),ns,nn1,6)
 50   continue
c... doublet/vacuum (ia)
      call natia(icr,cr(kcj),nstrt1,nlast1,1,nnst,8)
c... n-2/doublet (ia)
      call natia(icr,cr(kcj),nstrt2,nlast2,nstrt1,nlast1,9)
      call wrtsb(word,4)
c... doublet/doublet and n-2/n-2 (ab) spin natorb
      call natdub(icr,cr(kcj))
      call ertsb
      return
 6010 format (/1x,80('=')//' **** hamiltonian builder called at',f9.2,
     +        ' secs')
      end
**==ijkl.f
      subroutine ijkl(cr,iconf,cj,ck,cx,cy,nstart,nlast)
_IF(parallel)
c MPP
c MPP   parallel (partial) matrix-element generator over outer (ic) loop
c MPP   (that's good for genz as well and probably not too course)
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c     external ihint
      dimension cr(*),cj(*),ck(*),cx(*),cy(*)
      dimension xnumb(2)
c     compute the h-matrix between vacuum states
c     **     this subroutine computes the internal parts for     **
c     ** the (n-1/h/n-1) and (n-2/h/n-2) matrix elements as well **
      common/erg/potnuc,totnuc,energy
      common/spew/rmode(3)
      common/symcon/norbbb(7),nd,ii,jj,kk,ll,konte(14)
      common/three/n1l(1024),n2l(1024),minlhs(64),maxlhs(64),iocc(64),
     *n1r(1024),n2r(1024),minrhs(64),maxrhs(64),jocc(64),modelr
INCLUDE(common/sizes)
INCLUDE(common/orb)
INCLUDE(common/helpr)
      dimension iconf(5,*)
_IF(parallel)
c MPP
      logical oipsci
c MPP
_ENDIF
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      data xnumb/1.0d0,2.0d0/
c
      ihint(i,j) = iky(norbas(max(i,j))) + min(i,j)
     +             + nbase(iro(max(i,j)))
c
_IF(parallel)
c MPP
      iflop = iipsci()
c MPP
_ENDIF
      if (nlast.ge.nstart) then
        do 60 ic = nstart , nlast
_IF(parallel)
c MPP
         if (oipsci()) go to 60
c MPP
_ENDIF
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncols)
          do 50 ir = nstart , ic
            call expan2(iconf(1,ir))
c.... establish number of differences
c.... if more then 4 => no future in it
            if (nd.le.4) then
              call upack2(iconf(3,ir),nrowt,nrows)
              ndims = nrows*ncols
              ndimt = nrowt*ncolt
              call vclr(cx,1,ndims)
              call vclr(cy,1,ndimt)
              if (nd.lt.2) then
c.... no differences => diagonal element
c....  ( h is the common  coulomb-exchange/1-el contribution )
                h = potnuc
                do 30 i = 1 , nint
                  iocci = iocc(i)
                  if (iocci.ne.0) then
                    iocci2 = iocci - 2
                    pocci = xnumb(iocci)
                    je = i - 1
                    if (je.ne.0) then
                      do 20 j = 1 , je
                        ioccj = iocc(j)
                        if (ioccj.ne.0) then
                          coul = cr(irint(i,i,j,j))
                          exch = cr(irint(i,j,i,j))
                          ioccij = iocci2 + ioccj
                          if (ioccij.le.0) then
c.... coefficient  (both i and j are singly occupied)
                            call symb2(j,i)
                            h = h + coul
                            call symbtr(cj,ck,ifl,0,0)
                            if (ifl.ne.0) then
                              call daxpy(ndims,exch,ck,1,cx,1)
                            end if
                            call symbtr(cj,ck,ifl,1,1)
                            if (ifl.ne.0) then
                              call daxpy(ndimt,exch,ck,1,cy,1)
                            end if
                          else
c.... either one or both doubly occupied => (2j-k)*occ contribution
                            h = (coul+coul-exch)*xnumb(ioccij) + h
                          end if
                        end if
 20                   continue
                    end if
c.... (i/i) and (ii/ii) contributions
                    h = h + cr(ihint(i,i))*pocci
                    if (iocci2.eq.0) then
                      h = h + cr(irint(i,i,i,i))
                    end if
                  end if
 30             continue
c.... add coulomb/1-el contribution to diagonals
                if (ncols.ne.0) then
                  call ihk(cx,h,ncols)
_IF1()            call vsadd(cx,ncols+1,h,cx,ncols+1,ncols)
                end if
                if (ncolt.ne.0) then
                  call ihk(cy,h,ncolt)
_IF1()            call vsadd(cy,ncolt+1,h,cy,ncolt+1,ncolt)
                end if
              else if (nd.eq.2) then
c---- end of diagonal contributions --------------------------------
c.... single difference
                h = cr(ihint(ii,jj))
c.... first do the orbitals with differences
                if (jocc(ii).eq.2) h = h + cr(irint(ii,ii,ii,jj))
                if (iocc(jj).eq.2) h = h + cr(irint(ii,jj,jj,jj))
c.... fock-type contributions
                do 40 i = 1 , nint
                  iocci = iocc(i)
                  if (iocci.ne.0 .and. i.ne.ii .and. i.ne.jj) then
                    coul = cr(irint(ii,jj,i,i))
                    exch = cr(irint(ii,i,jj,i))
                    if (iocci.le.1) then
c.... i singly occupied => interesting coupling
                      call symb2(i,0)
                      call symbtr(cj,ck,ifl,0,0)
                      h = h + coul
                      if (ifl.ne.0) then
                        call daxpy(ndims,exch,ck,1,cx,1)
                      end if
                      call symbtr(cj,ck,ifl,1,1)
                      if (ifl.ne.0) then
                        call daxpy(ndimt,exch,ck,1,cy,1)
                      end if
                    else
c.... i doubly occupied => add contribution to h
                      h = h + coul + coul - exch
                    end if
                  end if
 40             continue
c.... now add accumulated h + 2j - k contributions
                call symb1
                call sym1tr(ck,ifl,0,0)
                if (ifl.ne.0) then
                  call daxpy(ndims,h,ck,1,cx,1)
                end if
                call sym1tr(ck,ifl,1,1)
                if (ifl.ne.0) then
                  call daxpy(ndimt,h,ck,1,cy,1)
                end if
              else
c---- end of delta = 1 -------------------------------------------------
c.... 2 orthogonalities
                call symb2(0,0)
                coul = cr(irint(ii,jj,kk,ll))
                exch = cr(irint(ii,kk,jj,ll))
                call symbtr(cj,ck,ifl,0,0)
                if (ifl.ne.0) then
                  call ihj(cx,cj,ck,coul,exch,ndims)
                end if
                call symbtr(cj,ck,ifl,1,1)
                if (ifl.ne.0) then
                  call ihj(cy,cj,ck,coul,exch,ndimt)
                end if
              end if
              call wrtsb(rmode(1),3)
              call wrtsb(cx,ndims)
              call wrtsb(cy,ndimt)
            end if
c---- end of delta = 2 -------------------------------------------------
 50       continue
 60     continue
      end if
      return
      end
**==ihj.f
      subroutine ihj(c,cx,cy,x,y,n)
      implicit REAL  (a-h,o-z),integer  (i-n)
c
c
      dimension c(*),cx(*),cy(*)
      call daxpy(n,x,cx,1,c,1)
      call daxpy(n,y,cy,1,c,1)
      return
      end
_IFN(fps,convex,alliant)
**==ihk.f
      subroutine ihk(c,x,n)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension c(*)
      n1 = n + 1
      i = 1
      do 20 l = 1 , n
        c(i) = c(i) + x
        i = i + n1
 20   continue
      return
      end
_ENDIF
**==ijka.f
      subroutine ijka(iconf,cj,mstart,mlast,nstart,nlast)
_IF(parallel)
c MPP
c MPP   parallelism should be ok with cijka
c MPP   i.e. for most external state (n-1/n-2) => outer loop
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,conner
_ELSEIF(convex,i8drct)
      integer *8 iconf,conner
_ELSE
      REAL iconf,conner
_ENDIF
_IF(parallel)
c MPP
      logical oipsci
c MPP
_ENDIF
      logical ncoltt,kltj
      dimension cj(*)
      common/dopey/ndimsd,ncols,ncold,isingl(62)
      common/spew/rmode(4)
      common/symcon/norbbb(7),ndiff,ii,jj,kk,ll,kount,kounte
INCLUDE(common/sizes)
INCLUDE(common/orb)
      common/expanc/conner(10),iocc(64)
INCLUDE(common/symci)
INCLUDE(common/table)
INCLUDE(common/helpr)
      dimension iconf(5,*)
      equivalence (rmode(1),model),(rmode(2),ic),
c... to evaluate spin coupling coeffs. for (ij/ka) and (ia) integrals
     *            (rmode(3),ir   ),(rmode(4),ninti)
      model = model + 1
      if (mlast.ge.mstart .and. nlast.ge.nstart) then
c... loop over internally least occupied states
_IF(parallel)
c MPP
      iflop = iipsci()
c MPP
_ENDIF
        do 50 ic = mstart , mlast
_IF(parallel)
c MPP
         if (oipsci()) go to 50
c MPP
_ENDIF
_IF(cray,t3e,convex,i8,i8drct)
          isym = iconf(4,ic)
_ELSE
          isym = ipad(iconf(4,ic))
_ENDIF
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncolss)
          ncols = ncolss + ncolt
          nintt = 0
          do 20 i = 1 , nint
            if (iocc(i).eq.1) then
              nintt = nintt + 1
              isingl(nintt) = i
            end if
 20       continue
c... above loop isolates singly occ. internal levels of 1st config.
c... loop over internally more occupied state
          do 40 ir = nstart , nlast
_IF(cray,t3e,convex,i8,i8drct)
            jsym = iconf(4,ir)
_ELSE
            jsym = ipad(iconf(4,ir))
_ENDIF
            ksym = mult(jsym,isym)
            nek = norbe(ksym)
            if (nek.ne.0) then
              call expan2(iconf(1,ir))
_IF(cray,t3e,convex,i8,i8drct)
              ncold = iconf(3,ir)
_ELSE
              ncold = ipad(iconf(3,ir))
_ENDIF
              ndimsd = ncold*ncols
              mubas = ncolss*ncold
              ncoltt = ncolt.eq.0
              if (ksym.eq.jsym .and. nek.eq.1) then
                if (ncolss.eq.0) go to 40
                ndimsd = mubas
                ncoltt = .true.
              end if
              kltj = ksym.lt.jsym .or. ncoltt
              nubas = ndimsd - mubas
c... if no. of differences gt 4 no future in it
              if (ndiff.lt.4) then
c... single orthogonality-----------------------------------------------
c... first evaluate fock contributions
                ninti = iap(ii)
                call wrtsb(rmode(1),4)
                call symb1
                call symbt(cj,0,kount)
                if (.not.(kltj)) then
                  call dscal(nubas,-1.0d0
     +                       ,cj(mubas+1),1)
                end if
                call wrtsb(cj,ndimsd)
c... now for contribution of singly occupied internal levels
                if (nintt.ne.0) then
                  do 30 i = 1 , nintt
                    j = isingl(i)
                    if (j.ne.ii) then
                      ninti = ijkap(j,ii,j)
                      call wrtsb(rmode(1),4)
                      call symb2(j,0)
                      call symbt(cj,0,kounte)
                      if (.not.(kltj)) then
                call dscal(nubas,-1.0d0
     +            ,cj(mubas+1),1)
                      end if
                      call wrtsb(cj,ndimsd)
                    end if
 30               continue
                end if
              else if (ndiff.eq.4) then
c... double orthogonality-----------------------------------------------
                call symb2(0,0)
                if (kk.le.nint) then
                  ll = kk
                  j = jj
                  jj = ii
                  ii = j
                end if
                ninti = ijkap(jj,ii,ll)
                call wrtsb(rmode(1),4)
                call symbt(cj,kounte,kount)
                if (.not.(kltj)) then
              call dscal(nubas,-1.0d0
     +        ,cj(mubas+1),1)
                end if
c... check for coincidence of ii and ll meaning no distinct exchange
                if (ii.ne.ll) then
                  call wrtsb(cj,ndimsd)
                  ninti = ijkap(jj,ll,ii)
                  call wrtsb(rmode(1),4)
                  call symbt(cj,0,kounte)
                  if (.not.(kltj)) then
                    call dscal(nubas,-1.0d0
     +                       ,cj(mubas+1),1)
                  end if
                end if
                call wrtsb(cj,ndimsd)
              end if
            end if
 40       continue
 50     continue
      end if
      return
      end
**==ijab.f
      subroutine ijab(cr,iconf,cj,ck,nstart,nlast,i1)
_IF(parallel)
c MPP
c MPP   parallel for fock-contributions in line with cijab
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
_IF(parallel)
c MPP
      logical oipsci
c MPP
_ENDIF
      dimension cj(*),ck(*)
      dimension iconf(5,*),cr(*)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
      common/spew/rmode(5)
      common/symcon/norbbb(7),ndiff,ii,jj,kk,ll,kount,kounte
     *,minmax(12),ndimax
INCLUDE(common/orb)
      common/three/n1l(1024),n2l(1024),minl(64),maxl(64),iocc(64),
     *n1r(1024),n2r(1024),minr(64),maxr(64),jocc(64),mode
c... to calc. spin coupling coeffs. for (ij/ab),(ia/jb),(ab) ints
c... for n-1/n-1 and n-2/n-2 interactions
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *            (rmode(4),kkkbas),(rmode(5),itrans)
      ndimax = 2
      nintp1 = nint + i1
      modelx = model
      modely = modelx + 2
_IF(parallel)
c MPP
      iflop = iipsci()
c MPP
_ENDIF
      do 40 ic = nstart , nlast
_IF(parallel)
c MPP
       if (oipsci()) go to 40
c MPP
_ENDIF
        call expan1(iconf(1,ic))
        call upack2(iconf(3,ic),ncolt,ncols)
        do 30 ir = nstart , ic
          call expan2(iconf(1,ir))
c... if more than 2 differences no future in it
          if (ndiff.le.2) then
            call upack2(iconf(3,ir),nrowt,nrows)
            ndimss = ncols*nrows
            ndimtt = ncolt*nrowt
            ndimst = ncols*nrowt
            if (ndiff.ne.0) then
c... single orthogonality-----------------------------------------------
              call symb2(nintp1,0)
c... must find out which state is more occ at level ii
              if (mode.lt.0) then
c... conf 2 has 1 more electron in level jj
                ninti = jj
                nintj = ii
              else
c... conf 2 has 1 more electron in level ii
                ninti = ii
                nintj = jj
              end if
              if (ndimss+ndimtt.ne.0) then
                kkkbas = ijabp(ninti,nintj)
                call wrtsb(rmode(1),4)
              end if
              model = modely
              kkkbas = iajbp(ninti,nintj)
              call wrtsc(cr,rmode,5)
              model = modelx
              if (ndimss.ne.0) then
                call symbtr(cj,ck,ifl,0,0)
                call wrtsb(cj,ndimss)
                call wrtsc(cr,ck,ndimss)
              end if
c... triplet - triplet
              if (ndimtt.ne.0) then
                call symbtr(cj,ck,ifl,1,1)
                call wrtsb(cj,ndimtt)
                call wrtsc(cr,ck,ndimtt)
              end if
c... singlet - triplet
              if (ndimst.ne.0) then
                call symbtr(cj,ck,ifl,0,1)
                call wrtsc(cr,ck,ndimst)
              end if
c... triplet - singlet
              ndimts = ncolt*nrows
              if (ndimts.ne.0) then
                call symbtr(cj,ck,ifl,1,0)
                call wrtsc(cr,ck,ndimts)
              end if
            else
c... no orthogonality---------------------------------------------------
              kkkbas = 0
c... output fock reference
              call wrtsb(rmode(1),4)
c... evaluate exchange interactions with singly occ internal levels
              do 20 ninti = 1 , nint
                if (iocc(ninti).eq.1) then
                  kkkbas = mbfock(ninti)
                  call wrtsb(rmode(1),4)
                  call symb2(ninti,nintp1)
c... singlet - singlet or doublet - doublet
                  if (ndimss.ne.0) then
                    call symbtr(cj,ck,ifl,0,0)
                    call wrtsb(ck,ndimss)
                  end if
c... triplet - triplet
                  if (ndimtt.ne.0) then
                    call symbtr(cj,ck,ifl,1,1)
                    call wrtsb(ck,ndimtt)
c... singlet - triplet
                    if (ndimst.ne.0) then
                      call symbtr(cj,ck,ifl,0,1)
                      call wrtsb(ck,ndimst)
                    end if
                  end if
                end if
c... no need to load triplet / singlet in diagonal case
c... as this block is equal to singlet / triplet transposed
 20           continue
            end if
          end if
 30     continue
 40   continue
c...     reset nintp1
      nintp1 = nint + 1
c...
      return
      end
**==iajb.f
      subroutine iajb(cr,iconf,cj)
_IF(parallel)
c MPP
c MPP   parallel on outer index
c MPP   no synchronisation need
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
_IF(parallel)
c MPP
       logical oipsci
c MPP
_ENDIF
      dimension cj(*)
      dimension cr(*),iconf(5,*)
      common/dopey/ndimsd,ncols,ncold,isingl(62)
      common/spew/rmode(6)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
      common/symcon/norb,norbbb(6),ndiff,ii,jj,kk,ll,kount,kounte
     *,minmax(12),ndimax
INCLUDE(common/ccntl)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),ninti),(rmode(5),nintj),(rmode(6),nintk)
c... singlet/triplet - vacuum (ia/jb) contributions
      if (nnst.ne.0) then
        ndimax = 4
c... loop over n-2 states
        do 30 ic = nstrt2 , nlast2
_IF(parallel)
c MPP
         if (oipsci()) go to 30
c MPP
_ENDIF
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncolss)
          ncols = ncolss + ncolt
c... loop over vacuum states
          do 20 ir = 1 , nnst
            call expan2(iconf(1,ir))
c... check if no. of diffs gt 4
            if (ndiff.le.4) then
_IF(cray,t3e,convex,i8,i8drct)
              ncold = iconf(3,ir)
_ELSE
              ncold = ipad(iconf(3,ir))
_ENDIF
              if (ii.ne.ll) then
                ninti = iajbp(ll,ii)
                model = 11
                ndimsd = ncold*ncols
                call wrtsb(rmode(1),4)
                call symb2(0,0)
                call symbt(cj,kounte,kount)
                call wrtsb(cj,ndimsd)
              else
c... if i=j only singlet/vacuum interactions
                ndimsd = ncolss*ncold
                if (ndimsd.ne.0) then
                  ninti = mbfock(ii)
                  model = 12
                  call wrtsc(cr,rmode,4)
                  call symb2(0,0)
                  call symbt(cj,kounte,kount)
                  call wrtsc(cr,cj,ndimsd)
                end if
              end if
            end if
 20       continue
 30     continue
      end if
      return
      end
**==iabc.f
      subroutine iabc(iconf,cj)
_IF(parallel)
c MPP
c MPP   parallel on outer index
c MPP   no synchronisation needed
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
_IF(parallel)
c MPP
      logical oipsci
c MPP
_ENDIF
      dimension cj(*)
      dimension iconf(5,*)
      common/dopey/ndimsd,ncols,ncold,mnumb(62),iibas(62)
      common/stoctl/khbill
INCLUDE(common/corctl)
      common/spew/rmode(6)
      common/symcon/norb,norbbb(6),ndiff,ii,jj,kk,ll,kount,kounte
     *,minmax(12),ndimax
INCLUDE(common/orb)
INCLUDE(common/ccntl)
_IF(ga)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
_IF(ga)
      logical opcivect
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
_ENDIF
_ENDIF
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),ninti),(rmode(5),nintj),(rmode(6),nintk)
      narsil = 0
_IF(ga)
      call setmynorbe
_ENDIF
c... singlet/triplet - doublet (ia/bc) contributions
      if (nmin1.ne.0) then
        ndimax = 2
        model = 13
        ngott = (ngot-khbill)/nint
        do 20 loop = 1 , nint
          mnumb(loop) = 0
          iibas(loop) = narsil
          narsil = narsil + ngott
 20     continue
_IF(ga)
        iofff=0
_ENDIF
c... loop over n-2 states
        do 40 ic = nstrt2 , nlast2
_IF(parallel)
c MPP
c       if (oipsci()) go to 40
c MPP
_ENDIF
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncolss)
          ncols = ncolss + ncolt
c... loop over doublet states
          do 30 ir = nstrt1 , nlast1
            call expan2(iconf(1,ir))
c... check if no. of diffs gt 2
            if (ndiff.le.2) then
_IF(cray,t3e,convex,i8,i8drct)
              ncold = iconf(3,ir)
_IF(ga)
              isr = iconf(4,ir)
_ENDIF
_ELSE
              ncold = ipad(iconf(3,ir))
_IF(ga)
              isr = ipad(iconf(4,ir))
_ENDIF
_ENDIF
_IF(ga)
              if (.not.opcivect) then
                if (mynorbe(ii,isr).eq.0) goto 30
              endif
_ENDIF
              ndimsd = ncold*ncols
              ndim4 = ndimsd + 4
_IF(ga)
              mars = iofff
_ELSE
              mars = mnumb(ii)
_ENDIF
              narsil = mars + ndim4
c             if (narsil.gt.ngott) call errors(991)
              if (narsil.gt.ngott) then
         print *,'remco narsil::',narsil,ngott
c    1 call caserr('memory problem in iabc => give it (a lot) more')
       call caserr('memory problem in iabc => give it (a lot) more')
              endif
_IF(ga)
              jjbas = mars
_ELSE
              jjbas = iibas(ii) + mars
_ENDIF
              ninti = ii
              call fmove(rmode,cj(jjbas+1),4)
              call symb1
              call symbt(cj(jjbas+5),0,kount)
_IF(ga)
              iofff = narsil
_ELSE
              mnumb(ii) = narsil
_ENDIF
            end if
 30       continue
 40     continue
_IF(ga)
        call wrtsb(cj(1),iofff)
_ELSE
        do 50 loop = 1 , nint
          call wrtsb(cj(iibas(loop)+1),mnumb(loop))
 50     continue
_ENDIF
      end if
      return
      end
**==natij.f
      subroutine natij(iconf,cj,nstart,nlast,nodel)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,con1,con2
_ELSEIF(convex,i8drct)
      integer *8 iconf,con1,con2
_ELSE
      REAL iconf,con1,con2
_ENDIF
      logical ls2,lt2
      dimension cj(*)
      dimension xnumb(2)
c... to compute coupling coefficients for (ij) integrals
c... for natural orbital generation
      common/spew/rmode(4)
      common/symcon/norbbb(7),nd,ii,jj
INCLUDE(common/sizes)
INCLUDE(common/helpr)
      common/expanc/con1(5),con2(5),iocc(64)
INCLUDE(common/orb)
INCLUDE(common/natorb)
      dimension iconf(5,*)
      equivalence (bilbo,ir)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *            (rmode(4),ninti)
      data xnumb/-2.0d0,-1.0d0/
      if (nlast.ge.nstart) then
        nodel1 = nodel + 1
        nodel8 = nodel + 8
        do 40 ic = nstart , nlast
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncols)
          ic1 = ic - 1
          if (ic1.ge.nstart) then
            do 20 ir = nstart , ic1
              call expan2(iconf(1,ir))
c.... establish number of differences
c.... if more then 2 => no future in it
              if (nd.eq.2) then
c... single difference
                ninti = iky(jj) + ii
                model = nodel
                call wrtsb(rmode(1),4)
                call upack2(iconf(3,ir),nrowt,nrows)
                nss = nrows*ncols
                ntt = nrowt*ncolt
                call symb1
                call sym1tr(cj,ifl,0,0)
                call wrtsb(cj,nss)
                call sym1tr(cj,ifl,1,1)
                call wrtsb(cj,ntt)
                if (.not.(mspin1)) then
                  model = nodel8
                  call wrtsb(rmode(1),4)
                  call spinn
                  call sym1tr(cj,ifl,0,0)
                  call wrtsb(cj,nss)
                  call sym1tr(cj,ifl,1,1)
                  call wrtsb(cj,ntt)
                end if
              end if
c---- end of delta = 1 -------------------------------------------------
 20         continue
          end if
c... no differences
          nss = ncols*ncols
          ntt = ncolt*ncolt
          ls2 = ncols.ge.2
          lt2 = ncolt.ge.2
          call expan2(iconf(1,ic))
          do 30 joke = 1 , nint
            ic1 = iocc(joke)
            if (ic1.ne.2) then
              model = nodel1
              bilbo = xnumb(ic1+1)
              ninti = iky(joke+1)
              call wrtsb(rmode(1),4)
              if (.not.(mspin1)) then
                if (ic1.ne.0) then
                  model = nodel8
                  ii = joke
                  jj = joke
                  ir = ic
                  call wrtsb(rmode(1),4)
                  call spinn
                  call sym1tr(cj,ifl,0,0)
                  if (ls2) call dubzer(cj,ncols)
                  call wrtsb(cj,nss)
                  call sym1tr(cj,ifl,1,1)
                  if (lt2) call dubzer(cj,ncolt)
                  call wrtsb(cj,ntt)
                end if
              end if
            end if
 30       continue
c---- end of delta = 0 -------------------------------------------------
 40     continue
      end if
      return
      end
**==dubzer.f
      subroutine dubzer(cj,n)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension cj(*)
      joop = n
      do 20 loop = 2 , n
        loop1 = loop - 1
        call vclr(cj(joop+1),1,loop1)
        call dscal(loop1,2.0d0,cj(loop),n)
        joop = joop + n
 20   continue
      return
      end
**==natia.f
      subroutine natia(iconf,cj,mstart,mlast,nstart,nlast,nodel)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      logical klej
      dimension cj(*)
      common/dopey/ndimsd,ncolu,ncold,isingl(62)
      common/spew/rmode(4)
      common/symcon/norbbb(7),ndiff,ii,jj,kk,ll,kount,kounte
INCLUDE(common/symci)
INCLUDE(common/table)
INCLUDE(common/natorb)
      dimension iconf(5,*)
      equivalence (rmode(1),model),(rmode(2),ic),
     *            (rmode(3),ir   ),(rmode(4),ninti)
c... to evaluate coupling coeffs. for (ia) integrals
c... for natural orbital generation
      if (mlast.ge.mstart .and. nlast.ge.nstart) then
        nodel7 = nodel + 7
c... loop over internally least occupied states
        do 30 ic = mstart , mlast
_IF(cray,t3e,convex,i8,i8drct)
          isym = iconf(4,ic)
_ELSE
          isym = ipad(iconf(4,ic))
_ENDIF
          call expan1(iconf(1,ic))
          call upack2(iconf(3,ic),ncolt,ncolss)
          ncols = ncolss + ncolt
c... loop over internally more occupied state
          do 20 ir = nstart , nlast
            call expan2(iconf(1,ir))
            if (ndiff.eq.2) then
_IF(cray,t3e,convex,i8,i8drct)
              jsym = iconf(4,ir)
_ELSE
              jsym = ipad(iconf(4,ir))
_ENDIF
              ksym = mult(jsym,isym)
              ne = norbe(ksym)
c... check that there is only one orthogonality
              if (ndiff.eq.2 .and. ne.ne.0) then
                ncolx = ncols
                if ((isym+ne).eq.2) ncolx = ncolss
                if (ncolx.ne.0) then
_IF(cray,t3e,convex,i8,i8drct)
                  ncold = iconf(3,ir)
_ELSE
                  ncold = ipad(iconf(3,ir))
_ENDIF
                  ndimsd = ncold*ncolx
                  ninti = ii
                  model = nodel
                  call wrtsb(rmode(1),4)
                  call symb1
                  call symbt(cj,0,kount)
                  klej = ksym.le.jsym
                  if (.not.(klej)) then
                    mubas = ncolss*ncold
                    ndltmu = ndimsd - mubas
                    call dscal(ndltmu,-1.0d0
     +                       ,cj(mubas+1),1)
                  end if
                  call wrtsb(cj,ndimsd)
                  if (.not.(mspin1)) then
                    model = nodel7
                    call wrtsb(rmode(1),4)
                    call spinn
                    call symbt(cj,0,kount)
                    if (.not.(klej)) then
                      call dscal(ndltmu,-1.0d0
     +                       ,cj(mubas+1),1)
                    end if
                    call wrtsb(cj,ndimsd)
                  end if
                end if
              end if
            end if
 20       continue
 30     continue
      end if
      return
      end
**==natdub.f
      subroutine natdub(iconf,cj)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cj(*)
      common/dopey/ndim,ncol,nrow
INCLUDE(common/ccntl)
      common/symcon/norb,norbb(7),ii,jj,kk,ll,kount
INCLUDE(common/natorb)
      dimension iconf(5,*)
c... to compute coupling coeffs for doublet/doublet and n-2/n-2
c... (ab) interactions for spin natorb
      if (.not.(mspin1)) then
        ii = norb
        jj = norb
        if (nmin1.ne.0) then
c... doublet/doublet
          do 20 ic = nstrt1 , nlast1
            call expan1(iconf(1,ic))
            call expan2(iconf(1,ic))
            call spinn
_IF(cray,t3e,convex,i8,i8drct)
            nrow = iconf(3,ic)
_ELSE
            nrow = ipad(iconf(3,ic))
_ENDIF
            ndim = nrow*nrow
            call symbt(cj,0,kount)
            call wrtsb(cj,ndim)
 20       continue
        end if
c... n-2/n-2
        if (nmin2.ne.0) then
          do 40 ic = nstrt2 , nlast2
            call upack2(iconf(3,ic),lamct,lamc)
c... if no triplets not interested
            if (lamct.ne.0) then
              call expan1(iconf(1,ic))
              call expan2(iconf(1,ic))
              call spinn
              call sym1tr(cj,ifl,1,1)
              call dscal(lamct,0.5d0,cj(1),lamct+1)
              if (lamct.ne.1) then
                junke = lamct
                do 30 loop = 2 , lamct
                  call vclr(cj(junke+1),1,loop-1)
                  junke = junke + lamct
 30             continue
              end if
              call wrtsb(cj,lamct*lamct)
              if (lamc.ne.0) then
                call sym1tr(cj,ifl,1,0)
                call wrtsb(cj,lamct*lamc)
              end if
            end if
 40       continue
        end if
      end if
      return
      end
**==irint.f
      function irint(ii,jj,kk,ll)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
c... iro i.r. index of internal orbitals in program order
c... norbi no. of internal orbitals of given i.r. index
c... npairi no. of internal pairs for given symmetry block
c... mbase base point for this  block of fully internal 2-e ints
c...
c.... to compute index of fully internal integral
c...
      i = max(ii,jj)
      j = min(ii,jj)
      k = max(kk,ll)
      l = min(kk,ll)
      if (i.lt.k) then
      else if (i.eq.k) then
        if (j.ge.l) go to 20
      else
        go to 20
      end if
      m = i
      i = k
      k = m
      m = j
      j = l
      l = m
 20   isi = iro(i)
      jsi = iro(j)
      ijsi = iky(isi) + jsi
      ksi = iro(k)
      lsi = iro(l)
      i = norbas(i)
      j = norbas(j)
      k = norbas(k)
      l = norbas(l)
      if (isi.ne.jsi) then
        ij = j + (i-1)*norbi(jsi)
      else
        ij = iky(i) + j
      end if
      if (ksi.ne.lsi) then
        kl = l + (k-1)*norbi(lsi)
      else
        kl = iky(k) + l
      end if
      if (ijsi.ne.(iky(ksi)+lsi)) then
        nasty = ij + (kl-1)*npairi(ijsi)
      else
        nasty = (ij*(ij-1))/2 + kl
      end if
      irint = nasty + mbase(niky(isi)+isi*jsi+ksi)
c... sub block order is i=1,nirr  j=1,i  k=1,i
      return
      end
_IF1()**==ihint.f
_IF1()     function ihint(ii,jj)
_IF1()     implicit REAL  (a-h,o-z),integer  (i-n)
_IF1()     common/helpr/iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
_IF1()    *m2exch,m2fock,m2psup,m2qsup,
_IF1()    *mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8)
_IF1()    *,iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36)
_IF1()    *,norbas(nd200)
_IF1()     if(ii.ge.jj)goto 2
_IF1()     i=jj
_IF1()     j=ii
_IF1()     goto 1
_IF1()  2  i=ii
_IF1()     j=jj
_IF1()  1  ihint=iky(norbas(i))+j+nbase(iro(i))
_IF1()     return
_IF1()     end
**==ijkap.f
      function ijkap(ii,jj,k)
c... to compute base point of (ij/ka) vector of integrals
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/table)
INCLUDE(common/symci)
INCLUDE(common/helpr)
      i = max(ii,jj)
      j = min(ii,jj)
      isi = iro(i)
      jsi = iro(j)
      i = norbas(i)
      j = norbas(j)
      if (isi.ne.jsi) then
        ij = j + (i-1)*norbi(jsi)
      else
        ij = iky(i) + j
      end if
      ksi = iro(k)
      nasi = mult(isi,jsi)
      nasi = mult(nasi,ksi)
      ijkap = (ij*norbi(ksi)+mnbasi(ksi)+k)*norbe(nasi)
     +        + mbase((iky(isi)+jsi-1)*nirr+ksi)
c... sub block ordering is i=1,nirr   j=1,i  k=1,nirr
      return
      end
**==basmak.f
      subroutine basmak
c... to generate constants to help in retrieval of (ijk/a)
c... (ij/ab) and (ia/jb) integrals
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/symci)
INCLUDE(common/table)
      common/stoctl/khbill,khbm2e,khiaib,khiajb,khiabc
INCLUDE(common/auxh)
      mbass = 0
      nbass = 0
      m2e = 0
      lmn = nint
      do 50 ioop = 1 , nirr
        nbase(ioop) = lmn
        lmn = lmn + norbe(ioop)
        klnorb = norbsi(ioop)
        item2 = norbtr(ioop)
        nsitot(ioop) = klnorb*nspips(ioop)
        ntrtot(ioop) = item2*nspipt(ioop)
        if (nubsi.lt.klnorb) nubsi = klnorb
        if (nubtr.lt.item2) nubtr = item2
        if (nubsq.lt.norbsq(ioop)) nubsq = norbsq(ioop)
        inorb = 0
        do 20 joop = 1 , nirr
          item2 = norbe(joop)
          mbassq(joop,ioop) = ibassq(joop,ioop) - item2
          ibass3(joop,ioop) = inorb
          inorb = inorb + item2*norbsq(mult(joop,ioop))
 20     continue
        niky(ioop) = inorb
        if (m2psup.lt.inorb) m2psup = inorb
        inorb = norbi(ioop)
        do 40 joop = 1 , ioop
          item2 = mult(joop,ioop)
          madsi(joop,ioop) = ibassi(joop,item2)
          madsq(joop,ioop) = ibassq(joop,item2)
          nadsq(joop,ioop) = ibassq(ioop,item2)
          if (item2.ne.1) then
            ijnorb = norbi(joop)*inorb
            klnorb = ijnorb
          else
            ijnorb = iky(inorb+1)
            klnorb = 0
            if (inorb.gt.1) klnorb = iky(inorb)
          end if
          ijsi = mult(ioop,joop)
          mbass = mbass + 1
          mbcoul(mbass) = m2coul
          mbexch(mbass) = m2exch
          m2coul = m2coul + klnorb*norbsi(ijsi)
          m2exch = m2exch + klnorb*norbsq(ijsi)
          do 30 koop = 1 , nirr
            nbass = nbass + 1
            mbase(nbass) = m2e
            m2e = m2e + ijnorb*norbi(koop)*norbe(mult(koop,ijsi))
 30       continue
 40     continue
 50   continue
      khbm2e = khbill + m2e
      khiaib = khbill
      do 60 ioop = 1 , nint
        iap(ioop) = iap(ioop) + khbm2e
        mbfock(ioop) = khiaib
        khiaib = khiaib + norbsi(1)
 60   continue
      m2fock = khiaib - khbill
      khiajb = khiaib + m2coul
      m2em1e = m2e + m1ia
      m2exab = m2fock + m2coul
      khiabc = khbill + m2psup
      do 70 ioop = 1 , mbass
        mbexch(ioop) = mbexch(ioop) + khbill
        mbcoul(ioop) = mbcoul(ioop) + khiaib
 70   continue
      do 80 ioop = 1 , nbass
        mbase(ioop) = mbase(ioop) + khbill
 80   continue
      return
      end
**==iajbp.f
      function iajbp(ii,jj)
      implicit REAL  (a-h,o-z),integer  (i-n)
       common/spew/rmode(6)
INCLUDE(common/sizes)
INCLUDE(common/table)
INCLUDE(common/symci)
INCLUDE(common/helpr)
      equivalence (rmode(5),itrans),(rmode(6),ijsi)
c... to compute base of (ia/jb) matrix
      itrans = ii - jj
      if (itrans.lt.0) then
        i = jj
        j = ii
      else
        i = ii
        j = jj
      end if
      isi = iro(i)
      jsi = iro(j)
      ijsi = mult(isi,jsi)
      i = norbas(i)
      if (isi.ne.jsi) then
        ix = norbi(jsi)*(i-1)
      else
        ix = iky(i-1)
      end if
      iajbp = mbexch(iky(isi)+jsi) + norbsq(ijsi)*(ix+norbas(j)-1)
      return
      end
**==ijabp.f
      function ijabp(ii,jj)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/spew/rmode(6)
INCLUDE(common/sizes)
INCLUDE(common/table)
INCLUDE(common/symci)
INCLUDE(common/helpr)
      equivalence (rmode(6),ijsi)
c... to compute base of (ij/ab) matrix
      i = max(ii,jj)
      j = min(ii,jj)
      isi = iro(i)
      jsi = iro(j)
      ijsi = mult(isi,jsi)
      i = norbas(i)
      if (isi.ne.jsi) then
        ix = norbi(jsi)*(i-1)
      else
        ix = iky(i-1)
      end if
      ijabp = mbcoul(iky(isi)+jsi) + norbsi(ijsi)*(ix+norbas(j)-1)
      return
      end
_EXTRACT(gijkl,sv1)
**==gijkl.f
      subroutine gijkl(cr,icr)
_IF(parallel)
c MPP
c MPP   inititate presort and gather ijkl-integrals (GGOP)
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF1()      external ihint
_IF(linux)
      external fget
_ENDIF
INCLUDE(common/sizes)
      dimension cr(*),icr(*)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
INCLUDE(common/machin)
INCLUDE(common/prnprn)
INCLUDE(common/ccntl)
INCLUDE(common/mapp)
      common/stoctl/khbill,khbm1e,khbm2e,ijklno
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/orb)
      common/blkin/value(510),mmword,mspace,mult8(8,8)
INCLUDE(common/iofile)
      common/bufb /nwbnwb,lnklnk,buf2(1)
_IF(ibm,vax)
      common/stak/btri,mlowww(2),iblock
_ELSE
      common/stak/btri,mlowww,nscrp,iblock
_ENDIF
INCLUDE(common/presrt)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/scra/vscr(3680),istscr(3680)
_IF(ibm,vax)
     *,iscr(3680),jscr(3680),kscr(3680),lscr(3680)
_ELSE
     *,ijscr(3680),klscr(3680)
      common/mapper/iky2(5,maxorb),i4096(1)
_ENDIF
      common/junk/icyb(3412),jcyb(3412),kcyb(3412),lcyb(3412)
INCLUDE(common/file12)
      common/lsort/occ(maxorb),potn,core,
     +    ilsort(4*maxorb+7),itype(nd200)
INCLUDE(common/discc)
INCLUDE(common/corctl)
      common/erg/potnuc,totnuc,energy
INCLUDE(common/prints)
INCLUDE(common/symci)
_IF(ga)
      dimension nactive(2*maxorb+6)
      data isecb/470/
_ENDIF
c
      ihint(i,j) = iky(norbas(max(i,j))) + min(i,j)
     +             + nbase(iro(max(i,j)))
c
_IF(parallel)
c MPP
c MPP   open seperate sortfiles for each node
c MPP
      call closbf(0)
      call setbfa(-1)
c MPP
_ENDIF
      joop = 0
      nin = 0
      nscrp = 0
      ivmax = 0
      nav = lenwrd()
c...
c... 1-e ij in core
c... 1-e ia bucket 1
c... 1-e ab bucket 2
c... 2-e ijkl bucket 8
c... 2-e ijka bucket 3
c... 2-e ijab bucket 4
c... 2-e iiab bucket 5
c... 2-e iajb bucket 6
c... 2-e iaib bucket 7
c... 2-e iabc bucket 8+nirr+i
c... 2-e abcd bucket 8+iksym
c...
      iblock = iblkp
      top = cpulft(1)
      khbm1e = khbill + m1e
      khbm2e = khbm1e + m2e
      nirr8 = nirr + 8
      nint9 = nint + nirr8
_IF(ibm,vax)
      nuse = khbm1e
_ELSE
      nuse = 0
      nusei = 0
_ENDIF
      do 20 loop = 1 , nint9
        lnumb(loop) = nuse
_IF(ibm,vax)
        nuse = nuse + nsz510
_ELSE
        nuse = nuse + nsz340
        ibasen(loop) = nusei
        nusei = nusei + nsz680
_ENDIF
 20   continue
      do 30 loop = 1 , ijklno
        mbase(loop) = mbase(loop) + khbm1e
 30   continue
      do 40 loop = 1 , nirr
        joop = iky(norbi(loop)+norbe(loop)+1) + joop
        nbase(loop) = nbase(loop) + khbill
 40   continue
_IF(ibm,vax)
      nuse = max(khbm2e,nint9*nsz510+khbm1e)
_ELSE
      nuse = max(khbm2e,nint9*nsz340*2+khbm1e)
_ENDIF
      call corfai('presort')
_IF(ibm,vax)
      call setstz(nint9*nsz510,0,cr(khbm1e+1))
_ELSE
      i10 = khbm1e + 1
      i20 = i10 + nint9*nsz340
      i20 = (i20-1)*nav + 1
_IF(cray,t3e,i8)
      call setsto(nint9*nsz340*2,0,cr(i10))
_ELSE
      call setstz(nint9*nsz340*2,0,cr(i10))
_ENDIF
_ENDIF
      call setsto(nd200,-1000,itype)
      call setsto(ninnex,1,itype)
      call setsto(nint,0,itype)
      call setsto(1360,0,i205)
_IF(ibm,vax)
      call setsto(7360,0,kscr)
_ELSE
      call setsto(3680,0,klscr)
_ENDIF
c...
c... sort 1-elec integrals
c...
      call secget(isec3,1004,jblkk)
      call rdedx(occ,mach(15),jblkk,idaf)

      if(.not.oprint(28)) write (iwr,6010) top , isec3 , ibl3d , 
     +     yed(idaf)
      call corfai('presort')
      potnuc = core
      totnuc = potn
      jblkk = jblkk + lensec(mach(15))
      call search(jblkk,idaf)
      valmax = 0.0d0
_IF(ga)
c
c     now restore =active= and =core= specifications
c
      call secget(isecb,1005,iblka)
      call readi(nactive,mach(14)*nav,iblka,idaf)
c
c     numberofints=(nint+next)*(nint+next+1)/2
      numberofints=nactive(1)*(nactive(1)+1)/2
      kword=min(511,numberofints)
      ibkk=jblkk
      call rdedx(value,kword,ibkk,idaf)
      numberofints=numberofints-kword
      ibkk=ibkk+1
_ELSE
      call fget(value,kword,idaf)
_ENDIF
      nwbnwb = nsz340
      do 60 iii = 1 , nd200
        iiii = mapei(iii)
        do 50 jjj = 1 , iii
          j = mapei(jjj)
          i = iiii
          if (i.lt.j) then
            i = j
            j = iiii
          end if
          nin = nin + 1
          if (nin.gt.kword) then
_IF(ga)
            kword=min(511,numberofints)
            call rdedx(value,kword,ibkk,idaf)
            numberofints=numberofints-kword
            ibkk=ibkk+1
_ELSE
            call fget(value,kword,idaf)
_ENDIF
            if(nscrp.ge.3170)
_IF(ibm,vax)
     *       call stack(cr,nscrp)
_ELSE
     *       call stack(cr(i10),icr(i20))
_ENDIF
            nin = 1
          end if
          val = value(nin)
          if (iro(j).eq.iro(i)) then
            istack = itype(i) + itype(j)
            if (istack.lt.0) go to 50
            if (istack.eq.0) then
c... (i/i)
              cr(ihint(i,j)) = val
            else
c... (e/i) or (e/e)
              nscrp = nscrp + 1
              vscr(nscrp) = val
              istscr(nscrp) = istack
_IF(ibm,vax)
              iscr(nscrp) = i
              jscr(nscrp) = j
_ELSE
              ijscr(nscrp) = j + i4096(i)
_ENDIF
            end if
            joop = joop - 1
            if (joop.eq.0) go to 70
c...   symmetry forbidden integral
          else if (dabs(val).gt.dabs(valmax)) then
            valmax = val
            ivmax = i
            jvmax = j
          end if
 50     continue
 60   continue
c...
c... sort 2-elec integrals
c...
 70   if (ivmax.ne.0) write (iwr,6020) ivmax , jvmax , valmax
c
      do 90 loop = 1 , nirr
        do 80 joop = 1 , nirr
          mult8(joop,loop) = mult(joop,loop) + 8
 80     continue
 90   continue
      do 180 ioopn = 1 , nsect2
        ibl = iblk2(ioopn)
        lbl = ibl - lblk2(ioopn)
        if (lbl.ne.0) then
          iunit = notap2(ioopn)
          call search(ibl,iunit)
 100      call fget(value,joop,iunit)
        if(nscrp.ge.3001)
_IF(ibm,vax)
     *      call stack(cr,nscrp)
_ELSE
     *      call stack(cr(i10),icr(i20))
_ENDIF
          if (joop.ne.0) then
c... process block
            if (modei) then
_IF(ibm,vax)
              call upak8v(value(num2e+1),i205)
_ELSE
              call unpack(value(num2e+1),lab816,i205,numlab)
_ENDIF
            else
              call upak8w(value(num2e+1),i205,mapei)
            end if
_IFN(ibm,vax)
            int4 = 1
_ENDIF
            do 170 iwword = 1 , mmword
_IF(ibm,vax)
              i = i205(iwword)
              j = j205(iwword)
              k = k205(iwword)
              l = l205(iwword)
_ELSEIF(littleendian)
              j = i205(int4  )
              i = i205(int4+1)
              l = i205(int4+2)
              k = i205(int4+3)
_ELSE
              i = i205(int4  )
              j = i205(int4+1)
              k = i205(int4+2)
              l = i205(int4+3)
_ENDIF
              val = value(iwword)
              iroi = iro(i)
              iroki = mult8(iro(k),iroi)
              irol = iro(l)
              if (mult8(iro(j),irol).ne.iroki) go to 160
              iti = itype(i) + itype(j) + itype(k) + itype(l)
              if (iti.lt.0) go to 160
              if (iti.eq.0) then
c... (ii/ii)
                istack = 8
                go to 150
              else
                go to (110,120,130,140) , iti
              end if
c... (ai/ii)
 110          if (nmin1.eq.0) go to 160
              istack = 3
              go to 150
c... (aa/ii) or (ai/ai)
 120          if (j.le.nint) then
c... (aj/ai)
                if (j.eq.l) then
c... (ai/ai)
                  istack = 7
                else
                  istack = 6
                end if
c... (aa/ij)
              else if (k.eq.l) then
c... (aa/ii)
                istack = 5
              else
                istack = 4
              end if
              go to 150
c... (aa/ai)-standard order   or   (ai/aa)
 130          if (nmin12.eq.0) go to 160
              if (l.gt.nint) then
                joop = i
                i = k
                k = joop
                joop = j
                j = l
                l = joop
              end if
              istack = l + nirr8
              go to 150
c... (aa/aa)
c... p(ki,lj) reference
 140          if (nmin2.eq.0) go to 160
              istack = iroki
              if (i.eq.j) then
                if (i.eq.k) val = val + val
              else if (k.ne.l) then
                bilbo = val
                if (i.eq.k .or. j.eq.l) val = bilbo + bilbo
                nscrp = nscrp + 1
                vscr(nscrp) = val
                istscr(nscrp) = istack
_IF(ibm,vax)
                iscr(nscrp) = i
                jscr(nscrp) = j
                kscr(nscrp) = k
                lscr(nscrp) = l
_ELSE
                ijscr(nscrp) = j + i4096(i)
                klscr(nscrp) = l + i4096(k)
_ENDIF
                val = bilbo
                if (j.eq.k) val = bilbo + bilbo
c... p(li,kj) reference
                istack = mult8(irol,iroi)
                joop = k
                k = l
                l = joop
              else
                if (k.eq.j) val = val + val
              end if
 150          nscrp = nscrp + 1
              vscr(nscrp) = val
              istscr(nscrp) = istack
_IF(ibm,vax)
              iscr(nscrp) = i
              jscr(nscrp) = j
              kscr(nscrp) = k
              lscr(nscrp) = l
 160          continue
_ELSE
              ijscr(nscrp) = j + i4096(i)
              klscr(nscrp) = l + i4096(k)
 160          int4 = int4 + 4
_ENDIF
 170        continue
            lbl = lbl + 1
            if (lbl.ne.0) go to 100
          end if
        end if
 180  continue
c...
      if (nscrp.ne.0)
_IF(ibm,vax)
     * call stack(cr,nscrp)
_ELSE
     * call stack(cr(i10),icr(i20))
_ENDIF
c...
c... clear up bucketing
c...
      lbl = nsz340 + 1
      do 190 joop = 1 , nint9
        nwb = nwbuf(joop)
        if (nwb.ne.0) then
          call stopbk
          nwbnwb = nwb
          lnklnk = link(joop)
          iiii = lnumb(joop)
          jjjj = ibasen(joop)
_IF(ibm,vax)
          call fmove(cr(iiii+1),buf2,nwb)
          call fmove(cr(iiii+lbl),buf2(lbl),min(nsz170,nwb))
_ELSE
          call dcopy(nwb,cr(i10+iiii),1,buf2,1)
          call pack(buf2(nsz341),lab1632,icr(i20+jjjj),nsz680)
_ENDIF
c     write(6,*) 'gijkl - buf2'
c     j = nsz340/2 + nsz341
c     write(6,99996) (buf2(i),i=nsz341,j)
c99996 format(1x,4(1x,z16))
          call sttout
          link(joop) = iblock
          iblock = iblock + nsz
        end if
 190  continue
      if(.not.oprint(28)) write (iwr,6030) iblock
      call stopbk
      iblkp = iblock
c...
c... get (ij/kl) integrals into core
c...
      nuse = khbm2e
      call corfai('presort')
      call vclr(cr(khbm1e+1),1,m2e)
      call setsto(nsz340,0,icyb)
      call setsto(nsz340,0,jcyb)
      call setsto(nsz340,0,kcyb)
      call setsto(nsz340,0,lcyb)
      iblok = link(8)
 200  if (iblok.ge.0) then
        call rdbak(iblok)
        call stopbk
        call upak8s(buf2,icyb)
        if (odebug(40)) call wrtsor('ijkl  ',icyb,jcyb,kcyb,lcyb,buf2,
     +                              nwbnwb)
        do 210 loop = 1 , nwbnwb
          cr(irint(icyb(loop),jcyb(loop),kcyb(loop),lcyb(loop)))
     +      = buf2(loop)
 210    continue
        iblok = lnklnk
        go to 200
      else
_IF(parallel)
c MPP
c MPP    all integrals read // gop (merge) the 2-electron ints
c MPP    and (therefore) synchronise
c MPP
        call pggop(7001,cr(khbm1e+1),m2e,'+',vscr,3680)
c MPP
_ENDIF
        return
      end if
 6010 format (/1x,80('-')
     +        //' **** transformed integral pre-sort called at',f9.2,
     +        ' secs'//
     +        ' transformed 1-electron integrals restored from section',
     +        i4,' of dumpfile starting at block',i8,' of ',a4/)
 6020 format (/' *** symmetry check / largest forbidden h-integral :',
     +        2i4,e12.5)
 6030 format (/' sort file',i8,' blocks long')
      end
**==stack.f
_IF(ibm,vax)
      subroutine stack(buf,nstak)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension buf(*)
INCLUDE(common/presrt)
      common/scra/v(3680),ibu(3680)
INCLUDE(common/blksiz)
     *,i(3680),j(3680),k(3680),l(3680)
      common/stak/btri,mlowww(2),iblock
      common/bufb/nwbnwb,lnklnk,buf2(1)
      mstak = 1
      call stackx(nstak)
 20   istak = ld340(mstak,nstak,nsz170,nwbuf,lnumb,v,ibu,i,buf)
      if (istak.ne.0) then
        call stopbk
        lnklnk = link(istak)
        ib = lnumb(istak) + 1
        call fmove(buf(ib),buf2,nsz510)
        call sttout
        link(istak) = iblock
        iblock = iblock + nsz
        mstak = mstak + 1
        if (mstak.le.nstak) go to 20
      end if
      nstak = 0
      return
      end
      subroutine stackx(nstak)
      implicit REAL  (a-h,o-z),integer (i-n)
      logical * 1 i,j,k,l
      common/scra/v(3680),ib(3680),
     * i(14720),j(14720),k(14720),l(14720)
c
_IF1(i)      ii4=4
_IF1(v)      ii4=1
      ii=1
      do 1 n=1,nstak
      i(ii  )=i(ii4)
      i(ii+1)=j(ii4)
      i(ii+2)=k(ii4)
      i(ii+3)=l(ii4)
      ii4=ii4+4
    1 ii =ii +4
      return
      end
      function ld340(mstack,nstak,nsz170,nwbuf,lnumb,
     * v,ib,i,buf)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension nwbuf(*),lnumb(*),v(*),ib(*)
      dimension i(*),buf(*)
      dimension inter(2)
      equivalence (inter(1),xxxx)
      mstak=mstack
      nsz340=nsz170+nsz170
    2 ld340=ib(mstak)
      vv=v(mstak)
      ijkl=i(mstak)
      nw=nwbuf(ld340)+1
      ln=lnumb(ld340)+nw
      buf(ln)=vv
      if(nw.gt.nsz170)go to   1
_IF1(i)      inter(1)=0
_IF1(i)      inter(2)=ijkl
_IF1(v)      inter(1)=ijkl
_IF1(v)      inter(2)=0
      buf(ln+nsz340)=xxxx
    3 nwbuf(ld340)=nw
      mstak=mstak+1
      if(mstak.le.nstak)go to   2
      ld340=0
      return
    1 continue
      xxxx=buf(ln+nsz170)
_IF1(i)      inter(1)=ijkl
_IF1(v)      inter(2)=ijkl
      buf(ln+nsz170)=xxxx
      if(nw.lt.nsz340)go to   3
      nwbuf(ld340)=0
      mstack=mstak
      return
      end
_ELSE
      subroutine stack(buf,ijklbuf)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension buf(*)
      dimension ijklbuf(*)
INCLUDE(common/presrt)
      common/scra/v(3680),ibu(3680)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/stak/btri,mlow,nstak,iblock,mstak
      common/bufb/nwbnwb,lnklnk,buf2(1)
      mstak = 1
 20   istak = ld340(buf,ijklbuf)
      if (istak.ne.0) then
        call stopbk
        lnklnk = link(istak)
        ibn=ibasen(istak)+1
        ib = lnumb(istak) + 1
_IFN1(iv)        call pack(buf2(nsz341),lab1632,ijklbuf(ibn),nsz680)
_IF1(iv)      call fmove(ijklbuf(ibn),gout(nsz341),nsz170)
        call dcopy(nsz340,buf(ib),1,buf2,1)
        call sttout
        link(istak) = iblock
        iblock = iblock + nsz
        mstak = mstak + 1
        if (mstak.le.nstak) go to 20
      end if
      nstak = 0
      return
      end
      function ld340(buf,ijklbuf)
      implicit REAL (a-h,o-z)
INCLUDE(common/stak)
INCLUDE(common/presrt)
      common/scra/v(3680),ib(3680),itx(3680),ktx(3680)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      dimension buf(*),ijklbuf(*)
      mstak=mstack
    2 ld340=ib(mstak)
      vv=v(mstak)
      nw=nwbuf(ld340)+1
      ln=lnumb(ld340)+nw
      buf(ln)  =vv
      ijklbuf(ln+ln-1)=itx(mstak)
      ijklbuf(ln+ln  )=ktx(mstak)
*     write(6,121)ln,buf(ln),ijklbuf(ln+ln-1),ijklbuf(ln+ln  )
*121  format('ld340: ',1x,i5,f15.5,2x,2i10)
      if(nw.ge.nsz340)go to 4
      nwbuf(ld340)=nw
      mstak=mstak+1
      if(mstak.le.nstack)go to 2
      ld340=0
      return
    4 nwbuf(ld340)=0
      mstack=mstak
      return
      end
_ENDIF
_ENDEXTRACT
**==gijka.f
      subroutine gijka(cr,iconf)
_IF(parallel)
c MPP
c MPP   backchain 2-electron integrals in parallel
c MPP   ggop 2-electron integrals (ijka)
c MPP   root writes ijka ints to common ED7 (numci)
c MPP   backchain 1-electron ints on each node (complete set)
c MPP   calculate FOCK-matrices in parallel (synchronised with ijka)
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
_IF(parallel)
c MPP
      logical oipsci,opg_root
c MPP
_ENDIF
INCLUDE(common/sizes)
      dimension iconf(5,*),cr(*),ibuf2(2)
INCLUDE(common/prnprn)
      common/stoctl/khbill,khbm2e
INCLUDE(common/symci)
INCLUDE(common/iofile)
INCLUDE(common/ccntl)
INCLUDE(common/helpr)
INCLUDE(common/presrt)
      common/scra  /rmode(3),result(nd200),xmult(125),iocc(62),jap(62)
INCLUDE(common/disktl)
INCLUDE(common/orb)
      common / bufb  / nwb,lnklnk,buf2(1)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
INCLUDE(common/blksiz)
INCLUDE(common/corctl)
INCLUDE(common/prints)
INCLUDE(common/table)
      common/xynumb/xnumb(3),ynumb(3)
      equivalence (rmode(1),ic),(rmode(2),integr),(rmode(3),ne)
      equivalence (buf2(1),ibuf2(1)),(xxs1,ms1)
      character*8 zblank
      data zblank / ' '/
c... to fetch (ij/ka) and (ia) integrals into store
c... and to form fock(ia) operators for doublets +  n-2 states
c... (ij/ka) integrals loaded to ed7
      ms1 = 99999999
      top = cpulft(1)
      if(.not.oprint(28)) write (iwr,6010) top
      if (nminv1+nmin12.ne.0) then
        nuse = nintwo*m1ia + khbm2e
        call corfai(zblank)
        call vclr(cr(khbill+1),1,m2e)
        call setsto(nsz340,0,i205)
        call setsto(nsz340,0,j205)
        call setsto(nsz340,0,k205)
        call setsto(nsz340,0,l205)
c...
c... read in (ij/ka) from bucket 3 of pre sort file
c...
        iblok = link(3)
 20     if (iblok.lt.0) then
c...
c... route (ij/ka) to numci
_IF(parallel)
c                          == common ed7
_ENDIF
c...
_IFN(parallel)
          index(2) = iposun(numci)
          call wrt3s(cr(khbill+1),m2e,numci)
_ENDIF
_IF(parallel)
c MPP
c MPP   ggop integrals all needed since we'll try to build
c MPP   Fock operators in parallel
c MPP
          call pggop(7002,cr(khbill+1),m2e,'+',i205,3412)
c         call search(1,numcic)
c MPP
c         if (opg_root()) then
c           call comed7
            index(2) = iposun(numci)
            call wrt3s(cr(khbill+1),m2e,numci)
c           call pried7
c         end if
c         call broadc(7003,index(2),1)
c MPP
_ENDIF
c...
c... get j and 2j-k operators neatly lined up
c...
          mbass = khbm2e
          do 50 loop = 1 , nint
            jap(loop) = iap(loop) - m1ia
            do 40 koop = 1 , nirr
              ni = norbi(koop)
              ne = norbe(koop)
              if (ni.ne.0 .and. ne.ne.0) then
                mi = mbasi(koop)
                do 30 joop = 1 , ni
                  ii = joop + mi
                  icoul = ijkap(loop,loop,ii)
                  iexch = ijkap(loop,ii,loop)
                  call dcopy(ne,cr(icoul+1),1,cr(mbass+1),1)
_IF(cray,t3e,ibm,vax,i8)
                  call gtrian(ne,2.0d0,cr(mbass+m1ia+1),
     +                        cr(iexch+1),cr(mbass+1))
_ELSE
                  call vsmsb(cr(mbass+1),1,2.0d0,cr(iexch+1),1,
     +                       cr(mbass+m1ia+1),1,ne)
_ENDIF
                  mbass = mbass + ne
 30             continue
              end if
 40         continue
            mbass = mbass + m1ia
 50       continue
c...
c... read in (ia) from bucket 1 of pre sort file
c...
          iblok = link(1)
 60       if (iblok.lt.0) then
c...
c... fock(ia) operators for doublets
c...
            xmult(1) = 1.0d0
            index(3) = iposun(numci)
            call brtsb(index(3))
            if (nnst.ne.0) then
_IF(parallel)
c MPP
c MPP    calculated fock-operators in parallel (need carefull sync)
c MPP
             iflop = iipsci()
c MPP
_ENDIF
              do 90 ic = nstrt1 , nlast1
_IF(parallel)
c MPP
              if (oipsci()) go to 90
c MPP
_ENDIF
                call upack(iconf(1,ic),nint,iocc)
_IF(cray,t3e,convex,i8,i8drct)
                jsym = iconf(4,ic)
_ELSE
                jsym = ipad(iconf(4,ic))
_ENDIF
                ni = norbi(jsym)
                if (ni.ne.0) then
                  mbass = 2
                  do 70 loop = 1 , nint
                    jocc = iocc(loop)
                    xmult(mbass) = xnumb(jocc+1)
                    xmult(mbass+1) = ynumb(jocc+1)
                    mbass = mbass + 2
 70               continue
                  ne = norbe(jsym)
                  ne3 = ne + 3
                  mb = mbasi(jsym)
                  do 80 loop = 1 , ni
                    intern = loop + mb
                    if (iocc(intern).ne.2 .and. iro(intern).eq.jsym)
     +                  then
                      call mxmd(cr(jap(intern)+1),1,m1ia,xmult,1,1,
     +                          result,1,1,ne,nint21,1)
                      integr = iap(intern)
                      call wrtsb(rmode(1),ne3)
                    end if
 80               continue
                end if
 90           continue
            end if
c...
c... fock(ia) operators for n-2 states
c...
_IF(parallel)
c MPP
            iflop = iipsci()
c MPP
_ENDIF
            if (nmin2.ne.0) then
              do 120 ic = nstrt2 , nlast2
_IF(parallel)
c MPP
              if (oipsci()) go to 120
c MPP
_ENDIF
                call upack(iconf(1,ic),nint,iocc)
                mbass = 2
                do 100 loop = 1 , nint
                  jocc = iocc(loop)
                  xmult(mbass) = xnumb(jocc+1)
                  xmult(mbass+1) = ynumb(jocc+1)
                  mbass = mbass + 2
 100            continue
                do 110 intern = 1 , nint
                  if (iocc(intern).ne.2) then
                    ne = norbe(iro(intern))
                    if (ne.ne.0) then
                      call mxmd(cr(jap(intern)+1),1,m1ia,xmult,1,1,
     +                          result,1,1,ne,nint21,1)
                      integr = iap(intern)
                      call wrtsb(rmode(1),ne+3)
                    end if
                  end if
 110            continue
 120          continue
            end if
            call wrtsb(xxs1,1)
            call ertsb
          else
            call rdbak(iblok)
            call stopbk
            call upak8s(ibuf2,i205)
            if (odebug(40)) call wrtsor('ia    ',i205,j205,k205,l205,
     +                                  buf2,nwb)
            do 130 io = 1 , nwb
              cr(jap(j205(io))+norbas(i205(io))) = buf2(io)
 130        continue
            iblok = lnklnk
            go to 60
          end if
        else
          call rdbak(iblok)
          call stopbk
          call upak8s(ibuf2,i205)
          if (odebug(40)) call wrtsor('ijka  ',i205,j205,k205,l205,buf2,
     +                                nwb)
          do 140 io = 1 , nwb
            cr(ijkap(k205(io),l205(io),j205(io))+norbas(i205(io)))
     +        = buf2(io)
 140      continue
          iblok = lnklnk
          go to 20
        end if
      end if
      return
 6010 format (//1x,80('-')//' **** transformed integral sort called at'
     +        ,f9.2,' secs')
      end
**==gijab.f
      subroutine gijab(cr,iconf)
_IF(parallel)
c MPP
c MPP   backchain 2-electron integrals in parallel
c MPP   ggop 2-electron integrals (ijab)
c MPP   root writes ijab ints to common ED7 (numci)
c MPP   backchain 1-electron ints on each node (complete set)
c MPP   calculate FOCK-matrices in parallel (synchronised with ijab)
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/sizes)
      dimension iconf(5,*),cr(*),ibuf2(2)
INCLUDE(common/prnprn)
INCLUDE(common/symci)
      common/spew/rmode(6)
INCLUDE(common/ccntl)
      common/stoctl/khbill,khbm2e,khiaib,khiajb
INCLUDE(common/helpr)
INCLUDE(common/presrt)
      common/scra  /xmult(125),iocc(62),ikym(nd200),ikyn(nd200),jap(62)
INCLUDE(common/disktl)
INCLUDE(common/orb)
      common /bufb / nwb,lnklnk,buf2(1)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
INCLUDE(common/blksiz)
INCLUDE(common/corctl)
INCLUDE(common/table)
      common/xynumb/xnumb(3),ynumb(3)
_IF(parallel)
c MPP
      logical oipsci,opg_root
c MPP
_ENDIF
      equivalence (rmode(5),itrans),(rmode(6),ijsi)
      equivalence (ibuf2(1),buf2(1))
      character*8 zblank
      data zblank / ' '/
c... to fetch (ij/ab) (ia/jb) and (ab) integrals into store
c... and to form fock(ab) operators for doublets +  n-2 states
c... and to form (ia/jb) p and q supermatrices
      call setsto(nsz340,0,i205)
      call setsto(nsz340,0,j205)
      call setsto(nsz340,0,k205)
      call setsto(nsz340,0,l205)
      khiiab = khiaib + m2fock
      khmug = khiiab + m2fock
      khab = khmug + norbsi(1)
      nuse = khab + norbsi(1)
      call corfai(zblank)
      do 20 loop = nintp1 , ninnex
        m = ibassi(iro(loop),1) + iky(norbas(loop))
        ikyn(loop) = m
        ikym(loop) = m + khmug
 20   continue
      call vclr(cr(khbill+1),1,m2fock+m2fock)
c...
c... read in (ia/ib) from bucket 7 of pre sort file
c...
      iblok = link(7)
 30   if (iblok.lt.0) then
c...
c... read in (ii/ab) from bucket 5 of pre sort file
c...
        do 40 loop = 1 , nint
          jap(loop) = mbfock(loop) + m2fock
 40     continue
        iblok = link(5)
 50     if (iblok.lt.0) then
c...
c... read in (ab) from bucket 2 of pre sort file
c...
_IF(parallel)
c MPP
c MPP   first ggop the (ia/ib) and (ii/ab) ints
c MPP
          call pggop(7004,cr(khbill+1),m2fock+m2fock,'+',i205,3412)
c MPP
_ENDIF
          iblok = link(2)
 60       if (iblok.ge.0) then
            call rdbak(iblok)
            call stopbk
            call upak8s(ibuf2,i205)
            if (odebug(40)) call wrtsor('  ab  ',i205,j205,k205,l205,
     +                                  buf2,nwb)
            do 70 io = 1 , nwb
              cr(ikym(i205(io))+norbas(j205(io))) = buf2(io)
 70         continue
            iblok = lnklnk
            go to 60
          else
c...
c... create (ia/ib) = (ii/ab)*2-(ia/ib)
c...
_IF(cray,t3e,ibm,vax,i8)
            call gtrian(m2fock,2.0d0,cr(khiiab+1),cr(khbill+1),
     *                  cr(khiaib+1))
_ELSE
            call vsmsb(cr(khiaib+1),1,2.0d0,cr(khbill+1),1,cr(khiiab+1),
     +           1,m2fock)
_ENDIF
c...
c... fock(ab) operators for doublets
c...
            xmult(nint21) = 1.0d0
            if (nmin1.ne.0) then
              index(7) = iposun(numci)
              call brtsb(index(7))
_IF(parallel)
c MPP
              iflop = iipsci()
c MPP
_ENDIF
              do 80 ic = nstrt1 , nlast1
_IF(parallel)
c MPP
                if (oipsci()) go to 80
c MPP
_ENDIF
                call upack(iconf(1,ic),nint,iocc)
_IF(cray,t3e,convex,i8,i8drct)
                jsi = iconf(4,ic)
_ELSE
                jsi = ipad(iconf(4,ic))
_ENDIF
_IF(cray)
                call gather(nint,xmult(nintp1),ynumb(2),iocc)
                call gather(nint,xmult,xnumb(2),iocc)
_ELSE
                call dgthr(nint,ynumb(2),xmult(nintp1),iocc)
                call dgthr(nint,xnumb(2),xmult,iocc)
_ENDIF
                nk = iky(norbe(jsi)+1)
                call mxmd(cr(ibassi(jsi,1)+khiaib+1),1,norbsi(1),xmult,
     +                    1,1,cr(khab+1),1,1,nk,nint21,1)
                call wrtsb(cr(khab+1),nk)
 80           continue
              call ertsb
            end if
c...
c... fock(ab) operators for n-2 states
c...
            if (nmin2.ne.0) then
              index(6) = iposun(numci)
              call brtsb(index(6))
_IF(parallel)
c MPP
              iflop = iipsci()
c MPP
_ENDIF
              do 90 ic = nstrt2 , nlast2
_IF(parallel)
c MPP
                if (oipsci()) go to 90
c MPP
_ENDIF
                call upack(iconf(1,ic),nint,iocc)
_IF(cray)
                call gather(nint,xmult(nintp1),ynumb(2),iocc)
                call gather(nint,xmult,xnumb(2),iocc)
_ELSE
                call dgthr(nint,ynumb(2),xmult(nintp1),iocc)
                call dgthr(nint,xnumb(2),xmult,iocc)
_ENDIF
                call mxmd(cr(khiaib+1),1,norbsi(1),xmult,1,1,cr(khab+1),
     +                    1,1,norbsi(1),nint21,1)
                call wrtsb(cr(khab+1),norbsi(1))
 90           continue
              call ertsb
            end if
c...
c... read in (ij/ab) from bucket 4 of pre sort file
c...
            nuse = khiajb
            call corfai(zblank)
            call vclr(cr(khiaib+1),1,m2coul)
            iblok = link(4)
 100        if (iblok.lt.0) then
c...
c... route (ia/ib)+(ij/ab) to ed7
c...
_IFN(parallel)
              index(4) = iposun(numci)
              call wrt3s(cr(khbill+1),m2fock,numci)
              call wrt3s(cr(khiaib+1),m2coul,numci)
_ENDIF
_IF(parallel)
c MPP
c             call comed7
              index(4) = iposun(numci)
c MPP
c MPP   only root writes / need to gop ijab (ia/jb) already done above??
c MPP
              call pggop(7005,cr(khiaib+1),m2coul,'+',i205,3412)
c MPP
c             if (opg_root()) then
                      call wrt3s(cr(khbill+1),m2fock,numci)
                      call wrt3s(cr(khiaib+1),m2coul,numci)
c             end if
c             call broadc(7006,index(4),1)
c             call pried7
c MPP
_ENDIF
c...
c... read in (ia/jb) from bucket 6 of pre sort file
c...
              khiajb = khbill + m2exch
              nuse = khiajb
              call corfai(zblank)
              call vclr(cr(khbill+1),1,m2exch)
              iblok = link(6)
 110          if (iblok.lt.0) then
c...
c... route (ia/jb) integrals to ed7
c...
_IFN(parallel)
                index(8) = iposun(numci)
                call wrt3s(cr(khbill+1),m2exch,numci)
_ENDIF
_IF(parallel)
c MPP
c               call comed7
                index(8) = iposun(numci)
                call pggop(7007,cr(khbill+1),m2exch,'+',i205,3412)
c               if (opg_root())
                call wrt3s(cr(khbill+1),m2exch,numci)
c               call broadc(7008,index(8),1)
c               call pried7
c MPP
c MPP   next bit only root  (it's ED7 is the common one)
c MPP   obviously not very parallel ....
c MPP   in genz one must be aware where the integrals should be read from
c MPP
_ENDIF
c...
c... form (ia/jb)    p+q     supermatrices
c...
                if (nminv2.ne.0) then
                  khbx = khiajb + mubsi
                  nuse = khbx + nubtr
                  call corfai(zblank)
_IF(parallel)
c MPP
c                 call comed7
_ENDIF
                  index(5) = iposun(numci)
_IF(parallel)
c MPP
c                 call broadc(7009,index(5),1)
c                 if (opg_root()) then
c MPP
_ENDIF
                  call brtsb(index(5))
                  kbas = khbill
                  do 190 isi = 1 , nirr
                    ni = norbi(isi)
                    if (ni.ne.0) then
                      do 180 jsi = 1 , isi
                        if (jsi.ne.isi) then
c...      ij not totally symmetric
                          nj = norbi(jsi)
                          if (nj.ne.0) then
                            ijsi = mult(jsi,isi)
                            nijsi = norbsq(ijsi)
                            nijtr = norbtr(ijsi)
                            do 140 ii = 1 , ni
                              do 130 jj = 1 , nj
                                do 120 ksi = 1 , nirr
                                  lsi = mult(ksi,ijsi)
                                  if (lsi.gt.ksi) then
                                    nk = norbe(ksi)
                                    nl = norbe(lsi)
                                    nkl = nk*nl
                                    if (nkl.ne.0) then
                                      kbas1 = kbas + ibassq(ksi,ijsi)
                                      kbas2 = kbas + ibassq(lsi,ijsi)
                                      call symm2(cr(khiajb+1),
     +                                  cr(kbas1+1),cr(kbas2+1),nk,nl)
                                      call wrtsb(cr(khiajb+1),nkl)
                                      call anti2
     +                                  (cr(khbx+ibassi(ksi,ijsi)+1),
     +                                  cr(kbas1+1),cr(kbas2+1),nk,nl)
                                    end if
                                  end if
 120                            continue
                                call wrtsb(cr(khbx+1),nijtr)
                                kbas = kbas + nijsi
 130                          continue
 140                        continue
                          end if
                        else
c... ij  totally  symmetric
                          if (ni.eq.1) go to 190
                          do 170 ii = 2 , ni
                            ii1 = ii - 1
                            do 160 jj = 1 , ii1
                              do 150 ksi = 1 , nirr
                                nk = norbe(ksi)
                                kbas1 = kbas + ibassq(ksi,1)
                                if (nk.lt.1) go to 150
                                if (nk.ne.1) then
                                  call anti1(cr(khbx+ibastr(ksi,1)+1),
     +                              cr(kbas1+1),nk)
                                end if
                                call symm1(cr(khiajb+1),cr(kbas1+1),nk)
                                call wrtsb(cr(khiajb+1),iky(nk+1))
 150                          continue
                              call wrtsb(cr(khbx+1),norbtr(1))
                              kbas = kbas + norbsq(1)
 160                        continue
 170                      continue
                        end if
 180                  continue
                    end if
 190              continue
                  call ertsb
                end if
_IF(parallel)
c               end if
c MPP
c               call pried7
c MPP
_ENDIF
                return
              else
                call rdbak(iblok)
                call stopbk
                call upak8s(ibuf2,i205)
                if (odebug(40)) call wrtsor('iajb  ',i205,j205,k205,
     +              l205,buf2,nwb)
                do 200 io = 1 , nwb
                  jbas = iajbp(j205(io),l205(io))
                  ioop = i205(io)
                  koop = k205(io)
                  if (itrans.lt.0) then
                    m = ioop
                    ioop = koop
                    koop = m
                  end if
                  isi = iro(koop)
                  cr(jbas+mbassq(isi,ijsi)+norbas(koop)+norbas(ioop)
     +              *norbe(isi)) = buf2(io)
 200            continue
                iblok = lnklnk
                go to 110
              end if
            else
              call rdbak(iblok)
              call stopbk
              call upak8s(ibuf2,i205)
              if (odebug(40)) call wrtsor('ijab  ',i205,j205,k205,l205,
     +            buf2,nwb)
              do 210 io = 1 , nwb
                jbas = ijabp(k205(io),l205(io))
                joop = j205(io)
                m = norbas(i205(io))
                isi = iro(joop)
                if (ijsi.ne.1) then
                  loop = mbassi(isi,ijsi) + norbe(isi)*m
                else
                  loop = ibassi(isi,1) + iky(m)
                end if
                cr(jbas+loop+norbas(joop)) = buf2(io)
 210          continue
              iblok = lnklnk
              go to 100
            end if
          end if
        else
          call rdbak(iblok)
          call stopbk
          call upak8s(ibuf2,i205)
          if (odebug(40)) call wrtsor('iiab  ',i205,j205,k205,l205,buf2,
     +                                nwb)
          do 220 io = 1 , nwb
            cr(jap(k205(io))+ikyn(i205(io))+norbas(j205(io))) = buf2(io)
 220      continue
          iblok = lnklnk
          go to 50
        end if
      else
        call rdbak(iblok)
        call stopbk
        call upak8s(ibuf2,i205)
        if (odebug(40)) call wrtsor('iaib  ',i205,j205,k205,l205,buf2,
     +                              nwb)
        do 230 io = 1 , nwb
          cr(mbfock(j205(io))+ikyn(i205(io))+norbas(k205(io)))
     +      = buf2(io)
 230    continue
        iblok = lnklnk
        go to 30
      end if
      end
**==giabc.f
      subroutine giabc(cr)
_IF(parallel)
c MPP
c MPP    read backchains in parallel and gop
c MPP    let root only do the writing
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
      dimension cr(*),ibuf2(2)
INCLUDE(common/prnprn)
      common/stoctl/khbill,khbm2e,khiaib,khiajb,khiabc
INCLUDE(common/symci)
INCLUDE(common/ccntl)
INCLUDE(common/helpr)
INCLUDE(common/presrt)
INCLUDE(common/disktl)
INCLUDE(common/orb)
      common /bufb / nwb,lnklnk,buf2(1)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
INCLUDE(common/corctl)
INCLUDE(common/table)
      common/scra  /lbas(nd200),multx(8),norbsx(8)
INCLUDE(common/blksiz)
_IF(ga)
c MPP
      logical opcivect
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
      logical opg_root,oipsci,odont
c MPP
_ENDIF
      equivalence (ibuf2(1),buf2(1))
      character*8 zblank
      data zblank / ' '/
      iii = 0
      call setsto(nsz340,0,i205)
      call setsto(nsz340,0,j205)
      call setsto(nsz340,0,k205)
      call setsto(nsz340,0,l205)
_IF(ga)
c     call setsto(62*nd200,0,myext)
      call setsto(nd200*8,0,myibass3)
      call setsto(8*nd200,0,myniky)
c     call setsto(8*nd200,0,mynorbe)
_ENDIF
c... to read back (ia/bc) from pre sort file
c... to produce (ia/bc) p and q supermatrices on ed7
      if (nmin12.ne.0) then
        khbx = khiabc + mubsi
        nuse = khbx + nubtr
        call corfai(zblank)
_IFN(parallel)
        index(9) = iposun(numci)
        call brtsb(index(9))
_ENDIF
_IF(parallel)
c MPP
c       if (opg_root()) then
c          call comed7
          index(9) = iposun(numci)
          call brtsb(index(9))
c       end if
c       call broadc(7010,index(9),1)
c MPP
        iflop=iipsci()
        llfock=0
_ENDIF
        do 90 lsi = 1 , nirr
          ni = norbi(lsi)
          if (ni.ne.0) then
            kbas = khbill
_IF(cray,t3e)
            call fmove(mult(1,lsi),multx,nirr)
_ELSE
            call icopy(nirr,mult(1,lsi),1,multx,1)
_ENDIF
_IF(cray)
            call gather(nirr,norbsx,norbsq,multx)
_ELSE
            call igthr(nirr,norbsq,norbsx,multx)
_ENDIF
            do 20 kk = nintp1 , ninnex
              lbas(kk) = kbas
              kbas = kbas + norbsx(iro(kk))
 20         continue
            l2iabc = niky(lsi)
            do 80 ll = 1 , ni
              iii = iii + 1
              call vclr(cr(khbill+1),1,l2iabc)
c... read in back chain of (ia/bc) integrals
              iblok = link(iii+nirr8)
 30           if (iblok.lt.0) then
_IF(parallel)
c MPP
              call pggop(7011,cr(khbill+1),l2iabc,'+',i205,3412)
c             if (.not.opg_root()) go to 80
c MPP
_ENDIF
                do 60 kk = nintp1 , ninnex
                  ksi = iro(kk)
_IF(ga)
                if (oipsci().and..not.opcivect) goto 60
                 if (myext(iii,kk).ne.1) call caserr('error myext')
c                 myext(iii,kk)=1
_ENDIF
                  klsi = multx(ksi)
                  kbas = lbas(kk)
                  if (klsi.ne.1) then
c... kl not totally symmetric
                    do 40 jsi = 1 , nirr
                      isi = mult(jsi,klsi)
                      if (isi.gt.jsi) then
                        nj = norbe(jsi)
                        nk = norbe(isi)
                        nij = nj*nk
                        if (nij.ne.0) then
                          kbas1 = kbas + ibassq(jsi,klsi)
                          kbas2 = kbas + ibassq(isi,klsi)
                          call symm2(cr(khiabc+1),cr(kbas1+1),
     +                               cr(kbas2+1),nj,nk)
                          call wrtsb(cr(khiabc+1),nij)
_IF(ga)
                  myniky(iii,lsi)=myniky(iii,lsi)+nij
                  myibass3(ksi,iii)=myibass3(ksi,iii)+nij
_ENDIF
                          call anti2(cr(khbx+ibassi(jsi,klsi)+1),
     +                               cr(kbas1+1),cr(kbas2+1),nj,nk)
                        end if
                      end if
 40                 continue
                  else
c... kl totally symmetric
                    do 50 jsi = 1 , nirr
                      nj = norbe(jsi)
                      kbas1 = kbas + ibassq(jsi,1)
                      if (nj.lt.1) go to 50
                      if (nj.ne.1) then
                        call anti1(cr(khbx+ibastr(jsi,1)+1),cr(kbas1+1),
     +                             nj)
                      end if
                      call symm1(cr(khiabc+1),cr(kbas1+1),nj)
                      call wrtsb(cr(khiabc+1),iky(nj+1))
_IF(ga)
                  myniky(iii,lsi)=myniky(iii,lsi)+iky(nj+1)
                  myibass3(ksi,iii)=myibass3(ksi,iii)+iky(nj+1)
_ENDIF
 50                 continue
                  end if
                  call wrtsb(cr(khbx+1),norbtr(klsi))
_IF(ga)
                  myniky(iii,lsi)=myniky(iii,lsi)+norbtr(klsi)
                  myibass3(ksi,iii)=myibass3(ksi,iii)+norbtr(klsi)
_ENDIF
 60             continue
              else
                call rdbak(iblok)
                call stopbk
                call upak8s(ibuf2,i205)
                if (odebug(40)) call wrtsor('iabc  ',i205,j205,k205,
     +              l205,buf2,nwb)
                do 70 ioo = 1 , nwb
                  io = i205(ioo)
                  jo = j205(ioo)
                  ko = k205(ioo)
                  ksi = iro(ko)
                  nk = norbe(ksi)
                  kk = norbas(ko)
                  valu = buf2(ioo)
                  if (io.ne.jo) then
c... process integral in its (ba/ci) form
                    cr(lbas(io)+mbassq(ksi,multx(iro(io)))+norbas(jo)
     +                *nk+kk) = valu
                  end if
c... process integral in its (ab/ci) form
                  cr(lbas(jo)+mbassq(ksi,multx(iro(jo)))+norbas(io)
     +              *nk+kk) = valu
 70             continue
                iblok = lnklnk
                go to 30
              end if
 80         continue
          end if
 90     continue
        call ertsb
      end if
_IF(parallel)
c MPP
      call pried7
      do  iij=1,nintp1-1
         ioff=0
         do kij=1,nirr
            myib=myibass3(kij,iij)
            myibass3(kij,iij)=ioff
            ioff=ioff+myib
         enddo
      enddo
c MPP
_ENDIF
      return
      end
_EXTRACT(gabcd,mips4)
**==gabcd.f
      subroutine gabcd(cr,i127)
_IF(parallel)
c MPP
c MPP   backchain in parallel for abcd integrals on presort file
c MPP   post-sort in parallel as well; then gop write by root only
c MPP   this seems easiest ..
c MPP
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer ibuf2
_ELSEIF(convex,i8drct)
      integer *8 ibuf2
_ELSE
      REAL ibuf2
_ENDIF
      dimension cr(*),i127(*),ibuf2(2)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
INCLUDE(common/prints)
      common/stoctl/khbill
INCLUDE(common/symci)
INCLUDE(common/ccntl)
INCLUDE(common/helpr)
INCLUDE(common/presrt)
      common/lsort/ilif(256),lbas(nd200),llbas(nd200,8)
      common/dopey/n1271,n127p1,mkk,maxblk,iksi,nine,mmm127,nbuck
     *,khb127,nrem,ndimsi(8),nremsi(8),m127si(8)
INCLUDE(common/orb)
      common /bufb / nwb,lnklnk,buf2(1)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
_IF(ibm,vax)
     * , ibm205(3412),jbm205(3412)
_ENDIF
INCLUDE(common/corctl)
INCLUDE(common/table)
      common/scra/mnumb(2016),mink(2016),mwbuf(2016)
      common/stak2/btri,mlowww(2),iblock
INCLUDE(common/blksiz)
INCLUDE(common/iofile)
_IF(parallel)
c MPP
      logical opg_root
c MPP
_ENDIF
      equivalence (ibuf2(1),buf2(1))
      character*8 zblank
      data zblank / ' '/
c

cjmht kbas & kkbas need to be initialised to 0 each time this routine is called
cjmht However, according to Marc a bug in gfortran means that the initialisation
cjmht statements are optimised away and so we also add additional data statements
cjmht to ensure they are set. These can be removed once the bug is fixed.
      
      data kkbas/0/
      data kbas/0/
      kkbas = 0
      kbas = 0

      nav = lenwrd() 
c...
_IF(parallel)
c MPP
c MPP  open seperate postsort-files (routines must be checked)
c MPP  check stack2 and fabcd
c MPP
c     call comed7
      call pried7
c MPP
_ENDIF
c...  to produce (ab/cd) p and q supermatrices by symmetry block
c...  each in sub-blocks of  n127 * n127 dimension
c...  the right hand border blocks have dimension  n127 * nrem
c...
      if (nmin2.ne.0) then
        lhb127 = khbill + nubsi
        khbil2 = khbill*nav
        nuse = lhb127 + 510
        call corfai(zblank)
c... determine blocking factor
        ngott = (ngot-lhb127)/510
        if (ngott.gt.2016) ngott = 2016
c... defines a limit of 63*(63+1)/2 sub blocks
        nk = nubsi - 1
        do 20 loop = 1 , 5
          if (n127.gt.maxsq) call caserr('dimensioning error in gabcd')
          mk = nk/n127 + 1
          if ((mk*(mk+1)/2).le.ngott) go to 30
          n127 = n127 + nincr
 20     continue
 30     if (n127.gt.nubsi) n127 = nubsi
        nl = n127 + 1
        n127p1 = nl
        n1271 = nl - 2
        if (nl.lt.mk) nl = mk
        do 40 loop = 1 , nl
          ilif(loop) = kbas
          kbas = kbas + n127
 40     continue
        n127sq = ilif(n127p1)
        nuse = lhb127 + n127sq
        call corfai(zblank)
        do 50 loop = 1 , nubsi
          i127(loop+khbil2) = (loop-1)/n127 + 1
 50     continue
        nbuck = iky(mk+1)
        khb127 = lhb127
        do 60 loop = 1 , nbuck
          mnumb(loop) = lhb127
          lhb127 = lhb127 + 510
 60     continue
        maxblk = 0
        call vclr(cr(khb127+1),1,nbuck*510)
c...
c... generate (ab/cd) p supermatrices
c...
        nine = 9
        do 130 iksi = 1 , nirr
          kbas = 0
          if (iksi.ne.1) then
            do 80 loop = 1 , nirr
              koop = mult(loop,iksi)
              if (koop.gt.loop) then
                nk = norbe(koop)
                if (nk.ne.0) then
                  nl = norbe(loop)
                  mk = nbase(koop)
                  do 70 koop = 1 , nk
                    mmm = koop + mk
                    llbas(mmm,iksi) = kbas
                    lbas(mmm) = kbas
                    kbas = kbas + nl
 70               continue
                end if
              end if
 80         continue
          else
            do 100 loop = 1 , nirr
              nk = norbe(loop)
              if (nk.ne.0) then
                mk = nbase(loop)
                do 90 koop = 1 , nk
                  lbas(koop+mk) = kbas
                  kbas = kbas + koop
 90             continue
              end if
 100        continue
          end if
          if (kbas.ne.0) then
            mkk = i127(kbas+khbil2)
            ndimsi(iksi) = mkk
            nrem = kbas - ilif(mkk)
            if (nrem.eq.0) nrem = n127
            nremsi(iksi) = nrem
            mmm127 = ilif(nrem+1)
            m127si(iksi) = mmm127
            nbuck = iky(mkk+1)
            call setsto(nbuck,-1,mink)
            call setsto(nbuck,0,mwbuf)
            call setsto(nsz340,0,i205)
            call setsto(nsz340,0,j205)
            call setsto(nsz340,0,k205)
            call setsto(nsz340,0,l205)
_IF(parallel)
c **** MPP
            call closbf2(0)
            call setbfb(-1)
c **** MPP
_ENDIF
            iblock = 0
c... process this symmetry back chain on pre-sort file
            iblok = link(iksi+8)
 110        if (iblok.lt.0) then
              call fabcd(cr)
            else
              call rdbak(iblok)
              call stopbk
              call upak8s(ibuf2,i205)
              if (odebug(40)) call wrtsor('abcd  ',i205,j205,k205,l205,
     +            buf2,nwb)
              do 120 ioo = 1 , nwb
_IF(ibm,vax)
               ibm205(ioo) = norbas(k205(ioo)) + lbas(i205(ioo))
_ELSE
               i205(ioo) = norbas(k205(ioo)) + lbas(i205(ioo))
_ENDIF
                jo = j205(ioo)
                lo = l205(ioo)
                j = norbas(min(jo,lo)) + lbas(max(jo,lo))
_IF(ibm,vax)
                jbm205(ioo) = j
_ELSE
                j205(ioo) = j
_ENDIF
 120          continue
              call stack2(nwb,cr,i127)
              iblok = lnklnk
              go to 110
            end if
          end if
 130    continue
c...
c... generate (ab/cd) q supermatrices
c...
_IF(parallel)
c **** MPP
      call closbf2(0)
      call setbfb(-1)
_ENDIF
        nine = 17
        do 190 iksi = 1 , nirr
          if (iksi.ne.1) then
            kkbas = norbtr(iksi)
            if (kkbas.eq.0) go to 190
            nrem = nremsi(iksi)
            mkk = ndimsi(iksi)
            mmm127 = m127si(iksi)
            do 140 loop = nintp1 , ninnex
              lbas(loop) = llbas(loop,iksi)
 140        continue
          else
            do 160 loop = 1 , nirr
              nk = norbe(loop)
              if (nk.gt.1) then
                mk = nbase(loop)
                do 150 koop = 2 , nk
                  lbas(koop+mk) = kkbas
                  kkbas = kkbas + koop - 1
 150            continue
              end if
 160        continue
            if (kkbas.eq.0) go to 190
            mkk = i127(kkbas+khbil2)
            nrem = kkbas - ilif(mkk)
            if (nrem.eq.0) nrem = n127
            mmm127 = ilif(nrem+1)
          end if
          nbuck = iky(mkk+1)
          call setsto(nbuck,-1,mink)
          call setsto(nbuck,0,mwbuf)
_IF(parallel)
c **** MPP
          call closbf2(0)
          call setbfb(-1)
c **** MPP
_ENDIF
          iblock = 0
c... process this symmetry back chain on pre-sort file
          iblok = link(iksi+8)
 170      if (iblok.lt.0) then
            call fabcd(cr)
          else
            call rdbak(iblok)
            call stopbk
            call upak8s(ibuf2,i205)
            if (odebug(40)) call wrtsor('abcd  ',i205,j205,k205,l205,
     +                                  buf2,nwb)
            nnnn = 0
            do 180 ioo = 1 , nwb
              io = i205(ioo)
              jo = j205(ioo)
              ko = k205(ioo)
              lo = l205(ioo)
              if (io.ne.ko) then
                val = buf2(ioo)
                if (lo.lt.jo) then
                  j = norbas(lo) + lbas(jo)
                else if (lo.eq.jo) then
                  go to 180
                else
                  j = norbas(jo) + lbas(lo)
                  val = -val
                end if
                i = norbas(ko) + lbas(io)
                nnnn = nnnn + 1
_IF(ibm,vax)
                ibm205(nnnn) = i
                jbm205(nnnn) = j
_ELSE
                i205(nnnn) = i
                j205(nnnn) = j
_ENDIF
                buf2(nnnn) = val
              end if
 180        continue
            if (nnnn.ne.0) call stack2(nnnn,cr,i127)
            iblok = lnklnk
            go to 170
          end if
 190    continue
        if(.not.oprint(28)) write (iwr,6010) n127 , maxblk
      end if
_IF(parallel)
c MPP
      call pried7
c MPP
_ENDIF
      return
 6010 format (/' blocking factor for external supermatrices=',i5,
     +        //' psort file',i8,' blocks long')
      end
_ENDEXTRACT
      subroutine stack2(nnn,buf,i127)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(convex,i8drct)
      integer *8 buf,buf2,ij
_ELSEIF(ibm,vax)
      integer *2 mmmm,nnnn
      dimension mmmm(4),nnnn(2)
_ELSEIF(i8)
      integer buf,buf2,ij
_ENDIF
      dimension buf(*),i127(*)
INCLUDE(common/sizes)
INCLUDE(common/cdcryi)
      common/stoctl/khbill
INCLUDE(common/helpr)
INCLUDE(common/presrt)
_IF(ibm,vax)
      common/junk/isp(4,3412),i205(3412),j205(3412)
_ELSE
      common/junk/i205(3412),j205(3412)
_ENDIF
      common/bufb/nword(2),buf2(1)
      common/lsort/ilif(256),lbas(nd200),llbas(nd200,8)
      common/scra/mnumb(2016),mink(2016),mwbuf(2016)
      common/bufc/nwbnwb,lnklnk,buf3(510)
      common/stak2/btri,mlowww(2),iblock
_IF(ibm,vax)
      equivalence (mmmm(1),xxxx),(nnnn(1),ij)
      data mmmm,nnnn/4*0,2*0/
_ENDIF
_IF(t3e)
_IFN(t3d)
      integer shift
      shift(ival,i) = ishftc(ival,i,-64)
_ENDIF
_ENDIF
      khbil2=khbill*lenwrd()
      do 99 loop=1,nnn
      i=max(i205(loop),j205(loop))
      j=min(i205(loop),j205(loop))
      ni=     i127(i+khbil2)
      nj=    i127(j+khbil2)
      nij=iky(ni)+nj
      ij=ilif(i-ilif(ni))+j-ilif(nj)
      nwb=mwbuf(nij)
      ni=nwb/4
      ibas=mnumb(nij)
      nwb=nwb+1
      buf(ibas+nwb)=buf2(loop)
_IF(ibm,vax)
      nj=ni*4
      nj=nwb-nj
      jbas=ibas+ni
      jbas49=jbas+409
      xxxx=buf(jbas49)
_IF(ibm)
      mmmm(nj)=nnnn(2)
_ELSE
      mmmm(nj)=nnnn(1)
_ENDIF
      buf(jbas49)=xxxx
_ELSEIF(cray,t3e,t3d)
      nj=shift(ni,2)
      nj=nwb-nj
      jbas49=ibas+ni+409
      go to (1,2,3,4), nj
_IF(cray,t3d)
    4 i127(jbas49)=ij.or.i127(jbas49)
_ELSE
    4 i127(jbas49) = or(ij,i127(jbas49))
_ENDIF
      go to 101
_IF(cray,t3d)
 3    i127(jbas49)=shift(ij,16).or.i127(jbas49)
_ELSE
 3    i127(jbas49) = or(shift(ij,16),i127(jbas49))
_ENDIF
      go to 101
_IF(cray,t3d)
   2  i127(jbas49)=shift(ij,32).or.i127(jbas49)
_ELSE
   2  i127(jbas49) = or(shift(ij,32),i127(jbas49))
_ENDIF
      go to 101
 1    i127(jbas49)=shift(ij,48)
 101  continue
_ELSEIF(bits8)
      nj=ni*4
      nj=nwb-nj
      jbas49=ibas+ni+409
      go to (1,2,3,4), nj
_IF(littleendian)
   4  call pad8(ij,padij)
      call shiftl8(padij,48,tmp8)
      call dor8(buf(jbas49),tmp8,buf(jbas49))
      go to 101
   3  call pad8(ij,padij)
      call shiftl8(padij,32,tmp8)
      call dor8(buf(jbas49),tmp8,buf(jbas49))
      go to 101
   2  call pad8(ij,padij)
      call shiftl8(padij,16,tmp8)
      call dor8(buf(jbas49),tmp8,buf(jbas49))
      go to 101
   1  call pad8(ij,buf(jbas49))
_ELSE
   4  call pad8(ij,padij)
      call dor8(buf(jbas49),padij,tmpij)
      buf(jbas49) = tmpij
      go to 101
   3  call pad8(ij,padij)
      call shiftl8(padij,16,tmp8)
      call dor8(buf(jbas49),tmp8,tmpij)
      buf(jbas49) = tmpij
      go to 101
   2  call pad8(ij,padij)
      call shiftl8(padij,32,tmp8)
      call dor8(buf(jbas49),tmp8,tmpij)
      buf(jbas49) = tmpij
      go to 101
   1  call pad8(ij,padij)
      call shiftl8(padij,48,buf(jbas49))
_ENDIF
  101   continue
_ELSE
      nj=ni*4
      nj=nwb-nj
      jbas49=ibas+ni+409
      go to (1,2,3,4), nj
_IF(i8,i8drct)
_IF(littleendian)
   4  buf(jbas49)=ior(buf(jbas49), ishft(ij,48) )
      go to 101
   3  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,32) )
      go to 101
   2  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,16) )
      go to 101
   1  buf(jbas49)=ij
_ELSE
   4  buf(jbas49)=ior(buf(jbas49), ij)
      go to 101
   3  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,16) )
      go to 101
   2  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,32) )
      go to 101
   1  buf(jbas49)=ishft(ij,48)
_ENDIF
_ELSE
_IF(convex)
   4  buf(jbas49)=ior(buf(jbas49), ij)
      go to 101
   3  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,16) )
      go to 101
   2  buf(jbas49)=ior(buf(jbas49) ,ishft(ij,32) )
      go to 101
   1  buf(jbas49)= ishft(ij,48)
_ELSEIF(littleendian)
   4  buf(jbas49)=dor(buf(jbas49), shiftl(pad(ij),48) )
      go to 101
   3  buf(jbas49)=dor(buf(jbas49) ,shiftl(pad(ij),32) )
      go to 101
   2  buf(jbas49)=dor(buf(jbas49) ,shiftl(pad(ij),16) )
      go to 101
   1  buf(jbas49)= pad(ij)
_ELSE
   4  buf(jbas49)=dor(buf(jbas49), pad(ij))
      go to 101
   3  buf(jbas49)=dor(buf(jbas49) ,shiftl(pad(ij),16) )
      go to 101
   2  buf(jbas49)=dor(buf(jbas49) ,shiftl(pad(ij),32) )
      go to 101
   1  buf(jbas49)= shiftl(pad(ij),48)
_ENDIF
_ENDIF
  101 continue
_ENDIF
      if(nwb.eq.408)goto 102
      mwbuf(nij)=nwb
      go to 99
102   call stopbk2
      lnklnk=mink(nij)
      call fmove(buf(ibas+1),buf3,510)
      call sttout2
      mink(nij)=iblock
      iblock=iblock+1
      mwbuf(nij)=0
 99   continue
      return
      end
      subroutine fabcd(cr)
_IF(parallel)
c MPP
c MPP   read backchain in parallel
c MPP   let only root write
c MPP
_ENDIF
      implicit REAL (a-h,o-z),integer  (i-n)
_IF(convex,i8drct)
      integer *8 ibuf3
_ELSEIF(ibm,vax)
      integer *2 ibuf3
_ELSE
      REAL ibuf3
_ENDIF
INCLUDE(common/sizes)
      dimension cr(*)
INCLUDE(common/prnprn)
INCLUDE(common/symci)
INCLUDE(common/presrt)
      common/lsort/ilif(256),lbas(nd200),llbas(nd200,8)
      common/dopey/n1271,n127p1,mkk,maxblk,iksi,nine,mmm127,nbuck
     *,khb127,nrem
      common/scra/mnumb(2016),mink(2016),mwbuf(2016)
_IFN1(iv)     * ,ij408(408)
INCLUDE(common/disktl)
INCLUDE(common/cdcryz)
INCLUDE(common/cdcryi)
      common/bufc/nwbnwb,lnklnk,buf3(408),ibuf3(1)
      common/stak2/btri,mlowww(2),iblock
_IF(parallel)
c MPP
      logical opg_root,odo,oipsci
c MPP
_ENDIF
c...
c... first clear up bucketing to sort file
c...
      do 7777 loop=1,nbuck
      nwb=mwbuf(loop)
      if(nwb)7776,7775,7776
7775  mwbuf(loop)=408
      goto 7777
7776  call stopbk2
      lnklnk=mink(loop)
_IF(cray,ibm,vax)
      call fmove(cr(mnumb(loop)+1),buf3,510)
_ELSE
      call vmov(cr(mnumb(loop)+1),1,buf3,1,510)
_ENDIF
      call sttout2
      mink(loop)=iblock
      iblock=iblock+1
7777  continue
      if(maxblk.lt.iblock)maxblk=iblock
      ndim=n127sq
      nsize=n127
      nsize1=n1271
      nbuk=0
      mwb=iposun(numci)
      index(nine+iksi)=mwb
_IF(parallel)
c MPP
c     call broadc(7012,index(nine+iksi),1)
c     if (opg_root())
c    *call brtsb(mwb)
      call brtsb(mwb)
c MPP
_ELSE
      call brtsb(mwb)
_ENDIF
      call stopbk2
_IF(parallel)
      iflop=iipsci()
_ENDIF
      do 11 ir=1,mkk
_IF(parallel)
         odo=.not.oipsci()
_ENDIF
      if(ir.ne.mkk)goto 2345
      ndim=mmm127
      nsize=nrem
      nsize1=nsize-1
2345  do 11 ic=1,ir
c...
c... read in back chain from sort file
c... and generate sub block of supermatrix
c...
      call vclr(cr(khb127+1),1,ndim)
      nbuk=nbuk+1
      nwb=mwbuf(nbuk)
      iblok=mink(nbuk)
777   if(iblok)999,888,888
888   call rdbak2(iblok)
      call stopbk2
_IF(ibm,vax)
      do 33 iwb=1,nwb
      ij=ibuf3(iwb)
      ij =ij +khb127
   33 cr(ij)=cr(ij)+buf3(iwb)
_ELSE
      call unpack(ibuf3,16,ij408,408)
      do 33 loop=1,nwb
      ijkl=ij408(loop)+khb127
  33  cr(ijkl)=cr(ijkl)+buf3(loop)
_ENDIF
      if(odebug(40)) then
       write(6,77776)nwb,iblok
77776  format(/1x,'fabcd, nwb, iblok = ',2i7)
       write(6,*)'buf3'
       write(6,77777) (buf3(mfg),mfg=1,nwb)
77777  format(1x,7f12.6)
       write(6,*)'ij408'
       write(6,77778) (ij408(mfg),mfg=1,nwb)
77778  format(1x,12i7)
      endif
      nwb=408
      iblok=lnklnk
      goto 777
999   continue
_IF(parallel)
c MPP
c MPP    gop / scratch buffer perhaps a bit small ??
c MPP
      call pggop(7013,cr(khb127+1),ndim,'+',buf3,500)
c     if (.not.opg_root()) go to 11
      if (.not.odo) go to 11
c MPP
_ENDIF
      if(ic.ne.ir)goto 1111
      if(nsize1)3456,1111,3456
c...
c... diagonal block   ---   complete the square
c...
3456  jwb=khb127
      kwb=khb127
      do 44 loop=1,nsize1
      jwb=jwb+n127p1
      call dcopy(nsize-loop,cr(jwb),n127,cr(kwb+loop+1),1)
44    kwb=kwb+n127
c...
c... output sub block to ed7
c...
1111  if(odebug(40).and.ndim.gt.0) then
       call dbgvec('fabcd   ',cr(khb127+1),ndim)
      endif
      call wrtsb(cr(khb127+1),ndim)
11    continue
_IF(parallel)
c     if (opg_root())
c    *call ertsb
      call ertsb
_ELSE
      call ertsb
_ENDIF
      return
      end
      subroutine brtsb(ibll)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/disktl)
c     initialise the process of writing a contiguous set
c     of data to numci
      common/disksa/bof(511),no,iblo
      iblo=ibll
      no=0
      call search(iblo,numci)
      return
      end
      subroutine wrtsb(q,nword)
      implicit REAL  (a-h,o-z),integer  (i-n)
c.... write full records (see brtsb)
      dimension q(*)
INCLUDE(common/disktl)
_IF(i8drct)
      integer*8 bof,q
_ENDIF
      common/disksa/bof(511),no,iblo
       if(nword)77,88,77
77    k = nword
      j = 1
c.... if no elements in bof ,try to write directly first
      if (no.eq.0) go to 20
c.... fill bof
10    nn = min(k,511-no)
_IF(i8drct)
      do loop =1,nn
       bof(no+loop) = q(j+loop-1)
      enddo
_ELSE
      call fmove(q(j),bof(no+1),nn)
_ENDIF
      no = no + nn
      if (no.lt.511) return
c.... buffer full => write out / check positioning
_IF(i8drct)
      call puti8(bof,511,numci)
_ELSE
      call put(bof,511,numci)
_ENDIF
      iblo=iblo+1
      no = 0
      k = k - nn
      j = j + nn
c.... check if it is worthwhile to write some complete records
20    if (k-511) 10,30,30
_IF(i8drct)
30    call puti8(q(j),511,numci)
_ELSE
30    call put(q(j),511,numci)
_ENDIF
      iblo=iblo+1
      k=k-511
      j=j+511
      if(k.ge.511)goto 30
      if(k)88,88,10
88    return
      end
      subroutine ertsb
      implicit REAL  (a-h,o-z),integer  (i-n)
c.... close writing of a contiguous array
_IF(i8drct)
      integer*8 bof
_ENDIF
      common/disksa/bof(511),no,iblo
INCLUDE(common/disktl)
      if (no) 2,1,2
c.... output last (partial) buffer
_IF(i8drct)
2     call puti8(bof,no,numci)
_ELSE
2     call put(bof,no,numci)
_ENDIF
1     return
      end
      subroutine whersb(m)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/spcoef/nwwwww(9),iwwwww(9)
_IF(i8drct)
      integer*8 bof
_ENDIF
      common/disksa/bof(511),no,iblo
      nwwwww(m)=no
      if(no)888,999,888
999   iwwwww(m)=iblo
      goto 777
888   iwwwww(m)=iblo+1
777   return
      end
      subroutine termsb(jconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
c
c...  only called from eblock (sorting of matrix-el for out store)
c...  so common lsort has to agree with that
c...  boldly changed ibang to integer and added other copy array
c
      integer ibang
_IF(cray,t3e,i8)
      integer   jconf
_ELSEIF(convex,i8drct)
      integer *8   jconf
_ELSE
      REAL  jconf
_ENDIF
INCLUDE(common/disktl)
_IF(i8drct)
      integer*8 scr,bof
_ENDIF
      common/disksa/bof(511),no,iblo
      common/lsort/ibang(255),nwr(256),icnow,ko,kblo,iw1
     1      ,idum,scr(511)
      dimension jconf(*)
_IF(linux)
      external fget
_ENDIF
c...
c... first position the file
c...
      no=ko-1
      iblo=kblo-1
      if(no)777,888,666
777   no=510
666   call search(iblo,numci)
c...   get block iblo
_IF(i8drct)
      call find(numci)
      call geti8(bof,m)
_ELSE
      call fget(bof,m,numci)
_ENDIF
888   call search(iblo,numci)
c...
c... now output the sorted coefficients
c... ibang is a integer*8 array preferably
c...  ibang can be straight integer now / scr for copying (jvl 93)
c...
       do 30 loop=1,icnow
       iiint = ibang(loop)+1
30     call wrtsb(jconf(iiint),nwr(loop))
c...
c...   clear up last block
c...
       if(no)333,999,333
_IF(i8drct)
333    call find(numci)
       call geti8(scr,m)
       do loop =1,m-no
        bof(no+loop) = scr(no+loop)
       enddo
       call wrt3i8(bof,m,iblo,numci)
_ELSE
333    call fget(scr,m,numci)
       call fmove(scr(no+1),bof(no+1),m-no)
       call wrt3(bof,m,iblo,numci)
_ENDIF
       iblo=iblo+1
c...    getsb is just for proper position
999    call getsb(scr,iw1)
       return
      end
      subroutine possb(mo,mblo)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(i8drct)
      integer*8 bof
_ENDIF
INCLUDE(common/disktl)
_IF(linux)
      external fget
_ENDIF
      common/disksa/bof(511),no,iblo
      no=mo
      if(iblo.eq.mblo)goto 999
      iblo=mblo
      if(no)888,999,888
888   if(iposun(numci).ne.iblo-1)call search(iblo-1,numci)
_IF(i8drct)
      call find(numci)
      call geti8(bof,m)
_ELSE
      call fget(bof,m,numci)
_ENDIF
999   return
      end
      subroutine skipsb(nw)
      implicit REAL  (a-h,o-z), integer  (i-n)
_IF(i8drct)
      integer*8 bof
_ENDIF
INCLUDE(common/disktl)
_IF(linux)
      external fget
_ENDIF
      common/disksa/bof(511),no,iblo
       if(no)4,3,4
3       no511=nw
        goto 1
4      no=no+nw
      no511=no-511
      if(no511)999,1,1
1     manx=no511/511
      iblo=iblo+manx
      no=no511-manx*511
      if(no)2,999,2
2     if(iposun(numci).ne.iblo)call search(iblo,numci)
_IF(i8drct)
      call find(numci)
      call geti8(bof,m)
_ELSE
      call fget(bof,m,numci)
_ENDIF
      iblo=iblo+1
999   return
      end
      subroutine getsb(q,nword)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(i8drct)
      integer *8 q, bof
_ENDIF
      dimension q(*)
      common/disksa/bof(511),no,iblo
INCLUDE(common/disktl)
_IF(linux)
      external fget
_ENDIF
      k=nword
      j=1
      if(no)30,20,30
30    nn=min(k,511-no)
_IF(i8drct)
      do loop =1, nn
       q(j+loop-1) = bof(no+loop)
      enddo
_ELSE
      call fmove(bof(no+1),q(j),nn)
_ENDIF
      k=k-nn
      j=j+nn
      no=no+nn
      if(no.eq.511)no=0
20    if(iposun(numci).ne.iblo)call search(iblo,numci)
       k511=k-511
      if(k511)10,40,40
_IF(i8drct)
40    call find(numci)
      call geti8(q(j),m)
_ELSE
40    call fget(q(j),m,numci)
_ENDIF
      iblo=iblo+1
      k=k511
      j=j+511
      goto 20
10    if(k)50,999,50
_IF(i8drct)
50    call find(numci)
      call geti8(bof,m)
_ELSE
50    call fget(bof,m,numci)
_ENDIF
      iblo=iblo+1
      no=k
_IF(i8drct)
      do loop =1, k
       q(j+loop-1) = bof(loop)
      enddo
_ELSE
      call fmove(bof,q(j),k)
_ENDIF
999   return
      end
      subroutine wrtsc(cr,a,nword)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension a(*)
      dimension cr(*)
        common/bufb /ifi,ifo,iiiblk
INCLUDE(common/corctl)
      if(nword)77,88,77
77    n=ifi+nword
      m=ngot-n
      if(m)888,999,999
888   jword=ifi-ifo
      nblk=jword/511
      if(nblk)666,777,666
777   call caserr('invalid block specified on cifile')
c... use ed8 as overspill area
666   iword=nblk*511
      call wrt3(cr(ifo+1),iword,iiiblk,9)
      jword=jword-iword
      iiiblk=iiiblk+nblk
      call fmove(cr(ifo+iword+1),cr(ifo+1),jword)
      ifi=ifo+jword
      goto 77
999   call fmove(a,cr(ifi+1),nword)
      ifi=n
88    return
      end
      subroutine ertsc(cr)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension cr(*)
      common/bufb /ifi,ifo,iiiblk,jjjblk,buf(511)
_IF(linux)
      external fget
_ENDIF
      nblk=iiiblk-1
      if(nblk)888,999,888
c... clear up overspill on ed8
888   call search(1,9)
      do 1 loop=1,nblk
      call fget(buf,n511,9)
1     call wrtsb(buf,n511)
999   call wrtsb(cr(ifo+1),ifi-ifo)
      return
      end
      subroutine brtsc(iblock)
      implicit REAL  (a-h,o-z),integer (i-n)
      common/bufb /bof(511),no,iblo
      iblo=iblock
      no=0
      return
      end
      subroutine getsc(q,nword)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension q(*)
INCLUDE(common/disktl)
      common/bufb /bof(511),no,iblo
_IF(linux)
      external fget
_ENDIF
      k=nword
      j=1
      if(no)30,20,30
30    nn=min(k,511-no)
      call fmove(bof(no+1),q(j),nn)
      k=k-nn
      j=j+nn
      no=no+nn
      if(no.eq.511)no=0
20    if(iposun(numci).ne.iblo)call search(iblo,numci)
      k511=k-511
      if(k511)10,40,40
40    call fget(q(j),m,numci)
      iblo=iblo+1
      k=k511
      j=j+511
      goto 20
10    if(k)50,999,50
50    call fget(bof,m,numci)
      iblo=iblo+1
      no=k
      call fmove(bof,q(j),k)
999   return
      end
      subroutine skipsc(nw)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/disktl)
      common/bufb /bof(511),no,iblo
_IF(linux)
      external fget
_ENDIF
       if(no)4,3,4
3       no511=nw
       goto 1
4      no=no+nw
      no511=no-511
      if(no511)999,1,1
1     manx=no511/511
      iblo=iblo+manx
      no=no511-manx*511
      if(no)2,999,2
2     if(iposun(numci).ne.iblo)call search(iblo,numci)
      call fget(bof,m,numci)
      iblo=iblo+1
999   return
      end
      subroutine symm2(r,a,b,nk,nl)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*),a(*),b(*)
_IF1(civ)      call dagger(nl,nk,b,nl,r,nk)
_IFN1(civ)      call rmtran(b,nl,r,nk,nl,nk)
      call vadd(r,1,a,1,r,1,nk*nl)
      return
      end
      subroutine anti2(r,a,b,nk,nl)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*),a(*),b(*)
_IF1(civ)      call dagger(nl,nk,b,nl,r,nk)
_IFN1(civ)      call rmtran(b,nl,r,nk,nl,nk)
      call vsub(a,1,r,1,r,1,nk*nl)
      return
      end
      subroutine ver_dircta(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dircta.m,v $
     +     "/
      data revision /"$Revision: 6183 $"/
      data date /"$Date: 2010-08-13 21:08:12 +0200 (Fri, 13 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
_IF(ga)
      subroutine switchnum(icr)
      implicit REAL  (a-h,o-z),integer  (i-n)  
INCLUDE(common/disktl) 
      save numciold
      if (icr.eq.0) then
         numciold=numci
         numci=numcz
      endif
      if (icr.eq.1) then
         numci=numciold
      endif
      return
      end

      subroutine setmynorbe
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/table)
INCLUDE(common/orb)
      logical opcivect
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
      logical oipsci
      iflop=iipsci()
      iii=0
      call setsto(62*nd200,0,myext)
      call setsto(8*nd200,0,mynorbe)
      do lsi = 1,nirr
         ni=norbi(lsi)
         if (ni.ne.0) then
            do ll = 1, ni
               iii=iii+1
               do kk = nintp1,ninnex
                   if (oipsci().and..not.opcivect) goto 60
                   myext(iii,kk)=1
 60                continue
               enddo
            enddo
         endif
      enddo
      do iext=nintp1,ninnex
          do int=1,nintp1-1
             if(myext(int,iext).eq.1) then
        mynorbe(int,iro(iext))=mynorbe(int,iro(iext))+1
             endif
          enddo
      enddo
      return
      end
_ENDIF
