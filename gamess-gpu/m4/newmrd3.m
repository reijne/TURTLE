_IF(hpux_parisc)
c Reduce optimization due to high memory consumption on PA-RISC
c$HP$ OPTIMIZE LEVEL2
_ENDIF
c
c
c  $Author: hvd $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd3.m,v $
c  $State: Exp $
c
c ******************************************************
c ******************************************************
c             =   parkwa = newmrd3    =
c ******************************************************
c ******************************************************
c
      subroutine bambiu(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +                  sx,lsx,
     +                  pey,acoul,aexc,ideks,ndeks,
     +                  peyo,acoulo,aexco,nf31o,
     +                  sdelta,seng,lbuffs,
     +                  vect,wect,coff,nvect,
     +                  ot,nt444,
     +                  jkan,njkan,nrunt)
c
c diagonalises the extended reference space
c
c writes geyser0-raum onto ft99
c writes geyser1-raum onto ft99
c
c  nrunt : integer nruntype, if =1, then ft99 doesn't get printed
c  also no geyser1-raum calculated anymore
c
c  ---- ci-coefficients get printed !
c
c ---------------------------------------------------------
c
c    iclasl = counter for configuration classification integer
c
c    specification :
c
c    iclas(iclasl) = 1       : main
c    iclas(iclasl) = 2       : main (ev. geyser1-raum)
c    iclas(iclasl) = 3       : single excitation
c    iclas(iclasl) = 4       < trwa
c    iclas(iclasl) = 5
c    iclas(iclasl) = 6
c    iclas(iclasl) = 7
c
c ---------------------------------------------------------

      implicit REAL (a-h,o-z)
c
      character*1 wgl(50)
      REAL sch
      REAL zero0,rmilli
      REAL ci2,rbga,rbgb
      REAL crus,crud,crub
      integer ncci(36),nnn,nc1o,nbel
      integer nrunt
      integer isp
c
*     dimension a(5292),e(ndimf)
      REAL hp, a, w, hy, sx
      integer nt64, mxref2, ndimf, lsx
      dimension hp(nt64), a(nt64), hy(nt64)
      dimension w(mxref2)
c
      REAL e, f
      dimension e(ndimf), f(ndimf)
      dimension sx(lsx)
c
      REAL sdelta, seng
      integer lbuffs
      dimension sdelta(5,lbuffs),seng(lbuffs)
c
      REAL vect,wect,coff
      integer nvect
      dimension vect(nvect),wect(nvect),coff(nvect)
c
      REAL peyo
      REAL acoulo,aexco
      integer nf31o
      dimension peyo(nf31o),acoulo(nf31o),aexco(nf31o)
c
      REAL pey,acoul,aexc
      integer ndeks,ideks
      dimension ideks(ndeks)
      dimension pey(ndeks),acoul(ndeks),aexc(ndeks)
c
      integer ot
      dimension ot(nt444)
c
      integer jkan, njkan
      dimension jkan(njkan)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
INCLUDE(common/prints)
c
      parameter (sch = 0.0005d0)
      parameter (zero0 = 0.0d0)
      parameter (rmilli = 1.0d-6)
      parameter (crus = 1000.0d0)
      parameter (crud =  940.0d0)
      parameter (crub =  999.0d0)
      parameter (rmaxe= -999.0d4)
c
c --- commons
c
      REAL egey1,tras0,tdel0
      integer negey1,nko2,nstarv
      common /cegey1/egey1,tras0,tdel0,negey1,nko2,nstarv
c
      common /cczero/ zero
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
      integer nod
      common /cnod/ nod(ndk5)
      integer nform
      common /cffor/ nform

c
      integer mxretain
      parameter (mxretain = 50)
      REAL energret,facret
      integer itest,nopret,ipret,jkonret,nretain
      integer iretain
      common /cski1/itest(maxshl),
     + energret(mxretain),facret,nopret(mxretain),
     + ipret(mxretain),jkonret(mxretain*maxshl),nretain,
     + iretain
c
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
c
      integer ij, nplu
      common /scrtch/ ij(maxsym),nplu(ndk5)
c
c --- data to generate compatibility with the old ft31 format
c rh
c -----------------------------------------------------
      integer ihog
      common /cihog/ ihog(48)
      common /ccore/ core,isc(nn3)
      common /csac/  sac(nopmax+1)
c
      integer nytl, ndub
      common /cny /  nytl(ndk5),ndub(ndk5)
c     common /einf/  mhe(maxref),imain(maxref),iselct
INCLUDE(common/comrjb2)
      common /g/ trwe,trwf,trwg,trwh,trwi,trwj,imsec
      integer nsc, mn
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +nconf(ndk5),mconf(ndk5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
_IF(notused)
c
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
_ELSE
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,ixex,i21,i2ex,i31,i3ex,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
_ENDIF
c
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
      common /bx/   jsec,kmax
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
_IF(notused)
      parameter (nt44i = nt44 / 4)
      parameter (nt44r = nt44 / 8)
      logical*1 pt
      integer iw
      dimension iw(nt44i)
      dimension pt(nt44)
      equivalence (ptiw(1),iw(1),pt(1))
_ENDIF
c
      REAL skoef, savect, devect, doff
      integer ilifq
_IF(notused)
      REAL ptiw
      common/junk/ ptiw(nt44r),
_ELSE
      common/junk/ 
_ENDIF
     + skoef(5292),
     + savect(mxroot,42),devect(mxroot,126),
     + doff(maxref),ilifq(maxref)
c
      common/ftap/ntapess(10),mstvt,nf01,nf62,nhead,
     +            nf99,mtapev,nf11
c
      integer icald, maxci2, ipt00, isym
      logical debugs
      common /tap/ icald(9),maxci2,ipt00,isym(mxroot),
     +             debugs
c
      REAL pey0
      dimension pey0(10)
c
      integer imap
      REAL trsum,ew
      integer itym,iper
      common /cshu1/ ew(mxroot),trsum(10,mxroot),
     +               istm(nopmax+1),itym(mxroot),
     +               imap(504),iper(4)
c
      REAL ht, hs
      REAL ea, senk
      REAL cbeitr, dbeitr
      integer jkon
      common /junk2/ ea(mxroot),senk(500),cbeitr(mxroot),dbeitr(126),
     +               ht(9),hs(maxref),jkon(ndcon)
c
      integer isuper,inumer,ichek1,nrec,jtrk,ktrk
      integer mkon,ikan,iqr,kref
      common /junk3/ isuper(ndimh),inumer(ndimh),ichek1(maxshl),
     +               nrec(ndk5),jtrk(maxref),ktrk(maxref),
     +               mkon(ndimh),ikan(ndimh),iqr(ndimh),kref(ndimh)
c
********************************************************************
* table-ci : *
**************
      REAL t, af, ae
      common /a/ t(nlca)
      integer j9
      dimension j9(49),ae(ndtab)
      integer idra,idrc
      integer ndet,nsac
      integer jdrc,iaz,iez,jan,jbn,isc
      dimension ndet(ndk5)
      dimension nsac(ndk5)
      dimension idra(ndk5),idrc(ndk5)
      dimension jdrc(ndk5),jan(7)
      dimension jbn(7),iaz(5),iez(5)
      dimension af(8000)
      equivalence (af(1),t(1))
      equivalence (t(1),j9(1)),(j9(1),ndet(1)),(j9(6),nsac(1)),
     1(j9(11),iaz(1)),(j9(16),iez(1)),(j9(21),idra(1)),
     2(j9(26),idrc(1)),(j9(31),jdrc(1)),(j9(36),jan(1)),
     3(j9(43),jbn(1))
_IF(cray,ksr,i8)
      equivalence (af(50),ae(1))
_ELSE
      equivalence (af(26),ae(1))
_ENDIF
********************************************************************
      REAL rmin,sm,
     2 zero,core,butzi,corea,eigv,trash,
     3 tdel,trwa,trwb,trwc,trwd,trwe,trwf,trwg,trwh,trwi,trwj,dif
     4,care,sac,exc,tim,coul,
     5 dif1
*     REAL xint,yint,factor,pmt
c
      character *3 itag, char3i
      dimension itag(maxshl)
c
      character*1000 remstring
_IF(notused)
_ELSE
      dimension lout(510)
      equivalence (gout(1),lout(1))
      dimension ifront(10)
      data ifront / 1, 2, 6, 20, 70, 1, 1, 2, 5, 14 /
_ENDIF
c
      if (nrunt.eq.1) then
         do i=1,ndk5
            nplu(i) = nzero
         enddo
      end if
      write(iwr,*) '============================================'
      write(iwr,*) '= Solve zero-order problem and select safs = '
      write(iwr,*) '============================================'
      write(iwr,*)
c     write(iwr,*) 'runtype =',nrunt
c     write(iwr,*) 'built configurations to'
c     write(iwr,*) 'kfile   =',kfile
c     write(iwr,*) 'generated configurations from'
c     write(iwr,*) 'mtape   =',mtape
*     write(iwr,*) 'nston   =',nston
c     write(iwr,*)
*     write(iwr,*) 'nytl :',nytl(1),nytl(2),nytl(3),nytl(4),nytl(5)
*     write(iwr,*) 'nplu :',nplu(1),nplu(2),nplu(3),nplu(4),nplu(5)
*     write(iwr,*) 'nod  :',nod(1),nod(2),nod(3),nod(4),nod(5)
*     write(iwr,*) 'ndub :',ndub(1),ndub(2),ndub(3),ndub(4),ndub(5)
*     write(iwr,*) 'nconf:',nconf(1),nconf(2),nconf(3),nconf(4),nconf(5)
c     write(iwr,*) 'mconf:',mconf(1),mconf(2),mconf(3),mconf(4),mconf(5)
c --- preliminaries
      rmin =0.5d0*sqrt(2.0d0)
      egey1  = egey1 * rmilli
      ngy10  = nzero
      negey1 = nzero
      nko2   = nzero
      ncount = nzero
c --- setting the parameters
      trash = tras0
      tdel  = tdel0
c --- set max. no of configurations
c
      if(maxci2.le.0) maxci2 = maxci
c
c --- initialising the counters
clear      iclasl = nzero
      do i=1,50
        wgl(i)=' '
      enddo
      do loop=1,nnid
       olab(loop) = 0
      enddo
      ci2 = zero0
      nbel = nzero
      do i=1,36
         ncci(i) = nzero
      enddo
      do i=1,10
      istm(i) = nzero
        do ii=1,mxroot
          trsum(i,ii) = zero0
        enddo
      enddo
c
      nretain = 0
      nretg = 0
c
      jfile = 18
      jdeli = 19
      call rewftn(ntape)
      call rewftn(ideli)
      call rewftn(kfile)
      if (ical.gt.1) call rewftn(jdeli)
      call rewftn(linf)
      if (ical.lt.3) go to 156
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call rewftn(jfile)
c ---- jkan : mains read from parke4
      read (jfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     * ipt0,nplu,ndub,nod,nir,loc,ideks,jkan
      write (kfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     * ipt0,nplu,ndub,nod,nir,loc,ideks,jkan
      call rewftn(kfile)
c --- entry-point incase ical<3
 156  continue
      nrmap =     9
c     allow for max no. of configurations to veride value in newmrd_parinc
      imsec = max(maxci,maxci2)
      nad   =  4000
      nad1  = -3999
      call rewftn(mtape)
      ix    =   100
c ---- jkan : mains read from parke4
      read (kfile) ij,mconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     * ipt0,nplu,ndub,nod,nir,loc,ideks,jkan
*     if (nrunt.ne.nzero) then
*       lulu = nko2
*     end if
c change due to modified selection scheme
c restoring of jkan because the information is used by eigen but
c afterwards it is destroyed
c --- is probably not always at position 1000 !!!!
      do  neng=1,ndimh
        kref(neng)=jkan(neng)
      enddo
c end of changing
      if (ical.gt.1) ipt0=nzero
c --- reading the mains jkon !!!
      if (nrunt.eq.nzero) then
        read (mtape)
     . imo,m,nmul,ispace,nconf,nytl,mxex,jkon,nko,iswh
      else
        read (mtape)
        nko = nko2
      end if
      mult = nplu(1) + 1
_IF(notused)
      call dropenp(ntab,nad)
_ELSE
      call setbfc
c
c ... now check the contents of the table data base
c ... must exhibit controlled abort if non-valid first record (usually
c ... due to failure to assign correctly) rather than floating point error!
      call stopbk3
      is = 0
      call rdbak3(is)
      call stopbk3
      do i = 1,10
      if (lout(i).ne.ifront(i)) then
       write (iwr,6613) i, lout(i), ifront(i)
       call caserr('table data base has incorrect format')
      endif
      enddo
      if(.not.oprint(29)) write(iwr,6614)
c
_ENDIF
      is = jerk(mult)
      il = jbnk(mult)
      jz = jbun(mult)
      ifr1 = 1 - irsf
      ix = ifr1
_IF(notused)
      do i=1,il
        is = is +    1
        ix = ix + irsf
        call dreadp(af(ix),is)
      enddo
_ELSE
      ill=2
      call rdbak3(is)
      call stopbk3
_IF(cray)
      call fmove(lout,ndet,49)
_ELSE
      call icopy(49,lout,1,ndet,1)
_ENDIF
      is=is+nsz
      ik=1
 5551 if (il.lt.ill) go to 5550
      call rdbak3(is)
      is=is+nsz
      call stopbk3
      call dcopy(nword,gout(1),1,ae(ik),1)
      ik=ik+nword
      ill=ill+1
      go to 5551
c
 5550 continue
_ENDIF
      jsec = nzero
      do  i=1,jswh
        isc(i) = jsec
        jsec   = jsec + mconf(i)*nsac(i)
      enddo
      jsac = jsec
clear if(ipt0) 680, 600,681
c --- ipt0 <> 0
c --- extremly mysterious construction!
*680  do 780 i=1,lulu
*     do 781 j=1,lulu
*       if(mn(j).eq.i) goto 780
*781  continue
*780  jtrk(i)=j
c --- end of the extremly mysterious construction!
c --- replacement !
c --- ipt0<0
c680  continue
      if (ipt0.lt.nzero) then
       do  i=1,lulu
        j = mn(i)
        jtrk(j)=i
       enddo
c --- end replacement !
       ix = nzero
       do  i=1,lulu
        ktrk(i) = ix
        jj      = jtrk(i)
        ns      = nsc(jj)
        ix      = ix + nsac(ns)
       enddo
       ix = nzero
c
       do i=1,nrootx
        read(ird,546)(doff(j),j=1,jsec)
        corea = zero0
        do j=1,jsec
          butzi = doff(j)
          corea = corea + butzi*butzi
        enddo
        corea = 1.0d0/dsqrt(corea)
        do j=1,jsec
          doff(j) = corea*doff(j)
        enddo
        do j=1,lulu
          jj = mn(j)
          iz = ktrk(jj)
          ns = nsc(j)
          ns = nsac(ns)
          do k=1,ns
            iz = iz + 1
            ix = ix + 1
            coff(ix) = doff(iz)
          enddo
        enddo
       enddo
c
       write(iwr,286) (coff(i),i=1,ix)
clear go to 600
c
c ---- ipt0>0
      else if (ipt0.gt.nzero) then
clear 681  continue
c      read(ird,682) (isym(i),i=1,nrootx)
       write(iwr,*) 'zero-order vectors for root-homing:'
       write(iwr,682) (isym(i),i=1,nrootx)
       write(iwr,*)
      end if
c
c ---- ipt0=0
clear 600  continue
c --- perparation for the call to skina
      jto = nzero
      mm  = nzero
      kft = nzero
c
c
      if (nrunt.eq.1)then
         call rewftn(ntape)
*****    read(ntape)
      end if
c
*     write (6,*) 'read from ntape !!!'
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
      read (mtype) h
      if (ical.eq.2) read (kfile) iqr
      if (nform.gt.nzero) then
        read (nston,err=61010,end=61001) pey
        read (nston,err=61020,end=61001) acoul
        read (nston,err=61030,end=61001) aexc
        read (nston,err=61040,end=61001) core
c       write(iwr,*) 'pey:',(pey(ii1),ii1=1,10)
c       write(iwr,*) 'acoul',(acoul(ii1),ii1=1,10)
c       write(iwr,*) 'aexc:',(aexc(ii1),ii1=1,10)
c       write(iwr,*) 'core:',core
c       write(iwr,*) 'pey:'
c       write(iwr,'(5(1x,f12.8))') (pey(ii1),ii1=1,100)
c       write(iwr,*) 'acoul:'
c       write(iwr,'(5(1x,f12.8))') (acoul(ii1),ii1=1,100)
c       write(iwr,*) 'aexc:'
c       write(iwr,'(5(1x,f12.8))') (aexc(ii1),ii1=1,100)
      else
        read (nston,err=61050,end=61001) 
     +        peyo,acoulo,aexco,core
        do ii=1,ndeks
           pey(ii)   = peyo(ii)
           acoul(ii) = acoulo(ii)
           aexc(ii)  = aexco(ii)
        enddo
        write(iwr,*) 'pey:',(pey(ii1),ii1=1,10)
        write(iwr,*) 'pey:',(pey(ii1),ii1=1,10)
        write(iwr,*) 'acoul',(acoul(ii1),ii1=1,10)
        write(iwr,*) 'aexc:',(aexc(ii1),ii1=1,10)
        write(iwr,*) 'pey:'
        write(iwr,'(5(1x,f12.8))') (pey(ii1),ii1=1,100)
        write(iwr,*) 'acoul:'
        write(iwr,'(5(1x,f12.8))') (acoul(ii1),ii1=1,100)
        write(iwr,*) 'aexc:'
        write(iwr,'(5(1x,f12.8))') (aexc(ii1),ii1=1,100)
      end if
      call rewftn(nston)
csut zero should be passed on!!!!!
clear read (nston) ibm,zero
      read (nston,err=61060,end=61001)
      if(debugs) then
        write(iwr,*) 'ezero vor skina call ezero=',zero
        write(iwr,*) 'core  vor skina call core =',core
      endif
c
      call skina(pey,acoul,aexc,ideks,ndeks,
     +           f,ndimf,w,mxref2,ot,nt444,
     +           jkan,njkan)
c     maq = ideks(jsec+1)
c     write(iwr,861)(w(k),k=1,maq)
c new expression of the w matrix
      kend = nzero
      if(oprint(31)) then
      write(iwr,862)
      call writel(w,jsec)
c     write(iwr,*)
c      do 1101 iaus=1,jsec
c        kanf = kend + 1
c        kend = kend + iaus
c        write(iwr,861)(w(kaus),kaus=kanf,kend)
c1101  continue
      endif
c
      if(jsec.le.1) then
       hy(1) = 1.0d0
      else
       call eigenp(w,hy,jsec,0)
      endif
c
c ------------------------------------------------------
c ---- results of eigen stored on unit mstvt (77)
       if (nstarv.gt.nzero) then
*       write(iwr,*)
*       write(iwr,*) '================================'
*       write(iwr,*) '== root eigenvectors (ft77)   =='
*       write(iwr,*) '================================'
*       write(iwr,*)
        call rewftn(mstvt)
        write(mstvt) jsec
        write(mstvt) hy
        write(mstvt) w
      end if
c ------------------------------------------------------
      ib = nzero
      do i=1,jsec
         ib   = ib + i
         hs(i)= w(ib)
         ilifq(i)=(jsec-i)*jsec
      enddo
      ib = jsec+1
      km = jsec*jsec
      ico = nzero
      do 1155 i=1,jsec
        ib = ib - 1
        eigv=hs(ib) + zero
        eigvalr(i) = eigv
c        write(iwr,77)i,eigv
        lm = km - jsec
        lq = lm + 1
c        write(iwr,861)(hy(j),j=lq,km)
c --- code for ipt0 < 0 :
c     check zero order sec. eq. test
        if (ipt0.lt.0) then
         ilc = nzero
         ico = nzero
         do j=1,nrootx
           corea = 0.0d0
           do k=lq,km
             ilc  = ilc   + 1
             corea= corea + coff(ilc)*hy(k)
           enddo
           corea=dabs(corea)
           if (corea.gt.rmin) ico=j
           ht(j)=corea
         enddo
c
         write(iwr,1555) (ht(j),j=1,nrootx)
         jtrk(i) = ico
        end if
c end check zero order sec. eq. test
        km = lm
        ico = ico + maxref
1155  continue
c

      if(.not.oprint(29)) then
*      write(6,*) (mn(i),i=1,jsec)
       do  i=1,lulu
        j = mn(i)
        id(j)=i
       enddo
*      write(6,*)(id(i),i=1,lulu)
       ix = nzero
       do  i=1,lulu
        ktrk(i) = ix
        jj  = id(i)
        ns  = nsc(jj)
        ix  = ix + nsac(ns)
       enddo
*      write(6,*) (ktrk(i),i=1,lulu)
       ix = nzero
       isafs = 0
       do j=1,lulu
         jj = mn(j)
         iz = ktrk(jj)
         ns = nsc(j)
         ns = nsac(ns)
         do k=1,ns
           iz = iz + 1
           ix = ix + 1
           isafs = isafs + 1
           jtrk(iz) = ix
           id(iz) = jj
         enddo
       enddo
      if (isafs.gt.maxref) then
       write(iwr,4059) lulu, isafs, maxref
       call caserr('too many SAFS in reference space')
      endif
*      write(6,*) (jtrk(i),i=1,jsec)
       write(iwr,9000)
       call writey(hy,ilifq,coff,jtrk,id,jsec,jsec,nrootx)
       write(iwr,9003)
       write(iwr,9002)(eigvalr(k),k=1,jsec)
      endif
c
      irom = 4
      if (ical.lt.2) go to 2533
c --- code for ical>=2
      read (jdeli) trwa,trwb,trwc,trash,tdel,irom,wect
      iroma = irom + 1
      trash = trash+ (4-irom)*tdel
      ilc   = jsec*jsec
c
      do i=1,jsec
        ilc  = ilc - jsec
        w(i) = nzero
        ix   = nzero
        do j=1,nrootx
          corea = 0.0d0
          kmax  = ilc
          do k=1,jsec
            kmax = kmax + 1
            ix   = ix   + 1
            corea= corea+ hy(kmax)*wect(ix)
          enddo
          w(i) = w(i) + corea*corea
        enddo
      enddo
c
      write (iwr,861) (w(i),i=1,jsec)
      do i=1,nrootx
        corea = zero0
        do j=1,jsec
          if (w(j).gt.corea) then
             corea = w(j)
             k     = j
          end if
        enddo
        w(k) = -2.0d0
        isym(i) = k
      enddo
c
      go to 684
c
c --- ende code for ical>=2 ipt0 : root certification , =0 (normally)
2533  if (ipt0) 683,643,684
c----- preferably this should be changed !
c ---- ipt0 = - 1
 683  do 556 i=1,nrootx
        do j=1,jsec
          if (jtrk(j).eq.i) goto 556
        enddo
        go to 759
556   isym(i)=j
c
c --- ipt0 = +1
684   continue
      do i=1,nrootx
        ilc = 10000
        do j=1,nrootx
          if (isym(j).lt.ilc) then
            ilc = isym(j)
            ix  = j
          end if
        enddo
        itym(i)  = ilc
        isym(ix) = 10001
      enddo
c
      do i=1,nrootx
        isym(i) = itym(i)
      enddo
      write(iwr,559)(isym(i),i=1,nrootx)
      go to 645
c --- ipt0 = 0
 643  continue
      do i=1,nrootx
         isym(i) = i
      enddo
 645  continue
      ilc = nzero
c
      do i=1,nrootx
        ix   = isym(i)
        ix   = jsec - ix + 1
        ew(i)= hs(ix)
        lm   = (ix-1)*jsec
        do j=1,jsec
          ilc = ilc + 1
          lm  = lm  + 1
          vect(ilc) = hy(lm)
        enddo
      enddo
c
      dmax = -50000.0d0
      do i=1,nrootx
        if (ew(i).gt.dmax) dmax=ew(i)
      enddo
c
      dmin = 50000.0d0
      do i=1,nrootx
        if (ew(i).lt.dmin) dmin=ew(i)
      enddo
cccc --- what is this ?
c  --- extension of the testspace for r--mains : input !!!
cccc
      dmin = dmin - 0.2d0
      dmax = dmax + 0.2d0
      if (ical.eq.2) go to 576
      write(linf) jsec,nrootx,nytl,nplu,ndub,vect,ew,nconf
c    this data appears never to initialised, and is not used.
c    *idprog,idvers,iddate,idtime,idplus
      if (ical.eq.3) go to 290
      go to 575
 576  read (linf)
      go to 290
ccccc 575  read (ird,75) trash,tdel
 575  continue
c changing due to testing of matrixelements
      if(trash.lt.0.0d0) then
        trash  = -trash
        idruck =      1
        read(ird,85)ilong
        read(ird,86)(ichek1(ibein),ibein=1,ilong)
        write(iwr,*) 
     +  'following configuration is examined in particular ',
     +  (ichek1(ibein),ibein=1,ilong)
      endif
c
c --- setting up of the various thresholds
c
      write (iwr,9300) trash,tdel
c
      trash = trash*rmilli
      tdel  = tdel*rmilli
c
      trwa  = 0.05d0*trash
      trwb  = 0.25d0*trash
      trwc  = 0.50d0*trash
  290 continue
      trwd = trash
      trwe = trash + tdel
      trwf = trwe  + tdel
      trwg = trwf  + tdel
      trwh = trwg  + tdel
      trwi = trwh  + tdel
      trwj = trwi  + tdel
c --- writing out the values
      if (tdel.gt.0.0d0) then
      write (iwr,1600)trash,trwa,trwb,trwc,trwd,trwe,
     +                trwf,trwg,trwh,trwi,trwj
      else
      write (iwr,1600) trash
      endif
_IF(notused)
c change to perform modified selection sheme
*     if(iselct.eq.1) then
*       read(ird,15) factor,iatom,iwelch
*       write(iwr,16) factor,iatom,iwelch
*       read(ird,17)test
*       read(ird,17)test1
*       write(iwr,18)test,test1
*     endif
c reading and storing of the needed integrals from ninte=60
*     ninte = 60
*     rewind ninte
*     iwe1 = nzero
*     do ieng=1,iatom
*       iwe = iwelch(ieng)
*       iwe1 = iwe - iwe1 - 1
*       if(iwe1.ne.0) then
*         do iweg=1,iwe1
*           read(ninte)
*         end do
*       end if
*       read(ninte)seng
*       do iein=1,(ndeks-1)
*         sdelta(ieng,iein) = seng(iein)
*       end do
*       iwe1 = iwe
*     end do
*     write(iwr,9999) (sdelta(1,iein),iein=1,100)
c end of changings
_ENDIF
      do i=1,10
        istm(i) = nzero
        do j=1,nrootx
          trsum(i,j) = zero0
        enddo
      enddo
      do j=1,maxshl
        itest(j) = nzero
      enddo
      i = 4
      if(ical.lt.2) goto 119
      i = irom
      trash = trwd + (i-4)*tdel
  119 write (ideli) trwa,trwb,trwc,trwd,tdel,i,vect,ew
c     write(iwr,*) 'ideli written , ideli=',ideli
c
      call setbfc
c
      jpag = 1
      ipag = 1
      ifr  = nzero
c new runner for the isuper,inumer array
      iandr = nzero
c ---- runs over the sk
      do 79 i=1,iswh
        nc = nconf(i)
        kml= nsac(i)
        ik = i + 3
        go to (603,605,606,607,608),i
c --- sk=1
 603    jsl = idra(1)
        ix  = nad1
        kss = 0
        if (jz-2) 481,482,483
 483    js = jsl
        ks = idra(3) + idrc(3) + jdrc(3) - js
_IF(notused)
 485    continue
        do j=1,ks
          ix = ix + nad
          js = js + 1
          call dreadp (pt(ix),js)
        enddo
        call unpk14(pt,ix+nad-1,ot)
        go to 484
_ELSE
 485    if (kss.eq.ks) go to 484
        ix=ix+nad
        call rdbak3(js)
        js=js+nsz
        call stopbk3
        call unpack(lout,8,ot(ix),nad)
        kss=kss+1
        go to 485
_ENDIF
c
 482    ks = idra(3) + idrc(3) - jsl
        js = jsl
        go to 485
 481    ks = idra(2) + idrc(2) - jsl
        js = jsl
_IF(notused)
        do j=1,ks
          ix = ix + nad
          js = js + 1
          call dreadp(pt(ix),js)
        enddo
        call unpk14(pt,ix+nad-1,ot)
_ELSE
 486    if (kss.eq.ks) go to 6620
        ix=ix+nad
        call rdbak3(js)
        js=js+nsz
        call stopbk3
        call unpack(lout,8,ot(ix),nad)
        kss=kss+1
        go to 486
6620    continue
_ENDIF
        js = jan(2)
        ks = jbn(2)
c        kss => 0 to allow high spin cases (JvL, rwah 2003 )
        kss = 0
_IF(notused)
        do j=1,ks
          ix = ix + nad
          js = js + 1
          call dreadp(pt(ix),js)
        enddo
        call unpk14(pt,ix+nad-1,ot)
_ELSE
 487    if (kss.eq.ks) go to 6621
        ix=ix+nad
        call rdbak3(js)
        js=js+nsz
        call stopbk3
        call unpack(lout,8,ot(ix),nad)
        kss=kss+1
        go to 487
 6621   continue
_ENDIF
c
        js = idra(3)
        ks = idrc(3)
_IF(notused)
        do j=1,ks
          ix = ix + nad
          js = js +   1
          call dreadp(pt(ix),js)
        enddo
        call unpk14(pt,ix+nad-1,ot)
_ELSE
 488    if (kss.eq.ks) go to 484
        ix=ix+nad
        call rdbak3(js)
        js=js+nsz
        call stopbk3
        call unpack(lout,8,ot(ix),nad)
        kss=kss+1
        go to 488
_ENDIF
 484    if (nc.eq.0) go to 79
        mcl = mconf(i)
        ihp = nzero
_IF(notused)
        ic = iw(1)
        i2 = iw(2)
        i3 = iw(3)
_ELSE
        call upackx(ot)
        ic7=ic
        i27=i2
        i37=i3
_ENDIF
        if (mcl.ne.nzero) then
_IF(notused)
          jqn = idrc(1)*iwod + 1
          jrn = idrc(1)*nad
          ica = iw(jqn)   + jrn
          i2a = iw(jqn+1) + jrn
          i3a = iw(jqn+2) + jrn
_ELSE
          jqn = idrc(1)*nad + 1
          jrn = idrc(1)*nad
          call upackx(ot(jqn))
          ica=ic+jrn
          i2a=i2+jrn
          i3a=i3+jrn
_ENDIF
        end if
        j8 = 3
        mc = mconf(j8)
        if(debugs) write(iwr,*) 'mc=',mc
        if (mc.ne.nzero) then
          kmj  = nsac(j8)
          iej  = iez(j8)
          nis  = -kmj - kmj + (jan(2)-jsl)*nad
          nzer = kmj*kml
          if(debugs) write(iwr,*) 'nzer=',nzer,kmj,kml
          iix  = kmj + kmj + 1
        end if
        j7 = 2
        kc = mconf(j7)
        if (kc.eq.0) go to 403
        kmk  = nsac(j7)
        iek  = iez(j7)
        kzer = kmk*kml
        kr   = nod(j7)
        mps  = nplu(j7)
        mms  = 1
        jca  = kr + mult
        j3a  = jca*kmk + 1
        jab  = jan(1) - jsl
_IF(notused)
        jaq  = jab*iwod + 1
        jab  = jab*nad
        ic1  = iw(jaq) + jab
        i21  = iw(jaq+1) + jab
        i31  = iw(jaq+2) + jab
_ELSE
        jaq  = jab*nad + 1
        jab  = jab*nad
        call upackx(ot(jaq))
        ic1 = ic+jab
        i21 = i2+jab
        i31 = i3+jab
_ENDIF
        iiq  = kmk + kmk + 1
        md1  = ndub(j7)
        go to 403
c ---- sk=2
 605    if (nc.eq.0) go to 79
        mcl = mconf(i)
        ihp = idra(2) - jsl
_IF(notused)
        jqn = ihp*iwod + 1
        ihp = ihp*nad
        ic = iw(jqn) + ihp
        i2 = iw(jqn+1) + ihp
        i3 = iw(jqn+2) + ihp
_ELSE
        jqn = ihp*nad + 1
        ihp = ihp*nad
        call upackx(ot(jqn))
        ic7 = ic+ihp
        i27 = i2+ihp
        i37 = 13+ihp
_ENDIF
        if (mcl.ne.nzero) then
          jqn = idra(2) + idrc(2) - jsl
_IF(notused)
          jrn = jqn*nad
          jqn = jqn*iwod + 1
          ica = iw(jqn) + jrn
          i2a = iw(jqn+1) + jrn
          i3a = iw(jqn+2) + jrn
_ELSE
          jrn = jqn*nad
          jqn = jqn*nad + 1
          call upackx(ot(jqn))
          ica = ic + jrn
          i2a = i2 + jrn
          i3a = i3 + jrn
_ENDIF
        end if
        j7 = 3
        kc = mconf(j7)
        if (kc.ne.nzero) then
          kmk = nsac(j7)
          iek = iez(j7)
          kzer= kmk*kml
          iiq = kmk + kmk + 1
          md1 = ndub(j7)
          kr  = nod(j7)
          jca = kr + mult
          j3a = jca*kmk + 1
          mps = nplu(j7)
          mms = 2
          jab = jan(3) - jsl
_IF(notused)
          jaq = jab*iwod + 1
          jab = jab*nad
          ic1 = iw(jaq) + jab
          i21 = iw(jaq+1) + jab
          i31 = iw(jaq+2) + jab
_ELSE
          jaq = jab*nad + 1
          jab = jab*nad
          call upackx(ot(jaq))
          ic1 = ic + jab
          i21 = i2 + jab
          i31 = i3 + jab
_ENDIF
        end if
        j6 = 1
        lc = mconf(j6)
        if (lc.eq.0) go to 403
        jab = jan(1) - jsl
 615  continue
      ndl = ndet(j6)
      km2 = nsac(j6)
      la  = iaz(j6)
      lzer= kml*km2
      lr  = nod(j6)
      jcb = nod(i)  + mult
      j3b = jcb*kml + 1
      lps = nplu(j6)
      lms = j6 - 1
_IF(notused)
      jaq = jab*iwod + 1
      jab = jab*nad
      ic2 = iw(jaq)  + jab
      i22 = iw(jaq+1)+ jab
      i32 = iw(jaq+2)+ jab
_ELSE
      jaq = jab*nad + 1
      jab = jab*nad
      call upackx(ot(jaq))
      ic2 = ic + jab
      i22 = i2 + jab
      i32 = i3 + jab
_ENDIF
      go to 403
c
c --- sk= i = 3 :
c
 606  if (nc.eq.0) go to 79
      mcl = mconf(i)
      ihp = idra(3) - jsl
      if (jz.eq.1) ihp=ihp-jdrc(2)-jbn(3)
_IF(notused)
      jqn = ihp*iwod + 1
      ihp = ihp*nad
      ic  = iw(jqn)   + ihp
      i2  = iw(jqn+1) + ihp
      i3  = iw(jqn+2) + ihp
_ELSE
      jqn = ihp*nad + 1
      ihp = ihp*nad
      call upackx(ot(jqn))
      ic7 = ic + ihp
      i27 = i2 + ihp
      i37 = i3 + ihp
_ENDIF
      if (mcl.eq.0) go to 611
      jqn = idra(3) + idrc(3) - jsl
_IF(notused)
      jrn = jqn*nad
      jqn = jqn*iwod  + 1
      ica = iw(jqn)   + jrn
      i2a = iw(jqn+1) + jrn
      i3a = iw(jqn+2) + jrn
_ELSE
      jrn = jqn*nad
      jqn = jqn*nad  + 1
      call upackx(ot(jqn))
      ica = ic + jrn
      i2a = i2 + jrn
      i3a = i3 + jrn
_ENDIF
 611  j5  = 1
      lcl = mconf(j5)
      if(lcl.eq.0) go to 612
      ndi = ndet(j5)
      kmi = nsac(j5)
      jai = iaz(j5)
      nzei= kmi*kml
      nisi=-kml - kml + (jan(2)-jsl)*nad
      if (jz.eq.1) nisi=nisi-(jdrc(2)*nad)
 612  continue
      j6 = 2
      lc = mconf(j6)
      if (lc.eq.0) go to 403
      jab = jan(3)-jsl
      go to 615
c
c ---- sk= i = 4 :
c
  607 if (nc.eq.0) go to 79
      j5 = 2
      lcl= mconf(j5)
      jsl= jan(4)
      if (jz.eq.3) go to 489
      ks  = jbn(4)
      jqn = ks
      js  = jsl
      ix  = nad1
_IF(notused)
      do j=1,ks
        ix = ix + nad
        js = js+1
        call dreadp(pt(ix),js)
      enddo
      call unpk14(pt,ix+nad-1,ot)
_ELSE
      kss=0
 490  if (kss.eq.ks) go to 6222
      ix = ix+nad
      call rdbak3(js)
      js = js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss = kss+1
      go to 490
6222  continue
_ENDIF
      ks = idrc(4)
      js = idra(4)
_IF(notused)
      do j=1,ks
        ix = ix + nad
        js = js + 1
        call dreadp(pt(ix),js)
      enddo
      call unpk14(pt,ix+nad-1,ot)
      go to 492
_ELSE
      kss=0
 491  if (kss.eq.ks) go to 492
      ix = ix+nad
      call rdbak3(js)
      js = js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss = kss+1
      go to 491
_ENDIF
 489  ks = idra(4)+idrc(4)-jsl
      js = jsl
      ix = nad1
_IF(notused)
      do j=1,ks
        ix = ix + nad
        js=js+1
        call dreadp(pt(ix),js)
      enddo
      call unpk14(pt,ix+nad-1,ot)
      jqn=jbn(4)+jbn(5)
  492 ihp=jqn*nad
      jqn=jqn*iwod+1
      ic = iw(jqn)   + ihp
      i2 = iw(jqn+1) + ihp
      i3 = iw(jqn+2) + ihp
_ELSE
      kss=0
 613  if (kss.eq.ks) go to 6223
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 613
 6223 jqn=jbn(4)+jbn(5)
  492 ihp=jqn*nad
      jqn = jqn*nad+1
      call upackx(ot(jqn))
      ic7 = ic+ihp
      i27 = i2+ihp
      i37 = i3+ihp
_ENDIF
      if (lcl.eq.0) go to 614
      ndi = ndet(j5)
      kmi = nsac(j5)
      jai = iaz(j5)
      nzei= kmi*kml
      nisi=-kml-kml
 614  j6 = 3
      lc=mconf(j6)
      if (lc.eq.0) go to 403
      jab=jbn(4)
      go to 615
c
c ---- sk = i = 5 :
c
  608 if (nc.eq.0) go to 79
      jsl = jan(6)
      ks  = jbn(6)
      js  = jsl
      ix  = nad1
_IF(notused)
      do j=1,ks
        ix=ix+nad
        js=js+1
        call dreadp(pt(ix),js)
      enddo
      call unpk14(pt,ix+nad-1,ot)
      ksl=idra(5)
      ks=idrc(5)
      js=ksl
      do j=1,ks
        ix=ix+nad
        js=js+1
        call dreadp(pt(ix),js)
      enddo
      call unpk14(pt,ix+nad-1,ot)
_ELSE
      kss=0
 616  if (kss.eq.ks) go to 6224
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 616
c
 6224 ksl=idra(5)
      ks=idrc(5)
      js=ksl
      kss=0
  617 if (kss.eq.ks) go to 6225
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 617
 6225 continue
_ENDIF
      j5=3
      lcl=mconf(j5)
      if (lcl.eq.0) go to 618
      ndi = ndet(j5)
      kmi = nsac(j5)
      jai = iaz(j5)
      nzei= kmi*kml
      nisi=-kml-kml
 618  ihp=jbn(6)
_IF(notused)
      jqn=ihp*iwod+1
      ihp=ihp*nad
      ic=iw(jqn)+ihp
      i2=iw(jqn+1)+ihp
      i3=iw(jqn+2)+ihp
 403  ndt=ndet(i)
_ELSE
      jqn=ihp*nad+1
      ihp=ihp*nad
      call upackx(ot(jqn))
      ic7=ic+ihp
      i27=i2+ihp
      i37=i3+ihp
 403  ndt=ndet(i)
      ic=ic7
      i2=i27
      i3=i37
_ENDIF
      ia=iaz(i)
      ie=iez(i)
      mzer=kml*kml
      nd=ndub(i)
      nps = nplu(i)
      nms = i - 1
      niw = nps*nms
      nr  = nod(i)
      nnx  = nytl(i)
      n1  = nr + 1
      n2  = nr + 2
      n3  = nr - 1
      np1 = nps - 1
      np2 = np1 + np1
      nd1 = nd - 1
      nm1 = i - 2
      nm2 = nm1 + nm1
      j2  = kml + 1
      jc  = nr + nr
      j3  = jc*kml + 1
      iiz = kml + kml + 1
      ii  = iiz + kml
c --- counter for the jkan--field
      lpix = nzero
c --- read the first mkon--field for sk
      read (mtape) mkon
      jg=nzero
      if (ical.eq.2) go to 619
      np = kml*ndt
      ix = ia + ndt
      do l=1,np
        ix       = ix + 1
        a(l)     = ae(ix)
        skoef(l) = ae(ix)
      enddo
c     write(iwr,*) 'the coefficients '
c     write(iwr,*)(a(iaus),iaus=1,np)
c     write(iwr,*)(skoef(iaus),iaus=1,np)
      write (linf) ndt,kml,a
      ix = ie
      do l=1,mzer
        ix = ix + 1
        e(l) = ae(ix)
      enddo
      ix = ihp
      do 1666 l=1,kml
        ix = ix + 1
        ihog(l) = ot(ix+n12)
1666  continue
      if (i.eq.1) go to 624
      if (i.gt.2) go to 625
      do l=1,ndt
        imap(l)=l
      enddo
      go to 624
 625  np=ndt*nms
      do l=1,nms
        iper(l)=l
        imap(l)=l
      enddo
      ib=nms
 628  jb=nms
      jm=nr
 629  km=iper(jb)
      if (km.ne.jm) go to 630
      jb=jb-1
      if (jb.eq.0) go to 631
      jm=jm-1
      go to 629
 630  ip=iper(jb)+1
      iper(jb)=ip
      if (jb.eq.nms) go to 432
      jb=jb+1
      do l=jb,nms
        ip=ip+1
        iper(l)=ip
      enddo
 432  continue
      do l=1,nms
        ib  =ib+1
        imap(ib)=iper(l)
      enddo
      go to 628
 631  if (ib.ne.np) go to 755
 624  write (linf) ihog,imap,e
      go to 436
 619  read(linf)
      read (linf)
 436  continue
c**********************************************************************
c  ---- main loop
      do 573 j=1,nc
        do k=1,maxshl
         itag(k)=' '
        enddo
       kk=0
       do k=1,nnx
        jg = jg + 1
        itest(k) = mkon(jg)
        if(itest(k).gt.0)then
          kk=kk+1
          itag(kk)=char3i(itest(k))
        endif
        if(jg.ge.iwod) then
          jg = nzero
          read (mtape) mkon
        end if
       enddo
c
      idruk1 = idruck
      mm = mm + 1
      kk = olab(mm)
      if(idruck.ne.1) go to 81
      if(nnx.ne.ilong) go to 81
      do 1726 ites=1,nnx
       if(itest(ites).ne.ichek1(ites)) go to 81
1726  continue
      write(iwr,*)'itest',(itest(iaus),iaus=1,9),'kk,mm ',kk,mm
      idruk1 = -1
  81  continue
      if (mm.lt.nnid) go to 437
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 437  continue
*     if(idruk1.eq.(-1)) then
*       write(iwr,*)'mm,otherwise olab',mm,(olab(iaus),iaus=mm+1,mm+7)
*       write(iwr,*)'jto,h(jto) ',jto,(h(iaus),iaus=jto,jto+5)
*       idrk2 = 2
*     else
*       idrk2 = nzero
*     endif
c  kk=6 : --> main
c  kk=7 : --> single excitation
      if (kk-6) 438,439,640
c  for the mains
 439  corea=crud
      nko2 = nko2 + 1
cold      write (iwr,703) ipag,jpag,i,kml,nr,(itest(ll),ll=1,nnx)
crh   write (iwr,*) ipag,jpag,i,kml,nr,(itest(ll),ll=1,nnx)
cold  write (nf99,1871) nr,(itest(ll),ll=1,nnx)
      write (nf99,*) nr,'   ,'

      do krem=1,1000
         remstring(krem:krem)=' '
      enddo
      krem=1
      do irem=1,nnx
         write(remstring(krem:krem+4),'(i4)')itest(irem)
         krem=krem+5
      enddo
      remstring(krem:krem)='/'
      write(remstring(krem+2:krem+12),'(f10.1)')rmaxe
      write(nf99,'(a)')remstring(1:krem+12)
      call flushn(nf99)

crh   write (nf99,*) (itest(ll),ll=1,nnx),'    /     ', rmaxe
cold703  format('m',2i5,i2,2i3,2x,40i3)
      go to 1
c  for the single excitation
 640  corea = crus
      if(oprint(32)) then
        write (iwr,*) ipag,jpag,i,kml,nr,(itag(ll),
     .                             ll=1,nnx)
      endif
      go to 1
c  configuration from the normal mrd-ci - space
 438  ifl  = nzero
      jbab = nzero
      do 850 k=1,jswh
c --- # configurations in the reference space per sk
       md=mconf(k)
       if (md.eq.0) go to 850
c  --- jump to the dk-cases : kfc
       kfc = ik - k
       go to (851,852,853,854,855), kfc
c ---  dk > 2; dk < -2
       mz=nsac(k)
       do 1786 ix=1,md
         kbab = jbab
         jbab = jbab + mz
         ig7  = kbab - jsec
         do 1785 l=1,kml
          ig7 = ig7 + jsec
          in3 = ig7
          do 1784 ll=1,mz
            in3    = in3 + 1
            hp(in3)= zero0
1784      continue
1785    continue
1786  continue
c
      go to 850
c --- dk=0
 853  continue
      do 130 lb=1,md
        kbab=jbab
        jbab=jbab+kml
        if (ifl.eq.1) go to 857
        ifl = 1
        go to 131
 857    mm = mm + 1
        kk = olab(mm)
        if (mm.lt.nnid) go to 131
        mm = nzero
        read (ntape) olab8
        call unpack(olab8,8,olab,nnid)
 131    go to (858,132,133,134,135),kk
 135    mm = mm + 1
        kk = olab(mm)
        if(mm.lt.nnid) go to 136
        mm = nzero
        read (ntape) olab8
        call unpack(olab8,8,olab,nnid)
 136    mm = mm + 1
        ll = olab(mm)
        if(mm.lt.nnid) go to 137
        mm = nzero
        read (ntape) olab8
        call unpack(olab8,8,olab,nnid)
 137    if(kk.lt.ll) go to 138
        kk=ideks(kk)+ll
        go to 139
 138    kk=ideks(ll)+kk
 139    continue
        do l=1,mzer
          f(l) = zero0
        enddo
        tim=aexc(kk)
        jx=-kml
        do l=1,kml
          kx = ot(l+ihp+n12)+ia
          jx = jx+1
          lx = jx
          do m=1,kml
            lx=lx+kml
            kx=kx+ndt
            f(lx)=f(lx)+tim*ae(kx)
          enddo
        enddo
 146    continue
        ig7 = kbab - jsec
      mx=ie
      do 145 l=1,kml
      ig7=ig7+jsec
      in3=ig7
      ly=0
      my=mx
      do 145 m=1,kml
      tim=0.0d0
      mx=my
      do 143 iff=1,kml
      ly=ly+1
      mx=mx+1
 143  tim = tim + f(ly)*ae(mx)
      in3 = in3 + 1
 145  hp(in3)=tim
      go to 130
 144  mx=ie
      ig7=kbab-jsec
      do 253 l=1,kml
       ly=0
       ig7=ig7+1
       my=mx
       in3=ig7
       do 253 m=1,kml
        in3=in3+jsec
        tim=0.0d0
        mx=my
         do 1869 iff=1,kml
          ly = ly + 1
          mx = mx + 1
          tim= tim+ f(ly)*ae(mx)
1869     continue
 253  hp(in3)=tim
      go to 130
 858  ig7=kbab-jsec
      do 1880 mz=1,kml
       ig7 = ig7 + jsec
       in3 = ig7
       do 1879 nz=1,kml
        in3 = in3 + 1
        hp(in3) = nzero
1879   continue
1880  continue
      go to 130
 132  mm = mm + 1
      ll = olab(mm)
      if(mm.lt.nnid) go to 147
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 147  mm = mm + 1
      jj = olab(mm)
      if(mm.ge.nnid) then
        mm = nzero
        read (ntape) olab8
        call unpack(olab8,8,olab,nnid)
      end if
      mm = mm + 1
      kk = olab(mm)
      if(mm.lt.nnid) go to 149
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 149  if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob = nzero
      go to 151
 150  jj=ideks(kk)+jj
      ibob = 1
 151  jj = jj*ii + ica
      jto= jto + 1
      ip = ot(jj)
      if(ip.eq.255) go to 830
      coul=h(jto)
      if(jto.lt.iwod) go to 152
      jto=0
      read(mtype) h
152   jto = jto + 1
      exc = -h(jto)
      if(jto.lt.iwod) go to 153
      jto = nzero
      read(mtype) h
      go to 153
830   coul = -h(jto)
      if (jto.lt.iwod) go to 831
      jto = nzero
      read (mtype) h
831   jto = jto + 1
      exc = h(jto)
      if(jto.lt.iwod) go to 153
      jto = nzero
      read (mtype) h
 153  if (ll-2) 832,833,834
 832  bob(1) = coul + exc
      bob(2) = exc
      bob(3) = coul
      go to 835
 833  bob(1) = exc
      bob(2) =-coul
      bob(3) = coul + exc
      go to 835
 834  bob(1) = -exc
      bob(2) = -coul-exc
      bob(3) = coul
 835  continue
      do 1945 l=1,mzer
        f(l) = zero0
1945  continue
      mx = -kml
      do 155 l=1,kml
      jj = jj + 1
      mx = mx + 1
      ip = ot(jj)
      if(ip.eq.2)  go to 159
      jj = jj + 1
      kx = ot(jj)
      kx = kx + ia
      jj = jj  +1
      tim= bob(1)
      lx = mx
      do 1962 m=1,kml
        kx = kx + ndt
        lx = lx + kml
        f(lx) = f(lx) + tim*ae(kx)
1962  continue
      go to 155
 159  jj = jj + 1
      kx = ot(jj)
      kx = kx + ia
      tim= bob(2)
      lx = mx
      do 1973 m=1,kml
        kx   = kx + ndt
        lx   = lx + kml
        f(lx)= f(lx) + tim*ae(kx)
1973  continue
      jj=jj+1
      kx=ot(jj)
      kx=kx+ia
      tim=bob(3)
      lx=mx
      do 1983 m=1,kml
        lx=lx+kml
        kx=kx+ndt
        f(lx)=f(lx)+tim*ae(kx)
1983  continue
 155  continue
      if (ibob.eq.1) go to 144
      go to 146
 133  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 168
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 168  mm = mm + 1
      jj = olab(mm)
      if(mm.lt.nnid) go to 169
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 169  if(jj.gt.ll) go to 170
      jj   = ideks(ll) + jj
      ibob = nzero
      go to 171
 170  jj=ideks(jj)+ll
      ibob = 1
 171  jto  = jto + 1
      tim  = h(jto)
      if(jto.lt.iwod) go to 172
      jto = nzero
      read(mtype)h
 172  jj = jj*j2 + i2a
      do 2013 l=1,mzer
        f(l) = zero0
2013  continue
      jx =-kml
      ip = ot(jj)
      if(ip.eq.255) tim=-tim
      do 2028 l=1,kml
        jx = jx + 1
        jj = jj + 1
        kx = ot(jj)
        kx = kx + ia
        lx = jx
        do 2027 m=1,kml
          kx   = kx + ndt
          lx   = lx + kml
          f(lx)= f(lx) + tim*ae(kx)
2027    continue
2028  continue
      if (ibob.eq.1) go to 144
      go to 146
 134  mm = mm + 1
*     if(idruk1.eq.(-1)) then
*       write(iwr,*)'start dk=0 p=3'
*     endif
      ll = olab(mm)
      if(mm.lt.nnid) go to 179
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 179  mm = mm + 1
      jj = olab(mm)
      if(mm.lt.nnid) go to 180
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 180  mm   = mm + 1
      ntar = olab(mm)
      if(mm.lt.nnid) go to 302
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 302  mm = mm + 1
      nb = olab(mm)
      if(mm.lt.nnid) go to 303
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 303  mm = mm + 1
      mb = olab(mm)
*     if(idruk1.eq.(-1)) then
*     write(iwr,*)'ntar,nb,ideks(nb),ij(ntar)',ntar,nb,ideks(nb),
*    1ij(ntar)
*     endif
      if (mm.lt.nnid) go to 48
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  48  nb = ideks(nb) + mb + ij(ntar)
      sm = pey(nb)
      if(oprint(32)) then
        if(idruk1.eq.(-1)) then
          write(iwr,*)'ll,jj,ntar,nb,sm',ll,jj,ntar,nb,sm
        endif
      endif
      if(ll.gt.128) go to 181
      if(ll.lt.jj) go to 182
      ibob = nzero
      jj   = ideks(ll) + jj
      go to 205
 182  jj   = ideks(jj) + ll
      ibob = 1
 205  jj   = jj*j3 + i3a
      if(nr.eq.1) go to 185
      do 186 l=1,n3
       jto = jto + 1
       sm  = sm  + h(jto)
       if(jto.lt.iwod) go to 187
       jto = nzero
       read(mtype) h
187    jto = jto + 1
       sac(l)=-h(jto)
       if(jto.lt.iwod) go to 186
       jto = nzero
       read(mtype) h
 186  continue
 185  if(nd.eq.0) go to 188
      do 189 l=1,nd
      jto = jto + 1
      tim = h(jto)
      sm  = sm  + tim + tim
      if(jto.lt.iwod) go to 190
      jto = nzero
      read(mtype) h
 190  jto = jto + 1
      sm  = sm  - h(jto)
      if(jto.lt.iwod) go to 189
      jto = nzero
      read(mtype) h
 189  continue
 188  ip = ot(jj)
      if(ip.eq.1) go to 191
      sm = -sm
      if (nr.eq.1) go to 191
      do 2114 l=1,n3
        sac(l)=-sac(l)
2114  continue
 191  continue
      jj = jj + 1
      do 2119 l=1,mzer
        f(l) = zero0
2119  continue
      jx = -kml
      do 194 l=1,kml
        in = jj
        im = ot(in)
        jx = jx + 1
        in = in + 1
        kx = ot(in)
        kx = kx+ia
        lx = jx
        tim= sm
        if(im.eq.2) go to 195
        if(nps.eq.1) go to 196
        do 2136 m=1,np1
          in = in + 2
          nn = ot(in)
          tim =tim+ sac(nn)
2136    continue
 196  continue
      do 2142 m=1,kml
        lx   = lx   + kml
        kx   = kx   + ndt
        f(lx)= f(lx)+ tim*ae(kx)
2142  continue
      if(nms.eq.0) go to 194
      do 2157 m=1,nms
        in = in + 1
        kx = ot(in)
        kx = kx + ia
        lx = jx
        in = in + 1
        mx = ot(in)
        tim=sac(mx)
        do 2156 iff=1,kml
          lx   = lx    + kml
          kx   = kx    + ndt
          f(lx)= f(lx) + tim*ae(kx)
2156    continue
2157  continue
      go to 194
 195  if(nms.ne.1) then
        do 2164 m=1,nm1
          in = in + 2
          nn = ot(in)
          tim= tim+ sac(nn)
2164    continue
      end if
      do 2170 m=1,kml
        lx   = lx   + kml
        kx   = kx   + ndt
        f(lx)= f(lx)+ tim*ae(kx)
2170  continue
c
      do 2185 m=1,nps
        in = in + 1
        kx = ot(in)
        kx = kx + ia
        lx = jx
        in = in + 1
        mx = ot(in)
        tim= sac(mx)
        do 2184 iff=1,kml
          lx   = lx   + kml
          kx   = kx   + ndt
          f(lx)= f(lx)+ tim*ae(kx)
2184    continue
2185  continue
 194  jj = jj + jc
      if (ibob.eq.1) go to 144
      go to 146
 181  ll = 256 - ll
      if(ll.lt.jj) go to 204
      ibob = nzero
      jj=ideks(ll)+jj
      go to 206
 204  jj = ideks(jj)+ll
      ibob = 1
 206  jj   = jj*j3 + i3a
      jto  = jto + 1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 207
      jto = nzero
      read(mtype) h
 207  jto = jto + 1
      sm  = sm  + h(jto)
      if(jto.lt.iwod) go to 208
      jto = nzero
      read(mtype) h
 208  if(nr.eq.1) go to 209
      do 2221 l=1,n3
        jto = jto + 1
        sm  = sm  + h(jto)
        if(jto.ge.iwod) then
          jto = nzero
          read(mtype) h
        end if
        jto   = jto + 1
        sac(l)=-h(jto)
        if(jto.ge.iwod) then
          jto = nzero
          read(mtype) h
        end if
2221  continue
c
 209  if(nd.eq.1) go to 212
      do 213 l=1,nd1
        jto = jto + 1
        tim = h(jto)
        sm  = sm + tim + tim
        if(jto.lt.iwod) go to 214
        jto = nzero
      read(mtype) h
 214  jto = jto + 1
      sm  = sm - h(jto)
      if (jto.lt.iwod) go to 213
      jto = nzero
      read(mtype) h
 213  continue
 212  ip = ot(jj)
      if(ip.eq.1) go to 215
      sm = -sm
      if(nr.eq.1) go to 215
      do 2243 l=1,n3
        sac(l) = -sac(l)
2243  continue
c
 215  jj = jj + 1
      do 2248 l=1,mzer
        f(l) = zero0
2248  continue
      jx = -kml
      do 218 l=1,kml
       in = jj
       im = ot(in)
       jx = jx + 1
       in = in + 1
       kx = ot(in)
       kx = kx + ia
       lx = jx
       tim=-sm
       if (im.eq.2) go to 219
       if(nms.eq.0) go to 220
       in = in + np2
       jn = in
c
       do 2268 m=1,nms
        jn = jn + 2
        nn = ot(jn)
        tim= tim - sac(nn)
2268   continue
 220   continue
       do 2274 m=1,kml
        lx = lx + kml
        kx = kx + ndt
        f(lx) = f(lx) + tim*ae(kx)
2274   continue
       if(nms.eq.0) go to 218
       do 2289 m=1,nms
        in  = in + 1
        kx  = ot(in)
        kx  = kx + ia
        lx  = jx
        in  = in + 1
        mx  = ot(in)
        tim = sac(mx)
        do 2288 iff =1,kml
          lx = lx + kml
          kx = kx + ndt
          f(lx) = f(lx) + tim*ae(kx)
2288    continue
2289   continue
       go to 218
 219  in = in + nm2
      jn = in
c
      do 2298 m=1,nps
        jn = jn + 2
        nn = ot(jn)
        tim= tim - sac(nn)
2298  continue
      do 2303 m=1,kml
        lx = lx + kml
        kx = kx + ndt
        f(lx)= f(lx) + tim*ae(kx)
2303  continue
      do 2317 m=1,nps
        in = in + 1
        kx = ot(in)
        kx = kx + ia
        lx = jx
        in = in + 1
        mx = ot(in)
        tim= sac(mx)
        do 2316 iff=1,kml
          lx   = lx + kml
          kx   = kx + ndt
          f(lx)= f(lx) + tim*ae(kx)
2316    continue
2317  continue
c
 218  jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
 130  continue
      go to 850
c --- dk=+2
 851  continue
      if(debugs) write(iwr,*) 'call to skin2a'
      call skin2a(hp,nt64,f,ndimf,ot,nt444)
      go to 850
c --- dk=-2
 855  continue
      if(debugs) write(iwr,*) 'call to skin2b'
      call skin2b(hp,nt64,f,ndimf,ot,nt444)
      go to 850
c --- dk=+1
 852  continue
      if(debugs) write(iwr,*) 'call to skin1a'
      call skin1a(pey,ideks,ndeks,
     +            hp,nt64,f,ndimf,ot,nt444)
      go to 850
c --- dk=-1
 854  call skin1b(pey,ideks,ndeks,
     +            hp,nt64,f,ndimf,ot,nt444)
 850  continue
      do 2344 l=1,mzer
        f(l) = zero0
2344  continue
      care = core
      if (nr.eq.0) go to 106
      do 2353 l=1,nr
        mal  = itest(l)
        ntar = nir(mal)
        mal  = loc(mal) + 1
        mal  = ideks(mal) + ij(ntar)
        care = care + pey(mal)
2353  continue
      if(nd.eq.0) go to 108
 106  continue
      do 2364 l=n1,nnx
        mal  = itest(l)
        ntar = nir(mal)
        nal  = ideks(mal+1)
        mal  = loc(mal)+1
        mal  = ideks(mal) + ij(ntar)
        sm   = pey(mal)
        care = care + sm + sm + acoul(nal)
2364  continue
 108  if(nr.lt.2) go to 110
      do 2374 l=2,nr
        mal = itest(l)
        mal = ideks(mal)
        it = l - 1
        do 2373 m=1,it
          nal = mal + itest(m)
          care= care+ acoul(nal)
2373    continue
2374  continue
 110  continue
      if(nd.lt.2) go to 112
      sm = nzero
      do 2387 l=n2,nnx
       mal = itest(l)
       mal = ideks(mal)
       it  = l - 1
       do 2386 m=n1,it
        nal = mal+itest(m)
        tim = acoul(nal)
        sm  = sm + tim + tim - aexc(nal)
2386   continue
2387  continue
      care = care + sm + sm
 112  if(nr.eq.0.or.nd.eq.0) go to 114
      do 115 l=1,nr
      mal=itest(l)
      do 115 m=n1,nnx
      nal=itest(m)
      if(nal.lt.mal) go to 116
      nal=ideks(nal)+mal
      go to 117
 116  nal=ideks(mal)+nal
 117  sm=acoul(nal)
 115  care=care+sm+sm-aexc(nal)
 114  kk = ic
      kn = i2
      ml=i3
      jx=-kml
      do 118 l=1,kml
      jx=jx+1
      lx=jx
      kx=ot(l+ihp+n12)+ia
      tim=care
      if(nps.lt.2) go to 122
      mx=kn+1
      do 2424 iff=2,nps
        mx=mx+1
        it=iff-1
        mal=ot(mx+n12)
        mal=itest(mal)
        mal=ideks(mal)
        mz=kn
        do 2423 in=1,it
          mz=mz+1
          nal=ot(mz+n12)
          nal=itest(nal)+mal
          tim=tim-aexc(nal)
2423    continue
2424  continue
      kn = mx
      if(nms.lt.2) go to 122
      mx = ml + 1
c
      do 2442 iff=2,nms
        mx  = mx + 1
        mal = ot(mx+n12)
        mal = itest(mal)
        mal = ideks(mal)
        mz  = ml
        it  = iff - 1
        do in=1,it
          mz  = mz + 1
          nal = ot(mz+n12)
          nal = itest(nal) + mal
          tim = tim - aexc(nal)
        enddo
2442  continue
      ml = mx
 122  continue
      do iff=1,kml
        lx    = lx + kml
        kx    = kx + ndt
        f(lx) = f(lx) + tim*ae(kx)
      enddo
      if(niw.eq.0) go to 118
      do 2468 m=1,niw
        kk = kk + 1
        kx = ot(kk+n12)+ia
        kk = kk + 1
        mal= ot(kk+n12)
        kk = kk + 1
        nal= ot(kk+n12)
        if (k.eq.0) write (6,241) kx,mal,nal,kk
        mal = itest(mal)
        nal = ideks(mal) + itest(nal)
        tim =-aexc(nal)
        lx  = jx
        do iff=1,kml
          lx   = lx + kml
          kx   = kx + ndt
          f(lx)= f(lx) + tim*ae(kx)
        enddo
2468  continue
      if (k.eq.0) write (6,239) lx,kx,mal,nal,tim,f(lx),ae(kx)
 118  continue
      ly  = ie
      mxq = nzero
      mz  = nzero
      do k=1,nrootx
        ea(k) = zero0
      enddo
      do 670 k=1,kml
      eigv=0.0d0
      do 2489 l=1,kml
        ly = ly + 1
        mz = mz + 1
*       if(idruk1.eq.(-1)) then
*         write(iwr,*)'eigv,mz,f(mz),ly,ae(ly)',eigv,mz,f(mz),ly,ae(ly)
*       endif
        eigv=eigv+f(mz)*ae(ly)
2489  continue
      if (k.gt.1) go to 184
      if (ical.ne.2) go to 305
      kft=kft+1
      kx=iqr(kft)
      if (kft.lt.iwod) go to 462
      kft=0
      read (kfile) iqr
 462  if (kx.eq.12) go to 381
 305  if (eigv.lt.dmin) go to 183
      if (eigv.lt.dmax) go to 382
      ilse=0
      go to 184
 183  ilse=1
      go to 184
 382  if (ical.ne.2)  go to 465
      ilse=0
      write (iwr,466) ipag,jpag,i,kml,nr,eigv,(itag(ll),ll=1,nnx)
 184  ilc=0
      do 2539 l=1,nrootx
        ll=mxq
        dif = zero0
c
        do 2520 mw=1,jsec
          ll  = ll  + 1
          ilc = ilc + 1
          if(debugs) then
          if(idruk1.eq.(-1)) then
             write(iwr,*) 'k,l,mw,ll,hp(ll),vect(ilc),eigv'
             write(iwr,*) k,l,mw,ll,hp(ll),vect(ilc),eigv
          endif
          endif
          dif = dif+hp(ll)*vect(ilc)
2520    continue
c
c ---- calculation of the ci--coefficients
        savect(l,k)=dif/(eigv-ew(l))
        if (abs(savect(l,k)).gt.1.0d-30) then
         ci2 = sqrt(savect(l,k)*savect(l,k))
         ci2 = dabs(ci2)
         nnn = 1 - int(1.302883443d0*(dlog(ci2)))
       end if
       if (nnn.lt.1) then
           nnn = 1
        end if
       if (nnn.gt.36) then
           nbel=nbel+1
        else
           ncci(nnn)=ncci(nnn)+1
        endif
c --- epstein-nesbet guess !
        ea(l)=ea(l)+savect(l,k)*dif
2539  continue
c --- epstein-nesbet guess !
c656  ea(l)=ea(l)+(dif*dif)/(eigv-ew(l))
c
 670  mxq=ll
*     if(idruk1.eq.(-1)) then
*      write(iwr,*) 'ea(1)',ea(1)
*     endif
      dif = zero0
         if(ilse.eq.0) go to 880
         do 2551 l=1,nrootx
           if(ea(l).lt.dif) dif=ea(l)
2551     continue
        dif=-dif
        go to 882
 880  continue
      do 2557 l=1,nrootx
        if (ea(l).gt.dif) dif=ea(l)
2557  continue
 882  dif=dif/kml
c changing to perform modified selection sheme
c     write(iwr,950)(itest(kaus),kaus=1,10),dif
      dif1=dif
_IF(notused)
      if(iselct.ne.1) go to 885
      if(dif1.gt.trwd) go to 885
      if((dif1*factor).gt.trwd) then
c     write(iwr,951)(itest(kaus),kaus=1,10),dif
c test of the diagonal contribution of testconfiguration
c       write(iwr,*)'savect'
c       write(iwr,*) (savect(1,iaus),iaus=1,kml)
c transformation of the saf coef to det coef
        do 2574 kana=1,nrootx
          do 2573 kanada=1,ndt
            devect(kana,kanada)=0.0d0
2573      continue
2574    continue
c
        do 2590 kana=1,nrootx
          istart = nzero
          do 2589 itrans=1,ndt
            istart=istart+1
            ilauf=istart
            do 2588 iform=1,kml
       devect(kana,itrans)=devect(kana,itrans)+savect(kana,iform)*
     .skoef(ilauf)
c           write(iwr,*)'savect',iform,savect(kana,iform)
c           write(iwr,*)'skoef',ilauf,skoef(ilauf)
c           write(iwr,*)'devect',itrans,devect(kana,itrans)
            ilauf=ilauf+ndt
2588        continue
2589      continue
2590    continue
c devect contains coefficients of determinants
c       write(iwr,*) 'devect  '
c       write(iwr,*)(devect(1,iaus),iaus=1,ndt)
c calculation of value with assuming only alpha spins
        do 887 nlauf=1,iatom
        xint=0.0d0
c       write(iwr,*) xint,nr
          do 2602 ntest=1,nr
            jd=itest(ntest)
            int0=jd*(jd+1)/2
            xint=xint+sdelta(nlauf,int0)
2602      continue
          if(dabs(xint).lt.0.01)go to 887
          do 2606 ialex=1,ndt
            dbeitr(ialex)=xint
2606      continue
          ausdr=dbeitr(1)
c         write(iwr,*)'dbeitr',(dbeitr(iaus),iaus=1,ndt)
c taking beta spins into account
          nbeta = nr-nps
          imapl = nzero
c
          do 2624 iber=1,ndt
c           write(iwr,*)'jd,kd,int0,yint'
            do 2623 imar=1,nbeta
              imapl=imapl+1
              jd=imap(imapl)
              kd=itest(jd)
              int0=kd*(kd+1)/2
              yint=sdelta(nlauf,int0)
c             write(iwr,*) jd,kd,int0,yint
              dbeitr(iber)=dbeitr(iber)-yint-yint
2623        continue
2624      continue
c dbeitr contains values for single determinants
c         write(iwr,*) 'dbeitr',(dbeitr(iaus),iaus=1,ndt)
c transforming from detvalues to configurationvalue
          do 2631 inull=1,nrootx
            cbeitr(inull)=0.0d0
2631      continue
          do 2635 kana=1,nrootx
           do 2634 ialex=1,ndt
            cbeitr(kana)=cbeitr(kana)+devect(kana,ialex)*dbeitr(ialex)
2634       continue
2635      continue
c cbeitr contains configuration value
          do 2644 ialex=1,nrootx
            if(dabs(cbeitr(ialex)).gt.test(nlauf)) then
              dif1=dif*factor
c     write(iwr,954)(itest(kaus),kaus=1,5),dif,dif1,xint
              write(iwr,*)'ausdr,cbeitr',ausdr,cbeitr(1)
              go to 885
            endif
2644      continue
 887    continue
c test wether testkonfiguration is important single excitation
c with respect to mains
c     write(iwr,952)(itest(kaus),kaus=1,10),dif
c952  format(1x,'ausserdiagonal-test  itest,dif ',10i3,2x,f10.5)
        nstart=1
        do 1500 ksk=1,iswh
          mnum=mconf(ksk)
          if(mnum.eq.0) go to 1500
          mnr=nod(ksk)
          mnd=ndub(ksk)
          mnx=nytl(ksk)
          if(iabs(ksk-i).gt.1) then
c at least dk=2 --- no one-electron interaction possible
            nstart=nstart+mnum*mnx
            go to 1500
          endif
          do 1650 klauf=1,mnum
            if(ksk-i)1000,1100,1200
c testkonf. and main in same sk
1100        continue
            call check0(nr,nd,itest(1),kref(nstart),icheck)
            go to 1700
c testkonf. in higher sk
1000        continue
            call check1(mnr,mnd,nr,nd,kref(nstart),itest(1),icheck)
            go to 1700
c testkonf. in lower sk
1200        continue
            call check1(nr,nd,mnr,mnd,itest(1),kref(nstart),icheck)
1700        if(icheck.ne.0)then
            kende=nstart+7
c           write(iwr,*)nstart,(kref(kaus),kaus=nstart,kende)
c           write(iwr,*)icheck,sdelta(1,icheck)
              do 2686 nlauf=1,iatom
                xint=sdelta(nlauf,icheck)
                if(dabs(xint).gt.test1(nlauf))then
                  dif1=dif*factor
c     write(iwr,955)(itest(kaus),kaus=1,5),dif,dif1,icheck,xint
                  go to 885
                endif
2686          continue
            endif
            nstart=nstart+mnx
1650      continue
1500    continue
      endif
 885  continue
_ENDIF
c end of changings
      if (ical.ne.2) go to 463
      if (kx.eq.11) go to 414
      istm(kx)=istm(kx)+kml
      do 2699 l=1,nrootx
         trsum(kx,l)=trsum(kx,l)+ea(l)
2699  continue
      if (kx.lt.iroma) go to 416
      go to 414
 463  kft=kft+1
c ---- writing the geyser1--space
      if (nrunt.eq.nzero) then
       if (egey1.gt.-998.9d0) then
        if (ngy10.lt.1) then
         write(iwr,*) 'build geyser1 - space :'
         write(iwr,*) '======================='
         write(iwr,*)
*****    write(nf99,*) 'geyser1-space :          egey1=',(egey1/rmilli)
         ngy10 = 2
        end if
        if (dif1.gt.egey1.and.nr.le.5) then
         dener = dif1 - egey1
         write(iwr,*) 'dif1-egey1',dener
c        write(nf99,1870)
c    .     nr,(itest(ll),ll=1,maxshl),dener
         write(nf99,*) nr

      do krem=1,1000
         remstring(krem:krem)=' '
      enddo
      krem=1
      do irem=1,maxshl
         write(remstring(krem:krem+4),'(i4)')itest(irem)
         krem=krem+5
      enddo
      remstring(krem:krem)='/'
      write(remstring(krem+2:krem+12),'(f10.1)')dener
      write(nf99,'(a)')remstring(1:krem+12)
      call flushn(nf99)
c        write(nf99,*) (itest(ll),ll=1,maxshl),'    /      ',dener
         nko2 = nko2 + 1
         write(iwr,1701)
     .     ipag,jpag,i,kml,nr,eigv,(itag(ll),ll=1,nnx)
         negey1 = negey1 + 1
        write (6,*) negey1,' geyser1 - safs built'
        end if
       end if
      end if
c to perform modified selection sheme dif was changed in dif1
c if this option was not choose dif1=dif is valid
      if (dif1.lt.trwa) go to 404
      if (dif1.lt.trwb) go to 405
      if (dif1.lt.trwc) go to 406
      if (dif1.lt.trwd) go to 407
      if (dif1.lt.trwe) go to 408
      if (dif1.lt.trwf) go to 409
      if (dif1.lt.trwg) go to 410
      if (dif1.lt.trwh) go to 411
      if (dif1.lt.trwi) go to 412
      if (dif1.lt.trwj) go to 413
      go to 5
c
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwa c
ccccccccccccccc
c
  404 continue
      istm(1)=istm(1)+kml
      do 2755 l=1,nrootx
        trsum(1,l)=trsum(1,l)+ea(l)
2755  continue
      iqr(kft)=1
clear iclasl = iclasl + 1
clear iclas(iclasl) = 4
      if (kft.lt.iwod) go to 416
      kft = nzero
      write (kfile) iqr
      go to 416
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwb c
ccccccccccccccc
c
  405 continue
      istm(2)=istm(2)+kml
      do 2733 l=1,nrootx
        trsum(2,l)=trsum(2,l)+ea(l)
2733  continue
      iqr(kft)=2
      if (kft.lt.iwod) go to 416
      kft = nzero
      write (kfile) iqr
      go to 416
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwc c
ccccccccccccccc
c
  406 continue
      istm(3)=istm(3)+kml
      do 2789 l=1,nrootx
        trsum(3,l)=trsum(3,l)+ea(l)
2789  continue
      iqr(kft)=3
      if (kft.lt.iwod) go to 416
      kft = nzero
      write (kfile) iqr
      go to 416
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwd c
ccccccccccccccc
c
  407 continue
      istm(4)=istm(4)+kml
      do 2805 l=1,nrootx
        trsum(4,l)=trsum(4,l)+ea(l)
2805  continue
      iqr(kft)=4
      if (kft.lt.iwod) go to 416
      kft = nzero
      write (kfile) iqr
      go to 416
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwe c
ccccccccccccccc
c
  408 continue
      istm(5)=istm(5)+kml
      do 2821 l=1,nrootx
        trsum(5,l)=trsum(5,l)+ea(l)
2821  continue
      iqr(kft)=5
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwe c
ccccccccccccccc
c
  409 continue
      istm(6)=istm(6)+kml
      do 2837 l=1,nrootx
        trsum(6,l)=trsum(6,l)+ea(l)
2837  continue
      iqr(kft)=6
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwf c
ccccccccccccccc
c
  410 continue
      istm(7)=istm(7)+kml
      do 2853 l=1,nrootx
        trsum(7,l)=trsum(7,l)+ea(l)
2853  continue
      iqr(kft)=7
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwg c
ccccccccccccccc
c
  411 continue
      istm(8)=istm(8)+kml
      do 2869 l=1,nrootx
        trsum(8,l)=trsum(8,l)+ea(l)
2869  continue
      iqr(kft)=8
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwh c
ccccccccccccccc
c
  412 continue
      istm(9)=istm(9)+kml
      do 2885 l=1,nrootx
        trsum(9,l)=trsum(9,l)+ea(l)
2885  continue
      iqr(kft)=9
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
c
ccccccccccccccc
c dif1 < trwi c
ccccccccccccccc
c
  413 continue
      istm(10)=istm(10)+kml
      do 2901 l=1,nrootx
        trsum(10,l)=trsum(10,l)+ea(l)
2901  continue
      iqr(kft)=10
      if (kft.lt.iwod) go to 414
      kft=0
      write (kfile) iqr
      go to 414
c------------------------------------------------------------------------
   5  iqr(kft)=11
      if (kft.lt.iwod) go to 414
      kft = nzero
      write (kfile) iqr
c------------------------------------------------------------------------
 414  continue
      corea=dif
      if(oprint(32)) then
      if (dif.ge.sch) then
      write(iwr,401)ipag,jpag,i,kml,nr,(itag(loop),loop=1,nnx)
      write(iwr,4401)(ea(mo),mo=1,nrootx)
      write(iwr,*)mm,jto,'mm und jto'
      endif
      endif
      if (ical.eq.2) go to 2
      go to 700
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc -----> jump at the beginning       ccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 416  if(.not.oprint(32)) go to 573
      write(iwr,125)ndt,kml,i,j,nr,(itag(loop),loop=1,nnx)
      write(iwr,4401)(ea(mo),mo=1,nrootx)
      go to 573
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   r-mains :    c
cccccccccccccccccccc
 465  kft=kft+1
      iqr(kft)=12
      if (kft.lt.iwod) go to 381
      kft = nzero
      write (kfile) iqr
 381  corea = crub
c rwah
      nxold=0
      if (nnx.gt.40) then
         write(iwr,'(a,i4,a)')'nnx greater than 40 ',nnx,
     *   ' reduced to 40'
         nxold=nnx
         nnx=40
      endif
      write(iwr,701) ipag,jpag,i,kml,nr,(itag(ll),ll=1,nnx)
      write(iwr,7701)eigv
c
      if (ifbuen) then
c
c     retain this configuration as a main configuration ?
c
       if (eigv-eigvalr(nrootx).le.facret.and.
     +     nretain.lt.mxretain.and.iretain.le.nbuenk) then
        nretain = nretain + 1
        energret(nretain) = eigv
        nopret(nretain) = nr
        ipret(nretain) = nytl(i)
        do loop = 1, ipret(nretain)
         jkonret(nretg+loop) = itest(loop)
        enddo
        write(iwr,12788)nretain, energret(nretain)
        write(iwr,12789)nopret(nretain)
        write(iwr,12790)(jkonret(nretg+loop),loop=1,ipret(nretain))
        write(iwr,12791)
        nretg = nretg + maxshl
       endif
      endif
c
      if (nxold.ne.0) then
         nnx=nxold
         nxold=0
      endif
c rwah
      if (ical.eq.2) go to 2
      go to 700
ccc -----> jump at the end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- it was a main !  ccc
ccccccccccccccccccccccccccccc
  1   kft = kft+1
      if (ical.eq.2) go to 2
      iqr(kft) = 11
      if (kft.lt.iwod) go to 700
      kft = nzero
      write (kfile) iqr
      go to 700
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc -----> jump at the end          --------------------------------\
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   |
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   v
  2   if (kft.lt.iwod) go to 686
      kft = nzero
      read (kfile) iqr
 686  lpix = lpix + nnx
      jpag = jpag + 1
      ipag = ipag + kml
      if (lpix.lt.iwod) go to 573
      lpix = lpix - iwod
      read (linf)
      go to 573
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc --- entry 700                                                     c
cccc     end                                                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
700   continue
      do 2981 ll=1,nnx
       lpix = lpix + 1
       jkan(lpix) = itest(ll)
       if (lpix.ge.iwod) then
        lpix = nzero
        write (linf)jkan
       end if
2981  continue
      jpag = jpag + 1
      ipag = ipag + kml
      do 705 ll=1,kml
       ifr=ifr+1
       senk(ifr)=corea
       if(ifr.lt.irsf) go to 705
       ifr = 0
       write(ideli) senk
       if(debugs) write(iwr,*) 'senk written!'
705   continue
c storing of sk and number within sk of selected konf
c will be used for wald
      iandr=iandr+1
      isuper(iandr)=i
      inumer(iandr)=j
      if(iandr.lt.iwod) go to 573
        iandr = nzero
        write(nf62) isuper,inumer
c ---- 573 : entry point for small thresholds and also
c            the end of the do--loops over #nconf/sk !
 573  continue
c end of main loop
c**********************************************************************
      if (ical.ne.2) go to 469
       read (linf)
      go to 79
 469  lpix=lpix+1
      jkan(lpix)=0
      write (linf) jkan
      if(debugs) then
       write(iwr,*)' auf file',linf
       write(iwr,*)(jkan(iii),iii=1,100)
      endif
  79  continue
c storing of sk and number within sk of selcted konf
c will be used for wald
      iandr=iandr+1
      isuper(iandr)=0
      inumer(iandr)=0
      write(nf62) isuper,inumer
      iqr(kft+1)=0
      if (ical.ne.2) write (kfile) iqr
      immm=ipag-1
c
      if(nretain.gt.0) write(iwr,12792) nretain
c
      write(iwr,*)
      write(iwr,*)
     .   'summary of energy lowerings'
      write(iwr,*)
     .   '==========================='
      write(iwr,*)
c
      pey0(1)=trwa*crus
      pey0(2)=trwb*crus
      pey0(3)=trwc*crus
      pey0(4)=trwd*crus
      pey0(5)=trwe*crus
      pey0(6)=trwf*crus
      pey0(7)=trwg*crus
      pey0(8)=trwh*crus
      pey0(9)=trwi*crus
      pey0(10)=trwj*crus
      write(iwr,430)
      do 440 i=1,irom
        immm=immm+istm(i)
 440  continue
      istm(1)=immm-istm(1)
      do 3048 i=2,10
       istm(i)=istm(i-1)-istm(i)
3048  continue
      do 442 i=1,nrootx
         pey(1)=trsum(1,i)*crus
         do 3055 j=2,10
           j1=j-1
           trsum(j,i)=trsum(j1,i)+trsum(j,i)
           pey(j)=trsum(j,i)*crus
3055     continue
c --- printing the results :
      write(iwr,'(2x,f14.9,2x,i12,2x,f12.5,11x,i2)')
     +            pey0(1),istm(1),pey(1),i
      if (tdel.gt.0.0d0) then
       do ii=2,10
         write(iwr,'(2x,f14.9,2x,i12,2x,f12.5,11x,i2)')
     +             pey0(ii),istm(ii),pey(ii),i
       enddo
      endif
      write(iwr,*)
      write(iwr,*)'from a total of ',immm,' safs'
      write(iwr,*)
c --- end of printing
  442 continue
c --- printing the guessed ci - coefficients
      nc1o   = nzero
      ncmax  = 1
      rbga = 1.0d0
      rbgb = zero0
      do 3074 ii=1,36
         nc1o   = nc1o + ncci(ii)
3074  continue
      write(iwr,*)
      write(iwr,*)
     . 'total number of guessed ci-coefficients:',nc1o
      write(iwr,*)
      write(iwr,*)
     . 'ci-coefficients below the guess :       ',nbel
c     write(iwr,*)
c    .'##############################################################'
c234567890123456789312345678941234567895123456789612345678971234567890
      write(iwr,*)
      write(iwr,*) 'guess of the ci - coefficients:'
      write(iwr,*) '*******************************'
      write(iwr,*)
      do 3091 i=1,36
        ntest = ncci(i)
        ncmax = max(ncmax,ntest)
3091  continue
c
      write(iwr,*) 'coefficients domain'
c     first decide on last finite coefficient domain
      l136 = 0
      do i=1,36
      if(ncci(i).gt.0) then
       l136 = i
      endif
      enddo
c
      do 3106 i=1,l136
        ncount = 50*float(ncci(i))/float(ncmax)
        do 3097 j=1,ncount
          wgl(j)='='
3097    continue
        rbgb = (10.0d0)**((1-i)*0.333333333d0)
        rbgb = (0.1d0)*rbgb
        write(iwr,999) rbga,'-',rbgb,' ',wgl,ncci(i)
        do 3102 j=1,50
          wgl(j)=' '
3102    continue
        rbga = rbgb
3106  continue
      write(iwr,*)
      write(iwr,*)
     .'##############################################################'
c234567890123456789312345678941234567895123456789612345678971234567890
      write(iwr,*)
c --- negey1 at the end of ft99
      if (nrunt.eq.nzero) then
       write(nf99,*) negey1,',          = negey1'
       write(nf99,*) 
     +      '=================================================='
       call flushn(nf99)
      end if
c
c--------------------------------------------------------------------
c --- rausschreiben der konfigurationsklassifizierung auf ft77
c --- writing the coefficients classification on ft77
c
clear write(77) iclas
c--------------------------------------------------------------------
      if (ical.eq.2) then
       irom=istm(irom)
       irom=(irom-1)/irsf+1
       do i=1,irom
         read (jdeli) senk
         write (ideli) senk
       enddo
       go to 428
      endif
c
      if (ifr.ne.0) then
       write(ideli) senk
      endif
      if (istm(4).gt.imsec) then
c
        call shuffl(a,nt64,e,ndimf,sx,lsx,vect,nvect,jkan,njkan)
        return
c
      endif
428   write (ideli) nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm
c     it appears that invoking drclos on the Origin causes
c     a segmentation violation buried in the system I/O routines ..
c     call drclos
      return
 759  write(iwr,765) i,nrmap
      call caserr(
     + 'desired root not found on basis of overlap criterion')
 755  write (iwr,756) ib,ndt,nms
      call caserr('error in imap generation')
61010 write (iwr,*)' ***** reading pey ****'
      go to 61000
61020 write (iwr,*)' ***** reading coul ****'
      go to 61000
61030 write (iwr,*)' ***** reading aexc ****'
      go to 61000
61040 write (iwr,*)' ***** reading core ****'
      go to 61000
61050 write (iwr,*)' ***** reading pey0 ****'
      go to 61000
61060 write (iwr,*)' ***** null read ****'
61000 write(iwr,61005)
      call caserr(
     + 'error in reading transformed integral interface')
61001 write(iwr,61006)
      call caserr(
     + 'unexpected end of transformed integral interface')
c --- formats
_IF(notused)
  15  format(f10.2,6i3)
  16  format(1x,'factor,iatom,iwelch ',f10.2,1x,6i3)
  17  format(5f8.1)
  18  format(1x,'test-treshold for single centers '/,
     *1x,'diagonal contribution      ',5f8.1/,
     *1x,'outerdiagonal contribution ',5f8.1/)
9999  format(10(1x,f11.2))
_ENDIF
12790 format(5x,'MOs: ',256i4)
12789 format(5x,'number of open shells = ',i2)
12788 format(5x,'Retained configuration ',i2,' energy = ', f12.5)
12791 format(/)
12792 format(/1x,'**** Number of Retained Configurations = ',i2)
61005 format(1x,'error in reading transformed integral interface')
61006 format(1x,'end of transformed integral interface')
 756  format(10x,'error in imap generation',3i6)
 765  format(/10x, 
     +'desired root not found on basis of overlap criterion',2i6)
c 75  format(6f12.8)
 999  format(1x,d12.6,1x,a,1x,d12.6,1x,a,50a,1x,i10)
6613  format(/1x,'*** table data base in error: word ',i3,' = ',i4,
     +' and should be ',i4/
     +' *** check assignment of table data set'//)
6614  format(/1x,'*** table data base assigned correctly')
 862  format(/20x,'root hamiltonian matrix'/20x,23('-')/)
4059  format(
     + 5x,'*****************************************'/
     + 5x,'TOO MANY SAFS in REFERENCE CONFIGURATIONS'/
     + 5x,'no. of reference functions = ', i4/
     + 5x,'no. of SAFS                = ', i4/
     + 5x,'no. of allowed SAFS        = ', i4/
     + 5x,'****************************************'/)
9000  format(/40x,50('=')/
     + 40x,'=            Root Eigenvectors                   ='/
     + 40x,'= (coefficients in input configuration ordering) ='/
     + 40x,50('='))
9003  format(40x,17('-')/
     +       40x,'root eigen values'/40x,17('-'))
9002  format(/10x,8f14.7)
 559  format(/1x,
     *'*** selection based on following roots of zero order problem:'/
     *10x,20i3)
9300  format(/' *** threshold specified ***'//
     *' minimal selection threshold ',f7.2,' microhartree'/
     *' threshold increment for use in selection ',f7.2,
     *' microhartree'/)
c1871 format(24i3,'     ','m')
1701  format('g',2i6,i2,2i3,2x,f12.8,2x,40(1x,a3))
 401  format('s',2i6,i2,2i3,2x,40(1x,a3))
4401  format(2x,4f12.8/2x,5f12.8)
 125  format(i4,i3,i2,i6,i2,5x,40(1x,a3))
 466  format(//'z',2i6,i2,2i3,2x,f12.8,2x,40(1x,a3)//)
 701  format('r',2i6,i2,2i3,5x,40(1x,a3))
7701  format(43x,f18.8)
c 77  format(20x,i5,f20.8)
  85  format(i3)
  86  format(80i2)
 239  format(2x,4i6,3f20.8)
 241  format(2x,21i6)
 286  format(2x,14f8.5)
 430  format(' threshold domain (mh)',4x,'safs',5x,
     .       'energy sum        root')
 546  format(14f5.2)
 682  format(14i5)
c702 format('t',2i5,i2,2i3,2x,20i3)
 861  format(2x,10f14.6)
c950  format(1x,'conf to be tested  itest,dif ',10i3,2x,f10.5)
c951  format(1x,'diagonal-test  itest,dif ',10i3,2x,f10.5)
c954  format(1x,'diagonal gross  itest,dif,dif1,xint',
c    *5i3,2(2x,f10.7),2x,f10.2)
c955  format(1x,'off-diagonal gross  itest,dif,dif1,icheck,xint',
c    *5i3,2(2x,f10.5),(2x,i3),(2x,f10.2))
 1555 format(/
     *' overlap between zero order eigenvector no.',i3,
     *' and input eigenvectors'/
     */6(/i3,5x,f8.4))
1600  format(/1x,'following threshold classes to be investigated'/
     +        1x,11f11.6/)
      end
      subroutine writey(w,ilifq,xcoff,iord,iroot,newbas,nnn,nroot)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
INCLUDE(common/prints)
      dimension w(*),ilifq(*),iord(*),iroot(*),xcoff(*)
100   format(/1x,'SAF',1x,'MAIN'/5x,8i14//)
101   format(/)
200   format(/1x,'SAF',1x,'MAIN'/6x,12i9//)
102   format(1x,i3,1x,i3,1x,8(f13.7,a1))
202   format(1x,i3,1x,i3,1x,12(f8.4,a1))
c
c     flag large coefficients for tagging
c
      do i = 1, newbas
       ij = ilifq(i)
       do j = 1, newbas
        xcoff(j+ij) = ' '
       enddo
       loop=idamax(newbas,w(ij+1),1)
       wmax = abs(w(ij+loop)) / 4.0d0
       do j = 1, newbas
        if(abs(w(ij+j)).ge.wmax) then
        xcoff(ij+j) = '*'
        endif
       enddo
      enddo
c
      lprnt = nnn
      if(.not.oprint(20)) lprnt=min(nroot+5,nnn)
      m=1
c
      m5=12
      if(oprint(20))m5=8
      n=m5
2106  if(lprnt.lt.m)return
      if(n.gt.lprnt)n=lprnt
      if(oprint(20))then
       write(iwr,100)(i,i=m,n)
      else
       write(iwr,200)(i,i=m,n)
      endif
      do 2111 j=1,newbas
      imain = iroot(j)
      if(oprint(20)) then
       write(iwr,102)j,imain,(w(iord(j)+ilifq(i)),
     +                    xcoff(iord(j)+ilifq(i)),i=m,n)
      else
       write(iwr,202)j,imain,(w(iord(j)+ilifq(i)),
     +                    xcoff(iord(j)+ilifq(i)),i=m,n)
      endif
2111  continue
      m=m+m5
      n=n+m5
      write(iwr,101)
      goto 2106
      end
      subroutine bearb(iot,niott,ideks,ndeks,jcon,nd8)
c
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndeks,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
c
      integer lj
      common /clj2/ lj(8)
c
      integer kc,kd,lab
      common /ckd/ kc(ndimh),kd(ndimh),lab(3)
      common /cnit/ nit(nitmax)
      common /a/ t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     +ndx,nston,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     +igmax,nr2,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +mconf(ndk5),nconf(ndk5),mh,jswh,
     +id(maxref),ie(maxref),nsc(maxref),
     +imo,mc,ihz,nd,nr,nnx,nt,ny,mt,mz,nzx,iht,nz,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr1,m,nad,nad1,ifr1,mx,md,mr,ipag
c
      l3=mz
      l4=mt
      do 105 k=1,mc
      if (nd.eq.0) go to 106
      jz=l3
      do 107 l=1,nd
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 108,109,107
 107  continue
      go to 106
 108  ip=l
      if (l.eq.nd) go to 110
      ia=l+1
      do 111 l=ia,nd
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj).ne.2) go to 112
  111 continue
      go to 110
  112 mm=mm+1
      olab(mm)=1
      if (mm.lt.nnid) go to 113
      mm=nzero
c
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      goto 113
 109  ip=l
      if (l.eq.nd) go to 114
      ia=l+1
      do 115 l=ia,nd
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 112,116,115
115   continue
      go to 114
116   if (l.eq.nd) go to 117
      ia=l+1
      do 118 l=ia,nd
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2)  go to 112
 118  continue
      go to 117
 110  if (nr.eq.0) go to 119
      jz=l4
      do 120 l=1,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,120,121
 120  continue
      go to 119
 121  ip=l
      if (l.eq.nr) go to 122
      ia=l+1
      do 123 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.1) go to 112
 123  continue
      go to 122
 119  mm=mm+1
      olab(mm)=5
      if (mm.lt.nnid) go to 124
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 124  jz=nz
      if (ip.eq.1) go to 125
      kz=l3
      ja=ip-1
      do 126 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 126  continue
 125  if (ip.eq.nd) go to 128
      ja=ip+1
      kz=l3+ip
      do 129 l=ja,nd
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 127
 129  continue
 128  kk=iot(jz+1)
  127 mm=mm+1
      olab(mm)=ll
      if (mm.lt.nnid) go to 130
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  130 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nnid) go to 113
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 113
  122 mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 131
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  131 jz=nt
      if (ip.eq.1) go to 132
      kz=l4
      ja=ip-1
      do 133 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 133  continue
 132  if (ip.eq.nr) go to 650
      ja=ip+1
      kz=l4+ip
      do 651 l=ja,nr
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 134
 651  continue
 650  jz=jz+1
      kk=iot(jz)
  134 mm=mm+1
      olab(mm)=jz-nt
      if (mm.lt.nnid) go to 135
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  135 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nnid) go to 136
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  136 nt1r=nir(ll)
      nl1r=loc(ll)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(kk)
   20 if (nt1r-nb1r) 23,22,99
   99 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq=icq+nm2r
      icq=icq+nm1r
   33 if (icq.lt.idq) go to 34
      icq = ideks(icq)+idq
      go to 29
   34 icq = ideks(idq)+icq
      go to 29
   23 iax=ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      go to 33
   22 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 46
      iay=ideks(nl1r)+nm1r
      go to 47
   46 iay= ideks(nm1r)+nl1r
   47 if (nl1r.lt.nm2r) go to 48
      iby=ideks(nl1r)+nm2r
      go to 50
   48 iby=ideks(nm2r)+nl1r
   50 icx=ideks(iax+1)
      if (iay.lt.iby) go to 53
      icq=ideks(iay)+iby
      go to 29
   53 icq=ideks(iby)+iay
 29   icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 113
cdel 55   format(3x,21i6)
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 113
 117  jz=l4
      do 137 l=1,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,137,138
  137 continue
  138 l1=l
      ia=l+1
      do 139 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,139,140
  139 continue
  140 l2=l
      if (l.eq.nr) go to 141
      ia=l+1
      do 142 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
 142  continue
 141  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 143
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  143 mm=mm+1
      olab(mm)=1
      if (mm.lt.nnid) go to 144
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 144  jz=nt
      do 145 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 146
 145  continue
 146  ip=l
      ia=l+1
      do 870 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.nn) go to 147
 870  continue
 147  mm=mm+1
      olab(mm)=ip+ideks(l-1)
      if (mm.lt.nnid) go to 148
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  148 mm=mm+1
      olab(mm)=ideks(l2-1)+l1
 540  if (mm.lt.nnid) go to 149
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 149  nt1r=nir(ll)
      nt2r=nir(nn)
      nl1r=loc(ll)
      nl2r=loc(nn)
 189  nb1r=nir(kk)
      nb2r=nir(jj)
      nm1r=loc(kk)
      nm2r=loc(jj)
      kix=3
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 224
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 225,226,227
  227 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  240 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 229
  226 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  241 idq = (nl2r-1)*ljn+nm2r
      go to 233
  225 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  235 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 229
  224 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 230,231,232
  232 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 229
  231 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  233 if (icq.lt.idq) go to 234
      icq = ideks(icq)+idq
      go to 229
  234 icq = ideks(idq)+icq
      go to 229
  230 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 235
  223 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 236
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 237,238,239
  239 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 240
  238 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 241
  237 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  245 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 229
  236 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 242,243,244
  244 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 229
  243 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 233
  242 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 245
  222 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
  247 if (nl2r.lt.nm2r) go to 248
      iby=ideks(nl2r)+nm2r
      go to 249
  248 iby=ideks(nm2r)+nl2r
  249 if (nt1r.eq.nt2r) go to 250
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 251
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 229
  251 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 229
  250 icx=ideks(iax+1)
      if (iay.lt.iby) go to 253
      icq=ideks(iay)+iby
      go to 229
  253 icq=ideks(iby)+iay
 229  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
260   if(kix.lt.0)  go to 113
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 220
 114  jz=l4
      do 150 l=1,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 151,150,152
  150 continue
  151 jp=l
      if (l.eq.nr) go to 153
      ia=l+1
      do 154 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,154,155
  154 continue
      go to 153
  155 kp=l
      if (l.eq.nr) go to 156
      ia=l+1
      do 157 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  157 continue
      go to 156
  152 jp=l
      if (l.eq.nr) go to 158
      ia=l+1
      do 159 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 160,159,112
  159 continue
      go to 158
  160 kp=l
      if (l.eq.nr) go to 161
      ia=l+1
      do 162 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  162 continue
      go to 161
  153 mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 163
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 163  jz=nt
      do 164 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 165
 164  continue
  165 mm=mm+1
      olab(mm)=l
      if (mm.lt.nnid) go to 166
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  166 mm=mm+1
      olab(mm)=jp
      if (mm.lt.nnid) go to 167
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  167 jz=nz
      if (ip.eq.1) go to 168
      kz=l3
      ja=ip-1
      do 169 l=1,ja
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  169 continue
  168 if (ip.eq.nd) go to 171
      ja=ip+1
      kz=l3+ip
      do 172 l=ja,nd
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  172 continue
  171 kk=iot(jz+1)
  170 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
      go to 20
  156 mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 173
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  173 jz=nt
      mm=mm+1
      kz=l4
      ig=0
      if (jp.eq.1) go to 174
      ka=1
      ja=jp-1
  176 do 175 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 175
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 800
 175  continue
      go to 174
 800  ka=l
      kz=kz-1
      go to 176
 174  la=kp-1
      if (jp.eq.la) go to 178
      ja=jp+1
      kz=kz+1
 180  do 179 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 179
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 801
 179  continue
      go to 178
 801  ja=l
      kz=kz-1
      go to 180
 178  if (kp.eq.nr) go to 181
      ia=kp+1
      kz=l4+kp
 182  do 183 l=ia,nr
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 183
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 802
 183  continue
      go to 181
 802  ia=l
      kz=kz-1
      go to 182
 181  lab(2)=nr
      if (ig.eq.1) go to 177
      lab(1)=nr1
 177  ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 184
      olab(mm)=2
 188  if (mm.lt.nnid) go to 185
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  185 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nnid) go to 186
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  186 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nnid) go to 187
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 187  nt1r=nir(nn)
      nt2r=nir(ll)
      nl1r=loc(nn)
      nl2r=loc(ll)
      go to 189
 184  olab(mm)=3
      jz=nt+ja
      jj=iot(jz)
      go to 188
  161 mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 190
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  190 jz=nt
      mm=mm+1
      kz=l4
      ig=0
      if (jp.eq.1) go to 191
      ka=1
      ja=jp-1
 192  do 193 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 193
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 803
 193  continue
      go to 191
 803  ka=l
      kz=kz-1
      go to 192
 191  la=kp-1
      if (jp.eq.la) go to 194
      ja=jp+1
      kz=kz+1
 195  do 196 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 196
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 804
196   continue
      go to 194
 804  ja=l
      kz=kz-1
      go to 195
 194  if (kp.eq.nr) go to 197
      ia=kp+1
      kz=l4+kp
 198  do 199 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 199
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 805
  199 continue
      go to 197
  805 ia=l
      kz=kz-1
      go to 198
  197 lab(2)=nr
      if (ig.eq.1) go to 806
      lab(1)=nr-1
  806 ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 501
      olab(mm)=3
 502  if (mm.lt.nnid) go to 503
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  503 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nnid) go to 504
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  504 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nnid) go to 505
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 505  nt1r=nir(nn)
      nt2r=nir(jj)
      nl1r=loc(nn)
      nl2r=loc(jj)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(ll)
      nm1r=loc(kk)
      nm2r=loc(ll)
      go to 220
 501  olab(mm)=2
      jz=nt+ja
      jj=iot(jz)
      go to 502
 158  mm=mm+1
      olab(mm)=4
      if (mm.lt.nnid) go to 506
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 506  jz=nt
      do 507 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 508
 507  continue
 508  mm=mm+1
      olab(mm)=-l
      if (mm.lt.nnid) go to 509
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 509  mm=mm+1
      olab(mm)=jp
      if (mm.lt.nnid) go to 510
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 510  mm=mm+1
      iax=nir(ll)
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      olab(mm)=iax
      if (mm.lt.nnid) go to 893
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 893  mm=mm+1
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      kl=km+kk
      nm=ideks(ii)
      nn=nm+ii
      if (ii.gt.kk) go to 513
      olab(mm)=kk
      if (mm.lt.nnid) go to 511
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 511  mm=mm+1
      olab(mm)=ii
      if (mm.lt.nnid) go to 512
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 512  jb=km+ii
      ib=ideks(jb)+nn+iat
      kb=ideks(kl)+jb+iat
      go to 514
 513  jb=nm+kk
      olab(mm)=ii
      if (mm.lt.nnid) go to 551
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 551  mm=mm+1
      olab(mm)=kk
      if (mm.lt.nnid) go to 552
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 552  ib=ideks(jb)+kl+iat
      kb=ideks(nn)+jb+iat
 514  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 515
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 515  idud=(kb-1)/igmax
      kb=kb-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=kb
      nix=1
      if (jto.lt.iwod) go to 516
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 516  if (nr.eq.1) go to 517
      jz=l4
      kix=1
      lp=jp
      ir=nr
 518  do 384 l=1,ir
      jz=jz+1
      if (l.eq.lp) go to 384
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  if(kb.lt.lb) go to 393
      ib=iat+ideks(kb) +lb
      go to 394
 393  ib=iat+ideks(lb) +kb
 394  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      if(kb.lt.lb) go to 398
      ib=ideks(kb)+lb+ibx
      go to 471
 398  ib=ideks(lb) +kb+ibx
 471  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      if(kb.lt.lb) go to 474
      ib=ideks(kb)+lb+ibx
      go to 475
 474  ib=ideks(lb)+kb+ibx
 475  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 113
      if (nix.lt.0) go to 519
 517  if (nd.eq.1) go to 113
      kix=-1
      jz=l3
      lp=ip
      ir=nd
      go to 518
  519 if (nd.eq.0) go to 113
      kix=-1
      jz=l3
      lp=0
      ir=nd
      go to 518
  106 jz=l4
      do 520 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 521,520,112
  520 continue
  521 ip=l
      if (l.eq.nr) go to 522
      ia=l+1
      do 523 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 524,523,112
  523 continue
      go to 522
  524 jp=l
      if (l.eq.nr) go to 525
      ia=l+1
      do 526 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  526 continue
  525 mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 527
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  527 mm=mm+1
      olab(mm)=1
      if (mm.lt.nnid) go to 528
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  528 mm=mm+1
      jz=nt
      kz=l4
      ig=0
      if (ip.eq.1) go to 529
      ka=1
      ja=ip-1
  531 do 530 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 530
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 807
  530 continue
      go to 529
  807 ka=l
      kz=kz-1
      go to 531
  529 la=jp-1
      if (ip.eq.la) go to 533
      ja=ip+1
      kz=kz+1
 534  do 535 l=ja,la
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 535
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 808
  535 continue
      go to 533
  808 ja=l
      kz=kz-1
      go to 534
  533 if (jp.eq.nr) go to 536
      ia=jp+1
      kz=l4+jp
  537 do 538 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 538
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 809
  538 continue
      go to 536
  809 ia=l
      kz=kz-1
      go to 537
  536 lab(2)=nr
      if (ig.eq.1) go to 532
      lab(1)=nr1
  532 l1=lab(1)
      jz=nt+l1
      kk=iot(jz)
      l2=lab(2)
      jz=nt+l2
      jj=iot(jz)
      olab(mm)=ideks(l2-1)+l1
      if (mm.lt.nnid) go to 539
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  539 mm=mm+1
      olab(mm)=ideks(jp-1)+ip
      go to 540
  522 mm=mm+1
      olab(mm)=4
      if (mm.lt.nnid) go to 541
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  541 mm=mm+1
      jz=nt
      if (nr.eq.1) go to 546
      kz=l4
      if (ip.eq.1) go to 542
      ja=ip-1
      do 543 l=1,ja
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  543 continue
      if (ip.eq.nr) go to 546
  542 ja=ip+1
      kz=kz+1
      do 547 l=ja,nr
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  547 continue
  546 nn=iot(jz+1)
      l1=nr
      go to 548
  544 l1=jz-nt
  548 olab(mm)=l1
      if (mm.lt.nnid) go to 549
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  549 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nnid) go to 550
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  550 mm=mm+1
      iax=nir(ll)
      olab(mm)=iax
      if (mm.lt.nnid) go to 897
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 897  mm=mm+1
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      nm=ideks(ii)
      if (ii.gt.kk) go to 560
      olab(mm)=kk
      if (mm.lt.nnid) go to 561
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  561 mm=mm+1
      olab(mm)=ii
      if (mm.lt.nnid) go to 562
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  562 jb=km+ii
      go to 565
  560 olab(mm)=ii
      if (mm.lt.nnid) go to 563
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  563 mm=mm+1
      olab(mm)=kk
      if (mm.lt.nnid) go to 564
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  564 jb=nm+kk
  565 nix=-1
      if (nr.eq.1) go to 519
      jz=l4
      kix=1
      lp=ip
      ir=nr
      go to 518
  113 l3=l3+nnx
  105 l4=l4+nnx
      return
cdel 54   format(4x,40i3)
      end
*************************************************************************
      subroutine bearz(iot,niott,ideks,ndeks)
c
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndeks
      dimension iot(niott),ideks(ndeks)
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      common /a/ t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
      parameter (nmmo=256)
c
c     nir(nmmo) : symmetry of mo's n : nir(n)
      integer loc, nir, jcon
      common /cloc/ loc(nmmo), nir(nmmo), jcon(nmmo)
c --- jcon : kodierte konfigurationen (0,1,2 - kode)
c
      integer lj
      common /clj2/ lj(8)
c
      common /cnit/   nit(nitmax)
      integer kc,kd,lab
      common /ckd/ kc(ndimh),kd(ndimh),lab(nn3)
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref),ie(maxref),nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
c
      if (nc.eq.1) return
cdebugdo i=1,9
cdebug write(6,*)'ideks',ideks(i),'jdeks',jdeks(i),kdeks(i)
cdebugend do
cdebugdo i=1,nmmo
cdebug write(6,*)'loc',loc(i),'  nir',nir(i),' jcon',jcon(i)
cdebugend do
cdebug
cdebug  write(6,*)'bearz'
cdebugwrite(6,*)'ntape =        ',ntape
cdebugwrite(6,*) 'lt=',lt,'  nnx=',nnx,'  nt=',nt
c
      nt = lt + nnx
      nz = lz + nnx
      j1 = nzero
c --- calculate jcon : coding with 0,1,2 as usual
      do 100 j=2,nc
       j1 = j1 + 1
       it = nt
       mt = lt
       mz = lz
       if (nr.eq.0) go to 102
       do 1191 k=1,nr
        it = it + 1
        ll = iot(it)
        jcon(ll) = 1
1191   continue
       if (nd.eq.0) go to 103
 102   continue
       do 1198 k=1,nd
        it = it + 1
        ll = iot(it)
        jcon(ll)=2
1198   continue
 103   continue
       do 105 k=1,j1
        if (nd.eq.0) go to 106
        jz = mz
        do 107 l=1,nd
         jz = jz + 1
         ll = iot(jz)
         if (jcon(ll)-1) 108,109,107
 107    continue
       go to 106
c --- jcon(ll) = 0 (ll not occupied)
 108   continue
       ip = l
       if (l.eq.nd) go to 110
       ia = l + 1
       do 111 l=ia,nd
        jz = jz + 1
        jj = iot(jz)
        if (jcon(jj).ne.2) go to 112
  111  continue
       go to 110
  112  mm = mm + 1
       olab(mm) = 1
       if (mm.lt.nnid) go to 113
       mm = nzero
c
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       goto 113
c --- jcon(ll) = 1 (ll singly occupied)
 109   continue
       ip = l
       if (l.eq.nd) go to 114
       ia = l + 1
       do 115 l=ia,nd
        jz = jz + 1
        nn = iot(jz)
        if (jcon(nn)-1) 112,116,115
115    continue
       go to 114
116    if (l.eq.nd) go to 117
       ia = l + 1
       do 118 l=ia,nd
        jz = jz + 1
        kk = iot(jz)
        if (jcon(kk).ne.2)  go to 112
 118   continue
       go to 117
 110   if (nr.eq.0) go to 119
       jz = mt
       do 120 l=1,nr
         jz = jz + 1
         jj = iot(jz)
         if (jcon(jj)-1) 112,120,121
 120   continue
       go to 119
 121   ip = l
       if (l.eq.nr) go to 122
       ia = l + 1
       do 123 l=ia,nr
        jz = jz + 1
        kk = iot(jz)
        if (jcon(kk).ne.1) go to 112
 123   continue
       go to 122
 119   mm = mm + 1
       olab(mm) = 5
       if (mm.lt.nnid) go to 124
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 124   jz = nz
       if (ip.eq.1) go to 125
       kz = mz
       ja = ip-1
       do 1279 l=1,ja
         kz = kz+1
         jz = jz+1
         kk = iot(jz)
         if (kk.ne.iot(kz)) go to 127
1279  continue
 125  if (ip.eq.nd) go to 128
      ja = ip + 1
      kz = mz + ip
      do 1288 l=ja,nd
       kz = kz + 1
       jz = jz + 1
       kk = iot(jz)
       if (kk.ne.iot(kz)) go to 127
1288  continue
 128  kk = iot(jz+1)
 127  mm = mm + 1
      olab(mm) = ll
      if (mm.ge.nnid) then
        mm = nzero
        call pack(olab8,8,olab,nnid)
        write (ntape) olab8
      end if
      mm = mm + 1
      olab(mm) = kk
      if (mm.lt.nnid) go to 113
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 113
 122  mm = mm + 1
      olab(mm) = 3
      if (mm.lt.nnid) go to 131
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  131 jz=nt
      if (ip.eq.1) go to 132
      kz = mt
      ja = ip - 1
      do 133 l=1,ja
       kz = kz +1
       jz = jz + 1
       kk = iot(jz)
       if (kk.ne.iot(kz)) go to 134
 133  continue
 132  if (ip.eq.nr) go to 650
      ja = ip + 1
      kz = mt + ip
      do 651 l=ja,nr
      kz = kz + 1
      jz = jz + 1
      kk = iot(jz)
      if (kk.ne.iot(kz)) go to 134
 651  continue
 650  jz = jz + 1
      kk = iot(jz)
  134 mm = mm + 1
      olab(mm) = jz - nt
      if (mm.lt.nnid) go to 135
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  135 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nnid) go to 136
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  136 nt1r = nir(ll)
      nl1r = loc(ll)
      nb1r = nir(jj)
      nm1r = loc(jj)
      nm2r = loc(kk)
   20 if (nt1r-nb1r) 23,22,99
   99 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq + nm2r
      icq = icq + nm1r
   33 if (icq.lt.idq) go to 34
      icq = ideks(icq) + idq
      go to 29
   34 icq = ideks(idq) + icq
      go to 29
   23 iax = ideks(nb1r) + nt1r
      icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      go to 33
   22 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 46
      iay=ideks(nl1r)+nm1r
      go to 47
   46 iay= ideks(nm1r)+nl1r
   47 if (nl1r.lt.nm2r) go to 48
      iby=ideks(nl1r)+nm2r
      go to 50
   48 iby=ideks(nm2r)+nl1r
   50 icx=ideks(iax+1)
      if (iay.lt.iby) go to 53
      icq=ideks(iay)+iby
      go to 29
   53 icq=ideks(iby)+iay
 29   icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 113
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 113
 117  jz=mt
      do 137 l=1,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,137,138
  137 continue
  138 l1=l
      ia=l+1
      do 139 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jcon(jj)-1) 112,139,140
  139 continue
  140 l2=l
      if (l.eq.nr) go to 141
      ia=l+1
      do 142 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
 142  continue
 141  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 143
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  143 mm=mm+1
      olab(mm)=1
      if (mm.lt.nnid) go to 144
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 144  jz = nt
      do 145 l=1,nr
       jz = jz + 1
       if (iot(jz).eq.ll) go to 146
 145  continue
 146  ip = l
      ia = l + 1
      do 1433 l=ia,nr
       jz = jz + 1
       if (iot(jz).eq.nn) go to 147
1433  continue
 147  mm = mm + 1
      olab(mm)=ip+ideks(l-1)
      if (mm.ge.nnid) then
        mm = nzero
        call pack(olab8,8,olab,nnid)
        write (ntape) olab8
      end if
      mm = mm + 1
      olab(mm)=ideks(l2-1)+l1
 540  if (mm.lt.nnid) go to 149
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 149  nt1r = nir(ll)
      nt2r = nir(nn)
      nl1r = loc(ll)
      nl2r = loc(nn)
 189  nb1r = nir(kk)
      nb2r = nir(jj)
      nm1r = loc(kk)
      nm2r = loc(jj)
      kix  = 3
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 224
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 225,226,227
  227 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  240 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 229
  226 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  241 idq = (nl2r-1)*ljn+nm2r
      go to 233
  225 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  235 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 229
  224 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 230,231,232
  232 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 229
  231 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  233 if (icq.lt.idq) go to 234
      icq = ideks(icq)+idq
      go to 229
  234 icq = ideks(idq)+icq
      go to 229
  230 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 235
  223 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 236
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 237,238,239
  239 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 240
  238 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 241
  237 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  245 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 229
  236 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 242,243,244
  244 icx = ideks(iax)+ibx
      icq = (nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 229
  243 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 233
  242 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 245
  222 iax = ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 246
      iay = ideks(nl1r)+nm1r
      go to 247
 246  iay = ideks(nm1r) + nl1r
 247  if (nl2r.lt.nm2r) go to 248
      iby = ideks(nl2r) + nm2r
      go to 249
 248  iby = ideks(nm2r) + nl2r
 249  if (nt1r.eq.nt2r) go to 250
      ibx = ideks(nt2r+1)
      if (iax.lt.ibx) go to 251
      icx = ideks(iax)+ibx
      ljn = lj(nt2r)
      icq = (iay-1)*ideks(ljn+1) + iby
      go to 229
 251  icx = ideks(ibx)+iax
      ljn = lj(nt1r)
      icq = (iby-1)*ideks(ljn+1)+iay
      go to 229
 250  icx = ideks(iax+1)
      if (iay.lt.iby) go to 253
      icq = ideks(iay)+iby
      go to 229
 253  icq  = ideks(iby)+iay
 229  icq  = nit(icx)  +  icq
      idud = (icq -1)/igmax
      jto  = jto + 1
      kc(jto) = idud  + 1
      icq  = icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto = nzero
      jblk= jblk  + 1
 260  if(kix.lt.0)  go to 113
      kix  = -3
      itr  = nb1r
      nb1r = nb2r
      nb2r = itr
      itr  = nm1r
      nm1r = nm2r
      nm2r = itr
      go to 220
 114  jz   = mt
      do 150 l=1,nr
        jz = jz + 1
        nn = iot(jz)
        if (jcon(nn)-1) 151,150,152
 150  continue
 151  jp = l
      if (l.eq.nr) go to 153
      ia=l+1
      do 154 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk)-1) 112,154,155
  154 continue
      go to 153
  155 kp = l
      if (l.eq.nr) go to 156
      ia = l + 1
      do 157 l=ia,nr
        jz = jz + 1
        ii = iot(jz)
        if (jcon(ii).ne.1) go to 112
  157 continue
      go to 156
  152 jp = l
      if (l.eq.nr) go to 158
      ia = l + 1
      do 159 l=ia,nr
        jz = jz + 1
        kk = iot(jz)
        if (jcon(kk)-1) 160,159,112
  159 continue
      go to 158
  160 kp = l
      if (l.eq.nr) go to 161
      ia = l + 1
      do 162 l=ia,nr
        jz = jz + 1
        ii = iot(jz)
        if (jcon(ii).ne.1) go to 112
  162 continue
      go to 161
  153 mm = mm + 1
      olab(mm) = 3
      if (mm.ge.nnid) then
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
      end if
      jz = nt
      do 164 l=1,nr
        jz = jz + 1
        if (iot(jz).eq.ll) go to 165
 164  continue
c
 165  mm = mm + 1
      olab(mm) = l
      if (mm.lt.nnid) go to 166
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 166  mm = mm + 1
      olab(mm)=jp
      if (mm.lt.nnid) go to 167
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  167 jz=nz
      if (ip.eq.1) go to 168
      kz=mz
      ja=ip-1
      do 169 l=1,ja
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  169 continue
  168 if (ip.eq.nd) go to 171
      ja=ip+1
      kz=mz+ip
      do 172 l=ja,nd
      jz=jz+1
      kz=kz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 170
  172 continue
  171 kk = iot(jz+1)
  170 nt1r = nir(kk)
      nl1r = loc(kk)
      nb1r = nir(ll)
      nm1r = loc(ll)
      nm2r = loc(nn)
      go to 20
  156 mm = mm + 1
      olab(mm)=2
      if (mm.lt.nnid) go to 173
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  173 jz = nt
      mm = mm + 1
      kz = mt
      ig = nzero
      if (jp.eq.1) go to 174
      ka=1
      ja=jp-1
  176 do 175 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 175
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 800
 175  continue
      go to 174
 800  ka=l
      kz=kz-1
      go to 176
 174  la=kp-1
      if (jp.eq.la) go to 178
      ja=jp+1
      kz=kz+1
 180  do 179 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 179
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 801
 179  continue
      go to 178
 801  ja=l
      kz=kz-1
      go to 180
 178  if (kp.eq.nr) go to 181
      ia=kp+1
      kz=mt+kp
 182  do 183 l=ia,nr
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 183
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 177
      go to 802
 183  continue
      go to 181
 802  ia=l
      kz=kz-1
      go to 182
 181  lab(2)=nr
      if (ig.eq.1) go to 177
      lab(1)=nr1
 177  ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 184
      olab(mm)=2
 188  if (mm.lt.nnid) go to 185
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  185 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nnid) go to 186
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  186 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nnid) go to 187
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 187  nt1r=nir(nn)
      nt2r=nir(ll)
      nl1r=loc(nn)
      nl2r=loc(ll)
      go to 189
 184  olab(mm)=3
      jz=nt+ja
      jj=iot(jz)
      go to 188
  161 mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 190
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  190 jz=nt
      mm=mm+1
      kz=mt
      ig=0
      if (jp.eq.1) go to 191
      ka=1
      ja=jp-1
 192  do 193 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 193
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 803
 193  continue
      go to 191
 803  ka=l
      kz=kz-1
      go to 192
 191  la=kp-1
      if (jp.eq.la) go to 194
      ja=jp+1
      kz=kz+1
 195  do 196 l=ja,la
      kz=kz+1
      jz=jz+1
      if (iot(kz).eq.iot(jz)) go to 196
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 804
196   continue
      go to 194
 804  ja=l
      kz=kz-1
      go to 195
 194  if (kp.eq.nr) go to 197
      ia=kp+1
      kz=mt+kp
 198  do 199 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 199
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 806
      go to 805
  199 continue
      go to 197
  805 ia=l
      kz=kz-1
      go to 198
  197 lab(2)=nr
      if (ig.eq.1) go to 806
      lab(1)=nr-1
  806 ia=lab(1)
      jz=ia+nt
      ja=lab(2)
      jj=iot(jz)
      if (jj.eq.ll) go to 501
      olab(mm)=3
 502  if (mm.lt.nnid) go to 503
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  503 mm=mm+1
      olab(mm)=ideks(ja-1)+ia
      if (mm.lt.nnid) go to 504
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  504 mm=mm+1
      olab(mm)=ideks(kp-1)+jp
      if (mm.lt.nnid) go to 505
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 505  nt1r=nir(nn)
      nt2r=nir(jj)
      nl1r=loc(nn)
      nl2r=loc(jj)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(ll)
      nm1r=loc(kk)
      nm2r=loc(ll)
      go to 220
 501  olab(mm)=2
      jz=nt+ja
      jj=iot(jz)
      go to 502
 158  mm=mm+1
      olab(mm)=4
      if (mm.lt.nnid) go to 506
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 506  jz=nt
      do 507 l=1,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 508
 507  continue
 508  mm=mm+1
      olab(mm)=-l
      if (mm.lt.nnid) go to 509
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 509  mm=mm+1
      olab(mm)=jp
      if (mm.lt.nnid) go to 510
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 510  mm=mm+1
      iax=nir(ll)
      iay=ideks(iax+1)
      iaz=ideks(iay+1)
      olab(mm)=iax
      if (mm.lt.nnid) go to 893
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 893  mm=mm+1
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(ll)
      kk=loc(nn)
      km=ideks(kk)
      kl=km+kk
      nm=ideks(ii)
      nn=nm+ii
      if (ii.gt.kk) go to 513
      olab(mm)=kk
      if (mm.lt.nnid) go to 511
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 511  mm=mm+1
      olab(mm)=ii
      if (mm.lt.nnid) go to 512
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 512  jb=km+ii
      ib=ideks(jb)+nn+iat
      kb=ideks(kl)+jb+iat
      go to 514
 513  jb=nm+kk
      olab(mm)=ii
      if (mm.lt.nnid) go to 551
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 551  mm=mm+1
      olab(mm)=kk
      if (mm.lt.nnid) go to 552
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 552  ib=ideks(jb)+kl+iat
      kb=ideks(nn)+jb+iat
 514  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 515
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 515  idud=(kb-1)/igmax
      kb=kb-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=kb
      nix=1
      if (jto.lt.iwod) go to 516
      jto=0
      write (ntype) kc,kd
      jblk=jblk+1
 516  if (nr.eq.1) go to 517
      jz=mt
      kix=1
      lp=jp
      ir=nr
 518  do 384 l=1,ir
      jz=jz+1
      if (l.eq.lp) go to 384
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  if(kb.lt.lb) go to 393
      ib=iat+ideks(kb) +lb
      go to 394
 393  ib=iat+ideks(lb) +kb
 394  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto     = jto+1
      kc(jto) = idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      if(kb.lt.lb) go to 398
      ib=ideks(kb)+lb+ibx
      go to 471
 398  ib=ideks(lb) +kb+ibx
 471  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      if(kb.lt.lb) go to 474
      ib=ideks(kb)+lb+ibx
      go to 475
 474  ib=ideks(lb)+kb+ibx
 475  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 113
      if (nix.lt.0) go to 519
 517  if (nd.eq.1) go to 113
      kix=-1
      jz=mz
      lp=ip
      ir=nd
      go to 518
  519 if (nd.eq.0) go to 113
      kix=-1
      jz=mz
      lp=0
      ir=nd
      go to 518
  106 jz=mt
      do 520 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 521,520,112
  520 continue
  521 ip=l
      if (l.eq.nr) go to 522
      ia=l+1
      do 523 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 524,523,112
  523 continue
      go to 522
  524 jp=l
      if (l.eq.nr) go to 525
      ia=l+1
      do 526 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 112
  526 continue
  525 mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 527
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      write(6,*)'olab8:',(olab(i),i=1,10)
  527 mm=mm+1
      olab(mm)=1
      if (mm.lt.nnid) go to 528
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  528 mm=mm+1
      jz=nt
      kz=mt
      ig=0
      if (ip.eq.1) go to 529
      ka=1
      ja=ip-1
  531 do 530 l=ka,ja
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 530
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 807
  530 continue
      go to 529
  807 ka=l
      kz=kz-1
      go to 531
  529 la=jp-1
      if (ip.eq.la) go to 533
      ja=ip+1
      kz=kz+1
 534  do 535 l=ja,la
      jz=jz+1
      kz=kz+1
      if (iot(jz).eq.iot(kz)) go to 535
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 808
  535 continue
      go to 533
  808 ja=l
      kz=kz-1
      go to 534
  533 if (jp.eq.nr) go to 536
      ia=jp+1
      kz=mt+jp
  537 do 538 l=ia,nr
      kz=kz+1
      jz=jz+1
      if (iot(jz).eq.iot(kz)) go to 538
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.2) go to 532
      go to 809
  538 continue
      go to 536
  809 ia=l
      kz=kz-1
      go to 537
  536 lab(2)=nr
      if (ig.eq.1) go to 532
      lab(1)=nr1
  532 l1=lab(1)
      jz=nt+l1
      kk=iot(jz)
      l2=lab(2)
      jz=nt+l2
      jj=iot(jz)
      olab(mm)=ideks(l2-1)+l1
      if (mm.lt.nnid) go to 539
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  539 mm=mm+1
      olab(mm)=ideks(jp-1)+ip
      go to 540
  522 mm=mm+1
      olab(mm)=4
      if (mm.lt.nnid) go to 541
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  541 mm=mm+1
      jz=nt
      if (nr.eq.1) go to 546
      kz=mt
      if (ip.eq.1) go to 542
      ja=ip-1
      do 543 l=1,ja
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  543 continue
      if (ip.eq.nr) go to 546
  542 ja=ip+1
      kz=kz+1
      do 547 l=ja,nr
      jz=jz+1
      kz=kz+1
      nn=iot(jz)
      if (nn.ne.iot(kz)) go to 544
  547 continue
  546 nn = iot(jz+1)
      l1 = nr
      go to 548
  544 l1=jz-nt
  548 olab(mm)=l1
      if (mm.lt.nnid) go to 549
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  549 mm=mm+1
      olab(mm)=ip
      if (mm.lt.nnid) go to 550
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  550 mm=mm+1
      iax=nir(ll)
      olab(mm)=iax
      if (mm.lt.nnid) go to 897
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 897  mm  = mm + 1
      iay = ideks(iax+1)
      iaz = ideks(iay+1)
      mq  = lj(iax)
      mv  = ideks(mq+1)
      iaq = iay - iax
      iaw = iaz - iay
      iat = nit(iaz)
      ii  = loc(ll)
      kk  = loc(nn)
      km  = ideks(kk)
      nm  = ideks(ii)
      if (ii.gt.kk) go to 560
      olab(mm) = kk
      if (mm.lt.nnid) go to 561
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  561 mm = mm + 1
      olab(mm) = ii
      if (mm.lt.nnid) go to 562
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  562 jb = km + ii
      go to 565
  560 olab(mm) = ii
      if (mm.lt.nnid) go to 563
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  563 mm = mm + 1
      olab(mm)=kk
      if (mm.lt.nnid) go to 564
      mm = nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  564 jb = nm + kk
  565 nix = -1
      if (nr.eq.1) go to 519
      jz  = mt
      kix = 1
      lp  = ip
      ir  = nr
      go to 518
  113 mz = mz + nnx
  105 mt = mt + nnx
      do 2291 k=1,nnx
       nt = nt + 1
       ll = iot(nt)
       jcon(ll) = nzero
2291  continue
  100 nz = nt + nr
c
      return
      end
*************************************************************************
      subroutine bruna(jdeks,ndeks)
c
      implicit REAL (a-h,o-z)
c
      integer jdeks,ndeks
      dimension jdeks(ndeks)
c
INCLUDE(common/newmrd_parinc)
c
       common /ftap/ ntape, mtape, mdisk,
     .               ideli, ltype, linf,  ntab,  kfile,
     .               kclab, ntype, mstvt, nf01, nf62,
     +               nhead, nf99,  mtapev, nf11
       common /tap/ ical,  m,     nko,
     .              mxex,  nmul,  ispace, nprin, ncorci,
     +              msng,  maxci2
c
      common /a/ b(nlca )
      common /cnbox/ nbox
c
      parameter (nmmo=256)
      integer jab, isym, mj, iab, kj, jcon, icuf, iduf
      integer kcon, ilee, ntil, nbal, jkon, nop
      integer isw, nyzl, nj, ijoe, ndub
      common /aeng/  nconf(ndk5),jab(n36),isym(maxsym),
     .               mj(maxsym),iab(nmmo),kj(maxsym),
     .               jcon(nmmo),icuf(n12),iduf(maxshl),
     .               kcon(ndimh),ilee(maxsym+1),
     .               ntil(maxsym),nbal(maxsym+1),
     .               jkon(maxref*maxshl),nop(maxref),
     .               isw(ndk5),    nyzl(maxref),  nj(maxsym),
     .               ijoe(maxref), ndub(maxref)
c
      common /b1/ ig,jpaw,klx,nshl,kly,ispin
c     integer jdeks,kdeks
c     common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
c
cccc replace ibl !!!
      ibl = ndimh
ccccc
      nbxa= nbox - 1
      nbxb= nbox - 2
      nbd = 1
      kar = nzero
414   kar = kar + 1
      if(kar.gt.nbox) go to 415
      if(jcon(kar).eq.2)  go to 414
      iq = iab(kar)
      if(iq.gt.jpaw) go to 416
      kpaw = jdeks(jpaw) + iq
      go to 417
416   kpaw = jdeks(iq) + jpaw
417   kpaw = jab(kpaw)
      kpaw = jdeks(kpaw)
      if(kpaw.eq.0) go to 414
      kas  = ilee(kpaw)
      npox = ilee(kpaw+ 1)
      jcon(kar) =jcon(kar)  + 1
418   kas = kas + 1
      if(kas.gt.npox) go to 419
      if(kas.eq.kar) go to 418
      if(jcon(kas).eq.0) go to 418
      jcon(kas) = jcon(kas) - 1
      lar=0
420   lar=lar + 1
      if(lar.eq.kar.or.lar.eq.kas) go to 420
      if(lar.gt.nbox) go to 421
      if(jcon(lar).lt.2) go to 420
      jcon(lar)=0
      las=0
422   las=las + 1
      if(las.eq.kar.or.las.eq.kas) go to 422
      if(las.eq.lar) go to 422
      if(las.gt.nbox) go to 423
      if(jcon(las).gt.0) go to 422
      jcon(las) = 2
      go to 33
423   jcon(lar) =2
      go to 420
421   jcon(kas) =jcon(kas)  + 1
      go to 418
419   jcon(kar)= jcon(kar) -1
      go to 414
400   jcon(las) = 0
      go to 422
415   nbd=2
      kar=0
424   kar = kar + 1
      if(kar.eq.nbxa) go to 425
      if(jcon(kar).eq.0 ) go to 424
      jcon(kar) = jcon(kar) - 1
      iq= iab(kar)
      if(iq.gt.jpaw) go to 426
      kpaw= jdeks(jpaw) + iq
      go to 427
426   kpaw=jdeks(iq)  + jpaw
427   kpaw= jab(kpaw)
      lar  = kar
428   lar = lar + 1
      if(lar.eq.nbox) go to 429
      if(jcon(lar).eq.0) go to 428
      jcon(lar) = jcon(lar) - 1
      iq=iab(lar)
      if(iq.gt.kpaw) go to 430
       lpaw=jdeks(kpaw)  + iq
      go to 431
430   lpaw= jdeks(iq)  + kpaw
431   lpaw= jab(lpaw)
      mar= lar
432   mar= mar + 1
      if(mar.gt.nbox) go to 433
      if(jcon(mar).eq.0) go to 432
      iq=iab(mar)
      if(iq.gt.lpaw) go to 434
      mpaw=jdeks(lpaw) + iq
      go to 435
434   mpaw= jdeks(iq)  + lpaw
435   mpaw=jab(mpaw)
      mpaw=jdeks(mpaw)
      if(mpaw.eq.0) go to 432
      kas=ilee(mpaw)
      npox=ilee(mpaw + 1)
      jcon(mar) =jcon(mar) - 1
436   kas=kas + 1
      if(kas.gt.npox) go to 437
      if(kas.eq.kar.or.kas.eq.lar) go to 436
      if(kas.eq.mar.or.jcon(kas).eq.2) go to 436
      jcon(kas) = jcon(kas) + 1
      las=0
438   las = las + 1
      if(las.gt.nbox) go to 439
      if(las.eq.kar.or.las.eq.lar) go to 438
      if(las.eq.mar.or.las.eq.kas) go to 438
      if(jcon(las).ne.0) go to 438
      jcon(las) =2
      go to 33
439   jcon(kas) =jcon(kas) - 1
      go to 436
401   jcon(las) =0
      go to 438
437   jcon(mar) = jcon(mar) + 1
      go to 432
433   jcon(lar) = jcon(lar) + 1
      go to 428
429   jcon(kar) =jcon(kar) + 1
      go to 424
425   nbd=3
      kar=0
440   kar=kar + 1
      if(kar.eq.nbxa) go to 441
      if(jcon(kar).eq.2) go to 440
      jcon(kar) = jcon(kar) + 1
      iq=iab(kar)
      if(iq.gt.jpaw) go to 442
      kpaw=jdeks(jpaw)  + iq
      go to 443
442   kpaw=jdeks(iq) + jpaw
443   kpaw=jab(kpaw)
      lar=kar
444   lar= lar + 1
      if(lar.eq.nbox) go to 445
      if(jcon(lar).eq.2) go to 444
      jcon(lar) = jcon(lar) + 1
      iq=iab(lar)
      if(iq.gt.kpaw) go to 446
      lpaw=jdeks(kpaw)  + iq
      go to 447
446   lpaw=jdeks(iq) + kpaw
447   lpaw=jab(lpaw)
      mar=lar
448   mar= mar + 1
      if(mar.gt.nbox) go to 449
      if(jcon(mar).eq.2) go to 448
      iq=iab(mar)
      if(iq.gt.lpaw) go to 450
      mpaw=jdeks(lpaw)  + iq
      go to 451
450   mpaw=jdeks(iq) + lpaw
451   mpaw=jab(mpaw)
      mpaw=jdeks(mpaw)
      if(mpaw.eq.0) go to 448
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      jcon(mar)=jcon(mar)  + 1
452   kas = kas + 1
      if(kas.gt.npox) go to 453
      if(kas.eq.kar.or.kas.eq.lar) go to 452
      if(kas.eq.mar.or.jcon(kas).eq.0) go to 452
      jcon(kas) = jcon(kas) - 1
      las = nzero
454   las = las + 1
      if(las.gt.nbox) go to 455
      if(las.eq.kar.or.las.eq.lar ) go to 454
      if(las.eq.mar.or.las.eq.kas) go to 454
      if(jcon(las).ne.2) go to 454
      jcon(las) = 0
      go to 33
455   jcon(kas) =jcon(kas) + 1
      go to 452
402   jcon(las) = 2
      go to 454
453   jcon(mar) = jcon(mar) - 1
      go to 448
449   jcon(lar) = jcon(lar) - 1
      go to 444
445   jcon(kar) = jcon(kar) - 1
      go to 440
441   nbd = 4
      kar=0
456   kar = kar + 1
      if(kar.eq.nbxa) go to 582
      if(jcon(kar).eq.2) go to 456
      jcon(kar) = jcon(kar)  + 1
      iq=iab(kar)
      if(iq.gt.jpaw) go to 458
      kpaw=jdeks(jpaw)  + iq
      go to 459
458   kpaw=jdeks(iq)  + jpaw
459   kpaw=jab(kpaw)
      lar= kar
460   lar=lar + 1
      if(lar.eq.nbox) go to 461
      if(jcon(lar).eq.2) go to 460
      jcon(lar) = jcon(lar ) + 1
      iq=iab(lar)
      if(iq.gt.kpaw) go to 462
      lpaw=jdeks(kpaw)  + iq
      go to 463
462   lpaw=jdeks(iq)  + kpaw
463   lpaw=jab(lpaw)
      mar = lar
464   mar = mar + 1
      if(mar.gt.nbox) go to 465
      if(jcon(mar).eq.2) go to 464
      iq=iab(mar)
      if(iq.gt.lpaw) go to 466
      mpaw=jdeks(lpaw)  + iq
      go to 467
466   mpaw=jdeks(iq)  + lpaw
467   mpaw= jab(mpaw)
      jcon(mar) = jcon(mar)  + 1
      kas =0
468   kas = kas + 1
      if(kas.eq.nbxa) go to 469
      if(kas.eq.kar.or.kas.eq.lar) go to 468
      if(kas.eq.mar.or.jcon(kas).eq.0) go to 468
      jcon(kas) = jcon(kas) -1
      iq=iab(kas)
      if(iq.gt.mpaw) go to 470
      npaw=jdeks(mpaw) + iq
      go to 471
470   npaw=jdeks(iq)  + mpaw
471   npaw=jab(npaw)
      las = kas
472   las = las + 1
      if(las.eq.nbox) go to 473
      if(las.eq.kar.or.las.eq.lar) go to 472
      if(las.eq.mar.or.jcon(las).eq.0) go to 472
      iq=iab(las)
      if(iq.gt.npaw) go to 474
      ipaw=jdeks(npaw)  + iq
      go to 475
474   ipaw=jdeks(iq)  + npaw
475   ipaw=jab(ipaw)
      ipaw=jdeks(ipaw)
      if(ipaw.eq.0) go to 472
      mom= ilee(ipaw)
      npox=ilee(ipaw + 1)
      if(mom.ge.las) go to 478
      if(las.ge.npox) go to 472
      mom = las
478   jcon(las) = jcon(las)  - 1
476   mom = mom + 1
      if(mom.gt.npox) go to 477
      if(mom.eq.kar.or.mom.eq.lar) go to 476
      if(mom.eq.mar.or.jcon(mom).eq.0) go to 476
      jcon(mom) = jcon(mom) - 1
      go to 33
477   jcon(las) = jcon(las) + 1
      go to 472
403   jcon(mom) = jcon(mom) + 1
      go to 476
473   jcon(kas) = jcon(kas) + 1
      go to 468
469   jcon(mar) = jcon(mar) - 1
      go to 464
465   jcon(lar) = jcon(lar) - 1
      go to 460
461   jcon(kar) = jcon(kar) - 1
      go to 456
582   if(mxex.eq.3) return
      nbd= 5
      if(jpaw.ne.1) go to 479
      kar=0
480   kar=kar + 1
      if(kar.eq.nbox) go to 479
      if(jcon(kar).ne.2) go to 480
      jcon(kar) = 0
      lar =kar
481   lar=lar + 1
      if(lar.gt.nbox) go to 482
      if(jcon(lar).ne.2) go to 481
      jcon(lar) = 0
      kas =0
483   kas = kas + 1
      if(kas.eq.nbox) go to 484
      if(kas.eq.kar.or.kas.eq.lar) go to 483
      if(jcon(kas).ne.0) go to 483
      jcon(kas) = 2
      las = kas
485   las = las + 1
      if(las.gt.nbox) go to 486
      if(las.eq.kar.or.las.eq.lar) go to 485
      if(jcon(las).ne.0)  go to 485
      jcon(las) = 2
      go to 33
486   jcon(kas) = nzero
      go to 483
404   jcon(las) = nzero
      go to 485
484   jcon(lar) = 2
      go to 481
482   jcon(kar) = 2
      go to 480
479   nbd= 6
      kar=0
487   kar=kar + 1
      if(kar.eq.nbox) go to 488
      if(jcon(kar).eq.0) go to 487
      iq=iab(kar)
      if(iq.gt.jpaw) go to 489
      kpaw=jdeks(jpaw)  + iq
      go to 490
489   kpaw=jdeks(iq)  + jpaw
490   kpaw=jab(kpaw)
      kpaw= jdeks(kpaw)
      if(kpaw.eq.0) go to 487
      lar=ilee(kpaw)
      npox= ilee(kpaw + 1)
      if(lar.ge.kar) go to 491
      if(kar.ge.npox)  go to 487
      lar =kar
491   jcon(kar) = jcon(kar) - 1
492   lar =lar + 1
      if(lar.gt.npox) go to 493
      if(jcon(lar).eq.0) go to 492
      jcon(lar) = jcon(lar) - 1
      mar = nzero
494   mar = mar + 1
      if(mar.gt.nbox) go to 495
      if(mar.eq.kar.or.mar.eq.lar) go to 494
      if(jcon(mar).ne.2) go to 494
      jcon(mar) = nzero
      kas = nzero
496   kas = kas + 1
      if(kas.eq.nbox) go to 497
      if(kas.eq.kar.or.kas.eq.lar) go to 496
      if(kas.eq.mar.or.jcon(kas).ne.0) go to 496
      jcon(kas) = 2
      las = kas
498   las=las + 1
      if(las.gt.nbox) go to 499
      if(las.eq.kar.or.las.eq.lar) go to 498
      if(las.eq.mar.or.jcon(las).ne.0) go to 498
      jcon(las) =2
      go to 33
499   jcon(kas) = nzero
      go to 496
405   jcon(las) = nzero
      go to 498
497   jcon(mar) = 2
      go to 494
495   jcon(lar) =jcon(lar)  + 1
      go to 492
493   jcon(kar) =jcon(kar) + 1
      go to 487
488   nbd=7
      kar = nzero
532   kar = kar + 1
      if(kar.eq.nbox) go to 533
      if(jcon(kar).eq.2) go to 532
      iq=iab(kar)
      if(iq.gt.jpaw) go to 534
      kpaw=jdeks(jpaw)  + iq
      go to 535
534   kpaw=jdeks(iq)  + jpaw
535   kpaw=jab(kpaw)
      kpaw= jdeks(kpaw)
      if(kpaw.eq.0) go to 532
      lar=ilee(kpaw)
      npox= ilee(kpaw + 1)
      if (lar.ge.kar) go to 536
      if (kar.ge.npox) go to 532
      lar = kar
536   jcon(kar) = jcon(kar) + 1
537   lar = lar + 1
      if(lar.gt.npox) go to 538
      if(jcon(lar).eq.2) go to 537
      jcon(lar) =jcon(lar) + 1
      mar = nzero
539   mar = mar + 1
      if(mar.gt.nbox) go to 540
      if(mar.eq.kar.or.mar.eq.lar) go to 539
      if(jcon(mar).ne.0) go to 539
      jcon(mar) =2
      kas = nzero
541   kas = kas + 1
      if(kas.eq.nbox) go to 542
      if(kas.eq.kar.or.kas.eq.lar) go to 541
      if(kas.eq.mar.or.jcon(kas).ne.2) go to 541
      jcon(kas) = nzero
      las = kas
543   las = las + 1
      if(las.gt.nbox) go to 544
      if(las.eq.kar.or.las.eq.lar) go to 543
      if (las.eq.mar.or.jcon(las).ne.2) go to 543
      jcon(las) = nzero
      go to 33
544   jcon(kas) = 2
      go to 541
406   jcon(las) = 2
      go to 543
542   jcon(mar) = nzero
      go to 539
540   jcon(lar) =jcon(lar) - 1
      go to 537
538   jcon(kar) =jcon(kar) -1
      go to 532
533   nbd = 8
      kar = nzero
545   kar = kar +1
      if (kar.eq.nbox) go to 546
      if (jcon(kar).eq.0) go to 545
      jcon(kar) =jcon (kar) - 1
      iq= iab(kar)
      if(iq.gt.jpaw) go to 800
      kpaw=jdeks(jpaw) + iq
      go to 547
800   kpaw=jdeks(iq)  + jpaw
547   kpaw=jab(kpaw)
      lar = kar
548   lar = lar + 1
      if(lar.gt.nbox) go to 549
      if(jcon(lar).eq.0 ) go to 548
      jcon(lar) =jcon(lar)  - 1
      iq = iab(lar)
      if(iq.gt.kpaw) go to 550
      lpaw=jdeks(kpaw)  + iq
      go to 551
550   lpaw=jdeks(iq)  + kpaw
551   lpaw= jab(lpaw)
      mar = nzero
552   mar = mar  + 1
      if(mar.eq.nbox)  go to 553
      if (mar.eq.kar.or.mar.eq.lar) go to 552
      if(jcon(mar).eq.2)  go to 552
      iq=iab(mar)
      if(iq.gt.lpaw) go to 554
      mpaw= jdeks(lpaw)   + iq
      go to 555
554   mpaw= jdeks(iq)   + lpaw
555   mpaw=jab(mpaw)
      mpaw=jdeks(mpaw)
      if(mpaw.eq.0)  go to 552
      kas = ilee(mpaw)
      npox= ilee(mpaw + 1)
      if (kas.ge.mar)   go to 556
      if (mar.ge.npox)   go to 552
      kas = mar
 556  jcon(mar)  = jcon(mar)  + 1
557   kas  = kas   + 1
      if(kas.gt.npox)   go to 558
      if(kas.eq.kar.or.kas.eq.lar)   go to 557
      if (jcon(kas).eq.2)   go to 557
      jcon(kas)  =jcon(kas)   + 1
      las =0
559   las  = las  + 1
      if(las.gt.nbox)  go to 560
      if(las.eq.kar.or.las.eq.lar)  go to 559
      if(las.eq.mar.or.las.eq.kas)   go to 559
      if (jcon(las).ne.2)   go to 559
      jcon(las)  =0
      mom = 0
561   mom = mom  + 1
      if(mom.gt.nbox)   go to 562
      if(mom.eq.kar.or.mom.eq.lar)   go to 561
      if(mom.eq.mar.or.mom.eq.kas)   go to 561
      if(mom.eq.las.or.jcon(mom).ne.0)   go to 561
      jcon(mom) = 2
      go to 33
562   jcon(las)   =2
      go to 559
407   jcon(mom)   = 0
      go to 561
560   jcon(kas)  =jcon(kas)  - 1
      go to 557
558   jcon(mar)  =jcon(mar)   - 1
      go to 552
553   jcon(lar)   =jcon(lar)  + 1
      go to 548
549   jcon(kar)   =jcon(kar)  + 1
      go to 545
546   nbd= 9
      kar =0
801   kar = kar  + 1
      if(kar.eq.nbxb)  go to 563
      if(jcon(kar).eq.0 )   go to 801
      jcon(kar) = jcon(kar)    - 1
      iq= iab(kar)
      if(iq.gt.jpaw)   go to 564
      kpaw=jdeks(jpaw)  + iq
      go to 565
564   kpaw= jdeks(iq)   + jpaw
565   kpaw=jab(kpaw)
      lar = kar
566   lar = lar + 1
      if (lar.eq.nbxa)  go to 567
      if (jcon(lar).eq.0)  go to 566
      jcon(lar)  = jcon(lar)   - 1
      iq=iab(lar)
      if(iq.gt.kpaw)  go to 568
      lpaw=jdeks(kpaw)  + iq
      go to 569
568   lpaw=jdeks(iq)   + kpaw
569   lpaw=jab(lpaw)
      mar = lar
570   mar = mar + 1
      if(mar.eq.nbox)  go to 571
      if (jcon(mar).eq. 0) go to 570
      iq=iab(mar)
      if(iq.gt.lpaw)  go to 572
      mpaw=jdeks(lpaw)   + iq
      go to 573
572   mpaw= jdeks(iq)   + lpaw
573   mpaw= jab(mpaw)
      mpaw=jdeks(mpaw)
      if(mpaw.eq.0)  go to 570
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      if(kas.ge.mar)   go to 574
      if(mar.ge.npox)   go to 570
      kas = mar
574   jcon(mar)  = jcon(mar)  - 1
575   kas = kas + 1
      if(kas.gt.npox) go to 576
      if(jcon(kas).eq.0)   go to 575
      jcon(kas)   = jcon(kas)  - 1
      las =0
577   las = las + 1
      if(las.eq.nbox)  go to 578
      if(las.eq.kar.or.las.eq.lar)   go to 577
      if (las.eq.mar.or.las.eq.kas)   go to 577
      if (jcon(las).ne.0) go to 577
      jcon(las)  = 2
      mom= las
579   mom=mom + 1
      if(mom.gt.nbox)  go to 580
      if(mom.eq.kar.or.mom.eq.lar)   go to 579
      if (mom.eq.mar.or.mom.eq.kas)  go to 579
      if (jcon(mom).ne.0)  go to 579
      jcon(mom)   = 2
      go to 33
580   jcon(las)   =0
      go to 577
408   jcon(mom)   =0
      go to 579
578   jcon(kas) = jcon(kas)  + 1
      go to 575
576   jcon(mar)   =jcon(mar)  + 1
      go to 570
571   jcon(lar)  = jcon(lar)  + 1
      go to 566
567   jcon(kar) = jcon(kar)  + 1
      go to 801
563   nbd=10
      kar =0
583   kar = kar + 1
      if(kar.eq.nbxb)  go to 584
      if(jcon(kar).eq.2)  go to 583
      jcon(kar) = jcon(kar)  + 1
      iq= iab(kar)
      if(iq.gt.jpaw)  go to 585
      kpaw=jdeks(jpaw)  + iq
      go to 586
585   kpaw=jdeks(iq) + jpaw
586   kpaw= jab(kpaw)
      lar = kar
587   lar = lar +1
      if(lar.eq.nbxa)  go to 588
      if(jcon(lar).eq.2)  go to 587
      jcon(lar) = jcon(lar)  + 1
      iq= iab(lar)
      if (iq.gt.kpaw)  go to 589
      lpaw=jdeks(kpaw)  + iq
      go to 590
589   lpaw=jdeks(iq)  + kpaw
590   lpaw=jab(lpaw)
      mar = lar
591   mar = mar + 1
      if(mar.eq.nbox)  go to 592
      if(jcon(mar).eq.2)  go to 591
      iq= iab(mar)
      if(iq.gt.lpaw)  go to 593
      mpaw=jdeks(lpaw)  + iq
      go to 594
593   mpaw= jdeks( iq)  + lpaw
594   mpaw= jab(mpaw)
      mpaw=jdeks(mpaw)
      if(mpaw.eq.0)    go to 591
      kas = ilee(mpaw)
      npox=ilee(mpaw + 1)
      if(kas.ge.mar)   go to 595
      if(mar.ge.npox)  go to 591
      kas = mar
595   jcon(mar) =jcon(mar) + 1
596   kas = kas + 1
      if(kas.gt.npox)  go to 597
      if(jcon(kas).eq.2)  go to 596
      jcon(kas) = jcon(kas) + 1
      las =0
598   las =las + 1
      if(las.eq.nbox)  go to 599
      if(las.eq.kar.or.las.eq.lar)  go to 598
      if(las.eq.mar.or.las.eq.kas)  go to 598
      if(jcon(las).ne.2)  go to 598
      jcon(las) = 0
       mom = las
600   mom = mom + 1
      if(mom.gt.nbox)  go to 601
      if(mom.eq.kar.or.mom.eq.lar)  go to 600
      if(mom.eq.mar.or.mom.eq.kas)   go to 600
      if(jcon(mom).ne.2)   go to 600
      jcon(mom) = 0
      go to 33
601   jcon(las) = 2
      go to 598
409   jcon(mom) = 2
      go to 600
599   jcon(kas) = jcon(kas) - 1
      go to 596
597   jcon(mar) = jcon(mar)  - 1
      go to 591
592   jcon(lar) =jcon(lar) - 1
      go to 587
588   jcon(kar) = jcon(kar) - 1
      go to 583
584   nbd=11
      kar=0
602   kar = kar  + 1
      if(kar.eq.nbxb)  go to 603
      if(jcon(kar).eq.0)  go to 602
      jcon(kar) = jcon(kar) - 1
      iq=iab(kar)
      if(iq.gt.jpaw) go to 604
      kpaw=jdeks(jpaw) + iq
      go to 605
604   kpaw=jdeks(iq)  + jpaw
605   kpaw=jab(kpaw)
      lar = kar
606   lar = lar + 1
      if(lar.eq.nbxa)  go to 607
      if(jcon(lar).eq.0)  go to 606
      jcon(lar) = jcon(lar) - 1
      iq=iab(lar)
      if(iq.gt.kpaw)  go to 608
      lpaw=jdeks(kpaw)  + iq
      go to 609
608   lpaw=jdeks(iq)   + kpaw
609   lpaw=jab(lpaw)
      mar = lar
610   mar = mar + 1
      if(mar.eq.nbox)  go to 611
      if(jcon(mar).eq.0)  go to 610
      jcon(mar) = jcon(mar)  - 1
      iq= iab(mar)
      if(iq.gt.lpaw)  go to 612
      mpaw=jdeks(lpaw)   + iq
      go to 613
612   mpaw= jdeks(iq)  + lpaw
613   mpaw= jab(mpaw)
      kas = mar
614   kas = kas + 1
      if(kas.gt.nbox)  go to 615
      if(jcon(kas).eq.0)   go to 614
      jcon(kas)= jcon(kas) - 1
      iq= iab(kas)
      if(iq.gt.mpaw)  go to 616
      npaw= jdeks(mpaw)  + iq
      go to 617
616   npaw= jdeks(iq)  + mpaw
617   npaw= jab(npaw)
      las =0
618   las = las + 1
      if(las.eq.nbox) go to 619
      if( las.eq.kar.or.las.eq.lar)   go to 618
      if (las.eq.mar.or.las.eq.kas)   go to 618
      if(jcon(las).eq.2)   go to 618
      iq= iab(las)
      if(iq.gt.npaw) go to 620
      ipaw=jdeks(npaw)  + iq
      go to 621
620   ipaw= jdeks(iq)  + npaw
621   ipaw= jab(ipaw)
      ipaw= jdeks(ipaw)
      if(ipaw.eq.0)   go to 618
      mas= ilee(ipaw)
      npox= ilee(ipaw + 1)
      if(mas.ge.las)   go to 622
      if(las.ge.npox)   go to 618
      mas = las
622   jcon(las) = jcon(las)  + 1
623   mas = mas + 1
      if(mas.gt.npox)   go to 624
      if(jcon(mas).eq.2)  go to 623
      if(mas.eq.kar.or.mas.eq.lar)  go to 623
      if(mas.eq.mar.or.mas.eq.kas)   go to 623
      jcon(mas) = jcon(mas) + 1
      mom =0
625   mom = mom + 1
      if(mom.gt.nbox)   go to 626
      if(jcon(mom).ne.0) go to 625
      if(mom.eq.kar.or.mom.eq.lar) go to 625
      if(mom.eq.mar.or.mom.eq.kas) go to 625
      if(mom.eq.las.or.mom.eq.mas) go to 625
      jcon(mom) =2
      go to 33
626   jcon(mas) = jcon(mas)  - 1
      go to 623
410   jcon(mom) = 0
      go to 625
624   jcon(las) = jcon(las) - 1
      go to 618
619   jcon(kas) = jcon(kas) + 1
      go to 614
615   jcon(mar) = jcon(mar)  + 1
      go to 610
611   jcon(lar) = jcon(lar) + 1
      go to 606
607   jcon(kar) =jcon(kar) + 1
      go to 602
603   nbd=12
      kar=0
627   kar = kar + 1
      if(kar.eq.nbxb) go to 628
      if(jcon(kar).eq.2)  go to 627
      jcon(kar) = jcon(kar) + 1
      iq=iab(kar)
      if(iq.gt.jpaw) go to 629
      kpaw=jdeks(jpaw)  + iq
      go to 630
629   kpaw= jdeks(iq)  + jpaw
630   kpaw= jab(kpaw)
      lar = kar
631   lar = lar + 1
      if(lar.eq.nbxa) go to 632
      if(jcon(lar).eq.2) go to 631
      jcon(lar) = jcon(lar) + 1
      iq= iab(lar)
      if(iq.gt.kpaw) go to 633
      lpaw= jdeks(kpaw) + iq
      go to 634
633   lpaw=jdeks(iq)  + kpaw
634   lpaw= jab(lpaw)
      mar = lar
635   mar = mar + 1
      if(mar.eq.nbox)  go to 636
      if(jcon(mar).eq.2)   go to 635
      jcon(mar) = jcon(mar)  + 1
      iq=iab(mar)
      if(iq.gt.lpaw)  go to 637
      mpaw=jdeks(lpaw)  + iq
      go to 638
637   mpaw= jdeks(iq)   + lpaw
638   mpaw= jab(mpaw)
      kas = mar
639   kas = kas + 1
      if(kas.gt.nbox)  go to 640
      if(jcon(kas).eq.2)   go to 639
      jcon(kas) = jcon(kas) + 1
      iq= iab(kas)
      if(iq.gt.mpaw)   go to 641
      npaw= jdeks(mpaw)  + iq
      go to 642
641   npaw= jdeks(iq)   + mpaw
642   npaw= jab(npaw)
      las=0
643   las = las + 1
      if(las.eq.nbox) go to 644
      if(las.eq.kar.or.las.eq.lar)   go to 643
      if(las.eq.mar.or.las.eq.kas)   go to 643
      if(jcon(las).eq.0)   go to 643
      iq= iab(las)
      if(iq.gt.npaw)   go to 645
      ipaw= jdeks(npaw)  + iq
      go to 646
645   ipaw=jdeks(iq)   + npaw
646   ipaw=jab(ipaw)
      ipaw= jdeks(ipaw)
      if(ipaw.eq.0)  go to 643
      mas = ilee(ipaw)
      npox= ilee(ipaw + 1)
      if(mas.ge.las)   go to 647
      if(las.ge.npox)   go to 643
      mas = las
647   jcon(las)=jcon(las)  - 1
648   mas = mas + 1
      if(mas.gt.npox)  go to 649
      if(jcon(mas).eq.0)   go to 648
      if(mas.eq.kar.or.mas.eq.lar)   go to 648
      if(mas.eq.mar.or.mas.eq.kas)   go to 648
      jcon(mas) = jcon(mas)   - 1
      mom = 0
650   mom = mom + 1
      if(mom.gt.nbox)   go to 651
      if(jcon(mom).ne.2)   go to 650
      if(mom.eq.kar.or.mom.eq.lar)   go to 650
      if(mom.eq.mar.or.mom.eq.kas)   go to 650
      if (mom.eq.las.or.mom.eq.mas)   go to 650
      jcon(mom)   = 0
      go to 33
651   jcon(mas)   = jcon(mas)   + 1
      go to 648
411   jcon(mom)   = 2
      go to 650
649   jcon(las) = jcon(las) + 1
      go to 643
644   jcon(kas) = jcon(kas)  - 1
      go to 639
640   jcon(mar) = jcon(mar)  - 1
      go to 635
636   jcon(lar) = jcon(lar)   - 1
      go to 631
632   jcon(kar) = jcon(kar)   -1
      go to 627
628   nbd=13
      kar = 0
652   kar = kar + 1
      if(kar.eq.nbxb)   return
      if(jcon(kar).eq.0)  go to 652
      jcon(kar) = jcon(kar) - 1
      iq= iab(kar)
      if(iq.gt.jpaw)  go to 654
      kpaw=jdeks(jpaw)  + iq
      go to 655
654   kpaw= jdeks(iq)  + jpaw
655   kpaw= jab(kpaw)
      lar = kar
656   lar = lar + 1
      if(lar.eq.nbxa)  go to 657
      if(jcon(lar).eq.0)   go to 656
      jcon(lar) = jcon(lar) -  1
      iq=iab(lar)
      if(iq.gt.kpaw)  go to 658
      lpaw= jdeks(kpaw)  +  iq
      go to 659
658   lpaw= jdeks(iq)   + kpaw
659   lpaw= jab(lpaw)
      mar = lar
660   mar = mar + 1
      if(mar.eq.nbox)   go to 661
       if(jcon(mar).eq.0)  go to 660
      jcon(mar) = jcon(mar)   -   1
      iq= iab(mar)
      if(iq.gt.lpaw)   go to 662
      mpaw= jdeks(lpaw)  + iq
      go to 663
662   mpaw= jdeks(iq)  + lpaw
663   mpaw= jab(mpaw)
      kas = mar
664   kas = kas + 1
      if(kas.gt.nbox)   go to 665
      if(jcon(kas).eq.0)   go to 664
      jcon(kas) = jcon(kas)   -1
      iq= iab(kas)
      if(iq.gt.mpaw)   go to 666
      npaw= jdeks(mpaw)   + iq
      go to 667
666   npaw= jdeks(iq)   + mpaw
667   npaw= jab(npaw)
      las =0
668   las = las + 1
      if(las.eq.nbxb)  go to 669
      if(las.eq.kar.or.las.eq.lar)   go to 668
      if(las.eq.mar.or.las.eq.kas)   go to 668
      if(jcon(las).eq.2)   go to 668
      jcon(las) = jcon(las)   +1
      iq= iab(las)
      if(iq.gt.npaw)   go to 802
      ipaw= jdeks(npaw)  + iq
      go to 671
802   ipaw= jdeks(iq)   + npaw
671   ipaw= jab(ipaw)
      mas = las
672   mas = mas + 1
      if(mas.eq.nbxa)   go to 673
      if(jcon(mas).eq.2)  go to 672
      if(mas.eq.kar.or.mas.eq.lar)   go to 672
      if(mas.eq.mar.or.mas.eq.kas)   go to 672
      jcon(mas) = jcon(mas)  +  1
      iq= iab(mas)
      if(iq.gt.ipaw)   go to 674
      ipop= jdeks( ipaw)  + iq
      go to 675
674   ipop= jdeks(iq)   + ipaw
675   ipop= jab(ipop)
      mom = mas
676   mom = mom + 1
      if(mom.eq.nbox)   go to 677
      if(jcon(mom).eq.2)   go to 676
      if(mom.eq.kar.or.mom.eq.lar)   go to 676
      if(mom.eq.mar.or.mom.eq.kas) go to 676
      iq= iab(mom)
      if(iq.gt.ipop)  go to 678
      idad= jdeks(ipop)  + iq
      go to 679
678   idad= jdeks(iq)   + ipop
679   idad= jab(idad)
      idad= jdeks(idad)
      if(idad.eq.0) go to 676
      mel= ilee(idad)
      npox= ilee(idad + 1)
      if(mel.ge.mom) go to 680
      if(mom.ge.npox)  go to 676
      mel= mom
680   jcon(mom) = jcon(mom)  +  1
681   mel= mel + 1
      if(mel.gt.npox)  go to 682
      if(jcon(mel).eq.2)   go to 681
      if(mel.eq.kar.or.mel.eq.lar)   go to 681
      if(mel.eq.mar.or.mel.eq.kas)   go to 681
      jcon(mel) = jcon(mel) + 1
      go to 33
682   jcon(mom) = jcon(mom)  -  1
      go to 676
412   jcon(mel) = jcon(mel) -  1
      go to 681
677   jcon(mas) = jcon(mas)     -1
      go to 672
673   jcon(las) = jcon(las)   -   1
      go to 668
669   jcon(kas) = jcon(kas)   +   1
      go to 664
665   jcon(mar) = jcon(mar)   +   1
      go to 660
661   jcon(lar) = jcon(lar)   +   1
      go to 656
657   jcon (kar) = jcon(kar)   +   1
      go to 652
33    if (klx.eq.1) go to 25
      lly=-nshl
      do 34 i=1,kly
      mx=0
      lly=lly+nshl
      np=nop(i)
      ja=ndub(i)
      llz=lly
      if (np.eq.0) go to 173
      do 172 j=1,np
      llz=llz+1
      jd=jkon(llz)
      if (jcon(jd).gt.0) go to 172
      mx=mx+1
      if (mx.gt.mxex) go to 34
 172  continue
      if (ja.eq.0) go to 36
 173  do 174 j=1,ja
      llz=llz+1
      jd=jkon(llz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 174
      mx=mx+jd
      if (mx.gt.mxex) go to 34
 174  continue
      go to 36
34    continue
25    nod=0
      ncl=0
      do 76 i=1,nbox
      jc=jcon(i)
      if (jc.eq.0) go to 76
      if (jc.eq.1) go to 77
      ncl=ncl+1
      iduf(ncl)=i
      go to 76
77    nod=nod+1
      icuf(nod)=i
76    continue
      nmns=nod/2+1-ispin
      if(nmns.le.0) goto 36
      ig=ig+1
      kcon(ig)=nmns
      nconf(nmns)=nconf(nmns)+1
      if (ig.lt.ibl) go to 85
      ig=0
      write(mdisk)kcon
 85   if (nod.eq.0) go to 80
      do 86 i=1,nod
      ig=ig+1
      kcon(ig)=icuf(i)
      if (ig.lt.ibl) go to 86
      ig=0
      write(mdisk) kcon
86    continue
      if (ncl.eq.0) go to 36
80    do 81 i=1,ncl
      ig=ig+1
      kcon(ig)=iduf(i)
      if (ig.lt.ibl) go to 81
      ig=0
      write(mdisk) kcon
81    continue
36    go to (400,401,402,403,404,405,406,407,408,409,410,411,
     *412), nbd
      return
      end
*************************************************************************
_IF(notused)
      subroutine check0 (nr,nd,itest,jtest,icheck)
      implicit REAL (a-h,o-z)
      integer itest,jtest
      dimension itest(*),jtest(*)
INCLUDE(common/iofile)
      nall=nd+nr
c     write(iwr,1021)nd,nr,nall
c1021 format(1x,'check0 was called '/3i3)
c     write(iwr,1122)(itest(kaus),kaus=1,nall)
c1122 format(1x,'itest ',10i3)
c     write(iwr,1121)(jtest(kaus),kaus=1,nall)
c1121 format(1x,'jtest ',10i3)
c possible interaction are a/b or abb/baa. the diagonal case
c is not considered
c if a single excitation is found icheck contains the integral
c address
      icheck=0
c check of the closed shells
      it=nr+1
      jt=nr+1
c     write(iwr,2) it,jt
    2 format(1x,'check of the closed shells starting with ',2i3)
      im=0
      jm=0
      nd1=nd-1
      do 100 i=1,nd1
c       write(iwr,10) it,itest(it),jt,jtest(jt)
c10     format(1x,'check doubly:it,itest(it),jt,jtest(jt)',4i3)
        if(itest(it).ne.jtest(jt)) then
          if(itest(it+1).eq.jtest(jt)) then
            if(im.ne.0) return
            im=itest(it)
            it=it+1
            go to 110
          endif
          if(itest(it).eq.jtest(jt+1)) then
            if(jm.ne.0) return
            jm=jtest(jt)
            jt=jt+1
            go to 110
          endif
          if(itest(it+1).eq.jtest(jt+1)) then
            if(jm.ne.0.or.im.ne.0) return
            jm=jtest(jt)
            im=itest(it)
            it=it+1
            jt=jt+1
            go to 110
          endif
c if he reaches that point at least two mo nonidentical
          return
        endif
  110   continue
        it=it+1
        jt=jt+1
  100 continue
      if(jm.eq.0.or.im.eq.0)then
c we have to look for that one which is zero
        if(im-jm) 500,600,700
c im.eq.0  jm.ne.0
  500   im=itest(it)
        go to 1000
c im.eq.0 and jm.eq.0
  600   if(itest(it).ne.jtest(jt)) then
          im=itest(it)
          jm=jtest(jt)
        endif
        go to 1000
c im.ne.0 and jm.eq.0
  700   jm=jtest(jt)
        go to 1000
      endif
 1000 continue
c     write(iwr,11)jm,im
   11 format(1x,'doubly finished   jm,im ',2i3)
      if(im.eq.0) then
c     write(iwr,3)
    3 format(1x,' closed shells are equal, a/b only ')
c here all closed shells are equal --->   a/b interaction only
        it=1
        jt=1
c you can run up to # open shells because the closed shells are
c identical
c     write(iwr,5) it,jt
    5 format(1x,'check of the open shells starting with ',2i3)
        do 200 i=1,nr
c         write(iwr,6) it,itest(it),jt,jtest(jt)
c   6     format(1x,'it,itest,jt,jtest ',i2,i4,i2,i4)
          if(itest(it).ne.jtest(jt)) then
            if(itest(it+1).eq.jtest(jt)) then
c             write(iwr,7)im
c7            format(1x,'itest(it+1)=jtest(jt)  im',i3)
              if(im.ne.0) return
              im=itest(it)
              it=it+1
              go to 210
            endif
            if(itest(it).eq.jtest(jt+1)) then
c             write(iwr,8) jm
c8            format(1x,'itest(it)=jtest(jt+1)   jm',i3)
              if(jm.ne.0) return
              jm=jtest(jt)
              jt=jt+1
              go to 210
            endif
            if(itest(it+1).eq.jtest(jt+1)) then
c             write(iwr,9) jm,im
c9            format(1x,'itest(it+1)=jtest(jt+1) jm,im',2i3)
              if(jm.ne.0.or.im.ne.0) return
              jm=jtest(jt)
              im=itest(it)
              it=it+1
              jt=jt+1
              go to 210
            endif
c if he reaches that point at least two mo nonidentical
            return
          endif
  210   continue
        it=it+1
        jt=jt+1
  200   continue
c if the program reaches this point an a/b interaction is found
c calculation of the integraladdress
        ma=max(im,jm)
        mi=min(im,jm)
        ma=ma-1
        ma=ma*(ma+1)/2
        icheck=ma+mi
        return
      else
c one missmatch in closed shells only abb/baa possible
c     write(iwr,4)
    4 format(1x,'one missmatch in closed shells abb/baa only')
        it=1
        jt=1
        nrend=nr-1
c running only up to nrend because there must be at least one
c missmatch
        im1=0
        do 300 i=1,nrend
c         write(iwr,15) it,jt
c15       format(1x,'it,jt ',2i3)
          if(itest(it).ne.jtest(jt)) then
            im1=1
            if(itest(it+1).eq.jtest(jt)) then
c             write(iwr,13)itest(it+1),it,jtest(jt),jt,im,jm
c13           format(1x,'itest(it+1)=jtest(jt) ',i3,i2,2(i4,i2))
              if(itest(it).ne.jm) return
              it=it+1
              go to 310
            endif
            if(itest(it).eq.jtest(jt+1)) then
c             write(iwr,16)itest(it),it,jtest(jt+1),jt,im,jm
c16           format(1x,'itest(it)=jtest(jt+1) ',i3,i2,2(i4,i2))
              if(jtest(jt).ne.im) return
              jt=jt+1
              go to 310
            endif
            return
          endif
  310     continue
          it=it+1
          jt=jt+1
  300   continue
c if im1=0 the last two open shells must be checked
        if(im1.eq.0) then
          if(itest(it).ne.jm.or.jtest(jt).ne.im)return
        endif
c if programm reachs this point an interaction abb/baa is found
c calculation of the integral address
        ma=max(im,jm)
        mi=min(im,jm)
        ma=ma-1
        ma=ma*(ma+1)/2
        icheck=ma+mi
      endif
      return
      end
*************************************************************************
      subroutine check1 (nr,nd,nr1,nd1,itest,jtest,icheck)
      implicit REAL (a-h,o-z)
      integer itest,jtest
      dimension itest(*),jtest(*),jm(2)
INCLUDE(common/iofile)
c     write(iwr,1)
c   1 format(1x,'check1 was called')
      nall1=nr+nd
c     write(iwr,3)nr,nd,nall1,(itest(kaus),kaus=1,nall1)
c   3 format(1x,'sk=i:nr,nd,nall1,itest ',3i3,2x,15i3)
      nall2=nr1+nd1
c     write(iwr,4)nr1,nd1,nall2,(jtest(kaus),kaus=1,nall2)
c   4 format(1x,'sk=i+1:nr1,nd1,nall2,jtest ',3i3,2x,15i3)
      icheck=0
      nall=nd+nr
      jm(1)=0
      jm(2)=0
c check of open shells of sk i, they must be contained in i+1
      it=0
      jt=0
      nix=0
      do 10 i=1,nr
        it=it+1
        jt=jt+1
 25     if(itest(it)-jtest(jt)) 20,10,30
c 20:if itest(it)<jtest(jt) it can not appear afterwards -- return
c if itest(it)>jtest(jt) it may appear a little later but jtest(jt)
c must be stored
 30     if(nix.eq.2) return
        nix=nix+1
        jm(nix)=jtest(jt)
        jt=jt+1
        go to 25
 10   continue
      jrest=nr1-jt
c     write(iwr,100) jrest,jt,nr1,jm,nix
c100  format(1x,'first part singles from jrest,jt,nr1,jm,nix ',6i3)
      if(jrest.eq.0) go to 40
      do 50 i=1,jrest
        if(nix.eq.2) return
        nix=nix+1
        jt=jt+1
 50   jm(nix)=jtest(jt)
 40   continue
c     write(iwr,2)jm
c 2   format(1x,'singles from  jm(1,2) ',2i3)
c test wether the doubly occupied in i+1 are contained in i
      niy=0
      do 60 i=1,nd1
        it=it+1
        jt=jt+1
        if(jtest(jt)-itest(it))20,60,70
c 20: jtest(jt)<itest(it) cannot appear afterwards
c 70: jtest(jt)<itest(it) may appear afterwards but must be equal
c to one of the missmatching single occupied shells
c 60: jtest(jt)=itest(it) peace, joy and teacakes
 70     if(niy.eq.1) return
        im=itest(it)
        if(im.ne.jm(1).and.im.ne.jm(2)) return
        niy=niy+1
        it=it+1
        if(jtest(jt)-itest(it))20,60,20
 60   continue
c if missmatching closed shell is the last
c     write(iwr,101) im,it,nd,nr,nall
c101  format(1x,'closed shells from im,it,nd,nr,nall ',5i3)
      if(it.lt.nall) im=itest(it+1)
      if(im.eq.jm(1)) then
        jim=jm(2)
      else if(im.eq.jm(2)) then
        jim=jm(1)
      else
        return
      endif
      ma=max(im,jim)
      mi=min(im,jim)
      ma=(ma-1)*ma/2
      icheck=ma+mi
 20   return
      end
************************************************************************
_ENDIF
************************************************************************
      subroutine coolez(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +                  sx,lsx,
     +                  pey,acoul,aexc,ideks,ndeks,
     +                  peyo,acoulo,aexco,nf31o,
     +                  sdelta,seng,lbuffs,
     +                  vect,wect,coff,nvect,
     +                  iot,niott,
     +                  ot,nt444,jcon1,nd8,jkan,njkan,
     +                  jtest)
c
c  --- computes the interaction with the geyser1-space and determines
c      the geysez1-space that will be stored on ft99.
c
c  mconf(5) : number of mains per super category
c  nconf(5) : number of generated configurations per super category
c
      implicit REAL (a-h,o-z)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
      REAL hy,hp,a,w,e,f,sx
      REAL sdelta,seng,vect,wect
      integer nt64, mxref2, ndimf, lsx, lbuffs, jcon1
      integer jtest
      integer jkan, njkan, nvect
      dimension hy(nt64),hp(nt64),w(mxref2),a(nt64)
      dimension e(ndimf),f(ndimf)
      dimension sx(lsx)
      dimension sdelta(5,lbuffs), seng(lbuffs)
      dimension vect(nvect),wect(nvect),coff(nvect)
      dimension jtest(nvect)
      dimension jcon1(nd8)
      dimension jkan(njkan)
c
      REAL pey,acoul,aexc
      integer ideks,ndeks
      dimension pey(ndeks),acoul(ndeks),aexc(ndeks)
      dimension ideks(ndeks)
c
      REAL peyo
      REAL acoulo,aexco
      integer nf31o
      dimension peyo(nf31o),acoulo(nf31o),aexco(nf31o)
c
      integer iot
      dimension iot(niott)
c
      integer ot
      dimension ot(nt444)
c
      REAL timeb,timea,rzero
c
      integer n,imo0,icisi,mmm
      integer nteste
      parameter (rzero = 0.0d0)
c
c --- commons
c
      integer nconf0
      integer jab0, isym0, mj0, iab0, kj0, jcon0
      integer icuf0, iduf0
      integer kcon0, ilee0, ntil0, nbal, nop
      integer isw, nyzl, nj0, ijoe0, ndub0
c
      parameter (nmmo=256)
c
      common /aeng/ nconf0(ndk5),jab0(n36),isym0(maxsym),
     .              mj0(maxsym),iab0(nmmo), kj0(maxsym),
     .              jcon0(nmmo),icuf0(n12), iduf0(maxshl),
     .              kcon0(ndimh),ilee0(maxsym+1),
     .              ntil0(maxsym),nbal(maxsym+1),
     .              jkon(maxref*maxshl),nop(maxref),
     .              isw(ndk5),nyzl(maxref),nj0(maxsym),
     .              ijoe0(maxref),ndub0(maxref)
c
      common /cimo/ imo0
      common /cisi/ icisi(maxref)
c
      REAL egey1,trash,tdel
      integer negey1,nko2,nstarv
      common /cegey1/ egey1,trash,tdel,negey1,nko2,nstarv
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      integer icon
      common /cicon/  icon(ndk5)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
      integer nytl,ndub
      common /cny /   nytl(ndk5),ndub(ndk5)
      integer nod
      common /cnod/   nod(ndk5)
c
      integer lj
      common /clj2/ lj(maxsym)
c
      integer ij
      common /ccij/   ij(maxsym)
      integer kc,kd,lab
      common /ckd/    kc(ndimh),kd(ndimh),lab(nn3)
      integer nhead, mtapev, nf11
      common /ftap/   nston, mtape, mdisk,
     .                ideli, ltype, linf,  ntab,  kfile,
     .                kclab, ntype, mstvt, nf01, nf62,
     +                nhead, nf99,  mtapev, nf11
      common /tap/    ical,  ma,    nko,
     .                mxex,  nmul,  ispace, nspace, ncorci,
     +                msng,  maxci2
c --- hier wird ein teil von tap auf b transferiert nston z.b.
      integer mn, nsc
      common /b/ iwod,lg,irsf,i11,i12,i3,i4,nc,lt,lz,nnx,nd,
     +i5,i8,i9,j1,iswh,nr,mm,jblk,
     +jto,igmax,nr1,j2,j3,j4,nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +mconf(ndk5),nconf(ndk5),mh,jswh,
     +nn(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx0,md,mr,ipag
c
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      common /einf/ mhe(maxref),imain(maxref),iselct
      common /ft31in/ n
c
      REAL ptiw
_IF(notused)
      common/junk/ ptiw(5500)
_ELSE
_ENDIF
c
      integer labs,mkon,nplu
      integer kj,mj,jconr
      REAL xemp
      common /miscop/ xemp(500),mkon(ndimh),labs(ndimh),nplu(ndk5),
     +                kj(maxsym),mj(maxsym),jconr(maxshl)
c 
c
c --------------------------------------------------------------
c --- copying of commons for calling the rumplx-routines and
c --- recasting data from i*2 to i*4
*     do i=1,nitmax
*        nit2(i) = (nit(i))
*     end do
*     do i=1,ideksm
*        ideks2(i) = (ideks(i))
*     end do
c --------------------------------------------------------------
      call rewftn(linf)
      call rewftn(mtape)
      call rewftn(nhead)
      call rewftn(ntype)
c
      write(iwr,*) 'coolez : ntype =',ntype
c
      call rewftn(kclab)
c -- kfile rewind ! because it gets overwritten !
      call rewftn(kfile)
c --- write header
      write(iwr,*) '==================================================='
      write(iwr,*) 'compute the discarded energy contributions (coolez)'
      write(iwr,*) '==================================================='
      write(iwr,*)
cdebugwrite(iwr,*) ' icon >>>>>'
cdebugdo i=1,ndk5
cdebug  write(iwr,*)icon(i)
cdebugend do
c
c --- egey changes sign for a normal run !!
c
      egey1 = - 999d+6
c changings to take more single excitations into account
      do 3710 ieng=1,maxref
        mhe(ieng)  = nzero
        imain(ieng)= nzero
3710  continue
c counter in mhe array
      imhe  = 1
c end of changings
c --- copying of /tap/ to /b/ is obsolete here ! leave out
      nirm  = nmmo
      igmax = n3zt
c --- restore common /b/
      nad = nzero
      nad1 = nzero
c --- new mtape
***   j1 = 97
c --- number of mo's : ndeks - 1
      im3  = ndeks - 1
      nshl = maxshl
c --- why nsec =80 ???????
      nsec = maxref
c ------
      read(mtape) nbox,m,nmul,ispace,mconf,nytl,mxex,jkon,nko,iswh,nyzl
c
c --- read the configuration from ft stream nf99 
      call rewftn(nf99)
      read(nf99,*)
      read(nf99,*)
      read(nf99,*)
      nt0= nzero
      im = nzero
c ---------------------------------------------------------
*     jblk = nzero
*     jto  = nzero
*     mm   = nzero
c ---- for the table
      irsf = 500
      lko  = maxref
      incr = 10
      msec = nsec + 1
      nawk = maxele/2+4
      nrmax= 4
c --- reading in the geyser-space:
      nt0= nzero
      im = nzero
      nops = 5
      write(iwr,*) 'nko2 = ',nko2
      do 9002 i=1,nko2
        nt0 = im + 1
        im  = im + nshl
c ---- reading the configurations
*       read(nf99,670)np,(jkon(j),j=nt0,im)
*       read(nf99,*)np,(jkon(j),j=nt0,im)
        read(nf99,*)np
        write(iwr,*) 'np=',np
        mmm = (ma+np)/2
        read(nf99,*)(jconr(j),j=1,mmm)
        do 3766 jj=nt0,im
           jkon(jj) = jconr(jj-nt0+1)
3766    continue
c ---- the configurations get reformatted and stored:
        write(iwr,*) 'nop :',np,'csf: ',i,' :',
     .   (jkon(j),j=nt0,im),' was read'
c
        if(np.gt.nops) then
          write(iwr,*) 'np=',np,' greater than ',nops
          call aborts('too many open shells')
        end if
        nyt = 1
        if(np.eq.0) go to  9994
        nyt = isw(np)
9994    nyt = nytl(nyt)
        nyzl(i) = nyt
        if(i.eq.1) go to 9503
        i1 = i - 1
        mx =-nshl
        do 9504 j=1,i1
         mx=mx+nshl
         if(nop(j).eq.np) then
          nv = mx
          la = nt0 - 1
          do 3792 k=1,nyt
            la = la + 1
            nv = nv + 1
            if(jkon(la).ne.jkon(nv)) go to 9504
3792      continue
          write(iwr,*) 'configuration was read incorrectly !!!'
         end if
9504    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*        write(iwr,*) 'half way coolez 1'
*        do iii=1,10
*           write(iwr,*) 'ideks(',iii,')=',ideks(iii)
*        end do
*        write(iwr,*) '-------'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
9503    na = nmul + np
        if(na-2*(na/2).eq.0) go to 773
        if(np.gt.m) go to 773
c
c     I dont understand the "np.gt.nmul+3" condition below given the
c     current setting for nopmax e.g. 6 open shells when nmul=1
c     will trigger it.
c     I've tried editing it out, but causes problems with
c     "not all mains generated"
c
      if(np.gt.nmul+3) go to 773
        nop(i) = np
        ndub0(i)= (m-np)/2
        if(np.lt.2) go to 9500
        nz0 = jkon(nt0)
        if(nz0.lt.1.or.nz0.gt.nbox) go to 773
        do 3818 j=2,np
         nt0 =nt0 + 1
         nv =jkon(nt0)
         if(nv.le.nz0) go to 773
         if(nv.gt.nbox) go to 773
         nz0=nv
3818    continue
        if(np.eq.m) go to 9002
9518    mx = np + im - nshl + 1
        nv = jkon(mx)
        if(nv.lt.1.or.nv.gt.nbox) go to 773
        lb = im - nshl
        do 9506 j=1,np
         lb = lb + 1
         if(jkon(lb)-nv) 9506,773,9507
9506    continue
9507    jm = np + 2
        if(jm.gt.nyt) go to 9002
        kp = mx
        do 9508 j=jm,nyt
         kg = mx
         mx = mx + 1
         nv = jkon(mx)
         if(nv.lt.1.or.nv.gt.nbox) go to 773
         do 9509 k=kp,kg
          nz0 = jkon(k)
          if(nv.le.nz0) go to 773
9509     continue
         kg=im-nshl
         do 9510 k=1,np
           kg = kg + 1
           if(jkon(kg)-nv) 9510,773,9508
9510     continue
9508    continue
        go to 9002
9500    if(np.eq.0) go to 9511
        nv = jkon(nt0)
        if(nv.lt.1.or.nv.gt.nbox) go to 773
        go to 9518
9511    nz0 = jkon(nt0)
        if(nz0.lt.1.or.nz0.gt.nbox) go to 773
        if(nyt.eq.1) go to 9002
        kp = nt0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 3867 j=2,nyt
         kg  = nt0
         nt0 = nt0 + 1
         nv  = jkon(nt0)
         if(nv.lt.1.or.nv.gt.nbox) go to 773
         do 3867 k=kp,kg
          nz0=jkon(k)
          if(nv.le.nz0) then
            write(iwr,*) 'nv=',nv,' nz=',nz0
            call aborts('error nv/nz!')
          end if
3867    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
9002   continue
       goto 9000
 773   continue
       call aborts('error !!!')
9000   continue
*      backspace (mtape)
*      write(mtape)
*    .            nbox,m,nmul,ispace,mconf,nytl,mxex,jkon,nko,iswh,nyzl
c --- preparations from parke4
*      kg = nbox*nko2
*      do i=1,kg
*         mkon(i) = nzero
*      end do
*      kg = -maxshl
*      kb = -nbox
*      do 9514 i=1,nko2
*         write(iwr,*) '9514 :i=',i
*         write(iwr,*)
*         kg = kg + nshl
*         kb = kb + nbox
*         np = nyzl(i)
*         if (np.eq.nzero) go to 9516
*         do j=1,np
*            nv = kb+jkon(j+kg)
*            mkon(nv) =1
*         end do
*         if (np.eq.m) goto 9514
*9516      np = np + 1
*         do jj=np,nyt
*            nv = kb + jkon(j+kg)
*            mkon(nv)=2
*         end do
*9514      continue
c
c --- nko enlarged !
      nko = nko2
      write(iwr,*) 'enlarged reference space :',nko,'csfs'
c --- recopying !
      imo = imo0
cdebugwrite(iwr,*) 'coolez : iwod=',iwod
c
      ix = nzero
c
      do 3921 i=1,n
       kap=lj(i)
       if (kap.ne.nzero) then
        do 3919 j=1,kap
         ix      = ix + 1
         nir(ix) = i
         loc(ix) = j
3919    continue
       end if
3921  continue
      nnx  = nmul - 2
      ndz = nytl(1) - nnx
      do 3930 i=1,iswh
       nnx = nnx + 1
       ndz= ndz - 1
       nplu(i) = nnx
       nod(i)  = nplu(i) + i - 1
       ndub0(i) = ndz
3930  continue
_IF(notused)
*     if (ical.gt.1) go to 520
*     read (ird,2) lulu,nrootx,lsng,nprin,ipt0,isec
*     write (iwr,2) lulu,nrootx,lsng,nprin,ipt0,isec
c change to perform modified selection scheme
c --- just temporarily fixed !!!!! (hus)
_ENDIF
      iselct = nzero
      if(lsng.lt.0) then
        iselct = 1
        lsng   =-lsng
      endif
c end of change
      if (mxex.eq.1) go to 500
      lulu = nko2
      go to 501
 500  continue
      if (lulu.eq.nko2) go to 501
      ix = nzero
c  ** check here, for jtest appears unset if lulu.ne.nko2
      do 502 i=1,lulu
       jt = jtest(i)
       if (jt.eq.i) go to 503
       jx = (jt-1)*nshl
       do 3957 j=1,nshl
        jx = jx + 1
        ix = ix + 1
        jkon(ix) = jkon(jx)
3957   continue
       go to 502
 503   ix = ix + nshl
 502  continue
 501  ny = nzero
      ix = nzero
c changing to take more single excitations into account
      lsng1 = nzero
      if(.not.(mxex.lt.2.or.lsng.eq.nzero)) then
        if(lsng.gt.25) then
          lsng =lulu
          lsng1=1
         do ii=1,lsng
            imain(ii) = (ii)
         end do
        else
         write(iwr,*) ' singles taken from geysez !'
cgeysez  read(ird,7)(imain(ieng),ieng=1,lsng)
         do 3974 ii=1,lsng
            imain(ii) = icisi(ii)
3974     continue
         write(iwr,9)(imain(ieng),ieng=1,lsng)
        endif
      endif
c     write(iwr,*) lsng,lsng1
c end of changings
c
      do 505 i=1,3
        nc = mconf(i)
        if (nc.eq.0) go to 505
        nt = nytl(i)
        if (ix.lt.lulu) go to 506
        nb = (nc*nt)/iwod+1
        do 3990 j=1,nb
          read (mtape)
3990    continue
        go to 505
 506    read (mtape) mkon
        ig=nzero
        do 508 j=1,nc
          do 4002 k=1,nt
           ig = ig + 1
           jtest(k) = mkon(ig)
           if (ig.ge.iwod) then
            ig = nzero
            read (mtape) mkon
           end if
4002    continue
        kg1 = -nshl
        do 510 k=1,lulu
          kg1 = kg1 + nshl
          if (nyzl(k).ne.nt) go to 510
          l7 = kg1
          do 4011 l=1,nt
            l7 = l7 + 1
            if (jtest(l).ne.jkon(l7)) go to 510
4011      continue
          ix = ix + 1
          mn(ix) = k
          nn(ix) = j
          nsc(ix) = i
c change to include all single excitations with respect to more than
c one reference conf.
          if(lsng.eq.nzero) go to 10
          if(lsng1.eq.1) then
           mhe(ix) = ix*imo
           write(iwr,*)mhe(ix),lsng1,lsng,imo
          else
           do 4029 ieng=1,lsng
            if(k.eq.imain(ieng)) then
             mhe(imhe) = ix*imo
             imhe = imhe + 1
             write(iwr,*)ix,mhe(imhe),lsng1,lsng,imo
            end if
4029       continue
          endif
  10    continue
c       if (k.eq.lsng) mh=ix*imo
c       write(iwr,*) mhe
c end of changings
        if (ix.ne.lulu) go to 508
        go to 512
 510  continue
 508  continue
      go to 505
 512  nb = (nc-j)*nt+ig
      nb = nb/iwod
      if (nb.eq.nzero) go to 505
      do 4045 j=1,nb
        read (mtape)
4045  continue
 505  continue
      if (ix.ne.lulu) then
        write (6,*)'ix=',ix,'    lulu=',lulu
        call aborts('not all reference configurations were generated')
      endif
      write (iwr,514)
      write (iwr,515)  (mn(i),i=1,lulu)
      do 4055 i=1,3
        nconf(i)=nzero
4055  continue
      do i=1,lulu
        ix = nsc(i)
        nconf(ix) = nconf(ix) + 1
      enddo
      jswh = ix
      ix   = nzero
      do 518 i=1,3
       icon(i) = ny
       nnx = nytl(i)
       nc = nconf(i)
       if (nc.eq.0) go to 518
       do 4076 j=1,nc
        ix = ix + 1
        jx = mn(ix)
        jx = (jx-1)*nshl
        do 4075 k=1,nnx
         jx = jx + 1
         ny = ny + 1
         if(ny.le.iwod) iot(ny) = jkon(jx)
4075    continue
4076   continue
 518  continue
c
      if(ny.gt.iwod) then
         call caserr('coolez : too many mains!')
      end if
c --- assigning the iot - field to jkan (in stead of equivalence)
      do 4085 ii=1,(maxshl*maxref)
         jkan(ii) = iot(ii)
4085  continue
      nko = nko2
c --- writing out the first 1000 iot's : why ????
      write(kfile) ij,nconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     x ipt0,nplu,ndub,nod,nir,loc,ideks,jkan
c--------------------------
c --- zeroing olab
      do 4094 iii=1,ndimh
         olab(iii) = nzero
4094  continue
c
*     write(iwr,*) 'iot - field :'
*     write(iwr,*) 'in coolez'
*     iii = nzero
*     iii = iii+1
*     do iii0=1,400
*       write(iwr,*) iot(iii),iot(iii+1),iot(iii+2),iot(iii+3)
*       iii = iii + 4
*     end do
      call cputim(timeb)
      call rumplz(iot,niott,ideks,ndeks)
      call cputim(timea)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*       write(iwr,*) 'end of rumplz'
*       do iii=1,10
*          write(iwr,*) 'ideks(',iii,')=',ideks(iii)
*       end do
*       write(iwr,*) '-------'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(iwr,*)
     +  'time used by rumplz :',timea-timeb,'      <--------'
*     write(iwr,*)'**  after rumplz  **'
*     write(iwr,*)'*****    kc   *******'
*     write(iwr,20)(kc(iaus),iaus=1,120)
*     write(iwr,*)'*****    kd   *******'
*     write(iwr,20)(kd(iaus),iaus=1,120)
*     write(iwr,*)'*****   olab  *******'
*     write(iwr,21)(olab(iaus),iaus=1,150)
      call cputim(timea)
      call rumpb(iot,niott,ideks,ndeks,jcon1,nd8)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*       write(iwr,*) 'end of rumpb'
*       do iii=1,10
*          write(iwr,*) 'ideks(',iii,')=',ideks(iii)
*       end do
*       write(iwr,*) '-------'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(iwr,*)
     +   'time used by rumpb  :',timea-timeb,'      <--------'
      call cputim(timeb)
*     write(iwr,*)'**  after rumpb   **'
*     write(iwr,*)'*****    kc   *******'
*     write(iwr,20)(kc(iaus),iaus=1,120)
*     write(iwr,*)'*****    kd   *******'
*     write(iwr,20)(kd(iaus),iaus=1,120)
*     write(iwr,*)'*****   olab  *******'
*     write(iwr,21)(olab(iaus),iaus=1,150)
      kc(jto+1) = nzero
      kc(iwod)  = nzero
      write (ntype) kc,kd
      call pack(olab8,8,olab,nnid)
      write (kclab) olab8
      call rewftn(ntype)
      call rewftn(kclab)
_IF(notused)
      go to 521
 520  jto = nzero
*      do 4152 jblk=1,nzmil
*       read (ntype) kc
*       if(kc(ndimh).eq.nzero) go to 887
*4152  continue
c --- endless loop over kc: to count.
c     must be changed !
 777  continue
       read (ntype) kc
       if(kc(ndimh).eq.nzero) go to 887
      go to 777
 887  continue
      do 4157 j=1,iwod
       jto = jto + 1
       if(kc(jto).eq.0) go to 886
4157  continue
 886  continue
      jto = jto - 1
      jblk= jblk - 1
      call rewftn(ntype)
c
 521  continue
_ENDIF
c --- changed call to determine the geyser-space
      if (egey1.lt.rzero) then
c ---
       call cputim(timeb)
       call stilet
       call cputim(timea)
       write(iwr,*)
     +   'time used by stilet  :',timea-timeb,'      <--------'
*
       call cputim(timeb)
       call bambiu(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     + sx,lsx,
     + pey,acoul,aexc,ideks,ndeks,
     + peyo,acoulo,aexco,nf31o,
     + sdelta,seng,lbuffs,vect,wect,coff,nvect,
     + ot,nt444,jkan,njkan,1)
       call cputim(timea)
       write(iwr,*)
     +   'time used by bambiu :',timea-timeb,'      <--------'
       return
      else
       call aborts('program error : not mendable !')
      end if
c --- format statements
c  2  format(14i5)
c  7  format(25i3)
   9  format(1x,'automatic inclusion of single excitations with ',
     .'respect to mains '/25i3)
c 20  format(1x,20(i6))
c 21  format(1x,30(i4))
c380  format (2x,20i6)
 514  format(5x,'permuted order of reference configurations')
 515  format(/20x,10i3/)
c670  format(24i3)
      end
*************************************************************************
*###### library routins for the mrdci-program
*
* v0.0a (21.8.1992)
*
*
      subroutine aborts(messag)
      implicit none
INCLUDE(common/iofile)
      character*(*) messag
      character*80 m1,m2
c
c     error message gets printed and the program execution terminated
c
      m1='=============================================='
      m2='!!!!   program termination :              !!!!'
      write(iwr,900) m1
      write(iwr,900) m2
      write(iwr,*) messag
      write(iwr,900) m1
      call flsout
c --- format
900   format(a80)
      call caserr('program termination !')
      return
      end

c-----------------------------------------------
      subroutine ft31sp
c
c --- the program splits ft31 in a 1 and 2-electron file and
c     on ft71 (1-electron-integrals) and
c        ft72 (2-electron-integrals)
c        ft72 (record1 and rewcord2) nhead
c
c     (1992) h.u. suter
c
      implicit REAL (a-h,o-z)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
c rh
      REAL e,c,vnuc,zero
c
      integer lg,ibal,itil,mcomp
      common /junk2/ lg(8),ibal(8),itil(8),mcomp(100)
c
      integer knu,n,nod,jod,ksum,nbox,nsym
      integer iorbs,nit,nn,i
      integer nfone, nftwo
      integer nsym0
      integer ncimo
      integer ndum,nform
c --- dimensions : to be removed
c --- commons
      integer ncomp
      common /rec2/ ncomp(100)
      common /cczero/ zero
      common /ft31in/ n
      common /cffor/  nform
c
      integer nston, mtape, mdisk
      integer ideli, ltype, linf,  ntab,  kfile
      integer kclab, ntype, mstvt, nf01, nf31
      integer nhead, nf99, mtapev, nf11
      common /ftap/ nston, mtape, mdisk,
     .             ideli, ltype, linf,  ntab,  kfile,
     .             kclab, ntype, mstvt, nf01, nf62,
     +             nhead, nf99, mtapev, nf11
c
      integer ical,  nele,  nko
      integer mxex,  nmul,  ispace, nprin, ncorci
      common /tap/ ical,  nele,  nko,
     .             mxex,  nmul,  ispace, nprin, ncorci
c --- first record from ft31
      integer nconf, jab
      integer mj, kj, nj, ntil, nbal, isym
      integer iab, jcon, icuf, iduf
      integer kcon, ilee, jkon, nop
      integer isw, nyzl, ijoe, ndub
c
      parameter (nmmo=256)
      common /aeng/ nconf(ndk5),jab(n36),isym(maxsym),
     .              mj(maxsym),iab(nmmo),kj(maxsym),
     .              jcon(nmmo),icuf(n12),iduf(maxshl),
     .              kcon(ndimh),ilee(maxsym+1),
     .              ntil(maxsym),nbal(maxsym+1),
     .              jkon(maxref*maxshl),nop(maxref),
     .              isw(ndk5),nyzl(maxref),nj(maxsym),
     .              ijoe(maxref),ndub(maxref)
***   common /efeld/ e(nmmo*(nmmo+1)/2)
      common /efeld/ e(2475)
      integer lsym
      common /clsym/ lsym(800)
c
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
c
      integer lj
      common /clj2/ lj(maxsym)
c
      common /cimo/ nbox
c --- second record of ft31
      common /cnit/ nit(nitmax)
      common /b/nod,nn
      integer ij
      common /ccij/ij(maxsym)
c
      nfone=72
      nftwo=71
c --- read first record
      call rewftn(nston)
      if (nform.gt.nzero) then
          read (nston,err=61000,end=61001) 
     +    n , nod , jod , ksum , nbox , mj,kj,lj,nj,nsym0,
     .    ntil, nbal, isym, jab, iorbs, knu,
     .    lsym, ncomp,vnuc,zero
      else
         read (nston,err=61000,end=61001) 
     +    n , nod , jod , ksum , nbox , mj,kj,lj,nj,nsym0,
     .    ntil, nbal, isym, jab, iorbs, knu,
     .    lsym, ncomp,e, vnuc,zero
      end if
cdebug
*     write(iwr,*)
c     write(iwr,*)'first record of ft31'
      write(iwr,'(A,i16)')' number of irreducible representations :'
     1 ,n
c     write(iwr,*)'nod=',nod
c     write(iwr,*)'jod=',jod
c     write(iwr,*)'ksum=',ksum
      write(iwr,'(A,i16)')' number of orbitals                    :'
     1 ,nbox
c     write(iwr,*)'mj=',mj
c     write(iwr,*)'kj=',kj
c     write(iwr,*)'lj=',lj
c     write(iwr,*)'nj=',nj
      write(iwr,'(A,f16.8)')' nuclear energy                        :'
     1 ,vnuc
c     write(iwr,*)'zero=',zero
*     write(iwr,*)
c
      call rewftn(nhead)
c
      write(nhead) n , nod , jod , ksum , nbox , mj,kj,lj,nj,nsym0,
     .             ntil, nbal, isym, jab, iorbs, knu,
     .             lsym, ncomp,e,vnuc,zero
c --- read second record
cdebugstop
***   
      read (nston,err=61000,end=61001) 
     +   nit,nn,lg,ibal,itil,mcomp
      if (nform.eq.1) then
        write(nhead) nit,nn,lg,ibal,itil,mcomp
      else
        write(nhead) nit,nn,lg,ibal,itil,mcomp
      end if
      call rewftn(nston)
      read(nston,err=61000,end=61001)
      read(nston,err=61000,end=61001)nit,nn,ij
cdebugwrite(iwr,*)'ij=',ij
cdebug
cdebugread (nston) stonn
cdebugwrite(iwr,*) (stonn(i),i=1,20)
cdebug
      write(iwr,'(A,i16)')' total number of integrals             :'
     1 ,nn
c --- read 2-electronen integrals upto record end
c     write(iwr,*) ' integralfiles not yet reformatted'
      call rewftn(nston)
      return
c --- format statements
61000 write(iwr,61005)
61005 format(1x,'error in reading transformed integral interface')
      call caserr(
     + 'error in reading transformed integral interface')
61001 write(iwr,61006)
61006 format(1x,'end of transformed integral interface')
      call caserr(
     + 'unexpected end of transformed integral interface')
      end
      subroutine geysez(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +                  sx,lsx,
     +                  pey,acoul,aexc,ideks,ndeks,
     +                  peyo,acoulo,aexco,nf31o,
     +                  sdelta,seng,lbuffs,
     +                  vect,wect,coff,nvect,
     +                  iot,niott,
     +                  ot,nt444,jcon0,nd8,jkan,njkan,
     +                  jtest)
c
c  --- computes the interaction with the geyser0-space and determines
c      the geysez1--space, which get stored on ft99
c
      implicit REAL (a-h,o-z)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
      REAL hy,hp,a,w,e,f,sx
      integer nt64, mxref2,ndimf, lsx
      dimension hy(nt64), hp(nt64), a(nt64), e(ndimf)
      dimension w(mxref2)
      dimension f(ndimf)
      dimension sx(lsx)
c
      REAL pey,acoul,aexc
      integer ndeks,ideks
      dimension ideks(ndeks)
      dimension pey(ndeks),acoul(ndeks),aexc(ndeks)
c
      REAL peyo
      REAL acoulo,aexco
      integer nf31o
      dimension peyo(nf31o),acoulo(nf31o),aexco(nf31o)
c
      REAL sdelta, seng
      integer lbuffs
      dimension sdelta(5,lbuffs), seng(lbuffs)
c
      REAL vect, wect, coff
      integer nvect,jtest
      dimension vect(nvect), wect(nvect), coff(nvect)
      dimension jtest(nvect)
c
      integer iot
      dimension iot(niot)
c
      integer ot
      dimension ot(nt444)
c
      integer jcon0, nd8
      dimension jcon0(nd8)
c
      integer jkan, njkan
      dimension jkan(njkan)
c
*     REAL timeb,timea,rzero
      integer n,imo0
*     parameter (rzero = 0.0d0)
c --- commons
c
      common /cimo/   imo0
      integer icisi
      common /cisi/   icisi(maxref)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      integer icon
      common /cicon/  icon(ndk5)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
      integer nytl, ndub
      common /cny /   nytl(ndk5),ndub(ndk5)
      integer nod
      common /cnod/   nod(ndk5)
c
* jcon feld wird auf imo*(lulu+1) ausgenullt,
* daher sollte es auf nmmo*maxref dimensioniert sein
      integer lj
      common /clj2/ lj(maxsym)
c
      integer ij
      common /ccij/   ij(maxsym)
      integer kc,kd,lab
      common /ckd/    kc(ndimh),kd(ndimh),lab(3)
      integer nhead, mtapev
      common /ftap/ nston,mtape,mdisk,ideli,ltype,linf,
     + ntab,kfile,kclab,ntype, mstvt, nf01, nf62,
     + nhead, nf99, mtapev, nf11
      common /tap/ ical,ma,nko,mxex,nmul,ispace,nspipr,ncorci,
     +             msng,
     +             maxci2,ipt0in,isymm(mxroot)
c --- hier wird ein teil von tap auf b transferiert nston z.b.
      integer nsc, mn
      common /b/ iwod,lg,irsf,i1,i2,i3,i4,nc,lt,lz,nnx,nd,i5,
     +i8,i9,j1,iswh,nr,mm,jblk,jto,igmax,nr1,j2,j3,j4,nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +nn(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
c
      REAL rtol,cptol,cptolcc,cptolm,hstor
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      common /einf/   mhe(maxref),imain(maxref),iselct
      common /ft31in/ n
      logical bypass, oprin, debugd, debugtab,debugci
      common /adlin/ hstor, rtol, cptol, cptolm, cptolcc,
     +                  nroot, ntch, kprin, iselec,
     1                  ndeci, icodea, konf, keps,
     2                  ioint, norhp, izus, issk, 
     4                  ifirst, ilast, istart, ndav,
     5                  iggey, nforma,
     5                  bypass(6),oprin(6,3),
     6                  debugd, debugtab,debugci,
     7                  nteint,mdi,iotm,nedim
c
      integer labs,mkon,nplu,nyzl
      REAL xemp
      common /miscop/ xemp(500),mkon(ndimh),labs(ndimh),nplu(ndk5),
     +                nyzl(maxref)
c
      integer jkon
_IF(notused)
      REAL ptiw
      common/junk/ ptiw(5500),jkon(ndcon)
_ELSE
      common/junk/ jkon(ndcon)
_ENDIF
c
*     REAL vnuc,zero,core,corea
c
      call rewftn(linf)
      call rewftn(mtape)
      call rewftn(nhead)
      call rewftn(ntype)
      call rewftn(kclab)
      call rewftn(kfile)
      lulu=nroot
      nrootx=nroot
      lsng = msng
      ipt0 = ipt0in
      nprin = nspipr
c --- write header
      write(iwr,*) '==================================================='
      write(iwr,*) 'Compute the discarded energy contributions (geysez)'
      write(iwr,*) '==================================================='
      write(iwr,*)

c changings to take more single excitations into account
      do 386 ieng=1,maxref
        mhe(ieng)  = nzero
        imain(ieng)= nzero
 386  continue
c counter in mhe array
      imhe=1
c end of changings
c --- umkpieren von /tap/ nach /b/
      i1 = kclab
      i2 = ntype
      i4 = ntab
      i3 = ical
      i5 = nston
      i8 = mdisk
      i9 = ideli
      j1 = mtape
      j2 = ltype
      j3 = linf
      j4 = kfile
c --- number of mo's
      nirm = nmmo
      igmax= n3zt
c --- number of mo's : ndeks - 1
      im3  = ndeks - 1
      nshl = maxshl
c --- why nsec =80 ???????
      nsec = maxref
c ---- jkon read in : configurations ????
      read(mtape) nbox,m,nmul,ispace,mconf,nytl,mxex,jkon,nko,iswh,nyzl
c
      jto=nzero
c ---- for the table
      irsf = 500
      lko  = maxref
      incr = 10
      msec = nsec+1
      nawk = maxele/2+4
      nrmax= 4
c --- compute binomial n above 2 iteratively !!!
      ideks(1)=nzero
      do 428 i=1,im3
       ideks(i+1)=ideks(i)+i
 428  continue
c --- compute binomial n above 3 iteratively !!!
      jdeks(1)=nzero
      do 433 i=1,(nopmax-1)
       jdeks(i+1)=jdeks(i)+ideks(i+1)
 433  continue
c --- compute binomial n above 4 iteratively !!!
      kdeks(1)=nzero
      do 436 i=1,(nopmax-2)
       kdeks(i+1)=kdeks(i)+jdeks(i+1)
 436  continue
c --- recopy !
      imo = imo0
c     write(iwr,*) 'geysez : iwod=',iwod
c
      ix = nzero
      do 453 i=1,n
       kap=lj(i)
       if (kap.ne.nzero) then
        do 451 j=1,kap
         ix      = ix + 1
         nir(ix) = i
         loc(ix) = j
 451    continue
       end if
 453  continue
      nnx  = nmul - 2
      ndz = nytl(1) - nnx
      do 462 i=1,iswh
       nnx = nnx + 1
       ndz= ndz - 1
       nplu(i) = nnx
       nod(i)  = nplu(i) + i - 1
       ndub(i) = ndz
 462  continue
      if (ical.gt.1) go to 520
c rwah
c     read (ird,2) lulu,nrootx,lsng,nprin,ipt0,isec
      write (iwr,2) lulu,nrootx,lsng,nprin,ipt0,isec
c change to perform modified selectionscheme
      iselct = nzero
c     if(lsng.lt.0) then
c       iselct = 1
c       lsng   =-lsng
c     endif
c end of change
      if (mxex.eq.1) go to 500
      lulu = nko
      go to 501
 500  continue
c     this line to fix existing problem with exci=1 (27/7)
      lulu = nko
      if (lulu.eq.nko) go to 501
c     read (ird,2) (jtest(i),i=1,lulu)
      do i=1,lulu
          jtest(i)=ifirst+i-1
      enddo
      ix = nzero
      do 502 i=1,lulu
       jt = jtest(i)
       if (jt.eq.i) go to 503
       jx = (jt-1)*nshl
       do 488 j=1,nshl
        jx = jx + 1
        ix = ix + 1
        jkon(ix) = jkon(jx)
 488   continue
       go to 502
 503   ix = ix + nshl
 502  continue
 501  ny = nzero
ccccccccccc changing to take more single excitations into account
      lsng1 = nzero
      if(.not.(mxex.lt.2.or.lsng.eq.nzero)) then
        if(abs(lsng).ge.lulu) then
          lsng =lulu
          lsng1=1
          write(iwr,*)
          write(iwr,*) 'all single excitations are taken into account'
          do ii=1,lsng
            icisi(ii) = (ii)
          end do
        else
c         read(ird,*)(imain(ieng),ieng=1,lsng)
c --- mains projected on common /cisi/
          imain(1) = lsng 
          lsng = 1
          do 507 ii=1,lsng
            icisi(ii) = imain (ii)
 507      continue
          write(iwr,*) 'the single excitations of the following mains'
          write(iwr,*)(imain(ieng),ieng=1,lsng)
          write(iwr,*) 'are taken into account.'
        endif
c --- all singles with respect to all mains lsng < 0
        if (lsng.lt.nzero) then
          write(iwr,*) 'all single excitations are taken into account'
          lsng = lulu
          lsng1=1
        end if
      endif
c  8  continue
c     write(iwr,*) lsng,lsng1
c end of changes
      ix = nzero
      do 505 i=1,3
        nc = mconf(i)
        if (nc.eq.0) go to 505
        nt = nytl(i)
        if (ix.lt.lulu) go to 506
        nb = (nc*nt)/iwod+1
        do 530 j=1,nb
          read (mtape)
 530    continue
        go to 505
 506    read (mtape) mkon
        ig=nzero
        do 508 j=1,nc
          do 509 k=1,nt
           ig = ig + 1
           jtest(k) = mkon(ig)
           if (ig.ge.iwod) then
            ig = nzero
            read (mtape) mkon
           end if
 509      continue
        kg = -nshl
        do 510 k=1,lulu
          kg = kg + nshl
          if (nyzl(k).ne.nt) go to 510
          l7 = kg
          do 551 l=1,nt
            l7 = l7 + 1
            if (jtest(l).ne.jkon(l7)) go to 510
 551      continue
          ix = ix + 1
          mn(ix) = k
          nn(ix) = j
          nsc(ix) = i
c change to include all single excitations with respect to more than
c one reference conf.
          if(lsng.eq.nzero) go to 10
          if(lsng1.eq.1) then
           mhe(ix) = ix*imo
           if(nprin.gt.1) write(iwr,*)mhe(ix),lsng1,lsng,imo
          else
           do 569 ieng=1,lsng
            if(k.eq.imain(ieng)) then
             mhe(imhe) = ix*imo
             imhe = imhe + 1
             if(nprin.gt.1) write(iwr,*)ix,mhe(imhe),lsng1,lsng,imo
            end if
 569       continue
          endif
  10    continue
c       if (k.eq.lsng) mh=ix*imo
c       write(iwr,*) mhe
c end of changings
        if (ix.ne.lulu) go to 508
        go to 512
 510  continue
 508  continue
      go to 505
 512  nb = (nc-j)*nt + ig
      nb = nb/iwod
      if (nb.eq.nzero) go to 505
      do 585 j=1,nb
        read (mtape)
 585  continue
 505  continue
      if (ix.ne.lulu) then
        write (iwr,*)'ix=',ix,'    lulu=',lulu
        call aborts('not all reference configurations generated')
      endif
      write (iwr,514)
      write (iwr,515)  (mn(i),i=1,lulu)
      do 516 i=1,3
       nconf(i) = nzero
 516  continue
      do i=1,lulu
        ix = nsc(i)
        nconf(ix) = nconf(ix) + 1
      enddo
      jswh = ix
      ix   = nzero
      do 518 i=1,3
       icon(i) = ny
       nnx = nytl(i)
       nc = nconf(i)
       if (nc.eq.0) go to 518
       do 616 j=1,nc
        ix = ix + 1
        jx = mn(ix)
        jx = (jx-1)*nshl
        do 615 k=1,nnx
         jx = jx + 1
         ny = ny + 1
         if(ny.le.niot) iot(ny) = jkon(jx)
 615     continue
 616   continue
 518  continue
c
cbe changed:
* the iot field is niot big and not iwod big
* this needs changing to enable handling really big 
* main sets.
* changed, b.e.    28. juni 1995
      if(ny.gt.niot) then
       write(iwr,*) 'ny,niot',ny,niot
       write(iwr,*) 'geysez : too many mains! ny greater than niot'
       call caserr('geysez : too many mains!')
      end if
c --- assignment of iot - field to jkan (in stead of equivalence)
      do 627 ii=1,(maxshl*maxref)
         jkan(ii) = iot(ii)
 627  continue
c --- are the mains on the iot-field ??????
*     write(iwr,*) 'iot'
*     do i=1,300
*       write(iwr,*) 'iot(',i,')=',iot(i)
*     end do
c --- writing out the first 1000 iot's : for mains !!!!
      write(kfile) ij,nconf,mn,jswh,lulu,nsc,nrootx,isec,nprin,
     . ipt0,nplu,ndub,nod,nir,loc,ideks,jkan
c--------------------------
c
cdebug
c     write(iwr,*)'call to rumplz'
cdebug
c     call cputim(timeb)
      call rumplz(iot,niott,ideks,ndeks)
c     call cputim(timea)
c     write(iwr,*)'time used by rumplz :',timea-timeb,'      <--------'
*     write(iwr,*)'**  after rumplz  **'
*     write(iwr,*)'*****    kc   *******'
*     write(iwr,20)(kc(iaus),iaus=1,120)
*     write(iwr,*)'*****    kd   *******'
*     write(iwr,20)(kd(iaus),iaus=1,120)
*     write(iwr,*)'*****   olab  *******'
*     write(iwr,21)(olab(iaus),iaus=1,150)
c     call cputim(timea)
      call rumpb(iot,niott,ideks,ndeks,jcon0,nd8)
c     write(iwr,*)'time used by rumbp  :',timea-timeb,'      <--------'
c     call cputim(timeb)
*     write(iwr,*)'**  after rumpb   **'
*     write(iwr,*)'*****    kc   *******'
*     write(iwr,20)(kc(iaus),iaus=1,120)
*     write(iwr,*)'*****    kd   *******'
*     write(iwr,20)(kd(iaus),iaus=1,120)
*     write(iwr,*)'*****   olab  *******'
*     write(iwr,21)(olab(iaus),iaus=1,150)
      kc(jto+1) = nzero
      kc(iwod)  = nzero
      write (ntype) kc,kd
      call pack(olab8,8,olab,nnid)
      write (kclab) olab8
      call rewftn(ntype)
      call rewftn(kclab)
      go to 521
 520  jto = nzero
csut      do 884 jblk=1,nzmil
csut       read (ntype) kc
csut       if(kc(ndimh).eq.nzero) go to 887
csut  884  continue
 777  continue
       read (ntype) kc
       if(kc(ndimh).eq.nzero) go to 887
      go to 777
csut   write(iwr,*) 'increase nzmil !',nzmil,' = nzmil'
csut   call aborts('ntype-file too small !')
  887 continue
      do 885 j=1,iwod
       jto = jto + 1
       if(kc(jto).eq.0) go to 886
 885  continue
 886  jto = jto - 1
      jblk= jblk - 1
      call rewftn(ntype)
c
 521  continue
c --- changed call for determining the geyser-space
*     if (egey1.lt.rzero) then
c      call cputim(timeb)
       call stilet
c      call cputim(timea)
c      write(iwr,*)'time used by stilet  :',timea-timeb,'      <--------'
c      call cputim(timeb)
c -- printing /b/
*      write(iwr,*) ' geysez:'
*      write(iwr,*) 'line 5:'
*      write(iwr,*)
*    5 ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr
       call bambiu(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +             sx,lsx,
     +             pey,acoul,aexc,ideks,ndeks,
     +             peyo,acoulo,aexco,nf31o,
     +             sdelta,seng,lbuffs,
     +             vect,wect,coff,nvect,ot,nt444,jkan,njkan,0)
c      call cputim(timea)
c      write(iwr,*)'time used by bambiu :',timea-timeb,'      <--------'
       return
*     else
*      call aborts('program error : not mendable!')
*     end if
*****   762 write(iwr,763) ny,iwod
*****   763 format(/10x,'storage of mains takes up too much space',2i6)
*****       stop
c --- format statements
   2  format(14i5)
***7  format(25i3)
c  9  format(1x,'automatic inclusion of single excitations with ',
c    *'respect to mains '/25i3)
c 20  format(1x,20(i6))
c 21  format(1x,30(i4))
c380  format (2x,20i4)
 514  format(/1x,'changed ordering of reference configurations')
 515  format(/10x,15i4/)
      end
*###### library routins for mrdci-program
*
* v0.0a (21.8.1992)
*
*

      subroutine eigenp(a,r,n,mv)
c
c primitive jacobi-diagonalisation
c
c
      implicit REAL (a-h,o-z)
      integer ij,j,i,n
      dimension a(n*(n+1)/2),r(n*n)
      REAL a,r,anorm,anrmx,thr,x,y,sinx,sinx2,cosx,cosx2,
     1 sincs
      if(mv-1) 10,25,10
  10  iq=-n
      do j=1,n
       iq = iq + n
       do i=1,n
        ij = iq + i
        r(ij)=0.0d0
        if(i.eq.j) r(ij)=1.0d0
       end do
      end do
  25  anorm=0.0d0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
  30  ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
  35  continue
      if(anorm) 165,165,40
  40  anorm=1.414d0*dsqrt(anorm)
      anrmx=anorm*1.0d-12/dfloat(n)
      ind=0
      thr=anorm
  45  thr=thr/dfloat(n)
  50  l=1
  55  m=l+1
  60  mq=((m*m)-m)/2
      lq=((l*l)-l)/2
      lm=l+mq
      if(dabs(a(lm))-thr) 130,65,65
  65  ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
  70  y=-y
  75  if(y .gt. 1.0d0) y=1.0d0
      if(y .lt. -1.0d0) y=-1.0d0
      sinx=y/dsqrt(2.0d0*(1.0d0+(dsqrt(1.0d0-y*y))))
      sinx2=sinx*sinx
      cosx=dsqrt(1.0d0-sinx2)
      cosx2=cosx*cosx
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
  80  if(i-m) 85,115,90
  85  im=i+mq
      go to 95
  90  im=m+iq
  95  if(i-l) 100,105,105
 100  il=i+lq
      go to 110
 105  il=l+iq
 110  x=a(il)*cosx-a(im)*sinx
      a(im)=(a(il)*sinx)+(a(im)*cosx)
      a(il)=x
 115  if(mv-1) 120,125,120
 120  ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=(r(ilr)*sinx)+(r(imr)*cosx)
      r(ilr)=x
 125  continue
      x=2.0d0*a(lm)*sincs
      y= (a(ll)*cosx2)+(a(mm)*sinx2)-x
      x=(a(ll)*sinx2)+(a(mm)*cosx2)+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
 130  if(m-n) 135,140,135
 135  m=m+1
      go to 60
 140  if(l-n+1) 145,150,145
 145  l=l+1
      go to 55
 150  if(ind-1) 160,155,160
 155  ind=0
      go to 50
 160  if(thr-anrmx) 165,165,45
 165  iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
 170  x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
 175  do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
 180  r(imr)=x
 185  continue
      return
      end
c-----------------------------------------------------------------------
cjj   subroutine flsout
c
c     system independent flush routine for the output file (i/o unit 6).
c     this routine explicitly closes the output file then
c     reopens it and skips to the end of file.
c
c     character*100 filn06
c     logical file6
c
c     ------------------------
c     if the file exists close
c     ------------------------
c     inquire(file=filn06,exist=file6)
c     if (file6) then
c        close(unit=6,status='keep')
c     end if
c     ------------------------
c     now open it back up and
c     the end-of-file.
c     ------------------------
c     open(unit=6,status='unknown',form='formatted')
c
c 100   continue
c     read(6,900,end=120)
c      go to 100
c120   continue
c900   format(1x)
c      return
cc...end of flsout
c      end
c*
c*
c*
      subroutine flsout
      implicit REAL (a-h,o-z)
INCLUDE(common/iofile)
      call flushn(iwr)
      return
      end
_IF(notused)
c---------------------------------------------------------------------
c     fortran77 dreadp/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     drclos only (for hand-linked new mrd-ci package)  (m.c. and g.j.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine drclos
      implicit REAL (a-h,o-z)
      common /directiop/ lu,nbytes,nwords,nrec,irc,iopt
      if (iopt.eq.1) then
      write (6,*) 'drclos'
      endif
      close (lu)
      return
      end
c---------------------------------------------------------------------
c     fortran77 dreadp/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c   dreadp only  (for hand-linked new mrd-ci package)  (g.j. and m.c.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dreadp(a,nr)
      implicit REAL (a-h,o-z)
      integer a(nwords)
      common /directiop/ lu,nbytes,nwords,nrec,irc,iopt
      if (iopt.eq.1) then
      write (6,*) 'dreadp nrec=',nr,' nbytes=',nbytes,
     *' nwords=',nwords
      endif
      read (lu,rec=nr,iostat=irc) a
      if (irc.ne.0) then
       write (6,100) lu,nr,nrec,nbytes,irc
 100   format(' ######  dreadp error on logical unit ',i10/,
     *        '        record number is             ',i10/,
     *        '        written by previous dwrite   ',i10/,
     *        '        record length (bytes)        ',i10/,
     *        '        returncode is                ',z8/)
        call aborts('error reading the da-file')
      endif
      return
      end
c---------------------------------------------------------------------
c     fortran77 dreadp/dwrite routine - bernd hess - 22.5.87
c      2.6.89 - ibm version: status='delete'
c     25.5.88 - drdel
c---------------------------------------------------------------------
c  delcx        convex version
c  delc@        debug
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c  dropenp only (for hand-linked new mrd-ci package)  (g.j. and m.c.)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dropenp(lun,nbyt)
      implicit REAL (a-h,o-z)
      common /directiop/ lu,nbytes,nwords,nrec,irc,iopt
c
c     machine dependent parameter: word length
c
      logical open
      data nwlen /4/
      if (iopt.eq.1) then
      write (6,*) 'dropenp lu=',lun,' nbyt=',nbyt
      endif
      nbytes=(nbyt/nwlen)*nwlen
      if (nbytes.ne.nbyt) then
       write (6,101) lun,nbyt,nwlen
 101   format(' ###### dropenp error on logical unit ',i10/,
     *        ' record length (in bytes)           : ',i10/,
     *        ' must be divisible by word length   : ',i10)
       call caserr('error 1 in opening tab file in dropenp')
      endif
      lu=lun
      nwords=nbyt/nwlen
      inquire (lun,opened=open)
cbe   if (.not.open) then
      open (lu,iostat=irc,access='direct',form='unformatted',
     +      status='unknown',file='table-ci',recl=nbyt)
      if (irc.ne.0) then
       write (6,100) lu,nbytes,irc
 100   format(' ###### dropenp error on logical unit ',i10/,
     *        '        record length (bytes)         ',i10/,
     *        '        returncode is                 ',z8/)
       call caserr('error 2 in opening tab file in dropenp')
      endif
c     endif
c@    write (6,*) ' dropenp lu',lu,' lrecl=',nbytes
      return
      end
_ENDIF
c******************************
      subroutine cputim(cpusec)
c******************************
c
c     obtain the elapsed cpu time (in seconds)
c
c                                                       -------
c                                                       vax/vms
c                                                       -------
cvax  integer*2 cptmcd,iffour,len2a(8)
cvax  integer cputm,cptmad,izero,izero1,sys$getjpi
cvax  double precision cpusec
cvax  equivalence (len2a(1),ifour ),(len2a(2),cptmcd),(len2a(3),cptmad),
cvax *            (len2a(5),izero ),(len2a(7),izero1)
cvax  data izero/0/, izero1/0/, ifour/4/, cptmcd/1031/
cvax  cptmad=%loc(cputm)
cvax  if(.not. sys$getjpi(,,,len2a,,,)) write(6,1)
cvax  cpusec=cputm/100.0d0
cvax1 format(1x,'error from sys$getjpi system service routine.')
c
c                                                      ---
c                                                      sun
c                                                      ---
csun  implicit double precision (a-h,o-z)
csun  real tarray(2),etime
csun  cpusec = dble(etime(tarray))
c                                                      ------
c                                                      convex
c                                                      ------
ccvx  implicit double precision (a-h,o-z)
ccvx  real tarray(2),etime
ccvx  cpusec = dble(etime(tarray))
c                                                      --------
c                                                      dec risc
c                                                      --------
cdrs  implicit double precision (a-h,o-z)
cdrs  real tarray(2),etime
cdrs  cpusec = dble(etime(tarray))
c                                                      -----------------
c                                                      ibm vs fortran
c                                                      call an assembler
c                                                      routine to return
c                                                      cpu time in micro
c                                                      the fortran manua
c                                                      suggests that clo
c                                                      returns cpu time
c                                                      has proven to be
c                                                      -----------------
cibm  double precision cpusec,xmicro
cibm  integer irtcod
cibm  call cputime(xmicro,irtcod)
cibm  cpusec = xmicro/1.0d+06
c
c                                                      -----------------
c                                                      microsoft fortran
c                                                      -----------------
cmsf  integer*2 ihr,imin,isec,i100th
cmsf  double precision cpusec
cmsf  call gettim(ihr,imin,isec,i100th)
cmsf  cpusec = 3600*ihr + 60*imin + isec + 0.01*i100th
c                                                       ---------------
c                                                       ibm aix fortran
c                                                       ---------------
      implicit REAL (a-h,o-z)

      cpusec=cpulft(1)

c     integer ibuf(4)
c     istat  = times(ibuf)
c     user   = ibuf(1)/100.0
c     system = ibuf(2)/100.0
c     cpusec = user+system
c                                                      -----------------
c                                                      hp-ux: note this
c                                                             a c routin
c                                                      -----------------
c
chpx  call second(cpusec)
c
c205  cpusec = second()
c
cray  call qtime(icpu,io,mem,fcpu,fio,fmem)
cray  cpusec = icpu*1.0e-6
c
      return
c...end of cputim
      end
*************************************************************************
_IF(notused)
      subroutine pack1(ii,n,pi)
      logical*1 ii(2,n)
      logical*1 pi(n)
      do 1 i=1,n
_IF(littleendian)
      pi(i) = ii(1,i)
_ELSE
      pi(i) = ii(2,i)
_ENDIF
    1 continue
*     write(6,*) 'pack1:n = ', n
*     write(6,*) (ii(loop),loop=1,n)
      return
      end
_ENDIF
      subroutine parke4(ideks,ndeks,lkon,nmax,mxref,oprinsym)
c
c   generating the csf's
c
cvp mains with up to 9 open shells are accepted.
cvp  configurations up to seventh sk are generated.
cvp  resulting configurations having more then 9 open shells will not
cvp  be taken into account any further.
      implicit REAL (a-h,o-z)
c
      integer ntotcf
      integer n
      integer ideks,ndeks,lkon,nmax,mxref
      logical oprinsym
      dimension ideks(ndeks),lkon(nmax)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
c NB 
INCLUDE(common/cepa_sn2)
INCLUDE(common/cepa_mrd1)
c NB
c
      parameter(mxset=10,mxdeg=100)
      character*8 char8
      logical coreci

      integer iret,nset
      common /inkurt/ coreci

      dimension nset(3)
      dimension ndeg(3,mxset),inocc(3,mxset)
      dimension nseto(3,mxdeg,mxset)
c
       common /indat/ ndeg,inocc,nseto
c
       common /ft31in/ n
       integer nhead
       common /ftap/ nston, mtape, mdisk,
     .               ideli, ltype, linf,  ntab,  kfile,
     .               kclab, ntype, mstvt, nf01, nf62,
     +               nhead, nf99, mtapev, nf11
       common /tap/  ical,  m,     nko,
     .               mxex,  nmul,  ispace, nprin, 
     .               ncorci
c
      common /a/ b(nlca )
c
      common /cnbox/ nbox
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer lj
      common /clj/  lj(maxsym)
      integer nytl,ndub0
      common /cny/  nytl(ndk5),ndub0(ndk5)
c
      integer nconf
      integer jab, isym, mj, iab, kj, jcon, icuf, iduf
      integer kcon, ilee, ntil, nbal, jkon, nop
      integer isw, nyzl, nj, ijoe, ndub
      parameter (nmmo=256)
      common /aeng/ nconf(ndk5),jab(n36),isym(maxsym),
     .              mj(maxsym),iab(nmmo),kj(maxsym),
     .              jcon(nmmo),icuf(n12),iduf(maxshl),
     .              kcon(ndimh),ilee(maxsym+1),
     .              ntil(maxsym),nbal(maxsym+1),
     .              jkon(ndcon),nop(maxref),
     .              isw(ndk5),nyzl(maxref),nj(maxsym),
     .              ijoe(maxref),ndub(maxref)
c
      common /b1/ ig,jpaw,klx,nshl,kly,ispin
c
      logical oconf,debugs
      common/parkin/egey1,trash,tdel,
     +       ical0,nele,nkodump(9),
     +       oconf,debugs,
     +       lsng,maxciin,ipt0,isymin(mxroot)
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
c     dimension nj(8)
c
      integer imj(8),ikj(8)
      integer jkonr(maxshl)
      integer kchelp
c
      integer mkon
      dimension mkon(ndimh)
      dimension kchelp(100),inver(8)
c
      integer isaf
      dimension isaf(0:9)
cvp
cvp
c new parameter
      parameter (ndkneu=ndk5+2,nopneu=nopmax+4)
c
c fields, that dimensioned anew
      integer nytln,nconfn,iswn,jdeksn,kdeksn
      common /momain/ nytln(ndkneu),nconfn(ndkneu),
     +      iswn(ndkneu),jdeksn(nopneu),kdeksn(nopneu)
      integer nirrep(8)
      character*100 remstring
      character*4 zsymx(8)
      logical oi3
c
      data zsymx/'s1','s2','s3','s4','s5','s6','s7',
     1           's8'/
c
c   iswh remains set to 5!
cvp
cvp ???
c --- initialisations
c --- open ft99 for geyser1 and geyser0-spaces
cmvs  open (unit=99,file='geyser.raum',status='unknown')
      write(nf99,*) 'geyser-raum'
      write(nf99,*)
     .   '==========================================================='
c
      write(nf99,*) nko,m,'    ,'
      ntotcf = nzero
      nshl   = maxshl
      nops   = 5
      ibl = ndimh
_IFN(cray)
      call setsto(nnid,0,olab)
_ENDIF
c
      iswh =  5
      lko  = maxref
      mmax = maxele
c
**** rewind nston
      call rewftn(nhead)
      call rewftn(mdisk)
c ---- read typical information from ft31
c     read (nston) n,nod,jod,ksum,nbox,mj,kj,lj,nj,nsym,ntil,nbal,isym,
c    1 jab,iorb
      read (nhead) innnn,nnnod,njod,
     +     iksum,nbox,
     +    ikj,imj,lj,nj,nsym,ntil,nbal,isym,jab
c
c-- changes by hutter
c   retain only
c   single excitations out of a specified set of "core"-orbitals
c   ("core/valence-correlations only")
c-- changes by hutter
      if (ncorci .eq. 0) then
       coreci = .false.
cc change: k.pfingst 8/10/92
       else if (((ncorci .eq. 1).or.(ncorci.eq.2.)).or.(ncorci.eq.3).or.
     *(ncorci.eq.4))
     *then
        coreci = .true.
        write (iwr,867)
 867    format
     * (/1x,'only core-valence excitations are retained for following'/)
c-- read input
        read (ird,*) (nset(i),i=1,3)
        do 820 j=1,3
        if(j.eq.1) then
         char8='at least'
        elseif(j.eq.2) then
         char8='at most '
        elseif(j.eq.3) then
         char8='exact   '
        endif
        do 810 i=1,nset(j)
        read (ird,*) ndeg(j,i),inocc(j,i),
     *  (nseto(j,mdeg,i),mdeg=1,ndeg(j,i))
        write (iwr,359)  ndeg(j,i),char8,inocc(j,i)
        write (iwr,360) (nseto(j,mdeg,i),mdeg=1,ndeg(j,i))
359   format(1x,i4,' orbital(s) with',1x,a8,1x,i3,' electron(s):')
360   format(5x,20i3)
 810    continue
 820    continue
        write (iwr,812) (nset(i),i=1,3)
 812    format (/,1x,3i3,' set(s) of orbitals are taken into account')
       else
       write (iwr,868)
 868   format (/1x,'error: ncorci must be 0,1,2,3 or 4'/)

       call caserr('invalid setting for ncorci')
      endif
      nbxa=nbox-1
      nbxb=nbox-2
      nyx=nsym*(nsym+1)/2
c     write (iwr,8600)
c-- changes by hutter
c8600 format(5x/,'nbox,no.electrons,no.mains,excitations,no.irs,multipli
c    +city,ispace,nprint,ncorci')
      nnnx=(m+nmul-3)/2
cvp   do 87 i=1,iswh
      do 87 i=1,ndkneu
      nnnx=nnnx+1
 87   nytln(i)=nnnx
      do 2011 kdoo=1,ndk5
        nytl(kdoo)=nytln(kdoo)
 2011 continue
      write(iwr,8497) nbox,m,nko,mxex,nsym,nmul,ispace,nprin,
     +                ncorci
8497  format(/,' number of orbitals ',i6/,
     &         ' no. of electrons   ',i6/,
     &         ' no. of mains       ',i6/,
     &         ' excitations        ',i6/,
     &         ' no. of irreps.     ',i6/,
     &         ' multiplicity       ',i6/,
     &         ' ispace             ',i6/,
     &         ' nprint             ',i6/,
     &         ' ncorci             ',i6)
      if(mxex.gt.4.or.mxex.lt.0) go to 767
      if(nmul.lt.1) go to 783
cvp   do 88 i=1,nops
      do 88 i=1,ndk5
88     isw(i)=0
      do 2088 i=1,ndkneu
 2088  iswn(i)=0
csut here something is wrong !!! nmul : multiplicity
cvp   if(nmul.gt.1) iswn(nmul-1)=1
cvp    iswn(nmul+1)=2
cvp   if(nmul.lt.3) then
cvp    iswn(nmul+3)=3
cvp   endif
      if (nmul.eq.1) then
        icount=1
        do 2051 kdoo=2,ndkneu,2
          icount=icount+1
          iswn(kdoo)=icount
 2051   continue
       else
        icount=0
        do 2052 kdoo=nmul-1,ndkneu,2
          icount=icount+1
          iswn(kdoo)=icount
 2052   continue
      endif
      if(m.lt.2.or.m.gt.mmax) go to 779
      if(ispace.lt.1.or.ispace.gt.8) go to 785
      if(nko.gt.lko) go to 766
c8498 format(1x,12i10/)
      do 95 i=1,nsym
 95   inver(i)=0
      do 96 i=1,n
      ix=isym(i)
      if (lj(i).eq.0) go to 96
      inver(ix)=i
 96   continue
      ix=0
      do 89 i=1,n
      ilee(i)=ix
      l=lj(i)
      if (l.eq.0) go to 89
      lx=isym(i)
      do 90 j=1,l
      ix=ix+1
   90 iab(ix)=lx
   89 continue
      ilee(n+1)=nbox
      im=0
c670  format(24i3)
csut      do 2 i=1,nko
csut      nt=im+1
csut      im=im+nshl
csut  read(ird,670)np,(jkon(j),j=nt,im)
c --- read the configurations (mains) from inputfile
      restape=ird
      if (oconf) call fileconfpa
      ird=nf11
      call rewftn(ird)
      do 2 i=1,nko
       nt = im +    1
       im = im + nshl
c ---- reading the configurations
cformated read(ird,670)np,(jkon(j),j=nt,im)
*      read(ird,*)np,(jkon(j),j=nt,im)
       do 9198 j=nt,im
          jkon(j) = 0
9198   continue
       read (ird,*)np
       mm = nint((m+np)/2.0)
       read (ird,*)(jkonr(j),j=1,mm)
CMR    do 9204 j=nt,im
       do 9204 j=1,mm
          jkon(j+nt-1) = jkonr(j)
9204   continue
cdebug
cdebug   write(iwr,*) 'csf',np,jkonr
cdebug
cvp   if(np.gt.nops) go to 775
cvp reference configurations with up to nopmax open shells are allowed
      if(np.gt.nopmax) go to 775
      nyt=1
      if(np.eq.0) go to  94
      nyt = iswn(np)
 94   nyt=nytln(nyt)
      nyzl(i)=nyt
      if(i.eq.1) go to 503
      i1=i-1
      mx=-nshl
      do 504 j=1,i1
      mx=mx+nshl
      if(nop(j).ne.np) go to 504
      nv=mx
      la=nt-1
      do 505 k=1,nyt
      la=la+1
      nv=nv+1
      if(jkon(la).ne.jkon(nv)) go to 504
 505  continue
      go to 777
 504  continue
 503  na=nmul+np
      if(na-2*(na/2).eq.0) go to 771
      if(np.gt.m) go to 771
c
c     I dont understand the "np.gt.nmul+3" condition below given the
c     current setting for nopmax e.g. 6 open shells when nmul=1
c     will trigger it.
c     I've tried editing it out, but causes problems with
c     "not all mains generated"
c
      if(np.gt.nmul+3) go to 771
      nop(i)=np
      ndub(i)=(m-np)/2
      if(np.lt.2) go to 500
      nz=jkon(nt)
      if(nz.lt.1.or.nz.gt.nbox) go to 773
      do 501 j=2,np
      nt=nt+1
      nv=jkon(nt)
      if(nv.le.nz) go to 769
      if(nv.gt.nbox) go to 773
501   nz=nv
      if(np.eq.m) go to 2
518   mx=np+im-nshl+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      lb=im-nshl
      do 506 j=1,np
      lb=lb+1
      if(jkon(lb)-nv) 506,769,507
506   continue
507   jm=np+2
      if(jm.gt.nyt) go to 2
      kp=mx
      do 508 j=jm,nyt
      kg=mx
      mx=mx+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      do 509 k=kp,kg
      nz=jkon(k)
      if(nv.le.nz) go to 769
 509   continue
      kg=im-nshl
      do 510 k=1,np
      kg=kg+1
      if(jkon(kg)-nv) 510,769,508
 510  continue
 508  continue
      go to 2
 500  if(np.eq.0) go to 511
      nv=jkon(nt)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      go to 518
 511  nz=jkon(nt)
      if(nz.lt.1.or.nz.gt.nbox) go to 773
      if(nyt.eq.1) go to 2
      kp=nt
      do 512 j=2,nyt
      kg=nt
      nt=nt+1
      nv=jkon(nt)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      do 512 k=kp,kg
      nz=jkon(k)
      if(nv.le.nz) go to 769
 512  continue
   2  continue
c
c     check requested size of zero order space
c
      do  j=(nmul-1),nopmax,2
       isaf(j)=numsaf(nmul,j)
      enddo
      k=0
      do j=1,nko
        k= k + isaf(nop(j))
      enddo
      write(iwr,553) nko,k
553   format(' The set of',i5,' main configurations leads to',
     + i5,' spin adapted functions.')
      if (k.gt.mxref)
     + call caserr('too many SAFs from specified mains')
c
cvp ? end of reading and testing the reference configurations
      kg=nbox*nko
      do 513 i=1,kg
 513  lkon(i)=0
      kg=-nshl
      kb=-nbox
      do 514 i=1,nko
         kg=kg+nshl
         kb=kb+nbox
         np=nop(i)
         nyt=nyzl(i)
         if(np.eq.0) go to 516
         do 515 j=1,np
           nv=kb+jkon(j+kg)
 515     lkon(nv)=1
         if(np.eq.m) go to 514
 516     np=np+1
         do j=np,nyt
           nv=kb+jkon(j+kg)
         lkon(nv)=2
         enddo
 514  continue
      write(iwr,520)
 520  format(/10x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,
     +'numbers of open shells and corresponding main configurations'/
     + 10x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + / 5x, 'input  number of     main'
     + / 5x, 'index  open shells   configuration'/)
      kg = -nshl
      do i=1,nko
         kg=kg+nshl
         np=nop(i)
         kb=kg+1
         nv=kg+nyzl(i)
         write(iwr,522)i, np,(jkon(j),j=kb,nv)
 522  format(5x,i3,3x,i3,9x,24i4)
      enddo
c
c NB select the vacuum space MOs
c
      if (cepai) then
c
      jjopen = nmul - 1
      nvacuum = 0
      do i=1,nyzl(1)
         movac(i)=jkon(i)
         nvacuum = nvacuum + 1
      end do
      kg=0
      do i=2,nko
         kg=kg+nshl
         kb=kg+1
         nv=kg+nyzl(i)
         do 542 j=kb,nv
            do k=1,nvacuum
               if(jkon(j).eq.movac(k)) go to 542
            end do
            nvacuum = nvacuum + 1
            movac(nvacuum) = jkon(j)
 542     continue 
      enddo
c
      if (ipcepa.ge.1) then
c
      write(iwr,'(//" MOs in vacuum space")')
      write(iwr,'(24(1x,i4))') (movac(i),i=1,nvacuum)
c532  format(2x,i4,1x)
c
      end if
c
c NB doc. MOs in ref. conf. selection
c
c i runs over doubly occupied MOs in the first ref.
c j runs over all the rest of ref. conf.
c k runs over doc. MOs within a ref. conf.
c more info in common/cepa_sn2
c
      np=nop(1)
      docnum = 0
      do i=np+1,nyzl(1)
         doctest = jkon(i)
         doccnt = 1
         kg = 0
         do 432 j=2,nko
            kg=kg + nshl
            np=nop(j)
            kb=kg+1
            kb=kb + nop(j)
            nv=kg+nyzl(j)
            do k=kb,nv
               if(doctest.eq.jkon(k)) then
                 doccnt = doccnt + 1
                 go to 432
               end if
            enddo
 432     continue
         if(doccnt.eq.nko) then
            docnum = docnum + 1
            modoc(docnum)=doctest
         end if
      enddo
c
      if(ipcepa.ge.1) then
c
      write(iwr,435)
 435  format(//'Doubly occupied MOs within main configurations')
c    +       //5x,15('='),/8x,'N',4x,'doc.MOs',/5x,15('='))
      write(iwr,'(24(1x,i4))') (modoc(i),i=1,docnum)
c
      end if
c
      end if
c
c NB end of vacuum and doc. MOs selection for cepa
c        
      if (oprinsym) then
      write(iwr,'(/a/)')'  print by symmetry'
      nirrep(1)=0
      oi3=.true.
      istep=3
      do iiz=2,8
         nirrep(iiz)=nirrep(iiz-1)+lj(iiz-1)
      enddo
      do iiz=1,8
         if (lj(iiz).gt.100) then
            oi3=.false.
            istep=4
         endif
      enddo
      kg=-nshl
      do i=1,nko
         kg=kg+nshl
         np=nop(i)
         kb=kg+1
         nv=kg+nyzl(i)
         krem=1
         write(remstring(krem:krem+4),'(i4)')np
         krem=krem+4
         nirr=1
         nsymm=0
         ist=4
         iopen=0
         do iiz=kb,nv
            iorb0=jkon(iiz)
            do iiy=1,8
               if (iorb0.lt.nirrep(iiy)) then
                  isymm=iiy-1
                  goto 6100
               endif
            enddo
           isymm=8
 6100      continue
           if (iopen.eq.np) then
              write(remstring(krem:krem+3),'(a3)')'   '
              krem=krem+3
              ist=krem
           endif
           if (isymm.ne.nsymm) then
              if (krem+4.ge.100) then
                  write(iwr,'(a)')remstring(1:krem)
                  do krem=1,100
                      remstring(krem:krem)=' '
                  enddo
                  krem=ist
              endif
              if (oi3) then
                 write(remstring(krem:krem+3),'(a1,a2)')' ',
     1           zsymx(isymm)
                 krem=krem+3
              else
                 write(remstring(krem:krem+4),'(a2,a2)')'  ',
     1           zsymx(isymm)
                 krem=krem+4
              endif
              nsymm=isymm
           endif
           iorb0=iorb0-nirrep(isymm)
           if (krem+istep.ge.100) then
              write(iwr,'(a)')remstring(1:krem)
              do krem=1,100
                 remstring(krem:krem)=' '
              enddo
              krem=ist
           endif
           if (oi3) then
              write(remstring(krem:krem+istep),'(i3)')iorb0
           else
              write(remstring(krem:krem+istep),'(i4)')iorb0
           endif
           iopen=iopen+1
           krem=krem+istep
         enddo
         write(iwr,'(a)')remstring(1:krem)
      enddo
c
      endif
c
 112  format(5x,'orbital   irrep')
 113  format(5x,i3,' - ',i3,i6)
      write(iwr,8499)
 8499 format(/2x,'symmetry data'/)
      write(iwr,111)
 111  format(5x,15('='))
      write(iwr,112)
      write(iwr,111)
      ist=1
      isymm=1
      do irem=1,nbox
         if (iab(irem).ne.isymm.or.irem.eq.nbox) then
            if (irem.eq.nbox) then
               iend=irem
            else 
               iend=irem-1
            endif
            write(iwr,113)ist,iend,iab(iend)
            ist=irem
            isymm=iab(irem)
         endif
      enddo
c     write(iwr,352)(iab(i),i=1,nbox)
c     write(iwr,352)(jab(i),i=1,nyx)
c352  format(5x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,
c    110i1,1x,10i1,1x,10i1)
      ispin=(nmul-1)/2
      fms=ispin
      if(m.ne.2*(m/2)) fms=fms+0.5d0
cvp ideks bereits in geysez berechnet und zwar von 1 bis ndeks-1
cvp ideks was calculated already in geysez from 1 upto ndeks-1
c --- compute binomial n above 2 iteratively !!!
      ideks(1)=0
cvp   do 7 i=1,nsym
      do 7 i=1,nopneu-1
   7  ideks(i+1)=ideks(i)+i
c --- compute binomial n above 3 iteratively !!!
c     jdeks(1)=0
c     do i=1,(nopneu-1)
c      jdeks(i+1)=jdeks(i)+ideks(i+1)
c     enddo
c --- compute binomial n above 4 iteratively !!!
c     kdeks(1)=0
c     do i=1,(nopneu-2)
c      kdeks(i+1)=kdeks(i)+jdeks(i+1)
c     enddo
cvp   write(6,*) ' nsym=',nsym
cvp   write(6,*) ' ideks=',(ideks(iii),iii=1,20)
cvp   write(6,*) ' jdeks=',(jdeks(iii),iii=1,nopmax)
cvp   write(6,*) ' kdeks=',(kdeks(iii),iii=1,nopmax)
      kg=-nshl
      do 524 i=1,nko
      np=nop(i)
      kg=kg+nshl
      if(np.gt.0) go to 525
      jg=1
      go to 530
 525  jg=jkon(kg+1)
      jg=iab(jg)
      if(np.eq.1) go to 530
      do 527 j=2,np
      nv=jkon(kg+j)
      nv=iab(nv)
      if(nv.lt.jg) go to 528
      jg=ideks(nv)+jg
      go to 527
 528   jg=ideks(jg)+nv
 527  jg=jab(jg)
 530  continue
      nmns=np/2+1-ispin
      if(nmns.le.0) go to 771
 524  ijoe(i)=jg
      do 8 i=1,ndk5
   8  nconf(i)=0
      do 2008 i=1,ndkneu
 2008 nconfn(i)=0
      ig=0
      klx=0
      llx=0
cvp
cvp ? pseudo-loop over all reference configurations
cvp
22    kly=klx
      klx=klx+1
      if (klx.gt.nko) go to 23
      do 24 i=1,nbox
      llx=llx+1
24    jcon(i)=lkon(llx)
      kar=0
      imel=ijoe(klx)
      nbd=1
      nod=nop(klx)
      ncl=ndub(klx)
      nmns=nod/2+1-ispin
      if (klx.eq.1) go to 26
      kg=-nshl
      do 162 i=1,kly
      mx=0
      kg=kg+nshl
      np=nop(i)
      ja=ndub(i)
      kz=kg
      if (np.eq.0) go to 163
      do 164 j=1,np
      kz=kz+1
      jd=jkon(kz)
      if (jcon(jd).gt.0) go to 164
      mx=mx+1
      if (mx.gt.mxex) go to 162
  164 continue
      if (ja.eq.0) go to 28
  163 do 165 j=1,ja
      kz=kz+1
      jd=jkon(kz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 165
      mx=mx+jd
      if (mx.gt.mxex) go to 162
  165 continue
      go to  28
 162  continue
  26  if (imel.ne.ispace) go to 28
      if (nmns.le.0) go to 28
      ig=ig+1
      kz=kly*nshl
      kcon(ig)=nmns
cvp   nconf(nmns)=nconf(nmns)+1
      nconfn(nmns)=nconfn(nmns)+1
      if (ig.lt.ibl) go to 166
      ig=0
      write (mdisk) kcon
 166  if (nod.eq.0) go to 167
      do 168 j=1,nod
      kz=kz+1
      ig=ig+1
      kcon(ig)=jkon(kz)
      if (ig.lt.ibl) go to 168
      ig=0
      write (mdisk) kcon
 168  continue
      if (ncl.eq.0) go to 28
 167  do 169 j=1,ncl
      kz=kz+1
      ig=ig+1
      kcon(ig)=jkon(kz)
      if (ig.lt.ibl) go to 169
      ig=0
      write (mdisk) kcon
 169  continue
  28  if (mxex.eq.0) go to 22
      if (imel.gt.ispace) go to 170
      jpaw=ideks(ispace)+imel
      go to 171
 170  jpaw=ideks(imel)+ispace
 171  jpaw=jab(jpaw)
38    kar=kar+1
      if (kar.gt.nbox) go to 39
      if (jcon(kar).eq.0) go to 38
      iq=iab(kar)
      if (iq.gt.jpaw) go to 175
      kpaw=ideks(jpaw)+iq
      go to 176
 175  kpaw=ideks(iq)+jpaw
 176  kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 38
      kas=ilee(kpaw)
      npox=ilee(kpaw+1)
      jcon(kar)=jcon(kar)-1
30    kas=kas+1
      if (kas.eq.kar) go to 30
      if (kas.gt.npox) go to 31
      if (jcon(kas).eq.2) go to 30
      jcon(kas)=jcon(kas)+1
33    if (klx.eq.1) go to 25
      lly=-nshl
      do 34 i=1,kly
      mx=0
      lly=lly+nshl
      np=nop(i)
      ja=ndub(i)
      llz=lly
      if (np.eq.0) go to 173
      do 172 j=1,np
      llz=llz+1
      jd=jkon(llz)
      if (jcon(jd).gt.0) go to 172
      mx=mx+1
      if (mx.gt.mxex) go to 34
 172  continue
      if (ja.eq.0) go to 36
 173  do 174 j=1,ja
      llz=llz+1
      jd=jkon(llz)
      jd=2-jcon(jd)
      if (jd.eq.0) go to 174
      mx=mx+jd
      if (mx.gt.mxex) go to 34
 174  continue
      go to 36
34    continue
      go to 25
31    jcon(kar)=jcon(kar)+1
      go to 38
 3    jcon(kas)=jcon(kas)-1
      go to 30
39    if(mxex.eq.1) go to 22
      kar=0
      nbd=2
      if (jpaw.ne.1) go to 47
46    kar=kar+1
      if (kar.gt.nbox) go to 47
      if (jcon(kar).lt.2) go to 46
      jcon(kar)=0
      kas=0
41    kas=kas+1
      if(kas.eq.kar) go to 41
      if (kas.gt.nbox) go to 44
      if (jcon(kas).gt. 0) go to 41
      jcon(kas) =2
      go to 33
44    jcon(kar)=2
      go to 46
 42   jcon(kas)=0
      go to 41
 47   nbd=3
      kar=0
53    kar=kar+1
      if (kar.gt.nbox) go to 54
      if (jcon(kar).lt.2) go to 53
      jcon(kar)=0
      las=0
55    las=las+1
      if (las.eq.nbox) go to 56
      if (las.eq.kar) go to 55
      if (jcon(las).eq.2)  go to 55
      iq=iab(las)
      if (iq.gt.jpaw) go to 180
      kpaw=ideks(jpaw)+iq
      go to 181
 180  kpaw=ideks(iq)+jpaw
 181  kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 55
      kas=ilee(kpaw)
      npox=ilee(kpaw+1)
      if (kas.ge.las) go to 184
      if (las.ge.npox) go to 55
      kas=las
 184  jcon(las)=jcon(las)+1
50    kas=kas+1
      if (kas.eq.kar) go to 50
      if (kas.gt.npox) go to 51
      if (jcon(kas).eq.2) go to 50
      jcon(kas)=jcon(kas)+1
      go to 33
 48   jcon(kas)=jcon(kas)-1
      go to 50
51    jcon(las)=jcon(las)-1
      go to 55
56    jcon(kar)=2
      go to 53
54    nbd=4
      kar=0
62    kar=kar+1
      if (kar.eq.nbox) go to 63
      if (jcon(kar).eq.0) go to 62
      iq=iab(kar)
      if (iq.gt.jpaw) go to 182
      kpaw=ideks(jpaw)+iq
      go to 183
 182  kpaw=ideks(iq)+jpaw
 183  kpaw=jab(kpaw)
      kpaw=inver(kpaw)
      if (kpaw.eq.0) go to 62
      lar=ilee(kpaw)
      npox=ilee(kpaw+1)
      if (lar.ge.kar) go to 185
      if (kar.ge.npox) go to 62
      lar=kar
 185  jcon(kar)=jcon(kar)-1
64    lar=lar+1
      if (lar.gt.npox) go to 65
      if (jcon(lar).eq.0) go to 64
      jcon(lar)=jcon(lar)-1
      kas=0
59    kas=kas+1
      if (kas.eq.lar) go to 59
      if (kas.eq.kar) go to 59
      if (kas.gt.nbox) go to 60
      if (jcon(kas).gt.0)  go to 59
      jcon(kas)=2
      go to 33
60    jcon(lar)=jcon(lar)+1
      go to 64
65    jcon(kar)=jcon(kar)+1
      go to 62
 57   jcon(kas)=0
      go to 59
63    nbd=5
      kar=0
71    kar=kar+1
      if (kar.eq.nbox) go to 581
      if (jcon(kar).eq.0) go to 71
      jcon(kar)=jcon(kar)-1
      iq=iab(kar)
      if (iq.gt.jpaw) go to 186
      kpaw=ideks(jpaw)+iq
      go to 187
186   kpaw=ideks(iq)+jpaw
187   kpaw=jab(kpaw)
      lar=kar
72    lar=lar+1
      if (lar.gt.nbox) go to 73
      if (jcon(lar).eq.0) go to 72
      jcon(lar)=jcon(lar)-1
      iq=iab(lar)
      if (iq.gt.kpaw) go to 188
      lpaw=ideks(kpaw)+iq
      go to 189
188   lpaw=ideks(iq)+kpaw
189   lpaw=jab(lpaw)
      las=0
74    las=las+1
      if (las.eq.nbox) go to 75
      if (las.eq.kar) go to 74
      if (las.eq.lar) go to 74
      if (jcon(las).eq.2) go to 74
      iq=iab(las)
      if (iq.gt.lpaw) go to 190
      mpaw=ideks(lpaw)+iq
      go to 191
 190  mpaw=ideks(iq)+lpaw
 191  mpaw=jab(mpaw)
      mpaw=inver(mpaw)
      if (mpaw.eq.0) go to 74
      kas=ilee(mpaw)
      npox=ilee(mpaw+1)
      if (kas.ge.las) go to 192
      if (las.ge.npox) go to 74
      kas=las
 192  jcon(las)=jcon(las)+1
68    kas=kas+1
      if (kas.eq.kar) go to 68
      if (kas.eq.lar) go to 68
      if (kas.gt.npox) go to 69
      if (jcon(kas).eq.2) go to 68
      jcon(kas)=jcon(kas)+1
      go to 33
66    jcon(kas)=jcon(kas)-1
      go to 68
69    jcon(las)=jcon(las)-1
      go to 74
75    jcon(lar)=jcon(lar)+1
      go to 72
73    jcon(kar)=jcon(kar)+1
      go to 71
581   if(mxex.eq.2) go to 22
cvp
cvp ? generation of configurations from the references
cvp call to bruna: might work incorrectly because jdeks isn't set in 
cvp  in bruna !!!
c     call caserr('call to bruna: possibly wrong results')
c     mfg - i think jdeks was being used to hold the higher
c     dimensionsed ideks .. try this ..
      call bruna(ideks,ndeks)
cvp
      go to 22
36    go to (3,42,48,57,66), nbd
25    nod=0
      ncl=0
      do 76 i=1,nbox
      jc=jcon(i)
      if (jc.eq.0) go to 76
      if (jc.eq.1) go to 77
      ncl=ncl+1
      iduf(ncl)=i
      go to 76
77    nod=nod+1
      icuf(nod)=i
76    continue
      nmns=nod/2+1-ispin
      if(nmns.le.0) goto 36
      ig=ig+1
      kcon(ig)=nmns
      nconfn(nmns)=nconfn(nmns)+1
      if (ig.lt.ibl) go to 85
        ig=0
        write(mdisk)kcon
 85   if (nod.eq.0) go to 80
      do 86 i=1,nod
      ig=ig+1
      kcon(ig)=icuf(i)
      if (ig.lt.ibl) go to 86
      ig=0
      write(mdisk) kcon
86    continue
      if (ncl.eq.0) go to 36
80    do 81 i=1,ncl
      ig=ig+1
      kcon(ig)=iduf(i)
      if (ig.lt.ibl) go to 81
      ig=0
      write(mdisk) kcon
81    continue
      go to 36
   23 write (iwr,8601)
cvp
cvp ? end of pseudo-loop
cvp
 8601 format(/1x,'no. of categories in each super-category'/)
      do j=1,ndkneu
        if(nconfn(j).ne.0)then
        write(iwr,9001)j,nconfn(j)
 9001   format(
     *  ' no. of configurations in super-category ',i2,' = ',i8)
        endif
      enddo
c     write(iwr,8498) nconfn
      write(mdisk) kcon
      if (.not.coreci) then
       call rewftn(mtape)
cvp if the number of open shells exceeds nopmax then this sk shouldn't
c    be taken into account any further.
      nodneu=nmul-3
      do 2013 kdoo=1,ndk5
        nodneu=nodneu+2
        nconf(kdoo)=nconfn(kdoo)
        if (nodneu.gt.nopmax) nconf(kdoo)=0
2013  continue
       write (mtape) nbox,m,nmul,ispace,nconf,nytl,mxex,jkon,nko,iswh,
     &     nyzl
      endif
      nod=nmul-3
c-- changes by hutter
c   change buffers
      if (coreci) then
       call rewftn(mtape)
       isave = mtape
       mtape = ltype
      endif
c--
cvp if the number of open shells exceeds nopmax then this sk shouldn't
c    be taken into account any further.
      nodneu=nmul-3
      do 2015 kdoo=1,ndk5
        nconf(kdoo)=nconfn(kdoo)
        nodneu=nodneu+2
        if (nodneu.gt.nopmax) nconf(kdoo)=0
 2015 continue
c-- loop over super categories
cvp loop just upto iswh=5, because only configurations upto the 5th sk
c    are taken into account
      do 100 iw=1,iswh
        nc=nconf(iw)
        nod=nod+2
        if (nc.eq.0) go to 100
        jg=0
        nl=nytln(iw)
        write (iwr,241) nc,nod
        call rewftn(mdisk)
        read(mdisk)kcon
        ig=0
        do 201 i=1,nc
202     ig=ig+1
        nad=kcon(ig)
        if (ig.lt.ibl) go to 203
        ig=0
        read(mdisk) kcon
203     if(nad.eq.iw) go to 220
        ip=nytln(nad)
        ig=ig+ip
        if(ig.lt.ibl) go to 202
        ig=ig-ibl
        read(mdisk)kcon
        go to 202
220     lx=ig+1
        ig=ig+nl
c-- changes by hutter
        if (coreci) then
c
c     make sure that buffer does not overflow
c
         if (ig.gt.ibl) then
          iig=0
          do 355 ii=lx,ibl
            iig=iig+1
  355     kchelp(iig)=kcon(ii)
          read (mdisk) kcon
          ig=ig-ibl
          do 356 ii=1,ig
           iig=iig+1
356       kchelp(iig)=kcon(ii)
          call singlep(kchelp(1),nl,m,nset,ndeg,inocc,nseto,iret)
          if (iret.eq.1) then
c$           call prk(kchelp(1),nl)
             nconf(iw) = nconf(iw)-1
c$           write(iwr,*) 'nconf(',iw,') = ',nconf(iw)
             iret = 0
             goto 201
          endif
          do 357 j=1,nl
            jg=jg+1
            mkon(jg)=kchelp(j)
            if (jg.lt.ibl) go to 357
            write(mtape)mkon
            jg=0
357       continue
          go to 201
         else
           call singlep(kcon(lx),nl,m,nset,ndeg,inocc,nseto,iret)
         endif
cvp
cvp ? end buffer-test
cvp
c-- if iret=1 remove configuration
        if (iret.eq.1) then
c$         call prk(kcon(lx),nl)
           nconf(iw) = nconf(iw)-1
c$         write(iwr,*) 'nconf(',iw,') = ',nconf(iw)
           if (ig.lt.ibl) goto 201
           read (mdisk) kcon
           ig=0
           goto 201
        endif
       endif
cvp
cvp ? end of ncorci-block
cvp
       if(ig.gt.ibl) go to 221
       do 222 j=lx,ig
       jg=jg+1
       mkon(jg)=kcon(j)
       if (jg.lt.ibl) go to 222
       write(mtape)mkon
       jg=0
222    continue
       if (ig.lt.ibl) go to 201
       ig=0
       read(mdisk)kcon
       go to 201
221    do 224 j=lx,ibl
          jg=jg+1
          mkon(jg)=kcon(j)
          if (jg.lt.ibl) go to 224
          write(mtape) mkon
          jg=0
224    continue
       ig=ig-ibl
       read(mdisk) kcon
       do 225 j=1,ig
          jg=jg+1
          mkon(jg)=kcon(j)
          if (jg.lt.ibl) go to 225
          write(mtape)mkon
          jg=0
225       continue
201       continue
       write(mtape)mkon
100   continue
cvp
cvp ? end of loop over the sk
c-- changes by hutter
c   bring cleaned configuration into old buffer
      if (coreci) then
       call rewftn(mtape)
       mtape = isave
       call rewftn(mtape)
       write (iwr,*) 'number of remaining configurations'
       nod1=nmul-3
       do 345 iw=1,iswh
        nod1=nod1+2
        if(nconf(iw).eq.0) goto 345
        write (iwr,241) nconf(iw),nod1
 345   continue
      write (mtape) nbox,m,nmul,ispace,nconf,nytl,mxex,jkon,nko,iswh,
     &    nyzl
 250   read(ltype,end=251) kcon
       write(mtape) kcon
       goto 250
 251   continue
       call rewftn(ltype)
      endif
      write(iwr,240)
240   format(/1x,'configuration generation complete'/)
241   format(10x,'transferred ',i9,' configurations with',i2,1x,
     +           'open shell(s)')
cvp
c  die alten felder werden werden wieder eingesetzt und mit den werten
c   der neu dimensionerten belegt
c  the old field are reused and set to the values of the newly 
c   dimensioned fields
      do 2010 kdoo=1,ndk5
        isw(kdoo)=iswn(kdoo)
 2010 continue
      do 2020 kdoo=1,nopmax
        jdeks(kdoo)=jdeksn(kdoo)
        kdeks(kdoo)=kdeksn(kdoo)
 2020 continue
      write(iwr,2033)
 2033 format(1x,'total number of configurations'/
     +       1x,'==============================')
      do kdoo=1,ndkneu
       write(iwr,2031) kdoo, nconfn(kdoo)
 2031  format(1x,'configurations in super category',i2,' = ',
     +         i8)
      enddo
cvp
      ird=restape
      return
766   write(iwr,760) nko,lko
760   format(5x,'too many mains',2i6)
      call caserr('too many reference configurations')
  767 write (iwr,768)
  768 format(5x, 'excitation class not allowed')
      call caserr('excitation class not allowed')
  779 write(iwr,780)
 780  format(5x,'requested number of electrons troublesome')
      call caserr('requested number of electrons troublesome')
 777  write(iwr,778)
 778  format(5x,'at least two of the main configurations are identical')
      call caserr('at least two identical main configurations')
 775  write(iwr,776)
 776  format(5x,'too many open shells in mains for dimensions')
      call caserr('too many open shells in mains for dimensions')
 773  write(iwr,774)
 774  format(5x,'orbital numbering is weird in mains')
      call caserr('orbital numbering is inconsistent in mains')
 771  write(iwr,772)
 772  format(5x,
     + 'open shell structure in mains inconsistent with multiplicity')
      call caserr(
     +  'open shell structure inconsistent with multiplicity')
 769  write(iwr,770) i,nv,nz,jkon(kg),kg,jkon(lb),lb
 770  format(5x,'pauli was right or maybe permutation error',7i5)
      call caserr('possible permutation error')
 783  write(iwr,784)
 784  format(5x,'multiplicity out of bounds')
      call caserr('multiplicity out of bounds')
 785  write(iwr,786)
 786  format(5x,'ispace parameter out of bounds')
      call caserr('ispace parameter out of bounds')
      return
      end
      subroutine fileconfpa
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
       common /ftap/ ntape, mtape, mdisk,
     +               ideli, ltype, linf,  ntab,  kfile,
     +               kclab, ntype, mstvt, nf01, nf62,
     +               nhead, nf99, mtapev, nf11
      integer maxroot
      parameter (maxroot=50)
      logical oconf,debugs
      common/parkin/egey1,trash,tdel,
     +       ical0,nele,nko,mxex,nmul,ispace,nprin,
     1       nft31,nstarv,ncorci,nform,oconf,debugs,
     +       lsng,maxci,ipt0,isym(maxroot)
_IF(notused)
      REAL ptiw
      common/junk/ptiw(5500),iclosed(256), iopen(256),ioff(8)
_ELSE
      common/junk/iclosed(256), iopen(256),ioff(8)
_ENDIF
      dimension parkconf(8)
      integer *4 ipark,mask,itest
      dimension ipark(16)
      logical odirect
      equivalence (ipark(1),parkconf(1))
      integer popmycnt
c let op 1024=maxorb!!
INCLUDE(common/newmrd_sort)
c read first record of ftn031
      integer mjl,kjl,ljl,njl,ntill,nball
      integer isyml,jabl,lsyml,ncompl
      integer ilecj2
c
      common/miscop/mjl(8),kjl(8),ljl(8),njl(8),ntill(8),
     +               nball(9),isyml(8),jabl(36),lsyml(800),
     +               ncompl(100),ilecj2(256,2)
c
      nf31=31
      call rewftn(nf31)
      read(nf31)n,nod,jod,nt,icmo,mjl,kjl,ljl,njl,nsel,
     +              ntill,nball,isyml,jabl,iorbs,knu,lsyml,
     +              ncompl,repel,etot,newlz
      call rewftn(nf31)
      call rewftn(nf01)
      read(nf01)odirect
      if (odirect) then
         do j=1,nko
            read(nf01) parkconf
            nopen=0
            iorbc=0
            nel=0
            iorb=0
            do i=16,1,-1
               mask=3
               do k=1,16
                   iorb=iorb+1
                   itest=iand(ipark(i),mask)
                   if (popmycnt(itest).eq.2) then
                      iorbc=iorbc+1
                      iclosed(iorbc)=newlz(iorb)
                   else if (popmycnt(itest).eq.1) then
                      nopen=nopen+1
                      iopen(nopen)=newlz(iorb)
                   end if
_IF(absoft,i8)
                   mask=ishft(mask,2)
_ELSE
                   mask=ISHFT(mask,2)
_ENDIF
               enddo
            enddo
            call mrdci_ssort(iclosed,iorbc)
            call mrdci_ssort(iopen,nopen)
            write(11,100)nopen,
     1      (iopen(ij),ij=1,nopen),(iclosed(ijk),ijk=1,iorbc)
         enddo
      else
         ioff(1)=0
         do i=2,8
            ioff(i)=ioff(i-1)+ljl(i-1)
         enddo
         do j=1,nko
            read(nf01)nopen,norbs,ilecj2
            do iiz=1,norbs
               iopen(iiz)=ilecj2(iiz,1)+ioff(ilecj2(iiz,2))
            enddo
            write(11,100)nopen,(iopen(ij),ij=1,norbs)
         enddo
      endif 
      return
 100  format (i4,',',/,256i4)
      end

      subroutine mrdci_ssort(ia,n)
      dimension ia(n)
      if (n.le.1) goto 20
      i=2
  10  continue
      if (ia(i-1).gt.ia(i)) then
         j=ia(i-1)
         ia(i-1)=ia(i)
         ia(i)=j
         if (i.gt.2) i=i-1
      else
         i=i+1
      endif
      if (i.eq.n+1) goto 20
      goto 10
  20  continue
      return
      end

      integer function popmycnt(itest)
      integer *4 itest,mask,m1
      m1 = 1
      mask=1
      ipop=0
      do i=1,32
         if (iand(itest,mask).ne.0)ipop=ipop+1
_IF(absoft,i8)
         mask=ishft(mask,m1)
_ELSE
         mask=ISHFT(mask,m1)
_ENDIF
      enddo
      popmycnt=ipop
      return
      end

      subroutine parkwa(q,debug,oprinsym)
      implicit REAL (a-h,o-z)
      REAL q
      logical debug,oprinsym
      dimension q(*)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c
      integer   nston,mtape,mdisk
      integer   ideli,ltype,linf,ntab,kfile,kclab
      integer   ntype, mstvt, nf01, nf62, nhead
      integer   nf99, mtapev, nf11
      common /ftap/ nston, mtape, mdisk,
     .              ideli, ltype, linf,  ntab,  kfile,
     .              kclab, ntype, mstvt, nf01, nf62,
     +              nhead, nf99, mtapev, nf11
c
      integer ihp, w, ie, if
      integer nt64, ndimf, maxref, mxref2, lsx
      integer nf31o
*  laenge der hamiltonmatrix fuer geyser-vektoren :
      parameter (maxref = 256)
      parameter (nt64 = maxref*maxref)
      parameter (mxref2 = maxref*(maxref+1)/2)
c
      parameter (maxshl = 50 )
      parameter (njkan = maxshl * maxref)
c
      parameter (ndimf = 2304)
*     parameter (ndimf = 48 * 48) !!!!
c
      parameter (lsx=132000)
*
*  anzahl der superkategorien
      parameter (ndk5 = 5 )
c
      integer niot, iot, nav, ot
      parameter (niot = 100 000)
c
      integer nt44
      parameter (nt44 =  44 000)
c
      parameter (nmmoo=256)
      parameter (nd8 = nmmoo*maxref)
c
c     now allocate integral arrays dynamically
c     determine space
c
      nav = lenwrd()
c
      nmmo = num
      ndeks = num*(num+1)/2 + 1
*  alte laenge des einelektronenrekords auf ft31
      nf31o = ndeks
      ndeks1 = ndeks-1
c
c hp, hy, a
c
      ihp = 1
      w = ihp + nt64
      ie = w + mxref2
      if = ie + ndimf
      isx = if + ndimf
      ipey = isx + lsx
      icoul = ipey + ndeks
      iexc = icoul + ndeks
      indeks= iexc + ndeks
      ipey1 = indeks + ndeks
      icoul1= ipey1 + nf31o
      iexc1= icoul1 + nf31o
      idelta = iexc1 + nf31o
      iseng = idelta + ndk5 * ndeks1
      iot = iseng + ndeks1
      ot = iot + niot / nav
      icon = ot + nt44 / nav
      jkan = icon + nd8 / nav
      need = jkan + njkan / nav
c
      ihp = igmem_alloc(need)
      w  = ihp + nt64
      ie = w + mxref2
      if = ie + ndimf
      isx = if + ndimf
      ipey = isx + lsx
      icoul = ipey + ndeks
      iexc = icoul + ndeks
      indeks= iexc + ndeks
      ipey1 = indeks + ndeks
      icoul1= ipey1+ nf31o
      iexc1= icoul1+ nf31o
      idelta = iexc1+  nf31o
      iseng = idelta + ndk5 * ndeks1
      iot = iseng + ndeks1
      ot = iot + niot / nav
      icon = ot + nt44 / nav
      jkan = icon + nd8 / nav
c
      call parkwa_dat
      call parkwa2(q,q(ihp),q(ihp),q(ihp),q(w),
     +             nt64,mxref2,
     +             q(ie),q(if),ndimf,q(isx),lsx,
     +             q(ipey),q(icoul),q(iexc),q(indeks),ndeks,
     +             q(ipey1),q(icoul1),q(iexc1),nf31o,
     +             q(idelta),q(iseng),ndeks1,q(iot),niot,
     +             q(ot),nt44,q(icon),nd8,q(jkan),njkan,maxref,
     +             oprinsym)
c
c      now free memory
c
      call gmem_free(ihp)
c
c     rewind units to be re-used by later modules
c
      call rewftn(mdisk)
      call rewftn(mtape)
      call rewftn(kfile)
      call rewftn(ntype)
c
      return
      end
      subroutine parkwa2(core,hy,hp,a,w,nt64,mxref2,e,f,ndimf,
     +                   sx,lsx,
     +                   pey,acoul,aexc,ideks,ndeks,
     +                   peyo,acoulo,aexco,nf31o,
     +                   sdelta,seng,lbuffs,
     +                   iot,niot,
     +                   ot,nt44,jcon,nd8,jkan,njkan,maxref,
     +                   oprinsym)
c
c -----------------------------------------------------
c
c dinge, welche zu erledigen sind :
c  -- eliminiere equivalenzen
c  -- input-format aendern
c  -- aritmetische ifs entfernen
c  -- loop struktur uebersichtlicher
c  -- goto's entfernen
c to do:
c  -- eliminate equivalences
c  -- change input-format
c  -- removing arithmic if's
c  -- more readable loop structure
c  -- remove goto's
c
c -----------------------------------------------------
c
c ------ evolved from 
c
c
c    the version of parkw that was send to Aachen
c
c          30.10.89 - store number and sk of selected konf      be
c                     on fi62
c                     will be used for wald,no write statements
c  sm77-6 - 28.9.89 - better modification of selection sheme    be
c  sm77-6 - 25.9.89 - if lsng<0 iselct is set to 1 and          be
c                     the selection procedure is modified
c                     see bambiu
c  sm77-5 - 21.9.89 - mhe(80) all single excitation with respect to be
c                     reference space are selected for lsng>0
c                     see geysez,rumpb
c  sm77-4 - 19.5.89 - nirm=99, ideks=4951
c  sm77-3 - 23.1.89 - additional goto2 in bambi to make merge work
c  sm77-2 - 18.12.88 - do not generate output if dim too large - bah
c  parkeu - fortran77 version - 23.10.87
c  mm=1..nroot bug not present in this version (destroyed mm if
c                                               lprint > 1)
c  hs(32)   bug in bambi not present (must read hs(80))
c  delc@    debug-version
c  delc&    activate timer routines
c $$$$$$$$$$$$$$$$$$$$$$$$ parkeub $$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit REAL (a-h,o-z)
c
      REAL hy,hp,a,w,e,f
      integer nt64, mxref2, ndimf
      dimension hy(nt64),hp(nt64),a(nt64)
      dimension w(mxref2)
      dimension e(ndimf), f(ndimf),core(*)
c
      REAL sx
      integer lsx
      dimension sx(lsx)
c
      REAL pey,acoul,aexc
      integer ndeks,ideks
      dimension pey(ndeks),acoul(ndeks),aexc(ndeks)
      dimension ideks(ndeks)
c
      REAL peyo,acoulo,aexco
      integer nf31o
      dimension peyo(nf31o),acoulo(nf31o),aexco(nf31o)
c
      REAL sdelta, seng
      integer lbuffs
      dimension sdelta(5,lbuffs), seng(lbuffs)
c
      integer iot
      dimension iot(niot)
      integer ot
      dimension ot(nt44)
      integer jcon, nd8
      dimension jcon(nd8)
      integer jkan, njkan
      dimension jkan(njkan)
c
      logical oprinsym
c
INCLUDE(common/iofile)
c
      integer  ical0,nft31
c
      integer   nzero
      parameter (nzero = 0)
      parameter (mxroot=50)
c
clear character*8 idprog,idfi32,idvers,iddate,idtime
c
      REAL egey1,trash,tdel
      integer ndum1,ndum2,nstarv
      common /cegey1/ egey1,trash,tdel,ndum1,ndum2,nstarv
      integer   nform
      common /cffor/  nform
c
      integer   nston,mtape,mdisk
      integer   ideli,ltype,linf,ntab,kfile,kclab
      integer   ntype, mstvt, nf01, nf62, nhead
      integer   nf99, mtapev, nf11
      common /ftap/ nston, mtape, mdisk,
     .              ideli, ltype, linf,  ntab,  kfile,
     .              kclab, ntype, mstvt, nf01, nf62,
     +              nhead, nf99, mtapev, nf11
c
      integer   ical,nele,nko,mxex,nmul,ispace,nprin
      integer   ncorci, lsng, maxci2, ipt0, isymm
      logical   debugs
      common /tap/  ical,  nele,  nko,
     .              mxex,  nmul,  ispace,nprin, ncorci,
     +              lsng,
     +              maxci2, ipt0,  isymm(mxroot),
     +              debugs
csut namelist - input
      REAL cpulft
c
      REAL egey1in, trashin, tdelin
      integer ical0in,nelein,nkoin,mxexin,nmulin,ispacein
      integer nprinin,nft31in,nstarvin
      integer ncorciin,nformin,maxciin
      integer ipt0in,isymin
      logical oconf,debugin
      common /parkin/ egey1in,trashin,tdelin,
     +   ical0in,nelein,nkoin,mxexin,nmulin,ispacein,
     +   nprinin,nft31in,nstarvin,ncorciin,nformin,oconf,
     +   debugin,lsngin,
     +   maxciin,ipt0in,isymin(mxroot)
      character*10 charwall

       ndum1  = nzero
       ndum2  = nzero
csut  set default-values
       egey1  = -999.0d7
       ical0  = ical0in
       nele   = nelein
       nko    = nkoin
       mxex   = mxexin
       nmul   = nmulin
       nstarv = nstarvin
       nform  = nformin
       ispace = ispacein
       nprin  = nprinin
       ncorci = ncorciin
       maxci2 = maxciin
       trash=trashin
       tdel=tdelin
       lsng = lsngin
       nft31  = 1
       ipt0 = ipt0in
       do i=1,mxroot
        isymm(i) = isymin(i)
       enddo
       debugs = debugin
c
       write(iwr,11)
 11    format(/1x,104('=')//30x,47('*')/
     + 30x,'*** MRD-CI V2.0: Configuration Selection Module'/
     + 30x,47('*')/)
  12   format(/5x,
     + '***  start of selection at ',f10.2,' seconds',a10,' wall'/)
       write (iwr,12) cpulft(1) ,charwall()
c
cclear       read (ird,1) ical,m,nko,mxex,nmul,ispace,nprin
c-- changes by hutter: read new input variable ncorci
c   ncorci = 0 => no deleting of configuration
c   ncorci = 1 => deleting of configurations with "from"-excitation of
c--               more than one electron
c   ncorci = 2 => deleting of configurations with "in"-excitation of
c--               more than one electron
c   ncorci = 3 => deleting of configurations with other then given
c--               occupation in chosen orbitals
cc aenderung: k.pfingst 8/10/92
c--               occupation in chosen orbitals
c   ncorci = 4 => combination of 1,2,3
c--
cclear   1   format(14i5)
c ---- namelist - input
c      read(ird,parkin,end=20)
c  20  continue
c ---- ende des namelist - inputs
       ical=ical0
c ---- print input-daten
  15   format(/104('-'))
  16   format(A37,6x,I10)
c 17   format(A37,5x,F16.4)
       write(iwr,15)
c      write(iwr,16)'ical                               :',ical
       write(iwr,16)'number of electrons                :',nele
       write(iwr,16)'number of reference configurations :',nko
       write(iwr,16)'excitation level                   :',mxex
       write(iwr,16)'multiplicity                       :',nmul
       write(iwr,16)'ncorci                             :',ncorci
       write(iwr,16)'spatial symmetry                   :',ispace
c      write(iwr,16)'lese von stoneyfile                :',nft31
c      write(iwr,17)'egey1 (geyser1-raum)               :',egey1
       write(iwr,9300) trash,tdel
 9300  format(/' *** threshold specified ***'/
     + ' minimal selection threshold ',f7.2,' microhartree'/
     + ' threshold increment for use in selection ',f7.2,
     + ' microhartree')
       write(iwr,15)

       if (nstarv.gt.nzero) then
          write(iwr,'(A)')'  *** zero-order vector as start vector:'
       end if
c
       if (nft31.eq.1) then
csut     write(iwr,*)'copy ft31-stoney into new files'
csut     write(iwr,*)'ft71 : ......'
csut     write(iwr,*)'ft72 : ......'
csut     write(iwr,*)'ft73 : ......'
         call ft31sp
         write(iwr,21)cpulft(1) ,charwall()
   21    format(/,5x,'***  end of reading ft31 at ',f10.2,' seconds',
     +          a10,' wall')
         write(iwr,15)
       end if
       if (nstarv.gt.nzero) then
cmvs
c         open(unit=mstvt,file='gstarv.dat',form='unformatted',
c    +    status='unknown')
c         write(iwr,*)'file gstarv.dat must be opened'
cmvs
cmvs      rewind (mstvt)
cmvs
       end if
       if (ical.eq.0) then
         ikon = igmem_alloc(nd8)
         call parke4(ideks,ndeks,core(ikon),nd8,maxref,oprinsym)
         call gmem_free(ikon)
         write(iwr,18)cpulft(1) ,charwall()
   18    format(/5x,'***  end of configuration generation at ',
     +          f10.2,' seconds',a10,' wall')
         write(iwr,15)
       end if
c
_IFN(cray)
       call setsto(nt44,0,ot)
_ENDIF
c
       mvect = maxref * mxroot
       mvect2 = mvect + mvect
       mvect3 = mvect2 + mvect
       mvect4 = mvect3 + mvect
       ivect = igmem_alloc(mvect4)
       call geysez(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +             sx,lsx,
     +             pey,acoul,aexc,ideks,ndeks,
     +             peyo,acoulo,aexco,nf31o,
     +             sdelta,seng,lbuffs,
     +             core(ivect),core(ivect+mvect),
     +             core(ivect+mvect2),mvect,
     +             iot,niot,
     +             ot,nt44,jcon,nd8,jkan,njkan,
     +             core(ivect+mvect3))
       if (egey1.gt.-998.0d0) then
          write(iwr,15)
          write(iwr,19)
   19     format(/5x,'*** guess according to larger space'/
     +            5x,'==================================='/)
          call coolez(hy,hp,a,nt64,w,mxref2,e,f,ndimf,
     +                sx,lsx,
     +                pey,acoul,aexc,ideks,ndeks,
     +                peyo,acoulo,aexco,nf31o,
     +                sdelta,seng,lbuffs,
     +                core(ivect),core(ivect+mvect),
     +                core(ivect+mvect2),mvect,
     +                iot,niot,
     +                ot,nt44,jcon,nd8,jkan,njkan,
     +                core(ivect+mvect3))
       end if
c
       call gmem_free(ivect)
c
c      remove large temporary files
c
       call shutftn(kclab)
       call shutftn(mtape)
       call shutftn(mdisk)

       write(iwr,22)cpulft(1) ,charwall()
   22  format(/5x,'***  end of configuration selection at ',f10.2,
     +         ' seconds',a10,' wall')
       return
       end
c
c-----------------------------------------------
      subroutine parkwa_dat
c
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
c     data jerk /0,29,116,155,244,270,316,0,0,0/
c     data jbnk /3,16,6,14,3,4,6,0,0,0/
c     data jbun /3,3,2,2,1,1,1,0,0,0/
c
      do loop =1,10
       jerk(loop) = 0
       jbnk(loop) = 0
       jbun(loop) = 0
      enddo
c
_IF(notused)
      jerk (2)  = 29
      jerk (3)  = 116
      jerk (4)  = 155
      jerk (5)  = 244
      jerk (6)  = 270
      jerk (7)  = 316
c
      jbnk (1)  = 3
      jbnk (2)  = 16
      jbnk (3)  = 6
      jbnk (4)  = 14
      jbnk (5)  = 3
      jbnk (6)  = 4
      jbnk (7)  = 6
_ELSE
      jerk (2)  = 30
      jerk (3)  = 118
      jerk (4)  = 158
      jerk (5)  = 248
      jerk (6)  = 275
      jerk (7)  = 322
c
      jbnk (1)  = 4
      jbnk (2)  = 17
      jbnk (3)  = 7
      jbnk (4)  = 15
      jbnk (5)  = 4
      jbnk (6)  = 5
      jbnk (7)  = 7
_ENDIF
c
      jbun (1) = 3
      jbun (2) = 3
      jbun (3) = 2
      jbun (4) = 2
      jbun (5) = 1
      jbun (6) = 1
      jbun (7) = 1
c
      return
      end
c-----------------------------------------------
      subroutine rumam(iot,niott,ideks,ndeks,jcon,nd8)
c
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndeks,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
c
      parameter (nmmo=256)
      integer lj
      common /clj2/ lj(8)
c
      common /cnit/ nit(nitmax)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
      integer kc, kd, lab
      common /ckd/  kc(ndimh),kd(ndimh),lab(3)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,mt2,
     +md2,nston,mdisk,ideli,mtape,iswh,mr2,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref), ie(maxref), nsc(maxref),
     +imo,mc,nzx,md,mr,mx,mt,ih2,iht,ihz,jh2,ny,mz,jy,
     +jx,nd1,nd2,md1,nd,nr3,nr4,mr1,nr,nt1,nt2,mt1,nnx,ih1,nt,ih3,
     +ih4,nz,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,ntx,ndx,nrx,
     +ipag
c
      l2=jy
      l3=nz
      l4=nt
      do 13 j=1,mc
      jz=mz
      do 19 l=1,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l2)-1) 20,21,19
   21 ia=l+1
      go to 22
   19 continue
   20 mm=mm+1
      olab(mm)=nzero
      if (mm.lt.nnid) go to 61
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=nzero
      go to 61
   22 do 23 l=ia,md
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn+l2)-1) 20,24,23
   24 ja=l+1
      go to 25
   23 continue
   25 if (ja.gt.md) go to 26
      do 9964 l=ja,md
        jz=jz+1
        ii=iot(jz)
        if (jcon(ii+l2).ne.2) go to 20
9964  continue
   26 if (mr.eq.0) go to 28
      jz=mt
      do 29 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii+l2).ne.1) go to 20
   29 continue
   28 jz=l4
      kz=mt+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 30 l=1,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.ll) go to 31
      if (jq.eq.mr.or.iq.ne.ii) go to 32
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   30 continue
   31 ip=l
      ia=l+1
      do 33 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.nn) go to 34
      if (jq.eq.mr.or.iq.ne.ii) go to 35
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   33 continue
   32 ip=l
      ia=l+1
      do 36 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.ll) go to 37
      if (jq.eq.mr.or.iq.ne.jj) go to 38
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   36 continue
   34 jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   39 continue
   40 kp=jdeks(l-2)
      ia=l+1
      do 41 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.iq.ne.jj) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   41 continue
   42 mm=mm+1
      olab(mm)=1
      lp=kdeks(l-3)
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   35 jp=ideks(l-1)
      ia=l+1
      do 44l=ia,nr
       jz=jz+1
       jj=iot(jz)
       if (jj.eq.nn) go to 45
       if (jq.eq.mr.or.jj.ne.iq) go to 46
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
   44 continue
   37 jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
       jz=jz+1
       jj=iot(jz)
       if (jj.eq.nn) go to 48
       if (jq.eq.mr.or.jj.ne.iq) go to 49
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.ll) go to 51
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp=jdeks(l-2)
      ia=l+1
      do 53 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
cx    mmr = nzero
      mm = nzero
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nnid) go to 170
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 if (icq.lt.idq) go to 134
      icq = ideks(icq)+idq
      go to 129
  134 icq = ideks(idq)+icq
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 146
      iay=ideks(nl1r)+nm1r
      go to 147
  146 iay= ideks(nm1r)+nl1r
  147 if (nl2r.lt.nm2r) go to 148
      iby=ideks(nl2r)+nm2r
      go to 149
  148 iby=ideks(nm2r)+nl2r
  149 if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      if (iay.lt.iby) go to 153
      icq=ideks(iay)+iby
      go to 129
  153 icq=ideks(iby)+iay
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
*     write(6,*) 'rumam: kc,kd=', kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 l4=l4+nnx
      l3=l3+nnx
  13  l2=l2+imo
      return
      end
      subroutine rumap(iot,niott,ideks,ndeks,jcon,nd8)
      implicit REAL (a-h,o-z)
c
c
      integer iot,niott,ideks,ndeks,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
c
      integer lj
      common /clj2/ lj(8)
c
      common /cnit/ nit(nitmax)
      integer kc,kd,lab
      common /ckd/  kc(ndimh),kd(ndimh),lab(3)
      common /a/ t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     +ndx,nston,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref),ie(maxref),nsc(maxref),
     +imo,mc,jh4,nd,nr,nnx,nt,ny,iht,ihz,nzx,ih4,nz,jy,
     +jx,nd1,md,md1,md2,nr3,mr,mr1,mr2,nt1,mx,mt1,mt2,ih1,ih2,ih3,
     +mt,jh2,jh3,mz,iab,nt3,iac,nr2,m,nad,nad1,ifr1,nt2,nd2,nr4
c
      l3=mz
      l4=mt
      do 17 k=1,mc
      jz=l3
      do 19 l=1,md
       jz=jz+1
       ll=iot(jz)
       if (jcon(ll)-1) 20,21,19
   21  ia=l+1
       go to 22
   19 continue
   20 mm=mm+1
      olab(mm)=nzero
      if (mm.lt.nnid) go to 61
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 61
   22 do 23 l=ia,md
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn)-1) 20,24,23
   24 ja=l+1
      go to 25
   23 continue
   25 if (ja.gt.md) go to 26
      do 27 l=ja,md
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.2) go to 20
   27 continue
   26 if (mr.eq.0) go to 28
      jz=l4
      do 29 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.1) go to 20
   29 continue
   28 jz=nt
      kz=l4+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 30 l=1,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.ll) go to 31
      if (jq.eq.mr.or.iq.ne.ii) go to 32
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   30 continue
   31 ip=l
      ia=l+1
      do 33 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (ii.eq.nn) go to 34
      if (jq.eq.mr.or.iq.ne.ii) go to 35
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   33 continue
   32 ip=l
      ia=l+1
      do 36 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.ll) go to 37
      if (jq.eq.mr.or.iq.ne.jj) go to 38
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   36 continue
   34 jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   39 continue
   40 kp=jdeks(l-2)
      ia=l+1
      do 41 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.iq.ne.jj) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   41 continue
   42 mm=mm+1
      olab(mm)=1
      lp=kdeks(l-3)
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   35 jp=ideks(l-1)
      ia=l+1
      do 44l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 45
      if (jq.eq.mr.or.jj.ne.iq) go to 46
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   44 continue
   37 jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jj.eq.nn) go to 48
      if (jq.eq.mr.or.jj.ne.iq) go to 49
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.ll) go to 51
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 42
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp=jdeks(l-2)
      ia=l+1
      do 53 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nnid) go to 170
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 if (icq.lt.idq) go to 134
      icq = ideks(icq)+idq
      go to 129
  134 icq = ideks(idq)+icq
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 146
      iay=ideks(nl1r)+nm1r
      go to 147
  146 iay= ideks(nm1r)+nl1r
  147 if (nl2r.lt.nm2r) go to 148
      iby=ideks(nl2r)+nm2r
      go to 149
  148 iby=ideks(nm2r)+nl2r
  149 if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      if (iay.lt.iby) go to 153
      icq=ideks(iay)+iby
      go to 129
  153 icq=ideks(iby)+iay
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 l3=l3+mx
   17 l4=l4+mx
      return
      end
      subroutine rumbm(iot,niott,ideks,ndeks,jcon,nd8)
      implicit REAL (a-h,o-z)
c
c
      integer iot,niott,ideks,ndek,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
c
      integer lj
      common /clj2/ lj(8)
c
      common /cnit/ nit(nitmax)
      integer kc,kd,lab
      common /ckd/  kc(ndimh),kd(ndimh),lab(3)
      common /a/ t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer ie, nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,mt1,
     +md1,nston,mdisk,ideli,mtape,iswh,mr1,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref),ie(maxref),nsc(maxref),
     +imo,mc,nzx,md,mr,mx,mt,ih1,iht,ihz,nz,ny,mz,jx,
     +jy,nd1,nd2,nd,md2,nr3,nr4,nr,mr2,nt1,nt2,nnx,mt2,nt,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,ntx,ndx,nrx
     +,ipag
c
c
      l3=jy
      l4=nt
      do 65 j=1,mc
      nt4=l4+1
      nz=l4+nr
      jz=mz
      do 70 l=1,md
       jz=jz+1
       jj=iot(jz)
       if (jcon(jj+l3)-1) 71,72,70
   70 continue
   71 if (l.eq.md) go to 73
      ia = l+1
      do 74 l=ia,md
       jz = jz + 1
       ll = iot(jz)
       if (jcon(ll+l3).ne.2) go to 75
   74 continue
      go to 73
  75  mm = mm + 1
      olab(mm) = 1
      if (mm.lt.nnid) go to 260
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm = nzero
      go to 260
   72 l1 = l
      if (l.eq.md) go to 76
      ia = l + 1
      do 78 l=ia,md
       jz = jz + 1
       ll = iot(jz)
       if (jcon(ll+l3)-1) 75,77,78
  78  continue
      go to 76
  77  l2 = l
      if (l.eq.md) go to 79
      ia = l + 1
      do 80 l=ia,md
       jz = jz + 1
       kk = iot(jz)
       if (jcon(kk+l3).ne.2) go to 75
   80 continue
      go to 79
   73 if (mr.eq.0) go to 81
      jz = mt
      do 82 l=1,mr
       jz = jz + 1
       ll = iot(jz)
       if (jcon(ll+l3).ne.1) go to 75
   82 continue
  81  mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 83
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm = nzero
  83  jz = l4
      kz = mt + 1
      jq = 0
      if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jq.eq.mr.or.iq.ne.ll) go to 85
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia=l+1
      ip=l
      do 86 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jq.eq.mr.or.iq.ne.nn) go to 87
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   86 continue
  87  mm=mm+1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nnid) go to 88
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm = nzero
  88  nt1r=nir(jj)
      nl1r=loc(jj)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 if (nl1r.lt.nm2r) go to 248
      iby=ideks(nl1r)+nm2r
      go to 253
  248 iby=ideks(nm2r)+nl1r
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=mt
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii+l3)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc+l3).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 94
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   94 jz=l4
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nnid) go to 99
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=mz
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=mz+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=mz+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 187
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 187  jz=l4
      kz=mt
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=mr
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+l4
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nnid) go to 196
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nnid) go to 197
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nnid) go to 198
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+l4
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+l4
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 if (icq.lt.idq) go to 334
      icq = ideks(icq)+idq
      go to 329
  334 icq = ideks(idq)+icq
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 346
      iay=ideks(nl1r)+nm1r
      go to 347
  346 iay= ideks(nm1r)+nl1r
  347 if (nl2r.lt.nm2r) go to 348
      iby=ideks(nl2r)+nm2r
      go to 349
  348 iby=ideks(nm2r)+nl2r
  349 if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      if (iay.lt.iby) go to 353
      icq=ideks(iay)+iby
      go to 329
  353 icq=ideks(iby)+iay
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 320
  76  if (mr.eq.0) go to 273
      jz=mt
      do 274 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll+l3)-1) 275,274,75
 274  continue
      go to 273
 275  l1=l
      if (l.eq.mr) go to 277
      ia=l+1
      do 276 l=ia,mr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn+l3).ne.1) go to 75
  276 continue
 277  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 278
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 278  jz=l4
      kz=mt
      ig=0
      if (l1.eq.1) go to 279
      ja=l1-1
      do 280 l=1,ja
      kz=kz+1
      jp=iot(kz)
 254  jz=jz+1
      if (jp.eq.iot(jz)) go to 280
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 281
      go to 254
 280  continue
 279  if (l1.eq.mr) go to 282
      ia=l1+1
      kz=kz+1
      do 283 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 255  jz=jz+1
      if (jp.eq.iot(jz)) go to 283
      ig=ig+1
      lab(ig)=jz-l4
      if (ig.eq.3) go to 281
      go to 255
  283 continue
  282 lab(3)=nr
      if (ig-1) 286,287,281
  286 lab(2)=nr1
      lab(1)=mr
      go to 281
  287 lab(2)=nr1
  281 jz=lab(1)+l4
      kk=iot(jz)
      mm=mm+1
      if (kk.ne.jj) go to 288
      olab(mm)=-1
      jz=lab(2)+l4
      kk=iot(jz)
      jz=lab(3)+l4
      nn=iot(jz)
 289  if(mm.lt.nnid) go to 290
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 290  mm=mm+1
      olab(mm)=l1
      if(mm.lt.nnid) go to 291
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 291  mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp) +jdeks(kp)
      if(mm.lt.nnid)go to 292
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 292
 288  jz=lab(2)+l4
      nn=iot(jz)
      if (nn.ne.jj) go to 293
      olab(mm)=-2
      jz=lab(3)+l4
      nn=iot(jz)
      go to 289
 293  olab(mm)=-3
      go to 289
 292  nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(nn)
      nm1r=loc(kk)
      nm2r=loc(nn)
      go to 320
 273  mm=mm+1
      olab(mm)=4
      if(mm.lt.nnid) go  to 294
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 294  jz=l4
      mm=mm+1
      if(mr.eq.0) go to 297
      kz=mt
      do 295 l=1,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 296
 295  continue
 297  ll=iot(jz+1)
      if(ll.eq.jj) go to 298
      olab(mm) =iab
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 298  ll=iot(jz+2)
      olab(mm)=-iab
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 296  l1 =jz-l4
      if(ll.eq.jj) go to 373
      do 374 l=1,nr
      jz=jz+1
      if(iot(jz).eq.jj) go to 375
 374  continue
 375  olab(mm)=l1+ideks(jz-nt4)
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 373  kz=kz-1
      do 376 ia=l,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 377
 376  continue
      jz=jz+1
      ll=iot(jz)
 377  olab(mm)=-l1-ideks(jz-nt4)
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 371  mm=mm+1
      iax=nir(jj)
      olab(mm)=iax
      if(mm.lt.nnid) go to 372
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  372 iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(jj)
      kk=loc(ll)
      km=ideks(kk)
      nm=ideks(ii)
      nn=nm+ii
      if(ii.gt.kk) go to 380
      mm=mm+1
      jb=ii+km
      olab(mm)=kk
      if(mm.lt.nnid) go to 378
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 378  mm=mm+1
      olab(mm)=ii
      if(mm.lt.nnid) go to 379
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 379  ib=ideks(jb) +nn+iat
      go to 381
 380  lb=kk-ii
      jb=nn+lb
      ib=ideks(nn+1)+lb+iat
      mm=mm+1
      olab(mm)=ii
      if(mm.lt.nnid) go to 900
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 900  mm=mm+1
      olab(mm)=kk
      if(mm.lt.nnid) go to 381
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 381  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 382
      jto=0
      write(ntype) kc,kd
      jblk=jblk+1
 382  if(mr.eq.0) go to 383
      jz=mt
      ir=mr
      kix=1
 472  do 384 l=1,ir
      jz=jz+1
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  if(kb.lt.lb) go to 393
      ib=iat+ideks(kb) +lb
      go to 394
 393  ib=iat+ideks(lb) +kb
 394  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      if(kb.lt.lb) go to 398
      ib=ideks(kb)+lb+ibx
      go to 471
 398  ib=ideks(lb) +kb+ibx
 471  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto = nzero
      jblk = jblk + 1
      write(ntype) kc,kd
 473  kb = (mlr-1)*mq
      lb = kb + ii
      kb = kb + kk
      if(kb.lt.lb) go to 474
      ib=ideks(kb)+lb+ibx
      go to 475
 474  ib=ideks(lb)+kb+ibx
 475  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 260
 383  if (nd.eq.0) go to 260
      kix=-1
      jz=nz
      ir=nd
      go to 472
 260  l4=l4+nnx
  65  l3=l3+imo
      return
      end
      subroutine rumbp(iot,niott,ideks,ndeks,jcon,nd8)
c
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndeks,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
c
      integer lj
      common /clj2/ lj(8)
c
      integer kc,kd,lab
      common /ckd/ kc(ndimh),kd(ndimh),lab(3)
      common /cnit/ nit(nitmax)
      common /a/ t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,ntx,
     +ndx,nston,mdisk,ideli,mtape,iswh,nrx,mm,jblk,jto,
     +igmax,nr2,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref),ie(maxref),nsc(maxref),
     +imo,mc,jh3,nd,nr,nnx,nt,ny,iht,ihz,nzx,ih3,nz,jy,
     +jx,md,nd2,md1,md2,mr,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,mt,
     +ih4,jh2,mz,jh4,iac,mx,iab,nr1,m,nad,nad1,ifr1,nt3,nd1,nr3,if
c
      l3=mz
      l4=mt
      do 69 k=1,mc
      jz=l3
      do 70 l=1,md
       jz=jz+1
       jj=iot(jz)
       if (jcon(jj)-1) 71,72,70
   70 continue
   71 if (l.eq.md) go to 73
      ia=l+1
      do 1579 l=ia,md
       jz = jz + 1
       ll = iot(jz)
       if (jcon(ll).ne.2) go to 75
1579  continue
      go to 73
  75  mm = mm + 1
      olab(mm) = 1
      if (mm.lt.nnid) go to 260
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm = nzero
      go to 260
   72 l1=l
      if (l.eq.md) go to 76
      ia=l+1
      do 78 l=ia,md
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 75,77,78
  78  continue
      go to 76
  77  l2=l
      if (l.eq.md) go to 79
      ia=l+1
      do 80 l=ia,md
      jz=jz+1
      kk=iot(jz)
      if (jcon(kk).ne.2) go to 75
   80 continue
      go to 79
   73 if (mr.eq.0) go to 81
      jz=l4
      do 82 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll).ne.1) go to 75
  82  continue
  81  mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 83
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
  83  jz=nt
      kz=l4+1
      jq=0
      if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
      jz=jz+1
      ll=iot(jz)
      if (jq.eq.mr.or.iq.ne.ll) go to 85
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia=l+1
      ip=l
      do 86 l=ia,nr
      jz=jz+1
      nn=iot(jz)
      if (jq.eq.mr.or.iq.ne.nn) go to 87
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   86 continue
  87  mm=mm+1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nnid) go to 88
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
  88  nt1r=nir(jj)
      nl1r=loc(jj)
      nb1r=nir(ll)
      nm1r=loc(ll)
      nm2r=loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 if (nl1r.lt.nm2r) go to 248
      iby=ideks(nl1r)+nm2r
      go to 253
  248 iby=ideks(nm2r)+nl1r
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=l4
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 94
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   94 jz=nt
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nnid) go to 99
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=l3
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=l3+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=l3+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 187
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 187  jz=nt
      kz=l4
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=mr
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+nt
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nnid) go to 196
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nnid) go to 197
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nnid) go to 198
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+nt
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+nt
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 if (icq.lt.idq) go to 334
      icq = ideks(icq)+idq
      go to 329
  334 icq = ideks(idq)+icq
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 346
      iay=ideks(nl1r)+nm1r
      go to 347
  346 iay= ideks(nm1r)+nl1r
  347 if (nl2r.lt.nm2r) go to 348
      iby=ideks(nl2r)+nm2r
      go to 349
  348 iby=ideks(nm2r)+nl2r
  349 if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      if (iay.lt.iby) go to 353
      icq=ideks(iay)+iby
      go to 329
  353 icq=ideks(iby)+iay
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 320
  76  if (mr.eq.0) go to 273
      jz=l4
      do 274 l=1,mr
      jz=jz+1
      ll=iot(jz)
      if (jcon(ll)-1) 275,274,75
 274  continue
      go to 273
 275  l1=l
      if (l.eq.mr) go to 277
      ia=l+1
      do 276 l=ia,mr
      jz=jz+1
      nn=iot(jz)
      if (jcon(nn).ne.1) go to 75
  276 continue
 277  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 278
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 278  jz=nt
      kz=l4
      ig=0
      if (l1.eq.1) go to 279
      ja=l1-1
      do 280 l=1,ja
      kz=kz+1
      jp=iot(kz)
 254  jz=jz+1
      if (jp.eq.iot(jz)) go to 280
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 254
 280  continue
 279  if (l1.eq.mr) go to 282
      ia=l1+1
      kz=kz+1
      do 283 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 255  jz=jz+1
      if (jp.eq.iot(jz)) go to 283
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 281
      go to 255
  283 continue
  282 lab(3)=nr
      if (ig-1) 286,287,281
  286 lab(2)=nr1
      lab(1)=mr
      go to 281
  287 lab(2)=nr1
  281 jz=lab(1)+nt
      kk=iot(jz)
      mm=mm+1
      if (kk.ne.jj) go to 288
      olab(mm)=-1
      jz=lab(2)+nt
      kk=iot(jz)
      jz=lab(3)+nt
      nn=iot(jz)
 289  if(mm.lt.nnid) go to 290
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 290  mm=mm+1
      olab(mm)=l1
      if(mm.lt.nnid) go to 291
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 291  mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp) +jdeks(kp)
      if(mm.lt.nnid)go to 292
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 292
 288  jz=lab(2)+nt
      nn=iot(jz)
      if (nn.ne.jj) go to 293
      olab(mm)=-2
      jz=lab(3)+nt
      nn=iot(jz)
      go to 289
 293  olab(mm)=-3
      go to 289
 292  nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(kk)
      nb2r=nir(nn)
      nm1r=loc(kk)
      nm2r=loc(nn)
      go to 320
 273  mm=mm+1
      olab(mm)=4
      if(mm.lt.nnid) go  to 294
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 294  jz=nt
      mm=mm+1
      if(mr.eq.0) go to 297
      kz=l4
      do 295 l=1,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 296
 295  continue
 297  ll=iot(jz+1)
      if(ll.eq.jj) go to 298
      olab(mm) =iab
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 298  ll=iot(jz+2)
      olab(mm)=-iab
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 296  l1 =jz-nt
      if(ll.eq.jj) go to 373
      do 374 l=1,nr
      jz=jz+1
      if(iot(jz).eq.jj) go to 375
 374  continue
 375  olab(mm)=l1+ideks(jz-nt1)
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 371
 373  kz=kz-1
      do 376 ia=l,mr
      kz=kz+1
      jz=jz+1
      ll=iot(jz)
      if(ll.ne.iot(kz)) go to 377
 376  continue
      jz=jz+1
      ll=iot(jz)
 377  olab(mm)=-l1-ideks(jz-nt1)
      if(mm.lt.nnid) go to 371
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 371  mm=mm+1
      iax=nir(jj)
      olab(mm)=iax
      if(mm.lt.nnid) go to 372
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  372 iay=ideks(iax+1)
      iaz=ideks(iay+1)
      mq=lj(iax)
      mv=ideks(mq+1)
      iaq=iay-iax
      iaw=iaz-iay
      iat=nit(iaz)
      ii=loc(jj)
      kk=loc(ll)
      km=ideks(kk)
      nm=ideks(ii)
      nn=nm+ii
      if(ii.gt.kk) go to 380
      mm=mm+1
      jb=ii+km
      olab(mm)=kk
      if(mm.lt.nnid) go to 378
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 378  mm=mm+1
      olab(mm)=ii
      if(mm.lt.nnid) go to 379
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 379  ib=ideks(jb) +nn+iat
      go to 381
 380  lb=kk-ii
      jb=nn+lb
      ib=ideks(nn+1)+lb+iat
      mm=mm+1
      olab(mm)=ii
      if(mm.lt.nnid) go to 900
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 900  mm=mm+1
      olab(mm)=kk
      if(mm.lt.nnid) go to 381
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 381  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 382
      jto=0
      write(ntype) kc,kd
      jblk=jblk+1
 382  if(mr.eq.0) go to 383
      jz=l4
      ir=mr
      kix=1
 472  do 384 l=1,ir
      jz=jz+1
      mi=iot(jz)
      mar=nir(mi)
      mlr=loc(mi)
      mlb=ideks(mlr)
      mla=mlb+mlr
      if (mar-iax) 385,395,396
 395  if(mla.lt.jb) go to 386
      ib=iat+ideks(mla) +jb
 388  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto = jto+1
      kc(jto)=idud+1
      kd(jto) =ib
      if(jto.lt.iwod) go to 387
      jto =0
      write(ntype) kc,kd
      jblk=jblk+1
      go to 387
 386  ib=iat+ideks(jb) +mla
      go to 388
 387  if(mlr.lt.ii) go to 389
      kb=mlb+ii
      go to 390
 389  kb=nm+mlr
 390  if(mlr.lt.kk) go to 391
      lb=mlb+kk
      go to 392
 391  lb=mlr+km
 392  if(kb.lt.lb) go to 393
      ib=iat+ideks(kb) +lb
      go to 394
 393  ib=iat+ideks(lb) +kb
 394  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
 385  iby=ideks(mar+1)
      iby=iaw+iby
      ibx=iaq+mar
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      iby=nit(iby)
      ml=lj(mar)
      ib=ideks(ml+1)*(jb-1) +mla+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 397
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 397  kb=(ii-1)*ml +mlr
      lb=(kk-1)*ml +mlr
      if(kb.lt.lb) go to 398
      ib=ideks(kb)+lb+ibx
      go to 471
 398  ib=ideks(lb) +kb+ibx
 471  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto) = ib
      if(jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
      go to 384
396   ibx=ideks(mar+1)
      iby=ideks(ibx) +iay
      iby=nit(iby)
      ibx=ibx-mar+iax
      ibx=ideks(ibx+1)
      ibx=nit(ibx)
      ib=(mla-1)*mv +jb+iby
      idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if(jto.lt.iwod) go to 473
      jto=0
      jblk=jblk+1
      write(ntype) kc,kd
 473  kb=(mlr-1)*mq
      lb=kb+ii
      kb=kb+kk
      if(kb.lt.lb) go to 474
      ib=ideks(kb)+lb+ibx
      go to 475
 474  ib=ideks(lb)+kb+ibx
 475  idud=(ib-1)/igmax
      ib=ib-idud*igmax
      jto=jto+1
      kc(jto)=idud+1
      kd(jto)=ib
      if (jto.lt.iwod) go to 384
      jto=0
      jblk=jblk+1
      write (ntype) kc,kd
 384  continue
      if (kix.lt.0) go to 260
 383  if (nd.eq.0) go to 260
      kix=-1
      jz=nz
      ir=nd
      go to 472
 260  l3=l3+mx
  69  l4=l4+mx
      return
      end
      subroutine rumpb(iot,niott,ideks,ndeks,jcon,nd8)
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndeks,jcon,nd8
      dimension iot(niott),ideks(ndeks),jcon(nd8)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      integer icon
      common /cicon/ icon(ndk5)
      integer nod
      common /cnod/ nod(ndk5)
      integer nytl,ndub
      common /cny/ nytl(ndk5),ndub(ndk5)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
c
      integer lj
      parameter (nmmo=256)
      common /clj2/ lj(8)
c
      common /einf/ mhe(maxref),imain(maxref),iselct
c
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,ltape,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     +nn(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,j
c
      integer ican
      dimension  ican(nn3)
      integer mkon
      dimension  mkon(ndimh)
      character*10 charwall
c
*     integer i11(ndk5),i22(ndk5),i33(ndk5),i44(ndk5)
*     integer i110(ndk5),i220(ndk5)
c
      call rewftn(mtape)
      read (mtape)
      ih=imo*(lulu+1)
c
      do 2360 i=1,ih
        jcon(i) = nzero
2360  continue
      ny = nzero
      nnx = nzero
c
c remove sk from main-set
      do 501 i=1,jswh
       nc=mconf(i)
       if (nc.eq.0) go to 501
       nr=nod(i)
       nd=ndub(i)
c --- the mains are recoded in jcon :
c ---     unoccupied     orbital  : 0
c ---     singly occupied besetzt : 1
c ---     doubly occupied besetzt : 2
       do 502 j=1,nc
        nnx = nnx + imo
        if (nr.eq.0) go to 504
        do 2386 k=1,nr
          ny = ny + 1
          jd = iot(ny) + nnx
          jcon(jd) = 1
2386   continue
c
        if (nd.eq.0) go to 502
 504    continue
        do 2394 k=1,nd
         ny = ny + 1
         jd = iot(ny) + nnx
         jcon(jd) = 2
2394    continue
 502   continue
 501  continue
c381  format(2x,25i5)
      ky=imo
      do 2402 i=1,jswh
        ican(i)=ky
        ky=ky+imo*mconf(i)
2402  continue
      nt1=ny+1
      lulb=1
      ksc=nsc(1)
      jnerd=nn(1)
c - dumb reordering arrays
*     i110(1)= 2
*     i110(2)= 3
*     i110(3)= 1
*     i110(4)= 2
*     i110(5)= 3
*     i220(1)= 3
*     i220(2)= 1
*     i220(3)= 2
*     i220(4)= 3
*     i220(5)= 1
c
*     i11(1) = 2
*     i11(2) = 3
*     i11(3) = 4
*     i11(4) = 5
*     i11(5) = 1
c
*     i22(1) = 3
*     i22(2) = 4
*     i22(3) = 5
*     i22(4) = 1
*     i22(5) = 2
c
*     i33(1) = 5
*     i33(2) = 1
*     i33(3) = 2
*     i33(4) = 3
*     i33(5) = 4
c
*     i44(1) = 4
*     i44(2) = 5
*     i44(3) = 1
*     i44(4) = 2
*     i44(5) = 3
c
      do 506 i=1,iswh
       nc = nconf(i)
       if (nc.eq.0) go to 506
       read (mtape) mkon
       ig = nzero
       ndx = ndub(i)
       nrx = nod(i)
       ntx = nytl(i)
       iht = icon(i)
       ihz = iht + nrx
       nzx = ny + nrx
       if ((i+2).le.nn3) jy  = ican(i+2)
       if ((i+1).le.nn3) jx  = ican(i+1)
       nd1 = ndx + 1
       nd2 = ndx + 2
       md1 = ndx - 1
       md2 = ndx - 2
       nr3 = nrx - 2
       nr4 = nrx - 4
       mr1 = nrx + 2
       mr2 = nrx + 4
       nt3 = ntx - 1
       nt2 = ntx - 2
       mt1 = ntx + 1
       mt2 = ntx + 2
       if ((i+1).le.ndk5) ih1 = icon(i+1)
       if ((i+2).le.ndk5) ih2 = icon(i+2)
       if (i.gt.1) ih3 = icon(i-1)
       if (i.gt.2) ih4 = icon(i-2)
       jh2 = ih2 + mr2
       jh3 = ih3 + nr3
       jh4 = ih4 + nr4
       iab = ideks(nrx+2)
       nr1 = nrx+1
       if (nrx.gt.0) iac = ideks(nrx)
       nr2 = nrx - 1
       do 507 j=1,nc
        ky = ny
        do 508 k=1,ntx
         ig = ig + 1
         ky = ky + 1
         iot(ky)=mkon(ig)
         if (ig.lt.iwod) go to 508
         ig = 0
         read (mtape) mkon
 508   continue
       if (ksc.ne.i.or.jnerd.ne.j) go to 525
c if program arrives here test configuration is reference conf.
      mm = mm + 1
      olab(mm)=6
      if (mm.lt.nnid) go to 527
      mm=nzero
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
 527  if (lulb.eq.lulu) go to 528
      lulb=lulb+1
      ksc=nsc(lulb)
      jnerd=nn(lulb)
      go to 507
 528  ksc=0
      go to 507
c end of reference configuration part
 525  if (lsng.eq.0) go to 529
c change to include all single excitations with respect to all
c reference configurations
      do 10 ieng=1,lsng
c       write(6,*) mhe(ieng),lsng
        ianf=mhe(ieng)+1
        iend=ianf+imo
c       write(6,*)(jcon(ie1),ie1=ianf,iend)
        nix=0
        ky=ny
        if (nrx.eq.0) go to 530
        do 531 k=1,nrx
          ky=ky+1
          jd=iot(ky)+mhe(ieng)
          if (jcon(jd).gt.0) go to 531
c if nix.eq.1,more than single excitation,next reference conf.
          if (nix.eq.1) go to 10
          nix=1
 531    continue
        if (ndx.eq.0) go to 532
 530    do 533 k=1,ndx
          ky=ky+1
          jd=iot(ky)+mhe(ieng)
          if (jcon(jd)-1) 10,534,533
c if nix.eq.1,more than single excitation,next reference conf.
 534      if (nix.eq.1) go to 10
          nix=1
 533    continue
c if arrives here test conf is single excitation
 532    mm=mm+1
        olab(mm)=7
        if (mm.lt.nnid) go to 507
        mm=0
        call pack(olab8,8,olab,nnid)
        write (ntape) olab8
c branch to 507, next test conf. no further inf. required because
c konf. is single exc.
        go to 507
  10  continue
c end of changes
 529  ky=ny
      if (nrx.eq.0) go to 509
      do 510 k=1,nrx
      ky=ky+1
      jd=iot(ky)
 510  jcon(jd)=1
      if (ndx.eq.0) go to 511
 509  do 512 k=1,ndx
      ky=ky+1
      jd=iot(ky)
 512  jcon(jd)=2
 511  do 513 k=1,jswh
      mc=mconf(k)
      if (mc.eq.0) go to 513
      kfc=i-k+3
      go to (520,521,522,523,524),kfc
      go to 513
 522  call bearb(iot,niott,ideks,ndeks,jcon,nd8)
      go to 513
 520  call rumam(iot,niott,ideks,ndeks,jcon,nd8)
      go to 513
 524  call rumap(iot,niott,ideks,ndeks,jcon,nd8)
      go to 513
 521  call rumbm(iot,niott,ideks,ndeks,jcon,nd8)
      go to 513
 523  call rumbp(iot,niott,ideks,ndeks,jcon,nd8)
 513  continue
      ky=ny
      do 2589 l=1,ntx
        ky  = ky + 1
        ll  = iot(ky)
        jcon(ll) = nzero
2589  continue
 507  continue
 506  continue
      write(iwr,431) cpulft(1) ,charwall()
 431  format(/1x,'end of label generation at ', f10.2,
     +           ' seconds',a10,' wall'/)
      return
      end
****************************************************
      subroutine rumplz(iot,niott,ideks,ndeks)
c
c   generates the matrix elements of the geyser-space
c
c  - the mains are read from iot-field
c
      implicit REAL (a-h,o-z)
c
      integer iot,niott,ideks,ndek
      dimension iot(niott),ideks(ndeks)
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
      common /a/ t(nlca)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      common /cnit/ nit(nitmax)
c
      parameter (nmmo=256)
      integer loc, nir, jcon
      common /cloc/ loc(nmmo), nir(nmmo), jcon(nmmo)
      integer nod
      common /cnod/ nod(ndk5)
      integer nytl,ndub
      common /cny/  nytl(ndk5),ndub(ndk5)
      integer icon
      common /cicon/ icon(ndk5)
c
c
      integer lj
      common /clj2/ lj(8)
c
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer kc,kd,lab
      common /ckd/ kc(ndimh),kd(ndimh),lab(3)
c
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,ltype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +mconf(ndk5),nconf(ndk5),mh,jswh,
     +id(maxref), ie(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
c
*     write(6,*)'rumplz : ntape = ',ntape
*     write(6,*)'rumplz : ntype = ',ntype
*     write(6,*) 'rumplz: id,ie'
*     do i=1,maxref
*        write(6,*) id(i),ie(i)
*     end do
c
c --- rewind on ntape ????
*     write(6,*) 'rewind ntape'
      call rewftn(ntape)
c -------------------------------------
c --- zero out arrays
      do 2656 ii=1,ndimh
        kc(ii) = nzero
        kd(ii) = nzero
2656  continue
      mm   = nzero
      jto  = nzero
      jblk = nzero
cdebugwrite(6,*) 'imo=',imo
cdebug       do 2664 ii=1,5
cdebug          write(6,*)'nod(',ii,')=',nod(ii)
cdebug 2664  continue
      do 2667 jj=1,imo
         jcon(jj) = nzero
2667  continue
c -------------------------------------
c
cdebugwrite(6,*) 'iswh=',iswh
      do 3870 i=1,iswh
*      write(6,*) 'i=',i
       nc = nconf(i)
*      write(6,*) 'nc=',nc
       if (nc.eq.0) go to 11
       lt  = icon(i)
       nnx  = nytl(i)
       nd  = ndub(i)
       nr  = nod(i)
c      write(6,*) 'nr=',nr
       if (nr.lt.1) then
c         write(6,*) 'warning !!!!'
c         irite(6,*) 'iab set to 1 '
          iab = 1
       else
          iab = ideks(nr)
       end if
c      write(6,*) 'iab=',iab
       nr1 = nr - 1
       nr2 = nr - 2
       lz  = lt + nr
*      write(6,*) 'i=',i,'lz=',lz
       if (i.lt.3) go to 12
c&     call cpu
c&     call setcpu
c
c --- dk=2
c
*      write(6,*) ' ---- dk=2 ----'
       j8 = i - 2
       mc = nconf(j8)
       if (mc.eq.0) go to 100
        mx = nnx - 2
        md = nd + 2
        mr = nr - 4
        nt = lt
        nz = lz
        jt = icon(j8)
        iz = mr + jt
        do 13 j=1,nc
         it = nt
         mt = jt
         mz = iz
         do 2718 k=1,nr
          it = it + 1
          ll = iot(it)
          jcon(ll) = 1
2718     continue
        if (nd.ne.nzero) then
         do 2724 k=1,nd
          it = it + 1
          ll = iot(it)
          jcon(ll) = 2
2724     continue
        end if
c  16   continue
        do 17  k=1,mc
         jz = mz
         do 19 l=1,md
           jz = jz + 1
           ll = iot(jz)
           if (jcon(ll)-1) 20,21,19
   21      ia = l + 1
         go to 22
   19   continue
   20   mm = mm + 1
        olab(mm) = nzero
        if (mm.lt.nnid) go to 61
        call pack(olab8,8,olab,nnid)
        write (ntape) olab8
        mm=0
        go to 61
   22   do 23 l=ia,md
          jz=jz+1
          nn=iot(jz)
          if (jcon(nn)-1) 20,24,23
   24     ja=l+1
          go to 25
   23   continue
   25 if (ja.gt.md) go to 26
      do 27 l=ja,md
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii).ne.2) go to 20
   27 continue
   26 if (mr.eq.0) go to 28
      jz=mt
      do 2762 l=1,mr
        jz = jz + 1
        ii = iot(jz)
        if (jcon(ii).ne.1) go to 20
2762  continue
   28 jz = nt
      kz = mt + 1
      jq = nzero
      if (mr.gt.0) iq=iot(kz)
      do 2777 l=1,nr
       jz = jz + 1
       ii = iot(jz)
       if (ii.eq.ll) go to 31
       if (jq.eq.mr.or.iq.ne.ii) go to 32
       jq =jq + 1
       kz =kz + 1
       if (jq.lt.mr) then
          iq = iot(kz)
       end if
2777  continue
  31  ip = l
      ia = l + 1
      do 2788 l=ia,nr
       jz = jz + 1
       ii = iot(jz)
       if (ii.eq.nn) go to 34
       if (jq.eq.mr.or.iq.ne.ii) go to 35
       jq = jq + 1
       kz = kz + 1
       if (jq.lt.mr) iq=iot(kz)
2788  continue
  32  ip=l
      ia=l+1
      do 2799 l=ia,nr
       jz = jz + 1
       jj = iot(jz)
       if (jj.eq.ll) go to 37
       if (jq.eq.mr.or.iq.ne.jj) go to 38
       jq = jq + 1
       kz = kz + 1
       if (jq.lt.mr) iq=iot(kz)
2799  continue
  34  jp=ideks(l-1)
      ia=l+1
      do 39 l=ia,nr
      jz=jz+1
      ii=iot(jz)
      if (jq.eq.mr.or.iq.ne.ii) go to 40
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
  39  continue
  40  kp = jdeks(l-2)
      ia = l+1
      do 2819 l=ia,nr
       jz = jz + 1
       jj = iot(jz)
       if (jq.eq.mr.or.iq.ne.jj) go to 42
       jq = jq + 1
       kz = kz + 1
       if (jq.lt.mr) iq=iot(kz)
2819  continue
c
  42  mm = mm + 1
      olab(mm) = 1
      lp = kdeks(l-3)
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm = nzero
      go to 43
  35  jp=ideks(l-1)
      ia = l + 1
      do 2839 l=ia,nr
       jz = jz + 1
       jj = iot(jz)
       if (jj.eq.nn) go to 45
       if (jq.eq.mr.or.jj.ne.iq) go to 46
       jq = jq + 1
       kz = kz + 1
       if (jq.lt.mr) iq=iot(kz)
2839  continue
  37  jp=ideks(l-1)
      ia=l+1
      do 47 l=ia,nr
       jz=jz+1
       jj=iot(jz)
       if (jj.eq.nn) go to 48
       if (jq.eq.mr.or.jj.ne.iq) go to 49
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
   47 continue
   38 jp=ideks(l-1)
      ia=l+1
      do 50 l=ia,nr
       jz=jz+1
       kk=iot(jz)
       if (kk.eq.ll) go to 51
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
   50 continue
   51 kp=jdeks(l-2)
      ia=l+1
      do 52 l=ia,nr
       jz=jz+1
       kk=iot(jz)
       if (kk.eq.nn) go to 42
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
   52 continue
   48 kp = jdeks(l-2)
      ia = l + 1
      do 53 l=ia,nr
      jz = jz + 1
      jj = iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 54
      jq = jq + 1
      kz = kz + 1
      if (jq.lt.mr) iq=iot(kz)
   53 continue
   54 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   49 kp=jdeks(l-2)
      ia=l+1
      do 55 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   55 continue
   56 lp=kdeks(l-3)
      mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 43
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
      go to 43
   45 kp=jdeks(l-2)
      ia=l+1
      do 57 l=ia,nr
      jz=jz+1
      jj=iot(jz)
      if (jq.eq.mr.or.jj.ne.iq) go to 56
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   57 continue
   46 kp=jdeks(l-2)
      ia=l+1
      do 58 l=ia,nr
      jz=jz+1
      kk=iot(jz)
      if (kk.eq.nn) go to 54
      jq=jq+1
      kz=kz+1
      if (jq.lt.mr) iq=iot(kz)
   58 continue
 43   mm=mm+1
      olab(mm)=ip+jp+kp+lp
      if (mm.lt.nnid) go to 170
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 170  nt1r = nir(ll)
      nl1r = loc (ll)
      nt2r = nir (nn)
      nl2r = loc (nn)
      kix=3
      nb1r = nir(ii)
      nm1r = loc (ii)
      nb2r = nir (jj)
      nm2r = loc (jj)
  120 if (nt1r-nb1r)123,122,199
  199 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 124
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 125,126,127
  127 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  140 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 129
  126 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  141 idq = (nl2r-1)*ljn+nm2r
      go to 133
  125 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  135 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 129
  124 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 130,131,132
  132 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 129
  131 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  133 if (icq.lt.idq) go to 134
      icq = ideks(icq)+idq
      go to 129
  134 icq = ideks(idq)+icq
      go to 129
  130 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 135
  123 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 136
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 137,138,139
  139 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 140
  138 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 141
  137 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  145 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 129
  136 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 142,143,144
  144 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 129
  143 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 133
  142 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 145
  122 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 146
      iay=ideks(nl1r)+nm1r
      go to 147
  146 iay= ideks(nm1r)+nl1r
  147 if (nl2r.lt.nm2r) go to 148
      iby=ideks(nl2r)+nm2r
      go to 149
  148 iby=ideks(nm2r)+nl2r
  149 if (nt1r.eq.nt2r) go to 150
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 151
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 129
  151 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 129
  150 icx=ideks(iax+1)
      if (iay.lt.iby) go to 153
      icq=ideks(iay)+iby
      go to 129
  153 icq=ideks(iby)+iay
 129  icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 160
      write(ntype)  kc,kd
clear    write(77,*)  kc,kd
      jto=0
      jblk=jblk  + 1
  160 if (kix.lt.0) go to 61
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
      go to 120
   61 mz=mz+mx
   17 mt=mt+mx
      do 62 k=1,nnx
      nt=nt+1
      ll=iot(nt)
   62 jcon(ll) = nzero
   13 nz=nt+nr
      go to 100
c
c --- dk=1
c
  12  continue
c     write(6,*) ' ---- dk=1 ---- '
      if (i.eq.1) go to 64
  100 j8 = i - 1
      mc = nconf(j8)
c&    call cpu
c&    call setcpu
      if (mc.eq.0) go to 64
      mx = nnx - 1
      md = nd + 1
      mr = nr - 2
      jt = icon(j8)
      iz = mr + jt
      nt = lt
c ---- big loop over configurations of k
      do 3862 j=1,nc
       it = nt
       nt1= nt + 1
       nz = nt +nr
       mt = jt
       mz = iz
c building jcon : 0,1,2
       do 66 k=1,nr
        it=it+1
        ll=iot(it)
   66  jcon(ll)=1
       if (nd.eq.0) go to 67
       do 68 k=1,nd
        it = it + 1
        ll = iot(it)
   68  jcon(ll)=2
c do loop over configurations of k-1
   67  do 69 k=1,mc
c check on place level
cdebug  write(6,*) 'j,k ',j,k
        jz = mz
        do 70 l=1,md
         jz=jz+1
         jj=iot(jz)
         if (jcon(jj)-1) 71,72,70
   70   continue
   71   if (l.eq.md) go to 73
        ia = l + 1
        do 74 l=ia,md
         jz=jz+1
         ll=iot(jz)
         if (jcon(ll).ne.2) go to 75
   74   continue
        go to 73
  75    mm = mm + 1
c no interaction 
        olab(mm) = 1
        if (mm.lt.nnid) go to 260
        call pack(olab8,8,olab,nnid)
        write (ntape) olab8
        mm = nzero
        go to 260
   72   l1 = l
        if (l.eq.md) go to 76
        ia = l + 1
        do 78 l=ia,md
         jz = jz + 1
         ll = iot(jz)
         if (jcon(ll)-1) 75,77,78
  78    continue
        go to 76
  77    l2 = l
        if (l.eq.md) go to 79
        ia = l + 1
        do 80 l=ia,md
         jz = jz + 1
         kk = iot(jz)
         if (jcon(kk).ne.2) go to 75
   80   continue
        go to 79
   73   if (mr.eq.0) go to 81
        jz = mt
        do 82 l=1,mr
          jz = jz + 1
          ll = iot(jz)
          if (jcon(ll).ne.1) go to 75
   82   continue
  81   mm=mm+1
c p=2
       olab(mm)=3
       if (mm.lt.nnid) go to 83
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       mm = nzero
  83   jz = nt
       kz = mt + 1
       jq = nzero
       if (mr.gt.0) iq=iot(kz)
      do 84 l=1,nr
       jz=jz+1
       ll=iot(jz)
       if (jq.eq.mr.or.iq.ne.ll) go to 85
       jq = jq + 1
       kz = kz + 1
       if (jq.lt.mr) iq=iot(kz)
   84 continue
   85 ia = l + 1
      ip = l
      do 86 l=ia,nr
       jz=jz+1
       nn=iot(jz)
       if (jq.eq.mr.or.iq.ne.nn) go to 87
       jq=jq+1
       kz=kz+1
       if (jq.lt.mr) iq=iot(kz)
  86  continue
  87  mm = mm + 1
      olab(mm)=ideks(l-1)+ip
      if (mm.lt.nnid) go to 88
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
c computing the integral addresses
  88  nt1r = nir(jj)
      nl1r = loc(jj)
      nb1r = nir(ll)
      nm1r = loc(ll)
      nm2r = loc(nn)
  220 if (nt1r-nb1r) 223,222,299
  299 iax = ideks(nt1r)+nb1r
      icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn
      idq = icq+nm2r
      icq=ideks(idq)+icq+nm1r
      go to 229
  223 iax = ideks(nb1r)+nt1r
      icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl1r
      icq=ideks(idq)+icq
      go to 229
  222 iax=ideks(nt1r+1)
      icx=ideks(iax+1)
      if (nl1r.lt.nm1r) go to 246
      iay=ideks(nl1r)+nm1r
      go to 247
  246 iay= ideks(nm1r)+nl1r
      iby=ideks(nm2r)+nl1r
      go to 253
  247 if (nl1r.lt.nm2r) go to 248
      iby=ideks(nl1r)+nm2r
      go to 253
  248 iby=ideks(nm2r)+nl1r
  253 icq=ideks(iby)+iay
  229 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
cdebugwrite(6,*) 'jto,kc,kd',jto,kc(jto),kd(jto)
      if(jto.lt.iwod)  go to 260
      write(ntype) kc,kd
      jto=0
      jblk=jblk  + 1
      go to 260
   79 if (mr.eq.0) go to 89
      jz=mt
      do 90 l=1,mr
      jz=jz+1
      ii=iot(jz)
      if (jcon(ii)-1) 75,90,91
  90  continue
      go to 89
  91  l1=l
      if (l.eq.mr) go to 93
      ia=l+1
      do 92 l=ia,mr
      jz=jz+1
      jc=iot(jz)
      if (jcon(jc).ne.1) go to 75
  92  continue
      go to 93
  89  mm=mm+1
      olab(mm)=3
      if (mm.lt.nnid) go to 94
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   94 jz=nt
      do 95 l=1,nr
      jz=jz+1
      if (iot(jz).eq.jj) go to 96
   95 continue
   96 ia=l+1
      ip=l
      do 97 l=ia,nr
      jz=jz+1
      if (iot(jz).eq.ll) go to 98
   97 continue
  98  mm=mm+1
      olab(mm)=-ip-ideks(l-1)
      if (mm.lt.nnid) go to 99
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
   99 jz=nz
      if (l1.eq.1) go to 180
      kz=mz
      ja=l1-1
      do 181 l=1,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  181 continue
  180 if (l2.eq.l1+1) go to 183
      kz=mz+l1
      ia=l1+1
      ja=l2-1
      do 184 l=ia,ja
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  184 continue
  183 if (l2.lt.md) go to 185
      kk=iot(jz+1)
      go to 182
  185 kz=mz+l2
      ia=l2+1
      do 186 l=ia,md
      kz=kz+1
      jz=jz+1
      kk=iot(jz)
      if (kk.ne.iot(kz)) go to 182
  186 continue
      kk=iot(jz+1)
  182 nt1r=nir(kk)
      nl1r=loc(kk)
      nb1r=nir(jj)
      nm1r=loc(jj)
      nm2r=loc(ll)
      go to 220
  93  mm=mm+1
      olab(mm)=2
      if (mm.lt.nnid) go to 187
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      mm=0
 187  jz=nt
      kz=mt
      ig=0
      if (l1.eq.1) go to 188
      ja=l1-1
      ka=1
      do 189 l=ka,ja
      kz=kz+1
      jp=iot(kz)
 802  jz=jz+1
      if (jp.eq.iot(jz)) go to 189
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 802
  189 continue
  188 if (l1.eq.mr) go to 191
      ia=l1+1
      kz=kz+1
      do 192 l=ia,mr
      kz=kz+1
      jp=iot(kz)
 803  jz=jz+1
      if (jp.eq.iot(jz)) go to 192
      ig=ig+1
      lab(ig)=jz-nt
      if (ig.eq.3) go to 190
      go to 803
  192 continue
  191 lab(3)=nr
      if (ig-1) 284,285,190
  284 lab(2)=nr1
      lab(1)=nr2
      go to 190
  285 lab(2)=nr1
  190 jz=lab(1)+nt
      nn=iot(jz)
      mm=mm+1
      if (nn.eq.jj) go to 194
      olab(mm)=1
  195 if (mm.lt.nnid) go to 196
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  196 mm=mm+1
      olab(mm)=l1
      if (mm.lt.nnid) go to 197
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
  197 mm=mm+1
      jp=lab(2)-1
      kp=lab(3)-2
      olab(mm)=lab(1)+ideks(jp)+jdeks(kp)
      if (mm.lt.nnid) go to 198
      mm=0
      call pack(olab8,8,olab,nnid)
      write (ntape) olab8
      go to 198
  194 jz=lab(2)+nt
      nn=iot(jz)
      if (nn.eq.ll) go to 271
      olab(mm)=2
      go to 195
  271 jz=lab(3)+nt
      nn=iot(jz)
      olab(mm)=3
      go to 195
  198 nt1r=nir(jj)
      nt2r=nir(ll)
      nl1r=loc(jj)
      nl2r=loc(ll)
      kix=3
      nb1r=nir(nn)
      nb2r=nir(ii)
      nm1r=loc(nn)
      nm2r=loc(ii)
  320 if (nt1r-nb1r) 323,322,399
  399 iax = ideks(nt1r)+nb1r
      if (nt2r.lt.nb2r) go to 324
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 325,326,327
  327 icx = ideks(iax) + ibx
      icq = (nl1r-1)*lj(nb1r) + nm1r
  340 ljn = lj(nb2r)
      idq = (nl2r-1)*ljn + nm2r
      icq = (icq-1)*ljn*lj(nt2r)+idq
      go to 329
  326 icx = ideks (iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
  341 idq = (nl2r-1)*ljn+nm2r
      go to 333
  325 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r)+nm2r
  335 ljn = lj(nb1r)
      idq = (nl1r-1)*ljn+nm1r
      icq =(icq-1)*ljn*lj(nt1r) +idq
      go to 329
  324 ibx = ideks(nb2r) + nt2r
      if (iax-ibx) 330,331,332
  332 icx = ideks(iax)+ibx
      icq = (nl1r-1)*lj(nb1r)+nm1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r) + idq
      go to 329
  331 icx = ideks(iax+1)
      ljn = lj(nb1r)
      icq = (nl1r-1)*ljn+nm1r
      idq = (nm2r-1)*ljn+nl2r
  333 if (icq.lt.idq) go to 334
      icq = ideks(icq)+idq
      go to 329
  334 icq = ideks(idq)+icq
      go to 329
  330 icx=ideks(ibx)+iax
      icq = (nm2r-1)*lj(nt2r)+nl2r
      go to 335
  323 iax = ideks(nb1r)+nt1r
      if (nt2r.lt.nb2r) go to 336
      ibx = ideks(nt2r)+nb2r
      if (iax-ibx) 337,338,339
  339 icx = ideks(iax)+ibx
      icq = (nm1r-1)*lj(nt1r)+nl1r
      go to 340
  338 icx = ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      go to 341
  337 icx = ideks(ibx) + iax
      icq = (nl2r-1)*lj(nb2r) +nm2r
  345 ljn = lj(nt1r)
      idq = (nm1r-1)*ljn         +nl1r
      icq = (icq-1)*ljn*lj(nb1r) +idq
      go to 329
  336 ibx = ideks(nb2r)+nt2r
      if (iax-ibx) 342,343,344
  344 icx = ideks(iax)+ibx
      icq =(nm1r-1) *lj(nt1r)+nl1r
      ljn = lj(nt2r)
      idq = (nm2r-1)*ljn+nl2r
      icq = (icq-1)*ljn*lj(nb2r)+idq
      go to 329
  343 icx=ideks(iax+1)
      ljn = lj(nt1r)
      icq = (nm1r-1)*ljn+nl1r
      idq = (nm2r-1)*ljn+nl2r
      go to 333
  342 icx = ideks(ibx)+iax
      icq =(nm2r-1)*lj(nt2r)+nl2r
      go to 345
  322 iax=ideks(nt1r+1)
      if (nl1r.lt.nm1r) go to 346
      iay=ideks(nl1r)+nm1r
      go to 347
  346 iay= ideks(nm1r)+nl1r
  347 if (nl2r.lt.nm2r) go to 348
      iby=ideks(nl2r)+nm2r
      go to 349
  348 iby=ideks(nm2r)+nl2r
  349 if (nt1r.eq.nt2r) go to 350
      ibx=ideks(nt2r+1)
      if (iax.lt.ibx) go to 351
      icx=ideks(iax)+ibx
      ljn =lj(nt2r)
      icq = (iay-1)*ideks(ljn+1)+iby
      go to 329
  351 icx=ideks(ibx)+iax
      ljn=lj(nt1r)
      icq =(iby-1)*ideks(ljn+1)+iay
      go to 329
  350 icx=ideks(iax+1)
      if (iay.lt.iby) go to 353
      icq=ideks(iay)+iby
      go to 329
  353 icq=ideks(iby)+iay
  329 icq=nit(icx)  +  icq
      idud=(icq -1)/igmax
      jto=jto + 1
      kc(jto) = idud  + 1
      icq=icq  -  idud*igmax
      kd(jto)  = icq
      if(jto.lt.iwod)  go to 360
      write (ntype) kc,kd
      jto = nzero
      jblk=jblk  + 1
 360  if (kix.lt.0) go to 260
      kix=-3
      itr=nb1r
      nb1r=nb2r
      nb2r=itr
      itr=nm1r
      nm1r=nm2r
      nm2r=itr
       go to 320
  76   if (mr.eq.0) go to 273
       jz = mt
       do 274 l=1,mr
        jz = jz + 1
        ll = iot(jz)
        if (jcon(ll)-1) 275,274,75
 274   continue
       go to 273
 275   l1=l
       if (l.eq.mr) go to 277
       ia=l+1
       do 276 l=ia,mr
        jz = jz + 1
        nn = iot(jz)
        if (jcon(nn).ne.1) go to 75
 276   continue
 277   mm=mm+1
       olab(mm)=2
       if (mm.lt.nnid) go to 278
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       mm = nzero
 278   jz = nt
       kz = mt
       ig = nzero
       if (l1.eq.1) go to 279
       ja=l1-1
       do 280 l=1,ja
        kz=kz+1
        jp=iot(kz)
 254    jz=jz+1
        if (jp.eq.iot(jz)) go to 280
        ig=ig+1
        lab(ig)=jz-nt
        if (ig.eq.3) go to 281
        go to 254
 280   continue
 279   if (l1.eq.mr) go to 282
       ia=l1+1
       kz=kz+1
       do 283 l=ia,mr
       kz=kz+1
       jp=iot(kz)
 255   jz=jz+1
       if (jp.eq.iot(jz)) go to 283
       ig=ig+1
       lab(ig)=jz-nt
       if (ig.eq.3) go to 281
       go to 255
  283  continue
  282  lab(3)=nr
       if (ig-1) 286,287,281
  286  lab(2)=nr1
       lab(1)=nr2
       go to 281
  287  lab(2)=nr1
  281  jz=lab(1)+nt
       kk=iot(jz)
       mm=mm+1
       if (kk.ne.jj) go to 288
       olab(mm)=-1
       jz=lab(2)+nt
       kk=iot(jz)
       jz=lab(3)+nt
       nn=iot(jz)
 289   if(mm.lt.nnid) go to 290
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 290   mm=mm+1
       olab(mm)=l1
       if(mm.lt.nnid) go to 291
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 291   mm = mm + 1
       jp = lab(2) - 1
       kp = lab(3) - 2
       olab(mm)=lab(1)+ideks(jp) +jdeks(kp)
       if(mm.lt.nnid)go to 292
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       go to 292
 288   jz = lab(2)+nt
       nn = iot(jz)
       if (nn.ne.jj) go to 293
       olab(mm) = -2
       jz = lab(3)+nt
       nn = iot(jz)
       go to 289
 293   olab(mm) = -3
       go to 289
 292   nt1r = nir(jj)
       nt2r = nir(ll)
       nl1r = loc(jj)
       nl2r = loc(ll)
       kix  = 3
       nb1r = nir(kk)
       nb2r = nir(nn)
       nm1r = loc(kk)
       nm2r = loc(nn)
       go to 320
 273   mm=mm+1
       olab(mm)=4
       if(mm.lt.nnid) go  to 294
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 294   jz=nt
       mm=mm+1
       if(mr.eq.0) go to 297
       kz=mt
       do 295 l=1,mr
        kz=kz+1
        jz=jz+1
        ll=iot(jz)
        if(ll.ne.iot(kz)) go to 296
 295   continue
 297   ll=iot(jz+1)
       if(ll.eq.jj) go to 298
       olab(mm) =iab
       if(mm.lt.nnid) go to 371
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       go to 371
 298   ll=iot(jz+2)
       olab(mm) = -iab
       if(mm.lt.nnid) go to 371
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       go to 371
 296   l1 =jz-nt
       if(ll.eq.jj) go to 373
       do 374 l=1,nr
        jz=jz+1
        if(iot(jz).eq.jj) go to 375
 374   continue
 375   olab(mm)=l1+ideks(jz-nt1)
       if(mm.lt.nnid) go to 371
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
       go to 371
 373   kz=kz-1
       do 376 ia=l,mr
       kz=kz+1
       jz=jz+1
       ll=iot(jz)
       if(ll.ne.iot(kz)) go to 377
 376   continue
       jz=jz+1
       ll=iot(jz)
 377   olab(mm)=-l1-ideks(jz-nt1)
       if(mm.lt.nnid) go to 371
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 371   mm=mm+1
       iax=nir(jj)
       olab(mm)=iax
       if(mm.lt.nnid) go to 372
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
  372  iay=ideks(iax+1)
       iaz=ideks(iay+1)
       mq=lj(iax)
       mv=ideks(mq+1)
       iaq=iay-iax
       iaw=iaz-iay
       iat=nit(iaz)
       ii=loc(jj)
       kk=loc(ll)
       km=ideks(kk)
       nm=ideks(ii)
       nn=nm+ii
       if(ii.gt.kk) go to 380
       mm=mm+1
       jb=ii+km
       olab(mm)=kk
       if(mm.lt.nnid) go to 378
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 378   mm=mm+1
       olab(mm)=ii
       if(mm.lt.nnid) go to 379
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 379   ib = ideks(jb) +nn+iat
       go to 381
 380   lb=kk-ii
       jb = nn + lb
       ib = ideks(nn+1)+lb+iat
       mm = mm + 1
       olab(mm) = ii
       if(mm.lt.nnid) go to 900
       mm = nzero
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 900   mm=mm+1
       olab(mm)=kk
       if(mm.lt.nnid) go to 381
       mm=0
       call pack(olab8,8,olab,nnid)
       write (ntape) olab8
 381   idud=(ib-1)/igmax
       ib=ib-idud*igmax
       jto=jto+1
       kc(jto)=idud+1
       kd(jto)=ib
       if(jto.lt.iwod) go to 382
       jto=0
       write(ntype) kc,kd
       jblk=jblk+1
 382   if(mr.eq.0) go to 383
       jz=mt
       ir=mr
       kix=1
 472   do 384 l=1,ir
       jz=jz+1
       mi=iot(jz)
       mar=nir(mi)
       mlr=loc(mi)
       mlb=ideks(mlr)
       mla=mlb+mlr
       if (mar-iax) 385,395,396
 395   if(mla.lt.jb) go to 386
       ib=iat+ideks(mla) +jb
 388   idud    = (ib-1)/igmax
       ib      = ib - idud*igmax
       jto     = jto + 1
       kc(jto) = idud + 1
       kd(jto) = ib
       if(jto.lt.iwod) go to 387
       jto = nzero
       write(ntype) kc,kd
       jblk = jblk + 1
       go to 387
 386   ib = iat + ideks(jb) +mla
       go to 388
 387   if(mlr.lt.ii) go to 389
       kb = mlb + ii
       go to 390
 389   kb = nm+mlr
 390   if(mlr.lt.kk) go to 391
       lb = mlb+kk
       go to 392
 391   lb = mlr+km
 392   if(kb.lt.lb) go to 393
       ib=iat+ideks(kb) +lb
       go to 394
 393   ib=iat+ideks(lb) +kb
 394   idud=(ib-1)/igmax
       ib=ib-idud*igmax
       jto=jto+1
       kc(jto)=idud+1
       kd(jto)=ib
       if(jto.lt.iwod) go to 384
       jto = nzero
       jblk=jblk+1
       write(ntype) kc,kd
       go to 384
 385   iby = ideks(mar+1)
       iby = iaw + iby
       ibx = iaq + mar
       ibx = ideks(ibx+1)
       ibx = nit(ibx)
       iby = nit(iby)
       ml  = lj(mar)
       ib=ideks(ml+1)*(jb-1) +mla+iby
       idud=(ib-1)/igmax
       ib=ib-idud*igmax
       jto=jto+1
       kc(jto)=idud+1
       kd(jto)=ib
       if(jto.lt.iwod) go to 397
       jto=0
       jblk=jblk+1
       write(ntype) kc,kd
 397   kb=(ii-1)*ml +mlr
       lb=(kk-1)*ml +mlr
       if(kb.lt.lb) go to 398
       ib=ideks(kb)+lb+ibx
       go to 471
 398   ib=ideks(lb) +kb+ibx
 471   idud=(ib-1)/igmax
       ib=ib-idud*igmax
       jto=jto+1
       kc(jto)=idud+1
       kd(jto) = ib
       if(jto.lt.iwod) go to 384
       jto = nzero
       jblk=jblk+1
       write(ntype) kc,kd
       go to 384
396    ibx=ideks(mar+1)
       iby=ideks(ibx) +iay
       iby=nit(iby)
       ibx=ibx-mar+iax
       ibx=ideks(ibx+1)
       ibx=nit(ibx)
       ib=(mla-1)*mv +jb+iby
       idud=(ib-1)/igmax
       ib=ib-idud*igmax
       jto = jto + 1
       kc(jto)=idud+1
       kd(jto)=ib
       if(jto.lt.iwod) go to 473
       jto = nzero
       jblk=jblk+1
       write(ntype) kc,kd
 473   kb=(mlr-1)*mq
       lb = kb + ii
       kb = kb + kk
       if(kb.lt.lb) go to 474
       ib =ideks(kb) + lb + ibx
       go to 475
 474   ib = ideks(lb) + kb + ibx
 475   idud = (ib-1)/igmax
       ib  = ib - idud*igmax
       jto = jto + 1
       kc(jto) = idud + 1
       kd(jto) = ib
       if (jto.lt.iwod) go to 384
       jto = nzero
       jblk=jblk+1
       write (ntype) kc,kd
 384   continue
       if (kix.lt.0) go to 260
 383   if (nd.eq.0) go to 260
       kix = -1
       jz  = nz
       ir  = nd
       go to 472
 260   mz=mz + mx
  69   mt = mt + mx
       do 3861 k=1,nnx
         nt = nt + 1
         ll = iot(nt)
         jcon(ll) = nzero
3861   continue
3862  continue
c---- end of big loop
c
   64 continue
c     write(6,*) 'rumplz calls bearz'
      call bearz(iot,niott,ideks,ndeks)
c
  11  continue
3870  continue
      return
      end
****************************************************
      subroutine shuffl(a,nt64,e,ndimf,sx,lsx,
     +                  vect,nvect,jkan,njkan)
c
      implicit REAL (a-h,o-z)
c
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
c
      REAL a,e,sx,vect
      integer ndimf,lsx,jkan,njkan,nvect
      dimension a(nt64),e(ndimf),sx(lsx),jkan(njkan),vect(nvect)
c
      REAL crus,sch
      parameter (crus = 1000.0d0)
c
      common /g/ trwe,trwf,trwg,trwh,trwi,trwj,imsec
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,nconf(5),mconf(5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
c
       common /d/ jbab,kml,ifl,ii,nis,iiz,nzer,kmj,ndj,jsac,ie,iej,
     * ja,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32,jcb,j3b
c
      common /a/ t(nlca)
      integer nytl,ndub
      common /cny/ nytl(ndk5),ndub(ndk5)
      integer ihog
      common /cihog/ ihog(48)
      integer imap
      REAL trsum,ew
      integer itym
      common /cshu1/ ew(mxroot),trsum(10,mxroot),
     +               istm(nopmax+1),itym(mxroot),
     +               imap(504)
c
      integer ij, nplu
      common /scrtch/ij(maxsym),nplu(ndk5)
c
      integer ikan, iqr, nrec
      REAL xemp
      common /miscop/ xemp(500), iqr(ndimh), ikan(ndimh), nrec(ndk5)
c
      REAL trash,pmt,trwe,trwf,trwg
      REAL trwh,trwi,trwj,trwa,trwb,trwc,trwd
      REAL tdel
c
cccc  crus=1000.0d0
c --- kfile
      call rewftn(kfile)
      read (kfile)
      call rewftn(ideli)
      ltype=mtype
      call rewftn(ltype)
      ltape=linf
      call rewftn(ltape)
      read (ltape)
c
      xemp(5)=trwe
      xemp(6)=trwf
      xemp(7)=trwg
      xemp(8)=trwh
      xemp(9)=trwi
      xemp(10)=trwj
c ----
      i=4
  449 i=i+1
      if (i.gt.10) go to 450
      if (istm(i).gt.imsec) go to 449
      trash=xemp(i)
      sch=trash*crus
      isp=i
      write(iwr,451) sch
      read(ideli) trwa,trwb,trwc,trwd,tdel
      read(ideli) xemp
      read (kfile) iqr
      kft=0
      imk=0
      imn=0
      do 452 i=1,iswh
      nc=nconf(i)
      if (nc.eq.0) go to 452
      nnx=nytl(i)
      read(ltape)ndt,kml,a
      write(ltype) ndt,kml,a
      read(ltape) ihog,imap,e
      write(ltype) ihog,imap,e
      read(ltape) jkan
      lpix=0
      jmk=0
      ibc=1
      do 453 j=1,nc
      kft=kft+1
      kx=iqr(kft)
      if (kx.gt.4) go to 708
      if (kft.lt.iwod) go to 453
      kft=0
      read (kfile) iqr
      go to 453
 708  imk=imk+kml
      if (imk.le.irsf) go to 3273
      imk=imk-irsf
      read(ideli) xemp
3273  pmt=xemp(imk)
      if (kft.lt.iwod) go to 707
      kft=0
      read (kfile) iqr
 707  if (pmt.lt.trash) go to 454
      do 455 k=1,kml
      imn=imn+1
455   sx(imn)=pmt
      do 456 k=1,nnx
      lpix=lpix+1
      jmk=jmk+1
      ikan(jmk)=jkan(lpix)
      if (lpix.lt.iwod) go to 457
      lpix=0
      read(ltape) jkan
457   if (jmk.lt.iwod) go to 456
      jmk=0
      ibc=ibc+1
      write(ltype) ikan
456   continue
      go to 453
454   lpix=lpix+nnx
      if (lpix.lt.iwod) go to 453
      lpix=lpix-iwod
      read(ltape) jkan
453   continue
      ikan(jmk+1)=0
      write(ltype) ikan
      nrec(i)=ibc
452   continue
      call rewftn(kfile)
      call rewftn(ltape)
      call rewftn(ltype)
      call rewftn(ideli)
      write(ideli)trwa,trwb,trwc,trash,tdel,isp,vect,ew
      i1=1-irsf
      i2=0
  464 i1=i1+irsf
      i2=i2+irsf
      write(ideli)(sx(i),i=i1,i2)
      if (i2.lt.imn) go to 464
      write (iwr,712) nrootx,imn,istm
      write(ideli)nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm
      read (ltape)
      do 460 i=1,iswh
       if (nconf(i).ne.nzero) then
        read (ltype) ndt,kml,a
        write(ltape) ndt,kml,a
        read(ltype) ihog,imap,e
        write(ltape) ihog,imap,e
        ibc=nrec(i)
        do 461 j=1,ibc
         read(ltype) ikan
         write(ltape) ikan
 461    continue
       end if
 460  continue
_IF(notused)
      call drclos
_ELSE
      call closbf3
_ENDIF
      return
 450   write(iwr,446) istm(4),imsec
 446   format(/10x,i6,'dimension exceeds maximum of',i7)
c
c     the following statements prevent rumple from getting into an
c     embarassing state ... it will thus just die with end of file
c
      call rewftn(ltape)
      write (ltape) istm(4),imsec
 451  format(/5x,'***** correct threshold is ',f11.5, '*****'/)
 712  format(5x,12i10)
      call caserr('selected no. of configurations invalid')
      return
      end
****************************************************
      subroutine singlep(kcon,nopcl,nelec,nset,ndeg,inocc,nseto,iret)
c-- created by hutter
c   proof with input given whether configuration should be hold or not
c--
      implicit REAL (a-h,o-z)
      parameter(mxset=10,mxdeg=100)
INCLUDE(common/iofile)
      dimension kcon(*)
      integer kcon
      logical    coreci
      dimension ndeg(3,*),inocc(3,*),nseto(3,mxdeg,*),nset(*)
c
       common /ftap/ ntape,mtape,mdisk,ideli,ltype,linf,
     *  ntab,kfile,kclab,ntype,mstvt, nf01, nf62,nhead,
     *  nf99, mtapev, nf11
      common /tap/ ical,m,nko,mxex,nmul,ispace,nprin,ncorci
c
      common /inkurt/ coreci
      integer nel,iret,nopcl,nelec,nset
      nclos = nelec-nopcl
      nopen = nopcl-nclos
c$    write(iwr,*)'nclos = ',nclos,'nopen = ',nopen
c$    write(iwr,*) 'single - your specified subroutine'
c-- loop over all required and specified orbitals
cc change: k.pfingst 8/10/92
      do 11 kl=1,3
      do 10 i=1,nset(kl)
       nel = 0
c$     write(iwr,*)'10a nel = ',nel
       do 70 j=1,ndeg(kl,i)
c$     write(iwr,*)'10b nel = ',nel
        iorb = nseto(kl,j,i)
c$      write(iwr,*)'nseto(',j,',',i,') = ',nseto(j,i)
c$      write(iwr,*)'iorb = nseto = ',iorb
c-- loop over open shells
        do 20 k=1,nopen
c$      write(iwr,*)'20 iorb = ',iorb
c$      write(iwr,*)'20 nel = ',nel
c$      write(iwr,*) 'kcon(',k,') = ',kcon(k)
        if (kcon(k).eq.iorb) then
c$       write(iwr,*)'if-case: nel = ',nel
         nel = nel+1
c$       write(iwr,*)'nel = ',nel
         goto 30
        endif
 20     continue
 30     continue
c-- loop over closed shells
        do 40 k=nopen+1,nopcl
c$      write(iwr,*) 'kcon(',k,') = ',kcon(k)
         if (kcon(k).eq.iorb) then
          nel = nel+2
c$        write(iwr,*)'nel = ',nel
          goto 50
         endif
 40     continue
 50     continue
 70    continue
       goto (1,2,3) kl
 1     if (nel.lt.inocc(kl,i)) then
c$     write(iwr,*)'inocc(',i,') = ',inocc(i)
        iret=1
c$      write(iwr,*)'iret = ',iret
        return
       endif
       goto 10
 2     if (nel.gt.inocc(kl,i)) then
c$     write(iwr,*)'inocc(',i,') = ',inocc(i)
        iret=1
c$      write(iwr,*)'iret = ',iret
        return
       endif
       goto 10
 3     if (nel.ne.inocc(kl,i)) then
c$     write(iwr,*)'inocc(',i,') = ',inocc(i)
        iret=1
c$      write(iwr,*)'iret = ',iret
        return
       endif
 10   continue
cc change: k.pfingst 8/10/92
 11   continue
      iret = 0
c$    write(iwr,*)'iret = ',iret
      return
c$    end
      end
****************************************************
      subroutine skin1a(pey,ideks,ndeks,
     +                  hp,nt64,f,ndimf,ot,nt444)
c
c  computes the matrix elements of the geyser-space for dk=+1
c
      implicit REAL (a-h,o-z)
      REAL sac
c
INCLUDE(common/newmrd_parinc)
c
      REAL pey
      integer ideks,ndeks
      dimension ideks(ndeks)
      dimension pey(ndeks)
c
      REAL hp, f
      integer nt64
      dimension hp(nt64),f(ndimf)
c
      integer ot
      dimension ot(nt444)
c
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      integer ij, nplu0
      common /scrtch/ ij(maxsym),nplu0(ndk5)
      common /csac/ sac(nopmax+1)
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
      common /a/ t(nlca)
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +md1,nston,mdisk,ideli,mtape,iswh,mr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,kc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,nd,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,nr,ipag
c
       common /d/ jbab,kmj,ifl,iiq,nis,iiz,kzer,kmk,ndt,jsec,iek,iej,
     * ia,mc,kml,ndk,kr,ii,ic,ic1,i2,i21,i3,i31,iix,
     * nzer,ndj,ie,j3a,j3,nps,mps,nms,mms,jca,jc,lcl,ndi,kmi,
     * ja,jai,nzei,nisi,lc,kk,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32
c
      dimension ae(ndtab)
c
c
*************************************************************************
* table-ci *
************
_IF(cray,ksr,i8)
      equivalence (t(50),ae(1))
_ELSE
      equivalence (t(26),ae(1))
_ENDIF
*************************************************************************
c
      do 29 l=1,mc
      kbab=jbab
      jbab=jbab+kml
      if (ifl.eq.1) go to 400
      ifl=1
      go to 30
 400  mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 30
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  30  go to (405,31,32,33), kk
  31  mm=mm+1
      ll=olab(mm)
c103  format(2x,20i6)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nnid) go to 34
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
34    mm=mm+1
      jj=olab(mm)
      if (mm.lt.nnid) go to 35
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 36
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto=nzero
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto=nzero
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=nzero
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=nzero
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
  49  continue
      do 5797 m=1,nzer
         f(m)=nzero
5797  continue
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
      kx=ja+my
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
43    lx=m-kml
      do 44 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
44    f(lx)=f(lx)+tim*ae(kx)
40    continue
c104  format(2x,10f12.8)
  47  ig7=kbab-jsec
      mx=ie
c244  format(2x,2i7,2f20.8)
      do 45 m=1,kml
      ly=0
      ig7=ig7+1
      in3=ig7
      my=mx
      do 45 if=1,kmj
      in3=in3+jsec
      tim=0.0d0
      mx=my
      do 5832 ig=1,kml
         ly=ly+1
         mx=mx+1
         tim=tim+f(ly)*ae(mx)
5832  continue
  45  hp(in3)=tim
      go to 29
 405  ig7=kbab-jsec
      do 406 m=1,kml
        ig7=ig7+1
        in3=ig7
        do 406 if=1,kmj
          in3=in3+jsec
          hp(in3)=0.0d0
 406  continue
      go to 29
  32  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 51
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
      kx=ja-my+256
      tim=-sm
      go to 59
58    tim=sm
      kx=my+ja
59    lx=m-kml
      do 60 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
60    f(lx)=f(lx)+tim*ae(kx)
56    continue
      go to 47
  33  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 61
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
61    mm=mm+1
      ntar=olab(mm)
      if (mm.lt.nnid) go to 62
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nnid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nnid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=pey(nb)
      jto=jto+1
      sm=sm+h(jto)
c247  format(3x,5i8,f20.8)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
c248  format(4x,4i8,f20.8)
c249  format(2x,30i4)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
  72  sac(m)=-sac(m)
  71  jj=jj+1
      do 5950 m=1,nzer
        f(m)=0.0d0
5950  continue
      jx=-kml
      do 73 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
c250  format(4x,6i8)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mps.eq.0) go to 78
      do 79 if=1,mps
      in=in+1
      im=ot(in)
79    tim=tim+sac(im)
78    lx=jx
c251  format(4x,i8,f20.8)
      do 81 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
81    f(lx)=f(lx)+tim*ae(kx)
      go to 73
 75   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
c252  format(4x,3i8,f20.8)
      if(mms.eq.0) go to 78
      do 83 if=1,mms
      in=in+1
      im=ot(in)
83    tim=tim-sac(im)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 84 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
84    f(lx)=f(lx)+tim*ae(kx)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 85 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
85    f(lx)=f(lx)+tim*ae(kx)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mms.eq.0) go to 92
      in=in+mps
      do 93 if=1,mms
      in=in+1
      im=ot(in)
93    tim=tim+sac(im)
92    lx=jx
      do 94 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
94    f(lx)=f(lx)+tim*ae(kx)
      go to 90
95    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mps.eq.0) go to 92
      in=in+mms
      do 99 if=1,mps
      in=in+1
      im=ot(in)
99    tim=tim-sac(im)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 100 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
100   f(lx)=f(lx)+tim*ae(kx)
      go to 90
97    do 6081 iif=1,nps
        in=in+1
        kx=ot(in)
        kx=kx+ja
        in=in+1
        im=ot(in)
        tim=-sac(im)
        lx=jx
        do 6080 ig=1,kmj
          kx=kx+ndj
          lx=lx+kml
          f(lx)=f(lx)+tim*ae(kx)
6080    continue
6081  continue
  90  jj=jj+jc
      go to 47
  29  continue
c246  format(2x,4i8)
      return
      end
****************************************************
      subroutine skin1b(pey,ideks,ndeks,
     +                  hp,nt64,f,ndimf,ot,nt444)
c
c ---- computes the matrix elements of the geyser-space for dk=-1
c
      implicit REAL (a-h,o-z)
c
INCLUDE(common/newmrd_parinc)
c
      REAL pey
      integer ndeks,ideks
      dimension ideks(ndeks)
      dimension pey(ndeks)
c
      REAL hp, f
      integer nt64, ndimf
      dimension hp(nt64),f(ndimf)
c
      integer ot
      dimension ot(nt444)
c --- commons
      common /a/    t(nlca)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
c
      integer ij, nplu
      common /scrtch/ ij(maxsym),nplu(ndk5)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      common /csac/ sac(nopmax+1)
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,lc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1, mx,md,lr,ipag
       common /d/ jbab,kml,ifl,iiz,nis,ii,lzer,km2,ndl,jsec,ie,iej,
     * la,kc,kmk,ndk,kr,iiq,ic1,ic2,i21,i22,i31,i32,iix,
     * kzer,ndt,iek,j3b,j3a,lps,nps,lms,nms,jcb,jca,lcl,ndi,kmi,
     * ia,jai,nzei,nisi,mc,kk,ndj,kmj,ja,nzer,mr,mps,mms,
     * ic,i2,i3,jc,j3
c
      dimension ae(ndtab)
c
c
************************************************************************
* table-ci  *
*************
_IF(cray,ksr,i8)
      equivalence (t(50),ae(1))
_ELSE
      equivalence (t(26),ae(1))
_ENDIF
************************************************************************
c
      do 29 l=1,mc
      kbab=jbab
      jbab=jbab+kmj
      if (ifl.eq.1) go to 400
      ifl=1
      go to 30
400   mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 30
c
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
 30   go to (405,31,32,33), kk
 31   mm=mm+1
      ll=olab(mm)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nnid) go to 34
      mm=nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
34    mm=mm+1
      jj=olab(mm)
      if (mm.lt.nnid) go to 35
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 36
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto=0
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=0
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
49    do 39 m=1,nzer
39    f(m)=0
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
      kx=ja+my
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
43    lx=m-kml
      do 44 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
44    f(lx)=f(lx)+tim*ae(kx)
40    continue
47    ig7=kbab-jsec
      mx=ie
      do 45 m=1,kml
      ly=0
      ig7=ig7+jsec
      in3=ig7
      my=mx
      do 45 if=1,kmj
      in3=in3+1
      tim=0
      mx=my
      do 46 ig=1,kml
      ly=ly+1
      mx=mx+1
46    tim=tim+f(ly)*ae(mx)
  45  hp(in3)=tim
      go to 29
405   ig7=kbab-jsec
      do 406 m=1,kml
      ig7=ig7+jsec
      in3=ig7
      do 406 if=1,kmj
      in3=in3+1
 406  hp(in3)=0
      go to 29
  32  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 51
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
      kx=ja-my+256
      tim=-sm
      go to 59
58    tim=sm
      kx=my+ja
59    lx=m-kml
      do 60 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
60    f(lx)=f(lx)+tim*ae(kx)
56    continue
      go to 47
  33  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 61
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
61    mm=mm+1
      ntar=olab(mm)
      if (mm.lt.nnid) go to 62
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nnid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nnid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=pey(nb)
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
      do 80 m=1,nzer
80    f(m)=0
      jx=-kml
      do 73 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mps.eq.0) go to 78
      do 79 if=1,mps
      in=in+1
      im=ot(in)
79    tim=tim+sac(im)
78    lx=jx
      do 81 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
81    f(lx)=f(lx)+tim*ae(kx)
      go to 73
 75   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mms.eq.0) go to 78
      do 83 if=1,mms
      in=in+1
      im=ot(in)
83    tim=tim-sac(im)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 84 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
84    f(lx)=f(lx)+tim*ae(kx)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 85 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
85    f(lx)=f(lx)+tim*ae(kx)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mms.eq.0) go to 92
      in=in+mps
      do 93 if=1,mms
      in=in+1
      im=ot(in)
93    tim=tim+sac(im)
92    lx=jx
      do 94 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
94    f(lx)=f(lx)+tim*ae(kx)
      go to 90
95    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mps.eq.0) go to 92
      in=in+mms
      do 99 if=1,mps
      in=in+1
      im=ot(in)
99    tim=tim-sac(im)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 100 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
100   f(lx)=f(lx)+tim*ae(kx)
      go to 90
97    do 101 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 101 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
101   f(lx)=f(lx)+tim*ae(kx)
90    jj=jj+jc
      go to 47
29    continue
      return
      end
*****************************************
      subroutine skin2a(hp,nt64,f,ndimf,ot,nt444)
c
c computes the matrix elements of the geyser-space for dk=+2
c
      implicit REAL (a-h,o-z)
c
      REAL hp, f
      integer nt64, ndimf
c
INCLUDE(common/newmrd_parinc)
c
      dimension hp(nt64), f(ndimf)
c
      integer ot
      dimension ot(nt444)
cdebug common /cdruck/idrk2
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
      common /a/ t(nlca)
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +mconf(ndk5),nconf(ndk5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
       common /d/ jbab,kmj,ifl,iiz,nis,iix,nzer,kml,ndt,jsec,iej,ie,
     . ia,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,ii,
     . kzer,ndj,iek,j3,j3a,mps,nps,mms,nms,jc,jca,lcl,ndi,kmi,
     . ja,jai,nzei,nisi,lc,ll,ndl,km2,la,lzer,lr,lps,lms,
     . ic2,i22,i32
c
      dimension ae(ndtab)
c
***********************************************************************
* table-ci ? *
**************
_IF(cray,ksr,i8)
      equivalence (t(50),ae(1))
_ELSE
      equivalence (t(26),ae(1))
_ENDIF
***********************************************************************
c
*     write(6,*) 'iwod=',iwod
cdebug
      do 13 l=1,mc
       kbab=jbab
       jbab=jbab+kml
       if (ifl.eq.1) go to 400
        ifl=1
        go to 14
 400   mm = mm + 1
       ll=olab(mm)
       if(mm.lt.nnid) go to 14
       mm=nzero
       read (ntape) olab8
       call unpack(olab8,8,olab,nnid)
14     if(ll.eq.0) go to 405
       mm=mm+1
       jj=olab(mm)
       if (mm.lt.nnid) go to 15
       mm=0
       read (ntape) olab8
       call unpack(olab8,8,olab,nnid)
15     jj=jj*ii+nis
       jto=jto+1
       ip=ot(jj)
       if(ip.eq.255) go to 803
       coul=h(jto)
       if (jto.lt.iwod) go to 16
       jto=0
       read(mtype) h
16     jto=jto+1
       exc=-h(jto)
       if (jto.lt.iwod) go to 805
       jto=0
       read(mtype)h
       go to 805
803    coul=-h(jto)
       if (jto.lt.iwod) go to 904
       jto=0
       read(mtype) h
904    jto=jto+1
       exc=h(jto)
       if(jto.lt.iwod) go to 805
       jto=0
       read(mtype) h
 805   if(ll-2) 800,801,802
 800   bob(1)=-coul-exc
       bob(2)=exc
       bob(3)=coul
       go to 804
 801   bob(1)=-exc
       bob(2)=coul+exc
       bob(3)=-coul
       go to 804
802    bob(1)=exc
       bob(2)=coul
       bob(3)=-coul-exc
 804   do 6599 m=1,nzer
         f(m)=0.0d0
6599   continue
c103   format(2x,10f12.8)
       do 18 m=1,kml
        jj=jj+1
        my=ot(jj)
        jj=jj+1
        if (my.eq.0) go to 18
        kx=ja+my
        mike=ot(jj)
        tim=bob(mike)
        lx=m-kml
        do 22 iff=1,kmj
         kx=kx+ndj
         lx=lx+kml
22      f(lx)=f(lx)+tim*ae(kx)
18     continue
       ig7=kbab-jsec
       mx=ie
       do 23 m=1,kml
        ly=0
        ig7=ig7+1
        my=mx
        in3=ig7
        do 23 iff=1,kmj
         in3=in3+jsec
         tim=0
         mx=my
         do 24 ig=1,kml
          ly=ly+1
          mx=mx+1
  24    tim=tim+f(ly)*ae(mx)
  23  hp(in3)=tim
cdebug  if (idrk2.gt.nzero) then
cdebug   write(6,*) 'skin2a'
cdebug   write(6,*) 'in3,hp',in3,hp(in3)
cdebug  end if
      go to 13
 405  ig7=kbab-jsec
      do 406 m=1,kml
       ig7=ig7+1
       in3=ig7
       do 406 if=1,kmj
        in3=in3+jsec
 406  hp(in3)=0.0d0
cdebug  if (idrk2.gt.nzero) then
cdebug   write(6,*) 'skin2a'
cdebug   write(6,*) 'in3,hp',in3,hp(in3)
cdebug  end if
  13  continue
c101  format(3x,20i6)
      return
      end
****************************************************
      subroutine skin2b(hp,nt64,f,ndimf,ot,nt444)
c
c   computes the matrix elements of the geyser-space for dk=-2
c
      implicit REAL (a-h,o-z)
      REAL hp, f
      integer nt64, ndimf
c
INCLUDE(common/newmrd_parinc)
c
      dimension hp(nt64), f(ndimf)
c
      integer ot
      dimension ot(nt444)
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
      common /a/  t(nlca)
c
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref), mn(maxref), nsc(maxref),
     +imo,lcl,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
c
       common /d/ jbab,kml,ifl,iiz,nisi,ii,nzei,kmi,ndi,jsec,ie,iej,
     * jai,kc,kmk,ndk,kr,iiq,ic1,ic,i21,i2,i31,i3,iix,
     * kzer,ndt,iek,j3,j3a,mps,nps,mms,nms,jc,jca,mc,ndj,kmj,
     * ia,ja,nzer,nis,lc,ll,ndl,km2,la,lzer,lr,lps,lms,
     * ic2,i22,i32
c
      dimension ae(ndtab)
c
c
*************************************************************************
* table-ci *
************
_IF(cray,ksr,i8)
      equivalence (t(50),ae(1))
_ELSE
      equivalence (t(26),ae(1))
_ENDIF
*************************************************************************
c
cdebug      if (idrk2.ge.1) then
cdebug       write(6,*)'iwod=',iwod
cdebug       write(6,*)'ae(1,5)',ae(1),ae(2),ae(3),ae(4),ae(5)
cdebug      end if
c
c
c
c
cdebug      if(idrk2.ge.1) then
cdebug        write(6,*)'jbab',jbab
cdebug        write(6,*)'jsec',jsec
cdebug        write(6,*)'nzer',nzer
cdebug        write(6,*)'nzei',nzei
cdebug        write(6,*)'kml ',kml
cdebug        write(6,*)'ndj ',ndj
cdebug        write(6,*)'kmj ',kmj
cdebug        write(6,*)'ia',ia
cdebug        write(6,*)'nis',nis
cdebug        write(6,*)'lc ',lc
cdebug        write(6,*)'ja',ja
cdebug        write(6,*)'ll',ll
cdebug        write(6,*)'ea',ae(1),ae(2),ae(3),ae(4),ae(5)
cdebug      end if
cdebug      if (idrk2.ge.1) then
cdebug        write(6,*) 'jbab    ',jbab
cdebug        write(6,*) 'kml     ',kml
cdebug        write(6,*) 'jsec    ',jsec
cdebug        write(6,*) 'nn      ',nn
cdebug        write(6,*) 'll=',ll,'    z-fall'
cdebug        write(6,*) 'nzer=',nzer , '# safs per sk'
cdebug      end if
cdebug
      do 13 l=1,mc
       kbab=jbab
       jbab=jbab+kmj
       if (ifl.ne.1) then
        ifl=1
        go to 14
      end if
      mm=mm+1
      ll=olab(mm)
cdebug      if (idrk2.ge.1) then
cdebug        write(6,*) 'll=',ll,'    z-fall'
cdebug        write(6,*) 'nzer=',nzer , '# safs per sk'
cdebug        write(6,*) 'kml =',kml  , '# safs per sk'
cdebug        write(6,*) 'kmj =',kmj  , '# safs per sk'
cdebug      end if
      if(mm.lt.nnid) go to 14
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
14    if(ll.eq.0) go to 405
      mm=mm+1
      jj=olab(mm)
      if (mm.lt.nnid) go to 15
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
15    jj=jj*ii+nis
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 803
      coul=h(jto)
      if (jto.lt.iwod) go to 16
      jto=nzero
      read(mtype) h
16    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 805
      jto=nzero
      read(mtype) h
      go to 805
803   coul=-h(jto)
      if (jto.lt.iwod) go to 904
      jto=nzero
      read(mtype) h
904   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 805
      jto=0
      read(mtype) h
c --- dk=2 case
 805  if(ll-2) 800,801,802
 800  bob(1)=-coul-exc
      bob(2)=exc
      bob(3)=coul
      go to 804
 801  bob(1)=-exc
      bob(2)=coul+exc
      bob(3)=-coul
      go to 804
 802  bob(1)=exc
      bob(2)=coul
      bob(3)=-coul-exc
 804  continue
*805  continue
c--- f77 : block if structure
*     if (ll.eq.1) then
*        bob(1)=-coul-exc
*        bob(2)=exc
*        bob(3)=coul
*     else if(ll.eq.2) then
*        bob(1)=-exc
*        bob(2)=coul+exc
*        bob(3)=-coul
*     else if(ll.eq.3) then
*        bob(1)=exc
*        bob(2)=coul
*        bob(3)=-coul-exc
*     else
*        write(6,*) ' error : ll=',ll
*     endif
* 804 continue
c --- f90 : case if structure
c90
*     case ll
*
*     end case
c ---
c
cdebug      do m=1,ndimf
      do 6825 m=1,nzer
        f(m) = 0.0d0
6825  continue
c
      do 18 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if (my.eq.0) go to 18
      kx=ja+my
      mike=ot(jj)
      tim=bob(mike)
cdebug      if (idrk2.ge.1) then
cdebug        write(6,*) 'mike=',mike,tim
cdebug      end if
      lx=m-kml
      do 6843 iff=1,kmj
        kx=kx+ndj
        lx=lx+kml
        f(lx)=f(lx)+tim*ae(kx)
6843  continue
  18  continue
      ig7=kbab-jsec
      mx=ie
      do 23 m=1,kml
       ly=nzero
       ig7=ig7+jsec
       my=mx
       in3=ig7
       do 23 iff=1,kmj
        in3=in3+1
        tim=0.0d0
        mx=my
c
        do 6861 ig=1,kml
         ly  = ly + 1
         mx  = mx + 1
         tim = tim+ f(ly)*ae(mx)
6861    continue
  23  hp(in3)=tim
cdebug      if (idrk2.gt.1) then
cdebug        write(6,*)'skin2b: vor goto 13'
cdebug        write(6,*)'in3, tim ',in3,tim
cdebug        write(6,*)'hp ',(hp(iii),iii=1,8)
cdebug      end if
       go to 13
 405  ig7=kbab-jsec
c
      do 6878 m=1,kml
        ig7=ig7+jsec
        in3=ig7
        do 6877 iff=1,kmj
          in3=in3+1
          hp(in3)=nzero
6877    continue
6878  continue
cdebug      if (idrk2.gt.1) then
cdebug        write(6,*)'skin2b: after goto 13'
cdebug        write(6,*)'in3, hp ',in3,hp(in3)
cdebug        write(6,*)'hp ',(hp(iii),iii=1,8)
cdebug      end if
  13  continue
c
      return
      end

****************************************************
      subroutine skina(pey,acoul,aexc,ideks,ndeks,
     +                 f,ndimf,w,mxref2,ot,nt444,
     +                 jkan,njkan)
c
c computing the matrix elements of the geyser-space
c
      implicit REAL (a-h,o-z)
c
      REAL pey,acoul,aexc
      integer ideks,ndeks
      dimension ideks(ndeks)
      dimension pey(ndeks),acoul(ndeks),aexc(ndeks)
c
      REAL f,w
      integer ndimf,mxref2
      dimension f(ndimf),w(mxref2)
      integer ot
      dimension ot(nt444)
      integer jkan, njkan
      dimension jkan(njkan)
c
      integer n
c
INCLUDE(common/newmrd_parinc)
c
c --- commons
c
_IF(notused)
      REAL ptiw
      common/junk/ ptiw(nt44r)
      parameter (nt44i = nt44 / 4)
      parameter (nt44r = nt44 / 8)
      logical*1 pt
      integer iw
      dimension iw(nt44i)
      dimension pt(nt44)
      equivalence (ptiw(1),iw(1),pt(1))
_ENDIF
c
      integer olab
_IF(cray,ksr,i8)
      integer olab8
_ENDIF
      common /clab/ olab(nnid),olab8(nnid8)
c
      integer nod
      common /cnod/ nod(ndk5)
      REAL h
      common /cbob/ bob(nn3),h(ndimh)
c
      integer ij, nplu
      common /scrtch/ij(maxsym),nplu(ndk5)
      integer itest
      common /cski1/itest(maxshl)
      common /ccore/core,isc(nn3)
c
      parameter (nmmo=256)
      integer loc, nir
      common /cloc/ loc(nmmo), nir(nmmo)
      integer nytl,ndub
      common /cny/  nytl(ndk5),ndub(ndk5)
      integer jdeks,kdeks
      common /cvt/ jdeks(nopmax),kdeks(nopmax)
      integer mn,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,
     +ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,jswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,
     +mconf(ndk5),nconf(ndk5),mh,iswh,
     +id(maxref), mn(maxref),nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
      common /bx/ jsec,kmax
c
c
********************************************************************
* table-ci *
**************
      REAL t, af, ae
      common /a/ t(nlca)
      integer j9
      dimension af(8010),j9(49),ae(ndtab)
      integer idra,idrc
      integer ndet,nsac
      integer jdrc,iaz,iez,jan,jbn,isc
      dimension ndet(ndk5)
      dimension nsac(ndk5)
      dimension idrc(ndk5)
      dimension idra(ndk5),jdrc(ndk5),jan(7)
      dimension jbn(7),iaz(ndk5),iez(ndk5)
      equivalence (af(1),t(1))
      equivalence (t(1),j9(1)),(j9(1),ndet(1)),(j9(6),nsac(1)),
     +(j9(11),iaz(1)),(j9(16),iez(1)),(j9(21),idra(1)),
     2(j9(26),idrc(1)),(j9(31),jdrc(1)),(j9(36),jan(1)),(j9(43),
     3jbn(1))
_IF(cray,ksr,i8)
      equivalence (af(50),ae(1))
_ELSE
      equivalence (af(26),ae(1))
_ENDIF
********************************************************************
      dimension sac(10)
      character*10 charwall
c
_IF(notused)
      equivalence (pt(1),ic)
      equivalence (pt(5),i2)
      equivalence (pt(9),i3)
_ELSE
      common/linkmr/ ic,i2,i3
      common/bufd/ gout(510),nword
      common/blksi3/nsz
      dimension lout(510)
      equivalence (gout(1), lout(1))
c
_ENDIF
c
      mult=nplu(1)+1
c
      ibl =nzero
      mq=ideks(jsec+1)
      call vclr(w,1,mq)
_IFN(cray)
      call setsto(nt444,0,ot)
_ENDIF
c
      ih = nzero
      do 6 i=1,iswh
      nc=nconf(i)
c@    write(6,*) 'aaaaaaaaaaaaaaaaaaaaaaaa'
      if (nc.gt.0) go to 7
      if (i.lt.3) go to 8
c@    write(6,*) 'bbbbbbbbbbbbbbbbbbbbbbbb'
      ibl=ibl+2
      go to 6
8     ibl=ibl+i-1
      go to 6
7     ndt=ndet(i)
c@    write(6,*) 'cccccccccccccccccccccccc'
      ia=iaz(i)
      kml=nsac(i)
      ie=iez(i)
      iz=isc(i)
      mzer=kml*kml
      nd=ndub(i)
      nps = nplu(i)
      nms = i - 1
      niw = nps*nms
      nr  = nod(i)
      nnx  = nytl(i)
      n1  = nr + 1
      n2  = nr + 2
      n3  = nr - 1
      np1 = nps- 1
      np2 = np1+ np1
      nd1 = nd - 1
      nm1 = i  - 2
      nm2 = nm1 + nm1
      ii  = kml + kml + 1
      nis = 1 - ii
      if (i.lt.3) go to 9
      ibl = ibl + 1
      j8  = i - 2
      mc  = nconf(j8)
      if(mc.eq.0) go to 10
c@    write(6,*) 'dddddddddddddddddddddddd'
      ndj=ndet(j8)
      kmj=nsac(j8)
c@    ne = ndj*kmj
c@    write(6,*) (ae(o),o=1,ne)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      js=jan(ibl)
      ks=jbn(ibl)
      ix=nad1
_IF(notused)
      do 11 j=1,ks
       ix = ix + nad
       js = js + 1
 11   call dreadp(pt(ix),js)
      call unpk14(pt,ix+nad-1,ot)
c --- test
c@    write(6,*) (ot(k),k=1,100)
_ELSE
      kss=0
  11  if (kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 11
      endif
_ENDIF
      kz=iz
      do 12 k=1,nc
      lz=jz
      do 13 l=1,mc
      mm=mm+1
      ll=olab(mm)
      if(mm.lt.nnid) go to 14
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
14    if(ll.eq.0) go to 13
      mm=mm+1
      jj=olab(mm)
      if (mm.lt.nnid) go to 15
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
15    jj=jj*ii+nis
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 803
      coul=h(jto)
      if (jto.lt.iwod) go to 16
      jto=0
      read(mtype) h
16    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 805
      jto=0
      read(mtype)h
      go to 805
803   coul=-h(jto)
      if (jto.lt.iwod) go to 904
      jto=0
      read(mtype) h
904   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 805
      jto = nzero
      read(mtype) h
 805  if(ll-2) 800,801,802
 800  bob(1)=-coul-exc
      bob(2)=exc
      bob(3)=coul
      go to 804
 801  bob(1)=-exc
      bob(2)=coul+exc
      bob(3)=-coul
      go to 804
 802  bob(1)=exc
      bob(2)=coul
      bob(3)=-coul-exc
804   do 905 m=1,nzer
        f(m)=0.0d0
905   continue
      do 18 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if (my.eq.0) go to 18
      kx=ja+my
      mike=ot(jj)
      tim=bob(mike)
      lx=m-kml
      do 4364 if0=1,kmj
        kx=kx+ndj
        lx=lx+kml
        f(lx)=f(lx)+tim*ae(kx)
4364  continue
  18  continue
      mz=kz
      mx=ie
      do 23 m=1,kml
      ly = nzero
      mz=mz+1
      my=mx
      nz=lz
      do 23 if0=1,kmj
      nz=nz+1
      tim=0
      mx=my
      do 24 ig=1,kml
      ly=ly+1
      mx=mx+1
24    tim=tim+f(ly)*ae(mx)
      if (dabs(tim).lt.1.0e-7) go to 23
      jg = ideks(mz)+nz
      w(jg)=tim
23    continue
13    lz=lz+kmj
12    kz=kz+kml
      go to 10
9     if (i.eq.1) go to 25
10    ibl=ibl+1
      j8=i-1
      mc=nconf(j8)
      if (mc.eq.0) go to 25
      js=jan(ibl)
      ks=jbn(ibl)
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      mps=nplu(j8)
      mms=j8-1
      ix=nad1
_IF(notused)
      do 27 j=1,ks
       ix=ix+nad
       js=js+1
 27   call dreadp(pt(ix),js)
      call unpk14(pt,ix+nad-1,ot)
_ELSE
      kss=0
 27   if (kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 27
      endif
      call upackx(ot)
_ENDIF
      jc=nr+mult
      j3=jc*kml+1
      kz=iz
      do 28 k=1,nc
      lz=jz
      do 29 l=1,mc
      mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 30
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
30    go to (29,31,32,33), kk
 31   mm=mm+1
      ll=olab(mm)
      if (ll.gt.128) ll=ll-256
      if (mm.lt.nnid) go to 34
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  34  mm=mm+1
      jj=olab(mm)
      if (mm.lt.nnid) go to 35
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
35    mm=mm+1
      kk=olab(mm)
      if (mm.lt.nnid) go to 36
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
36    jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      jto=jto+1
      if (ll.lt.0) go to 820
      coul = h(jto)
      if (jto.lt.iwod) go to 37
      jto = nzero
      read(mtype) h
37    jto=jto+1
      exc=-h(jto)
      if (jto.lt.iwod) go to 38
      jto = nzero
      read(mtype) h
      go to 38
820   coul=-h(jto)
      if (jto.lt.iwod) go to 821
      jto=0
      read(mtype) h
821   jto=jto+1
      ll=-ll
      exc=h(jto)
      if (jto.lt.iwod) go to 38
      jto=0
      read(mtype) h
38    if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
  49  continue
      do 4482 m=1,nzer
         f(m) = nzero
4482  continue
      do 40 m=1,kml
      jj=jj+1
      my=ot(jj)
      jj=jj+1
      if(my.eq.0) go to 40
      kx=ja+my
      mike=ot(jj)
      if(mike.gt.128) go to 42
      tim=bob(mike)
      go to 43
42    tim=-bob(256-mike)
43    lx=m-kml
c
      do 4501 if0=1,kmj
        kx    = kx + ndj
        lx    = lx + kml
        f(lx) = f(lx) + tim*ae(kx)
4501  continue
  40  continue
  47  mz = kz
      mx = ie
c244  format(2x,2i7,2f20.8)
      do 45 m=1,kml
        ly = nzero
        mz = mz + 1
        nz = lz
        my = mx
        do 45 if=1,kmj
          nz  = nz+1
          tim = 0.0d0
          mx  = my
          do 4518 ig=1,kml
            ly = ly + 1
            mx = mx + 1
            tim=tim+f(ly)*ae(mx)
4518      continue
        if (dabs(tim).lt.1.0e-7) go to 45
        jg = ideks(mz) + nz
        w(jg)=tim
  45  continue
      go to 29
  32  mm = mm + 1
      ll = olab(mm)
      if (mm.lt.nnid) go to 51
      mm = 0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
51    jto=jto+1
      sm=h(jto)
      if(jto.lt.iwod) go to 52
      jto=0
      read(mtype) h
52    if(ll.lt.128) go to 53
      sm=-sm
      ll=256-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0
      do 56 m=1,kml
      jj=jj+1
      my=ot(jj)
      if (my.eq.0) go to 56
      if(my.lt.128) go to 58
      kx=ja-my+256
      tim=-sm
      go to 59
58    tim=sm
      kx=my+ja
59    lx=m-kml
      do 60 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
60    f(lx)=f(lx)+tim*ae(kx)
56    continue
      go to 47
  33  mm = mm + 1
      ll = olab(mm)
      if (mm.lt.nnid) go to 61
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
61    mm = mm + 1
      ntar=olab(mm)
      if (mm.lt.nnid) go to 62
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
62    mm=mm+1
      nb=olab(mm)
      if(mm.lt.nnid) go to 63
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
63    mm=mm+1
      mb=olab(mm)
      if (mm.lt.nnid) go to 41
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  41  nb=ideks(nb)+mb+ij(ntar)
      sm=pey(nb)
      jto=jto+1
      sm=sm+h(jto)
c247  format(3x,5i8,f20.8)
      if(jto.lt.iwod) go to 64
      jto=0
      read(mtype) h
64    if (mr.eq.0 ) go to 65
      do 66 m=1,mr
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 68
      jto=0
      read(mtype) h
68    jto=jto+1
      sac(m)=-h(jto)
      if (jto.lt.iwod) go to 66
      jto=0
      read(mtype) h
66    continue
65    if(nd.eq.0) go to 67
      do 69 m=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 70
      jto=0
      read(mtype) h
70    jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 69
      jto=0
      read(mtype) h
69    continue
67    if(ll.gt.128) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
c248  format(4x,4i8,f20.8)
c249  format(2x,30i4)
      if(ip.eq.1) go to 71
      sm=-sm
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
      do 80 m=1,nzer
80    f(m)=0
      jx=-kml
      do 73 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
c250  format(4x,6i8)
      go to (74,75,76,77),im
 74   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mps.eq.0) go to 78
      do 79 if=1,mps
      in=in+1
      im=ot(in)
79    tim=tim+sac(im)
78    lx=jx
c251  format(4x,i8,f20.8)
      do 81 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
81    f(lx)=f(lx)+tim*ae(kx)
      go to 73
 75   in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
c252  format(4x,3i8,f20.8)
      if(mms.eq.0) go to 78
      do 83 if=1,mms
      in=in+1
      im=ot(in)
83    tim=tim-sac(im)
      go to 78
76    do 84 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 84 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
84    f(lx)=f(lx)+tim*ae(kx)
      go to 73
77    do 85 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 85 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
85    f(lx)=f(lx)+tim*ae(kx)
73    jj=jj+jc
      go to 47
86    jj=(256-ll)*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
      jx=jx+1
      in=jj
      im=ot(in)
      go to (91,95,96,97),im
91    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=sm
      if(mms.eq.0) go to 92
      in=in+mps
      do 93 if=1,mms
      in=in+1
      im=ot(in)
93    tim=tim+sac(im)
92    lx=jx
      do 94 if=1,kmj
      kx=kx+ndj
      lx=lx+kml
94    f(lx)=f(lx)+tim*ae(kx)
      go to 90
95    in=in+1
      kx=ot(in)
      kx=kx+ja
      tim=-sm
      if(mps.eq.0) go to 92
      in=in+mms
      do 99 if=1,mps
      in=in+1
      im=ot(in)
99    tim=tim-sac(im)
      go to 92
96    do 100 if=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=sac(im)
      lx=jx
      do 100 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
100   f(lx)=f(lx)+tim*ae(kx)
      go to 90
97    do 101 if=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ja
      in=in+1
      im=ot(in)
      tim=-sac(im)
      lx=jx
      do 101 ig=1,kmj
      kx=kx+ndj
      lx=lx+kml
101   f(lx)=f(lx)+tim*ae(kx)
90    jj=jj+jc
      go to 47
29    lz=lz+kmj
28    kz=kz+kml
25    js=idra(i)
      ks=idrc(i)
      ix=nad1
_IF(notused)
      do 102 j=1,ks
      ix=ix+nad
      js=js+1
 102  call dreadp(pt(ix),js)
      call unpk14(pt,ix+nad-1,ot)
_ELSE
      kss=0
  102 if (kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 102
      endif
      call upackx(ot)
_ENDIF
      kz=iz
      do 103 k=1,nc
      do 104 l=1,nnx
      ih=ih+1
      itest(l)=jkan(ih)
104   continue
      do 105 l=1,mzer
105   f(l)=0.0d0
      care=core
      if (nr.eq.0) go to 106
      do 107 l=1,nr
      mal=itest(l)
      ntar=nir(mal)
      mal=loc(mal)+1
      mal=ideks(mal)+ij(ntar)
107   care=care+pey(mal)
c235  format(10x,f20.8,2i6)
      if(nd.eq.0) go to 108
106   do 109 l=n1,nnx
        mal=itest(l)
        ntar=nir(mal)
        nal=ideks(mal+1)
        mal=loc(mal)+1
        mal=ideks(mal)+ij(ntar)
        sm=pey(mal)
109   care=care+sm+sm+acoul(nal)
108   if(nr.ge.2) then
        do 4806 l=2,nr
          mal = itest(l)
          mal = ideks(mal)
          it  = l - 1
          do 4805 m=1,it
            nal  = mal  + itest(m)
            care = care + acoul(nal)
4805      continue
4806    continue
      end if
      if(nd.lt.2) go to 112
      sm = 0
      do 113 l=n2,nnx
      mal=itest(l)
      mal=ideks(mal)
      it=l-1
      do 113 m=n1,it
      nal=mal+itest(m)
      tim=acoul(nal)
113   sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
112   if(nr.eq.0.or.nd.eq.0) go to 114
      do 115 l=1,nr
      mal=itest(l)
      do 115 m=n1,nnx
      nal=itest(m)
      if(nal.lt.mal) go to 116
      nal=ideks(nal)+mal
      go to 117
116   nal=ideks(mal)+nal
117   sm=acoul(nal)
115   care=care+sm+sm-aexc(nal)
c237  format(7x,f20.8,7i7)
114   kk=ic
      kn=i2
      ml=i3
      jx=-kml
c300  format(2x,6i10,f20.8)
      do 118 l=1,kml
c@    write(6,*) 'l,ohog(l),care'
c@    write(6,*) l,ohog(l),care
      jx=jx+1
      lx=jx
      kx=ot(l+n12)+ia
      tim=care
      if(nps.lt.2) go to 122
      mx=kn+1
      do 123 if=2,nps
      mx=mx+1
      it=if-1
      mal=ot(mx+n12)
c@    write(6,*) 'if,out(mx),mx'
c@    write(6,*) if,out(mx),mx
      mal=itest(mal)
      mal=ideks(mal)
      mz=kn
      do 123 in=1,it
      mz=mz+1
      nal=ot(mz+n12)
c@    write(6,*) 'in,mz,out(mz)'
c@    write(6,*) in,mz,out(mz)
      nal=itest(nal)+mal
123   tim=tim-aexc(nal)
      kn=mx
      if(nms.lt.2) go to 122
      mx=ml+1
      do 124 if=2,nms
      mx=mx+1
cs    mal=out(mx)
      mal=ot(mx+n12)
      mal=itest(mal)
      mal=ideks(mal)
      mz=ml
      it=if-1
      do 124 in=1,it
      mz=mz+1
csut  nal=out(mz)
      nal=ot(mz+n12)
      nal=itest(nal)+mal
124   tim=tim-aexc(nal)
      ml=mx
122   do 121 if=1,kml
      lx=lx+kml
      kx=kx+ndt
c@    write(6,*) 'if,f(lx),tim,ae(kx)'
c@    write(6,*) if,f(lx),tim,ae(kx)
121   f(lx)=f(lx)+tim*ae(kx)
c301  format(2x,3f20.8)
c906  format(2x,4f20.8)
c238  format(10x,2i6,3f20.8)
      if(niw.eq.0) go to 118
c240  format(2x,4i6,f20.8,2i6)
      do 120 m=1,niw
      kk=kk+1
csut  kx=out(kk)+ia
      kx=ot(kk+n12)+ia
      kk=kk+1
csut  mal=out(kk)
      mal=ot(kk+n12)
      kk=kk+1
csut  nal=out(kk)
      nal=ot(kk+n12)
      if (k.eq.0) write (6,241) kx,mal,nal,kk
      mal=itest(mal)
      nal=ideks(mal)+itest(nal)
      tim=-aexc(nal)
      lx=jx
 241  format(2x,21i6)
c242  format(2x,10f13.8)
      do 120 if=1,kml
      lx=lx+kml
      kx=kx+ndt
120   f(lx)=f(lx)+tim*ae(kx)
      if (k.eq.0) write (6,239) lx,kx,mal,nal,tim,f(lx),ae(kx)
239   format(2x,4i6,3f20.8)
118   continue
      mz=kz
      mx=ie
      do 125 l=1,kml
      ly=0
      nz=kz
      mz=mz+1
      my=mx
      do 125 m=1,l
      nz=nz+1
      tim=0
      mx=my
      do 126 if =1,kml
      ly=ly+1
      mx=mx+1
126   tim=tim+f(ly)*ae(mx)
      if(l.gt.m) go to 127
      jg=ideks(mz+1)
      w(jg)=tim
c178  format(3x,i10,f20.8,2i8)
      go to 125
127   if(dabs(tim).lt.1.0e-7) go to 125
      jg=ideks(mz)+nz
      w(jg)=tim
125   continue
103   kz=mz
      if(nc.eq.1) go to 6
      ks=jdrc(i)
      ix=nad1
_IF(notused)
      do 128 j=1,ks
      ix=ix+nad
      js=js+1
 128  call dreadp(pt(ix),js)
      call unpk14(pt,ix+nad-1,ot)
_ELSE
      kss = 0
 128  if (kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 128
      endif
      call upackx(ot)
_ENDIF
      ii=ii+kml
      j2=kml+1
      jc=nr+nr
      j3=jc*kml+1
      kz=iz+kml
      j8=0
      do 129 j=2,nc
      lz=iz
      j8=j8+1
      do 130 k=1,j8
      mm=mm+1
      kk=olab(mm)
c@    write(6,*) '******* kk ********'
c@    write(6,*) kk
      if (mm.lt.nnid) go to 131
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
c@    write(6,*) '******* mm ***131**'
c@    write(6,*) mm
 131  go to (130,132,133,134,135),kk
 135  mm=mm+1
      kk=olab(mm)
      if(mm.lt.nnid) go to 136
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
136   mm=mm+1
      ll=olab(mm)
      if(mm.lt.nnid) go to 137
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
137   if(kk.lt.ll) go to 138
      kk=ideks(kk)+ll
      go to 139
138   kk=ideks(ll)+kk
139   do 140 l=1,mzer
140   f(l)=0
      tim=aexc(kk)
      jx=-kml
      do 141 l=1,kml
csut  kx=ohog(l)+ia
      kx=ot(l+n12)+ia
      jx=jx+1
      lx=jx
      do 141 m=1,kml
      lx=lx+kml
      kx=kx+ndt
141   f(lx)=f(lx)+tim*ae(kx)
146   mz=kz
      mx=ie
      do 145 l=1,kml
      ly=0
      mz=mz+1
      my=mx
      nz=lz
      do 145 m=1,kml
      nz=nz+1
      tim=0
      mx=my
      do 143 if=1,kml
      ly=ly+1
      mx=mx+1
143   tim=tim+f(ly)*ae(mx)
      if(dabs(tim).lt.1.0e-7) go to 145
      jg=ideks(mz)+nz
      w(jg)=tim
145   continue
      go to 130
 144  mx=ie
      nzz=lz
      do 253 l=1,kml
      ly=0
      nzz=nzz+1
      my=mx
      mzz=kz
      do 253 m=1,kml
      mzz=mzz+1
      tim=0
      mx=my
      do 254 if=1,kml
      ly=ly+1
      mx=mx+1
 254  tim=tim+f(ly)*ae(mx)
      if (dabs(tim).lt.1.0e-7) go to 253
      jg=ideks(mzz)+nzz
      w(jg)=tim
 253  continue
      go to 130
 132  mm=mm+1
      ll=olab(mm)
      if(mm.lt.nnid) go to 147
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
147   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nnid) go to 148
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
148   mm=mm+1
      kk=olab(mm)
      if(mm.lt.nnid) go to 149
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
149   if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob=0
      go to 151
150   jj=ideks(kk)+jj
      ibob=1
151   jj=jj*ii+ic
      jto=jto+1
      ip=ot(jj)
      if(ip.eq.255) go to 830
      coul=h(jto)
      if(jto.lt.iwod) go to 152
      jto=0
      read(mtype) h
152   jto=jto+1
      exc=-h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read(mtype) h
      go to 153
830   coul=-h(jto)
      if (jto.lt.iwod) go to 831
      jto=0
      read (mtype) h
831   jto=jto+1
      exc=h(jto)
      if(jto.lt.iwod) go to 153
      jto=0
      read (mtype) h
  153 if (ll-2) 832,833,834
  832 bob(1)=coul+exc
      bob(2)=exc
      bob(3)=coul
      go to 835
  833 bob(1)=exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 835
 834  bob(1)=-exc
      bob(2)=-coul-exc
      bob(3)=coul
835   do 154 l=1,mzer
154   f(l)=0
      mx=-kml
      do 155 l=1,kml
      jj=jj+1
      mx=mx+1
      ip=ot(jj)
      if(ip.eq.2)  go to 159
      jj=jj+1
      kx=ot(jj)
      kx=kx+ia
      jj=jj+1
      tim=bob(1)
      lx=mx
      do 161 m=1,kml
      kx=kx+ndt
      lx=lx+kml
161   f(lx)=f(lx)+tim*ae(kx)
      go to 155
159   jj=jj+1
      kx=ot(jj)
      kx=kx+ia
      tim=bob(2)
      lx=mx
      do 160 m=1,kml
      kx=kx+ndt
      lx=lx+kml
160   f(lx)=f(lx)+tim*ae(kx)
      jj=jj+1
      kx=ot(jj)
      kx=kx+ia
      tim=bob(3)
      lx=mx
      do 167 m=1,kml
      lx=lx+kml
      kx=kx+ndt
167   f(lx)=f(lx)+tim*ae(kx)
155   continue
      if (ibob.eq.1) go to 144
      go to 146
 133  mm=mm+1
      ll=olab(mm)
      if (mm.lt.nnid) go to 168
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
168   mm=mm+1
      jj=olab(mm)
      if(mm.lt.nnid) go to 169
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
169   if(jj.gt.ll) go to 170
      jj=ideks(ll)+jj
      ibob=0
      go to 171
170   jj=ideks(jj)+ll
      ibob=1
171   jto=jto+1
      tim=h(jto)
      if(jto.lt.iwod) go to 172
      jto=0
      read(mtype)h
172   jj=jj*j2+i2
      do 177 l=1,mzer
177   f(l)=0
      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) tim=-tim
      do 173 l=1,kml
      jx=jx+1
      jj=jj+1
      kx=ot(jj)
      kx=kx+ia
      lx=jx
      do 173 m=1,kml
      kx=kx+ndt
      lx=lx+kml
173   f(lx)=f(lx)+tim*ae(kx)
      if (ibob.eq.1) go to 144
      go to 146
 134  mm=mm+1
c@    write(6,*) 'olab,mm,nnid 134'
c@    write(6,*) olab(mm),mm,nnid
      ll=olab(mm)
      if(mm.lt.nnid) go to 179
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
c@    write(6,*) '*******olab(1-nnid)*******'
c@    write(6,*) (olab(il),il=1,nnid)
179   mm=mm+1
c@    write(6,*) 'olab,mm,nnid 179'
c@    write(6,*) olab(mm),mm,nnid
      jj=olab(mm)
      if(mm.lt.nnid) go to 180
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
180   mm = mm + 1
c@    write(6,*) 'olab,mm,nnid 180'
c@    write(6,*) olab(mm),mm,nnid
      ntar = olab(mm)
      if(mm.lt.nnid) go to 183
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
183   mm = mm + 1
c@    write(6,*) 'olab,mm,nnid 183'
c@    write(6,*) olab(mm),mm,nnid
      nb = olab(mm)
      if(mm.lt.nnid) go to 184
      mm = nzero
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
184   mm=mm+1
c@    write(6,*) 'olab,mm,nnid 184'
c@    write(6,*) olab(mm),mm,nnid
      mb=olab(mm)
      if (mm.lt.nnid) go to 48
      mm=0
      read (ntape) olab8
      call unpack(olab8,8,olab,nnid)
  48  nb=ideks(nb)+mb+ij(ntar)
      sm=pey(nb)
c@    write(6,*) 'ideks(mb),mb,ij(ntar),ntar,pey(nb),ll,jj 48'
c@    write(6,*) ideks(mb),mb,ij(ntar),ntar,pey(nb),ll,jj
      if(ll.gt.128) go to 181
      if(ll.lt.jj) go to 182
      ibob=0
      jj=ideks(ll)+jj
      go to 205
182   jj=ideks(jj)+ll
      ibob=1
205   jj=jj*j3+i3
      if(nr.eq.1) go to 185
      do 186 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 187
      jto=0
      read(mtype) h
187   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 186
      jto=0
      read(mtype) h
186   continue
185   if(nd.eq.0) go to 188
      do 189 l=1,nd
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 190
      jto=0
      read(mtype) h
190   jto=jto+1
      sm=sm-h(jto)
      if(jto.lt.iwod) go to 189
      jto=0
      read(mtype) h
189   continue
188   ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      if (nr.eq.1) go to 191
      do 192 l=1,n3
192   sac(l)=-sac(l)
191   jj=jj+1
      do 193 l=1,mzer
193   f(l)=0
      jx=-kml
      do 194 l=1,kml
      in=jj
      im=ot(in)
      jx=jx+1
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      tim=sm
      if(im.eq.2) go to 195
      if(nps.eq.1) go to 196
      do 197 m=1,np1
      in=in+2
      nn=ot(in)
197   tim=tim+sac(nn)
196   do 198 m=1,kml
      lx=lx+kml
      kx=kx+ndt
198   f(lx)=f(lx)+tim*ae(kx)
      if(nms.eq.0) go to 194
      do 199 m=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
      do 199 if=1,kml
      lx=lx+kml
      kx=kx+ndt
199   f(lx)=f(lx)+tim*ae(kx)
      go to 194
195   if(nms.eq.1) go to 200
      do 201 m=1,nm1
      in=in+2
      nn=ot(in)
201   tim=tim+sac(nn)
200   do 202 m=1,kml
      lx=lx+kml
      kx=kx+ndt
202   f(lx)=f(lx)+tim*ae(kx)
      do 203 m=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
      do 203 if=1,kml
      lx=lx+kml
      kx=kx+ndt
203   f(lx)=f(lx)+tim*ae(kx)
194   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
181   ll=256-ll
      if(ll.lt.jj) go to 204
      ibob=0
      jj=ideks(ll)+jj
      go to 206
204   jj=ideks(jj)+ll
      ibob=1
206   jj=jj*j3+i3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 207
      jto=0
      read(mtype) h
207   jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 208
      jto=0
      read(mtype) h
208   if(nr.eq.1) go to 209
      do 211 l=1,n3
      jto=jto+1
      sm=sm+h(jto)
      if(jto.lt.iwod) go to 210
      jto=0
      read(mtype) h
210   jto=jto+1
      sac(l)=-h(jto)
      if(jto.lt.iwod) go to 211
      jto=0
      read(mtype) h
211   continue
209   if(nd.eq.1) go to 212
      do 213 l=1,nd1
      jto=jto+1
      tim=h(jto)
      sm=sm+tim+tim
      if(jto.lt.iwod) go to 214
      jto=0
      read(mtype) h
214   jto=jto+1
      sm=sm-h(jto)
      if (jto.lt.iwod) go to 213
      jto=0
      read(mtype) h
213   continue
212   ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      if(nr.eq.1) go to 215
      do 216 l=1,n3
216   sac(l)=-sac(l)
215   jj=jj+1
      do 217 l=1,mzer
217   f(l)=0
      jx=-kml
      do 218 l=1,kml
      in=jj
      im=ot(in)
      jx=jx+1
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      tim=-sm
      if (im.eq.2) go to 219
      if(nms.eq.0) go to 220
      in=in+np2
      jn=in
      do 221 m=1,nms
      jn=jn+2
      nn=ot(jn)
221   tim=tim-sac(nn)
220   do 222 m=1,kml
      lx=lx+kml
      kx=kx+ndt
222   f(lx)=f(lx)+tim*ae(kx)
      if(nms.eq.0) go to 218
      do 223 m=1,nms
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
      do 223 if =1,kml
      lx=lx+kml
      kx=kx+ndt
223   f(lx)=f(lx)+tim*ae(kx)
      go to 218
219   in=in+nm2
      jn=in
      do 224 m=1,nps
      jn=jn+2
      nn=ot(jn)
224   tim=tim-sac(nn)
      do 225 m=1,kml
      lx=lx+kml
      kx=kx+ndt
225   f(lx)=f(lx)+tim*ae(kx)
      do 226 m=1,nps
      in=in+1
      kx=ot(in)
      kx=kx+ia
      lx=jx
      in=in+1
      mx=ot(in)
      tim=sac(mx)
      do 226 if=1,kml
      lx=lx+kml
      kx=kx+ndt
226   f(lx)=f(lx)+tim*ae(kx)
218   jj=jj+jc
      if (ibob.eq.1) go to 144
      go to 146
130   lz = lz + kml
129   kz = kz + kml
6     continue
_IF(notused)
c     call drclos
_ELSE
      call clredx
      call closbf3
      cpu=cpulft(1)
       write (6,8000) cpu ,charwall()
 8000 format(//' end of zero order ci at ',f8.2,' seconds',a10,' wall')
_ENDIF
c 50  format(10x,25i5)
c245  format(2x,6i7,f20.8)
c246  format(2x,4i8)
      return
      end
****************************************************
      subroutine stilet
************************************************************************
      implicit REAL (a-h,q-z),logical*1(p)
************************************************************************
c
c this is a scalar optimised version of the stilx program,
c which essentially performs an integral sort.
c Further optimalisation opportunities:
c  - other sort algorithm
c  - adapting the file lengths according to the way the computer 
c    accesses the discs.
c
c  (1992) h.u. suter
c
************************************************************************
INCLUDE(common/newmrd_parinc)
INCLUDE(common/iofile)
************************************************************************
c
      integer ie,nsc
      common /b/ iwod,lg,irsf,ntape,ntype,ical,ntab,nc,lt,lz,nnx,
     +nd,nston,mdisk,ideli,mtape,iswh,nr,mm,jblk,jto,
     +igmax,nr1,mtype,linf,kfile,
     +nawk,isec,nshl,
     +nsec,lulu,nrootx,lsng,nprin,ipt0,mconf(5),nconf(5),mh,jswh,
     +id(maxref), ie(maxref), nsc(maxref),
     +imo,mc,mz,ndx,nrx,ntx,ny,nt,iht,ihz,nz,mt,nzx,jy,
     +jx,nd1,nd2,md1,md2,nr3,nr4,mr1,mr2,nt1,nt2,mt1,mt2,ih1,ih2,ih3,
     +ih4,jh2,jh3,jh4,iab,nt3,iac,nr2,m,nad,nad1,ifr1,mx,md,mr,ipag
************************************************************************
      REAL y,g,h
      integer lc,ld,kc,kd 
       common/junk/g(n3zt)  , h(n3zt),
     +             y(nnid)  ,
     +             lc(n3zt) , ld(n3zt),
     +             kc(ndimh), kd(ndimh)
************************************************************************
      write(iwr,5525) lg
5525  format(1x,'number of integrals :',i9)
      call rewftn(mtype)
      if (jto.gt.0) jblk=jblk+1
      lgh=(lg-1)/igmax+1
      ib8=(lg-1)/nnid+1
      igma=igmax/nnid
************************************************************************
cdebug      write(iwr,*) 'jblk=',jblk,'    jto=',jto
      igy  = igmax/iwod
      nres = ib8 - (lgh-1)*igma
      lx   = lg  - (ib8-1)*nnid
      iga  = (jblk-1)/igy
      jblk = jblk - iga*igy
      iga  = iga + 1
      igz  = igy
      iw   = iwod
      if(jto.eq.0) jto=iwod
************************************************************************
      do 5584 i=1,iga-1
        ig = nzero
c --- lade arrays lc und ld
        do 5523 j=1,igz-1
          read(ntype) kc,kd
          do 5522 k=1,iw
            ig    = ig + 1
            lc(ig)= kc(k)
            ld(ig)= kd(k)
5522      continue
5523    continue
        j = igz
        read(ntype) kc,kd
c ---------------
        do 5531 k=1,iw
          ig     = ig + 1
          lc(ig) = kc(k)
          ld(ig) = kd(k)
5531    continue
        call rewftn(nston)
        read (nston,err=61000,end=61001)
        read (nston,err=61000,end=61001)
        nrecs = 2
        nnx = igma
        lz = nnid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 5555 k=1,lgh-1
         il = nzero
         do 5546 l=1,nnx
          read(nston,err=61000,end=61001)y
          nrecs = nrecs + 1
          do 5545 ll=1,lz
           il=il+1
           g(il)=y(ll)
5545      continue
5546     continue
c
         do 5553 l=1,ig
            if (lc(l).eq.k) then
              kx=ld(l)
              h(l)=g(kx)
            end if
5553     continue
c
5555    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        k=lgh
        nnx=nres
        il=0
        do 5567 l=1,nnx
         read(nston,err=61000,end=61001)y
          nrecs = nrecs + 1
         if (l.eq.nres) lz=lx
         do 5566 ll=1,lz
           il=il+1
           g(il)=y(ll)
5566     continue
5567    continue
        do 5573 l=1,ig
           if (lc(l).eq.k) then
             kx=ld(l)
             h(l)=g(kx)
           end if
5573    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (ig.eq.igmax) then
          igy = (ig-1)/iwod+1
        end if
        la = nzero
        do 5583 k=1,igy
         lb=la+1
         la=la+iwod
         write(mtype)(h(l),l=lb,la)
5583    continue
5584  continue
************************************************************************
      i   = iga
      igz = jblk
      ig  = nzero
c --- load arrays lc and ld
      do 5595 j=1,igz-1
        read(ntype) kc,kd
        do 5594 k=1,iw
          ig = ig + 1
          lc(ig) = kc(k)
          ld(ig) = kd(k)
5594    continue
5595  continue
c
      j  =  igz
      iw =  jto
      read(ntype) kc,kd
      do 5604 k=1,iw
        ig=ig+1
        lc(ig)=kc(k)
        ld(ig)=kd(k)
5604  continue
c
      call rewftn(nston)
      nrecs = 2
      read(nston,err=61000,end=61001)
      read(nston,err=61000,end=61001)
      nnx=igma
      lz=nnid
*****************************************************************
      do 5628 k=1,lgh-1
        il=0
        do 5620 l=1,nnx
         read(nston,err=61000,end=61001)y
         nrecs = nrecs + 1
         do 5619 ll=1,lz
          il=il+1
          g(il)=y(ll)
5619     continue
5620    continue
c
        do 5627 l=1,ig
         if (lc(l).eq.k) then
           kx=ld(l)
           h(l)=g(kx)
         end if
5627    continue
5628  continue
*****************************************************************
      k  = lgh
      nnx = nres
      il = nzero
      do 5640 l=1,nnx
        read(nston,err=61000,end=61001)y
         nrecs = nrecs + 1
        if (l.eq.nres) lz=lx
        do 5639 ll=1,lz
         il=il+1
         g(il)=y(ll)
5639    continue
5640  continue
c
      do 5647 l=1,ig
         if (lc(l).eq.k) then
           kx=ld(l)
           h(l)=g(kx)
         end if
5647  continue
*****************************************************************
      if (ig.eq.igmax) then
        igy = (ig-1)/iwod+1
      end if
      la = nzero
c
      do 5658 k=1,igy
        lb = la + 1
        la = la + iwod
        write(mtype)(h(l),l=lb,la)
5658  continue
************************************************************************
      write(mtype) h
      write(mtype) h
      call rewftn(mtype)
************************************************************************
      if (lg-ib8*nnid.eq.0) then
        read (nston,err=61000,end=61001)
        nrecs = nrecs + 1
      endif
c
      write(iwr,61007) nston,nrecs
61007 format(/1x,
     + 'total number of 2e-integral records read from fort.',i2,
     + ' = ',i5/)
      return
61000 write(iwr,61005)
61005 format(1x,'error in reading transformed integral interface')
      call caserr(
     + 'error in reading transformed integral interface')
61001 write(iwr,61006)
61006 format(1x,'end of transformed integral interface')
      call caserr(
     + 'unexpected end of transformed integral interface')
      return
      end
      subroutine shutftn(ntape)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     close fortran unit ntape, delete the file
c     and reopen. This is to avoid keeping large
c     temporary files active during the MRDCI
c
      integer error_code
      character *132 filatt
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     *nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn),
     *iwhere(maxlfn),isync(maxlfn),
     *iposf(maxfrt),keepf(maxfrt),istatf(3,maxfrt),
     *oform(maxfrt)
INCLUDE(common/discc)
INCLUDE(common/utilc)
INCLUDE(common/iofile)
c
      if (ooprnt) then
       write (iwr,*) 'shut fortran stream', ntape
      endif
      if(keepf(ntape).gt.0) then
       close(unit=ntape, iostat=error_code, status='keep')
      else
       close(unit=ntape, iostat=error_code, status='delete')
      endif
      if(error_code.ne.0)then
         write(iwr,5)zftn(ntape)
 5      format(1x,'***** error closing file',1x,a8)
      endif
c
c     now re-open
c
      filatt=' '

c      call gtnv(zftn(ntape),filatt)
c      if (filatt.eq.' ') filatt = zftn(ntape)

cjmht psh change to get file routing to work with file directive
      call gtnv(zftfil(ntape),filatt)
      if (filatt.eq.' ') filatt = zftfil(ntape)

      do 10 loop = len(filatt),1,-1
        if(filatt(loop:loop).ne.' ') goto 20
 10   continue
 20   if(oform(ntape)) then
       open(unit=ntape,iostat=error_code,
     +      form='formatted',file=filatt,status='unknown')
      else
       open(unit=ntape,iostat=error_code,
     +      form='unformatted',file=filatt,status='unknown')
      endif
      if(error_code.ne.0)then
      write(iwr,3)yft(ntape),filatt
3     format(1x,'error opening fortran file',1x,a4,
     * ' : filename : ',a132)
      endif
c
      return
      end
_IF(notused)
      subroutine unpk14(in,n,out)
      logical*1 in(n)
      logical*1 out(4,*)
_IF(littleendian)
       do k =1, n
        out(1,k) = in(k)
       enddo
_ELSE
       do k =1, n
        out(4,k) = in(k)
       enddo
_ENDIF
      return
      end
_ENDIF
      subroutine ver_newmrd3(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd3.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
