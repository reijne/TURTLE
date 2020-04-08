c 
c  $Author: mrdj $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/casa.m,v $
c  $State: Exp $
c  
c     deck=casa
c ******************************************************
c ******************************************************
c             =   casa    =
c ******************************************************
c ******************************************************
      subroutine casscf(q)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/segm)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/machin)
INCLUDE(common/avstat)
INCLUDE(common/dm)
INCLUDE(common/harmon)
      common/scfopt/maxit(4),accdi(2),icoupl(4),dmpcut,
     * accc(2),etot
INCLUDE(common/finish)
      common/degen /ifdeg,idegen(100),ophase
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
      common/dims  /kkkk(12)
INCLUDE(common/exc)
      common/junke /maxt,ires,ipass,nteff,npass1,npass2,lentri,
     2                 nbuck,mloww,mhi,ntri,iacc
      common/atmol3/mina(2),mouta,moutb(3),
     *  numg(7),nsym,iorbsy(1)
      common/infoa/natt(3),num
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
      common/drtlnk/norbl,nsyml,ipassl,iprinl,labdrt(maxorb)
     * ,odrt, offmt0
      common/cone/m,na,nb,ispre(6),iperm(8,8)
INCLUDE(common/files)
      common/disc  /isel,iselr(2),irep,ichek,nt(maxlfn)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa
     1              ,nab,nst,nqe,isym,mults,macro,maxc,itypci
     2              ,itydrt
INCLUDE(common/caspar)
      common/l     /lock
INCLUDE(common/gjs)
      common/blkcore/corev(512),charge(10)
c
INCLUDE(common/runlab)
c
      dimension q(*),o1e(6)
      data m1/1/
      data m51,m500/51,500/
c
      nav = lenwrd()
c
c
      if(nprint.ne.-5)write(iwr,200)
 200  format(/40x,18('=')/
     = 40x,'casscf calculation'/40x,18('=')/)
c
      iblkq4=ibl3qa
      npass1=npas41
      npass2=npas42
      iacc=iacc4
c
      ncol4=num
      nsa4=num
      newb4=num
      nbas4=num
      ndump4=idaf
c     nblkq4=nbas4*ncol4 set to largest value (harmonic)
      nblkq4=newbas1*newbas1
      isecv4=mouta
c
      ma=newbas0
c
      call secget(isect(490),m51,ibloc)
c
      call readi(nirr,mach(13)*nav,ibloc,ndump4)
      odoubl=.true.
      ssz=1.3d0
c
      nsyml=nirr
      ntype=nirr
c
c     ----- core allocation
c
c     first determine required allocation
c
      i10=0
      i20=i10+nblkq4
      i30=i20+nblkq4
      last=i30+nblkq4
c
c     now determine offset and allocate
c
      i10 = igmem_alloc_all(last)
c
      i20=i10+nblkq4
      i30=i20+nblkq4
c
      call rdedx(q(i10),nblkq4,ibl3qa,idaf)
c
      if (newbas0.ne.newbas1) then
c
c...    compress vectors (harmonic)
c
         call comharm(q(i10),'vectors',ilifq)
         ncol4 = newbas0
         nsa4 = newbas0
         newb4 = newbas0
         nbas4 = newbas0
c
        
      end if
c
      do loop =1,6
       o1e(loop) = .false.
      enddo
      if(ophase) then
       o1e(2) = .true.
       call getmat(q(i20),q(i20),q(i20),q(i20),q(i20),q(i20),
     +             charge,nbas4,o1e,ions2)
       call mult2(q(i10),q(i30),q(i20),nbas4,nbas4,nbas4)
       call phase(q(i10),q(i30),ilifq,nbas4,nst,nprim,iwr,nprint)
       o1e(2) = .false.
      endif
c
c     restore s-matrix
c
      o1e(1) = .true.
      call getmat(q(i20),q(i20),q(i20),q(i20),q(i20),q(i20),
     +             charge,nbas4,o1e,ions2)
      call mult2(q(i10),q(i30),q(i20),nbas4,nbas4,nbas4)
      call orfog(q(i10),q(i10),q(i30),q(i20),
     * iky,ilifq,nbas4,nbas4,m1)
c
c     ----- orbsym invoked ?
c
      if(nsym.gt.0)call symorb(q(i10),iwr,nprint)
c
      call wrt3(q(i10),nblkq4,ibl3qa,idaf)
c
      call secput(isect(484),m500,m1,isecbl)
       call revind
      do 1320 i=1,nact
1320  icfcor(i+ncore)=iky(i)-ncore
      do 1330 i=1,maxorb
1330  ilifp(i)=(i-1)*nprim
      if(nprint.ne.-5)write(iwr,1370)
1370   format(/' 2-electron ao integral files'/1x,28('-'))
       if(nprint.ne.-5)call filprn(n1file,n1blk,n1last,n1tape)
       if(nprint.ne.-5)write(iwr,1380)
1380   format(/' secondary mainfile'/1x,18('-'))
       if(nprint.ne.-5)call filprn(n4file,n4blk,n4last,n4tape)
       if(nprint.ne.-5)write(iwr,1390)
1390   format(/' 2-electron mo integral files'/1x,28('-'))
       if(nprint.ne.-5)call filprn(n6file,n6blk,n6last,n6tape)
c
c     ----- set up base adress for casscf variable core
c     ----- ist free from preliminaries
c
      call gmem_free(i10)
c
c     now allocate all available memory
c
      ibase = igmem_alloc_all(lword4)
c
      if(lword4.lt.(nblkq4+lenb4+lenb4))call caserr(
     *'insufficient memory in casscf pre-processor')
      if(nwt.gt.kkkm)call caserr(
     + 'differing number of states in weight and strings')
      nmc=nprim
      ltemp=lock
      if(lock.eq.1)lock=nprim
      if(lock.eq.0)lock=ncore
      if(lock.eq.-1) lock = 0
      call symget(q(ibase),iwr,nprint)
      if(odrt)call guga(q(ibase),lword4)
c     ifclos=1
      call casdim(q(ibase),iwr)
c
      if(.not.oconv)etot=0.0d0
c
c     ----- reset for subsequent cycles
c
      if(zruntp.ne.'force') then
       odrt=.false.
       oconv=.false.
      else
       irestu=0
      endif
      lock=ltemp
      do 1420 i=1,maxorb
 1420 ilifm(i)=ilifq(i)
c
c     ----- reset core allocation
c
      call gmem_free(ibase)
c
      return
      end
      subroutine mastc(q,qfock,ccc,twoe,gam1,
     *gvec,mvec,km,nconf)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
c-----min. dimension = nact*(nact+1)/2
c     dimension ic2e(121),ic3e(120),ic4e(120),gam1(1770)
c
c-----min. dimension = total number of orbitals, ma
c     dimension ic1e(256,2),icf(256)
c
c-----dimension must be ntype - the dimension of symm group
c     dimension id(8),iperm(8,8),ims(8),nsymm(8),nblock(8,255)
c
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/popnos)
INCLUDE(common/avstat)
INCLUDE(common/simul)
      common/timez/timlim,date(5)
      common/disc  /isel,iselr,iselw,irep,icheck,ipos(maxlfn)
INCLUDE(common/dm)
INCLUDE(common/machin)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/cndx40)
INCLUDE(common/restar)
INCLUDE(common/iofile)
INCLUDE(common/ctrl)
      common/potn  /core,potnuc
      common/degen /ifdeg,idegen(100)
INCLUDE(common/exc)
      common/add   /iad1,iadf,iadd1,iadd2,iadci,iad2,
     1           iad105,iad11,iad12,iad13,iad14,iad15,iad16
     2              ,iad205,iad206,iad207,iad208
     3             ,iad311,iad312,iad313,iad314,iad315,iad316,iad317
     4             ,iad318,iad319,iad320
INCLUDE(common/caspar)
INCLUDE(common/finish)
      common/scfopt/maxit(4),accdi(2),icoupl(4),dmpcut
     *           ,acccc,enn,etotal,ehf
INCLUDE(common/prints)
      common/cone   /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),
     * ims(8)
      common/tab3/is1,is1p,is2,is2p
INCLUDE(common/gjs)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,ns
     1,n1e,isym,mults,macro,maxc,itypci,odrt
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/cipar/vlevc,accci
INCLUDE(common/ciconv)
      character*10 charwall
      dimension q(*),en(7),id(8)
      dimension qfock(*),ccc(nconf,*),twoe(*)
      dimension gam1(*),gvec(mvec,*)
c
      data m1,m10,m16/1,7,16/
      data timc,ifclos/0.0d0,1/
c
      timc=cpulft(1)
      nav = lenwrd()
c
      if(irest.ne.4)go to 11
c
      call secget(isect(494),m16,iblke)
      call rdedx(en,m10,iblke,ndump4)
      enn   =en(1)
      etotal=en(3)
      ehf   =en(2)
      go to 140
c
c.....start of iterations - 4-index first..
c
11    time=cpulft(1)
      if(nprint.ne.-5)write(iwr,90)m1,time
 10   call tractl(q,iwr,nprint)
c
c.....restore all-active integrals from disc
c
      call order(q(iad1),twoe,ccc,gam1,q(iad2),km,nconf)
c
c
      mark=1
      if(oprint(8))call fout(q(iad1),gam1(1),mark,iwr)
      if(oprint(8))write(iwr,20)
20    format(40x,'two electron integrals'/40x,22('-')/)
      if(oprint(8))write(iwr,30)(twoe(i),i=1,n2e)
30    format(8x,8f15.5)
c
c.....ci in active space
c
      idmgen=0
      call cicntl(q,gam1,twoe,ccc,
     *     gvec,mvec,idmgen,ifclos,km,nconf,iwr,nprint)
c
c
c.....ci hessian builder
c
      if(isimul(macro).ne.1)goto50
c
      call stabs(id,iperm,q(iad311),q(iad312))
      idumf = intnos(q(iad315),q(iad316))
      i=iad319-iad318+1
      call vclr(q(iad318),1,i)
      kk100=kk10+kk3-1
c
      call scanci(gam1,twoe,q(iad1),ccc,q(iad314),
     3q(iad317),q(iad318),q(iad311),q(iad312),
     4q(iad319),q(iad315),q(iad316),ifclos,q(iad320),nconf)
c
c
c.....orbital optimisation
c
50    call stabs(id,iperm,q(iad2),q(iad14))
      if(isimul(macro).ne.1)kk100=kk10
c
      if(mode.eq.0) then
      idum=nprim*nsec+ncore*nact+1
      idum=(idum-1)/nav+1
      ifffff = iad14 + idum
      call super(q,q(iad15),
     +q(iad11),q(iad2),q(iad14),q(iad1),q(iad12),q(iad13),
     +qfock,gam1,twoe,q(ifffff),iwr,nprint)
c
      else
c
      ireturn = 0
c
      call newton(q,q(iad2),q(iad14),q(iad1),q(iad15),
     + q(iad205),q(iad206),q(iad207),q(iad208),q(iad208),kk100,ccc,
     + qfock,gam1,twoe,nconf,ireturn)
c
c...   if ireturn is > 0 something went wrong => go to next cycle
c
        if (ireturn.gt.0) go to 10
c
      endif
c
      call gconv(nprim,gam1,twoe)
c
c
      if(isimul(macro).ne.1.or.oconv)goto80
c.....if simultaneous optimisation on this iteration, construct density
c     matrices for the revised ci vector
c
c
      idmgen=1
      call cicntl(q,gam1,twoe,ccc,
     *gvec,mvec,idmgen,ifclos,km,nconf,iwr,nprint)
c
80    mark=2
      if(oprint(8))call fout(q(iad1),gam1,mark,iwr)
c     mma=(nact+1)*nact/2
      time=cpulft(1)
      timel=timlim-time
      timc=(time-timc)*1.2d0
      call wrt3(q,nblkq4,iblkq4,ndump4)
      if(oconv)go to 120
      if(timc.gt.timel.or.macro.eq.maxc)go to 130
      timc=time
      macro=macro+1
      if(nprint.ne.-5)write(iwr,90)macro,time ,charwall()
 90   format(/1x,104('=')
     */40x,'start of cycle',i4,'  at time',f8.2,
     * ' seconds',a10,' wall'/40x,43('*')/)
      call clredx
      goto 10
 130  irest=3
      go to 100
c
c....analysis etc. of final wavefunction..
c
120   irest=0
100   call anal(q,ccc,q(iad1),qfock,gam1,twoe,nconf,nprint)
c
      call revise
      if(igrad.ne.1)return
 140  if(oconv)call grcntl(q,qfock,gam1,twoe,iwr)
      return
      end
      subroutine casdim(q,iwr)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb),len3(8),
     * ibext(8),lenext(8)
INCLUDE(common/machin)
INCLUDE(common/restar)
      integer  mxsurd
      parameter (mxsurd=3600)
      common/craypk/surd(mxsurd),nsurd,nbf,ncorf,nactf,ntypf,
     1itypf(maxorb),nwks,naaf,nabf,ispacf
INCLUDE(common/prints)
INCLUDE(common/mapper)
INCLUDE(common/simul)
INCLUDE(common/exc)
INCLUDE(common/caspar)
      common/add   /iadd1,iaddf,iaddd1,iaddd2,iaddci,iadd2,iad105,
     1iadd11,iadd12,iadd13,iadd14,iadd15,iadd16,
     2iad205,iad206,iad207,iad208
     3,iad311,iad312,iad313,iad314,iad315,iad316,iad317,iad318,iad319,
     4iad320
INCLUDE(common/cndx40)
      common/cone    /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8)
     *           , ims(8)
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
INCLUDE(common/qice)
INCLUDE(common/gjs)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,
     1nab,nst,n1e,isym,mults,macro,maxc,itypci
      common/block /length
INCLUDE(common/ctrl)
INCLUDE(common/harmon)
INCLUDE(common/discc)

      dimension q(*),id(8)
c
      nav = lenwrd()
c     nav16 = 16*nav
      na=naa
c
c
c....storage allocation for ci stage......
      nb=nab
      m=nact
      mm=m*(m+1)/2
      kk1=na*(m-na+1)
      kk1b=nb*(m-nb+1)
      if(kk1.lt.kk1b)kk1=kk1b
      kk2=na*(na-1)*(m-na+2)*(m-na+1)/4
      kk2b=nb*(nb-1)*(m-nb+2)*(m-nb+1)/4
      if(kk2.lt.kk2b)kk2=kk2b
      kk4=ma
      kk5=na+1
      kk7=nact+1
      kk8=ic2e(mm)+ic3e(mm)
      kk9=ic1e(nprimp,1)
      n8=na+1
      m=nact
      mf=m
      nm1a=na-1
      nm1b=nb-1
c
c
c
c     formula tape ci
c
      ipass=0
      ifillp=m9tape(1)
      iblk=m9blk(1)
      if(nprint.ne.-5)write(iwr,400)
 400  format(/1x,104('*')/ )
 260  ipass=ipass+1
      call rdedx(surd,mach(17),iblk,ifillp)
      if(nprint.ne.-5)write(iwr,270)yed(ifillp)
270   format(/' header block from formula tape ( ',a4,') read')
      if(nbf.ne.ma)call caserr('nbasis on formula tape is wrong')
      if(ncorf.ne.ncore)call caserr
     1                         ('ncore on formula tape is wrong')
      if(nactf.ne.nact)call caserr
     1                         ( 'nact on formula tape is wrong')
      if(ntypf.ne.ntype)then
       write(6,*)' ntypf, ntype = ', ntypf, ntype
       call caserr ('ntype on formula tape is wrong')
      endif
      do 290 i=1,nprim
      if(itypf(i).eq.isymmo(i))goto 290
      write(iwr,280)(j,isymmo(j),itypf(j),j=1,nprim)
280   format(' orbital',i3,' has symmetry',i2,' but formula tape'
     1,'was constructed with symmetry',i2)
      call caserr('orbital symmetries on formula tape are wrong')
290   continue
      if(ipass.eq.2)go to 300
      if(ifsimu.ne.1)go to 300
      ifillp=mttape(1)
      iblk =mtblk(1)
      go to 260
300   kk6=nwks
      kk3=nwks
c
      lenb4=ma*(ma+1)/2
c     nblkq4=ma*ma  leave at largest value (harmonic)
      nblkq4=newbas1*newbas1
      lensq=length*length
      lenb=nprim*nsec+ncore*nact
c
c
c ----
      maxvec=44
      maxv2=maxvec*maxvec
      maxv22=maxv2
      iadd1=nblkq4+1
      iaddf=iadd1+2*n1e
c...   max for case harmonic has reduced ma below newbas1 (num)
c...   this is accounted for in lenbn, used in assignment of iaddf below
      lenbn = max(lenb4,newbas1*(newbas1+1)/2)
      if(lenbn.gt.2*n1e)iaddf=iadd1+lenbn
c     if(lenb4.gt.2*n1e)iaddf=iadd1+lenb4
c ----- /fock/
      iaddd1=iaddf+nprim*nprim
c ----- /den1/
      iaddd2=iaddd1+kk9
c
      lenint=n2e+ic1e(nprimp,1)-ic1e(nst,1)
      if(lenint.lt.kk8)lenint=kk8
c ----- /junk/
      iaddci=iaddd2+lenint
c ----- /three/
      lencon=kk3*kkkm
      if(lencon.lt.maxv22)lencon=maxv22
      iadd2=iaddci+lencon
      iad105=iadd2
      nbasec=iad105+kk3
c
c
      lfree=lword4-nbasec
      ibuffc=lfree
      nbuffc=ibuffc/kk3
      mbuffc=(nbuffc-2)/2
      lfree=lfree-(kk3+kk8+kk9)
      if(lfree.le.0)go to 430
c
c......allocation for full hessian formation
      lenact=iky(nact+1)
      do 320 i=1,8
320   lena(i)=0
      do 330 i=1,lenact
330   lena(ic4e(i))=lena(ic4e(i))+1
      do 335 i=1,nact
      isymi=isymmo(i+ncore)
      len3(isymi)=0
      do 335 j=1,nact
 335  len3(isymi)=len3(isymi)+lena(iperm(isymmo(j+ncore),isymi))
      len8=0
      do 340 i=1,lenact
340   len8=len8+lena(ic4e(i))*2
      len8=len8-nact*lena(1)
      n12e=n2e+ic1e(nprimp,1)-ic1e(nst,1)
      if(len8.lt.n12e)len8=n12e
      iad311=nbasec+(kk2-1)/nav+1
      call stabs(id,iperm,q(iad311),q(iad311+lenb))
      iad312=iad311+(kk10-1)/nav+1
      kk10m=kk10-1
      if(nprint.ne.-5)write(iwr,350)lword4,kk3,kk10m
 350  format(/1x,'store allocation stastics'/1x,27('*')/
     1' size of available store =',i8/
     2' size of ci problem      =',i8/
     3' size of orbital problem =',i8/)
      kk100=kk10+kk3-1
      if(ifsimu.ne.1)go to 500
      iad313=iad312+(kk10-1)/nav+1
      iad314=iad313+kk9
      iad315=iad314+len8
      iad316=iad315+(nbas4*nprim-1)/nav+1
      iad317=iad316+iky(nprimp)
      iad318=iad317+n12e
      iad319=iad318+intnos(q(iad315),q(iad316))
      iad320=iad319+kk100
      lfreeh=lword4-(iad320+kk100)
      if(lfreeh.le.0.and.ifsimu.eq.1)goto 470
c
c.....allocation for super-ci
c
500   iadd14=iadd2+(kk10-1)/nav+1
      iadd15=iadd14+(kk10-1)/nav+1
      iadd11=iadd15+lensq
      iadd12=iadd11+nblkq4
      iadd13=iadd12+kk10+1
      iadd16=iadd13+kk10+1
      ibuffs=lword4-iadd16
      if(ibuffs.le.0)goto410
c
c....allocation for newton-raphson
c
      if(nrever.eq.0)go to 510
      iad205=iadd15+kk10+kk3+1
      iad206=iad205+kk10+kk3+1
      iad207=iad206+(kk10+kk3)/nav+1
      iad208=iad207+nblkq4+nbas4+1
      m6=iad207+lensq
      iad208=max(iad208,m6)
      ibuffn=lword4-iad208
      lfree=lword4-iad208-lensq
      m6=lword4-iad208-kk10*(kk10+1)/2
      lfree=min(lfree,m6)
      if(lfree.le.0)goto450
c
c
c
510   if(nprint.ne.-5)write(iwr,360)ibuffc,mbuffc
360   format(' amount of free store in ci stage =',i8,
     *1x,'( equivalent to',
     *i6,' vectors in davidson diagonalisation')
c
      if(ifsimu.eq.1.and.nprint.ne.-5)write(iwr,370)lfreeh
370   format(/' ci hessian builder will have',i8,' words ',
     1'more than it needs')
c
      if(nprint.ne.-5)write(iwr,380)ibuffs
380   format(/' super-ci will have',i8,' words more than it needs')
c
      if(nrever.ne.0.and.nprint.ne.-5)write(iwr,390)ibuffn
390   format(/' buffer size for newton-raphson routine is',i8,' words')
c
      if(irestu.eq.0)call cigues(q(iaddci),q(iaddd1),q(iaddd2),
     +                           kkkm,kk3)
      call mastc(q,q(iaddf),q(iaddci),q(iaddd2),q(iaddd1),q(iaddci),
     *maxvec,kkkm,kk3)
      return
c
c
c.....error messages
c
c
410   ibuffs=-ibuffs
      write(iwr,420)ibuffs
420   format(/' super-ci needs a further',i8,' words of core')
      goto 490
c
430   lfree=-lfree
      write(iwr,440)lfree
440   format(/' ci stage requires a further',i8,' words of core')
      go to 490
c
450   lfree=-lfree
      write(iwr,460)lfree
460   format(/' newton-raphson needs an extra',i8,' words of core')
      goto490
c
470   lfree=-lfreeh
      write(iwr,480)lfree
480   format(/' ci hessian builder needs a further',i8,' words of core')
      goto 490
c
490   call caserr('insufficient memory available')
      return
      end
      function intnos(ic7e,ic8e)
c
c     set up look up tables for single excitation integrals (ic7e)
c     and expanded form of 2pdm (ic8e) for use in construction of full
c     hessian. the size of the single excitation integral buffer is
c     returned via the function value.  lena(i) gives the number of
c     active-active orb pairs with symmetry i
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb),len3(8),ibext(8)
     +             ,lenext(8)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      common/cone   /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),ims(8)
INCLUDE(common/qice)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa
     1                  ,nab,nst,nqe,isym,mults,macros,maxc
      dimension ic7e(ma),ic8e(ma)
c
c.... set up the look-up for single-excitation integrals
c
c     len7=ma*nprim
c
      lenact=iky(nact+1)
      call setsto(8,0,lena)
      do 20 i=1,lenact
20    lena(ic4e(i))=lena(ic4e(i))+1
      intnos=0
      if(ncore.eq.0)goto 29
c.....(iu\vx)
      do 25 i=1,ncore
      ilifp(i)=(i-1)*nprim
      do 25 iu=nst,nprim
      ic7e(ilifp(i)+iu)=intnos
25    intnos=intnos+lena(iperm(isymmo(i),isymmo(iu)))
c.....(au\vx)
29    do 41 isyma=1,nirr
      ibext(isyma)=intnos+1
      lenext(isyma)=0
      do 40 i=nprimp,ma
      ilifp(i)=(i-1)*nprim
      if(isymmo(i).ne.isyma)goto 40
      lenext(isyma)=lenext(isyma)+1
      do 30 iu=nst,nprim
      ic7e(ilifp(i)+iu)=intnos
      isymiu=iperm(isymmo(i),isymmo(iu))
30    intnos=intnos+lena(isymiu)
40    continue
41    continue
      if(ncore.eq.0)goto 70
c.....a-i integrals for fock term
      do 50 i=1,ncore
      do 60 ia=nprimp,ma
      if(isymmo(i).ne.isymmo(ia))goto 60
      ic7e(ilifp(ia)+i)=intnos
      intnos=intnos+lena(1)
60    continue
50    continue
c
c.....set up expanded density matrix look-up
70    ic8e(1)=1
      do 80 i=1,nact
      ic8e(i+1)=ic8e(i)
      do 80 j=1,nact
80    ic8e(i+1)=ic8e(i+1)+lena(iperm(isymmo(i+ncore),isymmo(j+ncore)))
c
      return
      end
      subroutine cigues(ccc,gam1,gam2,km,nconf)
c
c....subroutine to generate a ci dump representing single determinant
c    trial wavefunctions using the information given by the strings
c    directive.
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/caspar)
INCLUDE(common/qice)
INCLUDE(common/exc)
INCLUDE(common/machin)
      common/atmol3/mina,minb,mouta,moutb
INCLUDE(common/mapper)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults,macro
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
     *, kk11,nroot
      common/blkin/kkbuf(12)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
      dimension kkout(12),ccc(nconf,*),gam1(*),gam2(*)
      equivalence (kkout(1),kk1)
      data m3,done/3,1.0d0/
c
c ..... allocate ci section on dumpfile
c
      nroot=km
      nav = lenwrd()
      nw=12
c...  take maximum size for section => newbas1
      lenv=1 + lensec(mach(8)) +lensec(mach(9)) + 
     1     lensec(newbas1*newbas1)
      kk3bl=lensec(kk3)
      kk9bl=lensec(kk9)
      kk8bl=lensec(kk8)
      len  = (kk3bl+kk8bl+kk9bl)*km + lenv + 1
      call secput(moutb,m3,len,iblk)
      iblkci=iblk+lenv
c
c ..... do we restore ci vectors ?
c
      if(minb.le.0)go to 800
      call secget(minb,m3,iblk)
      iblk=iblk+lenv
      call search(iblk,idaf)
      call readis(kkbuf,nw,idaf)
      if(kkbuf(3).ne.kk3) go to 200
      if(kkbuf(8).ne.kk8) go to 200
      if(kkbuf(9).ne.kk9) go to 200
      if(kkbuf(12).eq.km) go to 201
 200  call caserr(
     *' invalid ci section restored from dumpfile')
 201  do 90 kkk=1,km
      call reads(ccc(1,kkk),kk3,idaf)
      call reads(gam1,kk9,idaf)
 90   call reads(gam2,kk8,idaf)
c
c .... now load to moutb section
c
       call search(iblkci,idaf)
       call wrt3is(kkout,nw,idaf)
       iblkci=iblkci+1
       do 95 kkk=1,km
       call wrt3s(ccc(1,kkk),kk3,idaf)
       call wrt3s(gam1,kk9,idaf)
 95    call wrt3s(gam2,kk8,idaf)
       go to 100
 800  continue
      call vclr(ccc,1,lencon)
      call vclr(gam1,1,kk9)
      call vclr(gam2,1,kk8)
c
c.....guga trial vector
c
c
      call search(iblkci,idaf)
c
      call wrt3is(kkout,nw,idaf)
c
      iblkci=iblkci+1
c
      do 20 kkk=1,km
20    idom(kkk)=kkk
c
c
c
      do 70 kkk=1,km
c
      ccc(idom(kkk),kkk)=done
c
c.....construct 1-determinant density matrices
c
      do 60 i=1,nact
      m1=ic1e(i+ncore,1)+ic1e(i+ncore,2)
      itemp1=iocca(kkk,i)+ioccb(kkk,i)
      if(itemp1.eq.0)go to 60
      gam1(m1)= dfloat(itemp1)
c
      ii=iky(i)+i
      if(itemp1.eq.1)go to 40
      n1=ic2e(ii)+ic3e(ii)
      gam2(n1)=done
40    continue
c
      if(i.eq.1)go to 60
      im=i-1
      do 50 j=1,im
      ij=iky(i)+j
      jj=iky(j)+j
      itemp=iocca(kkk,j)+ioccb(kkk,j)
      if(itemp.eq.0)go to 50
      itemp=itemp1+itemp
      n1=ic2e(ii)+ic3e(jj)
      n2=ic2e(ij)+ic3e(ij)
      temp2=4.0d0
      temp3=-2.0d0
      if(itemp.lt.4)temp3=-1.0d0
      if(itemp.lt.4)temp2=2.0d0
      if(itemp.lt.3)temp2=1.0d0
      if(iocca(kkk,i).ne.iocca(kkk,j).and.itemp.eq.2)temp3=0.0d0
      gam2(n1)=temp2
      gam2(n2)=temp3
50    continue
60    continue
c
       call wrt3s(ccc(1,kkk),kk3,idaf)
       call wrt3s(gam1,kk9,idaf)
       call wrt3s(gam2,kk8,idaf)
70    continue
c
100   irestu=minb
      return
      end
      subroutine first
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
INCLUDE(common/prints)
INCLUDE(common/gjs)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,nst
     1,n1e,isym,mults,macro
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/cone   /m,n,mf,n8,nm1,nm2,iii(3),iperm(8,8),ims(1)
      dimension id(8)
c
      macro=1
      oprint(10)=.true.
c
c.....set up one-electron look-up table
      kt=0
      ic1e(1,1)=0
      call setsto(ntype,0,id)
      do 30 i=1,ma
      iz=isymmo(i)
      ic1e(i,1)=kt
      id(iz)=id(iz)+1
      ic1e(i,2)=id(iz)
      do 20 j=1,i
20    if (isymmo(i).eq.isymmo(j)) kt=kt+1
30    continue
      n1e=kt
c
      call setsto(ntype,0,id)
c
c-----ic4e to hold the symmetry of the kt'th ij pair,ic3e the number of
c     pairs so far with this symmetry
c
      kt=0
      do 50 i=1,nact
      do 50 j=1,i
      kt=kt+1
      iz=iperm(isymmo(i+ncore),isymmo(j+ncore))
      id(iz)=id(iz)+1
      ic4e(kt)=iz
50    ic3e(kt)=id(iz)
      kt1=kt
c     (kt is the total no.of ij pairs)
      mt=0
      do 70 i=1,kt1
      do 60 j=1,i
60    if (ic4e(i).eq.ic4e(j)) mt=mt+1
70    ic2e(i+1)=mt
      n2e=mt
      ic2e(1)=0
      return
      end
      subroutine fout(skv,gam1,mark,iwr)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
       common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,
     1nab,nst,n1e,isym,mults
INCLUDE(common/gjs)
      dimension skv(*),temp(8),itemp(8),gam1(*)
      if(mark.eq.0)write(iwr,10)
      if(mark.eq.1)write(iwr,20)
      if(mark.eq.2)write(iwr,30)
      if(mark.eq.3)write(iwr,40)
10    format(40x,'total  fock  matrix'            /40x,19('-'))
20    format(40x,'core  fock  matrix'             /40x,18('-'))
30    format(40x,'active  space  fock  matrix'    /40x,27('-'))
40    format(40x,'first  order  density  matrix'  /40x,29('-'))
      m=ma
      if(mark.eq.3)m=nprim
      mapl=m+1
      k=-8
50    k=k+8
      if(m-k.le.0)go to 51
      do 110 i=1,mapl
      l=8
      if(m-k.lt.8)l=m-k
      if(i.eq.1)go to 80
      do 60 j=1,l
      temp(j)=0.0d0
      if(isymmo(i-1).ne.isymmo(j+k)) go to 60
      m1=ic1e(max(i-1,j+k),1)+ic1e(min(i-1,j+k),2)
      if(mark.eq.0)temp(j)=skv(m1)+skv(n1e+m1)
      if(mark.eq.1)temp(j)=skv(m1)
      if(mark.eq.2)temp(j)=skv(n1e+m1)
      if(mark.eq.3)temp(j)=gam1(m1)
60    continue
      m1=i-1
      write(iwr,70)m1,(temp(j),j=1,l)
70    format(1x,i5,8f15.7)
      go to 110
80    do 90 j=1,l
90    itemp(j)=j+k
      write(iwr,100)(itemp(j),j=1,l)
100   format(//2x,8i15)
110   continue
      go to 50
 51   return
      end
      subroutine symorb(q,iwr,nprint)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gjs)
      common/atmol3/mina(6),numg(7),nsym,iorbsy(1)
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
      dimension q(*)
      data small/1.0d-6/
c
c.....get symmetries of trial mo's
      do 20 i=1,nbas4
      do 10 j=1,nbas4
      if( dabs(q(ilifq(i)+j)).gt.small)goto 20
10    continue
      call caserr('error in getting symmetries of trial mo''s')
20    isymmo(i)=isymao(j)
c
c.....process required symmetries
      if(nsym.eq.0)call caserr('no orbitals specified in orbsym')
      if(nsym.gt.nbas4)call caserr('too many orbitals given in orbsym')
      do 60 i=1,nsym
      is=iorbsy(i)
      if(is.eq.isymmo(i))goto 60
      ip=i+1
      if(i.eq.nbas4)goto 40
      do 30 j=ip,nbas4
      if(isymmo(j).eq.is)goto 50
30    continue
40    call caserr(
     + 'invalid no. of orbitals of given symmetry specified by orbsym')
50    continue
c
      call dswap(nbas4,q(ilifq(i)+1),1,q(ilifq(j)+1),1)
      isymmo(j)=isymmo(i)
      isymmo(i)=is
60    continue
      if(nprint.ne.-5)write(iwr,70)
70    format(/1x,'vectors processed by orbsym and written to dumpfile')
      return
      end
      subroutine dout(ma,mark,qmat,length,iwr)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qmat(length,*),itemp(8)
c
c.....write title.....
      if(mark.eq.0)write(iwr,10)
      if(mark.eq.1)write(iwr,11)
      if(mark.eq.2)write(iwr,12)
      if(mark.eq.4)write(iwr,14)
10    format(40x,'first order density matrix'/)
11    format(40x,'total fock matrix'/)
12    format(40x,'new molecular orbitals in basis of current mo''s'/)
c
14    format(/40x,'molecular orbital coefficients'/)
      mapl=ma+1
      k=-8
15    k=k+8
      if(ma-k.le.0)go to 30
      do 20 i=1,mapl
      l=8
      if(ma-k.lt.8)l=ma-k
      if(i.eq.1)go to 29
      m1=i-1
      lk=k+l
      kp=k+1
      write(iwr,201)m1,(qmat(m1,j),j=kp,lk)
201   format(i4,8f15.8)
      go to 20
29    do 27 j=1,l
27    itemp(j)=j+k
      write(iwr,200)(itemp(j),j=1,l)
200   format(//2x,8i15)
20    continue
      go to 15
30    write(iwr,300)
300   format(////1x)
      return
      end
      subroutine phase(q,smo,ilifq,nbasis,nst,nprim,iwr,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/degen/ifdeg,idegen(1)
      dimension q(*),smo(*),ilifq(*)
      data m0,m1,tol/0,1,1.0d-10/
c
c     subroutine to align the phases of orbitals of degenerate pairs
c     of symmetries.  smo  a 1-particle matrix over mo's.
c
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      if(ifdeg.ne.m1)return
      do 100 ix=nst,nprim
      if(idegen(ix).eq.m0)go to 100
      iy=idegen(ix)
c
      do 90 ixx=1,nprim
      indx=ind(ix,ixx)
      sx=smo(indx)
      absx= dabs(sx)
      if(absx.lt.tol)go to 90
c
      iyy=idegen(ixx)
      if(iyy.eq.m0)call caserr('invalid degeneracy spec. encountered')
      indy=ind(iy,iyy)
      sy=smo(indy)
      if(sx*sy.gt.0)go to 90
c
c.....we now have an out of phase set -  switch
c
      ilif=ilifq(iyy)
      do 85 j=1,nbasis
      ilif=ilif+m1
      q(ilif)=-q(ilif)
      indd=ind(j,iyy)
85    smo(indd)=-smo(indd)
c
90    continue
c
100   continue
c
      if(nprint.ne.-5)write(iwr,600)
600   format(' degenerate trial mos phase adjusted')
      return
      end
      subroutine stabs(id,iperm,ic5e,ic6e)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gjs)
INCLUDE(common/qice)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     +            naa,nab,nst,n1e,isym,mults
      common/dims/kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/tab3/is1,is1p,is2,is2p
      dimension iperm(8,*),id(*)
      dimension ic6e(*),ic5e(*)
c
      is1=1
      is2=1
      kk=1
      ic5e(1)=0
      ic6e(1)=0
      call setsto(ntype,0,id)
      if(ncore.eq.0)go to 30
      do 20 ip=1,ncore
      do 20 iq=nst,nprim
      if(iperm(isymmo(ip),isymmo(iq)).ne.1)go to 20
      kk=kk+1
      ic5e(kk)=ip
      ic6e(kk)=iq
20    continue
      is1=kk
30    do 50 ip=1,nprim
      do 40 iq=nprimp,ma
      if(iperm(isymmo(ip),isymmo(iq)).ne.1)go to 40
      kk=kk+1
      ic5e(kk)=ip
      ic6e(kk)=iq
40    continue
      if(ip.eq.ncore)is2=kk
50    continue
      kk10=kk
      is1p=is1+1
      is2p=is2+1
      return
      end
      subroutine symget(q,iwr,nprint)
c
c......read mo's from dumpfile and work out their symmetries....
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/exc)
      common/cone   /m,na,nb,ispre(6),iperm(8,8)
      common/block /length,nsymm(8),nblock(8,maxorb)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa
     1            ,nab,nst,nqe,isym,mults,macro,maxc
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
      dimension q(*)
      data small/1.d-6/
      call rdedx(q,nblkq4,iblkq4,ndump4)
      call setsto(8,-1,nfin)
      call setsto(8,0,mstart)
      do 30 i=1,nbas4
      do 20 j=1,nbas4
      if( dabs(q(ilifq(i)+j)).gt.small)go to 30
20    continue
      call caserr('error in assigning mo symmetries')
30    isymmo(i)=isymao(j)
      iprev=0
      do 50 i=1,nbas4
      if(isymao(i).eq.iprev)go to 50
      mstart(isymao(i))=i
      if(iprev.eq.0)go to 40
      nfin(iprev)=i-1
40    iprev=isymao(i)
50    continue
      nfin(iprev)=nbas4
c
c
      length=0
      do 60 i=1,nirr
      new=nfin(i)-mstart(i)+1
60    length=max(length,new)
      if(nprint.ne.-5)write(iwr,70)
70    format(/40x,'mo symmetries'/40x,13('*')/)
      if(nprint.ne.-5)write(iwr,80)(isymmo(i),i=1,ma)
80    format(10x,8i10)
c
      call setsto(ntype,0,nsymm)
      do 100 j=1,nbas4
      ibsym=isymmo(j)
      nsymm(ibsym)=nsymm(ibsym)+1
100   nblock(ibsym,nsymm(ibsym))=j
c
      call first
      do 120 kkk=1,kkkm
      ibsym=1
      do 110 i=1,nact
110   if(iocca(kkk,i)+ioccb(kkk,i).eq.1)ibsym=iperm(ibsym
     1,isymmo(i+ncore))
      if(ibsym.ne.isym.and.kkk.ne.1) call  caserr
     + ('states in strings directive have different space symmetries')
120   isym=ibsym
c
c     if(lock.eq.1)lock=nprim
c     if(lock.ne.nprim)lock=ncore
      if(nprint.ne.-5)write(iwr,130)isym
130   format(/' space symmetry of wavefunction is',i4)
      return
      end
      subroutine blockk
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
INCLUDE(common/cndx40)
      common/blkin/gout(510),mword
      common/disc/isel(5),ipos(maxlfn)
INCLUDE(common/restar)
      data m511/511/
      if(ipos(junits).ne.jblkas)call search(jblkas,junits)
      call put(gout(1),m511,junits)
      mword=0
      jblkas=jblkas+1
      jblkrs=jblkrs+1
      if(jblkrs)41,40,41
c...   change channel
 40    m4file=nfiles
       m4tape(m4file)=n4tape(nfiles)
       m4blk (m4file)= n4blk( nfiles)
       m4last(m4file)= jblkas
       nfiles=nfiles+1
      if(nfiles.gt.n4file)call caserr(
     *'insufficient no. of integral sections allocated')
      junits=n4tape(nfiles)
      jblkas=n4blk(nfiles)
      jblkrs=jblkas - n4last(nfiles)
  41  return
      end
      subroutine blockl
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
INCLUDE(common/cndx40)
      common/blkin/gout(510),mword
      common/disc/isel(5),ipos(maxlfn)
INCLUDE(common/restar)
      data m511/511/
      if(ipos(junitf).ne.jblkaf)call search(jblkaf,junitf)
      call put(gout(1),m511,junitf)
      mword=0
      jblkaf=jblkaf+1
      jblkrf=jblkrf+1
      if(jblkrf)41,40,41
 40   m6file=nfilef
      m6tape(m6file)=n6tape(nfilef)
      m6blk(m6file) =n6blk(  nfilef)
      m6last(m6file)=jblkaf
      nfilef=nfilef+1
      if(nfilef.gt.n6file)call caserr(
     *'insufficient no. of integral sections allocated')
      junitf=n6tape(nfilef)
      jblkaf=n6blk(nfilef)
      jblkrf=jblkaf-n6last(nfilef)
41    return
      end
      subroutine dpconv(nprim,gam1,gam2)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
INCLUDE(common/caspar)
      common/qpar/m8,n2e,ntype,ncore,nact
INCLUDE(common/mapper)
INCLUDE(common/gjs)
c
c.....routine converts density matrices gamma to d,p
c
      dimension gam1(*),gam2(*)
      if(.not.odoubl)return
      odoubl=.false.
      dumy=0.5d0
      go to 10
      entry gconv(nprim,gam1,gam2)
c
c.....converts d,p to gammas
c
      if(odoubl)return
      odoubl=.true.
      dumy=2.0d0
10    continue
      do 30 i=1,nact
      do 30 j=1,i
      ij=iky(i)+j
      do 30 k=1,i
      l1=k
      if(k.eq.i)l1=j
      do 20 l=1,l1
      kl=iky(k)+l
      if(ic4e(kl).ne.ic4e(ij))go to 20
      n=ic2e(ij)+ic3e(kl)
      dumx=1.0d0
      if(i.ne.j)dumx=dumy
      if(k.ne.l)dumx=dumx*dumy
      if(i.ne.k.or.j.ne.l)dumx=dumx*dumy
      gam2(n)=gam2(n)*dumx
20    continue
30    continue
c
      do 40 i=2,nprim
      im=i-1
      do 40 j=1,im
      if(isymmo(i).ne.isymmo(j))go to 40
      m=ic1e(i,1)+ic1e(j,2)
      gam1(m)=gam1(m)*dumy
40    continue
      return
      end
      subroutine trnsfo(q,work,dsx,iset,length)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/block /leng,nsymm(8),nblock(8,maxorb)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      dimension dsx(length,*),work(length,*),q(*)
c
      nprev=mstart(iset)-1
      nt=nsymm(iset)
      do 10 i=1,nt
      map=nprev+ilifq(nblock(iset,i))
      do 10 j=1,nt
10    work(i,j)=q(map+j)
c
      do 20 i=1,nt
      map=nprev+ilifq(nblock(iset,i))
      do 20 j=1,nt
      q(map+j)=ddot(nt,work(1,j),1,dsx(1,i),1)
 20   continue
      return
      end
      subroutine rdfock(vector,d,f,skv,gam1)
c
c      construct the 'active' fock matrix by transforming
c      the active space density matrix to the atomic orbital
c      basis, scanning the atomic integrals, and then transforming
c      the resultant fock matrix back to the mo basis.
c
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
INCLUDE(common/cndx40)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     1            naa,nab,nst,n1e,isym,mults
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      dimension vector(ma,*),d(*),f(*),skv(*)
      dimension gam1(*)
      lenb2=lenb4+lenb4
c
      call vclr(d,1,lenb4)
      call vclr(f,1,lenb2)
      do 30 it=nst,nprim
      do 30 iu=nst,it
      if(isymmo(iu).ne.isymmo(it)) go to 30
      gam=gam1(ic1e(it,1)+ic1e(iu,2))
      kl1=0
      do 20 k=1,ma
      do 20 l=1,k
      kl1=kl1+1
20    d(kl1)=d(kl1)+gam*(vector(k,it)*vector(l,iu)+
     1vector(k,iu)*vector(l,it))
c
30    continue
c
c
      do 4 i=1,nbas4
      mm=iky(i+1)
    4 d(mm)=d(mm)*0.5d0
      call scang(f,d)
      call mult2(vector,d,f,ncol4,ncol4,nbas4)
      nn=0
      do 40 i=1,nbas4
      do 40 j=1,i
      nn=nn+1
      if(isymmo(i).ne.isymmo(j))go to 40
      ll=ic1e(i,1)+ic1e(j,2)
      skv(n1e+ll)=d(nn)*0.5d0
40    continue
      return
      end
      subroutine avge(gam1,gam2,buffer,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
INCLUDE(common/popnos)
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9
INCLUDE(common/avstat)
INCLUDE(common/iofile)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      common/degen /ifdeg,idegen(1)
      common/cone  /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),ims(8)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults,macro,maxc
c
      dimension gam1(*),gam2(*),buffer(*)
      data m0/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine to average first and second order density matrices
c       over selected orbitals detailed in idegen
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nwt.eq.1)goto50
c....average over states
      call search(iblkci,idaf)
      call vclr(gam1,1,kk9)
      call vclr(gam2,1,kk8)
      do 30 istat=1,nwt
      wt=weight(istat)
      call reads(buffer,kk3,idaf)
      call reads(buffer,kk9,idaf)
      call reads(buffer(kk9+1),kk8,idaf)
      call daxpy(kk9,wt,buffer,1,gam1,1)
      call daxpy(kk8,wt,buffer(kk9+1),1,gam2,1)
30    continue
      if(nprint.ne.-5)write(iwr,40)
40    format(/' ***** density matrices averaged over states')
50    do 110 i=nst,nprim
      if(idegen(i).eq.m0)go to 110
      j=idegen(i)
      isymi=isymmo(i)
      isymj=isymmo(j)
c     isymij=iperm(isymi,isymj)
c
      do 100 it=nst,nprim
      iu=idegen(it)
      if(iu.eq.m0)go to 100
      isymt=isymmo(it)
      if(isymt.ne.isymi)go to 60
      iit=ic1e(max(it,i),1)+ic1e(min(it,i),2)
      iju=ic1e(max(iu,j),1)+ic1e(min(iu,j),2)
      gam1(iit)=(gam1(iit)+gam1(iju))*0.5d0
      gam1(iju)=gam1(iit)
c
c
60    continue
      isymu=isymmo(iu)
      isymit=iperm(isymt,isymi)
      isymju=iperm(isymj,isymu)
      if(isymju.ne.isymit)go to 100
      iit=iky(max(it,i)-ncore)+min(it,i)-ncore
      iju=iky(max(j,iu)-ncore)+min(j,iu)-ncore
c
      do 90 iv=nst,nprim
      isymv=isymmo(iv)
      isymtv=(iperm(isymt,isymv))
      isym2=iperm(isymtv,isymi)
      do 80 ix=nst,nprim
      if(isym2.ne.isymmo(ix))go to 80
c
      iy=iv
      iz=ix
      if(idegen(iv).eq.m0)go to 70
      iy=idegen(iv)
      iz=idegen(ix)
c
c.....average p(t,i,v,x) and p(u,j,y,z)
70    n1=iky(max(iv,ix)-ncore)+min(iv,ix)-ncore
      n2=ic2e(max(iit,n1))+ic3e(min(iit,n1))
      n4=iky(max(iy,iz)-ncore)+min(iy,iz)-ncore
      n3=ic2e(max(n4,iju))+ic3e(min(n4,iju))
      gam2(n2)=(gam2(n2)+gam2(n3))*0.5d0
      gam2(n3)=gam2(n2)
c
c
c.....average p(t,v,i,x) and p(u,y,j,z)
      n1=iky(max(it,iv)-ncore)+min(it,iv)-ncore
      n2=iky(max(i,ix)-ncore)+min(i,ix)-ncore
      n3=iky(max(iu,iy)-ncore)+min(iu,iy)-ncore
      n4=iky(max(j,iz)-ncore)+min(j,iz)-ncore
      n1=ic2e(max(n1,n2))+ic3e(min(n1,n2))
      n2=ic2e(max(n3,n4))+ic3e(min(n3,n4))
      gam2(n1)=(gam2(n1)+gam2(n2))*0.5d0
      gam2(n2)=gam2(n1)
c
c
80    continue
90    continue
100   continue
110   continue
      do 120 i=nst,nprim
      n1=ic1e(i,1)+ic1e(i,2)
120   pop(i)=gam1(n1)
      if(ifdeg.eq.1.and.nprint.ne.-5)write(iwr,130)
130   format(/' ***** density matrices symmetry averaged')
      return
      end
      subroutine cicntl(q,gam1,gam2,ccc,
     *gvec,maxv,idmgen,ifclos,km,nconf,iwr,nprint)
c
c
c     control routine for ci stage
c
c      limited option version for gamess ..... guga loop formulae tape
c
c       itypci = 5    guga, loop formula tape,        ---"---
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/cipar/vlev,accci
INCLUDE(common/caspar)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1              nst,n1e,isym,mults,macro,maxc,itypci
INCLUDE(common/avstat)
INCLUDE(common/simul)
INCLUDE(common/exc)
      common/degen /ifdeg,idegen(100)
      common/add   /iad1,iadf,iadd1,iadd2,iadci,iad2,
     1          iad105,iad11,iad12,iad13,iad14,iad15,iad16
INCLUDE(common/prints)
      common/cone  /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm,ims
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
INCLUDE(common/ciconv)
      common/atmol3/mina(3),moutb
      character*10 charwall
      dimension q(*),gam1(*),gam2(*),ccc(nconf,*),gvec(maxv,*)
c
      ifclos=-1
c
      t=cpulft(1)
      if(nprint.eq.-5)go to 1
      if(idmgen.ne.0)go to 2
      if(noci(macro).ne.1)write(iwr,10)
 10   format(/1x,53('*')/
     *1x,'loop driven guga-ci : davidson diagonalization method'/
     *1x,53('*')/)
      go to 1
  2   write(iwr,20)t ,charwall()
20    format(/' commence density matrix construction at ',f8.2,
     *' seconds',a10,' wall')
c
 1    if(.not.ojust)go to 3
      if(fmax.ge.0.1d0) accci=0.005d0
      if(fmax.ge.0.01d0.and.fmax.lt.0.1d0)accci=0.001d0
      if(fmax.ge.0.001d0.and.fmax.lt.0.01d0)accci=0.0001d0
      if(fmax.ge.0.0001d0.and.fmax.lt.0.001d0)accci=0.00005d0
      if(fmax.lt.0.0001d0)accci=0.00001d0
      go to 4
  3    if(macro.le.10)accci=fudgit(macro)
c
c....loop driven guga from loop formula tape
c
 4    if(macro.gt.n20cas.or.noci(macro).ne.1.or.idmgen.eq.1)
     * call gugaci(gam2,gam1,gvec,
     1maxv,q(iad1),ccc,q(nbasec),q(iad105),idmgen,km,nconf)
c
c
      irestu=moutb
c
      if(ifdeg.eq.1.or.nwt.gt.1)call avge(gam1,gam2,q(nbasec),nprint)
c
      if(nprint.eq.-5)return
      t=cpulft(1)
      if(idmgen.eq.0)write(iwr,70)t ,charwall()
70    format (/' end of guga ci at',f8.2,' seconds',a10,' wall')
      if(idmgen.ne.0)write(iwr,80)t ,charwall()
 80   format(/' density matrix construction complete at ',
     *f8.2,' seconds',a10,' wall')
      return
      end
      subroutine order(skv,twoe,ccc,gam1,buff,km,nconf)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/craypk/i205(1360)
      common/actlen/lena(8),icfcor(100)
INCLUDE(common/iofile)
INCLUDE(common/exc)
INCLUDE(common/prints)
INCLUDE(common/filel)
INCLUDE(common/restar)
      common/potn  /core,potnuc
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1           nst,n1e,isym,mults,macro
      common/dims/kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/disc/isel,iselr,iselw,irep,ichek,ipos(maxlfn)
      common/cone  /m,na,nb,juju(6),iperm(8,8)
INCLUDE(common/simul)
      common/blkin /g(510),nword
INCLUDE(common/qice)
INCLUDE(common/atmblk)
      dimension ccc(nconf,*),skv(*),buff(*),twoe(*)
      dimension gam1(*)
c
      mac=min(macro,n20cas)
      if(noci(mac).ne.0)go to 90
c
      call setsto(1360,0,i205)
c
      lfile=m6file
      do 40 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i)  =m6blk(i)
  40  llblk(i)  =liblk(i)-m6last(i)
c
      call vclr(twoe,1,n2e)
      oprint(10)=.false.
c
      do 50 ii=1,lfile
      idev=lotape(ii)
      call search(liblk(ii),idev)
      call find(idev)
      jj=llblk(ii)
20    jj=jj+1
      call get(g(1),nw)
      if(nw.le.0)goto 50
      if(jj.ne.0)call find(idev)
_IF1(iv)      call upak8v(g(num2e+1),i205)
_IFN1(iv)      call unpack(g(num2e+1),lab816,i205,numlab)
_IFN1(iv)      int4=1
      do 3 kk=1,nword
_IF(ibm,vax)
      i=i205(kk)
_ELSEIF(littleendian)
      i=i205(int4+1)
_ELSE
      i=i205(int4  )
_ENDIF

      if(i.gt.nprim)go to 3
      if(i.le.ncore)go to 3
_IF(ibm,vax)
       j=j205(kk)
_ELSEIF(littleendian)
       j=i205(int4  )
_ELSE
       j=i205(int4+1)
_ENDIF
      if(j.le.ncore)go to 3
_IF(ibm,vax)
      k=k205(kk)
_ELSEIF(littleendian)
      k=i205(int4+3)
_ELSE
      k=i205(int4+2)
_ENDIF
      if(k.le.ncore)go to 3
_IF(ibm,vax)
      l=l205(kk)
_ELSEIF(littleendian)
      l=i205(int4+2)
_ELSE
      l=i205(int4+3)
_ENDIF
      if(l.le.ncore)go to 3
      twoe(ic2e(icfcor(i)+j)+ic3e(icfcor(k)+l))=g(kk)
_IF1(iv)    3 continue
_IFN1(iv)    3 int4=int4+4

      if(jj)20,50,20
50    continue
      call search(iblkci,idaf)
      do 70 i=1,km
      call reads(ccc(1,i),kk3,idaf)
      call reads(buff(kk3+1),kk9,idaf)
      call reads(buff(kk3+kk9+1),kk8,idaf)
      vl(i)=
     *ddot(kk9,buff(kk3+1),1,skv,1)+ddot(kk8,buff(kk3+kk9+1),1,twoe,1)
 70   continue
      return
c
 90   call search(iblkci,idaf)
      do 100 i=1,km
      call reads(ccc(1,i),kk3,idaf)
      call reads(gam1,kk9,idaf)
      call reads(twoe,kk8,idaf)
 100  continue
c
      return
      end
      subroutine gugaci(twoe,gam1,g,maxv,
     *skv,ccc,buffer,indx,idmgen,km,nconf)
c
c
c     control routine for loop driven formula tape ci using
c     davidson/liu diagonaliser - excited states done simultaneously
c
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      integer ilifg(5)
      common/disc/isel,iselr,iselw,irep,icheck,ipos(maxlfn)
      integer    mxsurd
      parameter (mxsurd=3600)
      common/blkin /vsurd(mxsurd),nnsurd
INCLUDE(common/iofile)
      common/craypk/nuwks(340),nlwks(340),iuwks(340),juwks(340),
     1              ints(340),icoup(340),nword,icode,ispeer,nnword
INCLUDE(common/caspar)
INCLUDE(common/machin)
INCLUDE(common/popnos)
INCLUDE(common/qice)
INCLUDE(common/exc)
INCLUDE(common/prints)
      common/cipar /vlev,acc,time,nitc,nfudge,ncfl,ncfu,nprinc
INCLUDE(common/gjs)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults,macro,maxc,itypci
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,k7,kk8,kk9,kk10
      common/cone  /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),ims(8)
      common/potn  /core
INCLUDE(common/restar)
c
      common/junk/vlamda(44),scrap(44),surd(mxsurd),r(5),q(5),sout(5)
      dimension twoe(*),g(maxv,*),gam1(*),iy(2)
      dimension skv(*),buffer(nconf,*)
      dimension ccc(nconf,*),indx(*)
c
      equivalence (iy(1),vsurd(1))
c
      data ilifg/0,44,88,132,176/
      data dzero,sixty/0.0d0,60.0d0/
      data thresh/1.d-3/
      data m1,m3/1,3/
c
      nav = lenwrd()
      odmgen=idmgen.ne.0
      ifail=0
      ifillp=m9tape(1)
      iblk0=m9blk(1)
c
      dumf = cpulft(1)
c
c     ibd=nbuffc-1
      maxvec=(nbuffc-1)/2
      maxvec=min(maxvec,maxv,kk3)
c ## store diagonal elements in more sensible position
c ## this removes the dynamic error on C2
      ibd=maxvec+maxvec+1
      if(.not.odmgen.and.nprint.ne.-5)write(iwr,10)maxvec
10    format(' maximum of',i3,' expansion vectors')
      do 20 kkk=1,km
20     if(oprint(3))write(iwr,30)kkk,(ccc(ll,kkk),ll=1,kk3)
30    format(40x,'trial vector number',i2/100(/10f13.6)/)
c
c....read header block from loop tape
      call rdedx(vsurd,mach(17),iblk0,ifillp)
      nsurd=nnsurd
      call dcopy(nsurd,vsurd,1,surd,1)
c
c....read index vector
      k32=(kk3-1)/nav+1
      call readis(indx,k32*nav,ifillp)
      iblk = lensec(mach(17))
     *    +lensec(k32)   + iblk0
c
c.....copy 1-electron integrals into zint
      n11e=ic1e(nprimp,1)-ic1e(nst,1)
      call dcopy(n11e,skv(ic1e(nst,1)+1),1,twoe(n2e+1),1)
c
      nstart=1
      nend=km
c     nvec=km
c.....addressing of buffer is
c
c    |    c    |   sigma     | diag |
c
c....copy trial vectors up to buffer.
c
      do 60 kkk=1,km
c
c.....normalise trial vector
c
      call dcopy(kk3,ccc(1,kkk),1,buffer(1,kkk),1)
      s=1.0d0/dnrm2(kk3,buffer(1,kkk),1)
      call dscal(kk3,s,buffer(1,kkk),1)
 60   continue
c
c
      nits=nitc
      mtx=0
c
      if(odmgen)goto 610
c
c   ********************************************************************
c
c
c.....get the diagonal matrix elements
      call vclr(buffer(1,ibd),1,kk3)
      call search(iblk,ifillp)
      call find(ifillp)
c
70    call get(vsurd,nw)
      if(nw.eq.0)goto 150
      call find(ifillp)
c
c
c**************************************
c
c.....loop tape
_IF1(iv)      call upakin(iy,nuwks,511)
_IFN1(iv)      call unpack(iy,16,nuwks,2044)
      if(icode.eq.1.or.nword.eq.0)goto 70
c.....process block
      do 100 iword=1,nword
      nlwk=nlwks(iword)
      nuwk=nuwks(iword)
      iuwk=iuwks(iword)
      val=surd(icoup(iword))*twoe(ints(iword))
c.....loop over matrix elements to which this loop contributes
      do 90 i=1,nlwk
      ii=indx(iuwk)
      do 80 j=1,nuwk
      buffer(ii,ibd)=buffer(ii,ibd)+val
80    ii=ii+1
90    iuwk=iuwk+1
100   continue
      goto 70
c
c*********************************
c
150   call search(iblk,ifillp)
      call find(ifillp)
c
c
c
c  **************
c
c      code to gen. trial vectors
      if(irestu.ne.0)go to 210
c
      do 200 kkk=1,km
      idom(kkk)=0
      kkm=kkk-1
      accum=sixty
      do 180 ll=1,kk3
      if(kkk.eq.1)goto 170
      otest=.false.
      do 160 i=1,kkm
160   otest=otest.or.ll.eq.idom(kkm)
      if(otest)goto 180
170   if(accum.lt.buffer(ll,ibd))goto 180
      idom(kkk)=ll
      accum=buffer(ll,ibd)
180   continue
      smfg=accum+core
      if(nprint.ne.-5)write(iwr,190)kkk,idom(kkk),smfg
 190  format(/
     *' trial vector for state',i2,' re-selected as unit vector',
     1' number',i6,'       energy =',f15.8)
      call vclr(buffer(1,kkk),1,kk3)
200   buffer(idom(kkk),kkk)=1.0d0
210   continue
      if(nprint.ne.-5)write(iwr,370)
 370  format(/1x,43('*')/
     *' cycle    time   tester               energy'/
     *1x,43('*')/)
c
c
c  **************
c
c
c.....start of iteration
220   mtx=mtx+1
c
c********************************
c.....loop tape
c
      do 230 i=nstart,nend
      call vmul(buffer(1,ibd),1,buffer(1,i),1,buffer(1,maxvec+i),1,kk3)
 230  continue
c
240   call get(vsurd,nw)
      if(nw.eq.0) go to 350
      call find(ifillp)
_IF1(iv)      call upakin(iy,nuwks,511)
_IFN1(iv)      call unpack(iy,16,nuwks,2044)
      if(icode.eq.0.or.nword.eq.0)goto240
c.....process block
      do 280 iword=1,nword
      nlwk=nlwks(iword)
      nuwk=nuwks(iword)
      val=surd(icoup(iword))*twoe(ints(iword))
      do 270 istat=nstart,nend
      iuwk=iuwks(iword)
      juwk=juwks(iword)
      do 260 i=1,nlwk
      ii=indx(iuwk)
      jj=indx(juwk)
      do 250 j=1,nuwk
      buffer(jj,maxvec+istat)=buffer(jj,maxvec+istat)
     1             +val*buffer(ii,istat)
      buffer(ii,maxvec+istat)=buffer(ii,maxvec+istat)
     1             +val*buffer(jj,istat)
      ii=ii+1
250   jj=jj+1
      iuwk=iuwk+1
260   juwk=juwk+1
270   continue
280   continue
      goto 240
350   call search(iblk,ifillp)
      call find(ifillp)
c
c
c***********************************************************************
c
c......guts of davidson algorithm
c
c
c
c......form small (g) matrix
      do 360 i=1,nend
      do 360 j=1,i
      g(i,j)=ddot(kk3,buffer(1,maxvec+i),1,buffer(1,j),1)
 360  continue
c      call mxma (buffer(1,maxvec+1),kk3,1, buffer,1,kk3,
c     +            g,1,maxv,    nend,kk3,nend )
c
c.....diagonalise g
      call f02abf(g,maxv,nend,vlamda,g,maxv,scrap,ifail)
      do 371 i=1,km
 371  sout(i)=vlamda(i)+core
      if(oprint(3))call writem(g,ilifg,nend,km)
      len=nend-nstart+1
c
c.....convergence test
c
      do 380 i=1,km
      r(i)=dnrm2(len,g(nstart,i),1)
 380  continue
      tsec=cpulft(1)
      test=dzero
      do 400 i=1,km
400   if (test.lt.r(i))test=r(i)
c
cjvl    if what is added to us satures the complete space, we are 
cjvl    converged, no matter what the program might think
cjvl
        if (nend.eq.kk3) test = 0.0d0
c
c
      if(nprint.ne.-5)write(iwr,390)mtx,tsec,test,(sout(i),i=1,km)
 390  format(i6,f8.2,f10.6,5f20.7)
c
c.....switch to density matrix mode if converged or too many iterations
      if(test.lt.acc.or.mtx.ge.nits)goto550
c
c....now form the m correction vectors
      ibase = nend
      if(nend+km.gt.maxvec) goto 510
      call vclr(buffer(1,nend+1),1,kk3*km)
      do 420 k=1,km
      vlk=vlamda(k)
      if(mtx.eq.1)vlk=vlk+.0000001d0
      do 420 j=1,nend
      gjk=g(j,k)
      do 420 i=1,kk3
420   buffer(i,k+nend)=buffer(i,k+nend)+gjk*(buffer(i,j+maxvec)
     +   - vlk*buffer(i,j) ) / (vlk-buffer(i,ibd) )
c
c....normalise the correction vectors
c
      do 430 k=1,km
      s=1.0d0/dnrm2(kk3,buffer(1,k+nend),1)
      call dscal(kk3, s, buffer(1,k+nend),1)
 430  continue
c
c
c....now take each correction vector in turn; orthogonalise it to
c    all previous expansion functions. if its norm is then less than
c     a certain threshold, reject it. otherwise, normalise and add to th
c    set.
c
      do 460 k=1,km
      do 440 l=1,nend
      buff=ddot(kk3,buffer(1,l),1,buffer(1,ibase+k),1)
      call daxpy(kk3,-buff
     +  ,buffer(1,l),1,buffer(1,ibase+k),1)
 440  continue
      s=dnrm2(kk3,buffer(1,ibase+k),1)
      if(s.lt.thresh.and.k.ne.1)goto460
      nend=nend+1
      s=1.0d0/s
_IF1(civ)      call scaler(kk3,s,buffer(1,nend),buffer(1,ibase+k))
_IFN1(civ)      call vsmul(buffer(1,ibase+k),1,s,buffer(1,nend),1,kk3)
c
c
460   continue
c
c.....refresh the orthonormality of the new set of expansion functions
c     because of round-off
c
      nstart=ibase+1
c     nvec=nend-ibase
      do 490 k=nstart,nend
      km1=k-1
      do 470 l=1,km1
      buff=ddot(kk3,buffer(1,l),1,buffer(1,k),1)
      call daxpy(kk3,-buff
     +   ,buffer(1,l),1,buffer(1,k),1)
c
 470  continue
c
      s=1.0d0/dnrm2(kk3,buffer(1,k),1)
      call dscal(kk3,s,buffer(1,k),1)
c
490   if(oprint(3).and.nprint.ne.-5)write(iwr,500)k,
     *                                       (buffer(i,k),i=1,kk3)
500   format(' new expansion function number',i3/100(10f13.6/))
c
c     end of iteration; now return and calculate sigma for the new funct
c
      goto220
c
c.....no space left
c.....collect up the solution so far
c
 510  ilen=km*kk3
      call vclr(buffer(1,maxvec+1),1,ilen)
      call mxmb(buffer,1,kk3,g,1,maxv,buffer(1,maxvec+1),1,kk3,
     * kk3,ibase,km)
c
      call dcopy(ilen,buffer(1,maxvec+1),1,buffer,1)
      nstart=1
      nend=km
c     nvec=km
      if(nprint.ne.-5)write(iwr,540)
540   format(' insufficient store available for further expansion vector
     1s - diagonalisation restarting with solution so far')
      goto 220
c
c
c
c
c
c.....end of diagonalisation
c.....construct final vectors and then get density matrices
c
 550  ilen=km*kk3
      call vclr(buffer(1,maxvec+1),1,ilen)
      call mxmb(buffer,1,kk3,g,1,maxv,buffer(1,maxvec+1),1,kk3,
     * kk3,nend,km)
c550  do 5580 kkk=1,km
c     do 5580 ll=1,kk3
c     buff=zeroe
c     do 5570 k=1,nend
c5570  buff=buff+buffer(ll,k)*g(k,kkk)
c5580  buffer(ll,maxvec+kkk)=buff
      do 600 kkk=1,km
      call dcopy(kk3,buffer(1,maxvec+kkk),1,ccc(1,kkk),1)
      if(oprint(3))write(iwr,590)kkk,(ccc(ll,kkk),ll=1,kk3)
590   format(' eigenvector number',i2/100(10f13.6/))
600   vl(kkk)=vlamda(kkk)
      kkk=km
      goto 640
c
c***********************************************************************
c
c.....generation of density matrices after full newton-raphson
c
610   call search(iblk,ifillp)
      call find(ifillp)
      call dcopy(kk3,buffer(1,km),1,ccc(1,km),1)
c
c.....construct density matrices
640   n21e=n2e+n11e
      call vclr(twoe,1,n21e)
c
c
c*********************************
c.....loop formula tape
650   call get(vsurd,nw)
      if(nw.eq.0)goto 770
      call find(ifillp)
_IF1(iv)      call upakin(iy,nuwks,511)
_IFN1(iv)      call unpack(iy,16,nuwks,2044)
      if(nword.eq.0)goto650
c.....on-diagonal
      if(icode.eq.1)goto 690
      do 680 iword=1,nword
      nlwk=nlwks(iword)
      nuwk=nuwks(iword)
      iuwk=iuwks(iword)
      val =surd(icoup(iword))
      gam=dzero
      do 670 j=1,nlwk
      ii=indx(iuwk)
      do 660 i=1,nuwk
      gam=gam+ccc(ii,km)*ccc(ii,km)*val
660   ii=ii+1
670   iuwk=iuwk+1
680   twoe(ints(iword))=twoe(ints(iword))+gam
      goto 650
c
c......off-diagonal
690   do 720 iword=1,nword
      nlwk=nlwks(iword)
      nuwk=nuwks(iword)
      iuwk=iuwks(iword)
      juwk=juwks(iword)
      val =surd(icoup(iword))
      val=val+val
      gam=dzero
      do 710 i=1,nlwk
      ii=indx(iuwk)
      jj=indx(juwk)
      do 700 j=1,nuwk
      gam=gam+ccc(jj,km)*ccc(ii,km)*val
      ii=ii+1
700   jj=jj+1
      iuwk=iuwk+1
710   juwk=juwk+1
720   twoe(ints(iword))=gam+twoe(ints(iword))
      goto 650
c
c
c.....now separate the 1-electron bit
c
 770  continue
      call vclr(gam1,1,kk9)
      call dcopy(n11e,twoe(n2e+1),1,gam1(ic1e(nst,1)+1),1)
c
c
      if(odmgen)goto 810
c
c***********************************************************************
c
c
c
      if(nprint.ne.-5)write(iwr,780)
 780  format(/1x,64('*')/
     *' state number',10x,'energy',10x,'overlap with trial function'/
     *1x,64('*'))
      do 790 kkk=1,km
      s=ddot(kk3,buffer(1,kkk),1,buffer(1,maxvec+kkk),1)
      ve=vl(kkk)+core
790   if(nprint.ne.-5)write(iwr,800)kkk,ve,s
800   format(i8,f25.11,f25.6)
810   call search(iblkci,idaf)
      do 820 kkk=1,km
      call wrt3s(ccc(1,kkk),kk3,idaf)
      call wrt3s(gam1,kk9,idaf)
820   call wrt3s(twoe,kk8,idaf)
c
c
      do 840 ml=1,nprim
      m1=ic1e(ml,1)+ic1e(ml,2)
      pop(ml)=gam1(m1)
840   continue
c
      if(oprint(1))call fout(skv,gam1,m3,iwr)
c
      return
      end
      subroutine super(q,dsx,work,ic5e,ic6e,skv,brill,scc,
     * fock,den1,den2,ffff,iwr,nprint)
c
c.....this routine does the super-ci stage in the casscf program
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c for symmetry assignment
      common/fsymas/dddum,otsym
      common/csymas/popp2(maxorb)
c
      common/craypk /eps(maxorb),iyc1(maxorb),scrtch(maxorb)
INCLUDE(common/avstat)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/popnos)
      common/l     /lock
INCLUDE(common/exc)
INCLUDE(common/gjs)
INCLUDE(common/ciconv)
INCLUDE(common/mapper)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1             nst,n1e,isym,mults
      common/potn/core
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/blbpar/acc,vlev,time,nit,nprinb
INCLUDE(common/prints)
      common/tab3  /is1,is1p,is2,is2p
INCLUDE(common/dm)
INCLUDE(common/finish)
      common/junk/popp(maxorb),map(15),map1(15)
c
      character*10 charwall
      dimension q(*),skv(*),dsx(length,*),ffff(*)
      dimension work(length,*)
     2,ic5e(*),ic6e(*),scc(*),brill(*)
      dimension fock(*),den1(*),den2(*)
c
c
      tst=cpulft(1)
      if(nprint.ne.-5)write(iwr,430)tst ,charwall()
 430  format(/1x,37('*')/1x,
     * 'commence super-ci at ',f8.2,' seconds',a10,' wall'/1x,37('*')/)
      sq2= dsqrt(2.0d0)
      ss=1.0d0
      nits=nit
      mta=0
c
      do 20 i=nst,nprim
      if(pop(i).gt.1.999999999999d0)goto30
      popsqm(i)= dsqrt(2.0d0-pop(i))
      go to 10
 30   popsqm(i)=1.0d-6
 10   if(pop(i).lt.1.0d-09)goto40
      popsqn(i)= dsqrt(pop(i))
      go to 20
 40   popsqn(i)=1.0d-5
 20   continue
c
c.....scan through integral file to pick up brillouin matrix els.
c
      call vclr(skv(n1e+1),1,n1e)
      call getblb(dsx(kk10+1,1),kk10,q,ic5e,ic6e,skv,dsx,
     + fock,den1,den2,ffff,iwr)
      call dcopy(kk10,dsx,1,brill,1)
      call dpconv(nprim,den1,den2)
      call vclr(scc(2),1,kk10)
      scc(1)=1.0d0
c
c......get starting energy for super-ci
c
      ve=vl(kkkm)
      if(nwt.eq.1)goto 80
      ve=ddot(nwt,weight,1,vl,1)
80    ve0=ve
c
c
_IF1(v)      vmax=0.0d0
_IF1(v)      do 9 kk=2,kk10
_IF1(v)      if(vmax.lt. dabs(brill(kk)))vmax= dabs(brill(kk))
_IF1(v)    9 continue
c
      kk=idamax(kk10-1,brill(2),1)
      vmax=dabs(brill(kk+1))
c.....enter first cycle of super-ci diagonalisation..
100   mta=mta+1
      if(mta.eq.nits+1)go to 130
c.....scan each row constructing zs (sigma),then use it to
c     increment coefficient,energy and overlap
c
      test=-10.0d0
      do 110 kk=2,kk10
      call matrix(vlev,vs,vd,kk,ve0,ic5e,ic6e,skv,brill,scc)
      vb=vd-ve
      va=(vd*scc(kk)-vs)/ss
      vs=vs-ve*scc(kk)
      dc= dsqrt(vb*vb/4.0d0-va*vs)
      if(vb.le.0)dc=(dc-vb/2.0d0)/va
      if(vb.gt.0)dc=vs/(-dc-vb/2.0d0)
      ss=ss+(2.0d0*scc(kk)+dc)*dc
      scc(kk)=scc(kk)+dc
      ve=ve+(2.0d0*vs+dc*vb)*dc/ss
c
      dc= dabs(dc)
      if(dc.gt.test)test=dc
110   continue
c
      if(test.le.acc)nits=mta
      go to 100
c.....end of iterative diagonalisation
 130   t=cpulft(1)
       if(nprint.ne.-5)
     * write(iwr,120)mta,t,charwall(),test
 120   format(/1x,
     *'iterative diagonalization converged at cycle ',i3,
     *'  after ',f8.2,' seconds',a10,' wall'
     */1x,'tester = ',f12.7)
      ss=1.0d0/ dsqrt(ss)
      call dscal(kk10,ss,scc,1)
      smfg=ve+core
c
      do 180 kk=2,kk10
      if(ic5e(kk).gt.ncore)go to 160
      if(ic6e(kk).le.nprim)scc(kk)=scc(kk)/popsqm(ic6e(kk))
      if(ic6e(kk).gt.nprim)scc(kk)=scc(kk)/sq2
      go to 170
160   scc(kk)=scc(kk)/popsqn(ic5e(kk))
170   test= dabs(scc(kk))
180   continue
      fmax=vmax
      if(nprint.ne.-5)write(iwr,200)smfg,vmax
      if(oprint(8).and.nprint.ne.-5)write(iwr,190)(scc(ll),ll=1,kk10)
190   format(/40x,'super-ci coefficients'/
     1/100(10f13.8/)/)
 200  format(/1x,
     *'final eigenvalue = ',f19.11,5x,
     *'maximum brillouin matrix element  =',f12.8/)
c
c...  avoid superci convergence if dynamic shifting is on
c
      oconv=vmax.lt.cccnv
      if (osuped.and.swnr.ne.0.0d0) oconv = .false.
c
      call vclr(popp,1,ma)
c
c.....optimise orbitals for each symmetry type in turn
c
      do 410 iset=1,ntype
      nthis=nsymm(iset)
      if(nthis.eq.0)go to 410
      if(oprint(4))write(iwr,210)iset
210   format(/' symmetry',i2)
c
c.....set up density matrix dsx
c
      call dens(dsx,work,iset,nthis,ic5e,ic6e,scc,den1,den2)
c
c.....print density matrix....
c
      mark=0
      mmmm=nthis
      if(oprint(5))call dout(mmmm,mark,dsx,length,iwr)
c.....diagonalise density matrix..
_IF()
cjvl  replaced givens by a jacobi for max overlap with prev vectors
      ifail=0
      call f02abf(dsx,length,nthis,scrtch,dsx,length,eps,ifail)
_ENDIF
c
      call f02rep(dsx,length,nthis,work,eps,iky,scrtch)
      call dcopy(length,eps,1,scrtch,1)
c
c.....sort eigenvectors....
c
      mark=2
      if(oprint(5))call dout(mmmm,mark,dsx,length,iwr)
      call setsto(nthis,-100,iyc1)
c
c ..select locked orbitals first (either core or all primary)
c ..lock off yields lock = 0 here which is caused by lock=-1 in master
c ..( no max overlap sorting then)
c
      icount=0
      if(lock.eq.0)go to 260
      do 250 ii=1,lock
      if(isymmo(ii).ne.iset)go to 250
      icount=icount+1
      test=-1000.0d0
      do 230 j=1,nthis
      if(iyc1(j).gt.0)goto230
      test1= dabs(dsx(icount,j))
      if(test1.le.test)go to 230
      test=test1
      mmax=j
230   continue
      call dcopy(nthis,dsx(1,mmax),1,work(1,icount),1)
      iyc1(mmax)=10
      eps(icount)=scrtch(mmax)
250   continue
260   icount=icount+1
      itop=0
c
c.....loop over final vectors..
      do 290 i=icount,nthis
      test=-10000.0d0
      if(nblock(iset,i).le.nprim)itop=i
c
c.....loop over initial vectors
      do 270 j=1,nthis
      if(iyc1(j).gt.0)go to 270
      if(scrtch(j).le.test)go to 270
      test=scrtch(j)
      mmax=j
270   continue
c
      call dcopy(nthis,dsx(1,mmax),1,work(1,i),1)
      iyc1(mmax)=10
290   eps(i)=scrtch(mmax)
      if(itop.eq.0.or.lock.eq.0)go to 340
      if(oprint(5))call dout(mmmm,mark,work,length,iwr)
c
c.....re-sort in active space to restore original configs.
      do 300 i=icount,itop
      scrtch(i)=eps(i)
      call dcopy(nthis,work(1,i),1,dsx(1,i),1)
 300  continue
c
      do 310 i=icount,itop
      ic=nblock(iset,i)-ncore
310   map(i)=iocca(kkkm,ic)+ioccb(kkkm,ic)
c
      ic=icount-1
      do 320 i=1,3
      ii=3-i
      do 320 j=icount,itop
      if(map(j).ne.ii)go to 320
      ic=ic+1
      map1(j)=ic
320   continue
c
      do 330 i=icount,itop
      mapp=map1(i)
      eps(mapp)=scrtch(i)
      call dcopy(nthis,dsx(1,i),1,work(1,mapp),1)
330   continue
c
340   continue
c.....check phases of transform. matrix
      do 380 i=1,nthis
      test=0.0d0
      do 360 j=1,nthis
      iphase=1
      test1=work(j,i)
      if(test1.gt.0.0d0)go to 350
      test1=-test1
      iphase=-iphase
350   if(test1.le.test)go to 360
      test=test1
      k=iphase
360   continue
      if(k.gt.0)go to 380
      do 370 j=1,nthis
370   work(j,i)=-work(j,i)
380   continue
c
c.....write out eigenvectors of density matrix
      mark=2
      mmmm=nthis
      if(oprint(4)) call dout(mmmm,mark,work,length,iwr)
c.....write out eigenvalues of density matrix..
      if(oprint(4))write(iwr,390)
      if(oprint(4))write(iwr,400)(eps(i),i=1,nthis)
390   format('+',40x,'eigenvalues of density matrix'/)
400   format(6x,8f14.7)
      mark=0
      do 411 i=1,ma
      if(isymmo(i).ne.iset)go to 411
      mark=mark+1
      popp(i)=eps(mark)
 411  continue
      if(.not.oconv)
     1           call trnsfo(q,dsx,work,iset,length)
410   continue
c
      if(otsym) then
       do 1001,i=1,ma
1001   popp2(i)=popp(i)
      endif
c
      if(nprint.ne.-5)write(iwr,412)
 412  format(/40x,
     *'occupation numbers'/40x,18('-'))
      if(nprint.ne.-5)
     *write(iwr,413)(popp(i),i=1,ma)
 413  format(/10x,8f14.7)
      mark=4
      mmmm=ma
      if(oprint(2)) call dout(mmmm,mark,q,ma,iwr)
      t=cpulft(1)
      if(nprint.ne.-5)write(iwr,420)t ,charwall()
420   format(/1x,'super-ci procedure complete at ',f8.2,
     1' seconds',a10,' wall'/)
      return
      end
      subroutine f02rep(dsx,length,nthis,work,eps,iky,ilifq)
c
      implicit REAL (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...   a f02abf (symmetric diagonaliser replacement) in case 
c...   this givens routine is not appreciated 
c...   e.g. IT SCREWS UP THE SUPER CI in CASSCF (jvl 1995)
c...        because it mixes degenerate natural orbitals
c...   dsx is dimensioned dsx(length,nthis) (length > nthis we hope)
c...   called from super
c
       dimension dsx(length,nthis),work(*),eps(nthis),iky(*)
       dimension ilifq(nthis)
c
c...   bring matrix to triangular form
c
      kk = 0
      do 10 i=1,nthis
        do 10 j=1,i
          kk = kk + 1
          work(kk) = dsx(i,j)
10    continue
c
      do 20 i=1,nthis
         ilifq(i) = (i-1)*length
20    continue
c
c...  jacobi call parameter one before last determines eigenvalue order
c...  0 nop / 1 unsorted / 2 increasing / 3 decreasing / 4 lock
c
      call jacobi(work,iky,nthis,dsx,ilifq,nthis,eps,2,4,1.0d-12)
c
      return
      end
      subroutine dens(dsx,iwork,iset,nthis,ic5e,ic6e,scc,gam1,gam2)
c.....this routine sets the values of the (ip,iq) elements of
c      the super-ci density matrix. it is returned in dsx.
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
INCLUDE(common/popnos)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1           nst,n1e,isym,mults
      common/dims/kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/tab3  /is1,is1p,is2,is2p
      dimension scc(*),gam1(*),gam2(*)
      dimension dsx(length,*),iwork(ma,*),ic5e(*),ic6e(*)
c
      scc(kk10+1)=0.0d0
      len2=ma*ma
      call setsto(len2,kk10+1,iwork)
c
      len2=length*length
      call vclr(dsx,1,len2)
c
c     vmax=0.0d0
c     sq2= dsqrt(2.0d0)
c
      do 30 kk=2,kk10
30    iwork(ic5e(kk),ic6e(kk))=kk
c
      onocor=ncore.eq.0
      do 290 mmm=1,nthis
      ip=nblock(iset,mmm)
      do 290 nnn=1,mmm
      iq=nblock(iset,nnn)
      dsx(mmm,nnn)=0.0d0
      if(ip.le.ncore)go to 60
      if(ip.gt.nprim)go to 40
      if(iq-ncore)90,90,130
40    if(iq-ncore)120,120,50
50    if(iq-nprim)210,210,240
c
c---------------------------------------------------------------------
c.....core-core element
c     -----------------
60    if(ip.eq.iq)dsx(mmm,nnn)=2.0d0
      do 70 ia=nprimp,ma
70    dsx(mmm,nnn)=dsx(mmm,nnn)-scc(iwork(ip,ia))*scc(iwork(iq,ia))*2
c
      do 80 it=nst,nprim
      do 80 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(it))go to 80
      m1=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      temp=0.0d0
      if(it.eq.iu)temp=2.0d0
      dsx(mmm,nnn)=dsx(mmm,nnn)-scc(iwork(ip,it))*scc(iwork(iq,iu))*
     1(temp-gam1(m1))
80    continue
      dsx(nnn,mmm)=dsx(mmm,nnn)
      go to 280
c-----------------------------------------------------------------------
c.....core-as element
c     ---------------
90    do 110 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(ip))go to 110
      m1=ic1e(max(ip,iu),1)+ic1e(min(ip,iu),2)
      temp=-gam1(m1)
      if(ip.eq.iu)temp=temp+2.0d0
      temp=temp*scc(1)*scc(iwork(iq,iu))
      temp1=0.0d0
      do 100 ia=nprimp,ma
100   temp1=temp1+scc(iwork(iq,ia))*scc(iwork(iu,ia))
      dsx(mmm,nnn)=dsx(mmm,nnn)+temp-temp1*gam1(m1)
110   continue
      dsx(nnn,mmm)=dsx(mmm,nnn)
      go to 280
c-----------------------------------------------------------------------
c.....core-secondary element
c     ----------------------
120   dsx(mmm,nnn)=scc(1)*scc(iwork(iq,ip))*2.0d0
      dsx(nnn,mmm)=dsx(mmm,nnn)
      go to 280
c-----------------------------------------------------------------------
c.....as-as element
c     -------------
130   m1=ic1e(ip,1)+ic1e(iq,2)
      n1=iky(ip-ncore)+iq-ncore
      dsx(mmm,nnn)=gam1(m1)
      do 200 iv=nst,nprim
      do 170 ix=nst,nprim
      if(isymmo(ix).ne.isymmo(iv))go to 170
      n2=iky(max(iv,ix)-ncore)+min(iv,ix)-ncore
      m2=ic1e(max(iv,ix),1)+ic1e(min(iv,ix),2)
      n2=ic2e(max(n1,n2))+ic3e(min(n1,n2))
      temp1=2.0d0*gam2(n2)-gam1(m1)*gam1(m2)
      temp=0.0d0
      do 140 ia=nprimp,ma
140   temp=temp+scc(iwork(iv,ia))*scc(iwork(ix,ia))
      if(onocor)go to 160
      do 150 ii=1,ncore
150   temp=temp-scc(iwork(ii,iv))*scc(iwork(ii,ix))
160   dsx(mmm,nnn)=dsx(mmm,nnn)+temp*temp1
170   continue
c
c
      if(onocor)go to 200
      if(isymmo(iv).ne.isymmo(iq))go to 200
      temp=0.0d0
      do 180 ii=1,ncore
180   temp=temp+scc(iwork(ii,iv))*scc(iwork(ii,iq))
      m3=ic1e(max(ip,iv),1)+ic1e(min(ip,iv),2)
      temp1=-gam1(m3)
      if(ip.eq.iv)temp1=temp1+1.0d0
      dsx(mmm,nnn)=dsx(mmm,nnn)+temp1*temp
      temp=0.0d0
      do 190 ii=1,ncore
190   temp=temp+scc(iwork(ii,iv))*scc(iwork(ii,ip))
      m3=ic1e(max(iq,iv),1)+ic1e(min(iq,iv),2)
      temp1=-gam1(m3)
      if(iq.eq.iv)temp1=temp1+1.0d0
      dsx(mmm,nnn)=dsx(mmm,nnn)+temp1*temp
200   continue
      dsx(nnn,mmm)=dsx(mmm,nnn)
      go to 280
c---------------------------------------------------------------------
c.....as-secondary elements
c     ---------------------
c
210   do 230 iu=nst,nprim
      if(isymmo(iq).ne.isymmo(iu))go to 230
      m1=ic1e(max(iu,iq),1)+ic1e(min(iu,iq),2)
      dsx(mmm,nnn)=scc(1)*scc(iwork(iu,ip))*gam1(m1)+dsx(mmm,nnn)
      if(onocor)go to 230
      do 220 ii=1,ncore
      temp=-gam1(m1)
      if(iq.eq.iu)temp=temp+2.0d0
      dsx(mmm,nnn)=dsx(mmm,nnn)+temp*scc(iwork(ii,ip))*scc(iwork(ii,iu))
220   continue
230   continue
      dsx(nnn,mmm)=dsx(mmm,nnn)
      go to 280
c-----------------------------------------------------------------------
c.....secondary-secondary elements
c     ----------------------------
c
240   if(onocor)go to 260
      do 250 ii=1,ncore
250   dsx(mmm,nnn)=dsx(mmm,nnn)+scc(iwork(ii,ip))*scc(iwork(ii,iq))
      dsx(mmm,nnn)=dsx(mmm,nnn)*2.0d0
c
260   do 270 it=nst,nprim
      do 270 iu=nst,nprim
      if(isymmo(it).ne.isymmo(iu))go to 270
      m1=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      dsx(mmm,nnn)=dsx(mmm,nnn)+scc(iwork(it,ip))*scc(iwork(iu,iq))*
     1gam1(m1)
270   continue
      dsx(nnn,mmm)=dsx(mmm,nnn)
280   if(mmm.eq.nnn)go to 290
      dsx(nnn,mmm)=dsx(mmm,nnn)
290   continue
       return
      end
      subroutine getblb(work,kk10,q,ic5e,ic6e,skv,brill,qq,
     + gam1,gam2,ffff,iwr)
c
c.....scans the file of integrals with 1  index outside the active
c     space, and puts their contributions into
c     the brillouin matrix elements
c
c.....requires gamma-type density matrices
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/craypk/i205(1360)
INCLUDE(common/restar)
      common/cone  /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),ims(8)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/popnos)
INCLUDE(common/prints)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
INCLUDE(common/filel)
      common/disc  /isel,iselr,iselw,irep,ichek,ipos(maxlfn)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults
      common/blkin /g(510),nword
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
INCLUDE(common/gjs)
      dimension id(4),brill(*),work(ma,*),q(ma,*),ic5e(*),ic6e(*),
     +          skv(*),gam1(*),gam2(*),qq(*),ffff(*)
      data half/0.5d0/, two/2.0d0 /
c
      lentri=ma*(ma+1)/2
      len2=lentri+lentri
      call rdfock(q,ffff,ffff(len2+1),skv,gam1)
      sq2=dsqrt(2.0d0)
c
      kk=nprim*nprim
      call vclr(qq,1,kk)
      kk101=kk10-1
      call vclr(brill(2),1,kk101)
      do 10 kk=2,kk10
10    work(ic5e(kk),ic6e(kk))=kk+.1d0
c
      lfile=m6file
      do 180 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i)=m6blk(i)
 180  llblk(i)=liblk(i)-m6last(i)
c
      do 190 ii=1,lfile
      iunit=lotape(ii)
      call search(liblk(ii),iunit)
      call find(iunit)
      jj=llblk(ii)
20    jj=jj+1
      call get(g(1),nw)
      if(nw)190,190,30
30    if(jj.ne.0)call find(iunit)
_IF1(iv)      call upak8v(g(num2e+1),i205)
_IFN1(iv)      call unpack(g(num2e+1),lab816,i205,numlab)
_IFN1(iv)      int4=1
      do 3 kk=1,nword
_IF(ibm,vax)
      i=i205(kk)
      j=j205(kk)
      k=k205(kk)
      l=l205(kk)
_ELSEIF(littleendian)
      j=i205(int4  )
      i=i205(int4+1)
      l=i205(int4+2)
      k=i205(int4+3)
_ELSE
      i=i205(int4  )
      j=i205(int4+1)
      k=i205(int4+2)
      l=i205(int4+3)
_ENDIF
      ggg=g(kk)
      do 40 iii=1,4
40    id(iii)=0
c
      if(i.le.ncore.or.i.gt.nprim)id(1)=1
      if(j.le.ncore.or.j.gt.nprim)id(2)=1
      if(k.le.ncore.or.k.gt.nprim)id(3)=1
      if(l.le.ncore.or.l.gt.nprim)id(4)=1
c
      idtot=1
      do 50 iii=1,4
50    idtot=idtot+id(iii)
c
      go to (120,60,3,3,3),idtot
c
60    if(id(1).eq.1)go to 70
      if(id(2))100,100,90
c
c----------------------------------------------------------------------
c.....(au|vx) integrals
c
70    n2=iky(k-ncore)+l-ncore
c
      do 80 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))go to 80
      ll=work(it,i)
      n1=iky(max(it,j)-ncore)+min(it,j)-ncore
      n3=ic2e(max(n1,n2))+ic3e(min(n1,n2))
      temp=ggg*gam2(n3)
      if(j.ne.it)temp=temp*0.5d0
      if(n1.ne.n2)temp=temp*0.5d0
      brill(ll)=brill(ll)+temp
80    continue
      go to 3
c
c-----------------------------------------------------
c.....(vi|tu)---> (tu|vi)
c
90    k1=k
      l1=l
      k=i
      l=j
      i=k1
      j=l1
c
c.....all (tu|vi) integrals now
c
100   n1=iky(i-ncore)+j-ncore
      do 110 it=nst,nprim
      if(isymmo(it).ne.isymmo(l))go to 110
      ll=work(l,it)
      n2=iky(max(it,k)-ncore)+min(it,k)-ncore
      n3=ic2e(max(n1,n2))+ic3e(min(n1,n2))
      temp=ggg*gam2(n3)
      if(k.ne.it)temp=temp*0.5d0
      if(n1.ne.n2)temp=temp*0.5d0
      brill(ll)=brill(ll)-temp
110   continue
      goto 3
c.....(tu|vx)
120   iuv=iky(i-ncore)+j-ncore
      ixy=iky(k-ncore)+l-ncore
      do 160 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto130
      in1=nprim*(it-1)+i
      in2=iky(max(it,j)-ncore)+min(it,j)-ncore
      in3=ic2e(max(ixy,in2))+ic3e(min(ixy,in2))
      dumx=gam2(in3)*ggg
      if(it.ne.j)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
130   if(isymmo(it).ne.isymmo(j).or.i.eq.j)goto140
      in1=nprim*(it-1)+j
      in2=iky(max(it,i)-ncore)+min(it,i)-ncore
      in3=ic2e(max(in2,ixy))+ic3e(min(in2,ixy))
      dumx=gam2(in3)*ggg
      if(it.ne.i)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
140   if(iuv.eq.ixy)goto160
      if(isymmo(it).ne.isymmo(k))goto150
      in1=nprim*(it-1)+k
      in2=iky(max(it,l)-ncore)+min(it,l)-ncore
      in3=ic2e(max(in2,iuv))+ic3e(min(in2,iuv))
      dumx=gam2(in3)*ggg
      if(it.ne.l)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
150   if(isymmo(it).ne.isymmo(l).or.k.eq.l)goto160
      in1=nprim*(it-1)+l
      in2=iky(max(it,k)-ncore)+min(it,k)-ncore
      in3=ic2e(max(in2,iuv))+ic3e(min(in2,iuv))
      dumx=gam2(in3)*ggg
      if(it.ne.k)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
160   continue
_IF1(iv)    3 continue
_IFN1(iv)    3 int4=int4+4
      if(jj)20,190,20
190   continue
      call dscal(kk101,two,brill(2),1)
c------------------------------------------------------------
c.....now finish evaluating brill. elements
c------------------------------------------------------------
c.....<0|h|i-a> elements
c
      if(ncore.eq.0)go to 230
      do 220 i=1,ncore
      do 210 ia=nprimp,ma
      if(isymmo(i).ne.isymmo(ia))go to 210
      m1=ic1e(ia,1)+ic1e(i,2)
      kk=work(i,ia)
      brill(kk)=sq2*(skv(n1e+m1)+skv(m1))
210   continue
220   continue
c
c.....<0|h|t-a> elements
c
230   do 250 it=nst,nprim
      do 250 ia=nprimp,ma
      if(isymmo(it).ne.isymmo(ia))go to 250
      kk=work(it,ia)
      do 240 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(it)) go to 240
      m1=ic1e(ia,1)+ic1e(iu,2)
      m2=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      temp=skv(m1)*gam1(m2)
      if(it.ne.iu)temp=temp*0.5d0
      brill(kk)=brill(kk)+temp
240   continue
      brill(kk)=brill(kk)/popsqn(it)
250   continue
c
c.....<0|h|i-t> elements
c
      if(ncore.eq.0)go to 280
      do 270 i=1,ncore
      do 270 it=nst,nprim
      if(isymmo(i).ne.isymmo(it))go to 270
      kk=work(i,it)
      m1=ic1e(it,1)+ic1e(i,2)
      brill(kk)=brill(kk)+2.0d0*skv(n1e+m1)
      do 260 iu=nst,nprim
      if(isymmo(it).ne.isymmo(iu)) go to 260
      m1=ic1e(iu,1)+ic1e(i,2)
      m2=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      if(it.eq.iu)brill(kk)=brill(kk)+skv(m1)*(2.0d0-gam1(m2))
      if(it.ne.iu)brill(kk)=brill(kk)-skv(m1)*gam1(m2)*0.5d0
260   continue
      brill(kk)=brill(kk)/popsqm(it)
270   continue
280    continue
c.....fock operator
      do 300 i=nst,nprim
      i1=ilifp(i)
      isymmoi=isymmo(i)
      do 300 j=nst,nprim
      if(isymmoi.ne.isymmo(j)) goto 300
      buff=0.0d0
      do 290 k=nst,i
      if(isymmoi.eq.isymmo(k)) buff=buff +
     1 skv(ic1e(max(j,k),1)+ic1e(min(j,k),2)) *
     2 gam1(ic1e(i,1)+ic1e(k,2) )
290   continue
      do 291 k=i,nprim
      if(isymmoi.eq.isymmo(k)) buff=buff+
     + skv(ic1e(max(j,k),1)+ic1e(min(j,k),2) )*
     + gam1(ic1e(k,1)+ic1e(i,2) )
291   continue
      qq(i1+j)=qq(i1+j)+buff*0.25d0
300   continue
      if(.not.oprint(12))goto330
      if(nprint.eq.-5)go to 330
      write(iwr,310)
310   format(/50x,'brillouin elements'/50x,18('*')/)
      write(iwr,320)(brill(kk),kk=2,kk10)
320   format(1x,10f13.6)
330   call clredx
      return
      end
_EXTRACT(matrix,_AND(hp800,itanium))
      subroutine matrix(vlev,vs,vd,k2,ve0,ic5e,ic6e,skv,brill,scc)
c
c.....this routine works out zs (sigma vector) for row kk, and also
c     zd, the diagonal element <kk|h|kk>
c     in the super-ci matrix.
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gjs)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/popnos)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/tab3  /is1,is1p,is2,is2p
      dimension skv(*),ic5e(*),ic6e(*)
      dimension scc(*)
      dimension brill(*)
c
      sq2= dsqrt(2.0d0)
      kk=k2
      vs=0.0d0
c
c.....loop over all ll
c
c
      do 120 ll=1,kk10
c
      k1=max(kk,ll)
      l1=min(kk,ll)
      ip=ic5e(k1)
      iq=ic6e(k1)
      ir=ic5e(l1)
      is=ic6e(l1)
c
      sxhsx=0.0d0
c.....transfer control to the relevant formula
c
c
      if(ir.eq.0)go to 10
      if(ip.le.ncore)go to 20
      if(ir.gt.ncore)go to 100
      if(is-nprim)80,80,90
c
c-----------------------------------------------------------------
c
c.....<i-t|h|0> elements
c     ------------------
c
10    if(ip.gt.ncore)go to 70
      if(iq.gt.nprim)go to 40
      sxhsx=brill(k1)
      go to 110
c---------------------------------------------------------------------
c
c.....<i-t|h|j-u> elements
c     --------------------
20    if(is.gt.nprim)go to 60
      if(iq.gt.nprim)go to 50
      ii=ip
      it=iq
      ij=ir
      iu=is
      m1a=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      m1b=ic1e(max(ii,ij),1)+ic1e(min(ii,ij),2)
      temp=0
      if(ii.ne.ij)go to 30
      if(it.ne.iu)temp=(2.0d0-pop(it))*(2.0d0-pop(iu))*
     1(skv(m1a)+skv(n1e+m1a))/2.0d0
      if(it.eq.iu)temp=(2.0d0-pop(it))*(skv(m1a)+skv(n1e+m1a))
30    if(it.eq.iu)temp=temp-(2.0d0-pop(it))*(skv(m1b)+skv(n1e+m1b))
      sxhsx=temp/(popsqm(it)*popsqm(iu))
      if(kk.ne.ll)go to 110
      sxhsx=sxhsx+vlev+ve0
      vd=sxhsx
      go to 110
c---------------------------------------------------------------------
c
c.....<i-a|h|0> element
c     -----------------
40    ii=ip
      ia=iq
      if(isymmo(ia).ne.isymmo(ii)) go to 110
      sxhsx=brill(k1)
      go to 110
c----------------------------------------------------------------------
c
c.....<i-a|h|j-t> element
c     -------------------
50    ii=ip
      ia=iq
      ij=ir
      it=is
      if(ii.ne.ij.or.isymmo(ia).ne.isymmo(it))go to 110
      m1=ic1e(ia,1)+ic1e(it,2)
      sxhsx=(popsqm(it)/sq2)*(skv(m1)+skv(n1e+m1))
      go to 110
c-----------------------------------------------------------------------
c
c.....<i-a|h|j-b> element
c     -------------------
60    ii=ip
      ia=iq
      ij=ir
      ib=is
      if(isymmo(ia).ne.isymmo(ib).or.isymmo(ii).ne.isymmo(ij))
     +   go to 110
      m1=ic1e(max(ii,ij),1)+ic1e(min(ii,ij),2)
      m1a=ic1e(max(ia,ib),1)+ic1e(min(ia,ib),2)
      if(ii.eq.ij)sxhsx=skv(m1a)+skv(n1e+m1a)
      if(ia.eq.ib)sxhsx=sxhsx-skv(m1)-skv(n1e+m1)
      if(kk.ne.ll)go to 110
      sxhsx=sxhsx+vlev+ve0
      vd=sxhsx
      go to 110
c----------------------------------------------------------------------
c
c.....<t-a|h|0> element
c     -----------------
70    it=ip
      ia=iq
      if(isymmo(ia).ne.isymmo(it))go to 110
      sxhsx=brill(k1)
      go to 110
c-------------------------------------------------------------------
c
c.....<t-a|h|i-u>
80    go to 110
c-------------------------------------------------------------------
c
c.....<t-a|h|i-b> element
c     -------------------
90    it=ip
      ia=iq
      ii=ir
      ib=is
      if(isymmo(it).ne.isymmo(ii).or.ia.ne.ib)go to 110
      m1=ic1e(it,1)+ic1e(ii,2)
      sxhsx=-(popsqn(it)/sq2)*(skv(m1)+skv(n1e+m1))
      go to 110
c------------------------------------------------------------------
c
c.....<t-a|h|u-b> element
c     -------------------
100   it=ip
      ia=iq
      iu=ir
      ib=is
      sxhsx=0.0d0
      if(isymmo(ia).ne.isymmo(ib).or.isymmo(it).ne.isymmo(iu))
     +   go to 110
      m1=ic1e(max(ia,ib),1)+ic1e(min(ia,ib),2)
      m1a=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      if(it.eq.iu)sxhsx=sxhsx+skv(m1)+skv(n1e+m1)
      if(ia.ne.ib)go to 110
      temp=skv(m1a)+skv(n1e+m1a)
      if(it.ne.iu)temp=temp*popsqn(it)*popsqn(iu)*0.5d0
      sxhsx=sxhsx-temp
c.....test for diagonal element
      if(kk.ne.ll) go to 110
      sxhsx=sxhsx+vlev+ve0
      vd=sxhsx
      go to 110
c----------------------------------------------------------------------
110   vs=vs+sxhsx*scc(ll)
c
c
c
120   continue
      return
      end
_ENDEXTRACT
      subroutine tractl(q,iwr,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/ciconv)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/blksiz)
      common/cone   /m,na,nb,ispre(6),iperm(8,8)
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
INCLUDE(common/files)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsed,nprimp,naa
     +              ,nab,nst,n1e,isym,mults,macro,maxc
INCLUDE(common/gjs)
      character*10 charwall
       dimension q(*)
c...    keep track of # cycles
      save nnncyc
c
       data m500/500/
c
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,1)top ,charwall()
 1    format( 1x,
     *'commence integral transformation at ',f8.2,' seconds',
     *a10,' wall')
      nsa4=nbas4
c     ion=1
      nfiles=1
      nfilef=1
      indxi=1
      indxj=1
      irrr=1
      iss=1
      master=0
      junits=n4tape(1)
      jblkas=n4blk(1)
      jblkrs=jblkas-n4last(1)
      junitf=n6tape(1)
      jblkaf=n6blk(1)
      jblkrf=jblkaf-n6last(1)
      call secput(isect(484),m500,1,isecbl)
      icyc=macro
      if(macro.gt.n20cas)icyc=n20cas
cjvl
cjvl  reorthogonalise every 5th iteration
cjvl
c     if ((macro/5)*5.eq.macro) 
c    +  call casort(q,q(nblkq4+1),q(nblkq4+lenb4+1))
c
      if (osuped) then
c... dynamic switching of super to nr 
c... mode : 0 : superci ; 1 newton-raphson ; 2 : nr+hessian-construction
c... test on fmax
         mode = 0
         if (fmax.lt.swnr) mode = 2
         if (mode.eq.0) then
            nnncyc = 0
         else
            nnncyc = nnncyc + 1
            if (nnncyc.ne.1) mode = 1
            if (nnncyc.eq.nhesnr) nnncyc = 0
         end if
         if (fmax.lt.swsimu) then
            isimul(icyc) = 1
         else
            isimul(icyc) = 0
         end if
      else
         mode=isort(icyc)
      end if
      call index4(q,macro)
      i10=nblkq4+1
      i20=i10+lenb4
      i30=i20+lenb4
      i40=i30+lenb4
      call hambld(q(i10),q(i30),q(i40),q)
      call revind
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,2)top ,charwall()
 2    format(/1x,
     *'integral transformation complete at ',f8.2,' seconds',
     *a10,' wall')
      return
      end
      subroutine index4(q,macro)
c...   control routine for the 4-index transformation routines
c...   sort1,calc1,sort2,calc2
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/machin)
      common/craypk/kbcray(2044)
INCLUDE(common/gjs)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
      common/disc  /isel(5),ipos(maxlfn)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
INCLUDE(common/iofile)
      common/junke /maxt,ires,ipass,nteff,npass1,npass2,lentri,
     +              nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/block /length,nsymm(8)
INCLUDE(common/mapper)
INCLUDE(common/files)
      common/blkin/gout(511)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/trntim)
INCLUDE(common/restar)
INCLUDE(common/filel)
      character*10 charwall
      dimension q(*)
      data m1,m2/1,2/
c
      oprint=macro.ne.1.and.mode.lt.2
      slder=cpulft(1)
      call ibasgn(1500,0,nsz340,ibase)
      call ibasgn(1500,0,nsz680,ibasen)
      call setsto(2044,0,kbcray)
      niqq=nblkq4
      do  i=1,nirr
       niqq=niqq+nsymm(i)**2
      enddo
      maxt=(lword4-niqq)/lenb4
_IF(ibm,vax)
      nword=(lword4*2)/3
      nword=(nword/4)*4
      nword2=nword/4
_ELSE
      if(o255i) then
       nword=lword4/3
      else
       nword=lword4/2
      endif
      nword=(nword/2)*2
      nword2=nword
_ENDIF
      ires=nword/nsz340
      if(ires.lt.1.or.maxt.lt.1)call caserr('not enough store')
      if(ires.gt.1500)ires=1500
       if(master.eq.lenb4)goto 7766
c...   determine min. no. of passes for sort1/calc1
      i=lenb4-master-1
102   nteff=i/npass1+1
      if(((nteff-1)/ires).lt.maxt)goto 103
      npass1=npass1+1
      goto 102
103   if(npass1.gt.i)npass1=i+1
      if(nprint.eq.-5)go to 104
      write(iwr,205)slder ,charwall()
      if(macro.eq.1)write(iwr,204)m1,npass1
 104  lfile=m1file
      do 1000 i=1,lfile
      lotape(i)=m1tape(i)
      liblk(i) = m1blk(i)
 1000 llblk(i) = liblk(i)-m1last(i)
      do 1 ipass=1,npass1
      call setsto(1360,0,kbcray)
      call sort1(q,q(nword+1))
      call calc1(q,q(nblkq4+1),q(niqq+1))
1     call dumper(m1,ipass,npass1,jblkas,junits,iwr)
7766  if(indxi.gt.nsa4)goto 7755
      lfile=m4file
      do 9998 i=1,lfile
      lotape(i)=m4tape(i)
      liblk (i)=m4blk(i)
 9998 llblk(i) =liblk(i)-m4last(i)
      if(oprint)goto4008
      if(nprint.ne.-5)write(iwr,4007)
4007  format(/' status of secondary mainfile'/1x,28('-'))
      if(nprint.ne.-5)call filprn(m4file,m4blk,m4last,m4tape)
4008  lentri=iky(nsa4+1)
c      if(mode.lt.m2)lentri=iky(nmc+1)
c...   determine min. no. of passes for sort2/calc2
      i=lentri-iky(indxi)-indxj
      npasss=npass2
202   nteff=i/npass2+1
      if(((nteff-1)/ires).lt.maxt)goto 203
      npass2=npass2+1
      goto 202
203   if(npass2.gt.i)npass2=i+1
       if(macro.eq.1.and.nprint.ne.-5)write(iwr,204)m2,npass2
204   format(/' no. of sort',i1,' passes=',i4)
205   format(/1x,'commence 4-index transformation at ',f8.2,
     *' seconds',a10,' wall')
       do 2 ipass=1,npass2
       call sort2(q,q(nword+1))
       call calc2(q,q(nblkq4+1),q(niqq+1))
2      call dumper(m2,ipass,npass2,jblkaf,junitf,iwr)
      npass2=npasss
7755  slder=cpulft(1)
       call clredx
      if(nprint.ne.-5)write(iwr,300)slder ,charwall()
300    format(/' end of 4-index transformation at',
     *f8.2,' seconds',a10,' wall')
      if(oprint)go to 9977
      if(nprint.ne.-5)write(iwr,3001)
3001  format(/' status of mo integral file'/1x,26('-'))
      if(nprint.ne.-5)call filprn(m6file,m6blk,m6last,m6tape)
9977  call revise
        return
      end
      subroutine sort1(g,nijkl)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c...   sorts mainfile onto lfn sort --- so that for a
c...   given rs comb. al pq combs. available
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
      dimension g(*),nijkl(*)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     2nbuck,mloww,mhi,ntri
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      common/bufb/nwbnwb,lnklnk,gout(1)
INCLUDE(common/cndx40)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/filel)
      common/blkin/gin(1)
_IF(linux)
      external fget
_ENDIF
      data maxb/9999999/
c
c...   determine base and limit triangles for this pass
      mloww=master
      mhi=master+nteff
      mlow=mloww+1
      if(mhi.gt.lenb4)mhi=lenb4
      mtri=mhi-mloww
      nbuck=ires
c...  determine minimum number of buckets
10    ntri=mtri/nbuck
      if(ntri.ge.maxt)goto 20
      nbuck=nbuck-1
      if(nbuck)10,20,10
20    nbuck=nbuck+1
      ntri=mtri/nbuck+1
c...  ntri = maximum number of triangles controlled by 1 bucket
c...  nbuck = number of buckets
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
c...   start loop over mainfile blocks
      do 140 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)40,140,40
 40   iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
 50   call fget(gin,m,iunit)
      if(m)60,140,60
 60   nnn=isort1(itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.2721)
     *  call stackr(g,nijkl)
      lbl=lbl+1
      if(lbl)50,140,50
 140  continue
c...  mainfile now swept
      if(nstack.ne.0)
     *  call stackr(g,nijkl)
c...  clear up output
      do 180 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)150,180,150
150   ib=ibase(ibuck)
      ibn=ibasen(ibuck)
c...  output code
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g(ib+1),1,gout,1)
_IFN1(iv)      call pack(gout(nsz341),lab1632,nijkl(ibn+1),nsz680)
_IF1(iv)      call fmove(nijkl(ibn+1),gout(nsz341),nsz170)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
180   continue
      return
      end
_EXTRACT(calc1,mips4)
      subroutine calc1(q,qt,qq)
c
c     4-index transformation, first phase. machine independent version
c     mode.ne.2 -> (..|..) -> (..|vx)
c     mode.eq.2 -> (..|..) -> (..|vx),(..|vi),(..|ij),(..|va),(..|ia)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),qt(*),qq(*)
INCLUDE(common/sizes)
INCLUDE(common/restar)
      common /block/ length,nsymm(8),nblock(8,maxorb)
      common /lsort/ isymge(8,maxorb),mapto(maxorb),mapbak(maxorb)
     *,v(maxorb)
      common /qpar/ ii4(3),ncore
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/gjs)
      common /junke/ maxt,ires,ipass,nteff,npass1,npass2,lentri,
     +                nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/junk/nwbuck(1500),itx(6800),gtx(3400),
     +            p(maxorb),ncor(8),ini(8),nm(8),nmcp(8),las(8),
     +            nsk,nfk,nsl,nfl,lenk,lenl
      common /bufb/ nkk,mkk,gin(1)
      common /blkin/ gout(510),mword
      common /three/ mark(1500)
INCLUDE(common/mapper)
_IF1(iv)      common /craypk/ kbcray(340),macray(340)
_IFN1(iv)       common/craypk/kbcray(680)
INCLUDE(common/atmblk)
      data maxb/9999999/
c
      call rdedx(q,nblkq4,iblkq4,ndump4)
      nintmx1 = nintmx-1
      omode=mode.ne.2
      last=nsa4
      init=1
      if(omode)last=nmc
      if(omode)init=ncore+1
c
c.....symmetry block the mo's
      ilifm(1)=0
      do isym=1,nirr
      len=nsymm(isym)
      if(len.gt.0)then
       do i=1,len
        ii=nblock(isym,i)
        iii=i+mstart(isym)-1
        ilifm(iii+1)=ilifm(iii)+len
c.....ii-->iii
        call dcopy(len,q(ilifq(ii)+mstart(isym)),1,qq(ilifm(iii)+1),1)
        mapto(ii)=iii
        mapie(iii)=iky(ii)
        mapbak(iii)=ii
       enddo
      endif
      enddo
c
      mapto(nbas4+1)=nbas4+1
      call dcopy(ilifm(nbas4+1),qq,1,q,1)
c
c.....set up base & top for each symmetry
      do 50 l=1,nbas4
      do 50 isym=1,nirr
      if(nsymm(isym).le.0)go to 50
      do 40 k=l,nbas4
      if(isymmo(k).eq.isym)goto 150
40    continue
      k=nbas4+1
150   isymge(isym,mapto(l))=mapto(k)
 50   continue
      do 80 isym=1,nirr
      nm(isym)=0
      las(isym)=0
      ini(isym)=0
      if(nsymm(isym).le.0)go to 80
      nst=mstart(isym)
      nfi=nfin  (isym)
      do 60 l=nst,nfi
      if(mapbak(l).le.last) las(isym)=l
      if(mapbak(l).le.nmc )  nm(isym)=l
60    continue
      do 70 l=nst,nfi
      if(mapbak(l).ge.init)goto 180
70    continue
      go to 80
180   ini(isym)=l
80    continue
c
      mword=0
_IFN1(iv)      int4=-1
c     small=10.0d0**(-iacc-2)
      call stopbk
c
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
      do 260 ibuck=1,nbuck
      mhigh=min(mhi,mloww+ntri)
      mtri=mhigh-mloww
      mloww=mloww+1
c
c...   read in a core load
      call vclr(qq,1,mtri*lenb4)
      mkk=mark(ibuck)
      go to 90
10    mkkkk = mkk
      call rdbak(mkkkk)
      call stopbk
      call transc(qq)
 90   if(mkk.ne.maxb)go to 10
c
      map=0
c.... loop over triangles in this bucket
      do 250 itri=1,mtri
      isymtr=mult(isymao(irrr),isymao(iss))
      master=master+1
c.... loop over symmetry blocks in this triangle
      do 240 isymk=1,nirr
      isyml=mult(isymtr,isymk)
      lenk =nsymm(isymk)
      lenl =nsymm(isyml)
      if(lenk*lenl.le.0)goto 240
      nsk  =mstart(isymk)
      nfk  =nfin  (isymk)
      nsl  =mstart(isyml)
      nfl  =nfin  (isyml)
c.... copy block to rectangle
      mapp=map+iky(nsk)+nsl
      if(nsl.gt.nsk) mapp=map+iky(nsl)+nsk
      call rectri (qt,qq(mapp))
c
c..... now do the transformation
      jfin=las(isyml)
      inik=ini(isymk)
      nmk =nm (isymk)
      if(nmk.le.0.or.inik.le.0.or.inik.gt.nmk)go to 240
      do 230 k=inik,nmk
      jstart=isymge(isyml,k)
      if(jstart.gt.jfin) goto 230
      call mxm (qt,lenl,q(ilifm(k)+1),lenk,p,1)
      lenj=jfin-jstart+1
_IF1(c)      call mxma (q(ilifm(jstart)+1),lenl,1,p,1,lenl,v(jstart),1,lenj,
_IFN1(c)      call mxmaa(q(ilifm(jstart)+1),lenl,1,p,1,lenl,v(jstart),1,lenj,
     1                 lenj,lenl,1)
c.... put result vector y into output buffer
      mapk=mapbak(k)
      jfinb=min(jfin,nintmx1-mword+jstart)
      do j=jstart,jfinb
      mword=mword+1
      gout(mword)=v(j)
_IF1(iv)      kbcray(mword)=mapie(j)+mapk
_IF1(iv)      macray(mword)=master
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mapie(j)+mapk
_IFN1(iv)      kbcray(int4+1)=master
      enddo
      if(mword.ne.nintmx)goto230
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockk
      jfinb1=jfinb+1
      do j=jfinb1,jfin
      mword=mword+1
      gout(mword)=v(j)
_IF1(iv)      kbcray(mword)=mapie(j)+mapk
_IF1(iv)      macray(mword)=master
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mapie(j)+mapk
_IFN1(iv)      kbcray(int4+1)=master
      enddo
230   continue
240   continue
      iss=iss+1
      if(iss.le.irrr)go to 250
      irrr=irrr+1
      iss=1
250   map=map+lenb4
260   mloww=mhigh
      if(mword.lt.1)return
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockk
      return
      end
_ENDEXTRACT
      subroutine rectri (qt,qq)
c
c     copy the block of qq -> rectangel qt. qq & qt both entered at star
c     the block
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
       common/junk/nwbuck(1500),itx(6800),gtx(3400),
     *               p(maxorb),junk(40),nsk,nfk,nsl,nfl,lenk,lenl
      dimension qt(*),qq(*)
c
      if(nsl-nsk) 110,140,170
110   nsize=1
      map=1
      do 130 j=nsk,nfk
      call dcopy(lenl,qq(nsize),1,qt(map),1)
      nsize=nsize+j
130   map=map+lenl
      return
c
170   nsize=1
      map=1
      do 190 i=nsl,nfl
      call dcopy(lenk,qq(nsize),1,qt(map),lenl)
      map=map+1
190   nsize=nsize+i
      return
c
140   nsize=1-nsl
      map1=nsize
      map=nsize
      do 160 i=nsl,nfl
      map2=map+i
      do 150 j=nsl,i
      top=qq(nsize+j)
      qt(map1+j)=top
      qt(map2)=top
150   map2=map2+lenl
      map1=map1+lenl
160   nsize=nsize+i
      return
      end
      subroutine sort2(g,nijkl)
c...   sorts secondary mainfile onto lfn sort---so that
c...   for a given ij comb. all rs combs. available
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
      dimension g(*),nijkl(*)
INCLUDE(common/sizes)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
INCLUDE(common/stak)
INCLUDE(common/mapper)
      common/junke /maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     2nbuck,mloww,mhi,ntri,iacc
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      common/bufb/nwbnwb,lnklnk,gout(1)
INCLUDE(common/cndx40)
       common/blkin/gin(1)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/filel)
_IF(linux)
      external fget
_ENDIF
      data maxb/9999999/
c
c...   determine base and limit triangles for this pass
       mloww=iky(indxi)+indxj-1
      mlow=mloww+1
       mhi=mloww+nteff
      if(mhi.gt.lentri)mhi=lentri
      mtri=mhi-mloww
c...  determine minimum number of buckets
      nbuck=ires
10    ntri=mtri/nbuck
      if(ntri.ge.maxt)goto 20
      nbuck=nbuck-1
      if(nbuck)10,20,10
20    nbuck=nbuck+1
      ntri=mtri/nbuck+1
c...  ntri = max. number of triangles controlled by 1 bucket
c...  nbuck = number of buckets
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
c...  start loop over secondary mainfile blocks
      do 140 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      lbl=llblk(ifile)
40    call fget(gin,m,iunit)
      if(m)50,140,50
 50   nnn=isort2(itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.3061)
     * call stackr(g,nijkl)
      lbl=lbl+1
      if(lbl)40,140,40
 140  continue
c...  secondary mainfile now swept
      if(nstack.ne.0)
     *  call stackr(g,nijkl)
c...  clear up output
      do 180 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)150,180,150
150   ib=ibase(ibuck)
      ibn=ibasen(ibuck)
c...  output code
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g(ib+1),1,gout,1)
_IF(ibm,vax)
      call fmove(nijkl(ibn+1),gout(nsz341),nsz170)
_ELSE
      call pack(gout(nsz341),lab1632,nijkl(ibn+1),nsz680)
_ENDIF
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
180   continue
      return
      end
      subroutine calc2(q,qt,qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
      common /block / length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/gjs)
      common /junke/ maxt,ires,ipass,
     1               nteff,npass1,npass2,lentri,
     2               nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common /lsort/ isymle(8,maxorb),mapto(maxorb),mapbak(maxorb)
     *,v(maxorb)
      common /bufb/nkk,mkk,gin(1)
      common /junk/nwbuck(1500),itx(6800),gtx(3400),
     +            p(maxorb),ncor(8),ini(8),nm(8),nmcp(8),las(8),
     +            nsk,nfk,nsl,nfl,lenk,lenl
      common /qpar/ ii4(3),ncore
INCLUDE(common/mapper)
      common /blkin/ gout(510),mword
_IF1(iv)      common /craypk/ kbcray(340),macray(340)
_IFN1(iv)      common /craypk/ kbcray(680)
      common /three/ mark(1500)
INCLUDE(common/atmblk)
      dimension q(*),qt(*),qq(*)
      data maxb/9999999/
c
      omode=mode.ne.2
      init=1
      if(omode)init=ncore+1
      oomode=.not.omode.or.ncore.eq.0
c     nmc1=nmc+1
      acc1=10.0d0**(-iacc)
      call rdedx(q,nblkq4,iblkq4,ndump4)
c
c.....symmetry block the mo's
      ilifm(1)=0
      do 10 isym=1,nirr
      len=nsymm(isym)
      if(len.le.0)go to 10
      do 130 i=1,len
      ii=nblock(isym,i)
      iii=i+mstart(isym)-1
      ilifm(iii+1)=ilifm(iii)+len
c.....ii-->iii
      call dcopy(len,q(ilifq(ii)+mstart(isym)),1,qq(ilifm(iii)+1),1)
      mapto(ii)=iii
130   mapbak(iii)=ii
10    continue
      mapto(nbas4+1)=nbas4+1
      call dcopy(ilifm(nbas4+1),qq,1,q,1)
c
c.....set up base & top for each symmetry
      do 30 l=1,nbas4
      do 30 isym=1,nirr
      isymle(isym,mapto(l))=0
      do 30 k=1,l
      if(isymmo(k).eq.isym)isymle(isym,mapto(l))=mapto(k)
30    continue
      do 80 isym=1,nirr
      ncor(isym)=0
      nm(isym)=0
      las(isym)=0
      ini(isym)=0
      nmcp(isym)=0
      if(nsymm(isym).le.0)go to 80
      nst = mstart(isym)
      nfi = nfin  (isym)
      do 40 l=nst,nfi
      if(mapbak(l).le.ncore)ncor(isym)=l
      if(mapbak(l).le.nmc )   nm(isym)=l
      if(mapbak(l).le.nbas4)las(isym)=l
40    continue
      do 50 l=nst,nfi
      if(mapbak(l).ge.init)goto 60
50    continue
      go to 61
60    ini(isym)=l
 61   do 70 l=nst,nfi
      if(mapbak(l).gt.nmc)goto 180
70    continue
      go to 80
180   nmcp(isym)=l
80    continue
c
      mword=0
_IFN1(iv)      int4=-1
      call stopbk
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
      do 380 ibuck=1,nbuck
      mhigh=min(mhi,mloww+ntri)
      mtri=mhigh-mloww
      mloww=mloww+1
      call vclr(qq,1,mtri*lenb4)
c...   read in a core load
      mkk=mark(ibuck)
110   if(mkk.eq.maxb)go to 120
      mkkkk = mkk
      call rdbak(mkkkk)
      call stopbk
      call transc(qq)
      goto 110
c
120   map=0
      do 370 itri=1,mtri
      if((omode.and.indxi.gt.nmc).or.indxj.gt.nmc.or.indxj.lt.init)
     +            go to 360
      mword1=indxj+i4096(indxi)
      isymtr=mult(isymmo(indxi),isymmo(indxj))
      do 350 isymk=1,nirr
      isyml=mult(isymtr,isymk)
      lenk=nsymm(isymk)
      lenl=nsymm(isyml)
      if(lenk*lenl.eq.0)goto 350
      nsk=mstart(isymk)
      nfk=nfin(isymk)
      nsl=mstart(isyml)
      nfl=nfin(isyml)
c.... copy triangle block of qq into qt
      mapp=map+iky(nsk)+nsl
      if(nsl.gt.nsk) mapp=map+iky(nsl)+nsk
      call rectri (qt,qq(mapp))
c
      if(indxi.gt.nmc)go to 320
c
c.....(tu|vx) and (tu|va)/(tu|vi)
      jstart=mstart(isyml)
      jfirst=nmcp(isyml)
      jfin  =las(isyml)
      inik = ini(isymk)
      nmk  = nm (isymk)
      if(nmk.eq.0.or.inik.gt.nmk. or. inik.eq.0)go to 351
      do 290 k=inik,nmk
      mapbk=mapbak(k)
      call mxm(qt,lenl,q(ilifm(k)+1),lenk,p,1)
c
      last=isymle(isyml,k)
      if(mapbk-indxi)250,240,230
230   if(oomode)goto 270
      last=ncor(isyml)
      goto 250
240   last=isymle(isyml,mapto(indxj))
250   lenj=last-jstart+1
      if(lenj.lt.1)goto 270
      mapk=i4096(mapbk)
      if(jstart.gt.last. or. last.eq.0)go to 270
c     write(6,5000)inik,nmk,k,jstart,last
      do 260 j=jstart,last
c      top=y(j)
c.....lenj internal orbitals only, so this maximises vector length
      top=ddot(lenl,q(ilifm(j)+1),1,p,1)
      if( dabs(top).lt.acc1)go to 260
      mword=mword+1
c     write(6,5001)mword,jstart,last,j,indxi,indxj,mapbk,mapbak(j),top
      gout(mword)=top
_IF1(iv)      kbcray(mword)=mword1
_IF1(iv)      macray(mword)=mapbak(j)+mapk
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mword1
_IFN1(iv)      kbcray(int4+1)=mapbak(j)+mapk
      if(mword.ne.nintmx)goto260
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockl
260   continue
c.....(tu|va)
270   lenj=jfin+1-jfirst
      if(jfirst.gt.jfin.or.jfirst.eq.0)go to 290
_IF1(c)      if(lenj.ge.1)call mxma (q(ilifm(jfirst)+1),lenl,1,p,1,lenl,
_IFN1(c)      if(lenj.ge.1)call mxmaa(q(ilifm(jfirst)+1),lenl,1,p,1,lenl,
     1             v(jfirst),1,lenj,     lenj,lenl,1)
c     write(6,5002)inik,nmk,k,jfirst,jfin
      do 280 j=jfirst,jfin
      top=v(j)
      if( dabs(top).lt.acc1)go to 280
      mword=mword+1
c     write(6,5001)mword,jfirst,jfin,j,mapbak(j),mapbak(k),indxi,
c    * indxj,top
      gout(mword)=top
_IF1(iv)      macray(mword)=mword1
_IF1(iv)      kbcray(mword)=mapbk+i4096(mapbak(j))
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mapbk+i4096(mapbak(j))
_IFN1(iv)      kbcray(int4+1)=mword1
      if(mword.ne.nintmx)goto280
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockl
280   continue
290   continue
351   if(omode)goto350
c.....(tu|ab)
      nst = nmcp(isymk)
      nfi = las (isymk)
      if(nfi.eq.0.or.nst.gt.nfi.or.nst.eq.0)go to 350
      do 310 k=nst,nfi
      jend=isymle(isyml,k)
      if(jfirst.gt.jend.or.jend.eq.0.or.jfirst.eq.0) goto 310
      call mxm(qt,lenl,q(ilifm(k)+1),lenk,p,1)
      lenj=jend+1-jfirst
c     write(6,6000)nst,nfi,k,isymk,isyml,jend,jfirst,lenj
_IF1(c)      call mxma (q(ilifm(jfirst)+1),lenl,1,p,1,lenl,v(jfirst),1,lenj,
_IFN1(c)      call mxmaa(q(ilifm(jfirst)+1),lenl,1,p,1,lenl,v(jfirst),1,lenj,
     1           lenj,lenl,1)
      mapk=i4096(mapbak(k))
      do 300 j=jfirst,jend
      top=v(j)
      if( dabs(top).lt.acc1)go to 300
      mword=mword+1
      gout(mword)=top
c     write(6,6001)j,mapbak(k),mapbak(j),indxi,indxj,top
_IF1(iv)      macray(mword)=mword1
_IF1(iv)      kbcray(mword)=mapbak(j)+mapk
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mapbak(j)+mapk
_IFN1(iv)      kbcray(int4+1)=mword1
      if(mword.ne.nintmx)goto300
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockl
300   continue
310   continue
      goto 350
c....(av|bx)
320   if(omode)goto350
      last=isymle(isyml,mapto(indxi))
      lastm=isymle (isyml,mapto(indxi-1))
      jstart=nmcp(isyml)
      nmcl=nm(isyml)
      nst = mstart(isymk)
      nfi = nm(isymk)
      if(nst.gt.nfi.or.nfi.eq.0)go to 350
      do 340 k=nst,nfi
      mapk=mapbak(k)
c     write(6,6002)last,lastm,jstart,nmcl,nst,nfi,k,mapk
      if(mapk.gt.indxj)last=lastm
      if(last.eq.nmcl.or.jstart.gt.last.or.jstart.eq.0)goto 340
      call mxm (qt,lenl,q(ilifm(k)+1),lenk,p,1)
      lenj=last+1-jstart
_IF1(c)      call mxma (q(ilifm(jstart)+1),lenl,1,p,1,lenl,v(jstart),1,lenj,
_IFN1(c)      call mxmaa(q(ilifm(jstart)+1),lenl,1,p,1,lenl,v(jstart),1,lenj,
     1           lenj,lenl,1)
      do 330 j=jstart,last
      top=v(j)
      if( dabs(top).lt.acc1)go to 330
      mword=mword+1
      gout(mword)=top
c     write(6,6003)j,jstart,last,indxi,indxj,mapbak(j),mapk,top
_IF1(iv)      kbcray(mword)=mword1
_IF1(iv)      macray(mword)=mapk+i4096(mapbak(j))
_IFN1(iv)      int4=int4+2
_IFN1(iv)      kbcray(int4  )=mword1
_IFN1(iv)      kbcray(int4+1)=mapk+i4096(mapbak(j))
      if(mword.ne.nintmx)goto330
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockl
330   continue
340   continue
350   continue
360   indxj=indxj+1
      if(indxj.le.indxi)goto 370
      indxi=indxi+1
      indxj=1
370   map=map+lenb4
380   mloww=mhigh
      if(ipass.eq.npass2)call rdedx(q,nblkq4,iblkq4,ndump4)
      if(mword.lt.1)return
_IF1(iv)      call pak4v(kbcray,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,kbcray,numlabp)
_IFN1(iv)      int4=-1
      call blockl
      return
      end
      subroutine hambld(q1,q3,q4,q2)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/dump3)
      common/potn  /core,potnuc
INCLUDE(common/gjs)
INCLUDE(common/cndx40)
      common/blkin /corev(512),pp(4),v(508)
INCLUDE(common/machin)
INCLUDE(common/mapper)
INCLUDE(common/qice)
      common/qpar  /i4(3),ncore
      dimension q1(*),q3(*),q2(*),q4(*),o1e(6)
c
c     x3/x4  and  x1/x5  adjacent in core
c
      call rdedx(q2(1),nblkq4,iblkq4,ndump4)
c
      do loop =1,6
       o1e(loop) = .false.
      enddo
      o1e(3) = .true.
      call getmat(q3,q3,q3,q3,q3,q3,pp,nbas4,o1e,ions2)
c .. ao matrix loaded int q3
      call vclr(q4,1,lenb4)
      lenb2=lenb4+lenb4
      call vclr(q1,1,lenb2)
       if(ncore.eq.0)go to 50
      do 40 k=1,ncore
      mm=ilifq(k)
      m=0
      do 30 i=1,nbas4
      top=q2(mm+i)
      v(i)=top
      top=top+top
      if(top.ne.0.0d0)  then
      call daxpy(i,top,v,1,q4(m+1),1)
        endif
 30   m=m+i
40    continue
      do 3 i=1,nbas4
      mm=iky(i+1)
    3 q4(mm)=q4(mm)*0.5d0
c...   zero order rho in x1
c...   process mainfile
c....core energy
 50   potnuc=pp(1)
      core=potnuc+ddot(lenb4,q4,1,q3,1)
      if(ncore.le.0)go to 70
      call scang(q1,q4)
      call vadd(q3,1,q1,1,q3,1,lenb4)
 70   continue
      core=core+ddot(lenb4,q4,1,q3,1)
      call mult2(q2,q1,q3,ncol4,ncol4,nbas4)
c ... required ints in x1
c.....compress matrix
      nn=0
      do 60 i=1,nbas4
      do 60 j=1,i
      nn=nn+1
      if(isymmo(i).ne.isymmo(j))go to 60
      q1(ic1e(i,1)+ic1e(j,2))=q1(nn)
60    continue
      return
      end
      subroutine dumper(iphase,ipass,npass,jbl,jun,iwr)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
INCLUDE(common/files)
INCLUDE(common/cndx40)
      common/blkin /gout(511)
      data m0/0/
c
      if(ipass.ne.npass)goto 88
      call search(jbl,jun)
      call put(gout,m0,jun)
      jbl=jbl+1
 88   if(iphase-1)89,89,90
 89   m4file=nfiles
      m4tape(m4file)=n4tape(nfiles)
      m4blk( m4file)=n4blk (nfiles)
      m4last(m4file)=jbl
      go to 66
 90   m6file=nfilef
      m6tape(m6file)=n6tape(nfilef)
      m6blk( m6file)=n6blk (nfilef)
      m6last(m6file)=jbl
 66   call search(isecbl,ndump4)
      m12=12/lenwrd()
      call put(master,m12,ndump4)
      b=cpulft(1)
      call clredx
      if(nprint.ne.-5)write(iwr,100)iphase,ipass,b
100   format(/' end of sort',i1,' pass',i6,' at ',f8.2,' seconds')
        return
      end
_IF()
      subroutine casort(q,s,scra)
c
cjvl  orbital reorthogonalisation for casa
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*),s(*),scra(*),o1e(6)
c
INCLUDE(common/cndx40)
INCLUDE(common/dump3)
      common/blkin /corev(512),pp(4),v(508)
c
c...   read vectors
      call rdedx(q(1),nblkq4,iblkq4,ndump4)
c...   get s-matrix
      do loop =1,6
       o1e(loop) = .false.
      enddo
      o1e(1) = .true.
      call getmat(s,s,s,s,s,s,pp,nbas4,o1e,ions2)
c
c...   schmidt orthonormalise
c     
      call normvc(q,s,scra,nbas4,nbas4,1.0d-6)
c...  write vectors
      call wrt3(q(1),nblkq4,iblkq4,ndump4)
c
      return
      end
_ENDIF
      subroutine normvc( vc, s, w, ndim,nmdim, critlow)
c
c...  schmidt othogonalisation of nmdim vectors of length ndim
c
c..   vc  ndim*nmdim
c..   s   ndim*(ndim+1)/2
c..   w   ndim*2
c
      implicit REAL (a-h,o-z)
      dimension vc(*),s(*),w(*)
c
      if (nmdim.le.0) return
c
      if (ndim .gt. 1) go to 100
      if (vc(1) * vc(1) .le. critlow) goto 160
      vc(1) = 1.0d0/dsqrt(s(1))
      return
  100 np1 = ndim+1
      ivci = 1
      do 150 i=1,nmdim
      call cntrc(s,vc(ivci),w,ndim,critlow)
      t = 0
      tnorm = 0
      ivcj = 1
      ivcw = ndim
      do 120 j=1,i
      tnorm = tnorm - t*t
      t = 0
      do 110 jj=1,ndim
      t = vc(ivcj)*w(jj) + t
  110 ivcj = ivcj + 1
      ivcw = ivcw + 1
  120 w(ivcw) = -t
      if (tnorm + t .le. critlow) goto 160
      tnorm = 1/dsqrt(tnorm+t)
      w(ivcw) = 1
      do 140 k = 1,ndim
      t = 0
      ivcw = np1
      do 130 j=k,ivci,ndim
      t = vc(j)*w(ivcw) + t
  130 ivcw = ivcw + 1
      vc(ivci) = t*tnorm
  140 ivci = ivci + 1
  150 continue
      return
  160 call caserr(' orthogonalisation failure - normvc ')
      return
      end
      subroutine cntrc(s, c, v, n, eps)
c
c     cntrc computes the inproducts of the vector c with the rows
c     of the triangular matrix s
c     the results are stored in the vector v
c
      implicit REAL (a-h,o-z)
      dimension s(*),c(*),v(*)
      do 10 mu=1,n
10    v(mu) = 0.0d0
      ix=0
      do 50 nu=1,n
      cnu=c(nu)
      if (dabs(cnu) .ge. eps) then
        munu=ix
        do 20 mu=1,nu
        munu=munu+1
   20   v(mu)=v(mu)+s(munu)*cnu
        if(nu.eq.n) go to 50
        nu1=nu+1
        munu=munu+nu
        do 30 mu=nu1,n
        v(mu)=v(mu)+s(munu)*cnu
   30   munu=munu+mu
      endif
   50 ix=ix+nu
      return
      end
      subroutine ver_casa(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/casa.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
