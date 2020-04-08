c 
c  $Author: jmht $
c  $Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
c  $Locker:  $
c  $Revision: 6090 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util5.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   util5(ci) =
c ******************************************************
c ******************************************************
_IF1()      subroutine gathr(n,r,a,map)
_IF1()      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      dimension r(*),a(*)
_IF1()      if(n.le.0)go to 2
_IF1()      ii=1
_IF1()      do 1 loop=1,n
_IF1()      r(loop)=a(ii)
_IF1()   1  ii=ii+map
_IF1()   2  return
_IF1()      end
_IF1()      subroutine vivits(n,scalar,r,mr,a,ma)
_IF1()      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      dimension a(*),r(*)
_IF1()      iii=1
_IF1()      jjj=1
_IF1()      do 1 loop=1,n
_IF1()      r(iii)=a(jjj)*scalar
_IF1()      iii=iii+mr
_IF1()    1 jjj=jjj+ma
_IF1()    2 return
_IF1()      end
      subroutine pinver(ip,ipb,n)
c
c     invert permutation ip ==> ipb
c
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension ip(n),ipb(n)
c
      do 10 i=1,n
10    ipb(ip(i)) = i
c
      return
      end
      subroutine anti1(r,a,n)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*),a(n,*)
      m=1
      do 1 i=2,n
      im1=i-1
      do 1 j=1,im1
      r(m)=a(j,i)-a(i,j)
    1 m=m+1
      return
      end
_IF1()      subroutine imove(ia,ib,n)
_IF1()      implicit REAL  (a-h,o-z),integer  (i-n)
_IF1()      dimension ia(*),ib(*)
_IF1()      if(n.le.0)return
_IF1()      do 1 loop=1,n
_IF1()    1 ib(loop)=ia(loop)
_IF1()      return
_IF1()      end
_IFN(cray,t3d)
      subroutine zbasgn(n,ibegin,iadd,iarr)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(ksr,convex,i8,i8drct)
_IF(convex,i8drct)
      integer *8 iarr
_ELSE
      integer iarr
_ENDIF
_ELSE
      REAL  iarr
_ENDIF
      dimension iarr(*)
      incr=ibegin
      do 1 loop=1,n
_IF(ksr,convex,i8,i8drct)
      iarr(loop)=incr
_ELSEIF(bits8)
      call pad8(incr,iarr(loop))
_ELSE
      iarr(loop)=pad(incr)
_ENDIF
    1 incr=incr+iadd
      return
      end
_ENDIF
_IFN(cray,t3d)
_IF(i8)
      subroutine smask
      implicit real*8  (a-h,m-z),integer    (i-l)
      integer mask,z1
      common/maskc/mask(64)
      data z1/z'8000000000000000'/
      mask(1)=z1
      do 1 i=2,64
    1 mask(i)=ior(mask(i-1),ishft(z1,1-i))
      return
      end
_ELSE
      subroutine smask
      implicit REAL  (a-h,m-z),integer    (i-l)
c...   note default mask : real
_IF(convex,i8drct)
      integer *8 mask,z1
_ELSEIF(ksr)
      integer mask,z1
_ELSE
_ENDIF
      common/maskc/mask(64)
_IF(ksr,convex,i8drct)
      data z1/z'8000000000000000'/
_ELSEIF(ibm)
      data z1/z8000000000000000/
_ELSEIF(bits8)
      call pad8(1,tmp1)
      call shift8(tmp1,63,z1)
_ELSE
_IFN1(gh)      z1 = shift(pad(1),63)
_IF1(gh)      z1 = shiftc(pad(1),63)
_ENDIF
      mask(1)=z1
      do 1 i=2,64
_IF(convex,i8drct)
    1 mask(i)=ior(mask(i-1),ishft (z1,1-i))
_ELSEIF(bits8)
      call shiftr8(z1,i-1,tmp8)
    1 call dxor8(mask(i-1),tmp8,mask(i))
_ELSE
    1 mask(i)=dor(mask(i-1),shiftr(z1,i-1))
_ENDIF
      return
      end
_ENDIF
_ENDIF
      subroutine isquar(ir,ia,mrowr,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ir(*),ia(*)
c... convert triangle ia to square ir
      iii=0
      jjj=1
      k=1
      do 40 i=1,n
      jj=jjj
_IF1(t)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
      do 30 j=1,i
      ir(iii+j)=ia(k)
      ir(jj)=ia(k)
      k=k+1
30    jj=jj+mrowr
      iii=iii+mrowr
40    jjj=jjj+1
      return
      end
      subroutine sqsitr(r,s,t,nk)
      implicit REAL  (a-h,o-z), integer  (i-n)
      dimension s(*),t(*)
      dimension r(nk,*)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
      do 1 i=1,nk
    1 r(i,i)=s(iky(i+1))
      m=1
      do 2 j=2,nk
      l=iky(j)
      j1=j-1
_IF1(t)cdir$ ivdep
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
      do 2 i=1,j1
      r(i,j)=s(l+i)+t(m)
      r(j,i)=s(l+i)-t(m)
   2  m=m+1
      return
      end
      subroutine symm1a(r,a,nk)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=1,nk
      do 1 j=1,i
      r(m)=r(m)+a(j,i)+a(i,j)
    1 m=m+1
      return
      end
      subroutine anti1s(r,a,nk)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=2,nk
      i1=i-1
      do 1 j=1,i1
      r(m)=r(m)+a(i,j)-a(j,i)
    1 m=m+1
      return
      end
      subroutine anti1a(r,a,nk)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*)
      dimension a(nk,*)
      m=1
      do 1 i=2,nk
      i1=i-1
      do 1 j=1,i1
      r(m)=r(m)-a(i,j)+a(j,i)
    1 m=m+1
      return
      end
_IF(ibm,vax)
      subroutine transe(qq,imap,local,maxtrs)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/blksiz)
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      integer*2 lin
      dimension qq(*),imap(*),local(*)
      dimension lin(2)
      equivalence (lin(1),gin(1))
      ij=nszij
      kl=nszkl
      do 1 loop=1,nkk
      qq(((imap(lin(ij))-mloww)*maxtrs+local(lin(ij)))*lenbas
     * + lin(kl)) = gin(loop)
      ij=ij+1
    1 kl=kl+1
      return
      end
      subroutine transb(qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qq(*)
      integer*2 lin
INCLUDE(common/blksiz)
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      dimension lin(2)
      equivalence (lin(1),gin(1))
      ij=nszij
      kl=nszkl
      do loop=1,nkk
       qq((lin(ij)-mloww)*lenbas + lin(kl)) = gin(loop)
       ij=ij+1
       kl=kl+1
      enddo
      return
      end
_ELSE
      subroutine transe(qq,imap,local,maxtrs)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qq(*),imap(*),local(*)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      common/scra/ijklin(2,3412),index(3412)
c
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=((imap(ijklin(1,loop))-mloww)*maxtrs+
     +  local(ijklin(1,loop)))* lenbas+ijklin(2,loop)
      enddo
_IF(cray)
      call scatter(nkk,qq,index,gin)
_ELSE
      call dsctr(nkk,gin,index,qq)
_ENDIF
      return
      end
      subroutine transb(qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qq(*)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt(6),lenbas,nb,mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
      common/scra/ijklin(2,3412),index(3412)
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=(ijklin(1,loop)-mloww)*lenbas + 
     +              ijklin(2,loop)
      enddo
_IF(cray)
      call scatter(nkk,qq,index,gin)
_ELSE
      call dsctr(nkk,gin,index,qq)
_ENDIF
      return
      end
_ENDIF
      subroutine jsortp(map)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/scra/ibuk(3400)
_IFN1(iv)     * ,             index(3400)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      dimension map(*)
      mstack=1
_IF1(c)      call gather(nstack,index,map,itx)
_IFN1(civ)      call igthr(nstack,map,index,itx)
      do 1 i=1,nstack
_IF1(iv)      ibuk(i)=(map(itx(i))-mlow)*btri+1
_IFN1(iv)      ibuk(i)=(index(i)-mlow)*btri+1
 1    continue
      return
      end
_IF(ibm,vax)
      subroutine pakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer * 2 intout
      dimension intin(*),intout(*)
      mword4=mword*4
      do 1 i=1,mword4,4
      intout(i  )=intin(i  )
      intout(i+1)=intin(i+1)
      intout(i+2)=intin(i+2)
    1 intout(i+3)=intin(i+3)
      return
      end
      subroutine upakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *2  intin
      dimension intin(*),intout(*)
      dimension i4(4)
      data i4/0,0,0,0/
      mword4=mword*4
      do 1 i=1,mword4,4
      i4(1)=intin(i)
      intout(i)=i4(1)
      i4(2)=intin(i+1)
      intout(i+1)=i4(2)
      i4(3)=intin(i+2)
      intout(i+2)=i4(3)
      i4(4)=intin(i+3)
    1 intout(i+3)=i4(4)
      return
      end
_ENDIF
_IF(cray,ksr,apollo,t3d)
_IF1() used in mcscf // only 32 bits integers
_IF1() assume for now that all others are like apollo
      logical function btest(i,ibit)
      implicit REAL  (a-h,o-z)
_IF(cray,t3d)
      j = shiftl(1,ibit)
      k = i.and.j
      btest = k.ne.0
_ELSEIF(ksr)
      j = shiftl(1,ibit)
      k=and(i,j)
      btest = k.ne.0
_ELSE
      j = lshft(1,ibit)
      k = and(i,j)
      btest = k.ne.0
_ENDIF
      return
      end
      function ibset(i,ibit)
      implicit REAL  (a-h,o-z)
_IF(cray,t3d)
      ibset = i.or.shiftl(1,ibit)
_ELSEIF(ksr)
      ibset = or (i,shiftl(1,ibit))
_ELSE
      ibset = or(i,lshft(1,ibit))
_ENDIF
      return
      end
      function ibclr(i,ibit)
      implicit REAL  (a-h,o-z)
_IF(cray,t3d)
      j = .not.shiftl(1,ibit)
      ibclr = i.and.j
_ELSEIF(ksr)
      j = not(shiftl(1,ibit))
      ibclr = and(i,j)
_ELSE
      j = not(lshft(1,ibit))
      ibclr = and(i,j)
_ENDIF
      return
      end
_ENDIF
      subroutine mcsrto (g,nijkl)
c
c...  subroutine to handle sortfile i/o. hopefully most machine dependnt
c     features are isolated here.
      implicit REAL  (a-h,o-z)
_IF1(iv)      integer *2 lout
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
INCLUDE(common/stak)
      common /bufb/ nwbnwb,lnklnk,gout(5118)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      dimension g(nwb),nijkl(*)
_IF1(iv)      dimension lout(2)
_IF1(iv)      equivalence (lout(1),gout(1))
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
      call dcopy(nwb,g,1,gout,1)
_IF1(iv)      ista=nszij
_IF1(iv)      do 10 i=1,nwb
_IF1(iv)      lout(ista  )=nijkl(i+i-1)
_IF1(iv)      lout(ista+1)=nijkl(i+i)
_IF1(iv)   10 ista=ista+2
_IFN1(iv)      call pack(gout(nsz341),lab1632,nijkl,nsz680)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
      nwb=0
      nwbuck(ibuck)=0
      return
      end
      subroutine mcsrti(mk,nk)
c
c...  subroutine to handle sortfile input. hopefully most machine
c     dependent features are isolated here.
      implicit REAL  (a-h,o-z)
_IF1(iv)      integer * 2 lout
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
INCLUDE(common/stak)
      common /bufb  / nwbnwb,lnklnk,gout(5118)
      common /junk  /klin(2,3412)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
_IF(ibm,vax)
      dimension lout(2)
      equivalence (gout(1),lout(1))
_ENDIF
      iblock = mk
      call rdbak(iblock)
      call stopbk
      mk=lnklnk
      nk=nwbnwb
_IF(ibm,vax)
      ista=nszij
      do 20 int=1,nwbnwb
      klin(1,int) = lout(ista)
      klin(2,int) = lout(ista+1)
   20 ista     = ista+2
_ELSE
      call unpack(gout(nsz341),lab1632,klin,nsz680)
_ENDIF
      if(odebug(34)) then
       write(6,*)'mcsrti: gout, kin, lin'
       do int = 1,nwbnwb
        write(6,100)int,gout(int),klin(1,int),klin(2,int)
100     format(1x,i5,f20.9,2(1x,i7))
       enddo
      endif
c
      return
      end
      subroutine mult3 (q,h,r,p,nactiv,nbasis)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension q(*),h(*),r(*),p(*)
INCLUDE(common/mlngth)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
      do 1 loop=1,nbasis
      m=iky(loop+1)
    1 h(m)=h(m)*0.5d0
      nacnba=nactiv*nbasis
      if(nbasd2.ge.nactiv)go to 2
      call vclr(p,1,nacnba+nactiv*nactiv)
      call mxmm(q(nanb+1),nsa4,h,p,nactiv,nactiv,nbasis)
      call mxmb(p,1,nactiv,q,1,nbasis,p(nacnba+1),1,nactiv,
     * nactiv,nbasis,nactiv)
      call symm1(r,p(nacnba+1),nactiv)
      return
    2 continue
      call vclr(p,1,nacnba)
      n=1
      l=1
      do   4 loop=1,nactiv
      m=1
      do   5 moop=1,nbasis
      call daxpy(moop,q(n),h(m),1,p(l),1)
      n=n+1
   5  m=m+moop
   4  l=l+nbasis
      m=1
      n=1
      do 6 loop=1,nactiv
      l=1
      do   7 moop=1,loop
      temp1=ddot(nbasis,q(n),1,p(l),1)
      r(m)=temp1+ddot(nbasis,q(l),1,p(n),1)
      m=m+1
   7  l=l+nbasis
   6  n=n+nbasis
      return
      end
      subroutine mrgle(n,scalar,r)
      implicit REAL  (a-h,o-z), integer (i-n)
      dimension r(*)
      do 1 loop=1,n
      if(dabs(r(loop)).ge.scalar) go to 1
      r(loop)=dsign(scalar,r(loop))
    1 continue
      return
      end
      subroutine mrgge(n,scalar,r)
      implicit REAL  (a-h,o-z), integer (i-n)
      dimension r(*)
      do 1 loop=1,n
      if(dabs(r(loop)).lt.scalar) go to 1
      r(loop)=dsign(scalar,r(loop))
    1 continue
      return
      end
      subroutine stackr(g,nijkl)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
      dimension g(*)
      dimension nijkl(*)
INCLUDE(common/stak)
      common/junk/nwbuck(1500),
     * itx(3400),ktx(3400),gtx(3400)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/bufb/nwbnwb,lnklnk,gout(5118)
INCLUDE(common/three)
      common/scra/ibu(3400),itxktx(3400)
      call jsort1
    1 ibuck=ld340t(g,nijkl)
      if(ibuck)3,2,3
    3 call stopbk
      lnklnk=mark(ibuck)
      ib=ibase(ibuck)+1
      ibn=ibasen(ibuck)+1
      isz=nsz340
      call dcopy(isz,g(ib),1,gout,1)
_IFN1(iv)      call pack(gout(nsz341),lab1632,nijkl(ibn),nsz680)
_IF1(iv)      call fmove(nijkl(ibn),gout(nsz341),nsz170)
      call sttout
      mark(ibuck)=iblock
      iblock=iblock+nsz
      mstack=mstack+1
      if(mstack-nstack)1,1,2
   2  nstack=0
      return
      end
      function isort1(itri,ktri,gtri)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
      common/craypk/integ(1360)
      common/blkin/gin(510),mword
INCLUDE(common/mapper)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ipad
      dimension itri(*),ktri(*),gtri(*)
c
      call unpack(gin(num2e+1),lab816,integ,numlab)
      n=0
      iword=1
      do 1 loop=1,mword
_IF(littleendian)
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      itx=iky(i)+j
      ktx=iky(k)+l
      gtx=gin(loop)
      if(mloww.ge.itx.or.mhi.lt.itx)go to 2
      n=n+1
      itri(n)=itx
      ktri(n)=ktx
      gtri(n)=gtx
    2 if(mloww.ge.ktx.or.mhi.lt.ktx.or.itx.eq.ktx)go to 1
      n=n+1
      itri(n)=ktx
      ktri(n)=itx
      gtri(n)=gtx
    1 iword=iword+4
*     write(6,121)i,j,k,l,itri(n),ktri(n),gtri(n)
*121  format('isort1: i,j,k,l,itri,ktri,G = ', 4i3,1x,2i6,2x,f10.5)
      isort1=n
      return
      end
      function isort2(itri,ktri,gtri)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/blkin/gin(510),nword
INCLUDE(common/atmblk)
      common/junke/maxt,ires,ipas,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ntri
_IF1(iv)      common/craypk/ij205(340),kl205(340),ipad(680)
_IFN1(iv)      common/craypk/ij205(2,340)
      dimension itri(*),ktri(*),gtri(*)
_IF1(iv)      call upak4v(gin(num2ep+1),ij205)
_IFN1(iv)      call unpack(gin(num2ep+1),lab1632,ij205,numlabp)
      n=0
      do 1 loop=1,nword
_IF1(iv)      ij=ij205(loop)
_IFN1(iv)      ij=ij205(1,loop)
      if(mloww.ge.ij.or.mhi.lt.ij)go to 1
      n=n+1
      gtri(n)=gin(loop)
      itri(n)=ij
_IF1(iv)      ktri(n)=kl205(loop)
_IFN1(iv)      ktri(n)=ij205(2,loop)
    1 continue
      isort2=n
      return
      end
      function isort3(imap,itri,ktri,gtri)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension itri(*),ktri(*),gtri(*),imap(*)
      common/blkin/gin(510),nword
      common/junke/maxt,ires,ipas,nteff,npass1,npass2,lentri,
     * nbuck,mloww,mhi,ntri
INCLUDE(common/atmblk)
_IF1(iv)      common/craypk/ij205(340),kl205(340),ipad(680)
_IFN1(iv)      common/craypk/ij205(2,340)
c
_IF1(iv)      call upak4v(gin(num2ep+1),ij205)
_IFN1(iv)      call unpack(gin(num2ep+1),lab1632,ij205,numlabp)
      n=0
      do 1 loop=1,nword
_IF1(iv)      ij=ij205(loop)
_IFN1(iv)      ij=ij205(1,loop)
      itrish=imap(ij)
      if(mloww.ge.itrish.or.mhi.lt.itrish)go to 1
      n=n+1
      gtri(n)=gin(loop)
      itri(n)=ij
_IF1(iv)      ktri(n)=kl205(loop)
_IFN1(iv)      ktri(n)=ij205(2,loop)
    1 continue
      isort3=n
      return
      end
      subroutine jsort1
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/scra/ibuk(3400)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      mstack=1
      do 1 i=1,nstack
    1 ibuk(i)=(itx(i)-mlow)*btri+1
      return
      end
      function ld340t(buf,ijklbuf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF1(iv)      integer *2 ijklbuf
      common/bufb/nkk,mkk,gout(5118)
INCLUDE(common/blksiz)
INCLUDE(common/stak)
c     common/three/mark(1500),lnumb(1500)
INCLUDE(common/three)
      common/scra/ib(3400)
      common/junk/nwbuf(1500),itx(3400),ktx(3400),v(3400)
      dimension buf(*),ijklbuf(*)
      mstak=mstack
    2 ld340t=ib(mstak)
      vv=v(mstak)
      nw=nwbuf(ld340t)+1
      ln=ibase(ld340t)+nw
      buf(ln)=vv
      ijklbuf(ln+ln-1)=itx(mstak)
      ijklbuf(ln+ln  )=ktx(mstak)
*     write(6,121)ln,buf(ln),ijklbuf(ln+ln-1),ijklbuf(ln+ln  )
*121  format('ld340t: ',1x,i5,f15.5,2x,2i10)
      if(nw.ge.nsz340)go to 4
      nwbuf(ld340t)=nw
      mstak=mstak+1
      if(mstak.le.nstack)go to 2
      ld340t=0
      return
    4 nwbuf(ld340t)=0
      mstack=mstak
      return
      end
_IF(ibm,vax)
      subroutine transc(qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer*2 lin
INCLUDE(common/blksiz)
      common/junke/maxt(8),mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
INCLUDE(common/cndx40)
      dimension qq(*)
      dimension lin(2)
      equivalence (lin(1),gin(1))
      ij=nszij
      kl=nszkl
      do 1 loop=1,nkk
      qq((lin(ij)-mloww)*lenb4 + lin(kl)) = gin(loop)
      ij=ij+1
    1 kl=kl+1
      return
      end
_ELSE
      subroutine transc(qq)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/junke/maxt(8),mloww,ipad(2)
      common/bufb/nkk,mkk,gin(5118)
INCLUDE(common/cndx40)
      common/scra/ijklin(2,3412),index(3412)
      dimension qq(*)
      call unpack(gin(nsz341),lab1632,ijklin,nsz680)
      do loop=1,nkk
       index(loop)=(ijklin(1,loop)-mloww)*lenb4 + 
     +              ijklin(2,loop)
      enddo
_IF(cray)
      call scatter(nkk,qq,index,gin)
_ELSE
      call dsctr(nkk,gin,index,qq)
_ENDIF
      return
      end
_ENDIF
_IFN(cray)
      function locate(ref,nref,nw,pat)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(ksr,i8)
      integer ref,pat
_ELSEIF(convex,i8drct)
      integer *8 ref,pat
_ELSE
      REAL  ref,pat
      logical eq
_ENDIF
      dimension ref(nw,*)
      do 1 loop=1,nref
_IF(ksr,convex,i8,i8drct)
      if(ref(1,loop).eq.pat)goto2
_ELSE
      if(eq(ref(1,loop),pat))goto2
_ENDIF
    1 continue
      locate=0
      return
    2 locate=loop
      return
      end
_ENDIF
_IF(ibm,vax)
      subroutine upak8w(ijkl,iiii,mapper)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF1(i)      logical *1 mij,mkl,ii,jj,kk,ll
_IF1(v)      logical *1 mij,mkl
      common/blkin/gout(340),mij(680),mkl(680),nwb
     common/craypk/i205(340),j205(340),k205(340),l205(340)
      dimension mapper(*),ijkl(*),iiii(*)
_IF1(i)      dimension ii(4),jj(4),kk(4),ll(4)
_IF1(i)      equivalence (ii(1),i),(jj(1),j),(kk(1),k),(ll(1),l)
      data i,j,k,l/4*0/
      call igthr(nwb*4,mapper,icray,i205)
      do 1 iw=1,nwb
      il=iw+iw
_IF1(i)      ii(4)=mij(il-1)
_IF1(i)      jj(4)=mij(il  )
_IF1(i)      kk(4)=mkl(il-1)
_IF1(i)      ll(4)=mkl(il  )
_IF1(v)      j=mij(il-1)
_IF1(v)      i=mij(il  )
_IF1(v)      l=mkl(il-1)
_IF1(v)      k=mkl(il  )
      i=mapper(i)
      j=mapper(j)
      k=mapper(k)
      l=mapper(l)
      if(i.ge.j)goto 2
      m=i
      i=j
      j=m
2     if(k.ge.l)goto 3
      m=k
      k=l
      l=m
3      if((i*256+j).ge.(k*256+l))goto 4
      m=i
      i=k
      k=m
      m=j
      j=l
      l=m
   4  i205(iw)=i
      j205(iw)=j
      k205(iw)=k
      l205(iw)=l
 1    continue
      return
      end
_ELSE
      subroutine upak8w(gijkl,iiii,mapper)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/atmblk)
      common/blkin/gout(510),nwb
      common/craypk/i205(4,340)
      dimension icray(1360)
      dimension mapper(*),gijkl(*),iiii(*)
c
      call unpack(gijkl,lab816,i205,numlab)
_IF1(c)      call gather(nwb*4,icray,mapper,i205)
_IFN1(c)      call igthr(nwb*4,mapper,icray,i205)
      int4=1
      do 1 iw=1,nwb
_IF(littleendian)
      j=icray(int4 )
      i=icray(int4+1)
      l=icray(int4+2)
      k=icray(int4+3)
_ELSE
      i=icray(int4 )
      j=icray(int4+1)
      k=icray(int4+2)
      l=icray(int4+3)
_ENDIF
      int4=int4+4
      if(i.ge.j)goto 2
      m=i
      i=j
      j=m
2     if(k.ge.l)goto 3
      m=k
      k=l
      l=m
3      if((i*256+j).ge.(k*256+l))goto 4
      m=i
      i=k
      k=m
      m=j
      j=l
      l=m
_IF(littleendian)
  4   i205(1,iw)=j
      i205(2,iw)=i
      i205(3,iw)=l
      i205(4,iw)=k
_ELSE
  4   i205(1,iw)=i
      i205(2,iw)=j
      i205(3,iw)=k
      i205(4,iw)=l
_ENDIF
 1    continue
      end
_ENDIF
_IF(ibm,vax)
      subroutine upak8s(gijkl,iiii)
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical *1 mij,i205,j205,k205,l205
      common/junk/i205(13648),j205(13648),k205(13648),l205(13648)
      common/bufb/nint,nlink,mij(5118)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      dimension gijkl(*),iiii(*)
c
      nint0=min(nint,nsz170)
_IF1(i)      i4=4
_IF1(i)      left=nsz340*8+5
_IF1(v)      i4=1
_IF1(v)      left=nsz340*8+1
      do 1 i=1,nint0
      i205(i4)=mij(left)
      j205(i4)=mij(left+1)
      k205(i4)=mij(left+2)
      l205(i4)=mij(left+3)
      left=left+8
    1 i4=i4+4
_IF1(i)      left=nsz340*8+1
_IF1(v)      left=nsz340*8+5
      if(nint.le.nsz170)go to 2
      do 3 i=nsz170+1,nint
      i205(i4)=mij(left  )
      j205(i4)=mij(left+1)
      k205(i4)=mij(left+2)
      l205(i4)=mij(left+3)
      i4=i4+4
    3 left=left+8
    2 continue
      return
      end
_ELSE
      subroutine upak8s(gijkl,iiii)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/junk/i205(3412),j205(3412),k205(3412),l205(3412)
      common/junk2/index(4,3412)
      common/bufb/nint,nlink,mij(5118)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      dimension gijkl(*),iiii(*)
c     write(6,*) 'upak8s'
c     j = nsz340/2 + nsz341
c     write(6,99996) (gijkl(i),i=nsz341,j)
c99996 format(1x,4(1x,z16))
_IFN1(ck)      call setsto(13648,0,index)
      call unpack(gijkl(nsz341),lab816,index,nsz340*4)

_IF(littleendian)
_IF(i8)
      do 1 i=1,(nint+1)/2
      i205(2*i-1)=index(2,2*i-1)
      j205(2*i-1)=index(1,2*i-1)
      k205(2*i-1)=index(4,2*i-1)
      l205(2*i-1)=index(3,2*i-1)
      i205(2*i)=index(2,2*i)
      j205(2*i)=index(1,2*i)
      k205(2*i)=index(4,2*i)
   1  l205(2*i)=index(3,2*i)
_ELSE
      if (o255i) then
       do i=1,(nint+1)/2
        i205(2*i-1)=index(2,2*i-1)
        j205(2*i-1)=index(1,2*i-1)
        k205(2*i-1)=index(4,2*i-1)
        l205(2*i-1)=index(3,2*i-1)
        i205(2*i)=index(2,2*i)
        j205(2*i)=index(1,2*i)
        k205(2*i)=index(4,2*i)
        l205(2*i)=index(3,2*i)
       enddo
      else
       do i=1,(nint+1)/2
        i205(2*i-1)=index(2,2*i-1)
        j205(2*i-1)=index(1,2*i-1)
        k205(2*i-1)=index(4,2*i-1)
        l205(2*i-1)=index(3,2*i-1)
        i205(2*i)=index(2,2*i)
        j205(2*i)=index(1,2*i)
        k205(2*i)=index(4,2*i)
        l205(2*i)=index(3,2*i)
       enddo
      endif
c     do  i = 1, nint
c      write(6,99997) nint, i205(i),j205(i),k205(i),l205(i),
c    +  gijkl(i)
c     enddo
99997  format('upak8s: nint,i,j,k,l,g', i6,2x,4i4,2x,f20.10)
_ENDIF
_ELSE
      do 1 i=1,nint
      i205(i)=index(1,i)
      j205(i)=index(2,i)
      k205(i)=index(3,i)
   1  l205(i)=index(4,i)
_ENDIF
      return
      end
_ENDIF
_IF1()      subroutine iscatr(n,ir,map,ia)
_IF1()      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      dimension ia(*),ir(*),map(*)
_IF1()      if(n.le.0)go to 2
_IF1()      do 1 i=1,n
_IF1()    1 ir(map(i))=ia(i)
_IF1()    2 return
_IF1()      end
      subroutine vvdv(n,r,a,b)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      dimension a(*),b(*),r(*)
      do 1 loop=1,n
    1 r(loop)=a(loop)/b(loop)
      return
      end
_IFN(cray)
      subroutine gtrian(n,scalar,r,a,b)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),b(*)
      if(n.le.0)go to 2
      do 1 loop=1,n
    1 r(loop)=scalar*b(loop)-a(loop)
    2 return
      end
_ENDIF
      subroutine wto (buf,lbuf)
      character*(*) buf
      character*80 obuf
      obuf=buf(1:lbuf)
      call sptchk(obuf)
      return
      end
      subroutine fixorb(v,ndim,nrow,ncol)
      implicit none
c
c     This routine canonicalises the sign of the vectors.
c     This is useful in comparing vectors from different runs for
c     debugging purposes.
c     The canonicalisation is done such that the sum of all coefficients
c     in every vector is positive.
c
      integer ndim      ! the leading dimension of the vectors
      integer nrow      ! the length of the vectors
      integer ncol      ! the number of vectors
      REAL v(ndim,ncol) ! the vectors
c
      integer i
      REAL sum, dsum
c
      do i = 1, ncol
         sum = dsum(nrow,v(1,i),1)
         if (sum.lt.0.0d0) then
            call dscal(nrow,-1.0d0,v(1,i),1)
         endif
      enddo
      end
_IFN(ibm,vax,cray,ipsc,t3d)
      subroutine dchksm(n,x,ix,string)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      dimension x(*)
      character*(*) string
c
c real check sum routine for debug
c
      sum = 0.0d0
      asum =0.0d0
      if (n.gt.0) then
_IFN1(civ)         sum = dsum(n,x,ix)
         asum = dasum(n,x,ix)
_IF1()c         write(6,2) (x(i),i=1,n)
_IF1()c2        format(1x,4d19.11)
      endif
      write(iwr,1) string,n,ix,sum,asum
 1    format(1x,a,2i11,2(2x,f23.15))
      return
      end
      subroutine pchksm(n,k,ik,string)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
_IF1(x)      integer*8 k(*),ix,io,ia
_IF1(k)      integer   k(*),ix,io,ia
_IF1(kx)c this equivalence due to bug writing i8 in hex on alliant
_IF1(kx)      equivalence (rx,ix),(ro,io),(ra,ia)
_IFN1(kx)      dimension k(2,*),ix(2),io(2),ia(2)
_IFN1(kx)      dimension iix(2),ioo(2),iaa(2)
      character*(*) string

_IF(apollo)
      data iix / 16#00000000,16#00000000/
      data ioo / 16#00000000,16#00000000/
      data iaa / 16#ffffffff,16#ffffffff/
_ELSEIF(sgi,ipsc)
      data iix / $00000000,$00000000/
      data ioo / $00000000,$00000000/
      data iaa / $ffffffff,$ffffffff/
_ELSEIF(_NOT(linux))
_IF(rs6000,win95)
      data iix / z00000000,z00000000/
      data ioo / z00000000,z00000000/
      data iaa / zffffffff,zffffffff/
_ENDIF
_ELSEIF(alliant,hp700,hpux11,dec)
      data iix / '00000000'x,'00000000'x/
      data ioo / '00000000'x,'00000000'x/
      data iaa / 'ffffffff'x,'ffffffff'x/
_ELSEIF(ksr,convex)
_ELSE
c
c this is the Sun format - replace code as required
c new machines
c
      data iix / z'00000000',z'00000000'/
      data ioo / z'00000000',z'00000000'/
      data iaa / z'ffffffff',z'ffffffff'/
_ENDIF
c
c packed integer check sum routine for debug
c

_IF(ksr,convex)
      ix = '0000000000000000'x
      io = '0000000000000000'x
      ia = 'ffffffffffffffff'x
_ELSE
      do 20 loop=1,2
      ix(loop)=iix(loop)
      io(loop)=ioo(loop)
 20   ia(loop)=iaa(loop)
_ENDIF
      if (n.gt.0) then
         do 10 i = 1,1+(n-1)*ik,ik
_IF(ksr)
            ix = xor(ix,k(i))
            io = or(io,k(i))
            ia = and(ia,k(i))
_ELSEIF(convex)
            ix = ieor(ix,k(i))
            io = ior(io,k(i))
            ia = iand(ia,k(i))
_ELSE
            ix(1) = IXOR32(ix(1),k(1,i))
            ix(2) = IXOR32(ix(2),k(2,i))
            io(1) = IOR32(io(1),k(1,i))
            io(2) = IOR32(io(2),k(2,i))
            ia(1) = IAND32(ia(1),k(1,i))
            ia(2) = IAND32(ia(2),k(2,i))
_ENDIF
 10      continue
c         write(6,2) (k(i),i=1,n)
c2        format(4(1x,z16))
      endif
c
_IF(ksr,convex)
      write(iwr,1) string,n,ik,rx,ro,ra
 1    format(1x,a,2i6,3(2x,z16))
_ELSE
      write(iwr,1) string,n,ik,ix,io,ia
 1    format(1x,a,2i6,3(2x,2z8))
_ENDIF
      return
      end
      subroutine ichksm(n,k,ik,string)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      dimension k(*)
      character*(*) string
c
c integer check sum routine for debug
c
      isum = 0
      iasum = 0
      if (n.gt.0) then
         do 10 i = 1,1+(n-1)*ik,ik
            isum = isum + k(i)
            iasum = iasum + iabs(k(i))
 10      continue
c         write(6,2) (k(i),i=1,n)
c2        format(7(1x,i10))
      endif
c
      write(iwr,1) string,n,ik,isum,iasum
 1    format(1x,a,2i11,2(2x,i11))
      return
      end
_ENDIF
      subroutine wrtsor( text,i205,j205,k205,l205,ggg,nwb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      character *(*) text
      dimension i205(*),j205(*),k205(*),l205(*),ggg(*)
      dimension f(4),ia(4),ja(4),ka(4),la(4)
      write(iwr,1) text,nwb
1     format(/1x,a6,' sorted integrals = ',i5)
      iseq = 0
      n = 0
      do loop=1,nwb
      n = n + 1
      f(n)  = ggg(loop)
      ia(n) = i205(loop)
      ja(n) = j205(loop)
      ka(n) = k205(loop)
      la(n) = l205(loop)
      if(n.eq.4) then
       iseq = iseq + 1
       write(iwr,200)iseq,(ia(k),ja(k),ka(k),la(k),
     +               f(k),k=1,n)
       n = 0
      endif
      enddo
      if(n.ne.0) then
       iseq=iseq+1
       write(iwr,200)iseq,(ia(k),ja(k),ka(k),la(k),
     +               f(k),k=1,n)
      endif
      return
200   format(i5,2x,4(4i4,f11.6))
      end
      subroutine dbgvec(text,array,n)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      character *(*) text
      dimension array(*)
      write(iwr,200)text,n,(array(i),i=1,n)
200   format(1x,a8,i8/(8f12.5))
      return
      end
      subroutine dbgvecv(text,array,n)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      character *(*) text
      dimension array(*)
      dimension f(6),ia(6)
      data thresh /1.0d-5/
c
      write(iwr,200)text,n
      iseq = 0
      m = 0
      do loop=1,n
       if(abs(array(loop)).ge.thresh) then
        m = m + 1
        f(m)  = array(loop)
        ia(m) = loop
        if(m.eq.6) then
        iseq = iseq + 1
        write(iwr,201) iseq,(ia(k),f(k),k=1,m)
        m = 0
        endif
       endif
      enddo
      if(m.ne.0) then
       iseq=iseq+1
       write(iwr,201)iseq,(ia(k),f(k),k=1,m)
      endif
      return
200   format(/1x,'**** ',a10,2x,i8/)
201   format(i5,2x,6(i6,f12.5))
      end
      subroutine wtfort(q,nword,irec,ifort)
c
c     routines for input/output on fortran streams
c     if called with nword.eq.0 no i/o takes place, but the file is
c     positioned ready to read or write at record irec.
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
c
      dimension q(*)
      rewind ifort
      assign 10 to label
      if(irec-1)10,10,30
10    if(nword.gt.0)call wtfors(q,nword,ifort)
      return
c
      entry rdfort(q,nword,irec,ifort)
      rewind ifort
      assign 20 to label
      if(irec-1)20,20,30
20    if(nword.gt.0)call rdfors(q,nword,ifort)
      return
c
30    irecm=irec-1
      do 40 irec=1,irecm
40    read(ifort,err=60,end=50)
      goto label,(10,20)
60    write(iwr,70)ifort,irec
      call caserr('error on attempting to read fortran dataset.')
      return
50    write(iwr,70)ifort,irec
      call caserr('end of file reached on reading fortran dataset.')
      return
70    format(' fault on fortran stream',i3,' at record',i7)
      end
      subroutine wtfors(q,nword,ifort)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(nword)
      if(nword.gt.0)write(ifort)q
      return
      end
      subroutine rdfors(q,nword,ifort)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(nword)
      if(nword.gt.0)read(ifort,end=60,err=50)q
      return
60    call caserr('end of file reached on reading fortran dataset.')
50    call caserr(' error on reading fortran dataset.')
      return
      end
_IF(parallel)
      subroutine broadc(icode,ii,ni)
c
c     broadcast integers from root to others
c     (REALs can be done in obvious way)
c
      implicit REAL  (a-h,o-z), integer (i-n)
      dimension ii(ni)
c
      ni4 = ni * 8 /lenwrd()
      call pg_brdcst(icode,ii,ni4,0)
c
      return
      end
      subroutine comed7
c
c...  switch numci to the number of the common ed7
c...  this file-handling needs a bit mre work
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/disktl)
c 
      numci = numcic
c
      return
c
      entry pried7
c
c..   set numci to private ed7
c
      numci = numcip
c
      return
      end
      subroutine pggop(icode,arr,len,op,buf,lenbuf)
c MPP
c MPP   ggop (do in parts / use buf(lenbuf) as scratch)
c MPP
      implicit REAL  (a-h,o-z), integer (i-n)
      character*1 op
      dimension arr(len),buf(lenbuf)
c MPP
      do 100 i=1,len,lenbuf
         nn = min(lenbuf,len-i+1)
         call pg_dgop(icode,arr(i),nn,op)
100   continue
c MPP
      return
      end
_ENDIF
_IF(fps)
      subroutine pak340(intin,intout)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/blksiz)
      dimension intin(*),intout(*)
      ii=1
      do 1 i=1,nsz85
      intout(i)= insert(
     *           insert(intin(ii+3),32,16,intin(ii+2)),0,32,
     *           insert(intin(ii+1),32,16,intin(ii  ))      )
   1  ii=ii+4
      return
      end
      subroutine upk340(intin,intout)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/blksiz)
      dimension intin(*),intout(*)
      ii=1
      do 1 i=1,nsz85
      intout(ii  )=extract(intin(i), 0,16)
      intout(ii+1)=extract(intin(i),16,16)
      intout(ii+2)=extract(intin(i),32,16)
      intout(ii+3)=extract(intin(i),48,16)
   1  ii=ii+4
      return
      end
      subroutine pakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension intin(*),intout(*)
      ii=1
      do 1 i=1,mword
      intout(i)= insert(
     *           insert(intin(ii+3),32,16,intin(ii+2)),0,32,
     *           insert(intin(ii+1),32,16,intin(ii  ))      )
   1  ii=ii+4
      return
      end
      subroutine upakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension intin(*),intout(*)
      ii=1
      do 1 i=1,mword
      intout(ii  )=extract(intin(i), 0,16)
      intout(ii+1)=extract(intin(i),16,16)
      intout(ii+2)=extract(intin(i),32,16)
      intout(ii+3)=extract(intin(i),48,16)
   1  ii=ii+4
      return
      end
      subroutine stackx(nstak)
      implicit REAL  (a-h,o-z),integer (i-n)
      word i,j,k,l
      common/scra/v(3680),ib(3680),
     * i(3680),j(3680),k(3680),l(3680)
      do 1 n=1,nstak
   1  i(n) = insert(
     *              insert(l(n),48,8,k(n)),32,16,
     *              insert(j(n),48,8,i(n)) )
c
      return
      end
      function ld340(mstack,nstak,nsz170,nwbuf,lnumb,
     * v,ib,i,buf)
      implicit REAL  (a-h,o-z),integer (i-n)
      word buf,ijkl,i
      dimension nwbuf(*),lnumb(*),v(*),ib(*)
      dimension i(*),buf(*)
      mstak=mstack
      nsz340=nsz170+nsz170
    2 ld340=ib(mstak)
      vv=v(mstak)
      ijkl=i(mstak)
      nw=nwbuf(ld340)+1
      ln=lnumb(ld340)+nw
      buf(ln)=vv
      if(nw.gt.nsz170)go to   1
      buf(ln+nsz340)=insert(buf(ln+nsz340),32,32,ijkl)
    3 nwbuf(ld340)=nw
      mstak=mstak+1
      if(mstak.le.nstak)go to   2
      ld340=0
      return
    1 continue
      buf(ln+nsz170)=insert(buf(ln+nsz170),0,32,ijkl)
      if(nw.lt.nsz340)go to   3
      nwbuf(ld340)=0
      mstack=mstak
      return
      end
_ENDIF
_IF(cyber205)
       subroutine pak4f(iout,ij,kl,nw)
cjvl
cjvl      pack ij,kl in iout  / 16 bits/word
cjvl      only called in mcsrto/mcsrti   (jvl 1988)
cjvl
       dimension iout(*),ij(nw),kl(nw)
       data mask /x'ffff'/
c
       nn = nw/2
       do 10 i=1,nn
          iout(i) = or(shift(ij(i*2-1),48),
     1              or(shift(kl(i*2-1),32),
     2              or(shift(ij(i*2),16),kl(i*2))))
10     continue
       if (nw.ne.nn*2) then
          iout(nn+1) = or(shift(ij(nw),48),shift(kl(nw),32))
       end if
c
       return
c
       entry upak4f(iout,ij,kl,nw)
cjvl
cjvl      undo pak4f
cjvl
       nn = nw/2
       do 20 i=1,nn
          ij(i*2-1) = and(shift(iout(i),16),mask)
          kl(i*2-1) = and(shift(iout(i),32),mask)
          ij(i*2) = and(shift(iout(i),48),mask)
          kl(i*2) = and(iout(i),mask)
20     continue
       if (nw.ne.nn*2) then
          ij(nw) = and(shift(iout(nn+1),16),mask)
          kl(nw) = and(shift(iout(nn+1),32),mask)
       end if
cjvl
       return
       end
      subroutine mult3(q,h,r,p,nactiv,nbasis)
      implicit REAL  (a-h,o-z)
      dimension q(*),h(*),r(*),p(*)
INCLUDE(common/mlngth)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
      call vxvxts(nbasis,0.5,h,iky(2),h,iky(2))
      nacnba=nactiv*nbasis
      call szero(p,nacnba)
      if(nbasd2.ge.nactiv)goto 1
      call mxmm(q(nanb+1),nsa4,h,p,nactiv,nactiv,nbasis)
      call mymo(p,nactiv,q,1,nbasis,p(nacnba+1),nactiv,
     *nactiv,nbasis,nactiv)
    2 call symm1(r,p(nacnba+1),nactiv)
      return
    1 call mxmtr(h,q,p,nbasis,nactiv)
      call mymi(q,nbasis,p,nbasis,p(nacnba+1),1,nactiv,
     * nactiv,nbasis,nactiv)
      goto 2
      end
      subroutine martyn(text)
c..
c..    debug-help routine for cyber205
c..    print message text and flush the output buffer
c..
       character*(*) text
INCLUDE(common/iofile)
c
       write(iwr,*) text
       endfile iwr
       backspace iwr
c
       return
       end
      subroutine pakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension intin(*),intout(*)
      ii=1
      do 1 i=1,mword
      intout(i)=or(shift(intin(ii),48),
     *             or(shift(intin(ii+1),32),
     *                or(shift(intin(ii+2),16),
     *                      intin(ii+3))))
   1  ii=ii+4
      return
      end
      subroutine upakin(intin,intout,mword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension intin(*),intout(*)
      data mask/x'ffff'/
      do 1 i=1,mword
       intout(i*4-3)=and(shift(intin(i),16),mask)
       intout(i*4-2)=and(shift(intin(i),32),mask)
       intout(i*4-1)=and(shift(intin(i),48),mask)
    1  intout(i*4  )=and(      intin(i)    ,mask)
      return
      end
      subroutine sqsitr(r,s,t,nk)
      implicit REAL  (a-h,o-z), integer  (i-n)
      dimension s(*),t(*)
      dimension r(*)
      common/scr205/temp(1)
      call square(r,s,nk,nk)
      call sqtrip(temp,t,nk)
      call vadd(r,1,temp,1,r,1,nk*nk)
      return
      end
      subroutine symm1a(r,a,nk)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/scr205/temp(1)
INCLUDE(common/helpr)
      dimension r(*)
      dimension a(*)
      call symm1(temp,a,nk)
      call vadd(r,1,temp,1,r,1,iky(nk+1))
      return
      end
      subroutine anti1s(r,a,nk)
      implicit REAL  (a-h,o-z),integer (i-n)
      common/scr205/temp(1)
INCLUDE(common/helpr)
      dimension r(*)
      dimension a(*)
      call anti1(temp,a,nk)
      call subvec(r,r,temp,iky(nk))
      return
      end
      subroutine anti1a(r,a,nk)
      implicit REAL  (a-h,o-z),integer  (i-n)
      common/scr205/temp(1)
INCLUDE(common/helpr)
      dimension r(*)
      dimension a(*)
      call anti1(temp,a,nk)
      call vadd(r,1,temp,1,r,1,iky(nk))
      return
      end
_ENDIF
      subroutine ver_util5(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util5.m,v $
     +     "/
      data revision /"$Revision: 6090 $"/
      data date /"$Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
