_IF(hpux11)
c$HP$ OPTIMIZE ASSUME_NO_PARAMETERS_OVERLAPS OFF
_ENDIF
c 
c  $Author: mrdj $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/casb.m,v $
c  $State: Exp $
c  
c     deck=casb
c ******************************************************
c ******************************************************
c             =   casguga    =
c ******************************************************
c ******************************************************
      subroutine guga(q,lword)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
      common/drtlnk/norbl(4),labdrt(maxorb),offmt0(2)
INCLUDE(common/machin)
      dimension q(*)
      character*10 charwall
c
      write(iwr,2)
2     format(/1x,104('=')//
     *' distinct row table generator'/
     *1x,28('=')/)
      call drt(q,lword)
      write(iwr,3)
 3    format(/1x,104('=')//
     * ' loop-driven formula tape generator'/
     *1x,34('=')/)
      call bugme(q,lword)
      if(offmt0(2))go to 5
      write(iwr,4)
 4    format(/1x,104('=')//
     *' formula tape construction and sorting'/
     *1x,37('=')/)
      call citape(q,lword)
c
5     call clredx
c
      call revise
      call whtps
c
      top=cpulft(1)
      write(iwr,6)top,charwall()
 6    format(/
     *' end of guga-ci module at ',f8.2,' seconds',a10,' wall')
      return
      end
      subroutine bugme(ia,maxx)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/machin)
INCLUDE(common/iofile)
      common/lsort /norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb),
     1        levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),levfrm
     2 , jspace(512)
c  option 1 - printing option
c  option 2 - full spin space
c  option 3 - use additional tape4 for small matrix
c  option 4 - excitation option
c  option 5 - do only configuration generation
c  option 6 - use smaller output buffer size
c  option 7 -
      common/loops /ix(6),ax(4),ix2(6),intsad,nmax,acrcy,ngrps
INCLUDE(common/files)
INCLUDE(common/restar)
      character*10 charwall
      dimension ia(*)
c
      itape8 = ntapes(1)
      mfilep=1
      mainp=n9tape(1)
      iblkmp=n9blk(1)
      mblp=iblkmp-n9last(1)
      nav = lenwrd()
      max=maxx*nav
       write(iwr,60)max
60     format(/1x,'store available is',i8,' integer words')
      rewind itape8
      do 80 i=1,14
80    read(itape8)
      read(itape8)ngrps,nxt,nmax
       rewind itape8
      read(itape8) nbf,norbs,nsym,nrows,nwks,levfrm,idum,iop
      write(iwr,100)
100   format(/     ' read in distinct row table')
      write(iwr,120)nwks
      write(iwr,130)nbf
      write(iwr,140)norbs
120   format(/' the number of configurations is',i6)
130   format(' the number of basis functions is',i5)
140   format(' the number of allowed orbitals is',i4)
      ih2=4*nrows+1
      ih3=ih2+nrows
      ih4=ih3+nrows
      ih5=1
      ih6=ih4+nrows
      ih7=ih6+nrows
      ih8=ih7+nrows
      ih11=ih8+nrows
      ih12=ih11+nrows*4
      nkl=nsym*norbs
      nij=(norbs*(norbs+1))/2
      ih13=ih12+nij
      ih14=ih13+nij
      ih15=ih14+nkl
      ih16=ih15+nkl
      ih17=ih16+200
      ih18=ih17+200
      ih19=ih18+200
      lvfrm=levfrm
      if(lvfrm.lt.0)lvfrm=0
      nlng=lvfrm*nsym
      ih20=ih19+nlng
      ih21=ih20+nlng
      ih22=ih21+nlng
      ih23=ih22+nlng
      ih26=ih23+nlng
      ih9=ih26+1
      ih27=ih9+nwks
      len=max-ih27
      intsad=ih27
      if((nmax).gt.(len-3))goto 200
      nrows4=nrows*4
      call getdrt(itape8,nrows,nrows4,nkl,nij,ia(ih2),ia(ih3),ia(ih4),
     1ia(ih5),ia(ih6),ia(ih7),ia(ih8),ia(ih9),ia(ih11),ia(ih12),
     2ia(ih13),ia(ih14),ia(ih15),ia(ih16),ia(ih17),ia(ih18),nwks,
     3nxt,ngrps)
      ifront=(nwks-1)/(511*nav)+1+lensec(mach(17))
      iblkmp=iblkmp+ifront
      mblp=mblp+ifront
      call initmx(nwks)
      call initcs(ia(1))
      if(lvfrm.le.0)goto 160
      call initex(ia,ia(ih3),ia(ih6),ia(ih7),ia(ih8),ia(ih19)
     1,ia(ih20),ia(ih21),ia(ih22),ia(ih23),ia(ih12),ia(ih14),ia(ih15)
     2,ia(ih2),ia(ih4),nrows,nij,nkl,nlng)
160   continue
       write(iwr,170)
170    format( /' generate loops using the loop-driven algorithm')
      call loopy(ia,nrows,ia(ih3),ia(ih6),ia(ih7),ia(ih8),ia(ih11)
     1,ia(ih2),ia(ih12),ia(ih13),ia(ih14),ia(ih15),ia(ih16),ia(ih17)
     2,ia(ih18),nxt,ngrps,nrows4,nij,nkl,ia(ih19),ia(ih20),
     3ia(ih21),ia(ih22),ia(ih23),ia(ih4),nlng)
      call casout(nbf,nwks,ia(ih9),itape8)
      call clredx
      top=cpulft(1)
      write(iwr,180)top,charwall()
 180  format(/
     *' end of loop formula generation at ',f8.2,' seconds',a10,' wall')
      return
200   continue
      write(iwr,210)nmax,len
210   format(' nmax=',i10,' len=',i10)
      call caserr('problem with nmax')
      write(iwr,230)ih27,len
230   format(' ih27 len ',2i10)
      call caserr('not enough memory available in guga')
c
      return
      end
      subroutine getdrt(itape8,nrows,nrows4,nkl,nij,nabca,nabcb,nabcs,
     1iarc,nlwks,nuwks,ipuwk,indx,iwght,ijadd,ijgrp,kadd,ladd,inext,
     2jmnnxt,jmxnxt,nwks,nxt,ngrps)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),levfrm
      dimension nabca(nrows),nabcb(nrows),nabcs(nrows)
      dimension iarc(nrows4),nlwks(nrows),nuwks(nrows),ipuwk(nrows)
      dimension indx(nwks),iwght(nrows4),ijadd(nij),ijgrp(nij)
      dimension kadd(nkl),ladd(nkl),inext(200),jmnnxt(200),jmxnxt(200)
c
      read(itape8)iout
      read(itape8)isym
      read(itape8)levnr
      read(itape8)levpt
      read(itape8)nabca
      read(itape8)nabcb
      read(itape8)nabcs
      read(itape8)nlwks
      read(itape8)nuwks
      read(itape8)ipuwk
      read(itape8)indx
      read(itape8)iarc
      read(itape8)iwght
      read(itape8)ngrps,nxt,idumf
      read(itape8)ijadd
      read(itape8)ijgrp
      read(itape8)kadd
      read(itape8)ladd
      read(itape8)inext
      read(itape8)jmnnxt
      read(itape8)jmxnxt
c
      return
      end
      subroutine initcs(ints)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      common/loops /lad,nuwk,nlwk,iuwk,juwk,itrack,a,c,d,val,ix(6),
     >                intsad,nmax,acrcy,nblks
      common/form  /nform,int1,int2,int3,icoup1,icoup2,icoup3,icoupz
      common/casscg/sqrt2,nuwkw(81),
     * iblk,nlext,nlint,lab1,lab2,lab3,lab4
      dimension ints(*)
c
      itape8 = ntapes(1)
      call isurd1(1.0d0,-1.0d0)
      sqrt2= dsqrt(2.0d0)
      third=1.0d0/3.0d0
      sqrt34= dsqrt(0.75d0)
      sqrth= dsqrt(0.5d0)
      lab1=isurd(1.0d0)
      lab2=isurd(-1.0d0)
      idumf=isurd(-sqrth)
      idumf=isurd(sqrth)
      idumf=isurd(-0.5d0)
      idumf=isurd(-sqrt34)
      idumf=isurd(0.5d0)
      idumf=isurd(third)
      idumf=isurd(-third)
      idumf=isurd(sqrt34)
      lab3=isurd(sqrt2)
      lab4=isurd(2.0d0)
      icoupz=isurd(0.0d0)
      iblk=1
      nlext=0
      nlint=0
      imn=intsad+2
      write(iwr,6532)nmax,imn,intsad
6532  format(/'  nmax  ',i10,' imn',i10,' intsad',i10)
      imnmax=imn+nmax
      imn=imn+1
      read(itape8)
      read(itape8)(ints(i),i=imn,imnmax)
c      write(iwr,6534)(ints(imn+i),i=1,nmax)
c  nowall integrals of iblk are in core
      return
      end
      subroutine nxtblk(ints)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/iofile)
      common/loops /lad,nuwk,nlwk,iuwk,juwk,itrack,a,c,d,val,ix(6),
     >                intsad,nmax,acrcy,nblks
      common/form  /nform,int1,int2,int3,icoup1,icoup2,icoup3,icoupz
      common/casscg/sqrt2,nuwkw(81),
     * iblk,nlext,nlint,lab1,lab2,lab3,lab4
c
      dimension ints(*)
c
      itape8 = ntapes(1)
      if(iblk.gt.nblks)goto 910
      iblk=iblk+1
      imn=intsad+2
      write(iwr,6532)nmax,imn,intsad
6532  format(/'  nmax  ',i10,' imn',i10,' intsad',i10)
      imnmax=imn+nmax
      imn=imn+1
      read(itape8)
      read(itape8)(ints(i),i=imn,imnmax)
      return
 910  write(iwr,915)iblk,nblks
 915  format('  iblk= ',i6,' nblks= ',i6)
      call caserr('trouble with blocks in nxtblk')
      return
      end
      subroutine  putext(ints)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      common/loops /lad,nuwk,nlwk,iuwk,juwk,itrack,a,c,d,val,ix(6),
     >                intsad,nmax,acrcy,nblks
      common/form  /nform,int1,int2,int3,icoup1,icoup2,icoup3,icoupz
      common/casscg/sqrt2,nuwkw(81),
     * iblk,nlext,nlint,lab1,lab2,lab3,lab4
      dimension ints(*)
c     itape8 = ntapes(1)
      nlext=nlext+1
      if(nuwk.eq.0)return
      goto (261,262,263,264,265,266,267,268,269,270,271,272,273,274,
     *275,276,277,278),itrack
      go to 400
261   nform=1
      int1=ints(lad+1)
      icoup1=isurd(a)
      go to 300
262   nform=1
      int1=ints(lad+2)
      icoup1=isurd(a)
      go to 300
263   nform=1
      int1=ints(lad+3)
      icoup1=isurd(a)
      go to 300
264   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=lab1
      icoup2=lab1
      go to 300
265   nform=2
      int1=ints(lad+2)
      int2=ints(lad+1)
      icoup1=lab4
      icoup2=lab1
      go to 300
266   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=lab1
      icoup2=lab1
      go to 300
267   nform=2
      int1=ints(lad+3)
      int2=ints(lad+1)
      icoup1=isurd(a)
      icoup2=icoup1
      go to 300
268   nform=3
      int1=ints(lad+3)
      int2=ints(lad+2)
      int3=ints(lad+1)
      icoup1=isurd(a)
      icoup2=icoup1
      icoup3=icoup1
      go to 300
269   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
270   nform=2
      int1=ints(lad+2)
      int2=ints(lad+1)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
271   nform=2
      int1=ints(lad+1)
      int2=ints(lad+3)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
272   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
273   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=icoup1
      go to 300
274   nform=2
      int1=ints(lad+1)
      int2=ints(lad+3)
      icoup1=lab1
      icoup2=lab1
      go to 300
275   nform=2
      int1=ints(lad+2)
      int2=ints(lad+3)
      icoup1=lab1
      icoup2=lab2
      goto 300
276   nform=1
      int1=ints(lad+1)
      icoup1=lab3
      goto 300
277   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=icoup1
      goto 300
278   nform=2
      int1=ints(lad+2)
      int2=ints(lad+3)
      icoup1=isurd(a)
      icoup2=isurd(-a)
 300  continue
      call outmxe
      return
400   call caserr('computed goto out of range in putout')
c
      return
      end
      subroutine finout
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      common/casscg/sqrt2,nuwkw(81),
     * iblk,nlext,nlint,lab1,lab2,lab3,lab4
c
      nltot=nlint+nlext
      write(iwr,404)nltot
 404  format(/,i12,' is the total number of generated loops')
      write(iwr,405)nlint
 405  format(i12,' were created by the loop-driven algorithm')
      write(iwr,406)nlext
 406  format(i12,' were created implicitly')
      call finmxe
c
      return
      end
      subroutine  putout(ints)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/loops /lad,nuwk,nlwk,iuwk,juwk,itrack,a,c,d,val,ix(6),
     >                intsad,nmax,acrcy,nblks
      common/form  /nform,int1,int2,int3,icoup1,icoup2,icoup3,icoupz
      common/casscg/sqrt2,nuwkw(81),
     * iblk,nlext,nlint,lab1,lab2,lab3,lab4
      dimension ints(*)
c
      nlint=nlint+1
      goto (261,262,263,264,265,266,267,268,269,270,271,272,273,274,
     *275,276,277,278),itrack
      go to 400
261   nform=1
      int1=ints(lad+1)
      icoup1=isurd(a)
      go to 300
262   nform=1
      int1=ints(lad+2)
      icoup1=isurd(a)
      go to 300
263   nform=1
      int1=ints(lad+3)
      icoup1=isurd(a)
      go to 300
264   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=lab1
      icoup2=lab1
      go to 300
265   nform=2
      int1=ints(lad+2)
      int2=ints(lad+1)
      icoup1=lab4
      icoup2=lab1
      go to 300
266   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=lab1
      icoup2=lab1
      go to 300
267   nform=2
      int1=ints(lad+3)
      int2=ints(lad+1)
      icoup1=isurd(a)
      icoup2=icoup1
      go to 300
268   nform=3
      int1=ints(lad+3)
      int2=ints(lad+2)
      int3=ints(lad+1)
      icoup1=isurd(a)
      icoup2=icoup1
      icoup3=icoup1
      go to 300
269   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
270   nform=2
      int1=ints(lad+2)
      int2=ints(lad+1)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
271   nform=2
      int1=ints(lad+1)
      int2=ints(lad+3)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
272   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=isurd(a*d)
      go to 300
273   nform=2
      int1=ints(lad+3)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=icoup1
      go to 300
274   nform=2
      int1=ints(lad+1)
      int2=ints(lad+3)
      icoup1=lab1
      icoup2=lab1
      go to 300
275   nform=2
      int1=ints(lad+2)
      int2=ints(lad+3)
      icoup1=lab1
      icoup2=lab2
      goto 300
276   nform=1
      int1=ints(lad+1)
      icoup1=lab3
      goto 300
277   nform=2
      int1=ints(lad+1)
      int2=ints(lad+2)
      icoup1=isurd(a)
      icoup2=icoup1
      goto 300
278   nform=2
      int1=ints(lad+2)
      int2=ints(lad+3)
      icoup1=isurd(a)
      icoup2=isurd(-a)
 300  continue
      call outmxe
      return
400   call caserr('computed goto out of range in putout')
c
      return
      end
      subroutine initmx(nwk)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junk/
     *nuwkd(340),nlwkd(340),iuwkd(340),juwkd(340),intd(340),icoupd(340),
     *nelind,icoded,
     2nuwko(340),nlwko(340),iuwko(340),juwko(340),into(340),icoupo(340)
     3,nelino,icodeo,  ispr1,ispr2
      common/craypk/ijkl(1)
      common/form  /nform,iint(3),icoup(3),icoupz
      common/loops /lad,nuwk,nlwk,iuwk,juwk,ix,ax(3),val,ilft(8),dumy
     *,imj
      common/casscg/sqrt2,nuwkw(88),
     * nwks,isize,neld,nelo,nbufd,nbufo,nwor
c
      nwor=2042/lenwrd()
      call setsto(2044,0,ijkl)
      call setsto(2042,0,nuwkd)
      call setsto(2042,0,nuwko)
      ispr1=0
      ispr2=0
      icoded=0
      icodeo=1
      nelind=0
      nelino=0
      nwks=nwk
      isize=340
      neld =0
      nelo =0
      nbufd=0
      nbufo=0
c
      return
      end
      subroutine outmxe
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junk/
     *nuwkd(340),nlwkd(340),iuwkd(340),juwkd(340),intd(340),icoupd(340),
     *nelind,icoded,
     2nuwko(340),nlwko(340),iuwko(340),juwko(340),into(340),icoupo(340)
     3,nelino,icodeo,  ispr1,ispr2
      common/craypk/ijkl(1)
      common/form  /nform,integ(3),icoup(3),icoupz
      common/loops /lad,nuwk,nlwk,iuwk,juwk,ix,ax(3),val,ilft(8),dumy
     *,imj
      common/casscg/sqrt2,nuwkw(88),
     * nwks,isize,neld,nelo,nbufd,nbufo,nwor
c
_IFN1(c)       nav = lenwrd()
      do 30 iform=1,nform
      if(icoup(iform).eq.icoupz)goto30
      if(iuwk.eq.juwk) go to 20
c off-diagonal tape construction
      nelo=nelo+1
      nelino=nelino+1
      nuwko (nelino)=nuwk
      nlwko (nelino)=nlwk
      iuwko (nelino)=iuwk
      juwko (nelino)=juwk
      into  (nelino)=integ(iform)
      icoupo(nelino)=icoup(iform)
      if(nelino.lt.isize)go to 30
      nbufo=nbufo+1
_IF1(c)      call fmove(nuwko,ijkl,nwor)
_IFN1(c)      call icopy(nwor*nav,nuwko,1,ijkl,1)
      call ed0out
      nelino=0
      goto 30
c
20    continue
c diagonal tape construction
      nelind=nelind+1
      nuwkd (nelind)=nuwk
      nlwkd (nelind)=nlwk
      iuwkd (nelind)=iuwk
      juwkd (nelind)=juwk
      intd  (nelind)=integ(iform)
      icoupd(nelind)=icoup(iform)
      neld=neld+1
      if(nelind.lt.isize) goto 30
      nbufd=nbufd+1
_IF1(c)      call fmove(nuwkd,ijkl,nwor)
_IFN1(c)      call icopy(nwor*nav,nuwkd,1,ijkl,1)
      call ed0out
      nelind=0
30    continue
c
      return
      end
      subroutine finmxe
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/junk/
     *nuwkd(340),nlwkd(340),iuwkd(340),juwkd(340),intd(340),icoupd(340),
     *nelind,icoded,
     2nuwko(340),nlwko(340),iuwko(340),juwko(340),into(340),icoupo(340)
     3,nelino,icodeo,  ispr1,ispr2
INCLUDE(common/iofile)
      common/blkin /dum(511)
      common/craypk/ijkl(1)
      common/form  /nform,iint(3),icoup(3),icoupz
INCLUDE(common/files)
INCLUDE(common/restar)
      common/loops /lad,nuwk,nlwk,iuwk,juwk,ix,ax(3),val,ilft(8),dumy
     *,imj
      common/casscg/sqrt2,nuwkw(88),
     * nwks,isize,neld,nelo,nbufd,nbufo,nwor
      data m0/0/
_IFN1(c)      nav = lenwrd()
      if(nelino.le.0)goto 40
_IF1(c)      call fmove(nuwko,ijkl,nwor)
_IFN1(c)      call icopy(nwor*nav,nuwko,1,ijkl,1)
      call ed0out
      nbufo=nbufo+1
40    if(nelind.le.0)goto 50
_IF1(c)      call fmove(nuwkd,ijkl ,nwor)
_IFN1(c)      call icopy(nwor*nav,nuwkd,1,ijkl,1)
      call ed0out
      nbufd=nbufd+1
50    call put(dum,m0,n9tape(1))
      neltot=neld+nelo
      write(iwr,60)neltot
60    format(/,i12,' is the total number of processed loop formulae')
      write(iwr,70)neld,nbufd
70    format(i12,'     diagonal loop formulae are stored in',
     1  i6, ' blocks')
      write(iwr,80)nelo,nbufo
80    format(i12,' off diagonal loop formulae are stored in',
     1  i6, ' blocks')
      m9file=mfilep
      m9tape(m9file)=n9tape(mfilep)
      m9blk(m9file) =n9blk(  mfilep)
      m9last(m9file)=iblkmp+1
      write(iwr,10)
 10   format(/' status of loop formula tape'/1x,27('*')/)
      call filprn(m9file,m9blk,m9last,m9tape)
c
      return
      end
      function isurd(val)
c
c     locate the index of coupling coefficient val by searching a balanc
c     binary tree. if unsuccessful, tree is updated, and re-balanced if
c     necessary. (see knuth, vol3, p451).                 pjk/15-iii-82
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      integer    mxsurd
      parameter (mxsurd=3600)
      common/three /surd(mxsurd),nsurd ,link(mxsurd+1,2),ibal(mxsurd+1)
c
      data small,dzero/1.0d-6,0.0d0/
c
      linki=link(mxsurd+1,2)
      is=linki
      it=mxsurd+1
10    isurd=linki
      test=val-surd(isurd)
      if( dabs(test).lt.small)return
      lr=1
      if(test.gt.dzero)lr=2
      linki=link(isurd,lr)
      if(linki.eq.0)goto 20
      if(ibal(linki).eq.0)goto 10
      it=isurd
      is=linki
      goto 10
20    nsurd=nsurd+1
      if(nsurd.gt.mxsurd) then
         write (iwr,'(//,''nsurd = '', i10,5x,''mxsurd = '',i10)')
     +          nsurd,mxsurd
         call caserr('table of coefficients too big')
      end if
      link(isurd,lr)=nsurd
      surd(nsurd)=val
      link(nsurd,1)=0
      link(nsurd,2)=0
      ibal(nsurd)=0
c.....adjust balance factors
      lr=1
      if(val.gt.surd(is))lr=2
      ir=link(is,lr)
      isurd=ir
30    if(isurd.eq.nsurd)goto 40
      lr=1
      if(val.gt.surd(isurd))lr=2
      ibal(isurd)=lr
      isurd=link(isurd,lr)
      goto 30
40    lr=1
      if(val.gt.surd(is))lr=2
      ibals=ibal(is)
      if(ibals.eq.lr)goto 50
      ibal(is)=0
      if(ibals.ne.0)return
      ibal(is)=lr
      link(mxsurd+1,1)=link(mxsurd+1,1)+1
      return
c.....rebalance the tree
50    if(ibal(ir).ne.lr)goto 60
c.... single rotation
      isurd=ir
      link(is,lr)=link(ir,3-lr)
      link(ir,3-lr)=is
      ibal(is)=0
      ibal(ir)=0
      goto 70
c.....double rotation
60    isurd=link(ir,3-lr)
      link(ir,3-lr)=link(isurd,lr)
      link(isurd,lr)=ir
      link(is,lr)=link(isurd,3-lr)
      link(isurd,3-lr)=is
      ibalss=ibal(isurd)
      ibal(is)=0
      ibal(ir)=0
      if(ibalss.eq.lr)ibal(is)=3-lr
      if(ibalss.eq.3-lr)ibal(ir)=lr
      ibal(isurd)=0
70    lr=1
      if(is.eq.link(it,2))lr=2
      link(it,lr)=isurd
      isurd=nsurd
      return
      end
      subroutine isurd1(val1st,val2nd)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer    mxsurd
      parameter (mxsurd=3600)
      common/three /surd(mxsurd),nsurd ,link(mxsurd+1,2),ibal(mxsurd+1)
      val1=val1st
      val2=val2nd
      if(val1.lt.val2)goto 80
      val1=val2nd
      val2=val1st
80    link(1,1)=0
      link(1,2)=2
      link(2,1)=0
      link(2,2)=0
      ibal(1)=2
      ibal(2)=0
      link(mxsurd+1,1)=2
      link(mxsurd+1,2)=1
      surd(1)=val1
      surd(2)=val2
      nsurd=2
c
      return
      end
      subroutine initex(iaa,nabcb,nlwks,nuwks,ipuwk,nuwky,nlwkx,
     *ipuwky,nlwkw,nlwky,ijadd,kadd,ladd,nabca,nabcs,nrows,nij,
     *nkl,nlng)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/cassav/mfrm
      common/lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     *,levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),levfrm
      common/loops/lad,nuwk,nlwk,iuwk,juwk,itrack,acoef,ccoef,
     *dcoef,val,levi,jmin,jmax,iadt,jsmt,ksb,intsad,nmax,acrcy,nblks
      common/cone/kmna(9),lkupsm(64)
      common/casscg/sqrt2,nuwkw(8),nuwkx(8),ismoff(8),ipuwkw(8),
     * ipuwkx(8),msmoff(8),mfrmpt(25),mults(8)
      dimension iaa(*)
      dimension nlwkw(nlng),nlwkx(nlng),nlwky(nlng),nuwky(nlng),
     *ipuwky(nlng)
      dimension nabcb(nrows),nlwks(nrows),nuwks(nrows),ipuwk(nrows),
     *ijadd(nij),kadd(nkl),ladd(nkl),nabca(nrows),nabcs(nrows)
      ism=0
      do 11 i=1,8
      mults(i)=ism
 11   ism=ism+8
      do 20 ism=1,nsym
      ismoff(ism)=(ism-1)*levfrm
       msmoff(ism)=(ism-1)*norbs
      mfrm=levpt(levfrm)
      do 10 i=1,levfrm
      ipt=ismoff(ism)+i
      nlwkw(ipt)=0
      nlwkx(ipt)=0
      nlwky(ipt)=0
      nuwky(ipt)=0
      ipuwky(ipt)=0
10    continue
      nuwkw(ism)=0
      nuwkx(ism)=0
      ipuwkw(ism)=0
      ipuwkx(ism)=0
20    continue
      do 30 im=1,25
30    mfrmpt(im)=0
      do 100 lev=2,levfrm
c     levm=lev-1
      npt=levpt(lev)
      nr=levnr(lev)
      do 90 ir=1,nr
      npt=npt+1
      ia=nabca(npt)
      if(ia.ge.2)goto 930
      ib=nabcb(npt)
      ism=nabcs(npt)
      ismpt=ismoff(ism)+lev
      if(ia.eq.1)goto 70
      if(ib.eq.2)goto 80
      if(ib.eq.0)goto 90
c  this is a y point
      nuwky(ismpt)=nuwks(npt)
      nlwky(ismpt)=nlwks(npt)
      ipuwky(ismpt)=ipuwk(npt)
      if(lev.lt.levfrm)goto 90
      goto 90
70    continue
c  this is a w point
      nlwkw(ismpt)=nlwks(npt)
      if(lev.lt.levfrm)goto 90
      nuwkw(ism)=nuwks(npt)
      ipuwkw(ism)=ipuwk(npt)
      mfrmpt(ir)=ism
      goto 90
80    continue
c  this is a x point
      nlwkx(ismpt)=nlwks(npt)
      if(lev.lt.levfrm)goto 90
      nuwkx(ism)=nuwks(npt)
      ipuwkx(ism)=ipuwk(npt)
      mfrmpt(ir)=-ism
90    continue
100   continue
      sqrt2= dsqrt(2.0d0)
      return
      entry extern(iaa,nabcb,nlwks,nuwks,ipuwk,nuwky,nlwkx,
     *ipuwky,nlwkw,nlwky,ijadd,kadd,ladd,nabca,nabcs,nrows,nij,
     *nkl,nlng)
      i=levi-1
      nlwk=1
      iad=(i*(i-1))/2
      jsm=isym(i)
c     isoff=ismoff(jsm)
      lkupj=mults(jsm)
      jkind=1
      ccoef=1.0d0
      dcoef=-1.0d0
      nuwk5w=nuwkw(1)
      juwk5w=ipuwkw(1)
      iuwk5w=juwk5w+nlwkw(i)
      do 800 j=jmin,jmax
      ij=iad+j
      jad=ijadd(ij)
      isj=isym(j)
      lkupsj=mults(isj)
      lkup=lkupj+isj
      ihsm3=lkupsm(lkup)
      lkupk=mults(ihsm3)
      ksmptx=msmoff(ihsm3)
      jsoff=ismoff(isj)
      ih3off=ismoff(ihsm3)
      nuwk3w=nuwkw(ihsm3)
      nuwk3x=nuwkx(ihsm3)
      ipt1=ih3off+i
      ipt2=jsoff+j
      kuwk3w=ipuwkw(ihsm3)
      kuwk3x=ipuwkx(ihsm3)
      iuwk3w=kuwk3w+nlwkw(ipt1)+nlwky(ipt2)
      iuwk3x=kuwk3x+nlwkx(ipt1)+nlwky(ipt2)
      if(i.ne.j)goto 160
c  cases where i=j
      jkind=2
      iuwk9w=juwk5w+nlwkw(i+1)-1
160   continue
      kkind=jkind
      do 750 k=1,j
      ksmpt=k+ksmptx
      kad=jad+kadd(ksmpt)
      isk=isym(k)
      lkup=lkupk+isk
      lsm=lkupsm(lkup)
      lkupl=mults(lsm)
      lsoff=ismoff(lsm)
      lsmptx=msmoff(lsm)
      ksoff=ismoff(isk)
      acoef=1.0d0
      ipt1=ih3off+k
      juwk3w=kuwk3w+nlwkw(ipt1)
      juwk3x=kuwk3x+nlwkx(ipt1)
      if(k.eq.j)goto 170
      lkup=lkupj+isk
      ihsm1=lkupsm(lkup)
      ih1off=ismoff(ihsm1)
      nuwk1w=nuwkw(ihsm1)
      nuwk1x=nuwkx(ihsm1)
      kuwk1w=ipuwkw(ihsm1)
      kuwk1x=ipuwkx(ihsm1)
      ipt1=ih1off+i
      ipt2=ksoff+k
      iuwk1w=kuwk1w+nlwkw(ipt1)+nlwky(ipt2)
      iuwk1x=kuwk1x+nlwkx(ipt1)+nlwky(ipt2)
      ipt1=ih1off+j
      juwk1w=kuwk1w+nlwkw(ipt1)
      juwk1x=kuwk1x+nlwkx(ipt1)
      lkup=lkupsj+isk
      ihsm2=lkupsm(lkup)
      ih2off=ismoff(ihsm2)
      nuwk2w=nuwkw(ihsm2)
      nuwk2x=nuwkx(ihsm2)
      kuwk2w=ipuwkw(ihsm2)
      kuwk2x=ipuwkx(ihsm2)
      ipt1=ih2off+i
      iuwk2w=kuwk2w+nlwkw(ipt1)
      iuwk2x=kuwk2x+nlwkx(ipt1)
      ipt1=ih2off+j
      juwk2w=kuwk2w+nlwkw(ipt1)+nlwky(ipt2)
      juwk2x=kuwk2x+nlwkx(ipt1)+nlwky(ipt2)
      goto 180
c  cases where j=k
170   continue
      kkind=kkind+2
      juwk7w=juwk5w+nlwkw(k+1)-1
      juwk9w=juwk5w+nlwkw(k)
180   continue
      lkind=kkind
      do 700 l=1,k
      if(isym(l).ne.lsm)goto 700
      lsmpt=l+lsmptx
      lad=kad+ladd(lsmpt)
      ihl=lsoff+l
      nlwkyl=nlwky(ihl)
      if(l.eq.k)lkind=lkind+4
      goto (250,300,350,400,500,550,700,600),lkind
250   continue
c  cases where all indices are different
c  type 1w
      nuwk=nuwk1w
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=6
      call putext(iaa)
c  type 1x
      nuwk=nuwk1x
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack =15
      call putext(iaa)
c  type 2w
      nuwk=nuwk2w
      iuwk=iuwk2w+nlwkyl
      juwk=juwk2w
      itrack=4
      call putext(iaa)
c  type 2x
      nuwk=nuwk2x
      iuwk=iuwk2x+nlwkyl
      juwk=juwk2x
      itrack=10
      call putext(iaa)
c  type 3w
      nuwk=nuwk3w
      iuwk=iuwk3w
      juwk=juwk3w+nlwkyl
      itrack=14
      call putext(iaa)
c  type 3x
      nuwk=nuwk3x
      iuwk=iuwk3x
      juwk=juwk3x+nlwkyl
      itrack=11
      call putext(iaa)
      goto 700
300   continue
c  cases where i=j
c  type 8w
      nuwk=nuwk1w
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=4
      call putext(iaa)
c  type 8x
      nuwk=nuwk1x
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack=10
      call putext(iaa)
c  type 10w
      nuwk=nuwk5w
      iuwk=iuwk9w
       juwk=juwk3w+nlwkyl
      itrack=16
      call putext(iaa)
      goto 700
350   continue
c  cases where j=k (types 6,7)
c  type 6w
      nuwk=nuwk3w
      iuwk=iuwk3w
      juwk=juwk3w+nlwkyl
      itrack=4
      call putext(iaa)
c  type 6x
      nuwk=nuwk3x
      iuwk=iuwk3x
      juwk=juwk3x+nlwkyl
      itrack=9
      call putext(iaa)
c  type 7w
      nuwk=nuwk5w
      iuwk=iuwk5w+nlwkyl
      juwk=juwk7w
      itrack=16
      call putext(iaa)
      goto 700
400   continue
c  cases where i=j=k (type 12)
c  type 12(a+c)w
      nuwk=nuwk5w
      iuwk=iuwk5w+nlwkyl
      juwk=juwk5w+nlwkw(l+1)-1
      acoef=sqrt2
      itrack=13
      call putext(iaa)
c  type 12(b+c)w
      iuwk=iuwk9w
      juwk=juwk9w+nlwkyl
      itrack=7
      call putext(iaa)
c  type 12wy and 12xy
       acoef=1.0d0
       itrack=3
      ipt1=l
       do 410 iz=1,nsym
       nlwk=nlwky(ipt1)
       if(nlwk.eq.0) go to 410
      lkup=lkupl+iz
      ihsym=lkupsm(lkup)
       kpt1=ismoff(ihsym)+i
       kpt2=ismoff(ihsym)+l
       nuwk=nuwkw(ihsym)
       iuwk=ipuwkw(ihsym)+nlwkw(kpt1)
       juwk=ipuwkw(ihsym)+nlwkw(kpt2)
      call putext(iaa)
       nuwk=nuwkx(ihsym)
       iuwk=ipuwkx(ihsym)+nlwkx(kpt1)
       juwk=ipuwkx(ihsym)+nlwkx(kpt2)
      call putext(iaa)
 410   ipt1=ipt1+levfrm
c  type 12yz
       nlwk=1
       ipt1=lsoff+i
       nuwk=nuwky(ipt1+1)
      iuwk=ipuwky(ipt1+1)+nlwky(ipt1)
      juwk=ipuwky(ipt1+1)+nlwkyl
      call putext(iaa)
c  types 12wz and 12xz
       izmax=i-1
       izmin=l+1
      if(izmax.lt.izmin)goto 145
       do 440 iz=izmin,izmax
      izsym=isym(iz)
      lkup=lkupl+izsym
      ihsym=lkupsm(lkup)
       kpt1=ismoff(ihsym)+i
       kpt2=ismoff(ihsym)+iz
       ipt1=ismoff(izsym)+iz
       nuwk=nuwkw(ihsym)
       iuwk=ipuwkw(ihsym)+nlwkw(kpt1)+nlwky(ipt1)
       juwk=ipuwkw(ihsym)+nlwkw(kpt2)+nlwkyl
       acoef=1.0d0
      call putext(iaa)
       nuwk=nuwkx(ihsym)
       iuwk=ipuwkx(ihsym)+nlwkx(kpt1)+nlwky(ipt1)
       juwk=ipuwkx(ihsym)+nlwkx(kpt2)+nlwkyl
       acoef=-1.0d0
      call putext(iaa)
 440   continue
145   continue
      acoef=1.0d0
       go to 700
500   continue
c  cases where l=k (types 4,5)
c  type 4w
      nuwk=nuwk1w
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=4
      call putext(iaa)
c  type 4x
      nuwk=nuwk1x
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack=10
      call putext(iaa)
c  type 5w
      nuwk=nuwk5w
      iuwk=iuwk3w
      juwk=juwk5w+nlwkw(l+1)-1
      itrack=16
      call putext(iaa)
      goto 700
550   continue
c  cases where i=j and k=l (types 11,13)
c  type 11w
      nuwk=nuwk1w
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=4
      call putext(iaa)
c  type 11x
      nuwk=nuwk1x
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack=10
      call putext(iaa)
c  type 13w
      nuwk=nuwk5w
      iuwk=iuwk9w
      juwk=juwk5w+nlwkw(l+1)-1
      itrack=1
      call putext(iaa)
      goto 700
600   continue
c  cases where i=j=k=l (type 14)
c  type 14(a+b)w
      nuwk=nuwk5w
      iuwk=iuwk9w
      juwk=iuwk
      itrack=5
      call putext(iaa)
c  types 14wy and 14xy
       itrack=2
      ipt1=i
       do 630 iz=1,nsym
       nlwk=nlwky(ipt1)
       if(nlwk.eq.0) go to 630
      lkup=lkupl+iz
      ihsym=lkupsm(lkup)
       kpt1=ismoff(ihsym)+i
       nuwk=nuwkw(ihsym)
       iuwk=ipuwkw(ihsym)+nlwkw(kpt1)
       juwk=iuwk
      call putext(iaa)
       nuwk=nuwkx(ihsym)
       iuwk=ipuwkx(ihsym)+nlwkx(kpt1)
       juwk=iuwk
      call putext(iaa)
 630   ipt1=ipt1+levfrm
c  type 14yz
       nlwk=1
       nuwk=nuwky(ihl+1)
      iuwk=ipuwky(ihl+1)+nlwkyl
       juwk=iuwk
      call putext(iaa)
700   continue
750   continue
800   continue
      return
c threex does the three external loops
      entry threex(iaa,nabcb,nlwks,nuwks,ipuwk,nuwky,nlwkx,
     *ipuwky,nlwkw,nlwky,ijadd,kadd,ladd,nabca,nabcs,nrows,nij,
     *nkl,nlng)
      ccoef=1.0d0
      dcoef=-1.0d0
      acoeft=acoef
      iuwkt=iuwk
      juwkt=juwk
      nlwk=1
      jmn=jmin
      if(jmin.lt.2)jmn=2
      jmx=jmax
      if(jmax.ge.levfrm)jmx=levfrm-1
      ksbl=ksb-mfrm
      jfsym=mfrmpt(ksbl)
      if(jfsym)2010,940,1000
1000  continue
      ih1off=ismoff(jfsym)
      lkupj=mults(jsmt)
      lkup=lkupj+jfsym
      ifsym=lkupsm(lkup)
      do 2000 j=jmn,jmx
      kkind=1
      ij=iadt+j
      jad=ijadd(ij)
      isj=isym(j)
      lkupsj=mults(isj)
      lkup=lkupj+isj
      ksm=lkupsm(lkup)
      lkupk=mults(ksm)
      ipt1=ih1off+j
      juwk1w=juwkt+nlwkw(ipt1)
      ksmptx=msmoff(ksm)
      jsoff=ismoff(isj)
      ipt2=jsoff+j
      iuwk3w=iuwkt+nlwky(ipt2)
      do 1900 k=1,j
      ihsm3=1
      if(ifsym.eq.isj)ihsm3=0
      ksmpt=k+ksmptx
      kad=jad+kadd(ksmpt)
      isk=isym(k)
      lkup=lkupk+isk
      lsm=lkupsm(lkup)
      lsoff=ismoff(lsm)
      lsmptx=msmoff(lsm)
      ksoff=ismoff(isk)
      ipt1=ih1off+k
      juwk3w=juwkt+nlwkw(ipt1)
      ihsm1=1
      if(ifsym.eq.isk)ihsm1=0
      lkup=lkupsj+isk
      ihsym2=lkupsm(lkup)
      ihsm2=1
      if(ihsym2.eq.jfsym)ihsm2=0
      if(j.eq.k)goto 1200
      ipt2=ksoff+k
      iuwk1w=iuwkt+nlwky(ipt2)
      juwk2w=juwk1w+nlwky(ipt2)
      goto 1250
1200  continue
      kkind=kkind+1
      juwk7w=juwkt+nlwkw(k+1)-1
1250  continue
      lkind=kkind
      do 1800 l=1,k
      if(isym(l).ne.lsm)goto 1800
      lsmpt=l+lsmptx
      lad=kad+ladd(lsmpt)
      ihl=lsoff+l
      nlwkyl=nlwky(ihl)
      if(l.eq.k)lkind=lkind+2
      goto (1300,1400,1500,1800),lkind
1300  continue
c cases where no indices are the same (1-3)
c type 1w
      if(ihsm1.ne.0)goto 1310
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=13
      call putext(iaa)
c type 2w
1310  continue
      if(ihsm2.ne.0)goto 1320
      iuwk=iuwkt+nlwkyl
      juwk=juwk2w
      itrack=17
      call putext(iaa)
c type 3w
1320  continue
      if(ihsm3.ne.0)goto 1800
      iuwk=iuwk3w
      juwk=juwk3w+nlwkyl
      itrack=7
      call putext(iaa)
      goto 1800
1400  continue
c cases where j=k (6,7)
c type 6w
      if(ihsm3.ne.0)goto 1410
      iuwk=iuwk3w
      juwk=juwk3w+nlwkyl
      itrack=17
      call putext(iaa)
c type 7w
1410  continue
      if(jfsym.ne.1)goto 1800
      iuwk=iuwkt+nlwkyl
      juwk=juwk7w
      itrack=1
      acoef=acoef*sqrt2
      call putext(iaa)
      acoef=acoeft
      goto 1800
1500  continue
c cases where l=k (4,5)
c type 4w
      if(ihsm1.ne.0)goto 1510
      iuwk=iuwk1w
      juwk=juwk1w+nlwkyl
      itrack=17
      call putext(iaa)
c type 5w
1510  continue
      if(jfsym.ne.1)goto 1800
      iuwk=iuwk3w
      juwk=juwkt+nlwkw(l+1)-1
      itrack=1
      acoef=acoef*sqrt2
      call putext(iaa)
      acoef=acoeft
1800  continue
1900  continue
2000  continue
      return
2010  continue
      jfsym=-jfsym
      ih1off=ismoff(jfsym)
      lkupj=mults(jsmt)
      lkup=lkupj+jfsym
      ifsym=lkupsm(lkup)
      do 4000 j=jmn,jmx
      kkind=1
      ij=iadt+j
      jad=ijadd(ij)
      isj=isym(j)
      lkupsj=mults(isj)
      lkup=lkupj+isj
      ksm=lkupsm(lkup)
      lkupk=mults(ksm)
      ksmptx=msmoff(ksm)
      ipt1=ih1off+j
      juwk1x=juwkt+nlwkx(ipt1)
      jsoff=ismoff(isj)
      ipt2=jsoff+j
      iuwk3x=iuwkt+nlwky(ipt2)
      do 3900 k=1,j
      ihsm3=1
      if(ifsym.eq.isj)ihsm3=0
      ksmpt=k+ksmptx
      kad=jad+kadd(ksmpt)
      isk=isym(k)
      lkup=lkupk+isk
      lsm=lkupsm(lkup)
      lsoff=ismoff(lsm)
      lsmptx=msmoff(lsm)
      ksoff=ismoff(isk)
      ipt1=ih1off+k
      juwk3x=juwkt+nlwkx(ipt1)
      ihsm1=1
      if(ifsym.eq.isk)ihsm1=0
      lkup=lkupsj+isk
      ihsym2=lkupsm(lkup)
      ihsm2=1
      if(ihsym2.eq.jfsym)ihsm2=0
      if(j.eq.k)goto 3200
      ipt2=ksoff+k
      iuwk1x=iuwkt+nlwky(ipt2)
      juwk2x=juwk1x+nlwky(ipt2)
      goto 3250
3200  continue
      kkind=kkind+1
3250  continue
      lkind=kkind
      do 3800 l=1,k
      if(isym(l).ne.lsm)goto 3800
      lsmpt=l+lsmptx
      lad=kad+ladd(lsmpt)
      ihl=lsoff+l
      nlwkyl=nlwky(ihl)
      if(l.eq.k)lkind=lkind+2
      goto (3300,3400,3500,3800),lkind
3300  continue
c cases where no indices are the same (1-3)
c type 1x
      if(ihsm1.ne.0)goto 3310
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack=18
      call putext(iaa)
c type 2x
3310  continue
      if(ihsm2.ne.0)goto 3320
      iuwk=iuwkt+nlwkyl
      juwk=juwk2x
      itrack=10
      call putext(iaa)
c type 3x
3320  continue
      if(ihsm3.ne.0)goto 3800
      iuwk=iuwk3x
      juwk=juwk3x+nlwkyl
      itrack=11
      call putext(iaa)
      goto 3800
3400  continue
c cases where j=k (6,7)
c type 6x
      if(ihsm3.ne.0)goto 3800
      iuwk=iuwk3x
      juwk=juwk3x+nlwkyl
      itrack=9
      call putext(iaa)
      goto 3800
3500  continue
      if(ihsm1.ne.0)goto 3800
c cases where l=k (4,5)
c type 4x
      iuwk=iuwk1x
      juwk=juwk1x+nlwkyl
      itrack=10
      call putext(iaa)
3800  continue
3900  continue
4000  continue
      return
 930  call caserr('point not w,x,y,or z called in initex')
940   call caserr('entered a y point')
      return
      end
      subroutine loopy(iarc,nrows,nabc,nlwks,nuwks,ipuwk,iwght,nelec,
     *  ijadd,ijgrp,kadd,ladd,inext,jmnnxt,jmxnxt,nxt,ngrps,
     *nrows4,nij,nkl,iext1,iext2,iext3,iext4,iext5,iext6,iext7)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/iofile)
      common/cone/mna(9),lkupsm(64)
      common/junk/nuwkd(6,340),nelin(2),nuwko(6,340),melin(4),
     * coeffs(20,21),acoef(maxorb),bcoef(maxorb),
     * isegm(maxorb),jsegm(maxorb),imain(maxorb),isub(maxorb),
     * iuwkmn(maxorb),iuwksb(maxorb),
     * itrack(maxorb),lmin(maxorb),ishift(4)
      common/lsort /norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     *,levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),levfrm
      common/loops /lad,nuwk,nlwk,iuwk,juwk,itrak,acf,ccf,d  ,val,
     * levi,jmin,jmax,iad,jsmt,ksbt,intsad,nmax,acrcy,nblks
c
       dimension ismoff(8)
      dimension iarc(*)
      dimension nabc(nrows),nlwks(nrows),nuwks(nrows),
     *ipuwk(nrows),iwght(nrows4),ijadd(nij),ijgrp(nij),kadd(nkl),
     *ladd(nkl),inext(200),jmnnxt(200),jmxnxt(200)
      dimension nelec(nrows),nlcsmn(22)
       dimension lcond(8)
       dimension jsegnr(22),jsegpt(22),iarcmn(228),iarcsb(228)
       dimension itrk(228),jcond(228),kcond(228),nxtseg(228)
       dimension cfs(420)
       dimension jsegpx(3)
      dimension mults(8)
c
      equivalence (coeffs(1,1),cfs(1))
c
       data mults/0,8,16,24,32,40,48,56/
       data jsegnr/16,34,52,63,75,92,102,118,128,137,148,155,162,172,
     a 179,186,193,200,207,214,221,228/
       data jcond/12*1,4*-1,13*1,5*-1,13*1,5*-1,176*0/
       data kcond/7*1,9*0,5*1,13*0,5*1,13*0,7*1,4*0,7*1,5*0,11*1,6*0,
     a 5*1,5*0,10*1,6*0,5*1,5*0,5*1,4*0,5*1,6*0,1,1,5*0,1,1,5*0,4*1,
     b 62*0/
       data itrk/10,2,2,3,7,3,7,1,3*9,5*1,0,2*10,2,10,4*11,2*9,3,3,6*0,
     a  10,10,2,10,4*11,2*9,2*3,27*0,1,3*0,4*12,8*0,2*3,1,2*0,2*1,5*0,
     b  2*1,4*0,4*1,6*0,1,0,0,1,1,14*0,1,35*0,8,6*0,8,6*0,13,6*0,13,33*0
     c /
       data nxtseg/3*0,17,15,18,16,10,5,5,4,7,2,2,3,3,20,3*21,22,11,11,
     a 12,12,6,6,7,7,4*2,3,19,3*22,21,11,11,13,13,6,6,9,9,4*3,2,3*0,
     b 21,21,22,22,4*4,3*0,21,21,22,22,4*5,7,3*0,21,21,22,22,19,19,20,
     c 20,4*6,7,9,0,21,21,20,20,4*7,8,0,0,21,21,22,22,19,19,20,20,4*8,
     d 7,9,0,22,22,19,19,4*9,8,0,21,21,22,22,4*10,0,21,21,22,22,4*11,
     e 12,13,21,21,3*12,12,14,22,22,4*13,14,21,21,22,22,4*14,12,13,
     f 0,0,4*15,16,0,0,4*16,15,0,0,4*17     ,18,0,0,4*18,17,0,0,4*19,
     g 20,0,0,20,3*20,19,0,0,4*21,22,0,0,4*22,21/
      data iarcmn/4,3,2,3,4,2,4,4,2,3,4,3,3,4,2,4,1,3,4,3*2,4,3,4,1,2,1,
     a 3,1,2,3,4,2,1,2,4,3*3,4,2,4,1,3,1,2,1,2,3,4,3,2,3,4,3,4,2,4,
     b 1,2,3,4,2,3,4,3,4,2,4,1,2,3,4,3,2,3,4,3,4,2,4,1,2,1,3,1,2,3,4,
     c 3,2,2,2,4,1,2,1,2,3,4,2,2,3,3,4,2,4,1,2,1,3,1,2,3,4,3,2,3,3,4,
     d 1,3,1,2,3,4,3,1,1,3,1,2,1,2,3,4,1,1,3,1,2,1,2,3,4,3,2,1,2,1,2,
     e 3,4,2,1,3,1,2,3,4,3,1,3,1,2,1,2,3,4,3,2,1,2,1,2,3,4,2,1,3,1,2,
     f 3,4,3,1,2,1,2,3,4,2,1,3,1,2,3,4,3,3,4,1,2,3,4,3,2,4,1,2,3,4,2,
     g 1,2,1,2,3,4,2,1,3,1,2,3,4,3/
      data iarcsb/4,3,2,1,2,1,3,1,2,3,4,2,1,2,1,3,4,3,4,2,3,1,3,1,2,3,4,
     a 2,4,1,2,3,4,3,4,2,4,3,2,1,2,1,3,2,4,3,4,1,2,3,4,2,2,3,4,1,2,1,3,
     b 1,2,3,4,2,3,4,1,2,1,3,1,2,3,4,2,2,3,4,1,2,1,3,3,4,2,4,1,2,3,4,
     c 2,3,3,1,3,3,4,1,2,3,4,3,2,3,1,2,1,3,3,4,2,4,1,2,3,4,2,3,2,1,2,
     d 2,4,1,2,3,4,2,4,2,4,3,4,1,2,3,4,4,2,4,3,4,1,2,3,4,2,3,3,4,1,2,
     e 3,4,3,2,4,1,2,3,4,2,2,4,3,4,1,2,3,4,2,3,3,4,1,2,3,4,3,2,4,1,2,
     f 3,4,2,3,4,1,2,3,4,3,2,4,1,2,3,4,2,1,2,1,2,3,4,2,1,3,1,2,3,4,3,
     g 3,4,1,2,3,4,3,2,4,1,2,3,4,2/
       data jsegpx/12,29,47/
       data nlcsmn/6*1,0,1,14*0/
c
      lvfrm=levfrm
      icnt88=0
      crite = 0.00001d0
      root2 =  dsqrt(2.0d0)
      rootn2 = -root2
      toor2 = 1.0d0/ root2
      toorn2 = -toor2
      jsegpt(1)=0
      do 130 i=1,21
 130  jsegpt(i+1)=jsegnr(i)
      do 135 i=1,2
      do 135 j=1,21
 135  coeffs(i,j)=0.0d0
      do 137 i=1,nij
137   ijadd(i)=ijadd(i)+intsad+2
      do 140 i=3,20
      a =  dfloat(i-2)
      coeffs(i,1) =  dsqrt(a/(a+1.0d0))
      coeffs(i,2) = -coeffs(i,1)
      coeffs(i,3) = coeffs(i,1)/ dsqrt(2.0d0)
      coeffs(i,4) = -coeffs(i,3)
      coeffs(i,5) =  dsqrt((a+1.0d0)/a)
      coeffs(i,6) = -coeffs(i,5)
      coeffs(i,7) = coeffs(i,5)/ dsqrt(2.0d0)
      coeffs(i,8) = -coeffs(i,7)
      coeffs(i,9) =  dsqrt((a+2.0d0)/(a*2.0d0))
      coeffs(i,10) = -coeffs(i,9)
      coeffs(i,11) =  dsqrt(a/(2.0d0*(a+2.0d0)))
      coeffs(i,12) = -coeffs(i,11)
      coeffs(i,13) =  dsqrt(2.0d0/(a*(a+1.0d0)))
      coeffs(i,14) =  dsqrt(a*(a+2.0d0))/(a+1.0d0)
      coeffs(i,15) = - dsqrt(a*(a+2.0d0))/(a+1.0d0)
      coeffs(i,16) =  dsqrt((a-1.0d0)*(a+2.0d0)/(a*(a+1.0d0)))
      coeffs(i,17) = -coeffs(i,16)
      coeffs(i,18)=- dsqrt(2.0d0/(a*(a+2.0d0)))
      coeffs(i,19) = 1.0d0/a
      coeffs(i,20) = -1.0d0/a
      coeffs(i,21) = - dsqrt(2.0d0)/a
140   continue
      do 150 k=1,4
      ishift(k)=(k-1)*nrows
150   continue
      do 155 i=1,nrows
      nelec(i)=nelec(i)+nelec(i)+nabc(i)
155   continue
      do 160 i=1,nsym
      ismoff(i)=(i-1)*norbs
 160  lcond(i)=0
      i=isym(1)
      lcond(i)=1
      lcond(1)=1
      nsm=0
      do 170 iorb=2,norbs
      do 165 i=1,nsym
      if(lcond(i).eq.0) go to 165
      ism=i
      lkup=mults(ism)+isym(iorb)
      j=lkupsm(lkup)
      if(lcond(j).gt.0) go to 165
      lcond(j)=iorb
      nsm=nsm+1
      if(nsm.eq.nsym) go to 175
 165  continue
 170  continue
 175  continue
      do 180 i=1,nsym
      if(lcond(i).eq.0) lcond(i)=norbs+1
 180  continue
      inxt=1
      write(iwr,186)
 186  format(/)
      do 500 iblock=1,ngrps
      write(iwr,187)iblock
187   format(' process integrals from group',i4)
190   continue
      i=inext(inxt)
      levi=i+1
      jmax=jmxnxt(inxt)
      jmin=jmnnxt(inxt)
      iad=(i*(i-1))/2
      ij=iad+jmax
      if(iblock.ne.ijgrp(ij)) go to 490
      if(levi.gt.levfrm) go to 195
      call extern(iarc
     *           ,nabc,nlwks,nuwks,ipuwk,iext1,iext2,iext3,iext4,iext5,
     *ijadd,kadd,ladd,nelec,iext6,nrows,nij,nkl,iext7)
      go to 485
195   continue
      jsm=isym(i)
      lev=levi
      levm=lev-1
      nr=levnr(lev)
      npt=levpt(lev)
      do 480 irow=1,nr
      npt=npt+1
      isegm(lev)=1
      iseg=1
      imn=npt
      isb=npt
      kseg=0
      ksegmx=jsegnr(iseg)
      lmin(lev)=lcond(jsm)
      iuwkmn(lev)=ipuwk(npt)
      iuwksb(lev)=ipuwk(npt)
      imain(lev)=npt
      isub(lev)=npt
      nuwk=nuwks(npt)
      acoef(lev)=1.0d0
c  test next segment of group
200   kseg=kseg+1
      if(kseg.gt.ksegmx) go to 440
      kmn=iarcmn(kseg)
      iarpt=imn+ishift(kmn)
      kmn=iarc(iarpt)
      if(kmn.eq.0) go to 200
      ksb=iarcsb(kseg)
      jarpt=isb+ishift(ksb)
      ksb=iarc(jarpt)
      if(ksb.eq.0) go to 200
      jsegm(lev)=kseg
      iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
      iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
      lmin(levm)=lmin(lev)
      if(jcond(kseg))220,240,230
 220  continue
      if(levm.le.jmin)goto 440
      goto 240
 230  continue
      if(levm.gt.jmax)goto 420
      ij=iad+levm
      jad=ijadd(ij)
      lkup=mults(jsm)+isym(levm)
      ksm=lkupsm(lkup)
      lmin(levm)=lcond(ksm)
      ksmptx=ismoff(ksm)
 240  continue
      if(kcond(kseg).eq.0)goto 260
      ksmpt=levm+ksmptx
      kad=jad+kadd(ksmpt)
      lkup=mults(ksm)+isym(levm)
      lsm=lkupsm(lkup)
      lmin(levm)=lcond(lsm)
      lsmptx=ismoff(lsm)
 260  continue
      if(itrk(kseg))280,280,270
270   itrack(levm)=itrk(kseg)
      goto 290
280   itrack(levm)=itrack(lev)
290   continue
_IF(sgi,sun,ksr,ipsc,rs6000,hp700,hpux11)
      if(kseg.le.120) then
      go to(50,1,1,1,3,1,4,1,44,45,51,1,1,3,1,4,5,52,53,6,54,40,41,1,7,
_ELSE
       go to(50,1,1,1,3,1,4,1,44,45,51,1,1,3,1,4,5,52,53,6,54,40,41,1,7,
_ENDIF
     a 46,47,1,8,1,6,2,2,9,10,52,53,11,55,42,43,1,12,48,49,1,13,1,2,11,
     b 2,36,77,77,78,77,79,77,80,1,1,1,1,87,68,82,68,69,67,70,71,75,
     c 76,71,83,87,68,82,68,69,67,70,68,69,67,70,71,75,76,71,83,83,
     d 6,6,74,6,8,1,16,16,1,17,18,19,19,20,18,21,19,20,18,21,1,22,23,
_IF(sgi,sun,ksr,ipsc,rs6000,hp700,hpux11)
     e 1,24,24,11,11) kseg
      else
       kseg1=kseg-120
      go to( 81,11,13,1,27,27,1,28,1,3,2,4,2,1,2,2,1,29,63,64,
_ELSE
     e 1,24,24,11,11,81,11,13,1,27,27,1,28,1,3,2,4,2,1,2,2,1,29,63,64,
_ENDIF
     f 65,66,71,72,73,71,84,85,56,57,1,30,1,1,86,56,58,1,1,31,1,32,
     g 59,60,61,62,1,33,88,1,34,35,1,5,1,6,2,2,9,1,10,1,2,11,2,36,
     h 1,5,1,6,2,2,9,1,10,1,2,11,2,36,1,4,1,37,2,2,38,1,3,1,2,37,2,39,
_IF(sgi,sun,ksr,ipsc,rs6000,hp700,hpux11)
     i 1,5,1,6,2,2,9,1,10,1,2,11,2,36),kseg1
      endif
_ELSE
     i 1,5,1,6,2,2,9,1,10,1,2,11,2,36),kseg
_ENDIF
 1    acoef(levm)=acoef(lev)
      goto 120
  2   acoef(levm)=-acoef(lev)
      goto 120
  3   ia = nabc(imn) + 2
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  4   ia = nabc(imn) + 83
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  5   ia = nabc(imn) + 82
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  6   ia = nabc(imn) + 261
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  7   ia = nabc(imn) + 1
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  8   ia = nabc(imn) + 102
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  9   ia = nabc(imn) + 362
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 10   ia = nabc(imn) + 3
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 11   ia = nabc(imn) + 263
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 12   ia = nabc(imn) + 84
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 13   ia = nabc(imn) + 23
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 16   ia = nabc(imn) + 281
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 17   ia = nabc(imn) + 402
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 18   ia = nabc(imn) + 162
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 19   ia = nabc(imn) + 222
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 20   ia = nabc(imn) + 143
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 21   ia = nabc(imn) + 42
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 22   ia = nabc(imn) + 302
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 23   ia = nabc(imn) + 303
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 24   ia = nabc(imn) + 342
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 27   ia = nabc(imn) + 283
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 28   ia=nabc(imn)+404
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 29   acoef(levm) = acoef(lev) * root2
      go to 120
 30   ia = nabc(imn) + 301
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 31   ia = nabc(imn) + 304
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 32   ia = nabc(imn) + 244
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 33   ia = nabc(imn) + 322
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 34   ia = nabc(imn) + 243
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 35   ia = nabc(imn) + 242
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 36   ia = nabc(imn) + 384
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 37   ia = nabc(imn) + 262
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 38   ia = nabc(imn) + 363
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 39   ia = nabc(imn) + 383
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 86   ia = nabc(imn) + 241
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 40   ia = nabc(imn) + 122
      ib=ia-61
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 41   ib = nabc(imn) + 162
      acoef(levm) = acoef(lev) * toorn2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 42   ia = nabc(imn) + 43
      ib = ia + 81
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 43   ib = nabc(imn) + 222
      acoef(levm) = acoef(lev) * toorn2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 44   ib=nabc(imn)+221
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 45   ib = nabc(imn) + 163
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 46   ib = nabc(imn) + 162
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 47   ia = nabc(imn) + 122
      ib = ia - 81
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 48   ib = nabc(imn) + 222
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 49   ia = nabc(imn) + 43
      ib = ia + 101
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 50   acoef(levm) = acoef(lev) + acoef(lev)
      d=0.5d0
      go to 120
 51   acoef(levm)=acoef(lev)*root2
      go to 120
 52   acoef(levm) = -acoef(lev)
      d= -1.0d0
      go to 120
 53   acoef(levm) = -acoef(lev) - acoef(lev)
      d = -0.5d0
      go to 120
 54   ia=nabc(imn)+362
      d=1.0d0/cfs(ia)
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 55   ia = nabc(imn) + 384
      d=1.0d0/cfs(ia)
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 56   acoef(levm) = acoef(lev)
      d = -1.0d0
      go to 120
 57   ia = nabc(imn) + 82
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 58   ia = nabc(imn) + 3
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 59   ia = nabc(imn) + 123
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 60   ia = nabc(imn) + 222
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 61   ia = nabc(imn) + 62
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 62   ia = nabc(imn) + 162
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 63   ia = nabc(imn) + 42
      ib = ia + 81
      acof = acoef(lev) * cfs(ia)
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 64   ib = nabc(imn) + 222
      acof = acoef(lev) * toorn2
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 65   ia = nabc(imn) + 123
      ib = ia - 61
      acof = acoef(lev) * cfs(ia)
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 66   ib = nabc(imn) + 162
      acof = acoef(lev) * toorn2
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 67   ib = nabc(imn) + 162
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 68   ib = nabc(imn) + 222
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 69   ia = nabc(imn) + 62
      ib = ia + 81
      dx=acoef(lev)*cfs(ia)
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 70   ia = nabc(imn) + 143
      ib = ia - 101
      dx=acoef(lev)*cfs(ia)
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 87   ib = nabc(imn) + 162
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm)=d
      d=-(dx+dx)/d
      go to 120
 71   acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev)
      go to 120
 72   ib = nabc(imn) + 322
      acoef(levm) = -acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 73   ib = nabc(imn) + 323
      acoef(levm) = -acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 74   ia=nabc(imn)+21
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 75   ib = nabc(imn) + 302
      acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 76   ib = nabc(imn) + 303
      acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 77   acoef(levm)=acoef(lev)*toorn2
      d=-2.0d0
      go to 120
 78   acoef(levm)=acoef(lev)*rootn2
      d=-2.0d0
      go to 120
 79   ia=nabc(imn)+62
      acoef(levm)=acoef(lev)*cfs(ia)
      d=-2.0d0
      go to 120
 80   ia=nabc(imn)+143
      acoef(levm)=acoef(lev)*cfs(ia)
      d=-2.0d0
      go to 120
 81   ia=nabc(imn)+104
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 82   acoef(levm) = acoef(lev) * rootn2
      d=-2.0d0
      go to 120
 83   ia = nabc(imn) + 342
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 84   ia = nabc(imn) + 243
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 85   ia = nabc(imn) + 242
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 88   ia = nabc(imn) + 323
      acoef(levm) = acoef(lev) * cfs(ia)
      icnt88=icnt88+1
      go to 120
110   itrack(levm)=3
      acoef(levm)=acof-bcof
      go to 120
111   itrack(levm) = 2
      acoef(levm)=-(dx+dx)
120   continue
      if(nxtseg(kseg).gt.0) go to 400
      if(isym(levm).ne.lsm) go to 200
      lsmpt=levm+lsmptx
      lad=kad+ladd(lsmpt)
      if(kmn-ksb) 300,380,300
300   levl=levm
      ksegmx=4
310   lev=levm
      levm=lev-1
      if(levm.gt.0) go to 315
      call caserr('problems with partial space')
 315  continue
      kseg=0
      imain(lev)=kmn
      imn=kmn
      isub(lev)=ksb
      isb=ksb
320   kseg=kseg+1
      if(kseg.gt.ksegmx)goto 360
      iarpt=imn+ishift(kseg)
      kmn=iarc(iarpt)
      if(kmn.le.0)goto 320
      jarpt=isb+ishift(kseg)
      ksb=iarc(jarpt)
      if(ksb.le.0)goto 320
      jsegm(lev)=kseg
      iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
      iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
      if(kmn-ksb) 310,340,310
340   nlwk = nlwks(kmn)
      iuwk = iuwkmn(levm)
      juwk = iuwksb(levm)
      itrak = itrack(levl)
      acf = acoef(levl)
      call putout(iarc)
      go to 320
360   if(lev.eq.levl) go to 440
      levm=lev
      lev=levm+1
      imn=imain(lev)
      isb=isub(lev)
      kseg=jsegm(lev)
      go to 320
380   nlwk = nlwks(kmn)
      iuwk = iuwkmn(levm)
      juwk = iuwksb(levm)
      itrak = itrack(levm)
      acf = acoef(levm)
      call putout(iarc)
      go to 200
400   continue
      if(levm.le.lmin(levm))goto 200
      iseg=nxtseg(kseg)
      if(nlcsmn(iseg).gt.0.and.nelec(kmn).eq.0)goto 200
      if(lvfrm.eq.levm.and.iseg.lt.4)goto 460
      lev=levm
      levm=lev-1
      isegm(lev)=iseg
      kseg=jsegpt(iseg)
      imn=kmn
      imain(lev)=kmn
      isb=ksb
      isub(lev)=ksb
      ksegmx=jsegnr(iseg)
      goto 200
420   continue
      kseg=jsegpx(iseg)
      goto 200
440   continue
      if(lev.eq.levi)goto 480
      levm=lev
      lev=levm+1
      iseg=isegm(lev)
      imn=imain(lev)
      isb=isub(lev)
      kseg=jsegm(lev)
      ksegmx=jsegnr(iseg)
      goto 200
460   continue
      iuwk=iuwkmn(levm)
      juwk=iuwksb(levm)
      acf=acoef(levm)
      jsmt=jsm
      ksbt=ksb
      call threex(iarc
     *           ,nabc,nlwks,nuwks,ipuwk,iext1,iext2,iext3,iext4,iext5,
     *ijadd,kadd,ladd,nelec,iext6,nrows,nij,nkl,iext7)
      goto 200
480   continue
485   continue
      inxt=inxt+1
      if(inxt.le.nxt) go to 190
      go to 500
490   call nxtblk(iarc)
500   continue
      call finout
      if(icnt88.gt.0)write(iwr,510)icnt88
 510  format(' *******'/1h0,'jseg 169 executed',i8,' times'/' *******')
      return
      end
      subroutine casout(nbf,nwks,indx,itape8)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/files)
      integer    mxsurd
      parameter (mxsurd=3600)
      common/three /vsurd(mxsurd),nnsurd
INCLUDE(common/machin)
      common/craypk /surd(mxsurd),nsurd,nbasis,ncore,nact,
     1 ntype,itype(maxorb),nwk,naa,nbb
      common/lsort /norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb),
     1        levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),levfrm
c
      dimension indx(*)
c
      nav = lenwrd()
      nwk=nwks
      nact=norbs
      ntype=nsym
      nbasis=nbf
      nsurd=nnsurd
      read(itape8)itype
      read(itape8)ncore,nactt,naa,nbb
      if(nact.ne.nactt)call caserr('problem with nact')
      call dcopy(nsurd,vsurd,1,surd,1)
      call wrt3(surd,mach(17),n9blk(1),n9tape(1))
      nwksh=(nwks-1)/nav+1
      call wrt3is(indx,nwksh*nav,n9tape(1))
      return
      end
      subroutine drt(ia,max0)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
      common /lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     * ,icode(maxorb),nlcs(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),nbff
c
c.....max. number of walks allowed
c
      integer  mxwks
      character*10 charwall
      parameter (mxwks=80000)
c
c  option 1 - printing option
c  option 2 - excitation option
c  option 3 - input this number of excitation types
c  option 4 - include non-interacting space if not 0
c  option 5 - if not zero, set levfrm to this value
c
      dimension ia(*)
      write(iwr,60)
      itape8 = ntapes(1)
      rewind itape8
      call reord (levfrm,npflg)
      nwks=0
      nav = lenwrd()
      iwordo=2048
      max=max0*nav
      nrows=(max/9)-2
      nrowso=nrows
      ih1a=1
      ih1b=ih1a+nrows
      ih1s=ih1b+nrows
      ih1t=ih1s+nrows
      ih2=ih1t+nrows
      ih3=ih2+nrows*4
      write(iwr,80)
      nrows4=4*nrows
      call paldus (nrows,nwks,ia(ih1a),ia(ih1b),ia(ih1s),ia(ih1t),ia(ih2
     1),ia(ih3),nrows4)
      nbf=norbs
      norbs=na+nb+nc
      if (iop(5).ne.0) levfrm=iop(5)
      write (itape8) nbf,norbs,nsym,nrows,nwks,levfrm,iwordo,iop
      write (itape8) iout
      write (itape8) isym
      write (itape8) levnr
      write (itape8) levpt
      norbs=na+nb+nc
      ih4b=ih1a+nrows
      ih4s=ih4b+nrows
      ih4t=ih4s+nrows
      ih4=ih4t+nrows
      ih5=ih4+nrows*4
      ih6=ih5+nrows
      ih7=ih6+nrows
      ih8=ih7+nrows
      nwksmx=max-ih8
cmax
      if (nwksmx.gt.mxwks) nwksmx=mxwks
c
      if (nwksmx.lt.nwks) go to 10
      if (nwksmx.lt.nrows*4) go to 10
      write(iwr,90)
      nrowsf=4*nrows
      call gindex(nrows,nwks,ia(ih1a),ia(ih1b),ia(ih1s),ia(ih1t),ia(ih2)
     1,ia(ih3),ia(ih4b),ia(ih4s),ia(ih4t),ia(ih4),ia(ih5),ia(ih6),ia(ih7
     2),ia(ih8),nrowso,itape8,ia(ih8),nrowsf,nrows4,npflg)
      nkl=nsym*norbs
      nij=(norbs*(norbs+1))/2
      write(iwr,100)
      ih9=1
      ih10=ih9+nij
      ih11=ih10+nij
      ih12=ih11+nkl
      ih13=ih12+nkl
      ih14=ih13+200
      ih15=ih14+200
      nlng=levfrm*nsym
      ih33=ih15+200
      ih34=ih33+nlng
      ih35=ih34+nlng
      ih36=ih35+nlng
      ih37=ih36+nlng
      ih16=ih37+nlng
      len=max-ih16-ih8-nrowsf
      nmax=len
      if (nmax.gt.25000) nmax=25000
      ih17=ih16+nij
      call symtri (ia(ih9),ia(ih10),ia(ih11),ia(ih12),ia(ih13),ia(ih14),
     1ia(ih15),ia(ih16),ia(ih17),nmax,nxt,ngrps,itape8,nij,nkl)
      top=cpulft(1)
      write(iwr,20)top,charwall()
 20   format(/' end of distinct row table generation at ',f8.2,
     *' seconds',a10,' wall' )
      return
   10 write(iwr,130) nwks,nwksmx
      call caserr('dimensioning problem in drt')
   80 format (/24x,'construct basic distinct row table')
   90 format (/24x,'construct indexing array and other drt arrays')
  100 format (/ 24x,'construct integral addressing arrays')
  130 format (25h   too many walks.  nwks=,i8,10h   nwksmx=,i8)
   60 format(/28x,'orbital information input'/28x,25('*')/)
      return
      end
      subroutine gindex(nrows,nwks,nabca,nabcbo,nabcso,nabcto,iarco,nlwk
     1so,nabcb,nabcs,nabct,iarc,nlwks,nuwks,ipuwk,indx,nrowso,itape8,nwg
     2ht,nrowsf,nrows4,npflg)
c
c   at end of subroutine, nelec is returned in nabca
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      common /lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     * ,icode(maxorb),nlcs(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8)
INCLUDE(common/iofile)
INCLUDE(common/prnprn)
c
      dimension iarco(nrows4), nlwkso(nrowso), iarc(nrowsf), nlwks(nrows
     1), nuwks(nrows), ipuwk(nrows), indx(nwks)
      dimension nabca(nrows), nabcb(nrows), nabcs(nrows), nabct(nrows),
     1nabcbo(nrowso), nabcso(nrowso), nabcto(nrowso)
c
      common/junk/iwght(maxorb),itest(maxorb),lwght(maxorb),levir(maxorb
     *), icase(128), ieconf(128),  iwksw(5), ishfto(4),
     *  ishift(4)  , jarcs(4)
      common/junkc/xtext(5,15)
      dimension xtag(3),yref(20)
      dimension nwght(nrowsf)
c
      data yref/'doc','uoc','fzc','fzv','alp','aos','bos','gva','gvb',
     1'cor','vir','val','vbe','vuo','vdo',5*'zzz'/
       data xtag/'0','1','2'/
       data xblank/' '/
       data dzero /0.0d0/
c
c  initialization and reduction of core space
c
      do 10 i=1,nrows
   10 nabcb(i)=nabcbo(i)
      do 20 i=1,nrows
   20 nabcs(i)=nabcso(i)
      do 30 i=1,nrows
   30 nabct(i)=nabcto(i)
      do 40 k=1,4
      ishift(k)=(k-1)*nrows
      ishfto(k)=(k-1)*nrowso
   40 continue
      do 50 k=1,4
      do 50 i=1,nrows
      iarpt=i+ishift(k)
      iarpto=i+ishfto(k)
      iarc(iarpt)=iarco(iarpto)
   50 continue
      do 60 i=1,nrows
   60 nlwks(i)=nlwkso(i)
      do 70 i=1,nwks
   70 indx(i)=0
      do 80 i=1,nrows
      nuwks(i)=0
      ipuwk(i)=0
   80 continue
c
c  generate number of upper walks array (nuwks)
c  generate primary upper walk array  (puwk)
c
      levm=nlevs
      nuwks(1)=1
      ipuwk(1)=1
   90 lev=levm
      levm=lev-1
      if (levm.eq.0) go to 150
      nr=levnr(lev)
      nrm=levnr(levm)
      nptmx=levpt(levm)
      do 100 i=1,nr
  100 iwght(i)=0
      do 140 k=1,4
      nptm=levpt(levm)
      do 110 i=1,nrm
      nptm=nptm+1
      itest(i)=0
      if (ipuwk(nptm).gt.0) itest(i)=-1
  110 continue
      npt=levpt(lev)
      do 130 i=1,nr
      npt=npt+1
      iarpt=npt+ishift(k)
      jarc=iarc(iarpt)
      if (jarc.eq.0) go to 130
      ipt=jarc+nptmx
      nuwks(ipt)=nuwks(ipt)+nuwks(npt)
      if (itest(jarc).lt.0) go to 120
      itest(jarc)=1
      ipuwk(ipt)=ipuwk(npt)+iwght(i)
  120 iwght(i)=iwght(i)+nlwks(ipt)
  130 continue
  140 continue
      go to 90
  150 continue
c
c  write out arrays
c
c     nrows8=nrows*8
c
      write (itape8) nabca
      write (itape8) nabcb
      write (itape8) nabcs
      write (itape8) nlwks
      write (itape8) nuwks
      write (itape8) ipuwk
      write(iwr,350)
      write(iwr,360) nrows
      write(iwr,370) nwks
      if ((iop(1)/2)*2.eq.iop(1)) go to 200
      iop(1)=iop(1)-1
c
      npuflg = 0
      if (opunch(16)) then
       open(ipu,file='civecs.ascii',form='formatted',status='unknown')
       rewind ipu
       npuflg = 1
      endif
c
      if(npflg.ne.0)write(iwr,380)
      lev=nlevs
  160 continue
      nr=levnr(lev)
      npt=levpt(lev)
      levm=lev-1
      if(npflg.ne.0)write(iwr,390)
      if(npflg.ne.0)write(iwr,400)
      do 180 i=1,nr
      npt=npt+1
      do 170 k=1,4
      iarpt=npt+ishift(k)
      jarcs(k)=iarc(iarpt)
  170 continue
      ia=nabca(npt)
      ib=nabcb(npt)
      ic=levm-ia-ib
      itp=nabct(npt)
      ism=nabcs(npt)
      if(npflg.ne.0)write(iwr,410)
     *           levm,i,ia,ib,ic,ism,itp,jarcs,nlwks(npt),nuwks(npt),
     1ipuwk(npt)
  180 continue
      if(npflg.ne.0)write(iwr,390)
      lev=lev-1
      if (lev.eq.0) go to 200
      icd=icode(lev)
      do 190 i=1,maxorb
      if (iout(i).eq.lev) iorb=i
  190 continue
      ism=isym(lev)
c
c  generate the indexing array  (indx)
c
      if(npflg.ne.0)write(iwr,420) lev,iorb,ism,yref(icd)
      go to 160
  200 continue
      lev=1
      iwk=0
      levm=1
      ir=1
      lwght(1)=0
      iw=1
  210 continue
      lwght(lev)=lwght(levm)+iw
      if (lev.eq.nlevs) go to 250
      levir(lev)=ir
      levm=lev
      lev=levm+1
      levir(lev)=levnr(lev)+1
  220 continue
      ir=levir(lev)
      nptx=levpt(lev)
      nptmx=levpt(levm)
      irm=levir(levm)
  230 continue
      ir=ir-1
      if (ir.eq.0) go to 260
      npt=ir+nptx
      iw=0
      do 240 k=1,4
      iarpt=npt+ishift(k)
      jarc=iarc(iarpt)
      if (jarc.eq.0) go to 240
      if (irm.eq.jarc) go to 210
      ipt=jarc+nptmx
      iw=iw+nlwks(ipt)
  240 continue
      go to 230
  250 continue
      iwk=iwk+1
      iw=lwght(nlevs)
      indx(iw)=iwk
  260 continue
      lev=levm
      levm=lev-1
      if (levm.gt.0) go to 220
      if ((iop(1)/4)*4.eq.iop(1)) go to 280
      iop(1)=iop(1)-2
      if(npflg.ne.0)write(iwr,430)
      ipti=1
  270 continue
      iptf=ipti+4
      if (iptf.gt.nwks) iptf=nwks
      if (iptf.lt.ipti) go to 280
      if(npflg.ne.0)write(iwr,440) (ipt,indx(ipt),ipt=ipti,iptf)
      ipti=iptf+1
      if (iptf.lt.nwks) go to 270
  280 continue
      write (itape8) indx
      do 290 k=1,4
      do 290 i=1,nrows
      iarpt=i+ishift(k)
      nwght(iarpt)=0
  290 continue
      do 330 lev=2,nlevs
      levm=lev-1
      nr=levnr(lev)
      npt=levpt(lev)
      nptm=levpt(levm)
      do 320 i=1,nr
      npt=npt+1
      do 310 k=1,4
      iarpt=npt+ishift(k)
      j=iarc(iarpt)
      iarptx=iarpt+nrows
      if (j.gt.0) go to 300
      if (k.eq.4) go to 310
      nwght(iarptx)=nwght(iarpt)
      go to 310
  300 continue
      iptm=j+nptm
      iarc(iarpt)=iptm
      if (k.eq.4) go to 310
      nwght(iarptx)=nwght(iarpt)+nlwks(iptm)
  310 continue
  320 continue
  330 continue
      write (itape8) iarc
      write (itape8) nwght
c
c     ----- print electronic configurations -----
c
      if (npflg .eq. 0 .and. npuflg.eq.0) go to 2004
      do i=1,5
       do j=1,15
        xtext(i,j)=xblank
       enddo
      enddo
      kw=0
      if (npflg .ne. 0) write (iwr,9148)
c     if (npuflg .ne. 0) write (ipu,9148)
      ncorbs = 0
      do i = 1,128
       if (iout(i) .lt. 0) ncorbs = ncorbs+1
      enddo
      iwks = 0
      lev = 1
      levm = 1
      ir0 = 1
  780 continue
      if (lev .eq. nlevs) go to 860
      levir(lev) = ir0
      levm = lev
      lev = levm+1
      levir(lev) = levnr(lev)+1
  800 continue
      ir0 = levir(lev)
      nptx = levpt(lev)
      nptm = levpt(levm)
      irm = levir(levm)
  820 continue
      ir0 = ir0-1
      if (ir0 .eq. 0) go to 1000
      npt = ir0+nptx
      do 840 k = 1,4
      iarpt = npt+ishift(k)
      jarc = iarc(iarpt)
      if (jarc .eq. 0) go to 840
      jarc = jarc-nptm
      icase(levm) = k
      if (irm .eq. jarc) go to 780
  840 continue
      go to 820
  860 continue
c
      do ilev = 1,norbs
       icas = icase(ilev)
       go to (880,900,900,920), icas
  880  iocc = 0
       go to 940
  900  iocc = 1
       go to 940
  920  iocc = 2
  940  continue
       do i = 1,128
        if (iout(i) .eq. ilev) iorb = i
       enddo
       ieconf(iorb) = iocc
      enddo
c
      iwks = iwks+1
      kw=kw+1
      do ilev=1,norbs
       iocc=ieconf(ncorbs+ilev)
       xtext(kw,ilev)=xtag(iocc+1)
      enddo
      iwksw(kw)=iwks
      if(kw.ne.5)go to 1000
      if (npflg .ne. 0) then
       write(iwr,5)(iwksw(kw),(xtext(kw,ilev),ilev=1,15),kw=1,5)
   5  format(/5(1x,i6,3x,15a1))
      endif
      if (npuflg .ne. 0) then
       do kw =1,5
       write(ipu,7)iwksw(kw),(xtext(kw,ilev),ilev=1,15),dzero
       enddo
   7  format(1x,i6,3x,15a1,5x,f20.10)
      endif
      kw=0
 1000 continue
      lev = levm
      levm = lev -1
      if (levm .gt. 0) go to 800
 9148 format(//20x,'list of electronic configurations',
     +     ' - occupation numbers for active orbitals'/ 20x,74('*'))
      if(kw.eq.0)go to 2004
      if (npflg .ne. 0) then
       write(iwr,5)(iwksw(k),(xtext(k,ilev),ilev=1,15),k=1,kw)
      endif
      if (npuflg .ne. 0) then
       do k =1, kw
        write(ipu,7)iwksw(k),(xtext(k,ilev),ilev=1,15),dzero
       enddo
       close(ipu,status='keep')
      endif
2004  return
c
  350 format (/)
  360 format (20h  number of rows  is,i8)
  370 format (20h  number of walks is,i8)
  380 format (/1x,104('*')//20x,'the distinct row table'/
     *20x,22('=')/)
  390 format (1x)
  400 format (7x,' i   j      a  b  c   sym t      k0 k1 k2 k3',
     +        '      nlwks    nuwks     puwk')
  410 format (5x,2i4,4x,3i3,3x,i2,i3,5x,4i3,2x,3i9)
  420 format (6h level,i3,9h  orbital,i3,6h  sym ,i2,4x,a4)
  430 format (1h1,17h   indexing array,/)
  440 format (1x,5(1h(,i5,1h,,i5,1h)))
      end
      subroutine paldus (nrows,nwks,nabca,nabcb,nabcs,nabct,iarc,nlwks,
     +                   nrows4)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
c
      common/cone/kmna(9),lkupsm(64)
      common /lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     * ,icode(maxorb),nlcs(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8)
      common /valen/ ivv,ifermi
c
      dimension nabca(nrows), nabcb(nrows), nabcs(nrows), nabct(nrows)
      dimension mult8(8)
      dimension iarc(nrows4), nlwks(nrows)
      dimension ishift(4), ired(200)
c
      data mult8/0,8,16,24,32,40,48,56/
c
c  generate the paldus array  (nabc)
c
      nrefs=iop(8)
      iexct=2
      if (iop(2).ne.0) iexct=iop(2)
      iop(2)=iexct
      if (ivv.eq.1) iexct=10
      iact=idocc/256
      idocc=idocc-iact*256
      nbmax=nb+iexct+2
      if (nbmax.gt.16.or.nb.gt.14)go to 600
      do 10 k=1,4
      ishift(k)=(k-1)*nrows
   10 continue
      do 20 i=1,nrows
      nlwks(i)=0
      nabca(i)=0
      nabcb(i)=0
      nabcs(i)=0
      nabct(i)=0
      do 20 k=1,4
      iarpt=i+ishift(k)
      iarc(iarpt)=0
   20 continue
      levm=nlevs
      levpt(levm)=0
      levnr(levm)=1
      nabca(1)=na
      nabcb(1)=nb
      nabcs(1)=msym
   30 lev=levm
      levm=lev-1
      if (levm.eq.0) go to 270
      nr=levnr(lev)
      npt=levpt(lev)
      levpt(levm)=npt+nr
      nrm=0
      icd=icode(levm)
      ibm=nbmax+1
      iam=na+1
      do 260 ias=1,iam
      ia=na-ias+1
      do 260 ibs=1,ibm
      ib=nbmax-ibs+1
      npt=levpt(lev)
      do 250 i=1,nr
      npt=npt+1
      do 250 k=1,4
      iabca=nabca(npt)
      iabcb=nabcb(npt)
      iabcs=nabcs(npt)
      iabct=nabct(npt)
      go to (40,50,60,70),k
   40 iabcc=levm-iabca-iabcb
      if (iabcc.lt.0) go to 250
      if (ia.ne.iabca) go to 250
      if (ib.ne.iabcb) go to 250
      go to 90
   50 iabcb=iabcb-1
      if (iabcb.lt.0) go to 250
      if (ia.ne.iabca) go to 250
      if (ib.ne.iabcb) go to 250
      go to 80
   60 iabca=iabca-1
      iabcb=iabcb+1
      iabcc=levm-iabca-iabcb
      if ((iabca.lt.0).or.(iabcc.lt.0).or.(iabcb.gt.nbmax)) go to 250
      if (ia.ne.iabca) go to 250
      if (ib.ne.iabcb) go to 250
      go to 80
   70 iabca=iabca-1
      if (iabca.lt.0) go to 250
      if (ia.ne.iabca) go to 250
      if (ib.ne.iabcb) go to 250
      go to 90
   80 lkup=mult8(iabcs)+isym(levm)
      iabcs=lkupsm(lkup)
   90 nele=ib+ia+ia
      ifut=0
      it=iabct
      if (it.ge.8) go to 100
      ifut=it
      go to 110
  100 it=it-10
      ifut=it+1
  110 continue
      ipast=it
      if (lev.eq.ifermi+1) iexct=iop(2)
      if (nele.gt.(nlcs(levm)-it+iexct)) go to 250
      go to (220,220,220,220,120,150,200,140,130,180,190,220),icd
  120 continue
c
c  alpha orbital uncoupled
c
      if(k.gt.2) it=it+1
      go to 210
  130 continue
      if (k.eq.4) it=it+1
      if (k.gt.1) it=it+1
      go to 210
  140 continue
      if (k.eq.1) it=it-1
      if (k.lt.4) it=it-1
      if (ib.eq.2.and.k.eq.3) it=it+1
      ifut=0
      if ((ib.eq.0).and.(k.eq.2).and.(nrefs.le.2)) it=it+1
      if (it.lt.0) it=0
      go to 210
  150 continue
c
c  alpha orbital to be singlet coupled
c
      it=ifut
      if (k.lt.4) go to 160
      ipast=ipast+1
      it=it+1
      ifut=ifut+1
  160 if (k.le.1) go to 170
      ipast=ipast+1
      ifut=ifut+1
  170 continue
      if (k.eq.3.and.(ib-1-it).gt.nb) it=it+1
      if (ipast.lt.it) it=ipast
      go to 210
  180 if (k.ne.4) go to 250
      go to 210
  190 if (k.ne.1) go to 250
      go to 210
  200 continue
c
c  singlet coupled beta orbital
c
      if (k.eq.4) it=it+1
      if (k.eq.2) it=it+1
      ifut=it
  210 continue
      if (it.gt.iexct) go to 250
      iabct=it
      ifut=ifut-it
      if (ifut.lt.0) ifut=0
      if ((ifut.eq.1).and.(iabct.lt.8)) iabct=iabct+10
  220 continue
      if (levm.eq.(idocc+1).and.(iabct.lt.8)) iabct=0
      if (levm.eq.(idocc+1).and.(iabct.ge.8)) iabct=10
      if (levm.eq.(iact+1).and.(iabct.ge.8)) iabct=iabct-10
      nptm=levpt(levm)
      do 230 krm=1,nrm
      nptm=nptm+1
      if ((iabca.eq.nabca(nptm)).and.(iabcb.eq.nabcb(nptm)).and.(iabcs.e
     1q.nabcs(nptm)).and.(iabct.eq.nabct(nptm))) go to 240
  230 continue
      nrm=nrm+1
      if (nrm.gt.200) go to 380
      nptm=nrm+levpt(levm)
      if (nptm.gt.nrows) go to 390
      nabca(nptm)=iabca
      nabcb(nptm)=iabcb
      nabcs(nptm)=iabcs
      nabct(nptm)=iabct
      if (iabct.ge.8) nabct(nptm)=iabct-1
      krm=nrm
  240 continue
      iarpt=npt+ishift(k)
      if (iarc(iarpt).ne.0) go to 430
      iarc(iarpt)=krm
  250 continue
  260 continue
      levnr(levm)=nrm
      go to 30
  270 continue
c
c  set weight of array bottom to unity. let all other bottoms equal zero
c
      itest=0
      npt=levpt(1)
      nr=levnr(1)
      nrows=nr+npt
      do 280 i=1,nr
      npt=npt+1
      if ((nabca(npt).ne.0).or.(nabcb(npt).ne.0).or.(nabcs(npt).ne.1).or
     1.(nabct(npt).ne.0)) go to 280
      if (itest.ne.0) go to 400
      itest=1
      nlwks(npt)=1
  280 continue
      if (itest.eq.0) go to 410
c
c  generate weights of each row  (nlwks)
c   remove all rows with zero weight
c
      do 330 lev=2,nlevs
      levm=lev-1
      nrm=levnr(levm)
      nptm=levpt(levm)
      nrmx=0
      do 300 i=1,nrm
      nptm=nptm+1
      ired(i)=0
      if (nlwks(nptm).eq.0) go to 300
      nrmx=nrmx+1
      ired(i)=nrmx
      nptx=nrmx+levpt(levm)
      nlwks(nptx)=nlwks(nptm)
      nabca(nptx)=nabca(nptm)
      nabcb(nptx)=nabcb(nptm)
      nabcs(nptx)=nabcs(nptm)
      nabct(nptx)=nabct(nptm)
      do 290 k=1,4
      iarptx=nptx+ishift(k)
      iarpt=nptm+ishift(k)
      iarc(iarptx)=iarc(iarpt)
  290 continue
  300 continue
      levnr(levm)=nrmx
      npt=levpt(lev)
      nr=levnr(lev)
      nptm=levpt(levm)
      do 320 i=1,nr
      npt=npt+1
      do 310 k=1,4
      iarpt=npt+ishift(k)
      j=iarc(iarpt)
      if (j.eq.0) go to 310
      j=ired(j)
      iarc(iarpt)=j
      if (j.eq.0) go to 310
      iptm=j+nptm
      nlwks(npt)=nlwks(npt)+nlwks(iptm)
  310 continue
  320 continue
  330 continue
c
c  contract the paldus and weight arrays
c
      ipt=0
      lev=nlevs
  340 continue
      if (lev.eq.0) go to 360
      nr=levnr(lev)
      npt=levpt(lev)
      levpt(lev)=ipt
      do 350 i=1,nr
      ipt=ipt+1
      npt=npt+1
      nabca(ipt)=nabca(npt)
      nabcb(ipt)=nabcb(npt)
      nabcs(ipt)=nabcs(npt)
      nabct(ipt)=nabct(npt)
      nlwks(ipt)=nlwks(npt)
      do 350 k=1,4
      iarpt=ipt+ishift(k)
      iarptx=npt+ishift(k)
      iarc(iarpt)=iarc(iarptx)
  350 continue
      lev=lev-1
      go to 340
  360 continue
      nrows=ipt
      nwks=nlwks(1)
      go to 610
  380 write(iwr,450) levm
  600 call caserr(
     *'dimensioning error in drt initialisation')
  390 write(iwr,460) nrows
      go to 600
  400 write(iwr,470)
      go to 600
  410 write(iwr,480) nr
      npt=levpt(1)
      do 420 i=1,nr
      npt=npt+1
      write(iwr,490) nabca(npt),nabcb(npt),nabcs(npt),nabct(npt)
  420 continue
  430 continue
      write(iwr,500)
      go to 600
 610  return
c
  450 format (24h  too many rows at level,i5,13h  (max is 63))
  460 format (28h  too many rows.  maximum is,i8)
  470 format (23h  duplicate graph tails)
  480 format (27h  no valid tail.  levnr(1)=,i5)
  490 format (3x,7h nabc =/1h ,4i15)
  500 format (28h in paldus, iarc(iarpt).ne.0)
      end
      subroutine rdrsym (icd,nsymtp,ntr,nsym)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common/cone/kmna(9),lkupsm(64)
c
      dimension ntp(8), ntr(8), nsymtp(10,8)
      dimension mult8(8)
c
      data mult8/0,8,16,24,32,40,48,56/
c
      do 10 i=1,nsym
      ntp(i)=nsymtp(icd,i)
   10 continue
      ntr(1)=1
      ntp(1)=0
      do 50 itr=2,nsym
      if (itr.eq.4) go to 30
      iin=1
      imax=0
      do 20 i=2,nsym
      if (ntp(i).lt.imax) go to 20
      iin=i
      imax=ntp(i)
   20 continue
      go to 40
   30 continue
      lkupa=ntr(2)
      lkup=mult8(lkupa)+ntr(3)
      iin=lkupsm(lkup)
   40 continue
      ntr(itr)=iin
      ntp(iin)=-1
   50 continue
      return
      end
      subroutine reord (levfrm,npflg)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/gjs)
INCLUDE(common/harmon)
INCLUDE(common/iofile)
      common/junk/iorbtp(12,maxorb),ismtp(maxorb)
      common/lsort /norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     * ,jcode(maxorb),nlcs(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8),nbf
      common/valen / ivv,ifermi
      common/casscg/ncore,nact,naa,nbb,itypec(maxorb)
      common/drtlnk/norbl,nsyml,lenblo,iprinl,labdrt(maxorb)
INCLUDE(common/runlab)
      common/cone/kmna(9),lkupsm(64)
c
c  input
c       (1)label
c       (2)iops
c       (3)nsym,norbs,nrefs
c       (4)iref(i),isym(i),i=1,norbs
c   iref codes       format(16(a3,i1,1x))
c   doc-doubly occupied
c   uoc-unoccupied
c   fzc-frozen core
c   fzv-frozen virtual
c   alp-alpha occupied(high spin states)
c   aos-alpha of open shell singlet
c   bos-beta of  open shell singlet
c   gva-first of gbv pair
c   gvb-second of gvb pair
_IF1()      dimension iref(16)
      dimension norbtp(12), mult8(8)
      dimension nsymtp(12,8)
      dimension isymtr(8)
      dimension ylbl(7),yref(12)
c
      data mult8/0,8,16,24,32,40,48,56/
      data yref/'doc','uoc','fzc','fzv','alp','aos','bos','gva','gvb',
     1'cor','vir','val'/
      data ylbl/'fzc','doc','alp','aos','bos','uoc','fzv'/
c
      ncore=0
      nact=0
      naa=0
      nbb=0
c     jblank=0
      ntps=12
      ivv=0
      do 20 i=1,ntps
      do 10 j=1,8
10    nsymtp(i,j)=0
20    norbtp(i)=0
      ispin=0
      msym=1
      npflg=0
      call setsto(128,0,iout)
c
c ----- default specifications for iop arrays
c
      do 410 i=1,8
 410  iop(i)=0
      iop(1)=1
      iop(5)=1
      iop(2)=10
      nrefs=1
c
c ----- link with casscf though /drtlnk/
c
      norbl = newbas0
      norb=norbl
      nsym=nsyml
c     lenb=lenblo
      if(iprinl.eq.1)npflg=1
      write(iwr,400)ztitle
      if (nrefs.gt.3) ivv=1
      write(iwr,420) nsym,norb,nrefs
      if(npflg.ne.0)write(iwr,550)
 550  format(/' print of drt and configuration list requested'/)
      iop(8)=nrefs
      write(iwr,440)
      if (norb.gt.maxorb.or.norb.le.0) go to 380
      nbf=norb
      norbs=newbas0
      if (nsym.gt.8.or.nsym.le.0) go to 370
      k=0
      do 90 i=1,maxorb
      icode=locatc(yref,ntps,ylbl(labdrt(i)))
      if(icode.eq.0)go to 350
      k=k+1
      ism=isymmo(i)
      itypec(k)=ism
      if(icode.eq.1.or.icode.eq.5.or.icode.eq.6)naa=naa+1
      if(icode.eq.1.or.icode.eq.7)nbb=nbb+1
      if(icode.eq.3)ncore=ncore+1
      if(icode.lt.3.or.(icode.ge.5.and.icode.le.7))nact=nact+1
      write(iwr,460) k,ism,yref(icode)
      ntp=norbtp(icode)+1
      norbtp(icode)=ntp
      iorbtp(icode,ntp)=k
      nsymtp(icode,ism)=nsymtp(icode,ism)+1
      if (ism.le.0.or.ism.gt.nsym) go to 340
      ismtp(k)=ism
      if (k.ge.norbs) go to 100
90    continue
100   continue
      do 110 ism=1,nsym
      do 110 icode=5,9
      nsymtp(1,ism)=nsymtp(1,ism)+nsymtp(icode,ism)*256
110   continue
c
c  frozen core orbitals
c
      icode=3
      ntp=norbtp(icode)
      if (ntp.lt.1) go to 130
      do 120 j=1,ntp
      k=iorbtp(icode,j)
      iout(k)=-1
120   continue
130   iorb=0
      na=0
      nb=0
      nc=0
      idocc=0
      nele=0
      iact=0
c
c  unoccupied orbitals
c
      icode=2
      ntp=norbtp(icode)
      call rdrsym (icode,nsymtp,isymtr,nsym)
      if (ntp.lt.1) go to 150
      do 140 iism=1,nsym
      ism=isymtr(iism)
      do 140 j=1,ntp
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 140
      iorb=iorb+1
      iout(k)=iorb
      idocc=idocc+1
      iact=iact+1
      nc=nc+1
      nlcs(iorb)=nele
      isym(iorb)=ismtp(k)
      jcode(iorb)=icode
140   continue
150   levfrm=iorb+1
c
c  restricted virtual orbitals.
c
      icode=11
      ntp=norbtp(icode)
      if (ntp.lt.1) go to 170
      do 160 ism=1,nsym
      do 160 j=1,ntp
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 160
      iorb=iorb+1
      iout(k)=iorb
      nc=nc+1
      nlcs(iorb)=nele
      isym(iorb)=ism
      jcode(iorb)=icode
160   continue
170   continue
c
c      valance virtuals
c
      icode=12
      ntp=norbtp(icode)
      call rdrsym (icode,nsymtp,isymtr,nsym)
      if (ntp.lt.1) go to 190
      ivv=1
      do 180 iism=1,nsym
      iiism=nsym-iism+1
      ism=isymtr(iiism)
      do 180 jj=1,ntp
      j=ntp-jj+1
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 180
      iorb=iorb+1
      iout(k)=iorb
      iact=iact+1
      nc=nc+1
      nlcs(iorb)=nele
      isym(iorb)=ismtp(k)
      jcode(iorb)=icode
180   continue
190   continue
c
c  doubly occupied orbitals
c
      icode=1
      ntp=norbtp(icode)
      call rdrsym (icode,nsymtp,isymtr,nsym)
      if (ntp.lt.1) go to 210
      do 200 iism=1,nsym
      iiism=nsym-iism+1
      ism=isymtr(iiism)
      do 200 jj=1,ntp
      j=ntp-jj+1
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 200
      iorb=iorb+1
      iout(k)=iorb
      iact=iact+1
      na=na+1
      nlcs(iorb)=nele
      nele=nele+2
      isym(iorb)=ismtp(k)
      jcode(iorb)=icode
200   continue
c
c  frozen virtual orbitals
c
210   icode=4
      ntp=norbtp(icode)
      if (ntp.lt.1) go to 230
      do 220 j=1,ntp
      k=iorbtp(icode,j)
      iout(k)=0
220   continue
230   continue
      do 300 icode=5,9
      ntp=norbtp(icode)
      call rdrsym (icode,nsymtp,isymtr,nsym)
      if (ntp.lt.1) go to 300
      do 290 iism=1,nsym
      iiism=nsym-iism+1
      ism=isymtr(iiism)
      do 290 jj=1,ntp
      j=ntp-jj+1
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 290
      iorb=iorb+1
      iout(k)=iorb
      nlcs(iorb)=nele
      isym(iorb)=ismtp(k)
      jcode(iorb)=icode
      go to (330,330,330,330,240,250,260,270,280),icode
240   continue
      nb=nb+1
      nele=nele+1
      lkup=mult8(msym)+isym(iorb)
      msym=lkupsm(lkup)
      ispin=ispin+1
      go to 290
250   continue
      nc=nc+1
      if (j.ne.ntp) go to 290
      lkup=mult8(msym)+isym(iorb)
      msym=lkupsm(lkup)
      go to 290
260   continue
      nele=nele+2
      na=na+1
      lkup=mult8(msym)+isym(iorb)
      msym=lkupsm(lkup)
      go to 290
270   continue
      nele=nele+1
      nc=nc+1
      go to 290
280   continue
      nele=nele+1
      na=na+1
290   continue
300   continue
c
c  restricted doubly occupied orbitals
c
      icode=10
      ntp=norbtp(icode)
      if (ntp.lt.1) go to 320
      do 310 ism=1,nsym
      do 310 jj=1,ntp
      j=ntp-jj+1
      k=iorbtp(icode,j)
      if (ism.ne.ismtp(k)) go to 310
      iorb=iorb+1
      iout(k)=iorb
      na=na+1
      nlcs(iorb)=nele
      nele=nele+2
      isym(iorb)=ismtp(k)
      jcode(iorb)=icode
310   continue
320   continue
      if (iorb.gt.127) go to 360
      if (iorb.ne.na+nb+nc) go to 330
      idocc=idocc+iact*256
      nlevs=na+nb+nc+1
      ifermi=levfrm
      return
c
330   write(iwr,470)
      write(iwr,480)
      go to 610
340   write(iwr,490) i
       go to 610
350   write(iwr,500) k
       go to 610
360   write(iwr,510) iorb
       go to 610
 370  write(iwr,520) nsym
       go to 610
380   write(iwr,530) norb
      write(iwr,540)
 610   call caserr('error in guga preprocessor')
       return
c
400   format (' title for this job - ',10a8)
420   format (/1x,'number of symmetry types   is',i4/
     +         1x,'number of basis functions  is',i4/
     +         1x,'number of reference states is',i4)
440   format (1x,'orbital information'/)
460   format (8h orbital,i3,6h  sym ,i2,4x,a4)
470   format (20h orbitals dont match)
480   format (20h too many references)
490   format (8h orbital,i4,29h  has an illegal symetry type)
500   format (8h orbital,i4,29h  has an illegal active space)
510   format (i5,51h  is too many unfrozen orbitals. the maximum is  63)
520   format (i8,41h  is not a valid number of symmetry types)
530   format (i8,35h  is not a valid number of orbitals)
540   format (25h too many active orbitals)
      end
      subroutine symtri (ijadd,ijgrp,kadd,ladd,inext,jmnnxt,jmxnxt,ijnin
     1t,ningrp,nmax,nxt,ngrps,itape8,nij,nkl)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/iofile)
      common/cone/kmna(9),lkupsm(64)
      common /lsort/norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb)
     * ,icode(maxorb),nlcs(maxorb),
     1levpt(maxorb),levnr(maxorb),iout(maxorb),iop(8)
c
      dimension ijadd(nij), ijgrp(nij), kadd(nkl), ladd(nkl), ijnint(nij
     1)
      dimension inext(200), jmnnxt(200), jmxnxt(200)
      dimension ningrp(200)
      dimension inint(maxorb)
      dimension ismoff(8)
      dimension ybsym(2)
      dimension mult8(8)
c
      data mult8/0,8,16,24,32,40,48,56/
      data ybsym/'kadd','ladd'/
      data yasym/'sym '/
c
      n=norbs
c     irecp8=0
      do 10 ism=1,nsym
      ismoff(ism)=(ism-1)*norbs
      do 10 i=1,n
      ismpt=ismoff(ism)+i
      kadd(ismpt)=0
      ladd(ismpt)=0
   10 continue
c
c  set nmaxmn based on core size neede later
c
      nmaxmn=10000
      inmx=0
c
c  count the integrals
c
      ij=0
      ijad=0
      do 50 i=1,n
      jsm=isym(i)
      jad=0
      do 40 j=1,i
      ij=ij+1
      lkup=mult8(jsm)+isym(j)
      ksm=lkupsm(lkup)
      kad=0
      do 30 k=1,j
      ksmpt=k+ismoff(ksm)
      if (kadd(ksmpt).ne.0.and.kadd(ksmpt).ne.kad) go to 210
      kadd(ksmpt)=kad
      lkup=mult8(ksm)+isym(k)
      lsm=lkupsm(lkup)
      lad=0
      do 20 l=1,k
      lsmpt=l+ismoff(lsm)
      if (ladd(lsmpt).ne.0.and.ladd(lsmpt).ne.lad) go to 220
      ladd(lsmpt)=lad
      if (lsm.ne.isym(l)) go to 20
      if (l.eq.j.and.l.lt.i) go to 20
      lad=lad+3
      if (k.eq.l) lad=lad-1
   20 continue
      kad=kad+lad
   30 continue
      ijad=ijad+kad
      jad=jad+kad
      ijnint(ij)=kad
      if (kad.gt.nmax) go to 260
   40 continue
      inint(i)=jad
      if (jad.gt.inmx) inmx=jad
   50 continue
c
c  change nmax if it is too large
c
      if (inmx.gt.nmaxmn) nmaxmn=inmx
      if (nmaxmn.lt.nmax) nmax=nmaxmn
c
c  generate the addressing  (ijadd and ijgrp)
c
      ngrps=0
      i=n
      jmin=1
   60 continue
      k=0
      mleft=nmax
      if (ngrps.lt.1) go to 80
      do 70 igrp=1,ngrps
      nleft=nmax-ningrp(igrp)-inint(i)
      if (nleft.lt.0) go to 70
      if (nleft.gt.mleft) go to 70
      mleft=nleft
      k=igrp
   70 continue
   80 continue
      if (k.gt.0) go to 90
      ngrps=ngrps+1
      if (ngrps.gt.200) go to 250
      k=ngrps
      ningrp(k)=0
   90 continue
      ij=(i*(i-1))/2+jmin
      do 100 j=jmin,i
      if ((ningrp(k)+ijnint(ij)).gt.nmax) go to 110
      ijgrp(ij)=k
      ijadd(ij)=ningrp(k)
      ningrp(k)=ningrp(k)+ijnint(ij)
      ij=ij+1
  100 continue
      i=i-1
      jmin=1
      if (i.eq.0) go to 120
      go to 60
  110 continue
      jmin=j
      inint(i)=inint(i)-ningrp(k)
      go to 60
  120 continue
c
c  generate  inext and jnext arrays for formula generation
c
      nxt=0
      do 160 igrp=1,ngrps
      i=n
  130 continue
      itest=0
      ij=(i*(i-1))/2
      do 150 j=1,i
      ij=ij+1
      if (itest.ne.0) go to 140
      if (ijgrp(ij).ne.igrp) go to 150
      itest=1
      nxt=nxt+1
      if (nxt.gt.199) go to 240
      inext(nxt)=i
      jmnnxt(nxt)=j
      jmxnxt(nxt)=i
      go to 150
  140 continue
      if (ijgrp(ij).eq.igrp) go to 150
      jmxnxt(nxt)=j-1
      itest=0
  150 continue
      i=i-1
      if (i.gt.0) go to 130
  160 continue
      ijlst=(n*(n+1))/2
      nout=2*nsym*norbs+2*ijlst+1
c
c  if there is only one group reset nmax to its size
c
      if (ngrps.eq.1) nmax=ijad
      write (itape8) ngrps,nxt,nmax
      write (itape8) ijadd
      write (itape8) ijgrp
      write (itape8) kadd
      write (itape8) ladd
      write (itape8) inext
      write (itape8) jmnnxt
      write (itape8) jmxnxt
      write (itape8) icode
      write(iwr,270) ngrps,ijad
      write(iwr,280) nmax
      if ((iop(1)/8)*8.eq.iop(1)) go to 200
      write(iwr,290)
      write(iwr,300)
      do 180 k=1,nxt
      i=inext(k)
      jmn=jmnnxt(k)
      jmx=jmxnxt(k)
      ij=(i*(i-1))/2+jmn
      lad=ijadd(ij)
      igrp=ijgrp(ij)
      jad=0
      do 170 l=jmn,jmx
      jad=jad+ijnint(ij)
      ij=ij+1
  170 continue
      write(iwr,310) k,igrp,i,jmn,jmx,jad,lad
  180 continue
      write(iwr,320) (yasym,k,k=1,nsym)
      write(iwr,330) (ybsym,k=1,nsym)
      nout=nsym*norbs
      do 190 k=1,n
      write(iwr,340) k,(kadd(l),ladd(l),l=k,nout,norbs)
  190 continue
  200 continue
      call reardr(ijgrp,ijadd,kadd,ladd,nij,nkl,nmax,itape8)
      end file itape8
      return
  210 write(iwr,350)
      go to 400
 220  write(iwr,360)
      go to 400
 240  write(iwr,370)
      go to 400
 250  write(iwr,380)
      go to 400
  260 write(iwr,390)
  400 call caserr('error in drt construction')
c
  270 format(/  ' number of groups is',i5/  '  number of integrals is',
     *i9)
  280 format (/' nmax is',1x,i6)
  290 format (1h1,27hintegral grouping breakdown)
  300 format (/,50h   nxt  block    i   jmn   jmx     nints     ijadd)
  310 format (5i6,4i10)
  320 format (/,5x,8(a4,2x,i2,7x))
  330 format (5h iorb,8(a4,2x,a4,5x))
  340 format (i5,8(2i6,3x))
  350 format (25h symmetry problems with k)
  360 format (25h symmetry problems with l)
  370 format (17h nxt is too large)
  380 format (34h there are too many integral bloks)
  390 format (18h nmax is too small)
      return
      end
      subroutine calcg(iqq,nele,lin,lenb,ntri,mt0)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common /lsort/ lout(3,681),iout
      common /junkg / ibl5,maxt,ires,ipass,nteff,npass,nbuck,mloww
     +                ,mhi,jtri,master,lenbas,nwks,ilow,niqq
     +                ,nwb,ibuck,indblk
      common /three/ mark(1500),nwbuck(1500)
      common/bufb/nkk,mkk,gin(1)
c
      dimension nele(*),iqq(lenb,4,ntri),lin(4,*)
c
      data maxb/9999999/
c
      iout=0
c
      call stopbk
c
c...  loop over buckets
c
      do 140 ibuck=1,nbuck
      mhigh=min(mloww+ntri,mhi)
      mtri=mhigh-mloww
c     mlow=mloww+1
      do 20 j=1,mtri
      nele(j)=0
      do 20 i=1,4
      call setsto(lenbas,0,iqq(1,i,j))
 20   continue
c
c...  load up
      mkk=mark(ibuck)
      go to 60
 40   call srtin (lin)
      do 50 iword=1,nkk
      itri=lin(1,iword)-mloww
      nel=nele(itri)+1
      iqq(nel,1,itri)=lin(2,iword)
      iqq(nel,2,itri)=lin(3,iword)
      iqq(nel,3,itri)=lin(4,iword)
50    nele(itri)=nel
 60   if(mkk.ne.maxb)go to 40
c
c...  process each row in this bucket
      do 90 itri=1,mtri
      master=master+1
      nel=nele(itri)
      iout=iout+1
      lout(1,iout)=master
      lout(2,iout)=nel
      lout(3,iout)=0
      if (iout.eq.681) call putblk(mt0)
      do 70 i=1,nel
      iout=iout+1
      lout(1,iout)=iqq(i,1,itri)
      lout(2,iout)=iqq(i,2,itri)
      lout(3,iout)=iqq(i,3,itri)
      if (iout.eq.681) call putblk(mt0)
70    continue
90    continue
140   mloww=mhigh
      if (iout.gt.0) call putblk(mt0)
      return
      end
      subroutine citape(q,lword)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      integer   ed0,mt0
c
INCLUDE(common/machin)
INCLUDE(common/restar)
INCLUDE(common/iofile)
       common/drtlnk/norbl(2),npassi
INCLUDE(common/files)
      common /junkg / ibl5,maxt,ires,ipass,nteff,npass,nbuck,mloww
     >    ,mhi,ntri,master,lenbas,nwks,ilow,niqq
     >    ,nwb,ibuck,indblk
INCLUDE(common/blksiz)
c
      integer    mxsurd, mxwks
      parameter (mxsurd=3600, mxwks=80000)
      common/lsort/surd(mxsurd),nsurd,nbasis,ncore,nact,ntype,
     1              itype(maxorb),nwk,na,nb,isppp
c
      character*10 charwall
      dimension q(*)
c
      data m0/0/
c
      npass = max (npassi,1)
c
      ed0=m9tape(1)
      mt0=nttape(1)
      nav = lenwrd()
c
c.....get header block and copy to output stream
c
      indblk = lensec(mach(17))+m9blk(1)
      call rdedx(surd,mach(17),m9blk(1),ed0)
      call wrt3(surd,mach(17),ntblk(1),mt0)
      nwks=nwk
      niqq=(nwks-1)/nav+ 1
c
c     nav2=nav+nav
c
c.... ibl5 is the number of 8 byte words which can be got into a
c     sortfile buffer (nsz * block)
      ibl5=nsz*512-(2-nav/2)
c
c      buffer space
c
      ibuff = (ibl5*4)/nav
c      core
      i10 = 1
      i20 = i10 + niqq
      i30 = i20 + ibuff
c
      master=0
c
      call stopbk
c
      call diagel (q(i10),q(i20)  ,ed0,mt0)
      t=cpulft(1)
      write (iwr,55) t,charwall()
55    format(/' end of diagonal formulae at ',f8.2,' seconds',
     1        a10,' wall')
      write(iwr,50)lenbas
50    format(/' length of sort buffer =',i8)
c
c...  use long  integers for buffers in both sort and calc
c... for short integers(*2) replace nav with nav2 in next 4 lines
c... and modify integer *2 argument calls etc
      nwords=lword-niqq -ibuff
      ires = ((nwords*nav)/4)/ibl5
      nwordc = lword-niqq - ibuff
      maxt = ((nwordc*nav)/4)/lenbas
      if(ires.lt.1.or.maxt.lt.1)callcaserr('not enough store')
      if(ires.gt.1500)ires=1500
c
c...   determine min. no. of passes
      i=nwks-master
70    nteff=i/npass+1
      if(((nteff-1)/ires).lt.maxt)goto 80
      npass=npass+1
      goto 70
80    if(npass.gt.i)npass=i+1
      write(iwr,90)npass
90    format(' number of passes =',i4/)
c
      do 100 ipass=1,npass
      call sortg(q(i20),q(i10),ed0,ibl5)
      call calcg(q(i30),q(i10),q(i20),lenbas,ntri,mt0)
      t=cpulft(1)
100   write (iwr,101) ipass,t,charwall()
101   format(/' end of pass',i3,' at ',f8.2,' seconds',a10,' wall')
c
      call put(q,m0,mt0)
c
c ----- at present 1-section only permitted.
c
      mtfile=ntfile
      mttape(mtfile)=nttape(1)
      mtblk( mtfile)=ntblk(1)
      mtlast(mtfile)=iposun(mttape(1))
c
      write(iwr,102)
 102  format(/1x,'status of reordered formula tape'/
     * 1x,32('-')/)
      call filprn(mtfile,mtblk,mtlast,mttape)
c
      top=cpulft(1)
      write(iwr,10)top,charwall()
 10   format(/
     *' end of formulae tape construction and sorting at ',
     *f8.2,' seconds',a10,' wall'/)
      return
      end
      subroutine diagel (index,nele,ed0,mt0)
c
c.... scan ed0 i) accumulating no.   of matrix els for each row
c             ii) resolving diag els & writing to ed7
c   b) now ndiag is known, copy ed7 to mt0
c   c) evaluate sorting length lenbas
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
       integer ed0,ed7,mt0
c
INCLUDE(common/iofile)
INCLUDE(common/scra7)
      common /junk / nuwks(340),nlwks(340),iuwks(340),juwks(340)
     >               ,integs(340),icoup(340),mword,icode ,isppp(2)
      common /lsort/ lout(3,681),iout
      common /junkg / ibl5,maxt,ires,ipass,nteff,npass,nbuck,mloww
     >    ,mhi,ntri,master,lenbas,nwks,ilow,niqq
     >    ,nwb,ibuck,indblk
c
      dimension       index(*),nele(*)
c
      data ed7/8/
c
      nav = lenwrd()
      call readi(index,niqq*nav,indblk,ed0)
      call setsto(nwks,0,nele)
      nblks=0
      iout=1
      call search (ibl7la,ed7)
c
50    call getblk (ed0,nw)
      if (nw.eq.0) goto 101
      if (mword.eq.0) goto 50
      if (icode.eq.1) goto 51
      do 120 iword=1,mword
      nlwk = nlwks(iword)
      nuwk = nuwks(iword)
      iuwk = iuwks(iword)
      do 110 i=1,nlwk
      iii=index(iuwk)
      do 100 j=1,nuwk
      iout=iout+1
      lout(1,iout)=iii
      lout(2,iout)=integs(iword)
      lout(3,iout)=icoup (iword)
      nele(iii)=nele(iii)+1
      if (iout.lt.681) goto 100
      call putblk(ed7)
      iout=0
      nblks=nblks+1
100   iii=iii+1
110   iuwk=iuwk+1
120   continue
      goto 50
c
51    do 320 iword=1,mword
      nlwk = nlwks(iword)
      nuwk = nuwks(iword)
      iuwk = iuwks(iword)
      juwk = juwks(iword)
      do 310 i=1,nlwk
      iii=index(iuwk)
      jjj=index(juwk)
      do 300 j=1,nuwk
      nele(iii)=nele(iii)+1
      nele(jjj)=nele(jjj)+1
      iii=iii+1
300   jjj=jjj+1
      iuwk=iuwk+1
310   juwk=juwk+1
320   continue
      goto 50
c
101   ndiag=nblks*681+iout-1
      if (iout.eq.0) goto 102
      call putblk(ed7)
      nblks=nblks+1
102   call clredx
      write (iwr,103) ndiag,nblks
103   format(/' there are',i8,' diagonal formulae occupying'
     < ,i6,' blocks')
c
      call search (ibl7la,ed7)
      i2 = 2044/nav
      do 200 i=1,nblks
      call getblk (ed7,nw)
_IF1(c)      call fmove (nuwks,lout,i2)
_IFN1(c)      call icopy(i2*nav,nuwks,1,lout,1)
      if (i.ne.1) goto 200
      ndiag1 = ndiag/65536
      ndiag2 = ndiag - ndiag1*65536
      lout(1,1) = ndiag1
      lout(2,1) = ndiag2
      lout(3,1) = 0
200   call putblk(mt0)
c
      lenbas=0
      do 400 iii=1,nwks
400   lenbas=max(lenbas,nele(iii))
      return
      end
      subroutine sortg(ibuff,index,ed0,ibl5)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer ed0
c
      common /junk / nuwks(340),nlwks(340),iuwks(340),juwks(340)
     >               ,ints (340),icoup(340),mword,icode
      common/junkg/iblsrt ,maxt,ires,ipass,nteff,npass,nbuck,mloww
     >                ,mhi,ntri,master,lenbas,nwks,ilow,niqq
     >                ,nwb,ibuck,indblk
      common /three/ mark(1500),nwbuck(1500)
      common /bufb /nwbnwb,lnklnk,gout(1)
      common /stak/ btri,mlowww,nstack,iblock,mstack
INCLUDE(common/blksiz)
c
      dimension ibuff(4,ibl5,*),index(*)
c
      data maxb/9999999/
c
      nav = lenwrd()
      call readi(index,niqq*nav,indblk,ed0)
c
c...  base and limit rows for this pass
      mloww=master
      mhi=min(master+nteff,nwks)
      mtri=mhi-mloww
      mlow=mloww+1
c
c...  number of buckets needed%
c
      nbuck=ires
10    ntri=(mtri-1)/nbuck+1
      if (ntri.gt.maxt) goto 20
      nbuck = nbuck-1
      if (nbuck.gt.0) goto 10
20    nbuck=nbuck+1
      ntri=(mtri-1)/nbuck+1
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0   ,nwbuck)
      call stopbk
      iblock=0
c
50    call getblk(ed0,nw)
      if (nw.eq.0) goto 101
      if (mword.eq.0) goto 50
      do 230 iword=1,mword
      nlwk=nlwks(iword)
      nuwk=nuwks(iword)
      iuwk=iuwks(iword)
      juwk=juwks(iword)
      if (icode.ne.1) juwk=iuwk
      do 230 i=1,nlwk
      iii=index(iuwk)
      jjj=index(juwk)
      do 220 j=1,nuwk
      if (iii.lt.mlow .or. iii.gt.mhi) goto 170
      ibuck=(iii-mlow)*btri+1
      nwb=nwbuck(ibuck)+1
      ibuff(1,nwb,ibuck)=iii
      ibuff(2,nwb,ibuck)=jjj
      ibuff(3,nwb,ibuck)=ints(iword)
      ibuff(4,nwb,ibuck)=icoup(iword)
      if (nwb.eq.ibl5) call srtout(ibuff(1,1,ibuck))
      nwbuck(ibuck)=nwb
170   if (iii.eq.jjj .or. jjj.gt.mhi .or. jjj.lt.mlow) goto 210
      ibuck=(jjj-mlow)*btri+1
      nwb=nwbuck(ibuck)+1
      ibuff(2,nwb,ibuck)=iii
      ibuff(1,nwb,ibuck)=jjj
      ibuff(3,nwb,ibuck)=ints(iword)
      ibuff(4,nwb,ibuck)=icoup(iword)
      if (nwb.eq.ibl5) call srtout(ibuff(1,1,ibuck))
      nwbuck(ibuck)=nwb
210   iii=iii+1
220   jjj=jjj+1
      iuwk=iuwk+1
230   juwk=juwk+1
      goto 50
c
101   do 280 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if (nwb.gt.0) call srtout(ibuff(1,1,ibuck))
280   continue
      return
      end
      subroutine getblk (iunit,nw)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common /junk / iin(2044)
      common /blkin / iy(1)
c
      call find (iunit)
      call get  (iy,nw)
      if (nw.eq.0) return
_IF(ibm,vax)
      call upakin(iy,iin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,iin,2044)
_ELSE
      call upakft(iy,16,iin,2044)
_ENDIF
      return
      end
_IFN(cray,t3e,i8,ibm,vax)
      subroutine upakft(in4,nbits,iout4,nw)
      implicit REAL (a-h,o-z)
c
      integer *2 in4,iout4
      dimension iout4(2,*),in4(*)
c
      if(nbits.ne.16)
     +  call caserr('invalid call to upakft')
      call setsto(nw,0,iout4)
      ind = 1
      do loop=1,nw/4

_IF(littleendian)
       iout4(1,ind  )=in4(ind  )
       iout4(1,ind+1)=in4(ind+1)
       iout4(1,ind+2)=in4(ind+2)
       iout4(1,ind+3)=in4(ind+3)
_ELSE
       iout4(2,ind  )=in4(ind  )
       iout4(2,ind+1)=in4(ind+1)
       iout4(2,ind+2)=in4(ind+2)
       iout4(2,ind+3)=in4(ind+3)
_ENDIF
       ind =ind +4
      enddo
      return
      end
_ENDIF
      subroutine putblk (iunit)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common/lsort/iout(2044)
      common/blkin/iy(1)
c
      data m511/511/
c
_IF1(iv)      call pakin(iout,iy,511)
_IFN1(iv)      call pack(iy,16,iout,2044)
      call put (iy,m511,iunit)
      iout(2044)=0
      return
      end

      subroutine srtout (ibuff)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common /junkg / ibl5,maxt,ires,ipass,nteff,npass,nbuck,mloww
     >                ,mhi,ntri,master,lenbas,nwks,ilow,niqq
     >                ,nwb,ibuck,indblk
      common /three/ mark(1500),nwbuck(1500)
INCLUDE(common/blksiz)
      common /bufb  / nwbnwb,lnklnk,iyout(1)
      common /stak / btri,mlowww(2),iblock,mstack
c
      dimension ibuff(*)
c
      call stopbk
      nwbnwb=nwb
      lnklnk=mark(ibuck)
_IF1(iv)      call pakin(ibuff,iyout,ibl5)
_IFN1(iv)      call pack(iyout,16,ibuff,ibl5*4)
      call sttout
      mark(ibuck)=iblock
      nwb=0
      iblock=iblock+nsz
      return
      end
      subroutine srtin (lin)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common/junkg/ibl5
      common /bufb  / nkk,mkk,iyout(1)
c
      dimension lin(*)
      mkkkk = mkk
      call rdbak(mkkkk)
      call stopbk
_IF1(iv)      call upakin(iyout,lin,ibl5)
_IFN1(iv)      call unpack(iyout,16,lin,ibl5*4)
      return
      end
      subroutine ed0out
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
      common/blkin/dum(511)
      common/craypk/ijkl(1)
      common/disc/isel(5),ipos(maxlfn)
INCLUDE(common/restar)
c
      dimension iyy(2)
      equivalence (iyy(1),dum(1))
      data m511/511/
c
      if(ipos(mainp).ne.iblkmp)callsearch(iblkmp,mainp)
_IF1(iv)      call pakin(ijkl,iyy,m511)
_IFN1(iv)      call pack(iyy,16,ijkl,2044)
      call put(dum,m511,mainp)
      iblkmp=iblkmp+1
      mblp=mblp+1
      if(mblp)41,40,41
c
c...   change channel
c
 40    m9file=mfilep
       m9tape(m9file)=n9tape(mfilep)
       m9blk(m9file)=n9blk(mfilep)
       m9last(m9file)=iblkmp
       mfilep=mfilep+1
      if(mfilep.gt.n9file)call caserr(
     *'insufficient sections allocated to loop-formula tape')
      mainp=n9tape(mfilep)
      iblkmp=n9blk(mfilep)
      mblp=iblkmp-n9last(mfilep)
  41  return
      end
c ******************************************************
c ******************************************************
c             =   casnr    =
c ******************************************************
c ******************************************************
      subroutine newton(q,ic5e,ic6e,skv,brill,v, icf,work,u, h,kk100
     * ,ccc , fock, gam1, gam2, nconf,ireturn)
c                        2     14  1    15 205 206 207 208 208
c
c.....this routine does the newton-raphson stage in the casscf program
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
INCLUDE(common/iofile)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/popnos)
INCLUDE(common/exc)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/gjs)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     +            nst,n1e,isym,mults,macro,iduqpr(4)
INCLUDE(common/ciconv)
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10ps
INCLUDE(common/prints)
      common/tab3  /is1,is1p,is2,is2p
INCLUDE(common/finish)
c
      dimension temp(maxorb)
      dimension skv(*),icf(*),u(length,*),ccc(nconf,*),q(*)
      dimension fock(*), gam1(*), gam2(*), work(ma,*), h(*),
     +          ic5e(*),ic6e(*),v(*),brill(*)
      character*10 charwall
c
      data vk20/.113d0/
c
      iunit1 = ntapes(1)
      iunit2 = ntapes(2)
c     sq2= dsqrt(2.0d0)
      kk10p=kk100
c     kkk10=kk10p-1
      t=cpulft(1)
      if(nprint.ne.-5)write(iwr,191)t ,charwall()
 191  format(/1x,53('*')//
     *' commence newton raphson procedure at ',f8.2,' seconds'
     *,a10,' wall')
      kk10s=kk10ps-1
      kk10=kk10p-1
c
      call gethes(h,kk10p,q,work,ic5e,ic6e,skv,brill,icf,kk10ps,work
     * , fock, gam1, gam2)
c
      call bldhes(h,kk10p,work,ic5e,ic6e,skv,brill,icf,kk10s,
     *   fock, gam1, gam2,nprint)
c
      if(isimul(macro).ne.0)call rdedx(brill,kk10,m4blk(1),m4tape(1))
      if(mode.lt.2)goto 50
      rewind iunit1
      ibase=1
      do 10 i=1,kk10s
      call wtfors(h(ibase),i,iunit1)
10    ibase=ibase+i
c
      if(isimul(macro).eq.0)goto50
c
c.....copy ci vector into hessian form
c
      nnn=kk10ps
      do 20 nn=1,kk3
      if(nn.eq.idom(kkkm))goto20
      h(nnn)=ccc(nn,kkkm)
      nnn=nnn+1
20    continue
c
c.....get gradient from disc
c
      rewind iunit2
      do 40 i=kk10ps,kk10
      call rdfors(v,i,iunit2)
      do 30 j=kk10ps,i
30    v(j)=v(j)-2.0d0*(brill(i)*h(j)+brill(j)*h(i))
40    call wtfors(v,i,iunit1)
      endfile iunit1
c
c
 50   continue
_IF1(v)      test=0.0d0
_IF1(v)      do 60 kk=1,kk10
_IF1(v)      tester= dabs(brill(kk))
_IF1(v)   60 if(test.lt.tester)test=tester
_IFN1(v)      kk=idamax(kk10,brill,1)
_IFN1(v)      test=dabs(brill(kk))
      s=dnrm2(kk10,brill,1)
      if(nprint.ne.-5)write(iwr,70)test,s
70    format(/' maximum first derivative =',f10.7
     1         ,'      norm of gradient =',f10.7/)
      oconv=test.lt.cccnv
      fmax = test
c
      if(oconv)go to 180
      if(iam(macro).ne.0)go to 121
      if(ipople.ne.0) go to 100
      if(ifhes.eq.1)
     1callcleski(h,kk10,v,v(kk10p),ibuffn,iunit1,iunit2,nprint)
      if(ifhes.eq.1.and.oprint(8))call writel(h,kk10)
c
      call baksub(h,kk10,brill,v,ibuffn,iunit2,nprint)
c
      go to 111
c
 100  call pople (h,ibuffn,brill,kk10,v,v(kk10p),nprint)
c
 111  continue
      s=dnrm2(kk10,brill,1)
      if(ifhes.eq.1.and.s.gt.vk20)
     +      call modhes(h,ibuffn,brill,kk10,v,v(kk10p),kk10s,nprint,
     +      ireturn)
      if (ireturn.gt.0) return
      go to 131
 121  call augmnt(h,ibuffn,brill,kk10,v,v(kk10p+1),nprint)
c
c.....x the solution of the n-r eqns - find u=exp(x) the orbital rotns.
c.....optimise orbitals for each symmetry type in turn
c
cc      if(cont(8))write(iwr,*)(brill(i),i=1,kk10)
 131  continue
      s=dnrm2(kk10,brill,1)
      if(nprint.ne.-5)write(iwr,110)s
110   format(/' norm of newton raphson solution =',f10.7)
      do 140 iset=1,ntype
      nthis=nsymm(iset)
      if(nthis.eq.0)go to 140
      if(oprint(4))write(iwr,120)iset
120   format(' symmetry',i2)
c
      lensq=length*length
      call vclr(u,1,lensq)
      do 130 kk=1,kk10s
      ip=ic5e(kk)
      if(isymmo(ip).ne.iset)goto130
      ip=mapsym(ip,iset,nblock)
      iq=mapsym(ic6e(kk),iset,nblock)
      u(ip,iq)=brill(kk)
      u(iq,ip)=-brill(kk)
130   continue
      nterm=4
      mark=2
      mmmm=nthis
      call expmat(u,length,nthis,work,temp,nterm)
c
      if(oprint(4)) call dout(mmmm,mark,u,length,iwr)
      if(.not.oconv)call trnsfo(q,work,u,iset,length)
c
140   continue
      if(isimul(macro).ne.1.or.oconv)goto180
      im=kk10s
      do150i=1,kk3
      if(i.eq.idom(kkkm))goto150
      im=im+1
      ccc(i,kkkm)=ccc(i,kkkm)-brill(im)
150   continue
      if(oprint(8))write(iwr,160)
160   format(//' improved ci vector'/1x,18('-')/)
      if(oprint(8))write(iwr,170)(ccc(i,kkkm),i=1,kk3)
170   format(5f23.15)
180   t=cpulft(1)
      if(nprint.ne.-5)write(iwr,190)t ,charwall()
190   format(/' end of newton-raphson procedure at',f8.2,
     1' seconds',a10,' wall')
      return
      end
      subroutine bldhes(h,kk10p,iwork,ic5e,ic6e,skv,brill,icf,kk10s,
     * qq,gam1,gam2,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/qice)
      common/block /length
INCLUDE(common/exc)
INCLUDE(common/simul)
      common/dims  /kk1,kk2,kk3
INCLUDE(common/iofile)
INCLUDE(common/prints)
INCLUDE(common/ctrl)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1              nst,n1e,isym,mults,macro
INCLUDE(common/gjs)
c
      character*10 charwall
      dimension brill(*),h(*),ic5e(*),ic6e(*),skv(*),icf(*),iwork(ma,*)
      dimension qq(*),gam1(*),gam2(*)
c
      data half,two/0.5d0,2.0d0/
c
      iunit = ntapes(1)
c
      do 10 i=1,kk10p
10    icf(i)=i*(i-1)/2
      icf(kk10p+1)=kk10p*(kk10p+1)/2
c     lentri=icf(kk10p)
      kk10=kk10p-1
c
c.....now finish off gradients, then hessian
c.....<0|h|i-a> elements
c
      if(ncore.eq.0)go to 40
      do 30 i=1,ncore
      do 20 ia=nprimp,ma
      if(isymmo(i).ne.isymmo(ia))go to 20
      m1=ic1e(ia,1)+ic1e(i,2)
      kk=iwork(i,ia)
      brill(kk)=(skv(n1e+m1)+skv(m1))*two
20    continue
30    continue
c
c.....<0|h|t-a> elements
c
40    do 60 it=nst,nprim
      do 60 ia=nprimp,ma
      if(isymmo(it).ne.isymmo(ia))go to 60
      kk=iwork(it,ia)
      do 50 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(it)) go to 50
      m1=ic1e(ia,1)+ic1e(iu,2)
      m2=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      temp=skv(m1)*gam1(m2)
      if(it.ne.iu)temp=temp*half
      brill(kk)=brill(kk)+temp
50    continue
60    continue
c
c.....<0|h|i-t> elements
c
      if(ncore.eq.0)go to 90
      do 80 i=1,ncore
      do 80 it=nst,nprim
      if(isymmo(i).ne.isymmo(it))go to 80
      kk=iwork(i,it)
      m1a=ic1e(it,1)+ic1e(i,2)
      brill(kk)=brill(kk)+two*skv(n1e+m1a)
      do 70 iu=nst,nprim
      if(isymmo(it).ne.isymmo(iu)) go to 70
      m1=ic1e(iu,1)+ic1e(i,2)
      m2=ic1e(max(it,iu),1)+ic1e(min(it,iu),2)
      if(it.eq.iu)brill(kk)=brill(kk)+skv(m1)*(two-gam1(m2))
      if(it.ne.iu)brill(kk)=brill(kk)-skv(m1)*gam1(m2)*half
70    continue
80    continue
90     continue
      if(nprint.eq.-5)go to 120
      if(.not.oprint(12))goto120
      write(iwr,100)
100   format(//1x,'first derivatives'/1x,17('-')/)
      write(iwr,110)(brill(kk),kk=1,kk10)
110   format(1x,10f13.6)
120   continue
c-----------------------------------------------------------------------
c
c.....additions to hessian now
      if(mode.ne.2)goto290
      do 260 kk=1,kk10s
      ip=ic5e(kk)
      iq=ic6e(kk)
      ind=icf(kk)
      if(ip.gt.ncore)goto170
      if(iq.gt.nprim)goto230
c
c......h(it,ju)..
c
      do 160 ll=1,kk
      ir=ic5e(ll)
      is=ic6e(ll)
      ind=ind+1
      if(isymmo(ip).ne.isymmo(ir))goto160
      in1=ic1e(max(iq,is),1)+ic1e(min(iq,is),2)
      in2=ic1e(ip,1)+ic1e(ir,2)
      dumx=skv(in2)*gam1(in1)
      if(iq.ne.is)dumx=dumx*half
      if(ip.ne.ir)goto140
      in3=(iq-1)*nprim+is
      in4=(is-1)*nprim+iq
      dumx=(skv(in1)+skv(in1+n1e))*two+dumx-qq(in3)-qq(in4)
      do 130 iv=nst,nprim
      if(isymmo(iv).ne.isymmo(iq))goto130
      in3=ic1e(max(iq,iv),1)+ic1e(min(iq,iv),2)
      in4=ic1e(max(is,iv),1)+ic1e(min(is,iv),2)
      dumy=gam1(in4)*skv(in3)*half
      if(iv.ne.is)dumy=dumy*half
      dumx=dumx-dumy
      dumy=gam1(in3)*skv(in4)*half
      if(iv.ne.iq)dumy=dumy*half
      dumx=dumx-dumy
130   continue
140   if(iq.ne.is)goto150
      dumx=dumx-two*(skv(in2)+skv(in2+n1e))
150   h(ind)=h(ind)+dumx
160   continue
      goto260
c
c
c
170   do 220 ll=1,kk
      ir=ic5e(ll)
      is=ic6e(ll)
      ind=ind+1
      if(ir.gt.ncore)goto190
      if(is.gt.nprim)goto180
      if(isymmo(ip).ne.isymmo(is))goto220
      in1=ic1e(iq,1)+ic1e(ir,2)
      in2=ic1e(max(ip,is),1)+ic1e(min(ip,is),2)
      dumx=-gam1(in2)*skv(in1)*half
      if(ip.eq.is)dumx=dumx*two+skv(in1)+skv(in1+n1e)
      h(ind)=h(ind)+dumx
      goto220
c
c.....h(ta,jb)...
c
180   if(iq.ne.is)goto220
      in1=ic1e(ip,1)+ic1e(ir,2)
      in2=iwork(ir,ip)
      dumx=half*brill(in2)-(skv(in1)+skv(in1+n1e))*two
      h(ind)=h(ind)+dumx
      goto220
c
c
c.....h(ta,ub)..
c
190   if(isymmo(ip).ne.isymmo(ir))goto220
      in1=ic1e(ip,1)+ic1e(ir,2)
      in2=ic1e(max(iq,is),1)+ic1e(min(iq,is),2)
      dumx=-gam1(in1)*skv(in2)
      if(ip.eq.ir)dumx=dumx*two
      if(iq.ne.is)goto210
      in1=(ip-1)*nprim+ir
      in3=(ir-1)*nprim+ip
      dumx=dumx+two*(qq(in1)+qq(in3))
      do 200 iv=nst,nprim
      if(isymmo(iv).ne.isymmo(ip))goto200
      in1=ic1e(max(iv,ip),1)+ic1e(min(iv,ip),2)
      in2=ic1e(max(iv,ir),1)+ic1e(min(iv,ir),2)
      dumx=dumx+skv(in1)*gam1(in2)+skv(in2)*gam1(in1)
      if(iv.ne.ip)dumx=dumx-skv(in2)*gam1(in1)*half
      if(iv.ne.ir)dumx=dumx-skv(in1)*gam1(in2)*half
200   continue
210   continue
      h(ind)=h(ind)-half*dumx
220   continue
      goto260
c
c
c
230   do 250 ll=1,kk
      ind=ind+1
      ir=ic5e(ll)
      is=ic6e(ll)
      if(is.gt.nprim)goto240
c
c.....h(ia,jt)
c
      if(ip.ne.ir)goto250
      in1=ic1e(iq,1)+ic1e(is,2)
      in2=iwork(is,iq)
      h(ind)=h(ind)+two*(skv(in1)+skv(in1+n1e))-half*brill(in2)
      goto250
c
c.....h(ia,jb)..
c
240   if(isymmo(ip).ne.isymmo(ir))goto250
      in1=ic1e(max(iq,is),1)+ic1e(min(iq,is),2)
      in2=ic1e(ip,1)+ic1e(ir,2)
      if(ip.eq.ir)h(ind)=h(ind)+two*(skv(in1+n1e)+skv(in1))
      if(iq.eq.is)h(ind)=h(ind)-two*(skv(in2+n1e)+skv(in2))
250   continue
c
260   continue
c
c.....fock operator
c
 290  call dpconv (nprim,gam1,gam2)
      do 280 i=nst,nprim
      i1=ilifp(i)
      isymmoi=isymmo(i)
      do 280 j=nst,nprim
      if(isymmoi.ne.isymmo(j)) goto 280
      buff=0.0d0
      do 270 k=nst,nprim
      if(isymmoi.eq.isymmo(k)) buff=buff +
     1 skv(ic1e(max(j,k),1)+ic1e(min(j,k),2)) *
     2 gam1(ic1e(max(i,k),1)+ic1e(min(i,k),2) )
270   continue
      qq(i1+j)=qq(i1+j)+buff*0.5d0
280   continue
      if(oprint(8).and.mode.eq.2)call writel( h,kk10)
      call clredx
      if(nprint.eq.-5)go to 330
      t=cpulft(1)
      if(mode.eq.1)write(iwr,300)t ,charwall()
      if(mode.eq.2)write(iwr,310)t ,charwall()
300   format(' first derivative evaluation complete at',f8.2,
     *' seconds',a10,' wall')
310   format(/' first and second derivative evaluation complete at'
     1,f8.2,' seconds',a10,' wall')
 330  ifhes=0
c
c     kkkk=nprim*nprim
      call dpconv(nprim,gam1,gam2)
      if(mode.eq.1)return
      if(oprint(8))write(iwr,320)
320   format(/40x,'hessian matrix'/40x,14('-')/)
      if(oprint(8))call writel(h,kk10)
      ifhes=1
      rewind iunit
c
      return
      end
_EXTRACT(gethes,gethes,pclinux)
      subroutine gethes(h,kk10p,q,iwork,
     +     ic5e,ic6e,skv,brill,iky,kk10ps,work,qq,gam1,gam2)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/craypk/i205(1360)
INCLUDE(common/iofile)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
     1                          ,lenext(8)
INCLUDE(common/prints)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/filel)
INCLUDE(common/restar)
      common/disc/isel,iselr,iselw,irep,ichek,ipos(maxlfn)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     +            nst,n1e,isym,mults,iduqpar(5)
      common/blkin /g(510),nword
INCLUDE(common/atmblk)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/gjs)
c
      character*10 charwall
      dimension work(*),iky(*),q(*),ij(2),qq(*),gam1(*),gam2(*)
      dimension brill(*),h(*),ic5e(*),ic6e(*),skv(*),iwork(ma,*)
c
c
      data half,donem,two,fourm/.5d0,-1.0d0,2.0d0,-4.0d0/
      data eight,quartm/8.0d0,-0.25d0/
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      t=cpulft(1)
      if(mode.ne.2.and.nprint.ne.-5)write(iwr,10)t ,charwall()
      if(mode.eq.2.and.nprint.ne.-5)write(iwr,20)t ,charwall()
c     nblkq=ma*ma
      kk10=kk10p-1
      kk10s=kk10ps-1
      lentri=ma*(ma+1)/2
      len2=lentri+lentri
      call rdfock(q,work,work(len2+1),skv,gam1)
      m2=2
      if(oprint(8))call fout(skv,gam1,m2,iwr)
      lentri=kk10s*(kk10s+1)/2
      do 30 in1=1,kk10p
30    iky(in1)=in1*(in1-1)/2
      iky(kk10p+1)=kk10p*(kk10p+1)/2
_IF1(civ)      call szero(h,lentri)
_IFN1(civ)      call vclr(h,1,lentri)
c
      do 40 kk=1,kk10s
      ic5e(kk)=ic5e(kk+1)
      ic6e(kk)=ic6e(kk+1)
40    iwork(ic5e(kk),ic6e(kk))=kk
      call vclr(brill,1,kk10)
c
      nprim2=nprim*nprim
      call vclr(qq,1,nprim2)
      lfile=m6file
      do 670 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i)  =m6blk(i)
 670  llblk(i)  =liblk(i)-m6last(i)
c
      do 680 ii=1,lfile
      idev=lotape(ii)
      call search(liblk(ii),idev)
      call find(idev)
      jj=llblk(ii)
50    jj=jj+1
      call get(g(1),nw)
      if(nw)680,680,60
60    if(jj.ne.0)call find(idev)
_IF1(iv)      call upak8v(g(num2e+1),i205)
_IFN1(iv)      call unpack(g(num2e+1),lab816,i205,numlab)
_IFN1(iv)      int4=1
      do 660 kk=1,nword
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
_IFN1(iv)      int4=int4+4
      ggg=g(kk)
c
      if(k.gt.nprim)goto70
      if(j.gt.nprim)goto160
      if(i.gt.nprim)goto260
      if(l.gt.ncore)goto440
      if(k.le.ncore)goto470
      if(j.gt.ncore)goto450
      goto540
70    if(j.gt.ncore.and.l.gt.ncore)goto 90
      if(j.gt.ncore)goto 120
      if(l.gt.ncore)goto 130
c
c.....(ai|bj) ints.
c
      if(isymmo(j).ne.isymmo(i).or.isymmo(k).ne.isymmo(l))goto80
      in=ind(iwork(j,i),iwork(l,k))
      h(in)=h(in)+eight*ggg
80    if(isymmo(j).ne.isymmo(k).or.isymmo(i).ne.isymmo(l))goto660
      in=ind(iwork(l,i),iwork(j,k))
      h(in)=h(in)-(ggg+ggg)
      goto660
c
c.....(at|bu) ints.
90    isymvx=mult(isymmo(j),isymmo(l))
      do 110 it=nst,nprim
      oitl=it.eq.l
      oitj=it.eq.j
      oitii=isymmo(it).ne.isymmo(i)
      oitik=isymmo(it).ne.isymmo(k)
      do 110 iu=nst,it
      oiui=isymmo(iu).ne.isymmo(i)
      oiuk=isymmo(iu).ne.isymmo(k)
      if(mult(isymmo(it),isymmo(iu)).ne.isymvx)goto110
      if(oitii) go to 100
      in1=icfcor(max(it,j))+min(it,j)
      in2=icfcor(max(iu,l))+min(iu,l)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      dumy=dumx
      if(oitj   )dumx=dumx+dumx
      if(iu.ne.l)dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      if(oiuk) go to 100
      in3=iky(iwork(it,i))+iwork(iu,k)
      h(in3)=h(in3)+dumx
c
100   if(it.eq.iu.or.oitik.or.i.eq.k)go to 110
      in1=icfcor(max(it,l))+min(it,l)
      in2=icfcor(max(iu,j))+min(iu,j)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      dumy=dumx
      if(oitl   )dumx=dumx+dumx
      if(iu.ne.j)dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      if(oiui) go to 110
      in3=iky(iwork(it,k))+iwork(iu,i)
      h(in3)=h(in3)+dumx
110   continue
c
      if(i.ne.k.or.j.le.l)goto660
c
c......recycle (av|ax) as (ax|av)
c
      in1=j
      j=l
      l=in1
      goto90
c
c
c.....(at|bi) ints.
120   k1=i
      l1=j
      i=k
      j=l
      k=k1
      l=l1
c
c.....(ai|bt)ints.
c
130   oij=isymmo(i).ne.isymmo(j)
      ojk=isymmo(j).ne.isymmo(k)
      do 150 it=nst,nprim
      if(isymmo(it).ne.isymmo(l))goto150
      dumx=gam1(ic1e(max(it,l),1)+ic1e(min(it,l),2))*(ggg+ggg)
      if(it.eq.l)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(k).or.oij)goto140
      in=iky(iwork(it,k))+iwork(j,i)
      h(in)=h(in)+dumx
140   if(isymmo(it).ne.isymmo(i).or.ojk)goto150
      in=iky(iwork(it,i))+iwork(j,k)
      h(in)=h(in)+dumx*quartm
150   continue
      goto660
c
160   if(l.gt.ncore)goto180
      if(k.gt.ncore)goto210
c
c.....(ab|ij) integrals
c
      ojk=isymmo(j).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      if(isymmo(i).ne.isymmo(k).or.ojl)goto170
      in1=iky(iwork(k,i))+iwork(l,j)
      h(in1)=h(in1)-(ggg+ggg)
170   if(isymmo(i).ne.isymmo(l).or.k.eq.l.or.i.eq.j.or.ojk)
     +   go to 660
c
c
c.....(ab|tu) integrals
      in2=iky(iwork(k,j))+iwork(l,i)
      if(in1.ne.in2)h(in2)=h(in2)-(ggg+ggg)
      goto660
180   isymvx=mult(isymmo(k),isymmo(l))
      ivx=iky(k-ncore)+l-ncore
      do 200 it=nst,nprim
      do 200 iu=nst,it
      isymtu=mult(isymmo(it),isymmo(iu))
      if(isymtu.ne.isymvx)goto200
      oiui=isymmo(iu).ne.isymmo(i)
      oiuj=isymmo(iu).ne.isymmo(j)
      itu=iky(it-ncore)+iu-ncore
      in=ic2e(max(itu,ivx))+ic3e(min(itu,ivx))
      dumx=gam2(in)*ggg
      dumy=dumx
      if(itu.ne.ivx)dumx=dumx*half
      if(it.eq.iu)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(i).or.oiuj)goto190
      in=ind(iwork(iu,j),iwork(it,i))
      h(in)=h(in)+dumx
190   if(isymmo(it).ne.isymmo(j).or.i.eq.j.or.it.eq.iu.or.oiui)
     +   go to 200
      in=ind(iwork(it,j),iwork(iu,i))
      h(in)=h(in)+dumx
200   continue
      goto660
c
c
c.....(ab|ti) integrals
210   continue
      if(isymmo(l).ne.isymmo(i))goto230
      in1=iwork(l,i)
      do 220 it=nst,nprim
      if(isymmo(it).ne.isymmo(j))goto220
      in2=iwork(it,j)
      in=iky(in2)+in1
      dumx=ggg*gam1(ic1e(max(k,it),1)+ic1e(min(k,it),2))
      if(it.ne.k)dumx=dumx*half
      h(in)=h(in)-dumx
220   continue
230   if(i.eq.j.or.isymmo(l).ne.isymmo(j))goto250
      in1=iwork(l,j)
      do 240 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto240
      in2=iwork(it,i)
      in=iky(in2)+in1
      dumx=ggg*gam1(ic1e(max(k,it),1)+ic1e(min(k,it),2))
      if(it.ne.k)dumx=dumx*half
      h(in)=h(in)-dumx
240   continue
250   continue
      goto660
c
260   if(l.gt.ncore)goto 370
      if(j.le.ncore)goto310
      if(k.le.ncore)goto340
c
c.....(at|ui) integrals
c
      isymvx=mult(isymmo(j),isymmo(k))
      do 280 it=nst,nprim
      oitk = it .ne. k
      if(isymmo(it).ne.isymmo(l))goto280
      do 270 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(i))goto270
      in1=icfcor(max(it,k))+min(it,k)
c
c......(ai|tj) integrals
      in2=icfcor(max(iu,j))+min(iu,j)
      dumx=ggg*gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))
      if(oitk   )dumx=dumx*half
      if(iu.eq.j)dumx=dumx+dumx
      if(in1.eq.in2)dumx=dumx+dumx
      in1=iky(iwork(iu,i))+iwork(l,it)
      h(in1)=h(in1)-dumx
270   continue
280   continue
c
      olj=isymmo(j).ne.isymmo(l)
      olk=isymmo(k).ne.isymmo(l)
      do 300 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto300
      in1=iky(iwork(it,i))
      if(isymmo(it).ne.isymmo(k).or.olj)goto290
      in2=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      in3=iwork(l,j)+in1
      dumx=gam1(in2)*ggg
      if(it.ne.k)dumx=dumx*half
      h(in3)=h(in3)-dumx
c
290   if(isymmo(it).ne.isymmo(j).or.olk)goto300
      dumx=gam1(ic1e(max(it,j),1)+ic1e(min(it,j),2))*(ggg+ggg)
      in3=iwork(l,k)+in1
      if(it.eq.j)dumx=dumx+dumx
      h(in3)=h(in3)+dumx
300   continue
      goto660
310   if(k.le.ncore)goto660
      oij=isymmo(i).ne.isymmo(j)
      oil=isymmo(i).ne.isymmo(l)
      do 330 it=nst,nprim
      if(isymmo(it).ne.isymmo(k))goto330
      in1=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      dumx=gam1(in1)*half
      if(it.eq.k)dumx=(gam1(in1)-two)
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(j).or.oil)goto320
      in3=ind(iwork(l,i),iwork(j,it))
      h(in3)=h(in3)+dumx
320   if(isymmo(it).ne.isymmo(l).or.oij)goto330
      in3=ind(iwork(l,it),iwork(j,i))
      h(in3)=h(in3)+dumx*fourm
330   continue
      goto660
c
c
c......(at|ij) integrals
 340  oik=isymmo(i).ne.isymmo(k)
      oil=isymmo(i).ne.isymmo(l)
      do 360 it=nst,nprim
      if(isymmo(it).ne.isymmo(j))goto360
      in1=ic1e(max(it,j),1)+ic1e(min(it,j),2)
      dumx=gam1(in1)*ggg*half
      if(it.eq.j)dumx=(gam1(in1)-two)*ggg
      if(isymmo(it).ne.isymmo(k).or.oil)goto350
      in2=ind(iwork(l,i),iwork(k,it))
      h(in2)=h(in2)+dumx
350   if(isymmo(it).ne.isymmo(l).or.k.eq.l.or.oik)goto360
      in2=iky(iwork(k,i))+iwork(l,it)
      h(in2)=h(in2)+dumx
360   continue
      goto660
c
370   if(j.gt.ncore)goto 420
c
c.....(ai|tu) integrals
c
      ivx=icfcor(k)+l
      isymvx=ic4e(ivx)
      do 390 it=nst,nprim
      oiti=isymmo(it).ne.isymmo(i)
      oitj=isymmo(it).ne.isymmo(j)
      do 390 iu=nst,it
      itu=icfcor(it)+iu
      if(ic4e(itu).ne.isymvx)goto390
      oiui=isymmo(iu).ne.isymmo(i)
      oiuj=isymmo(iu).ne.isymmo(j)
      dumx=ggg*gam2(ic2e(max(itu,ivx))+ic3e(min(itu,ivx)))
      if(itu.ne.ivx)dumx=dumx*half
      if(oitj.or.oiui) go to 380
      in1=iky(iwork(iu,i))+iwork(j,it)
      h(in1)=h(in1)-dumx
380   if(oiti.or.oiuj) go to 390
      in1=iky(iwork(it,i))+iwork(j,iu)
      h(in1)=h(in1)-dumx
390   continue
c
      ojk=isymmo(j).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      do 410 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto410
      in3=iky(iwork(it,i))
      if(isymmo(it).ne.isymmo(l).or.ojk)goto400
      dumx=gam1(ic1e(max(it,l),1)+ic1e(min(it,l),2))*ggg
      if(it.ne.l)dumx=dumx*half
      in2=in3+iwork(j,k)
      h(in2)=h(in2)-dumx
400   if(isymmo(it).ne.isymmo(k).or.k.eq.l.or.ojl)goto410
      dumx=gam1(ic1e(max(it,k),1)+ic1e(min(it,k),2))*ggg
      if(it.ne.k)dumx=dumx*half
      in2=in3+iwork(j,l)
      h(in2)=h(in2)-dumx
410   continue
c
      goto660
c
c.....(au/vx) integrals
c
420   n2=icfcor(k)+l
c
      do 430 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))go to 430
      ll=iwork(it,i)
      n1=icfcor(max(it,j))+min(it,j)
      temp=ggg*gam2(ic2e(max(n1,n2))+ic3e(min(n1,n2)))
      if(j.ne.it)temp=temp*half
      if(n1.ne.n2)temp=temp*half
      brill(ll)=brill(ll)+temp
430    continue
      goto660
c
c-----------------------------------------------------
c.....(vi/tu)---> (tu/vi)
c
440   if(j.gt.ncore)goto610
       k1=k
      l1=l
      k=i
      l=j
      i=k1
      j=l1
c
c
c.....(tu/ij) integrals
450   n1=icfcor(i)+j
      do 460 it=nst,nprim
      if(isymmo(it).ne.isymmo(l))go to 460
      ll=iwork(l,it)
      n2=icfcor(max(it,k))+min(it,k)
      temp=ggg*gam2(ic2e(max(n1,n2))+ic3e(min(n1,n2)))
      if(k.ne.it)temp=temp*half
      if(n1.ne.n2)temp=temp*half
      brill(ll)=brill(ll)-temp
460   continue
      goto660
470   if(j.le.ncore)goto660
      ivx=icfcor(i)+j
      isymvx=ic4e(ivx)
      okl=k.eq.l
      do 490 it=nst,nprim
      oitk= isymmo(it).ne.isymmo(k)
      oitl= isymmo(it).ne.isymmo(l)
      do 490 iu=nst,it
      in1=icfcor(it)+iu
      if(ic4e(in1).ne.isymvx)goto490
      oiul=isymmo(l).ne.isymmo(iu)
      oiuk=isymmo(k).ne.isymmo(iu)
      dumx=gam2(ic2e(max(ivx,in1))+ic3e(min(ivx,in1)))*ggg
      if(it.ne.iu)dumx=dumx*half
      if(ivx.eq.in1)dumx=dumx+dumx
      if(oitk.or.oiul) go to 480
      in=iky(iwork(k,it))+iwork(l,iu)
      h(in)=h(in)+dumx
480   if(oiuk.or.it.eq.iu.or.okl.or.oitl)goto490
      in=iky(iwork(k,iu))+iwork(l,it)
      h(in)=h(in)+dumx
 490  continue
      oil=isymmo(i).ne.isymmo(l)
      oik=isymmo(i).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      ojk=isymmo(j).ne.isymmo(k)
      do 530 it=nst,nprim
      if(isymmo(it).ne.isymmo(i).or.i.eq.j)goto510
      in1=ic1e(max(it,i),1)+ic1e(min(it,i),2)
      dumx=gam1(in1)*half
      if(it.eq.i)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(k).or.ojl.or.(k.eq.l.and.j.gt.it))
     +   go to 500
      in=ind(iwork(k,it),iwork(l,j))
      h(in)=h(in)+dumx
500   if(isymmo(it).ne.isymmo(l).or.ojk.or.(k.eq.l.and.it.gt.j))
     +   go to 510
      in=iky(iwork(k,j))+iwork(l,it)
      h(in)=h(in)+dumx
510   if(isymmo(it).ne.isymmo(j))goto530
      in1=ic1e(max(it,j),1)+ic1e(min(it,j),2)
      dumx=gam1(in1)*half
      if(it.eq.j)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(k).or.oil.or.(k.eq.l.and.i.gt.it))
     +   go to 520
      in=ind(iwork(k,it),iwork(l,i))
      h(in)=h(in)+dumx
520   if(isymmo(it).ne.isymmo(l).or.oik.or.(k.eq.l.and.it.gt.i))
     +   go to 530
      in=iky(iwork(k,i))+iwork(l,it)
      h(in)=h(in)+dumx
530   continue
      goto660
c
c.....(ti|uj) integrals
c
540   ivx=icfcor(i)+k
      isymvx=ic4e(ivx)
      do 560 it=nst,nprim
      oiti=it.ne.i
      oitk=it.ne.k
      oitj=isymmo(it).ne.isymmo(j)
      oitl=isymmo(it).ne.isymmo(l)
      do 560 iu=nst,it
      if(mult(isymmo(it),isymmo(iu)).ne.isymvx)goto560
      oiuj=isymmo(iu).ne.isymmo(j)
      oiul=isymmo(iu).ne.isymmo(l)
      in1=icfcor(max(iu,i))+min(iu,i)
      in2=icfcor(max(it,k))+min(it,k)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      if(iu.eq.i)dumx=dumx+dumx
      if(oitk   )dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      in1=icfcor(max(iu,k))+min(iu,k)
      in2=icfcor(max(it,i))+min(it,i)
      dumy=ggg*gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))
      if(oiti   )dumy=dumy*half
      if(iu.eq.k)dumy=dumy+dumy
      if(in1.eq.in2)dumy=dumy+dumy
      if(oitj.or.oiul         )goto550
      in3=ind(iwork(j,it),iwork(l,iu))
      h(in3)=h(in3)+dumy
550   if((j.eq.l.and.i.eq.k).or.(it.eq.iu.and.j.ne.l))goto560
      if(oitl.or.oiuj        )goto560
      in3=ind(iwork(l,it),iwork(j,iu))
      h(in3)=h(in3)+dumx
560   continue
c
      okj=isymmo(j).ne.isymmo(k)
      okl=isymmo(k).ne.isymmo(l)
      oil=isymmo(i).ne.isymmo(l)
      oij=isymmo(i).ne.isymmo(j)
      do 600 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto580
      in1=ic1e(max(it,i),1)+ic1e(min(it,i),2)
      dumx=gam1(in1)*half
      if(it.eq.i)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(j.eq.l.and.it.eq.k)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(l).or.okj)goto570
      in3=ind(iwork(l,it),iwork(j,k))
      h(in3)=h(in3)+dumx
570   if(isymmo(it).ne.isymmo(j).or.okl)goto580
      in3=ind(iwork(j,it),iwork(l,k))
      h(in3)=h(in3)+dumx*fourm
580   if(isymmo(it).ne.isymmo(k).or.(i.eq.k.and.j.eq.l))goto600
      in1=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      dumx=gam1(in1)*half
      if(it.eq.k)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(j.eq.l.and.it.eq.i)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(j).or.oil)goto590
      in3=ind(iwork(j,it),iwork(l,i))
      h(in3)=h(in3)+dumx
590   if(isymmo(it).ne.isymmo(l).or.oij)goto600
      in3=ind(iwork(l,it),iwork(j,i))
      h(in3)=h(in3)+dumx*fourm
600   continue
      goto660
610   iuv=icfcor(i)+j
      ixy=icfcor(k)+l
      do 650 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto620
      in1=ilifp(it)+i
      in2=icfcor(max(it,j))+min(it,j)
      dumx=gam2(ic2e(max(ixy,in2))+ic3e(min(ixy,in2)))*ggg
      if(it.ne.j)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
620   if(isymmo(it).ne.isymmo(j).or.i.eq.j)goto630
      in1=ilifp(it)+j
      in2=icfcor(max(it,i))+min(it,i)
      dumx=gam2(ic2e(max(in2,ixy))+ic3e(min(in2,ixy)))*ggg
      if(it.ne.i)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
630   if(iuv.eq.ixy)goto650
      if(isymmo(it).ne.isymmo(k))goto640
      in1=ilifp(it)+k
      in2=icfcor(max(it,l))+min(it,l)
      dumx=ggg*gam2(ic2e(max(in2,iuv))+ic3e(min(in2,iuv)))
      if(it.ne.l)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
640   if(isymmo(it).ne.isymmo(l).or.k.eq.l)goto650
      in1=ilifp(it)+l
      in2=icfcor(max(it,k))+min(it,k)
      dumx=gam2(ic2e(max(in2,iuv))+ic3e(min(in2,iuv)))*ggg
      if(it.ne.k)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
650   continue
c
660   continue
      if(jj)50,680,50
680   continue
      call dscal(kk10,two,brill,1)
      return
10    format(' commence first derivative evaluation at ',f8.2,
     * ' seconds',a10,' wall')
20    format(' commence first and second derivative',
     1' evaluation at ',f8.2,' seconds',a10,' wall')
      end
_ENDEXTRACT
      subroutine scanci(gam11,gam22,skv,ccc,gam2,
     +   twoe,excint,ic5e,ic6e,h,ic7e,ic8e,ifclos,grad,nconf)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
INCLUDE(common/restar)
c
      integer    mxsurd
      parameter (mxsurd=3600)
      common/blkin /vsurd(mxsurd),nnsurd
      common/junk /surd(mxsurd),gam1(401)
c
      common/craypk/lin(2043),nnword
INCLUDE(common/iofile)
INCLUDE(common/qice)
INCLUDE(common/machin)
      common/potn/core,potnuc
INCLUDE(common/simul)
INCLUDE(common/filel)
INCLUDE(common/exc)
INCLUDE(common/prints)
INCLUDE(common/gjs)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,ns
     1t,n1e,isym,mults,macro,maxc,itypci
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,k7,kk8,kk9,kk10p
      common/cone  /m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n8,iperm(8,8),ims(8)
c
      character*10 charwall
      dimension gam11(*),gam22(*),ccc(nconf,*),iy(2)
      dimension skv(*),ic7e(*),ic8e(*),twoe(*),h(kk3),
     +          gam2(n2e),excint(*),ic5e(*),ic6e(*),
     +          grad(*)
c
      equivalence (iy(1),vsurd(1))
c
      data m1/1/
c
c     iunit1 = ntapes(1)
      iunit2 = ntapes(2)
      rewind iunit2
c
c     iu=kk6+1
      idomkm=idom(kkkm)
      kk3m=kk3-1
      n1st=ic1e(nst,1)
      n11e=ic1e(nprimp,1)-n1st
      n12e=n2e+n11e
c
c.....read in integrals
      call vclr(twoe,1,n12e)
      call setsto(1360,0,lin)
      lfile=m6file
      do 600 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i)  =m6blk(i)
 600  llblk(i)  =liblk(i)-m6last(i)
      do 1 i=1,lfile
      iunit=lotape(i)
      call search(liblk(i),iunit)
      call find(iunit)
      j=llblk(i)
 2    j=j+1
      call get(vsurd(1),mword)
      if(mword)3,1,3
 3    if(j.ne.0)call find(iunit)
      call rdxint(twoe,excint,ic7e)
      if(j)2,1,2
 1    continue
      kk10=kk10p-1
      ifclos=-1
_IF1()      oclose=.false.
_IF1()      open=.not.oclose
c
c....construct orbital gradient
c
      kk100=kk10+kk3-1
_IF1(civ)      call szero(grad,kk100)
_IFN1(civ)      call vclr(grad,1,kk100)
      call dcopy(kk9,gam11,1,gam1,1)
      call wrt3(gam22,kk8,m4blk(1),m4tape(1))
c
      call grget(
     *grad,gam1,gam2,gam22,excint,ic5e,ic6e,ic7e,ic8e,skv)
      do 10 kk=1,kk10
      ip=ic5e(kk+1)
      if(ip.le.ncore)grad(kk)=grad(kk)
     1        +skv(ic1e(ic6e(kk+1),1)+ic1e(ip,2))*4.0d0
10    grad(kk)=grad(kk)*0.5d0
c
      ve=ddot(kk9,skv,1,gam11,1)+ddot(kk8,twoe,1,gam22,1)
      if(noci(macro).eq.0)goto50
      vet=ve+core
      vl(kkkm)=ve
      if(nprint.ne.-5)write(iwr,40)vet,ve,core
40    format(/' ci stage omitted - total energy =',f16.7,' ',
     1'electronic energy =',f16.7,' core energy =',f16.7)
50    continue
c     iu=kk3+kk10p
      tst=cpulft(1)
      if(nprint.ne.-5)write(iwr,60)tst ,charwall()
60    format(/' start of ci hessian construction'
     1,' at ',f8.2,'  seconds',a10,' wall')
      kkk=kkkm
      ll=0
c
c..... formula tape
      iblkk=mtblk(1)
      ifillp=mttape(1)
      call rdedx(vsurd,mach(17),iblkk,ifillp)
      nsurd=nnsurd
      call dcopy(nsurd,vsurd,1,surd,1)
      call find(ifillp)
      call get(vsurd,nw)
_IF(ibm,vax)
      call upakin(iy,lin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,lin,2044)
_ELSE
      call upakft(iy,16,lin,2044)
_ENDIF
      iblkk=lensec(mach(17))   + iblkk
      iblkk=iblkk + (lin(1)*65536+lin(2)-1)/681+1
      call search(iblkk,ifillp)
      call find(ifillp)
      call get(vsurd,nw)
      call find(ifillp)
_IF(ibm,vax)
      call upakin(iy,lin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,lin,2044)
_ELSE
      call upakft(iy,16,lin,2044)
_ENDIF
      map=1
      m2043=nnword*3
c
c.....copy up 1-electron integrals
      call dcopy(n11e,skv(n1st+1),1,twoe(n2e+1),1)
c
260   ll=ll+1
      if(ll.gt.kk3)goto 490
c
c
c....omit row of dom. config..
      if(ll.ne.idomkm)go to 280
c.....pass over formulae for the omitted row
c
      nelo=lin(map+1)+1
      do 270 i=1,nelo
      map=map+3
      if(map.lt.m2043)goto 270
      call get(vsurd,nw)
      call find (ifillp)
_IF(ibm,vax)
      call upakin(iy,lin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,lin,2044)
_ELSE
      call upakft(iy,16,lin,2044)
_ENDIF
      map=1
      m2043=nnword*3
270   continue
      goto 260
c
280   lll=ll+kk10
      call vclr(h,1,kk100)
      if(ll.gt.idomkm)lll=lll-1
      vt=1.0d0
_IF1()      if(mta.ne.mtb.and.oclose)vt=2.0d0
c     ir=1
      ved=0.0d0
      call vclr(gam1,1,kk9)
      call vclr(gam22,1,n12e)
c
c
c.....formula tape
c
      nelo=lin(map+1)
      if(ll.ne.lin(map))call caserr('formula tape screwed')
      map=map+3
      if(map.lt.m2043)goto 360
      call get(vsurd,nw)
      call find(ifillp)
_IF(ibm,vax)
      call upakin(iy,iin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,lin,2044)
_ELSE
      call upakft(iy,16,lin,2044)
_ENDIF
      map=1
      m2043=nnword*3
360   do 410 i=1,nelo
      nn=lin(map)
      intlab=lin(map+1)
      val=surd(lin(map+2))
      map=map+3
      if(map.lt.m2043)goto 370
      call get(vsurd,nw)
      if(nw.eq.0)go to 370
      call find(ifillp)
_IF(ibm,vax)
      call upakin(iy,lin,511)
_ELSEIF(cray,ksr,i8)
      call unpack(iy,16,lin,2044)
_ELSE
      call upakft(iy,16,lin,2044)
_ENDIF
      map=1
      m2043=nnword*3
370   h(nn)=h(nn)+twoe(intlab)*val
410   gam22(intlab)=gam22(intlab)+val*ccc(nn,kkkm)
      if(oprint(8))write(iwr,371)(h(nn),nn=1,kk3)
371   format(5f16.8)
c
c....push this row of the ci matrix into its place in the hessian
      ved=h(idomkm)
c
      is=idomkm+1
      ip=kk3+is
      iq=ip+kk10-1
      do 411 i=is,kk3
411   h(iq-i)=h(ip-i)
      is=idomkm-1
      ip=idomkm+kk10
      do 412 i=1,is
412   h(ip-i)=h(idomkm-i)
c
c.....copy down 1-particle density matrix
      call dcopy(n11e,gam22(n2e+1),1,gam1(n1st+1),1)
c
      h(lll)=h(lll)-ve*vt
      grad(lll)=grad(lll)+ved*ccc(idomkm,kkkm)
      vc=ccc(ll,kkkm)
      call daxpy(kk3m,vc,h(kk10p),1,grad(kk10p),1)
c
c.....get the contributions to the mixed hessian
c
      vc=vc*vt
      vc2m=-vc*2.0d0
_IF1(civ)      call scaler(kk10,vc2m,h,grad)
_IFN1(civ)      call vsmul(grad,1,vc2m,h,1,kk10)
      call grget (h,gam1,gam2,gam22,excint,ic5e,ic6e,ic7e,ic8e,skv)
      call daxpy(kk10,vc2m,grad,1,h,1)
      vc=vc*4.0d0
      do 460 kk=1,kk10
      ip=ic5e(kk+1)
      if(ip.le.ncore)h(kk)=h(kk)+vc*
     +      skv(ic1e(ic6e(kk+1),1)+ic1e(ip,2))
460   continue
c
      if(oprint(8))write(iwr,371)(h(kk),kk=1,ll)
c
c.....dump this row of the hessian onto fortran stream 2
ccray buffer out(2,0)(h(1),h(lll))
      call wtfors(h,lll,iunit2)
c
      go to 260
c
490   iq=2
      call rdedx(gam22,kk8,m4blk(1),m4tape(1))
      call wrt3(grad,kk100,m4blk(1),m4tape(1))
      vm=cpulft(1)
      if(nprint.ne.-5)write(iwr,500)vm ,charwall()
500   format('   end of ci hessian construction at ',f8.2,
     *' seconds',a10,' wall')
      return
      end
      subroutine grget(brill,gam1,gam2,gam2p,excint,ic5e,ic6e,ic7e,
     +                 ic8e,skv)
c
c
c     form brillouin matrix elements, given density matrices -
c       the routine will work for full and derivative density
c     matrices.
c
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb),len3(8) ,ibext(8)
     >                           ,lenext(8)
      common/tab3  /is1,is1p,is2,is2p
      common/dims  /kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10p
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/qice)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     +            naa,nab,nst,n1e,isym,mults,macro
c
      dimension gam1(*),gam2(*),brill(*),excint(*),ic7e(*),ic8e(*),
     +          ic5e(*),ic6e(*),skv(*),gam2p(*)
c
c
c....expand density matrix
c.....one-particle .. get rid of core-active 'holes'
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      iuv=1
      do 10 iu=nst,nprim
      do 10 iv=nst,iu
      if(isymmo(iu).ne.isymmo(iv))goto 10
      gam1(iuv)=gam1(ic1e(iu,1)+ic1e(iv,2))
      iuv=iuv+1
10    continue
c
c.....two-particle .. can be looked up both at ic8e(tu)+ic3e(vx) or
c     ic8e(vx)+ic3e(tu) .. also non standard factors of 2 for speed
      lenact=iky(nact+1)
      do 44 j=1,nact
44    ic8e(j+nact)=ic8e(j)
      ijp=1
      do 45 i=1,nact
      map=ic8e(i)-1
      do 45 j=1,i
      ij=ijp
      ijp=ij+1
      isymij=ic4e(ij)
      ic2eij=ic2e(ij)
      ic3eij=ic3e(ij)
      do 20 kl=1,ic3eij
20    gam2(map+kl)=gam2p(ic2eij+kl)
      map=map+ic3eij
      gam2(map)=gam2(map)+gam2(map)
      do 30 kl=ijp,lenact
      if(ic4e(kl).ne.isymij)goto 30
      map=map+1
      gam2(map)=gam2p(ic2e(kl)+ic3eij)
30    continue
      lenij=lena(isymij)
      if(i.ne.j)goto 41
      mapij = map-lenij+1
      call dscal(lenij,2.0d0,gam2(mapij),1)
      goto 45
41    call dcopy(lenij,gam2(map-lenij+1),1,gam2(ic8e(j+nact)),1)
      ic8e(j+nact)=ic8e(j+nact)+lenij
45    ic8e(i+nact)=ic8e(i+nact)+lenij
c
      if(ncore.eq.0)goto 79
c
c.....i-->a terms
c
      lia=is2-is1
      la=lena(1)
      ilifpp=ilifp(ic6e(is1p))+ic5e(is1p)
      ilifpp=ic7e(ilifpp)+1
      call mxmb (excint(ilifpp),la,1,
     >      gam1,1,la, brill(is1),1,lia,  lia,la,1)
      ijunk=is1-1
      do 100 kk=1,ijunk
      ip=ic5e(kk+1)
      iq=ic6e(kk+1)
c
c.....two particle term
c
      ilifpp=ilifp(ip)
      iqc=iq-ncore
      isymq=isymmo(iq)
      ilifpq=ic7e(ilifpp+nst)+1
      ilifpr=ic8e(iqc)
      buff=0.0d0
      if(len3(isymq).ne.0)buff=
     >   -ddot(len3(isymq),excint(ilifpq),1,gam2(ilifpr),1)
c
c.....fock and 1 particle term
c
      iuv=1
      map=ic7e(ilifpp+iq)
      itv=iky(iqc)
      do 50 iv=nst,nprim
      itv=itv+1
      if(iv.gt.iq)itv=icfcor(iv)+iq
      ic3tv=ic3e(itv)
      if(isymmo(iv).ne.isymq)goto 49
      dumx=gam1(ic3tv)*skv(ic1e(iv,1)+ic1e(ip,2))
      if(iq.eq.iv)dumx=dumx+dumx
      buff=buff-dumx
49    ic7pv=ic7e(ilifpp+iv)
      do 50 iu=nst,iv
      if(isymmo(iu).ne.isymmo(iv)) goto 50
      iuc=iu-ncore
      buff=buff+gam1(iuv)*(4.0d0*excint(map+iuv)
     1   -  excint(ic7pv+ic3e(ind(iuc,iqc)))
     2   -  excint(ic7e(ilifpp+iu)+ic3tv))
      iuv=iuv+1
50    continue
100   brill(kk)=brill(kk)+buff
c
c....t-->a    no fock term
c.....two-electron bit
c
79    kk=is2
      do 900 ipc=1,nact
      isymp=isymmo(ipc+ncore)
      l3=len3(isymp)
      le=lenext(isymp)
      iaaaaa=ibext(isymp)
      ibbbbb=ic8e(ipc)
      call mxmb (excint(iaaaaa      ),l3,1,  gam2(ibbbbb   ),1,l3,
     >           brill(kk),1,le,    le,l3,1)
900   kk=kk+le
      ijunk=kk10p-1
      do 1000 kk=is2,ijunk
      ip=ic5e(kk+1)
      iq=ic6e(kk+1)
      iqp=ic1e(iq,1)
      isymp=isymmo(iq)
      buff=0.0d0
      itv=icfcor(ip)
      do 90 iv=nst,ip
      if(isymp.eq.isymmo(iv)) buff = buff +
     >                     gam1(ic3e(itv+iv))*skv(iqp+ic1e(iv,2))
90    continue
      do 95 iv=ip,nprim
      if(isymp.eq.isymmo(iv)) buff = buff +
     >                    gam1(ic3e(icfcor(iv)+ip))*skv(iqp+ic1e(iv,2))
95    continue
1000  brill(kk)=buff+brill(kk)
      return
      end
      subroutine augmnt (q,max,brill,n4,buffio,diag,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      common /blkin / vlamda(maxorb),work(maxorb)
INCLUDE(common/iofile)
INCLUDE(common/exc)
INCLUDE(common/simul)
      common /qpar  / j4(3),j2(11),macro
c
      character*10 charwall
      dimension       q(max),brill(*),buffio(*),diag(*),iguess(5)
c
      data  thresh,huge /1.0d-3,1.0d16/
      data  nits/100/,dzero,done/0.0d0,1.0d0/
      data  acclo,acchi /1.0d-9,1.0d-4/,m1/1/
c
      iunith = ntapes(1)
      nroot=1
      ifail = 0
      if (isimul(macro).ne.0) nroot=kkkm
      nstart=1
      nend=nroot
      n=n4
      np=n+1
      np4=np
_IF1(civ)      call szero(q,max)
_IFN1(civ)      call vclr(q,1,max)
      acc=1.0d-2*ddot(n,brill,1,brill,1)
      if(acc.gt.acchi) acc=acchi
      if(acc.lt.acclo) acc=acclo
c
c.... ==storage==
c
c     | expansion set | sigma set  | small matrix |
c     ic              is           ig
c
      buff= dfloat(np4)
      maxvec =  dsqrt(buff**2+max) - buff
      maxvec = min(maxvec,maxorb,n4)
      if (maxvec.lt.nend*2) call caserr
     + ('nothing like enough store for augmented hessian diag.')
c
      ic    = 1
      is    = ic + np4*maxvec
      ig    = is + np4*maxvec
c
      tst = cpulft(1)
      if(nprint.ne.-5)write(iwr,10) tst ,charwall()
10    format(/' commence augmented hessian diagonalisation at',
     1f8.2,' seconds',a10,' wall')
c
c.....read h & get diag elements
c
      rewind iunith
      diag(1)=dzero
      do 20 irow=1,n
      call rdfors (buffio,irow,iunith)
20    diag (irow+1) = buffio (irow)
      rewind iunith
c
c.....trial solutions
      vlast=-huge
      do 40 i=1,nroot
      accum=huge
      do 30 j=1,np
      val = diag(j)
      if (val.le.vlast.or.val.gt.accum) goto 30
      iguess(i)=j
      accum=val
30    continue
      q (ic-1+(i-1)*np4+iguess(i)) = done
40    vlast=accum
c
c.... special code for 1st iter when vectors are (0..00100..0)
      mtx = 1
      if (iguess(1).eq.1.and.nroot.eq.1) goto 90
      do 80 irow = 1,n
      isig = is
      buffio(1) = brill (irow)
      call rdfors (buffio(2),irow,iunith)
      do 70 ivec = 1,nroot
      k = iguess(ivec)-1
      if (irow-k) 70,50,60
50    call dcopy(k,buffio,1,q(isig),1)
60    q(isig+irow)=buffio(k+1)
70    isig = isig + np4
80    continue
      rewind iunith
      goto 100
90    continue
c
      call dcopy(n,brill,1,q(is+1),1)
100   do 110 i=1,nend
      im=i-1
      do 110 j=1,i
      jm=j-1
      q(ig+jm*maxvec+im)=ddot(np,q(ic+im*np4),1,q(is+jm*np4),1)
 110  continue
c
      call f02abf (q(ig),maxvec,nend,vlamda,q(ig),maxvec,work,
     * ifail)
c
      test=dzero
      nvec = nend - nstart + 1
      igg = ig + (nstart-1)
      do 120 i=1,nroot
      r = dnrm2(nvec,q(igg),1)
      if (r.lt.test) goto 120
      iworst = i
      test = r
120   igg = igg+maxvec
      if (test.lt.acc) goto 230
      if (mtx .ge. nits) call caserr
     + ('augmented matrix diagonalisation has failed to converge')
      if (nend+nroot.gt.maxvec) goto 190
c
      ibase = nend
      i4 = ic+ibase*np4-1
      do 140 k=1,nroot
      vlk=vlamda(k)
      if (mtx.eq.1) vlk=vlk+1.0d-8
      call vclr(q(i4+1),1,np)
      do 130 j=1,ibase
      gjk=q(ig+(k-1)*maxvec+j-1)
      i44 = (j-1)*np4-1
      do 130 i=1,np4
130   q(i4+i)=q(i4+i)+gjk*(q(is+i44+i)-vlk*q(ic+i44+i))/(vlk-diag(i))
      s=done/dnrm2(np,q(i4),1)
      call dscal(np,s,q(i4+1),1)
c
c....project off components of previous vecs & see if worth keeping
c
      nendp = nend+1
      call schmid (q(ic),np,nendp,m1,s,ncol)
c
      if (s.lt.thresh .and. k.ne.iworst) goto 140
      nend = nendp
      i4 = i4 + np4
140   continue
c
c.....re-do orthonormality
150   nstart=ibase+1
      nvec=nend-ibase
      call schmid (q(ic),np,nstart,nvec,s,ncol)
c
c..... get sigma vector for current new vectors
      mtx=mtx+1
      i4=(nstart-1)*np4
      icc=ic+i4
      isig=is+i4
      do 170 ivec=nstart,nend
      q(isig)=ddot(n,brill,1,q(icc+1),1)
c     do 160 irow=1,n
c160  q(isig+irow)=brill(irow)*q(icc)
_IFN1(civ)      call vsmul(brill,1,q(icc),q(isig+1),1,n)
_IF1(civ)      call scaler(n,q(icc),q(isig+1),brill)
      icc=icc+np4
170   isig=isig+np4
c
      do 180 irow=1,n
c     irow2=irow
      irowm=irow-1
      icc=ic+i4
      isig=is+i4
      call rdfors (buffio,irow,iunith)
      do 180 ivec=nstart,nend
      q(isig+irow)=0.0d0
      if(irowm.gt.0)q(isig+irow) =
     *ddot(irowm,buffio(1),1,q(icc+1),1)
      call daxpy(irow,q(icc+irow),buffio,1,q(isig+1),1)
      icc=icc+np4
180   isig=isig+np4
      rewind iunith
      goto 100
c
c
c.... code to collect up when you run out of core
c     calculate solutions so far with sigma vectors, dump to disc, then
c     reload
190   igg=ig-1
      i4 = maxvec*np4
      call vclr(q(is),1,i4)
      call mxmb (q(ic),1,np, q(ig),1,maxvec, q(is),1,np, np,nend,nroot)
      call vclr(q(ic),1,i4)
      do 220 kkkk=1,nroot
      call dcopy(np,q(is+(kkkk-1)*np4),1,q(ic+(kkkk-1)*np4),1)
 220  continue
      call vclr(q(is),1,i4)
      ibase = 0
      nend=nroot
      goto 150
c
c
c.....convergence reached; analyse solution
230   t = cpulft(1) - tst
      if(nprint.ne.-5)write(iwr,240) test,mtx,t
240    format(/' diagonalisation converged to accuracy',e8.1,' in',
     1       i4,' iterations taking',
     2    f8.2,' seconds')
      if(nprint.ne.-5)write(iwr,250) (vlamda(i),i=1,nroot)
250   format(/' lowest eigenvalues of augmented hessian:',5e12.4)
      if(nprint.ne.-5)write(iwr,260)
260   format(/' the next few eigenvalues, and the accuracy to which '
     1  ,' they are known are:')
      jm = nroot
      if (nroot.eq.1) jm = 2
      if (nend.lt.nroot*2) jm = nend-nroot
      do 280 i=1,jm
      r = dnrm2(nvec,q(igg),1)
      if (i.gt.2.and.r.gt.acchi) goto 290
      if(nprint.ne.-5)write(iwr,270) vlamda(i+nroot),r
270   format(2e12.4)
280   igg = igg + maxvec
290   continue
c
      call vclr(buffio,1,np)
      call mxmb (q(ic),1,np, q(ig+(nroot-1)*maxvec),1,maxvec,
     <      buffio,1,np,  np,nend,1)
c     do 310 i=1,n
c     brill(i)=buffio(i+1)/(-buffio(1))
      test=-1.0d0/buffio(1)
_IFN1(civ)      call vsmul(buffio(2),1,test,brill,1,n)
_IF1(civ)      call scaler(n,test,brill,buffio(2))
      if(nprint.eq.-5)return
      top=cpulft(1)
      write(iwr,320)top ,charwall()
 320  format(/1x,
     *'augmented hessian diagonalization complete at ',f8.2,' seconds'
     *,a10,' wall')
      return
c
      end
      subroutine modhes(q,maxq,brill,n4,buffio,diag,norb,nprint,
     +                  ireturn)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
      common /blkin / vlamda(maxorb),work(maxorb)
INCLUDE(common/exc)
INCLUDE(common/simul)
INCLUDE(common/ciconv)
      common /qpar  / j4(3),j2(11),macro
      dimension q(maxq),brill(*),buffio(*),diag(*),iguess(5)
      character*10 charwall
      data  thresh,huge /1.0d-3,1.0d6/,nits/100/,
     +      dzero,done /0.0d0,1.0d0/ ,acclo,acchi/1.0d-6,1.0d-3/,
     +      m1/1/,vk10,vk20,vk3 /.426d0, .113d0, .2d0 /
      nroot=1
      iunith = ntapes(1)
      ifail=0
      if (isimul(macro).ne.0) nroot=kkkm
      nstate = nroot
      nstart=1
      nend=nroot+1
      n=n4
      call vclr(q,1,maxq)
      acc=dnrm2(n,brill,1)
      if(acc.gt.acchi) acc=acchi
      if(acc.lt.acclo) acc=acclo
      buff= dfloat(n4)
      maxvec =  dsqrt(buff**2+maxq) - buff
      maxvec = min(maxvec,maxorb,n4)
      if(maxvec.lt.nend*2) call caserr 
     + ('nothing like enough store available for hessian diag.')
      ic    = 1
      is    = ic + n4*maxvec
      ig    = is + n4*maxvec
      tst = cpulft(1)
      if(nprint.ne.-5)write(iwr,10) tst ,charwall()
10    format(/' commence hessian diagonalisation at ',
     1f8.2,' seconds',a10,' wall')
      rewind iunith
      do 20 irow=1,n
      call rdfors (buffio,irow,iunith)
20    diag (irow) = buffio (irow)
      rewind iunith
      vlast=-huge
      do 40 i=1,nroot
      accum=huge
      do 30 j=1,n
      val = diag(j)
      if (val.le.vlast.or.val.gt.accum) goto 30
      iguess(i)=j
      accum=val
30    continue
      q (ic-1+(i-1)*n4+iguess(i)) = done
40    vlast=accum
      i4 = ic + nroot*n4
      call dcopy(n,brill,1,q(i4),1)
      do 50 i=1,nroot
50    q(i4+iguess(i)-1) = dzero
      s = done/dnrm2(n,q(i4),1)
      call dscal(n,s,q(i4),1)
      mtx = 1
      do 90 irow = 1,n
      irowm = irow - 1
      isig = is
      call rdfors (buffio,irow,iunith)
      do 80 ivec = 1,nroot
      k = iguess(ivec)
      if (irow-k) 80,60,70
60    call dcopy(k,buffio,1,q(isig),1)
70    q(isig+irow-1)=buffio(k)
80    isig = isig + n4
      q(isig+irowm) = ddot(irowm,buffio,1,q(i4),1)
      call daxpy(irow,q(i4+irowm),buffio,1,q(isig),1)
90    continue
      rewind iunith
c     goto 100
100   do 110 i=1,nend
      im=i-1
      do 110 j=1,i
      jm=j-1
      q(ig+jm*maxvec+im)=ddot(n,q(ic+im*n4),1,q(is+jm*n4),1)
110   continue
c     call square(q(ic),q(ig),maxvec,nend)
c     call eigrs(maxvec,nend,q(ic),vlamda,work,q(ig),ifail)
      call f02abf (q(ig),maxvec,nend,vlamda,q(ig),maxvec,work,ifail)
      test=dzero
      nvec = nend - nstart + 1
      igg = ig + (nstart-1)
      do 120 i=1,nroot
      r = dnrm2(nvec,q(igg),1)
      if (r.lt.test) goto 120
      iworst = i
      test = r
120   igg = igg + maxvec
      if (test.lt.acc) goto 210
      if (mtx .ge. nits) then 
        if (osuped) then
         ireturn = 1
         fmax = 100.0d0
         swnr = swnr/2.0d0
         swsimu = dmin1(swsimu,swnr)
         write(iwr,121) swnr
121      format(' ** hessian diagonalisation has failed to converge'/
     +          ' switch to newton-raphson reset to',f14.8)
         return
        else
         call caserr(
     +   'hessian diagonalisation has failed to converge')
        end if
      end if
      if (nend+nroot.gt.maxvec) goto 170
      ibase = nend
      i4 = ic+ibase*n4-1
      do 140 k=1,nroot
      vlk=vlamda(k)
      if (mtx.eq.1) vlk=vlk+1.0d-8
      call vclr(q(i4+1),1,n)
      do 130 j=1,ibase
      gjk=q(ig+(k-1)*maxvec+j-1)
      i44 = (j-1)*n4-1
      do 130 i=1,n4
130   q(i4+i)=q(i4+i)+gjk*(q(is+i44+i)-vlk*q(ic+i44+i))/(vlk-diag(i))
      s=done/dnrm2(n,q(i4+1),1)
      call dscal(n,s,q(i4+1),1)
      nendp = nend+1
      call schmid (q(ic),n,nendp,m1,s,ncol)
      if (s.lt.thresh .and. k.ne.iworst) goto 140
      nend = nendp
      i4 = i4 + n4
140   continue
      assign 150 to iexit
      if (nend.gt.maxvec) goto 170
150   nstart=ibase+1
      nvec=nend-ibase
      call schmid (q(ic),n,nstart,nvec,s,ncol)
      mtx=mtx+1
      i4=(nstart-1)*n4
      do 160 irow=1,n
      irowm=irow-1
      icc=ic+i4
      isig=is+i4
      call rdfors (buffio,irow,iunith)
      do 160 ivec=nstart,nend
      q(isig+irowm)=ddot(irowm,buffio,1,q(icc),1)
      call daxpy(irow,q(icc+irowm),buffio,1,q(isig),1)
      icc=icc+n4
160   isig=isig+n4
      rewind iunith
      goto 100
 170  continue
      call vclr(q(is),1,n*nroot)
      call mxmb (q(ic),1,n, q(ig),1,maxvec, q(is),1,n, n,nend,nroot)
      i4 = maxvec*n4
      call vclr(q(ic),1,i4)
      do 200 kkkk=1,nroot
      call dcopy(n,q(is+(kkkk-1)*n4),1,q(ic+(kkkk-1)*n4),1)
 200  continue
      call vclr(q(is),1,i4)
      ibase = 0
      nend=nroot
      goto iexit,(150,220)
210   assign 220 to iexit
      nroot = nroot+1
      goto 170
220   nroot = nroot -1
      vll = ddot(n,brill,1,q(ic+(nroot-1)*n4),1)
      if ( dabs(vll).le.vk20 .and. vlamda(nroot).gt.dzero) goto 230
      if (nroot .ge. maxvec/2) goto 230
      nroot = nroot+1
      goto 150
230   t = cpulft(1) - tst
      if(nprint.ne.-5)write(iwr,240) mtx,t
240    format(/' diagonalisation converged in',
     1       i4,' iterations taking',
     2    f8.2,' seconds')
      i4 = nroot - nstate
      i4 = max(i4,1)
      vk1  = vk10 / dsqrt( dfloat(i4))
      vk2  = vk20 / dsqrt( dfloat(i4))
      jm = nroot*2
      if (nroot .eq. 1) jm = jm+1
      if (jm .gt. nend) jm = nend
      i4 = ic +norb
      i44 = ic
      igg = ig + (nstart-1)
      nci = n-norb
      if(nprint.ne.-5)write(iwr,250)
250   format(/5x,'mode',8x,'eigenvalue',6x,'accuracy',6x,
     1'config char',7x, 'nr step',6x,'trunc step'/1x,87('*')/)
      do 300 imode=1,jm
      tau = dnrm2(nci,q(i4),1)
      r = dnrm2(nvec,q(igg),1)
      if (r.gt.acchi .and. imode.gt.nroot+2) goto 310
      if(nprint.ne.-5)write(iwr,260) imode,vlamda(imode),r,tau
260   format(i8,e18.4,f14.7,f15.5)
      if (imode.gt.nroot) goto 290
      vll = ddot(n,brill,1,q(i44),1)
      vvll =  dabs(vll)
      vln = vll
      stepsz = vk1
      if (nci.ne.0) stepsz = dmin1(stepsz,vk2/ dsqrt(tau))
      if ( dabs(vlamda(imode)).lt.2d-3) stepsz = dmin1(stepsz,vk3)
      if (vvll .gt. stepsz) vln = stepsz * (vll/vvll)
      if (vlamda(imode).le.dzero .and. imode.ge.nstate) vln = - vln
      if(nprint.ne.-5)write(iwr,270) vll
270   format('+',e73.4)
      if (vln .eq. vll) goto 290
      if(nprint.ne.-5)write(iwr,280) vln
280   format('+',e88.4)
      vln = vln - vll
      call daxpy(n,vln,q(i44),1,brill,1)
290   igg = igg + maxvec
      i4 = i4 + n4
      i44 = i44 + n4
300   continue
      if(nprint.ne.-5)write(iwr,400)
400   format(1x,87('-')//)
310   continue
      return
      end
      subroutine expmat(u,ndim,n,w1,vec,nterm)
c
c      form an approximation to the exponential of a skew-symmetric
c     matrix, by expanding to order nterm, and then applying a lowdin
c     orthogonalisation
c
      implicit REAL  (a-h,p-x),integer (i-n),logical  (o)
      implicit character *8 (z)
      implicit character *4 (y)
      dimension u(ndim,n),w1(n,n),vec(n)
      data done,halfm,half3,big,ten/1.0d0,-0.5d0,1.5d0,1.0d5,1.0d1/
c
      vn= dfloat(n)
_IF(c90)
      small=vn*x02ajf(vn)*ten
_ELSE
      small=vn*x02aaf(vn)*ten
_ENDIF
      n2=n
      do 10 j=1,n
      do 10 i=1,n
10    w1(i,j)=u(i,j)
c.....expand exponential
      if(nterm.lt.2)goto 40
      do 30 iterm=2,nterm
      scaler=done/ dfloat(nterm+2-iterm)
      do 30 i=1,n
_IF1(civ)      do 20 j=1,n
_IF1(civ)   20 vec(j)=u(i,j)*scaler
_IFN1(civ)      call vsmul(u(i,1),ndim,scaler,vec,1,n)
      vec(i)=vec(i)+done
      do 30 j=1,n
      u(i,j)=ddot(n2,vec,1,w1(1,j),1)
 30   continue
40    do 50 i=1,n
50    u(i,i)=u(i,i)+done
c
c.....now apply iterative symmetric orthogonalisation
      ncyc=0
60    c=0.0d0
      do 70 l=1,n
      do 70 j=1,n
      buff=ddot(n2,u(1,l),1,u(1,j),1)
      c=c+ dabs(buff)
70    w1(l,j)=buff
c
      c= dabs(c-vn)
      ncyc=ncyc+1
      if(c.lt.small.or.(ncyc.gt.30.and.c.lt.done))return
      if(c.gt.big.or.ncyc.gt.30)
     1callcaserr('lowdin orthogonalisation won''t converge')
c
      do 110 i=1,n
      do 90 j=1,n
_IF1(iv)      buff=0.0d0
_IF1(iv)      do 1 l=1,n
_IF1(iv)    1 buff=buff+u(i,l)*w1(l,j)
_IF1(iv)      vec(j)=buff
_IFN1(iv)      vec(j)=ddot(n,u(i,1),ndim,w1(1,j),1)
 90   continue
      do 100 j=1,n
100   u(i,j)=half3*u(i,j)+halfm*vec(j)
110   continue
      goto60
c
      end
      function mapsym(i,iset,nblock)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension nblock(8,*)
c
      do 10 j=1,maxorb
      if(i.eq.nblock(iset,j))goto 20
10    continue
      callcaserr(
     + 'failure in construction of newton-raphson solution matrix')
20    mapsym=j
c
      return
      end
      subroutine pople (q,max,brill,n4,buffio,diag,nprint)
c
c     pople's linear equation solver
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common /blkin / rhs(128),soln(128),wks1(128),wks2(128)
INCLUDE(common/iofile)
      dimension       q(max),brill(*),buffio(*),diag(*)
      character*10 charwall
_IF1()      dimension iguess(5)
      data  nits/500/ ,acclo,acchi/1.0d-9,1.0d-4/
      data  dzero,done,two,four /0.0d0,1.0d0,2.0d0,4.0d0/
      data  m1/1/
c
      iunith = ntapes(1)
      mtx = 0
      ifail = 0
      ivec = 1
      call vclr(rhs,1,128)
      n=n4
      call vclr(q,1,max)
c
      acc=ddot(n,brill,1,brill,1)
      if(acc.gt.acchi) acc=acchi
      if(acc.lt.acclo) acc=acclo
c
c.... ==storage==
c
c     | expansion set | small matrix | work space |
c     ic              ig             iw
c
      buff= dfloat(n4)/four
      maxvec =  dsqrt(buff**2+( dfloat(max)/two)) - buff
      i4 = n + 1
      maxvec = min(maxvec,128,i4)
      if(maxvec.lt.3) call caserr
     + ('nothing like enough store for pople equation solver')
c
      ic    = 1
      ig    = ic + n4*maxvec
      iw    = ig + maxvec*maxvec
c
      tst = cpulft(1)
      if(nprint.ne.-5)write (iwr,10) tst ,charwall()
10    format(/1x,'start of pople equation solver at ',f8.2,
     *' seconds',a10,' wall'/)
c
c.....read h & get diag elements
c
      rewind iunith
      do 20 irow=1,n
      call rdfors (buffio,irow,iunith)
20    diag (irow) = done / buffio (irow)
      rewind iunith
c
c.....trial solution
_IF1(ivc)      do 3  i=1,n
_IF1(ivc)    3 q(ic-1+i) = brill(i) * diag(i)
_IFN1(civ)      call vmul(brill,1,diag,1,q(ic),1,n)
      call schmid (q(ic),n,m1,m1,rhs(1),ncol)
      icc = ic-1
c
40    mtx=mtx+1
      ico = icc
      icc = icc + n4
      call rdfors (buffio,1,iunith)
      do 50 irow=2,n
      call rdfors (buffio,irow,iunith)
      buff = dzero
      irowm = irow - 1
      do 1234 j=1,irowm
      q(icc+j) = q(icc+j) + buffio(j)*q(ico+irow)*diag(j)
1234  continue
      q(icc+irow) = ddot(irowm,buffio,1,q(ico+1),1) * diag(irow)
50    continue
      rewind iunith
c
      igg = ig + (ivec-1)*maxvec - 1
      do 60 i=1,ivec
      q(igg+i) = ddot(n,q(ic+(i-1)*n4),1,q(icc+1),1)
60    continue
      q(igg+ivec) = q(igg+ivec) + done
c
      i4 = ivec
      call f04atf(q(ig),maxvec,rhs,i4,soln,q(iw),maxvec,wks1,wks2,
     * ifail)
c
      if ( dabs(soln(ivec)).lt.acc) goto 90
      if (mtx .ge. nits) call caserr
     + ('pople method has failed to converge')
      ivec = ivec + 1
      if (ivec.ge.maxvec) goto 70
      call schmid (q(ic),n,ivec,m1,s,ncol)
      q(ig+maxvec*(ivec-2)+ivec-1) = s
      goto 40
 70   continue
c
c.....convergence
      call vclr(buffio,1,n)
      ivec = ivec - 1
      do 80 k=1,ivec
      call daxpy(n,soln(k),q(ic+(k-1)*n4),1,buffio,1)
 80   continue
      call schmid (buffio,n,m1,m1,s,ncol)
      ivec = 1
      i4 = maxvec*n4
      call vclr(q(ic),1,i4)
      call dcopy(n,buffio,1,q(ic),1)
      goto 40
c
c.....convergence
 90   continue
      call vclr(brill,1,n)
      do 100 k=1,ivec
      call daxpy(n,soln(k),q(ic+(k-1)*n4),1,brill,1)
 100  continue
      t = cpulft(1)     - tst
      if(nprint.ne.-5)write (iwr,110) mtx,t ,charwall()
110    format(/' pople method converged in',i3,' iterations taking',
     1    f8.2,' seconds',a10,' wall')
c
      return
      end
      subroutine rdxint(twoe,excint,ic7e)
c
c     routine to read through mo integral file and store all integrals
c     needed for ci hessian construction. dummy ep for parameter list.
c     designed to be called through scan.
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension excint(*),twoe(*),ic7e(*)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/craypk/i205(1360)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
      common/potn  /core,potnuc
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1               nst,n1e,isym,mults
      common/blkin/g(510),mword
INCLUDE(common/atmblk)
INCLUDE(common/qice)
INCLUDE(common/gjs)
c
_IF1(iv)      call upak8v(g(num2e+1),i205)
_IFN1(iv)      int4=1
_IFN1(iv)      call unpack(g(num2e+1),lab816,i205,numlab)
      do 1  iword=1,mword
_IF(ibm,vax)
      i=i205(iword)
      j=j205(iword)
_ELSEIF(littleendian)
      j=i205(int4  )
      i=i205(int4+1)
_ELSE
      i=i205(int4  )
      j=i205(int4+1)
_ENDIF
      if(j.gt.nprim)goto 1
_IF(ibm,vax)
      k=k205(iword)
_ELSEIF(littleendian)
      k=i205(int4+3)
_ELSE
      k=i205(int4+2)
_ENDIF
      if(k.gt.nprim)goto 1
      if(i.le.ncore)goto 1
_IF(ibm,vax)
      l=l205(iword)
_ELSEIF(littleendian)
      l=i205(int4+2)
_ELSE
      l=i205(int4+3)
_ENDIF
      if(l.le.ncore.and.j.le.ncore)goto1
      if(k.le.ncore)goto 1
      if(i.gt.nprim)goto 30
      if(j.le.ncore)goto 10
      if(l.le.ncore)goto 20
c....(tu|vx)
      twoe(ic2e(icfcor(i)+j)+ic3e(icfcor(k)+l))=g(iword)
      goto1
c....(ti|uv)
10    excint(ic7e(ilifp(j)+i)+ic3e(icfcor(k)+l))=g(iword)
      goto 1
c.....(tu|vi)
20    excint(ic7e(ilifp(l)+k)+ic3e(icfcor(i)+j))=g(iword)
      goto 1
30    if(j.le.ncore)goto 40
      if(l.le.ncore)goto 50
c.....(au|vx)
      excint(ic7e(ilifp(i)+j)+ic3e(icfcor(k)+l))=g(iword)
      goto 1
c.....(ai|tu)
40    if(isymmo(i).ne.isymmo(j))goto 1
      map=ic7e(ilifp(i)+j)+ic3e(icfcor(k)+l)
      excint(map)=excint(map)+g(iword)*4.0d0
      goto 1
50    if(isymmo(i).ne.isymmo(l))goto 1
      if(k.lt.j)goto 60
c.....(at|ui)
      map=ic7e(ilifp(i)+l)+ic3e(icfcor(k)+j)
      excint(map)=excint(map)-g(iword)
60    if(k.gt.j)goto 1
      map=ic7e(ilifp(i)+l)+ic3e(icfcor(j)+k)
c.....(au|ti)
      excint(map)=excint(map)-g(iword)
_IF1(iv)    1 continue
_IFN1(iv)    1 int4=int4+4
c
      return
      end
c ******************************************************
c ******************************************************
c             =   casdm12    =
c ******************************************************
c ******************************************************
      subroutine anal(q,ccc,skv,qq,gam1,gam2,nconf,nprint)
c
c    analysis of casscf wavefunction
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/fsymas/dddum,otsym
      common/csymas/popp(maxorb)
INCLUDE(common/popnos)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/prints)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/machin)
      common/potn  /core,potnuc,corejj(509)
      common/blkin /h1(1)
INCLUDE(common/simul)
INCLUDE(common/ctrl)
INCLUDE(common/mapper)
INCLUDE(common/dm)
INCLUDE(common/cndx40)
INCLUDE(common/exc)
      common/atmol3/mina(3),moutb
      common/add   /iad1(5),iad2
INCLUDE(common/finish)
INCLUDE(common/gjs)
INCLUDE(common/harmon)
INCLUDE(common/infoa)
      common/tran/ilifc(maxorb),ntranc(maxorb),itran(mxorb3),
     * ctran(mxorb3),otran
      common/scfopt/maxit(4),accdi(2),icoupl(4),dmpcut,
     *            accc,enuc,etotal,ehf,ehf0(2),iterc
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,nst
     +            ,n1e,isym,mults,macro,maxc
      common/dims/kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8,kk9,kk10
      common/junkc/zjob(29)
      common/junk/eig(maxorb),vocc(maxorb),
     +            energ,nbas,newb,nco,ivalue,iocc,ispare(3)
      dimension skv(*),q(*),energy(3)
      dimension qq(*),gam1(*),gam2(*)
      dimension ccc(nconf,*),ibuf(6),bufx(6)
c     move to dynamic memory?
      character *15 ccctxt
      character*10 charwall
      dimension ccctxt(2000)
c
      data m3,m29/3,29/
      data dzero,half,two/0.0d0,0.50d0,2.0d0/
c
      data cipunch /1.0d-10/
c
      nav = lenwrd()
      lensq=length*length
      iad22=iad2+n1e
      iad3=iad22+lensq
      iad4=iad3+lensq
      iad5=iad4+length
      iad6=iad5+length
      iad7=iad6+length
c
      if(nprint.ne.-5)write(iwr,10)macro
10    format(/10x,'***** mcscf iterations stopped after'
     1,' cycle',i4,' *****')
      if(oconv.and.nprint.ne.-5)write(iwr,20)
 20   format(/40x,29('=')/ 40x,
     *'convergence threshold reached' /40x,29('=')  )
c
CMR
CMR   write(*,*) 'kkkm= ',kkkm,' vl= ',vl(kkkm),' core =',core
CMR   write(*,*) 'potnuc= ',potnuc
      ve=vl(kkkm)+core
      if(nprint.ne.-5)write(iwr,30)ve
      iterc=macro
30    format(/10x, '***********************************'/
     1        10x, 'final energy =',f20.11/
     2        10x, '***********************************'/)
c
      if(oprint(6))calldmprin(q(iad22),gam2,iwr,nprint)
c
      callsecget(isecv4,m3,iblk)
      call rdchr(zjob,m29,iblk,ndump4)
      call reads(eig,mach(8),ndump4)
      len2=1+ lensec(mach(8))     +  lensec(mach(9))
      iblkv=iblk+len2
      call rdedx(q,nblkq4,iblkv,ndump4)
      calldippol(q(iad22),q(iad22+lenb4),q,gam1,nprint)
      ivalue=-1
      iocc=1
      call vclr(vocc,1,nbas4)
      do 70 i=1,nprim
      vocc(i)=two
      if(i.gt.ncore)vocc(i)=pop(i)
70    continue
      energy(1)=vl(kkkm)+core
      energy(2)=vl(kkkm)
      energy(3)=core
      enuc=potnuc
      etotal = energy(1)
      ehf=etotal-enuc
      energ=energy(1)
c
c...   harmonic
c
      if (newbas1.ne.newbas0) call expharm(q,'ctrans',flop)
c
      call wrtc(zjob,m29,iblk,ndump4)
      call wrt3s(eig,mach(8),ndump4)
      call wrt3is(ilifc,mach(9)*nav,ndump4)
c
c ..... now canonicalise over f(a) + f(i)
c ..... load up to section moutb
c
c     if(icanon.eq.0.and.mode.ne.0)goto120
c
      call secget(moutb,m3,iblk)
      iblkv=iblk+len2
      if(icanon.eq.0)goto100
c.....canonicalise active and/or inactive spaces over specified
c      hamiltonian
      call vclr(eig,1,nbas4)
      call vclr(h1,1,121)
      call vclr(q(iad2),1,n1e)
      call dpconv(nprim,gam1,gam2)
      call vadd(skv,1,skv(n1e+1),1,q(iad2),1,n1e)
      if(icanon.eq.1) then
        call dcopy(kk9,gam1,1,h1,1)
      endif
      if(icanon.eq.3) then
        call dcopy(kk9,q(iad2),1,h1,1)
      endif
      if(icanon.ne.2)goto80
c.........icaono.eq.2 ==> primary orbs over mcscf fock operator
      ic=1
      do 60 i=1,nprim
      do 60 j=1,i
      if(isymmo(i).ne.isymmo(j))goto60
      if(j.eq.i)h1(ic)=q(iad2-1+ic)
      if(j.le.ncore)goto50
      h1(ic)=half*(qq((i-1)*nprim+j)+qq((j-1)*nprim+i))
50    ic=ic+1
60    continue
      mark=1
      if(oprint(8))call fout(h1,gam1,mark,iwr)
80    call gconv(nprim,gam1,gam2)
      call conon(q(iad22),q,q(iad3),h1,q(iad2),
     1eig,q(iad4),q(iad5),q(iad6),iwr,nprint)
      ivalue=1
      if(icanon.ne.1)goto120
c
      do 90 i=nst,nprim
      pop(i)=eig(i)
90    eig(i)=dzero
c
100   iocc=1
      call vclr(vocc,1,nbas4)
      do 110 i=1,nprim
      vocc(i)=two
      if(i.gt.ncore)vocc(i)=pop(i)
110   continue
c
120   continue
c
c...  recover from harmonic
c
      l0 = newbas0
      if (newbas0.ne.newbas1) then
c
c...    expand vectors 
c
         call expharm(q,'vectors',ilifq)
         ncol4 = num
         nsa4 = num
         newb4 = num
         nbas4 = num
         ma = num
      end if
c
      call wrt3(q,nblkq4,iblkv,ndump4)
c
      cpu=cpulft(1)
      write(iwr,9368)iterc,cpu,charwall(),ehf,enuc,etotal
 9368 format(/20x,14('=')/
     *20x,'final energies  after ',i4,' cycles at ',f8.2,' seconds'
     *,a10,' wall'/
     *20x,14('=')//
     *20x,'electronic energy ',f16.10/
     *20x,'nuclear energy    ',f16.10/
     *20x,'total energy      ',f16.10/)
      call wrtc(zjob,m29,iblk,ndump4)
      call wrt3s(eig,mach(8),ndump4)
      call wrt3is(ilifc,mach(9)*nav,ndump4)
      if(nprint.eq.-5)go to 140
      l=0
c     print ci coefficients
      write(iwr,201)
      do 200 i=1,kk3
      if( dabs(ccc(i,kkkm)).lt.ciprnt)go to 200
      l=l+1
      ibuf(l)=i
      bufx(l)=ccc(i,kkkm)
      if(l.lt.6)go to 200
      write(iwr,202)(ibuf(l),bufx(l),l=1,6)
      l=0
 200  continue
      if(l.ne.0)write(iwr,202)
     *(ibuf(i),bufx(i),i=1,l)
 201  format(/ 1x,104('*')/ /40x,'ci coefficients'/
     *40x,15('=')/)
 202  format( 1x,6(i5,f11.6,5x))
c
c     punch ci coefficients and configurations
c     retrieve coefficient list and then add coefficients
      if(opunch(16))then
       open(ipu,file='civecs.ascii',form='formatted',status='unknown')
       rewind ipu
       do i=1,kk3
        read(ipu,203)itmp,ccctxt(i)
203     format(1x,i6,3x,a15)
       enddo
       rewind ipu
c     now punch with ci coefficients
       do i=1,kk3
        if(dabs(ccc(i,kkkm)).ge.cipunch) then
         write(ipu,204) i,ccctxt(i),ccc(i,kkkm)
204      format(1x,i6,3x,a15,5x,f20.10)
        endif
       enddo
       close(ipu,status='keep')
      endif
c
      if(oprint(11).and.(.not.otsym))goto140
      call tdown(q,ilifq,q,ilifq,l0)
      if(otsym)call symass(q,eig,popp,q(iad7))
      if(oprint(11))goto140
      write(iwr,130)
130   format(/1x,104('=')//50x,'casscf mos (canonicalised)'/
     * 50x,26('='))
      call prsql(q,l0,nbas4,nbas4)
140   continue
c
c...   expand the vectors on ilbkq4 as well (for density)
c
      if (newbas0.ne.newbas1) then
         call rdedx(q,nblkq4,iblkq4,ndump4)
         call expharm(q,'vectors',ilifq)
         call wrt3(q,nblkq4,iblkq4,ndump4)
      end if
c
      call p1out(q(iad22),gam1,gam2,iwr,nprint)
      if(igrad.eq.0) call p1den(q,q(nblkq4+1),gam1,nprint)
      if(oconv) then 
       if(ifout2.ne.0)call p2out(gam1,gam2,iwr)
      endif
      call revind
c
      return
      end
      subroutine conon(f,q,work,skv1,skv,eig,eig1,eig2,iyc1,
     +                 iwr,nprint)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/qice)
      common/block /length,nsymm(8),nblock(8,maxorb)
INCLUDE(common/finish)
INCLUDE(common/gjs)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     +            naa,nab,nst,n1e,isym,mults,macro
      dimension skv(*),eig1(*),eig2(*),f(length,*),q(ma,*),skv1(*),
     +           work(length,*),iyc1(*),eig(*)
c
c.....canonicalise secondary orbitals over 1-particle hamiltonian
c     supplied (packed) in skv, & primary orbitals over skv1
c
      if(nprint.ne.-5)write(iwr,10)icanon
10    format(/' canonicalisation of virtual orbitals over fock'
     1,' matrix ',i2,' starting')
      if(icanon.eq.1.and.nprint.ne.-5)write(iwr,230)
      if(icanon.eq.2.and.nprint.ne.-5)write(iwr,240)
      if(icanon.eq.3.and.nprint.ne.-5)write(iwr,250)
      do 220 iset=1,ntype
      nthis=nsymm(iset)
      if(nthis.eq.0)go to 220
      do 30 i=1,nthis
      if(nblock(iset,i).gt.nprim)go to 40
30    continue
40    ip=i-1
      do 60 i=1,nthis
      ii=nblock(iset,i)
      do 60 j=1,i
      jj=nblock(iset,j)
      if(isymmo(ii).ne.isymmo(jj))callcaserr
     1('some sort of error with symmetry in canonicalisation')
      k=ic1e(ii,1)+ic1e(jj,2)
      f(i,j)=skv(k)
      if(j.le.ip)f(i,j)=0.0d0
      if(i.le.ip)f(i,j)=skv1(k)
      if(jj.le.ncore.and.ii.ne.jj)f(i,j)=0.0d0
60    f(j,i)=f(i,j)
c     nthiss=nthis
      ifail=0
      call f02abf(f,length,nthis,eig2,work,length,eig1,ifail)
c.....sort vectors
      call setsto(nthis,-100,iyc1)
c......select primary orbitals
      if(ip.eq.0)go to 110
      do 100 i=1,ip
      test=-100.0d0
      do 80 j=1,nthis
      if(iyc1(j).gt.0)goto80
      test1= dabs(work(i,j))
      if(test1.le.test)go to 80
      test=test1
      max=j
80    continue
      call dcopy(nthis,work(1,max),1,f(1,i),1)
      iyc1(max)=10
      eig1(i)=eig2(max)
100   continue
c.....select secondary orbitals
110   ip=ip+1
      do 140 i=ip,nthis
      test=1.0d10
      do 120 j=1,nthis
      if(iyc1(j).gt.0)go to 120
      if(eig2(j).ge.test)go to 120
      test=eig2(j)
      max=j
120   continue
      call dcopy(nthis,work(1,max),1,f(1,i),1)
      iyc1(max)=10
140   eig1(i)=eig2(max)
c.....check phases
      do 180 i=1,nthis
      test=0.0d0
      do 160 j=1,nthis
      iphase=1
      test1=f(j,i)
      if(test1.gt.0.0d0)go to 150
      test1=-test1
      iphase=-iphase
150   if(test1.le.test)go to 160
      test=test1
      k=iphase
160   continue
      if(k.gt.0)go to 180
      do 170 j=1,nthis
170   f(j,i)=-f(j,i)
180   continue
c
      calltrnsfo (q,work,f,iset,length)
      do 210 i=1,nthis
210   eig(nblock(iset,i))=eig1(i)
220   continue
      if(nprint.ne.-5)write(iwr,190)
190   format(/40x,'eigenvalues of fock matrix'/40x,26('-'))
      if(nprint.ne.-5)write(iwr,200)(eig(i),i=1,ma)
200   format(/10x,8f14.7)
      return
230   format(/' active orbitals to diagonalise density matrix'
     1/       ' secondary orbitals to diagonalise f(a)+f(i) matrix'/)
240   format(/' active orbitals to diagonalise lagrangian matrix'
     1/       ' secondary orbitals to diagonalise f(a)+f(i) matrix'/)
250   format(/' active orbitals to diagonalise f(a)+f(i) matrix'
     1/       ' secondary orbitals to diagonalise f(a)+f(i) matrix'/)
c
      end
      subroutine dippol(d,p,q,gam1,nprint)
c
c     compute  dipole moment of casscf wavefunction
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/prints)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab
     1                 ,nst,n1e,isym,mults,macro
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/qice)
      common/blkin /corev(512),potnuc,dx(3),nato,numo
      dimension d(*),p(*),q(*),xx(3),gam1(*),o1e(6)
      dimension dxt(3)
      data half,two,small/0.5d0,2.0d0,1.0d-16/
      data xx/'x','y','z'/
      data o1e/ .false., .false., .false., .true., .true., .true. /
c
      nbasp=nbas4+1
      lenb2=lenb4+lenb4
      call vclr(d,1,lenb4)
      do 20 it=1,nprim
      ilift=ilifq(it)
      do 20 iu=1,it
      if((iu.le.ncore.and.it.ne.iu).or.isymmo(it).ne.isymmo(iu))
     +    go to 20
      qqq=gam1(ic1e(it,1)+ic1e(iu,2))
      if(it.le.ncore)qqq=two
      ilifu=ilifq(iu)
      kl1=0
      do 10 k=1,ma
      kt=ilift+k
      ku=ilifu+k
      do 10 l=1,k
      kl1=kl1+1
10    d(kl1)=d(kl1)+qqq*(q(kt)*q(ilifu+l)+q(ilift+l)*q(ku))
20    continue
      do 3 i=2,nbasp
    3 d(iky(i))=d(iky(i))*half
c
c.....ao density matrix in d
c
      if(oprint(8))call writel(d,nbas4)
c
c     restore dipole integrals
c
      call getmat(p(1),p(1),p(1),p(1),p(lenb4+1),p(lenb2+1),
     +            potnuc,nbas4,o1e,ions2)
c
      if(nprint.ne.-5)write(iwr,80)
      m=1
      do 70 i=1,3
      dxe=ddot(lenb4,d,1,p(m),1)
      dxt(i)=dx(i)+dxe
      if(nprint.ne.-5)write(iwr,90)xx(i),dx(i),dxe,dxt(i)
70    m=m+lenb4
      dtot=ddot(3,dxt,1,dxt,1)
      if(dtot.gt.small)dtot= dsqrt(dtot)
      if(nprint.ne.-5)write(iwr,100)dtot
      return
80    format(/17x,'dipole moments'//11x,'nuclear',6x,'electronic',
     111x,'total')
90    format(1x,a1,3f16.7)
100   format(/' total dipole moment(a.u.)=',f16.7/)
c
      end
      subroutine dmprin(p2,gam2,iwr,nprint)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa
     +            ,nab,nst,n1e,isym,mults,macros,maxc
INCLUDE(common/qice)
INCLUDE(common/gjs)
      dimension p2(*),gam2(*)
c
      if(nprint.ne.-5)write(iwr,10)ncore
10    format(' print of two-particle density matrix',
     a' -- ncore =',i3)
      nptri=iky(nact+1)
      do 50 ii=1,nact
c     i=ii+ncore
      do 50 jj=1,ii
c     j=jj+ncore
      ij=iky(ii)+jj
      if(nprint.ne.-5)write(iwr,20)ii,jj,ic4e(ij)
20    format(' orbital pair',i3,',',i2,' - symmetry =',i2)
      mapout=0
      call vclr(p2,1,nptri)
      do 40 kk=1,ii
      k=kk+ncore
      isymk=mult(isymmo(k),ic4e(ij))
      map=iky(kk)
      do 30 ll=1,kk
      l=ll+ncore
      map=map+1
      mapout=mapout+1
      if(isymmo(l).ne.isymk)goto 30
      ic=ic3e(map)+ic2e(ij)
      if(map.gt.ij)ic=ic2e(map)+ic3e(ij)
      p2(mapout)=gam2(ic)
30    continue
40    continue
50    if(nprint.ne.-5)call writel(p2,ii)
c
      return
      end
      subroutine grcntl(q,fock,gam1,gam2,iwr)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/timez)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
      common/qpar/ma
INCLUDE(common/restar)
      dimension q(*),gam1(*),gam2(*),fock(*)
c
      call rdedx(q,nblkq4,iblkq4,ndump4)
      call tdown(q,ilifq,q,ilifq,ma)
      if(irest.eq.4)go to 1
      call moden(q,fock,gam1,gam2,iwr,nprint)
      call pi2(q,q(nblkq4+1),gam1,nprint)
1     call pi4(q)
      if(tim.gt.timlim)irest=4
c
      return
      end
      subroutine moden(q,qq,gam1,gam2,iwr,nprint)
c
c      construct the lagrangian matrix (matrix elements of the
c       mcscf fock operator), and then transform to the hondo
c       basis
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
      common/add/iad1(5),iad2
INCLUDE(common/prints)
INCLUDE(common/qice)
INCLUDE(common/dm)
INCLUDE(common/gjs)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab
     1            ,nst,n1e,isym,mults,macro
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
      dimension q(*),gam1(*),gam2(*),qq(*)
      data m401,m991/401,991/
      data dzero/0.0d0/
c
      nblkq4=ma*ma
      ione=iad2
      itwo=ione+nblkq4
c
c     nprim2=nprim*nprim
c
      call dpconv(nprim,gam1,gam2)
c
      if(ncore.eq.0)goto30
      do 20 i=1,nprim
      i1=ilifp(i)
      isymi=isymmo(i)
      indi=ic1e(i,1)
      max=ncore
      if(i.lt.ncore)max=i
      do 20 j=1,max
      j1=ilifp(j)
      qq(i1+j)=dzero
      if(isymmo(j).ne.isymi)goto20
      ind=ic1e(j,2)+indi
      qq(i1+j)=(q(nblkq4+ind)+q(nblkq4+n1e+ind))*2.0d0
20    qq(j1+i)=qq(i1+j)
c
30    do 50 i=nst,nprim
      i1=ilifp(i)
      isymi=isymmo(i)
      do 50 j=nst,nprim
      if(isymmo(j).ne.isymi)goto50
      qq(i1+j)=qq(i1+j)*2.0d0
50    continue
      do 60 i=nst,nprim
      i1=ilifp(i)
      do 60 j=nst,i
      j1=ilifp(j)
      qq(i1+j)=0.5d0*(qq(i1+j)+qq(j1+i))
60    qq(j1+i)=qq(i1+j)
c
c....fock matrix over mo's now constructed
c
      if(ifmola.eq.0)goto100
      map=ione-1
      do70i=1,nprim
      i1=ilifp(i)
      do70j=1,i
      map=map+1
70    q(map)=qq(i1+j)
      if(nprint.ne.-5)write(iwr,80)
80    format(40x,' lagrangian matrix over mo''s')
      if(nprint.ne.-5)call writel(q(ione),nprim)
      lent=iky(nprim+1)
      len=lensec(lent)
      call secput(iseclm,m401,len,iblk)
      call wrt3(q(ione),lent,iblk,ndump4)
      if(nprint.ne.-5)write(iwr,90)iseclm
90    format(/' symmetric lagrangian matrix (mo basis) dumped to'
     1,' section ',i3/)
100   continue
c
c.....now express in ao basis
c
      if(iflagr.ne.1)return
      do 120 i=1,nprim
      isymi=isymmo(i)
      i1=ilifp(i)
      i2=ilifq(i)
      map=itwo-1+i2
      do 120 il=1,nbas4
      map=map+1
      dum=dzero
      do 110 j=1,nprim
      if(isymmo(j).ne.isymi)goto110
      dum=dum+q(ilifq(j)+il)*qq(i1+j)
110   continue
120   q(map)=-dum
c
      do 140 il=1,nbas4
      map1=il+itwo-1
      map=iky(il)+ione-1
      do 140 im=1,il
      map=map+1
      dum=dzero
      do 130 i=1,nprim
      map2=map1+ilifq(i)
130   dum=dum+q(map2)*q(ilifq(i)+im)
140   q(map)=dum
c
      call secput(isecla,m991,lensec(lenb4),iblk)
      call wrt3(q(ione),lenb4,iblk,ndump4)
c     iend=ione+lenb4-1
       if(nprint.ne.-5.and.oprint(16))write(iwr,6)
 6     format(/40x,' symmetric lagrangian matrix (ao basis)'/)
       if(nprint.ne.-5.and.oprint(16))call writel(q(ione),nbas4)
      call revind
      if(nprint.ne.-5)write(iwr,150)isecla
150   format(//' symmetric lagrangian matrix (ao basis) dumped to'
     1,' section ',i3/)
      call gconv(nprim,gam1,gam2)
c
      return
      end
      subroutine p1out(d,gam1,gam2,iwr,nprint)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/cndx40)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     +              nst,n1e,isym,mults,macro,maxc
INCLUDE(common/qice)
INCLUDE(common/dm)
      dimension d(*),gam1(*),gam2(*)
      data m400/400/
c
c      output one-particle density matrix of casscf wavefunction
c      to atmol dumpfile
c
c
      call dpconv(nprim,gam1,gam2)
      lenb4=iky(nbas4+1)
      call secput(isecdm,m400,lensec(lenb4),iblk)
      call vclr(d,1,lenb4)
      ind=0
      do 10 i=1,nprim
      do 10 j=1,i
      ind=ind+1
      if(isymmo(j).ne.isymmo(i))goto10
      in1=ic1e(i,1)+ic1e(j,2)
      d(ind)=gam1(in1)
      if(i.eq.j.and.i.le.ncore)d(ind)=2.0d0
10    continue
      call wrt3(d,lenb4,iblk,ndump4)
      if(nprint.ne.-5)write(iwr,20)isecdm
20    format(/' one particle density matrix (mo basis) dumped to section
     1',i4,' of dumpfile')
      call gconv(nprim,gam1,gam2)
      return
c
      end
c
      subroutine p2out(gam1,gam2,iwr)
c
c      output two-particle density matrix of casscf wavefunction
c      to mainfile (the one given by the mofile directive
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc  /isel(5),name(maxlfn)
_IF1(iv)      common/craypk/ij205(340),kl205(340)
_IFN1(iv)      common/craypk/ij205(680)
      common/blkin /gout(510),mword
INCLUDE(common/qice)
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/files)
INCLUDE(common/restar)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     1            nst,n1e,isym,mults,macro,maxc
INCLUDE(common/machin)
      dimension gam1(*),gam2(*)
      data m0,two/0,2.0d0/
c
      mfilep=1
      mainp=n6tape(1)
      iblkmp=n6blk(1)
      mblp=iblkmp-n6last(1)
c
      mword=0
c
_IFN1(iv)      int4=1
      call dpconv(nprim,gam1,gam2)
      if(ncore.eq.0)goto60
c
c.....core-core contributions...
c
      do 20 i=1,ncore
      do 20 j=1,i
      do 20 k=1,i
      last=k
      if(k.eq.i)last=j
      do 20 l=1,last
      gam=0.0d0
      if(i.eq.j.and.k.eq.l)gam=two
      if(i.eq.k.and.j.eq.l)gam=-0.5d0
      if(i.eq.j.and.k.eq.l.and.i.eq.k)gam=1.0d0
      if( dabs(gam).lt.0.1d0)goto20
      mword=mword+1
_IF1(iv)      ij205(mword)=i4096(i)+j
_IF1(iv)      kl205(mword)=i4096(k)+l
_IFN1(iv)      ij205(int4  )=i4096(i)+j
_IFN1(iv)      ij205(int4+1)=i4096(k)+l
_IFN1(iv)      int4=int4+2
      gout(mword)=gam
      if(mword.lt.nintmx)goto 20
      call block2
_IFN1(iv)      int4=1
20    continue
c
c....core-active
c
      do 50 it=nst,nprim
      do 50 iu=nst,it
      if(isymmo(it).ne.isymmo(iu))goto50
      gam=gam1(ic1e(it,1)+ic1e(iu,2))
      do 40 i=1,ncore
      mword=mword+1
      gout(mword)=gam
_IF1(iv)      ij205(mword)=i4096(it)+iu
_IF1(iv)      kl205(mword)=i4096(i)+i
_IFN1(iv)      ij205(int4  )=i4096(it)+iu
_IFN1(iv)      ij205(int4+1)=i4096(i)+i
_IFN1(iv)      int4=int4+2
      if(mword.lt.nintmx)goto 30
      call block2
_IFN1(iv)      int4=1
30    mword=mword+1
      gout(mword)=-gam*0.25d0
_IF1(iv)      ij205(mword)=i4096(it)+i
_IF1(iv)      kl205(mword)=i4096(iu)+i
_IFN1(iv)      ij205(int4  )=i4096(it)+i
_IFN1(iv)      ij205(int4+1)=i4096(iu)+i
_IFN1(iv)      int4=int4+2
      if(mword.lt.nintmx) goto 40
      call block2
_IFN1(iv)      int4=1
40    continue
50    continue
c
c.....all active density matrix..
c
60    do 70 i=nst,nprim
      do 70 j=nst,i
      iij=iky(i-ncore)+j-ncore
      isymij=ic4e(iij)
      indij=ic2e(iij)
      do 70 k=nst,i
      last=k
      if(k.eq.i)last=j
      do 70 l=nst,last
      ikl=iky(k-ncore)+l-ncore
      if(ic4e(ikl).ne.isymij)goto70
      ind=indij+ic3e(ikl)
      mword=mword+1
      gout(mword)=gam2(ind)
_IF1(iv)      ij205(mword)=i4096(i)+j
_IF1(iv)      kl205(mword)=i4096(k)+l
_IFN1(iv)      ij205(int4  )=i4096(i)+j
_IFN1(iv)      ij205(int4+1)=i4096(k)+l
_IFN1(iv)      int4=int4+2
      if(mword.lt.nintmx)goto 70
      call block2
_IFN1(iv)      int4=1
70    continue
      if(mword.eq.0)go to 80
      call block2
80    call put(gout,m0,mainp)
      m6file=mfilep
      m6tape(m6file)=n6tape(mfilep)
      m6blk (m6file)=n6blk(mfilep)
      m6last(m6file)=iblkmp+1
      call gconv(nprim,gam1,gam2)
      call revise
      call clredx
      if(nprint.ne.-5)write(iwr,10)
 10   format(//
     *' status of 2-particle mo-density file'/1x,36('-')/)
      if(nprint.ne.-5)
     * call filprn(m6file,m6blk,m6last,m6tape)
      return
c
      end
c
      subroutine block2
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/files)
INCLUDE(common/atmblk)
      common/blkin/gout(510),mword
      common/craypk/ij205(680)
      common/disc/isel(5),ipos(maxlfn)
INCLUDE(common/restar)
      if(ipos(mainp).ne.iblkmp)callsearch(iblkmp,mainp)
_IF1(iv)      call pak4v(ij205,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,ij205,numlabp)
      call put(gout(1),m511,mainp)
      mword=0
      iblkmp=iblkmp+1
      mblp=mblp+1
      if(mblp)41,40,41
c
c...   change channel
c
 40    m6file=mfilep
       m6tape(m6file)=n6tape(mfilep)
       m6blk(m6file)=n6blk(mfilep)
       m6last(m6file)=iblkmp
       mfilep=mfilep+1
      if(mfilep.gt.n6file)call errors(890)
       mainp=n6tape(mfilep)
      iblkmp=n6blk(mfilep)
      mblp=iblkmp-n6last(mfilep)
  41  return
      end
c
      subroutine pcalc1(q,iqq,qq)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
c
c...   processes sorted mainfile on lfn sort
c...   to produce integrals of the form (i j/r s)
c...   on secondary mainfile
c
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
      dimension iqq(*)
       dimension q(*),qq(*)
INCLUDE(common/gjs)
      common/junke/maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/junk/nwbuck(1500),inum(maxorb),jnum(maxorb),
     * p(maxorb),v(maxorb)
      common/bufb/nkk,mkk,gin(1)
_IF1(iv)      common/craypk/ij205(340),kl205(340)
_IFN1(iv)       common/craypk/ij205(680)
      common/blkin/gout(510),mword
      common/three/mark(1500)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      data maxb/9999999/
c
      call rdedx(q(1),nblkq4,iblkq4,ndump4)
      call tdown(q,ilifq,q,ilifq,nbas4)
      nmc1=nmc+1
      length=iky(nmc1)
c
c......transpose q
c
      call dcopy(nblkq4,q,1,q(nblkq4+1),1)
      do 90 i=1,nbas4
      map=ilifq(i)
      do 90 j=1,nbas4
90    q(map+j)=q(nblkq4+ilifq(j)+i)
      call vclr(p,1,nbas4)
      mword=0
_IFN1(iv)      int4=-1
      small=10.0d0**(-iacc-2)
      call stopbk
      do 1 ibuck=1,nbuck
      mhigh=mloww+ntri
      if(mhigh.gt.mhi)mhigh=mhi
      mtri=mhigh-mloww
      mloww=mloww+1
      nsize=mtri*lentri
      call vclr(qq,1,nsize)
c
c...   read in a core load
c
      mkk=mark(ibuck)
      go to 999
 888  iblock=mkk
      call rdbak(iblock)
      call stopbk
      call transb(qq)
 999  if(mkk.ne.maxb)go to 888
      map=0
      do 777 itri=1,mtri
      master=master+1
      if(master.gt.length)goto777
c
c...   master=iky(r)+s
c...   compute integrals of the form (i j/r s)
c... evaluate sparcity of triangle
c
      do 4000 i=1,nmc
      m=iky(i)
      nsize=m+map
      nlim=0
_IFN1(cx)      do 3 j=1,i
_IFN1(cx)      if(qq(nsize+j))4,3,4
_IFN1(cx)    4 nlim=nlim+1
_IFN1(cx)      iqq(m+nlim)=j
_IFN1(cx)    3 continue
_IF1(c)      call whenne(i,qq(nsize+1),1,0.0d0,iqq(m+1),nlim)
_IF1(x)      call dlstne(i,qq(nsize+1),1,0.0d0,nlim,iqq(m+1))
4000  inum(i)=nlim
      do 1001 k=1,nbas4
      call dcopy(nbas4,q(ilifq(k)+1),1,v(1),1)
      do 1003 i=1,nmc
      nlim=inum(i)
      top=0.0d0
      if(nlim)1003,1003,1002
1002  m=iky(i)
      nsize=map+m
      vv=v(i)
      if(vv)6000,6001,6000
6001   do 6002 jj=1,nlim
       j=iqq(jj+m)
6002   top=qq(nsize+j)*v(j)+top
      goto 1003
6000  do 1004 jj=1,nlim
      j=iqq(jj+m)
      bot=qq(nsize+j)
      top=bot*v(j)+top
1004  p(j)=bot*vv+p(j)
1003  p(i)=top
c
c... evaluate sparcity of x
c
      nlim=0
      do 1006 i=1,nmc
      if( dabs(p(i)).lt.small)goto1006
      nlim=nlim+1
      jnum(nlim)=i
1006  continue
      if(nlim)1001,1001,1007
1007  m=iky(k)
      do 1009 j=1,k
      top=0.0d0
      jj=ilifq(j)
      do 1008 i=1,nlim
      l=jnum(i)
1008  top=p(l)*q(jj+l)+top
c
c...    (i j/r s) now in top
c
      if( dabs(top).lt.small)goto1009
      mword=mword+1
      gout(mword)=top
_IF1(iv)      ij205(mword) = m+j
_IF1(iv)      kl205(mword) =master
_IFN1(iv)      int4=int4+2
_IFN1(iv)      ij205(int4  )=m+j
_IFN1(iv)      ij205(int4+1)=master
      if(mword.lt.nintmx)go to 1009
_IF1(iv)      call pak4v(ij205,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,ij205,numlabp)
_IFN1(iv)      int4=-1
      call blockk
1009  continue
1001  continue
777   map=map+lentri
1      mloww=mhigh
      if(mword.lt.1)return
_IF1(iv)      call pak4v(ij205,gout(num2ep+1))
_IFN1(iv)      call pack(gout(num2ep+1),lab1632,ij205,numlabp)
_IFN1(iv)      int4=-1
      call blockk
c
      return
c
      end
c
      subroutine pcalc2(imap1,q,qq,qqq,lokal,mxtyp)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/restar)
      dimension mtemp0(48),mtemp1(48),m2(48),m3(48)
      dimension lokal(*)
      dimension imap1(*)
      dimension q(*),qq(*),qqq(*)
c
c...  processes sorted secondary mainfile on lfn sort
c...  modified to produce density matrices over ao's in shells
c
      common/qpar/ma,n2e,ntypee
INCLUDE(common/symtry)
      common/lsort/iso(1)
      common/junk/nwbuck(1500),llin(6800),gtx(3400),
     +            lengt(mxshel),ish,jsh,maxtrs,lentrs,ishtri,
     +            icentr(mxshel)
INCLUDE(common/gjs)
      common/junke/maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/bufb/nkk,mkk,gin(1)
INCLUDE(common/mapper)
      common/blkin/gout(510),mword
_IF1(iv)      common/craypk/ij205(340),kl205(340)
_IFN1(iv)      common/craypk/ij205(680)
      common/three/mark(3000),p(maxorb),v(maxorb)
      data maxb/9999999/
      nm=maxtrs*nmc*mxtyp
      lentri=iky(nmc+1)
      acc1=10.0d0**(-iacc)
      call rdedx(q(1),nblkq4,iblkq4,ndump4)
      call tdown(q,ilifq,q,ilifq,nbas4)
c
c...... transpose q
c
      call dcopy(nblkq4,q,1,q(nblkq4+1),1)
      do 90 i=1,nbas4
      map=ilifq(i)
      do 90 j=1,nbas4
90    q(map+j)=q(nblkq4+ilifq(j)+i)
      call vclr(p,1,nbas4)
      mword=0
_IFN1(iv)      int4=-1
       call stopbk
      do 1 ibuck=1,nbuck
c
c.....lokal(ij) to give the position of real orbital-pair ij
c     in its shell-pair
c
      call setsto(lentrs,0,llin)
      do 10 i=1,nbas4
      ii=iky(i)
      ishell=iky(ilifm(i))
      do 10 j=1,i
      jshell=ishell+ilifm(j)
      ii=ii+1
      lokal(ii)=llin(jshell)
10    llin(jshell)=llin(jshell)+1
      mhigh=mloww+ntri
      if(mhigh.gt.mhi)mhigh=mhi
      mtri=mhigh-mloww
      mloww=mloww+1
      nsize=mtri*lentri*maxtrs
      call vclr(qq,1,nsize)
c...   read in a core load
c
      mkk=mark(ibuck)
c
      go to 999
888   iblock=mkk
       call rdbak(iblock)
      call stopbk
c
c...   move block to processing area
c
      call transe(qq,imap1,lokal,maxtrs)
 999  if(mkk.ne.maxb)go to 888
      nmc6=nmc*mxtyp
      mappp=1
      do 777 ishtr=1,mtri
      ishtri=ishtri+1
      do920 it=1,nt
      id=iso(indxi+iliso(it))
      if(id.gt.indxi)goto1002
920   mtemp0(it)=id
      do 950 it=1,nt
      id=mtemp0(it)
      jd=iso(indxj+iliso(it))
      if(jd.gt.indxi)goto1002
      if(id.ge.jd)goto940
      nd=id
      id=jd
      jd=nd
940   if(id.eq.indxi.and.jd.gt.indxj)goto1002
      mtemp1(it)=id
950   m2(it)=jd
      do 1001 ksh=1,indxi
      do970 it=1,nt
      kd=iso(ksh+iliso(it))
      if(kd.gt.indxi)goto1001
970   m3(it)=kd
      ocomcn=icentr(indxi).eq.icentr(indxj).and.
     +       icentr(indxi).eq.icentr(ksh)
      mapp=mappp
      call vclr(qqq,1,nm)
      mapit=1
      kbot=mapie(ksh)
      ktop=kbot+lengt(ksh)-1
      do 500 i=1,nbas4
      itri=iky(i)
      do 500 j=1,i
      itri=itri+1
      if(imap1(itri).ne.ishtri)goto500
      maput=mapit
      do 400 k=kbot,ktop
      call dcopy(nmc,q(ilifq(k)+1),1,v(1),1)
      map=mapp
      mapotr=maput
      do 360 ir=1,nmc
      top=0.0d0
      mapout=maput
      vv=v(ir)
      irm=ir-1
      if(irm)361,361,351
351   do 350 is=1,irm
      top=top+qq(map)*v(is)
      qqq(mapout)=qq(map)*vv+qqq(mapout)
      map=map+1
350   mapout=mapout+1
361   top=top+qq(map)*vv
      qqq(mapotr)=top
      map=map+1
      mapotr=mapotr+1
360   continue
      maput=maput+nmc
400   continue
      mapit=mapit+nmc6
      mapp=mapp+lentri
500   continue
c
c.....qqq holds alk (ij|ks) in shell triple ijk
c
      do 200 lsh=1,ksh
      if(ocomcn.and.icentr(lsh).eq.icentr(ksh))goto200
      do110 it=1,nt
      ld=iso(lsh+iliso(it))
      kd=m3(it)
      if(kd.ge.ld)goto990
      nd=kd
      kd=ld
      ld=nd
990   id=mtemp1(it)
      jd=m2(it)
c...    since ii is rubbish this is  correctly a do-nothing (mfg/jvl 88)
      if(id.ne.ii.and.kd.ne.ii)goto110
      if(kd.lt.id.or.(kd.eq.id.and.ld.le.jd))goto9100
      nd=id
      id=kd
      kd=nd
      nd=jd
      jd=ld
      ld=nd
9100  if(jd.gt.indxj.or.kd.gt.ksh.or.ld.gt.lsh)goto200
110   continue
      lbot=mapie(lsh)
      ltop=lbot+lengt(lsh)-1
      mapit=1
      do 250 i=1,nbas4
      do 250 j=1,i
      itri=iky(i)+j
      if(imap1(itri).ne.ishtri)goto250
      mword1=i4096(i)+j
      maput=mapit
      do 240 k=kbot,ktop
      last=ltop
      if(ksh.eq.lsh)last=k
      if(k-i)219,218,219
218   if(lsh.gt.indxj)goto240
      if(lsh.eq.indxj)last=j
219   do 230 l=lbot,last
      call dcopy(nmc,q(ilifq(l)+1),1,v(1),1)
      top=ddot(nmc,qqq(maput),1,v(1),1)
       if( dabs(top).lt.acc1)goto230
      mword=mword+1
      gout(mword)=top
_IF1(iv)      ij205(mword) =mword1
_IF1(iv)      kl205(mword) =l+i4096(k)
_IFN1(iv)      int4=int4+2
_IFN1(iv)      ij205(int4  )=mword1
_IFN1(iv)      ij205(int4+1)=l+i4096(k)
      if(mword.lt.nintmx)go to 230
      call blockp
_IFN1(iv)      int4=-1
230   continue
240   maput=maput+nmc
      mapit=mapit+nmc6
250   continue
200   continue
1001  continue
1002  indxj=indxj+1
      if(indxj.le.indxi)goto 777
      indxi=indxi+1
      indxj=1
777   mappp=mappp+lentri*maxtrs
1      mloww=mhigh
      if(mword.lt.1)return
      call blockp
      return
c
      end
c
      subroutine pdumpr(iphase,ipass,npass,jbl,jun,iwr)
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/timez)
INCLUDE(common/restar)
INCLUDE(common/files)
      common/trntim/told,alimin
INCLUDE(common/cndx40)
      common/blkin/gout(511)
      character*10 charwall
      data m0/0/
      if(ipass.ne.npass)goto 88
      call search(jbl,jun)
      call put(gout,m0,jun)
      jbl=jbl+1
 88    if(iphase-1)89,89,90
 89    m4file=nfiles
       m4tape(m4file)=n4tape(nfiles)
       m4blk (m4file)= n4blk( nfiles)
       m4last(m4file)= jbl
       go to 66
 90    m11fil=nfilef
       m11tap(m11fil)=n11tap(nfilef)
       m11bl (m11fil)= n11bl (nfilef)
       m11lst(m11fil)= jbl
66     call search(isecbl,ndump4)
       m12=12/lenwrd()
      call put(master,m12,ndump4)
      b=cpulft(1)
      call clredx
      if(nprint.ne.-5)write(iwr,100)iphase,ipass,b ,charwall()
100   format(/' job dumped in sort',i1,' pass',i6,' at ',f8.2,
     *' seconds',a10,' wall')
       if(ipass.eq.npass.and.iphase.eq.2)goto 99
      if(((b-told)*2.0d0+b).lt.alimin)goto 77
      write(iwr,101)
 101  format(/10x,50('=')/
     *10x,'density matrix transformation incomplete - restart'/
     *10x,50('=')/)
       call secsum
       call whtps
       tim=timlim+0.1d0
 77    told=b
99    return
c
      end
c
      subroutine pi2(vector,d,gam1,nprint)
c
c.....transforms mo 1-particle density matrix to ao basis
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension vector(ma,*),d(*),gam1(*)
INCLUDE(common/qice)
INCLUDE(common/cndx40)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     +           naa,nab,nst,n1e,isym,mults
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/dm)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      data half/0.5d0/
      data m990/990/
c
      lenb4=ma*(ma+1)/2
      call vclr(d,1,lenb4)
c
      do 20 it=1,nprim
      do 20 iu=1,it
      if(isymmo(iu).ne.isymmo(it)) go to 20
      if(iu.le.ncore.and.it.ne.iu)goto20
      m1=ic1e(it,1)+ic1e(iu,2)
      dum=gam1(m1)*half
      if(it.le.ncore)dum=1.0d0
      do 10 k=1,ma
      kl1=iky(k)
      do 10 l=1,k
      kl1=kl1+1
10    d(kl1)=d(kl1)+dum*(vector(k,it)*vector(l,iu)+
     1vector(k,iu)*vector(l,it))
c
20    continue
c
      call secput(isecda,m990,lensec(lenb4),iblk)
      call wrt3(d,lenb4,iblk,ndump4)
      call revind
      if(nprint.ne.-5)write(iwr,30)isecda
30    format(/' one-particle density matrix (ao basis ) dumped'
     1,' to section ',i3)
c
c ----- add to standard section
c
      call wrt3(d,lenb4,ibl3pa,ndump4)
      return
c
      end
c
      subroutine pi4(q)
c
c..... control routine for 4-index transformation of two-particle
c      density matrix from mo basis to gamess ao basis. calls
c      psort1, pcalc1, psort2, pcalc2.  call shellz gets the
c      gamess shell structure
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/junk/nwb(1500),itx(6800),gtx(3400),
     +            lengt(mxshel),ish,jsh,maxtrs,lentrs,ishtri,
     +            icentr(mxshel)
INCLUDE(common/blksiz)
INCLUDE(common/iofile)
INCLUDE(common/gjs)
INCLUDE(common/timez)
INCLUDE(common/restar)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/filel)
INCLUDE(common/machin)
      common/disc/isel,iselr,iselw,irep,ichek,ipos(maxlfn)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,lentri,
     *nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
      common/craypk/ijkl(1360)
INCLUDE(common/files)
INCLUDE(common/atmblk)
      common/blkin/gout(511)
      common/three/mark(1500),ibase(1500),ibasen(1500)
      common/trntim/top,alimin
      character*10 charwall
      dimension q(*)
      data m500/500/
      data m1,m2/1,2/
c
      nav = lenwrd()
      ishtri=0
      lentri=iky(nmc+1)
      top=cpulft(1)
      if(nprint.ne.-5)write(iwr,401)top ,charwall()
401   format(/' commence 2-particle density transformation at ',
     1f8.2,' seconds',a10,' wall'/1x,68('=')/)
c
      if(irest.eq.4)go to 3
c
      call secput(isect(485),m500,m1,isecbl)
      if(nprint.ne.-5)write(iwr,4)isect(485)
4     format(' 4-index restart information dumped to section',i4)
      nfiles=1
      nfilef=1
      indxi=1
      indxj=1
      irrr=1
      iss=1
      master=0
      junits=n4tape(1)
      jblkas=n4blk(1)
      jblkaf=jblkas-n4last(1)
      junitf=n11tap(1)
      jblkaf=n11bl(1)
      jblkrf=jblkaf-n11lst(1)
      goto5
3     callsecget(isect(485),m500,isecbl)
      m12=12/nav
      call readi(master,m12*nav,isecbl,ndump4)
      if(nprint.ne.-5)write(iwr,6)isect(485)
6     format(/' 4-index information recovered from section',i4)
5     continue
      call ibasgn(1500,0,nsz340,ibase)
      call ibasgn(1500,0,nsz680,ibasen)
      niqq=lenb4+nblkq4
c     len=iky(nmc+1)
      maxt=(lword4-niqq)/lentri
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
      if(ires.lt.m1.or.maxt.lt.m1)go to 3000
      if(ires.gt.1500)ires=1500
       if(master.eq.lentri)goto 7766
c
c...   determine min. no. of passes for sort1/calc1
c
      i=lentri-master-1
      npass1=1
102   nteff=i/npass1+1
      if(((nteff-1)/ires).lt.maxt)goto 103
      npass1=npass1+1
      goto 102
103   if(npass1.gt.i)npass1=i+1
       if(nprint.ne.-5)write(iwr,204)m1,npass1
      lfile=m6file
      do 1000 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i) =m6blk(i)
 1000 llblk(i)=liblk(i)-m6last(i)
      do 1 ipass=1,npass1
      call setsto(1360,0,ijkl)
      call psort1(q,q(nword+1))
       call pcalc1(q,q(nblkq4+1),q(niqq+1))
       call pdumpr(m1,ipass,npass1,jblkas,junits,iwr)
       if(tim.gt.timlim)go to 2223
 1     continue
7766  continue
      lfile=m4file
      do 9998 i=1,lfile
      lotape(i)=m4tape(i)
      liblk (i)=m4blk (i)
9998  llblk(i)=liblk(i)-m4last(i)
      if(nprint.ne.-5)write(iwr,4007)
4007  format(/' status of secondary mainfile')
      if(nprint.ne.-5)call filprn(m4file,m4blk,m4last,m4tape)
      callshellz (q,mxtyp)
      niqqq=lenb4*maxtrs+nblkq4
      niqq=niqqq+max(maxtrs*mxtyp*nmc,lenb4)
      maxt=(lword4-lenb4-niqq)/(lentri*maxtrs)
_IF(ibm,vax)
      nword=((lword4-lenb4)*2)/3
      nword=(nword/4)*4
      nword2=nword/4
_ELSE
      if(o255i) then
       nword=(lword4-lenb4)/3
      else
       nword=(lword4-lenb4)/2
      endif
      nword=(nword/2)*2
      nword2=nword
_ENDIF
      maxt1=nword/(lentri*maxtrs)
      if(maxt.gt.maxt1)maxt=maxt1
      ires=nword/nsz340
      if(ires.lt.m1.or.maxt.lt.m1)go to 3000
      if(ires.gt.1500)ires=1500
      if(maxt.ge.1)go to 3010
 3000 call caserr(
     *'insufficient memory for 2p-density transformation')
 3010 i=lentrs-iky(indxi)-indxj
      npass2=1
202   nteff=i/npass2+1
      if(((nteff-1)/ires).lt.maxt)goto 203
      npass2=npass2+1
      goto 202
203   if(npass2.gt.i)npass2=i+1
       if(nprint.ne.-5)write(iwr,204)m2,npass2
204   format(/' no. of sort',i1,' passes=',i4)
       do 2 ipass=1,npass2
       call psort2(q,q(lenb4+1),q(nword+lenb4+1))
       call pcalc2(q,q(lenb4+1),q(niqq+lenb4+1)
     1,q(niqqq+lenb4+1),q(niqqq+lenb4+1),mxtyp)
       call pdumpr(m2,ipass,npass2,jblkaf,junitf,iwr)
       if(tim.gt.timlim)go to 2223
 2     continue
      top=cpulft(1)
       call clredx
      if(nprint.ne.-5)write(iwr,300)top ,charwall()
300    format(/' 2-particle density transformation complete at ',
     * f8.2,' seconds',a10,' wall'//
     *' status of 2-particle density file'/1x,33('*'))
      if(nprint.ne.-5)call filprn(m11fil,m11bl,m11lst,m11tap)
      irest=0
2223  call revise
      return
      end
      subroutine psort1(g,nijkl)
c
c...   sorts mainfile onto lfn sort --- so that for a
c...   given rs comb. al pq combs. available
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
      dimension g(*),nijkl(*)
INCLUDE(common/sizes)
      common/junke /maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     2nbuck,mloww,mhi,ntri
INCLUDE(common/cndx40)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
      common/bufb/nwbnwb,lnklnk,gout(1)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/mapper)
INCLUDE(common/filel)
      common/blkin/gin(1)
_IF(linux)
      external fget
_ENDIF
      data maxb/9999999/
c
c...   determine base and limit triangles for this pass
c
      mloww=master
      mhi=master+nteff
c
c...   determine minimum no. of bucks.
c
      if(mhi.gt.lentri)mhi=lentri
      mtri=mhi-mloww
      nbuck=ires
10    ntri=(mtri-1)/nbuck+1
      if(ntri.gt.maxt)goto 20
      nbuck=nbuck-1
      if(nbuck)10,20,10
c
c...   ntri=max. no. of triangles controlled by 1 bucket
c...   nbuck=number of buckets
c
20    nbuck=nbuck+1
      ntri=(mtri-1)/nbuck+1
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
      mlow=mloww+1
c
c...   start loop over mainfile blocks
c
      do 140 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)40,140,40
 40   iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
50    call fget(gin,m,iunit)
      if(m)60,140,60
 60   nnn=isort1(itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.2721)
     * call stackr(g,nijkl)
      lbl=lbl+1
      if(lbl)50,140,50
 140  continue
c
c...   mainfile now swept
c...   clear up output
c
       if(nstack.ne.0)
     * call stackr(g,nijkl)
      do 180 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)150,180,150
150   ib=ibase(ibuck)
      ibn=ibasen(ibuck)
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
c
      subroutine psort2(imap1,g,nijkl)
c
c...   sorts secondary mainfile onto lfn sort---so that
c...   for a given ij comb. all rs combs. available
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
INCLUDE(common/sizes)
      dimension g(*),nijkl(*)
      common/junk/nwbuck(1500),itx(3400),ktx(3400),gtx(3400),
     *lengt(mxshel),ish,jsh,maxtrs,lentrs
INCLUDE(common/mapper)
      common/junke /maxt,ires,ipass,
     1nteff,npass1,npass2,lentri,
     2nbuck,mloww,mhi,ntri,iacc
INCLUDE(common/cndx40)
      common/bufb/nwbnwb,lnklnk,gout(1)
INCLUDE(common/blksiz)
INCLUDE(common/stak)
      common/blkin/gin(1)
      common/three/mark(1500),ibase(1500),ibasen(1500)
INCLUDE(common/atmblk)
INCLUDE(common/filel)
      dimension imap1(*)
_IF(linux)
      external fget
_ENDIF
      data maxb/9999999/
c
c...   determine base and limit shell-triangles for this pass
c
      mloww=iky(indxi)+indxj-1
      mtri=(lentrs-mloww)/(npass2-ipass+1)
      mhi=mloww+mtri
c
c...   determine minimum no. of bucks.
c
      nbuck=ires
10    ntri=(mtri-1)/nbuck+1
      if(ntri.gt.maxt)goto 20
      nbuck=nbuck-1
      if(nbuck)10,20,10
20    nbuck=nbuck+1
      ntri=(mtri-1)/nbuck+1
c
c...   ntri=max. no. of shell-triangles controlled by 1 bucket
c...   nbuck=number of buckets
c
      btri=ntri
      btri=1.000000001d0/btri
      call setsto(nbuck,maxb,mark)
      call setsto(nbuck,0,nwbuck)
      nstack=0
      nwbnwb=nsz340
      iblock=0
      mlow=mloww+1
c...   start loop over secondary mainfile blocks
      do 130 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      lbl=llblk(ifile)
 140  call fget(gin,m,iunit)
      if(m)150,130,150
 150  nnn=isort3(imap1,itx(nstack+1),ktx(nstack+1),gtx(nstack+1))
      nstack=nstack+nnn
      if(nstack.ge.3061)
     * call stackp(g,nijkl,imap1)
      lbl=lbl+1
      if(lbl)140,130,140
 130  continue
c
c...   mainfile now swept
c...   clear up output
c
      if(nstack.ne.0)
     * call stackp(g,nijkl,imap1)
c
      do 180 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if(nwb)160,180,160
160   ib=ibase(ibuck)
      ibn=ibasen(ibuck)
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
c
      subroutine shellz(imap,mxtyp)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/junk/nwbuck(1500),itx(6800),gtx(3400),
     *lengt(mxshel),ish,jsh,maxtrs,lentrs,ishtri,
     1icent(mxshel)
INCLUDE(common/cndx40)
INCLUDE(common/mapper)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
      common/lsort/iso(1)
c
      dimension imap(*)
c
c ----- lengt  ..  no. of basis fumctions in shell i
c ----- mapie  ..  1st basis function in shell i
c ----- icent  ..  nucleus of shell i
c ----- ilifm  ..  shell no. of basis function i
c ----- maxtrs ..  largest no. of triangles in a shell pair
c ----- mxtyp  ..  largest no. of orbitals  in a shell
c
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ndump4)
      do 1 i=1,nshell
      icent(i)=katom(i)
c
      mapie(i)=kloc(i)
      mini=kmin(i)
      maxi=kmax(i)
      loci=kloc(i)-mini
      do 1 j=mini,maxi
      ii=loci+j
c
 1    ilifm(ii)=i
      maxtrs=0
      do 2 i=1,nshell
      j=kmax(i)-kmin(i)+1
      lengt(i)=j
      if(maxtrs.lt.j)maxtrs=j
 2    continue
      mxtyp = maxtrs
      if( mxtyp.lt.6) mxtyp = 6
      j=0
      do 3 i=1,nshell
      mini=kmax(i)-kmin(i)-1
      if(mini.lt.maxtrs)go to 3
      j=j+1
 3    continue
      if(j.eq.1)maxtrs=iky(maxtrs+1)
      if(j.ne.1)maxtrs=maxtrs*maxtrs
c
      lentrs=iky(nshell+1)
c
      do 4 i=1,nbas4
      map=iky(i)
      mapp=iky(ilifm(i))
      do 4 j=1,i
      map=map+1
 4    imap(map)=mapp+ilifm(j)
      return
c
      end
c
      subroutine blockp
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      common/craypk/integ(680)
      common/blkin/gout(510),mword
INCLUDE(common/files)
INCLUDE(common/cndx40)
INCLUDE(common/restar)
INCLUDE(common/atmblk)
      common/disc/isel(5),ipos(maxlfn)
c
_IF(ibm,vax)
      call pak4v(integ,gout(num2ep+1)))
_ELSE
      call pack(gout(num2ep+1),lab1632,integ,numlabp)
_ENDIF
      if(ipos(junitf).ne.jblkaf)call search(jblkaf,junitf)
      call put(gout(1),m511,junitf)
      mword=0
      jblkaf=jblkaf+1
      jblkrf=jblkrf+1
      if(jblkrf)41,40,41
c
c     change channel
c
 40   m11fil=nfilef
      m11tap(m11fil)=n11tap(nfilef)
      m11bl (m11fil)=n11bl  (nfilef)
      m11lst(m11fil)=jblkaf
      nfilef=nfilef+1
      if(nfilef.gt.n11fil)call errors(889)
      junitf=n11tap(nfilef)
      jblkaf=n11bl(nfilef)
      jblkrf=jblkaf-n11lst(nfilef)
41    return
      end
      subroutine p1den(vector,d,gam1,nprint)
c
c.....transforms mo 1-particle density matrix to ao basis
c
      implicit REAL (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension vector(ma,*),d(*),gam1(*)
INCLUDE(common/qice)
INCLUDE(common/cndx40)
      common/qpar  /ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,
     1            naa,nab,nst,n1e,isym,mults
INCLUDE(common/gjs)
INCLUDE(common/mapper)
INCLUDE(common/dm)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      data half/0.5d0/
      data m990/990/
c
      call rdedx(vector,nblkq4,iblkq4,ndump4)
      call tdown(vector,ilifq,vector,ilifq,ma)
      lenb4=ma*(ma+1)/2
      call vclr(d,1,lenb4)
c
      do 20 it=1,nprim
      do 20 iu=1,it
      if(isymmo(iu).ne.isymmo(it)) go to 20
      if(iu.le.ncore.and.it.ne.iu)goto20
      m1=ic1e(it,1)+ic1e(iu,2)
      dum=gam1(m1)*half
      if(it.le.ncore)dum=1.0d0
      do 10 k=1,ma
      kl1=iky(k)
      do 10 l=1,k
      kl1=kl1+1
10    d(kl1)=d(kl1)+dum*(vector(k,it)*vector(l,iu)+
     1vector(k,iu)*vector(l,it))
c
20    continue
c
      call secput(isecda,m990,lensec(lenb4),ibl990)
      call wrt3(d,lenb4,ibl990,ndump4)
      call revind
      if(nprint.ne.-5)write(iwr,30)isecda
30    format(/' one-particle density matrix (ao basis ) dumped'
     1,' to section ',i3)
c
c ----- add to standard section
c
      call wrt3(d,lenb4,ibl3pa,ndump4)
      return
c
      end
      subroutine schmid (v,n,n1,nvec,smin,mincol)
      implicit REAL  (a-h,p-x),integer (i-n),logical  (o)
      implicit character *8 (z)
      implicit character *4 (y)
      dimension v(*)
cibm  data small/z0210000000000001/
      smin = 1.0d30
_IF(c90,convex)
      small=x02amf(dum)
_ELSE
      small=x02agf(dum)
_ENDIF
      ia = (n1-1)*n + 1
      n2 = n1 + nvec -1
      do 40 k=n1,n2
      km = k-1
      if (km.le.0) goto 20
      ib = 1
      do 10 i=1,km
      s =-ddot(n,v(ib),1,v(ia),1)
      call daxpy(n,s,v(ib),1,v(ia),1)
10    ib = ib +n
20    continue
      s = dnrm2(n,v(ia),1)
      if (s.lt.small) call caserr
     + ('schmidt orthogonaliser has hit a singularity')
      if (s.gt.smin) goto 30
      smin = s
      mincol = k
30    s = 1.0d0/s
      call dscal(n,s,v(ia),1)
40    ia = ia +n
      return
      end
      subroutine reardr(ijgrp,ijadd,kadd,ladd,nij,nkl,nmax,itape8)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
      common/lsort /norbs,nlevs,na,nb,nc,nsym,msym,idocc,isym(maxorb),
     *icode(maxorb),
     1nlcs(maxorb),levpt(maxorb),levnr(maxorb),iout(maxorb)
     2,iop(8),nbf
      common/casscg/ncore,nact,naa,nbb,itypec(maxorb)
c....should be enough for 256 bfns, 15 active orbitals
      common/three /ic1e1(maxorb),ic1e2(maxorb),icf(maxorb),ic2e(121),
     1 ic3e(121),    ic4e(121)
      common/junk  /itabb(7200)
_IF(linux)
      integer xor
_ENDIF
      dimension ijadd(nij),ijgrp(nij),kadd(nkl),ladd(nkl)
      dimension id(8),iperm(8,8)
      data iperm/1,2,3,4,5,6,7,8,
     1           2,1,4,3,6,5,8,7,
     2           3,4,1,2,7,8,5,6,
     3           4,3,2,1,8,7,6,5,
     4           5,6,7,8,1,2,3,4,
     5           6,5,8,7,2,1,4,3,
     6           7,8,5,6,3,4,1,2,
     7           8,7,6,5,4,3,2,1/
c
      ii=1
      nst=ncore+1
      nprim=nact+ncore
c
      write(iwr,10)
10    format(//23x,'construct casscf integral map'/23x,29('-'))
c
      do 20 i=1,maxorb
20    icf(i)=(i*(i-1))/2
c
c.....set up the casscf integral look-ups
c.....1-electron
      n1e=0
      do 30 i=1,8
30    id(i)=0
      do 50 i=1,nprim
      iz=itypec(i)
      ic1e1(i)=n1e
      id(iz)=id(iz)+1
      ic1e2(i)=id(iz)
      do 40 j=1,i
40    if(itypec(j).eq.iz)n1e=n1e+1
50    continue
c
      ic1e1(nprim+1)=n1e
c
c.....2-electron
      do 60 i=1,8
60    id(i)=0
      kt=0
      do 70 i=1,nact
      do 70 j=1,i
      kt=kt+1
      iz=iperm(itypec(i+ncore),itypec(j+ncore))
      id(iz)=id(iz)+1
      ic4e(kt)=iz
70    ic3e(kt)=id(iz)
      n2e=0
      do 90 i=1,kt
      do 80 j=1,i
80    if(ic4e(i).eq.ic4e(j))n2e=n2e+1
90    ic2e(i+1)=n2e
      ic2e(1)=0
c
c
c.....now adjust ic1e to look up active 1-e ints at the end of the 2-e
      ioffst=n2e-ic1e1(nst)
      do 100 i=nst,nprim
100   ic1e1(i)=ic1e1(i)+ioffst
c
c
      n12e=ic1e1(nprim+1)+ioffst
      write(iwr,110)n2e,n1e
110   format(' number of active  2-electron integrals =',i6/
     1       ' number of primary 1-electron integrals =',i6)
      if(n12e.le.7200)goto130
      write(iwr,120)
120   format(' too many integrals')
      call caserr('dimensioning problem in reardr')
130   continue
c      do 754 i=1,norbs
      do 140 i=1,nmax
140   itabb(i)=-1
c
c.....set up mapping of 1-e integrals
      do 170 i=1,nbf
      do 160 j=1,i
      if (iout(i).lt.1.or.iout(j).lt.1) go to 160
      i1=iout(i)
      j1=iout(j)
      if(isym(i1).ne.isym(j1))goto160
      if (i1.ge.j1) go to 150
      m1=i1
      i1=j1
      j1=m1
150   if (ijgrp(i1*(i1+1)/2).ne.ii) go to 160
      m1=2
      if (i1.ne.j1) m1=3
      is=isym(i1)-1
      js=isym(j1)-1
      if (is.ne.js) go to 160
_IF(cray)
      kss=is.xor.is
_ELSE
      kss=IXOR32(is,is)
_ENDIF
      ksh=kss*norbs
_IF(cray)
      lsh=(kss.xor.is)*norbs
_ELSE
      lsh=IXOR32(kss,is)*norbs
_ENDIF
      ijkl=ijadd(i1*(i1-1)/2+i1)+kadd(i1+ksh)+ladd(j1+lsh)+m1
      itabb(ijkl)=ic1e1(i)+ic1e2(j)
160   continue
170   continue
c
c......set up mapping of two-electron integrals
      ij=0
      do 300 ia=nst,nprim
      do 300 ja=nst,ia
      ij=ij+1
      kl=0
      do 300 ka=nst,ia
      do 300 la=nst,ka
      kl=kl+1
      if(ic4e(kl).ne.ic4e(ij).or.kl.gt.ij)goto 300
      i=iout(ia)
      j=iout(ja)
      k=iout(ka)
      l=iout(la)
      if (i.lt.1.or.j.lt.1.or.k.lt.1.or.l.lt.1) go to 300
      iii=max(i,j,k,l)
      mgrp=ijgrp(iii*(iii+1)/2)
      if (mgrp.ne.ii) go to 300
      if (i.ge.j) go to 180
      m=i
      i=j
      j=m
180   if (k.ge.l) go to 190
      m=k
      k=l
      l=m
190   if (i.ge.k) go to 200
      m=i
      i=k
      k=m
      m=j
      j=l
      l=m
200   if (i.gt.k) go to 210
      if (j.ge.l) go to 210
      m=j
      j=l
      l=m
210   i1=i
      j1=j
      k1=k
      l1=l
      if (j1.ge.k1) go to 220
      m=j1
      j1=k1
      k1=m
220   if (j1.ge.l1) go to 230
      m=j1
      j1=l1
      l1=m
230   if (k1.ge.l1) go to 240
      m=k1
      k1=l1
      l1=m
240   if (i1.ne.j1.and.j1.eq.k1.and.j1.eq.l1) go to 250
      go to 260
250   j1=i1
      k1=i1
260   continue
      is=isym(i1)-1
      js=isym(j1)-1
      ks=isym(k1)-1
      ls=isym(l1)-1

_IF(cray)
      kss=is.xor.js.xor.ks.xor.ls
_ELSE
      kss=IXOR32(is,js)
      kss=IXOR32(kss,ks)
      kss=IXOR32(kss,ls)
_ENDIF
      if (kss.ne.0) go to 300
_IF(cray)
      kss=is.xor.js
_ELSE
      kss=IXOR32(is,js)
_ENDIF
      ksh=kss*norbs
_IF(cray)
      lsh=(kss.xor.ks)*norbs
_ELSE
      lsh=IXOR32(kss,ks)*norbs
_ENDIF
      m=0
      if (i.gt.j) go to 270
      m=1
      if (k.ge.l.and.j.ne.k) m=2
      go to 290
270   if (k.gt.l) go to 280
      m=2
      go to 290
280   if (j.eq.k.or.j.eq.l.or.i.eq.k) m=1
      if (m.eq.1) go to 290
      if (j.gt.k) m=2
      if (j.lt.k.and.j.gt.l) m=1
      if (j.lt.k.and.j.lt.l) m=3
290   ijkl=ijadd(icf(i1)+j1)+kadd(k1+ksh)+ladd(l1+lsh)+m
      itabb(ijkl)=ic2e(ij)+ic3e(kl)
300   continue
      write (itape8) (itabb(i),i=1,nmax)
      write (itape8) itypec
      write (itape8) ncore,nact,naa,nbb
      if((iop(1)/8)*8.eq.iop(1))return
      write(iwr,310)
310   format(' map - (guga address,casscf address)'/)
      write (iwr,320) (i,itabb(i),i=1,nmax)
320   format(10(2i5,3x))
      return
      end
      subroutine baksub(a,n,vec,d,maxcor,ifort,nprint)
c
c      routine to perform forward and backward substitution in order
c     to solve a set of linear equations, given the choleski decompositi
c     of the matrix, which may be too large to hold in core.
c
c      the ldl(transpose) decomposition of the matrix is assumed to be o
c     fortran stream ifort, in triangular form with 1 row per record.
c     the diagonal (unit) elements of l should be written over by the
c     corresponding elements of d, i.e. in the same format as produced
c     by routine cleski.  on entry, vec should hold the right hand
c     side of the equations, and on exit it contains the solution vector
c
c      workspace required:   an array a of size  maxcor which should be
c                               as large as possible.
c                            a vector d of length n*2.
c
c      external routines required:   caserr - error handling of casscf p
c                                             can be replaced by stop.
c                                    sdot      cray library
c                                    cray buffered i/o
c
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      common/junk/ibase(200),maplen(200)
      dimension a(*),d(*),vec(*)
      data mxpass/200/
c
c
      if(n.gt.maxcor)call caserr(
     + 'store available for backsubstitution less than 1 record long')
c
      np=n+1
      npass=1
      ibase(1)=n
      itop=0
10    maplen(npass)=0
20    itop=itop+1
      maplen(npass)=maplen(npass)+itop
      if(itop.gt.n)goto30
      if(maplen(npass).le.maxcor)goto20
      maplen(npass)=maplen(npass)-itop
      itop=itop-1
      npass=npass+1
      if(npass.ge.mxpass)callcaserr
     + ('too many passes needed in backsubstution routine')
      ibase(npass)=n-itop
      goto10
30    maplen(npass)=maplen(npass)-itop
      if(nprint.ne.-5)write(iwr,40)npass
40    format(' backsubstitution in',i4,' passes')
      ibase(npass+1)=0
c
c....first phase of backsubstitution - forward read.
      rewind ifort
      im=0
      do 50 i=1,n
      call rdfors(a,i,ifort)
      d(i)=a(i)
      vec(i)=vec(i)-ddot(im,a,1,vec,1)
50    im=i
c
c.... reverse order for second phase
      do 60 i=1,n
60    a(np-i)=vec(i)/d(i)
      do 70 i=1,n
70    d(i)=a(i)
c
c.....sort decomposed matrix on ifort into reverse order andidomnext pha
c.....loop over output blocks
c
      do 110 ii=1,npass
      imax=ibase(ii)
      iminm=ibase(ii+1)
      leni=imax-iminm
      imin=iminm+1
c
c.....now read through the file and select the required matrix elements
c
c.....position file at start of the current block
      irec=imin
      call rdfort(vec,0,irec,ifort)
c
c.....loop over blocks jj.ge.ii
      mapp=maplen(ii)+np
      do 90 j=imin,n
      call rdfors(vec,j,ifort)
      map=mapp-j
      itop=min(j,imax)
      do 80 i=imin,itop
      map=map+i-np
80    a(map)=vec(i)
90    continue
c
c.....buffer should now be full - process it
c
      ibot=n-ibase(ii)
      map=1
      do 100 i=1,leni
      im=ibot
      ibot=ibot+1
      d(ibot)=d(ibot)-ddot(im,a(map),1,d,1)
100   map=map+ibot
110   continue
c
c.....now revert to forward ordering
      do 120 i=1,n
120   vec(np-i)=d(i)
      rewind ifort
      return
      end
      subroutine cleski(a,n,vec,d,maxcor,iin,iout,nprint)
c
c
c      subroutine to calculated the ldl(transpose) decomposition of
c     a symmetric matrix which may be too large to hold in core.
c
c      matrix is assumed to be on fortran stream iin, in triangular
c     form, with 1 row per record - i.e. the i'th record contains the
c     i matrix elements a(i,j), j.le.i.  the matrix l, with the diagonal
c     matrix d written onto its diagonal, is output to fortran stream
c     iout in the same format.
c
c      workspace required:  an array a  of size maxcor, which should
c                               be as large as possible.
c                           2 vectors vec,d each of length n. on exit,
c                               d contains the diagonal elements of the
c                               matrix d.
c
c      external routines required:   caserr - error handling of casscf
c                                             program. can be replaced
c                                             by stop.
c                                    sdot - cray library
c                                    cray buffered i/o
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      common/junk/ ibase(200)
      dimension a(*),d(*),vec(*)
c
c
      if(n.lt.3)callcaserr
     +  ('n < 3 for choleski decomposition not allowed')
c
      if(maxcor.eq.0)call caserr('in-core not supported')
c
c....calculate how many passes are needed
c
      middle=maxcor/2
      if(n.gt.middle)callcaserr(
     1'store available for choleski decomp. less than 1 record long')
      npass=1
      ibase(1)=1
      itop=0
10    len=0
20    itop=itop+1
      len=len+itop
      if(len.le.middle)goto20
      if(itop.gt.n)goto30
      npass=npass+1
      if(npass.ge.200)callcaserr
     + ('choleski decomposition needed > 200 passes')
      ibase(npass)=itop
      itop=itop-1
      goto10
30    if(npass.eq.2)npass=1
      if(nprint.ne.-5)write(iwr,40)npass
40    format(' choleski decomposition to be performed',
     1' in',i4,' passes')
      ibase(npass+1)=n+1
c
c
c.....loop over output blocks
      rewind iin
      do 150 ii=1,npass
      imin=ibase(ii)
      iminm=imin-1
      leni=ibase(ii+1)-imin
      imax=iminm+leni
c
c
c
      rewind iout
      if(ii.eq.1)goto110
      iim=ii-1
c
c.....loop over all the jj.lt.ii blocks of the output matrix to get
c       contributions
c
      do 100 jj=1,iim
      jmin=ibase(jj)
      jminm=jmin-1
      lenj=ibase(jj+1)-jmin
      jmax=jminm+lenj
c
c.....loop over the range of i in block ii to get contributions from
c     the rows j of the processed matrix in block jj
c
      mapi=0
      i=imin-1
      do 90 ip=1,leni
      i=i+1
      if(jj.eq.1)goto 60
      do 50 j=1,jminm
50    vec(j)=a(mapi+j)*d(j)
      goto 70
60    call rdfors(a(mapi+1),i,iin)
70    jm=jminm
      mapj=1+middle
      do 80 j=jmin,jmax
      mapij=mapi+j
      if(ip.eq.1)call rdfors(a(mapj),j,iout)
      vec(j)=a(mapij)-ddot(jm,a(mapj),1,vec,1)
      jm=j
      a(mapij)=vec(j)/d(j)
80    mapj=mapj+j
90    mapi=mapi+i
100   continue
c
c.....on-diagonal block jj=ii
c
110   mapi=0
      do 150 i=imin,imax
      if(ii.eq.1)call rdfors(a(mapi+1),i,iin)
      do 120 j=1,iminm
120   vec(j)=a(mapi+j)*d(j)
      im=i-1
      if(i.eq.imin)goto140
      jm=iminm
      mapj=1
      do 130 j=imin,im
      mapij=mapi+j
      vec(j)=a(mapij)-ddot(jm,a(mapj),1,vec,1)
      jm=j
      mapj=mapj+j
130   a(mapij)=vec(j)/d(j)
c
140   mapi=mapi+i
      d(i)=a(mapi)-ddot(im,vec,1,a(mapi-im),1)
      a(mapi)=d(i)
150   call wtfors(a(mapi-im),i,iout)
      endfile iout
      return
      end
      subroutine stackp(g,nijkl,map)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(iv)      integer *2 nijkl
      dimension g(*),map(*)
      dimension nijkl(*)
INCLUDE(common/stak)
      common/junk/nwbuck(1500),
     * itx(3400),ktx(3400),gtx(3400)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
      common/bufb/nwbnwb,lnklnk,gout(1)
      common/three/mark(1500),ibase(1500),ibasen(1500)
c
      call jsortp(map)
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
    2 nstack=0
      return
      end
      subroutine scang(gmat,pmat)
c
c     subroutine to construct 2-electron contributions to fock matrix
c     using 2-electron integral files
c
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/cndx40)
INCLUDE(common/filel)
INCLUDE(common/restar)
INCLUDE(common/atmblk)
      common/blkin/g(511)
      common/craypk/ijkl(1360)
      dimension gmat(*),pmat(*)
c
c.....2-electron integrals
c
      call setsto(1360,0,ijkl)
      lfile=m1file
      do 1 i=1,lfile
      lotape(i)=m1tape(i)
      liblk(i)= m1blk(i)
 1    llblk(i)=liblk(i)-m1last(i)
c
      call dscal(lenb4,0.5d0,pmat,1)
c
_IF1(iv)      do 2 loop=1,nbas4
_IF1(iv)      ii=iky(loop+1)
_IF1(iv)    2 pmat(ii)=pmat(ii)+pmat(ii)
      do 1000 i=1,lfile
      iunit=lotape(i)
      call search(liblk(i),iunit)
      call find(iunit)
      j=llblk(i)
 40   j=j+1
      call get(g(1),nw)
      if(nw)50,1000,50
 50   if(j.ne.0)call find(iunit)
_IF(ibm,vax)
      call gmake(gmat,pmat,iky)
_ELSE
      if(o255i)then
       call sgmata(gmat,pmat)
      else
       call sgmata_255(gmat,pmat)
      endif
_ENDIF
      if(j)40,1000,40
1000  continue
c
c
_IF1(iv)      do 3 loop=1,nbas4
_IF1(iv)      ii=iky(loop+1)
_IF1(iv)    3 pmat(ii)=pmat(ii)*0.5d0
c
      call dscal(lenb4,2.0d0,pmat,1)
c
      return
      end
      subroutine ver_casb(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/casb.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
_IF(hpux11)
c$HP$ OPTIMIZE ASSUME_NO_PARAMETERS_OVERLAPS ON
_ENDIF
