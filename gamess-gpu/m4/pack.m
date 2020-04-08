      subroutine pack(ip,nbits,iu,nw)
      implicit REAL  (a-h,o-z)
      REAL  ip
c
c...   takes  8 16 or 32 rightmost bits of u and packs them
c...   closely in p / nbits = # bits / nw = # u -elements
c
c
      dimension ip(*),iu(nw)
c
c***   not the right number of words does occurr !!!
c
      if (nbits.eq.8) then
         call vipk8(iu,ip,nw)
      else if (nbits.eq.16) then
         call vipk16(iu,ip,nw)
      else if (nbits.eq.32) then
         call vipk32(iu,ip,nw)
      else
         call caserr(' *** pack called with nbits ne 8,16,32 ***')
      end if
c
      return
      end
      subroutine unpack(ip,nbits,iu,nw)
      implicit REAL  (a-h,o-z)
      REAL  ip
c
c...   reverses action of pack / see there
c
c...   *** 8 / 16 / 32 are probably only ones for gamess ***
c
      dimension ip(*),iu(nw)
c
c        ((nw/nbits)*nbits.ne.nw)  does happen
c
      if (nbits.eq.8) then
c *** this redundant array zeroing to cure the problem in
c *** mrdci/skiny2 ... must be resolved
         call setsto(nw,0,iu)
c ***  above line to be removed once problem resolved 
         call viup8(ip,iu,nw)
      else if (nbits.eq.16) then
         call viup16(ip,iu,nw)
      else if (nbits.eq.32) then
         call viup32(ip,iu,nw)
      else
         call caserr(' *** pack called with nbits ne 8,16,32 ***')
      end if
c
      return
      end
      subroutine vipk8(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
      logical *1 iout4,in4
_IFN1(k)      dimension in4(4,*),iout4(*)
_IF1(k)      dimension in4(8,*),iout4(*)
c
      ind = 1
      do 1 loop=1,nw/4
_IF(ksr)
      iout4(ind  )=in4(8,ind  )
      iout4(ind+1)=in4(8,ind+1)
      iout4(ind+2)=in4(8,ind+2)
      iout4(ind+3)=in4(8,ind+3)
_ELSEIF(littleendian)
      iout4(ind  )=in4(1,ind  )
      iout4(ind+1)=in4(1,ind+1)
      iout4(ind+2)=in4(1,ind+2)
      iout4(ind+3)=in4(1,ind+3)
_ELSE
      iout4(ind  )=in4(4,ind  )
      iout4(ind+1)=in4(4,ind+1)
      iout4(ind+2)=in4(4,ind+2)
      iout4(ind+3)=in4(4,ind+3)
_ENDIF
 1    ind =ind+4
      return
      end
      subroutine vipk16(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
      integer *2 iout4
      dimension in4(*),iout4(*)
c
      ind = 1
      do 1 loop=1,nw/4
      iout4(ind  )=in4(ind  )
      iout4(ind+1)=in4(ind+1)
      iout4(ind+2)=in4(ind+2)
      iout4(ind+3)=in4(ind+3)
 1    ind =ind +4
      return
      end
      subroutine vipk32(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
_IF1(k)      integer *4 iout4
      dimension in4(*),iout4(*)
c
      ind = 1
      do 1 loop=1,nw/4
      iout4(ind  )=in4(ind  )
      iout4(ind+1)=in4(ind+1)
      iout4(ind+2)=in4(ind+2)
      iout4(ind+3)=in4(ind+3)
 1    ind =ind +4
      return
      end
      subroutine viup16(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
      integer *2 in4
      dimension iout4(*),in4(*)
c
      ind = 1
      do 1 loop=1,nw/4
      iout4(ind  )=in4(ind  )
      iout4(ind+1)=in4(ind+1)
      iout4(ind+2)=in4(ind+2)
      iout4(ind+3)=in4(ind+3)
 1    ind =ind +4
      return
      end
      subroutine viup32(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
_IF1(k)      integer *4 in4
      dimension in4(*),iout4(*)
c
      ind = 1
      do 1 loop=1,nw/4
      iout4(ind  )=in4(ind  )
      iout4(ind+1)=in4(ind+1)
      iout4(ind+2)=in4(ind+2)
      iout4(ind+3)=in4(ind+3)
 1    ind =ind +4
      return
      end
      subroutine viup8(in4,iout4,nw)
      implicit REAL  (a-h,o-z)
      logical *1 in4,iout4
_IFN1(k)      dimension in4(*),iout4(4,*)
_IF1(k)      dimension in4(*),iout4(8,*)
c
      ind = 1
      do 1 loop=1,nw/4
_IF(ksr)
      iout4(8,ind  )=in4(ind  )
      iout4(8,ind+1)=in4(ind+1)
      iout4(8,ind+2)=in4(ind+2)
      iout4(8,ind+3)=in4(ind+3)
_ELSEIF(littleendian)
      iout4(1,ind  )=in4(ind  )
      iout4(1,ind+1)=in4(ind+1)
      iout4(1,ind+2)=in4(ind+2)
      iout4(1,ind+3)=in4(ind+3)
_ELSE
      iout4(4,ind  )=in4(ind  )
      iout4(4,ind+1)=in4(ind+1)
      iout4(4,ind+2)=in4(ind+2)
      iout4(4,ind+3)=in4(ind+3)
_ENDIF
 1    ind =ind +4
      return
      end
c>>>>>>>>>>>c
c this stuff from util3
      function pack2(nw,m)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(convex)
      integer*8 pack,pack2
      dimension  i4(2)
      equivalence (i4(1),pack)
      i4(1)=nw
      i4(2)=m
      pack2 = pack
_ELSEIF(t3d)
      real pack2
      pack2=shiftl(nw,32).or.m
_ELSEIF(cray)
      integer pack2
      pack2=shiftl(nw,32).or.m
_ELSEIF(ksr)
      integer pack2,pack
      integer *4 i4(2)
      equivalence (i4(1),pack)
      i4(1)=nw
      i4(2)=m
      pack2 = pack
_ELSE
c
c these version for 32-bit integer machines
c
      dimension  i4(2)
      equivalence (i4(1),pack)
_IF(littleendian)
      i4(2)=nw
      i4(1)=m
_ELSE
      i4(1)=nw
      i4(2)=m
_ENDIF
      pack2 = pack
_ENDIF
_IF1()       pack2=or(shift(nw,32),m)
      return
      end
      subroutine upack2(b,left,iright)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(convex)
      dimension i4(2)
      equivalence (i4(1),unpack)
      unpack=b
      left=i4(1)
      iright=i4(2)
_ELSEIF(cray,t3d)
      left=shiftr(b,32)
      iright=b.and.37777777777b
_ELSEIF(ksr)
      integer *4 i4(2)
      equivalence (i4(1),unpack)
      unpack=b
      left=i4(1)
      iright=i4(2)
_ELSE
      dimension i4(2)
      equivalence (i4(1),unpack)
      unpack=b
_IF(littleendian)
      left=i4(2)
      iright=i4(1)
_ELSE
      left=i4(1)
      iright=i4(2)
_ENDIF
_ENDIF
_IF1()       data mask/x'ffffffff'/
_IF1()       left=and(shift(b,32),mask)
_IF1()       iright=and(b,mask)
      return
      end
      subroutine upack3(packed,i,j,k)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(t3d,cray)
      i=shiftr(packed,32)
      j=shiftr(packed,16).and.177777b
      k=packed.and.177777b
_ELSE
      integer * 2 i2(4)
      equivalence (i2(1),pack)
      pack=packed
      i = i2(2)
      j = i2(3)
      k = i2(4)
_ENDIF
      return
      end
      function pack3(i,j,k)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(t3d,cray)
      pack3=shiftl(i,32).or.shiftl(j,16).or.k
_ELSE
      integer * 2 i2(4)
      equivalence (i2(1),pack)
      i2(1) = 0
      i2(2) = i
      i2(3) = j
      i2(4) = k
      pack3=pack
_ENDIF
      return
      end
c
c
_IF1(ivf)       subroutine upak8v(ix,intij)
_IF1(ivf)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)      implicit character *8 (z),character *1 (x)
_IF1(ivf)      implicit character *4 (y)
_IF1(i)      logical * 1 mij,mkl
_IF1(i)      logical*  1 i205,j205,k205,l205
_IF1(i)      dimension ix(*),intij(*)
_IF1(i)      common/blkin/gout(340),mij(680),mkl(680),nint
_IF1(i)      common/craypk/i205(1360),j205(1360),k205(1360),l205(1360)
_IF1(i)      i4=4
_IF1(i)      do 1 i=1,nint
_IF1(i)      il=i+i
_IF1(i)      i205(i4)=mij(il-1)
_IF1(i)      j205(i4)=mij(il  )
_IF1(i)      k205(i4)=mkl(il-1)
_IF1(i)      l205(i4)=mkl(il  )
_IF1(i)    1 i4=i4+4
_IF1(v)      logical * 1 mij,mkl
_IF1(v)      dimension ix(*),intij(*)
_IF1(v)      common/blkin/gout(340),mij(680),mkl(680),nint
_IF1(v)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IF1(v)      do 1 i=1,nint
_IF1(v)      il=i+i
_IF1(v)      j205(i)=mij(il-1)
_IF1(v)      i205(i)=mij(il  )
_IF1(v)      l205(i)=mkl(il-1)
_IF1(v) 1    k205(i)=mkl(il  )
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      ii=1
_IF1(f)      do 1 i=1,85
_IF1(f)      intij(ii  )=extract(ix(i),0,8)
_IF1(f)      intij(ii+1)=extract(ix(i),16,8)
_IF1(f)      intij(ii+2)=extract(ix(i),32,8)
_IF1(f)      intij(ii+3)=extract(ix(i),48,8)
_IF1(f)      jj=ii+340
_IF1(f)      intij(jj  )=extract(ix(i),8,8)
_IF1(f)      intij(jj+1)=extract(ix(i),24,8)
_IF1(f)      intij(jj+2)=extract(ix(i),40,8)
_IF1(f)      intij(jj+3)=extract(ix(i),56,8)
_IF1(f)      kk=jj+340
_IF1(f)      j=i+85
_IF1(f)      intij(kk  )=extract(ix(j),0,8)
_IF1(f)      intij(kk+1)=extract(ix(j),16,8)
_IF1(f)      intij(kk+2)=extract(ix(j),32,8)
_IF1(f)      intij(kk+3)=extract(ix(j),48,8)
_IF1(f)      ll=kk+340
_IF1(f)      intij(ll  )=extract(ix(j),8,8)
_IF1(f)      intij(ll+1)=extract(ix(j),24,8)
_IF1(f)      intij(ll+2)=extract(ix(j),40,8)
_IF1(f)      intij(ll+3)=extract(ix(j),56,8)
_IF1(f)   1  ii=ii+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(ivf)       subroutine pak4v(intij,ix)
_IF1(ivf)       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)       implicit character *8 (z),character *1 (x)
_IF1(ivf)       implicit character *4 (y)
_IF1(iv)      integer*2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      j=341
_IF1(iv)      do 1 i=1,340
_IF1(iv)      ix(i  )= intij(i)
_IF1(iv)      ix(j  )= intij(j)
_IF1(iv)    1 j=j+1
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      int=1
_IF1(f)      do 1 i=1,170
_IF1(f)      ix(i)=insert(
_IF1(f)     *             insert(intij(int+3),32,16,intij(int+2)),0,32,
_IF1(f)     *             insert(intij(int+1),32,16,intij(int)) )
_IF1(f)   1  int=int+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(iv)      subroutine upak4v(ix,intij)
_IF1(iv)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(iv)      implicit character *8 (z),character *1 (x)
_IF1(iv)      implicit character *4 (y)
_IF1(iv)      integer * 2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      j=341
_IF1(iv)      do 1 i=1,340
_IF1(iv)      intij(i)=ix(i)
_IF1(iv)      intij(j)=ix(j)
_IF1(iv)    1 j=j+1
_IF1(iv)      return
_IF1(iv)      end
_IF1(ivf)       subroutine pak6v(intij,ix)
_IF1(ivf)       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)       implicit character *8 (z),character *1 (x)
_IF1(ivf)       implicit character *4 (y)
_IF1(iv)      integer*2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      common/blkin/gout(510),mword
_IF1(iv)      j=205
_IF1(iv)      do 1 i=1,mword
_IF1(iv)      ix(i  )= intij(i)
_IF1(iv)      ix(j  )= intij(j)
_IF1(iv)    1 j=j+1
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      int=1
_IF1(f)      do 1 i=1,102
_IF1(f)      ix(i)=insert(
_IF1(f)     *             insert(intij(int+3),32,16,intij(int+2)),0,32,
_IF1(f)     *             insert(intij(int+1),32,16,intij(int)) )
_IF1(f)    1 int=int+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(iv)      subroutine upak6v(ix,intij)
_IF1(iv)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(iv)      implicit character *8 (z),character *1 (x)
_IF1(iv)      implicit character *4 (y)
_IF1(iv)      integer * 2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      common/blkin/gout(510),mword
_IF1(iv)      j=205
_IF1(iv)      do 1 i=1,mword
_IF1(iv)      intij(i)=ix(i)
_IF1(iv)      intij(j)=ix(j)
_IF1(iv)    1 j=j+1
_IF1(iv)      return
_IF1(iv)      end

_IF(convex,titan)
      subroutine upak8z(nword,iii,i,j,k,l)
      implicit REAL  (a-h,o-z)
      logical *1 iii,i,j,k,l
      dimension iii(*),i(4,*),j(4,*),k(4,*),l(4,*)
      kk=1
      do 1 loop=1,nword
      i(4,loop)=iii(kk  )
      j(4,loop)=iii(kk+1)
      k(4,loop)=iii(kk+2)
      l(4,loop)=iii(kk+3)
  1   kk=kk+4
      return
      end
_ENDIF

c
c   isoin/isoout
c
      subroutine isoin(nt)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(cray,t3d)
      common/isopac/ indin(48),indout(12)
      dimension ibit(4)
      data ibit /0,16,32,48/
      iout=0
      j=0
      i=0
      k=1
    1 i=i+1
      if(i.gt.nt)go to 3
      ii=indin(i)
      j=j+1
      if(j.lt.5)go to 2
      indout(k)=iout
      k=k+1
      iout=0
      j=1
    2 its=shiftl(ii,ibit(j))
      iout=iout.or.its
      go to 1
    3 indout(k)=iout
_ELSEIF(ksr)
      integer *2 log,kog
      common/isopac/in(48),it(12)
      dimension log(4),kog(4)
      equivalence (log(1),ind),(kog(1),jnd)
c
      imax=0
      i=0
    1 imin=imax+1
      imax=imax+4
      if(imax.gt.nt)imax=nt
      do 2 ii=1,4
      n=imin-1+ii
      if(n.gt.nt)go to 3
      jnd=in(n)
    2 log(ii)=kog(4)
    3 i=i+1
      it(i)=ind
      if(imax.lt.nt)go to 1
_ELSE
      logical *1 log,kog
      common/isopac/in(48),it(12)
      dimension log(4),kog(4)
      equivalence (log(1),ind),(kog(1),jnd)
      imax=0
      i=0
    1 imin=imax+1
      imax=imax+4
      if(imax.gt.nt)imax=nt
      do 2 ii=1,4
      n=imin-1+ii
      if(n.gt.nt)go to 3
      jnd=in(n)
_IF(littleendian)
    2 log(ii)=kog(1)
_ELSE
    2 log(ii)=kog(4)
_ENDIF
    3 i=i+1
      it(i)=ind
      if(imax.lt.nt)go to 1
_ENDIF
      return
      end
      subroutine isoout(nt)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(cray,t3d)
      common/isopac/indin(48),indout(12)
      dimension mask(4),ibit(4)
      data ibit /0,16,32,48/
      data mask /
     *                177777b,
     *           37777600000b,
     *      7777740000000000b,
     *1777770000000000000000b /
c
      ij=0
      iunp=4
      ntwd=(nt+3)/4
      ilast=nt - (ntwd-1)* 4
      ii=0
    1 ii=ii+1
      if(ii-ntwd)2,3,4
    3 iunp=ilast
    2 do 5 i=1,iunp
      iq=indout(ii).and.mask(i)
      ij=ij+1
    5 indin(ij)=shiftr(iq,ibit(i))
      go to 1
    4 continue
_ELSEIF(ksr)
      integer *2 log,kog
      common/isopac/in(48),it(12)
      dimension log(4),kog(4)
      equivalence (log(1),ind),(kog(1),jnd)
      data jnd/0/
c
      do 5 i=1,nt
    5 in(i)=0
c
      i=0
      imax=0
   10 imin=imax+1
      imax=imax+4
      if(imax.gt.nt)imax=nt
      i=i+1
      ind=it(i)
      do 20 ii=1,4
      n=imin-1+ii
      if(n.gt.nt)go to 30
      kog(4)=log(ii)
   20 in(n)=jnd
   30 if(imax.lt.nt)go to 10
c
_ELSE
      logical *1 log,kog
      common/isopac/in(48),it(12)
      dimension log(4),kog(4)
      equivalence (log(1),ind),(kog(1),jnd)

      data jnd/0/
c
      do 5 i=1,nt
    5 in(i)=0
c
      i=0
      imax=0
   10 imin=imax+1
      imax=imax+4
      if(imax.gt.nt)imax=nt
      i=i+1
      ind=it(i)
      do 20 ii=1,4
      n=imin-1+ii
      if(n.gt.nt)go to 30
_IF(littleendian)
      kog(1)=log(ii)
_ELSE
      kog(4)=log(ii)
_ENDIF
   20 in(n)=jnd
   30 if(imax.lt.nt)go to 10
_ENDIF
      return
      end
      subroutine ver_pack(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/pack.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
