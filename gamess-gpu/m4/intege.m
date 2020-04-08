c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/intege.m,v $
c  $State: Exp $
c  
      subroutine cangt1(buf,ibuf,lbuf,inbrel,map,nzero,q,ifort)
      implicit REAL  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
c
c     reads in canonical list, version 1 - produces
c     map and condensed triangle
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real read from buf/ibuf
c     inbint = position of last integer read from buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      if (inbrel.eq.-1) then
         call reads(buf,lbuf,ifort)
         inbrel = 0
      end if
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     get value of nzero in buffer
c
 20   nzero = ibuf(lenrel(inbrel)+1)
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be read
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            map(nfirst+i) = ibuf(inbint+i)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            map(nfirst+i) = ibuf(inbint+i)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            q(nfirst+i) = buf(inbrel+i)
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            q(nfirst+i) = buf(inbrel+i)
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call reads(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60) , iret
c
 100  return
      end
      subroutine cangt2(buf,ibuf,lbuf,inbrel,map,nzero,q,ifort)
      implicit REAL  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
c
c     reads in canonical list
c     version 2 - produces map and expanded list of triangle
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real read from buf/ibuf
c     inbint = position of last integer read from buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      if (inbrel.eq.-1) then
         call reads(buf,lbuf,ifort)
         inbrel = 0
      end if
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     get value of nzero in buffer
c
 20   nzero = ibuf(lenrel(inbrel)+1)
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be read
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            map(nfirst+i) = ibuf(inbint+i)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            map(nfirst+i) = ibuf(inbint+i)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            q(map(nfirst+i)) = buf(inbrel+i)
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            q(map(nfirst+i)) = buf(inbrel+i)
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call reads(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60) , iret
c
 100  return
      end
      subroutine canonc(q,maxq,num,iblki,ifili,ifort,ltin)
c
c-------------------------------------------------------------
c     put list of integrals into canonical order
c-------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      logical ltin
      dimension q(maxq)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
INCLUDE(common/atmblk)
INCLUDE(common/blksiz)
      common/craypk/labs(1360)
c
c     ibl5 = no of integrals in block of sortfile
c     only machine dependent feature should be structure of /bufa/
c
      iwrt = 6
      ibl5 = nsz340
      nav = lenwrd()
      iilen = nsz340*nav/2
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = num*(num+1)/2
      n2 = nij
c
      maxt = (maxq-n2-lenint(n2))/nij
      nword = (maxq/(1+lenrel(1)))*lenrel(1)
c
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
c
c     nbuck is the number of buckets required
c
      nbuck = (nij/maxt) + 1
      nadd = min(maxt,nij)
      maxa = nbuck*(ibl5+ibl5i)
c
c
      if (nbuck.gt.maxb) then
         write (iwrt,6010) maxq , maxa
         call caserr('stop')
      end if
c
c     read through original file producing sorted file
c
      call vclr(q,1,maxa)
      call setsto(numlab,0,labs)
      call setbfa
      call canrd1(q,q,maxa,iblki,ifili,ltin)
c
c     read through the sort file to give final result
c
      maxqq = nij*nadd
      lbuf = nsz*512
      i1 = 1
      i2 = i1 + lbuf
      i3 = i2 + lenint(n2)
      call canwt1(q(i1),lbuf,q(i2),q(i3),maxqq,ifort)
c
      call closbf
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine canput(buf,ibuf,lbuf,inbrel,map,nzero,q,last,ifort)
      implicit REAL  (a-h,o-z)
      dimension buf(lbuf),ibuf(*),map(*),q(*)
      logical last
c
c     writes out canonical list
c
c     note : buf and ibuf ( = output buffer ) are
c     physically the same array.
c     inbrel = position of last real written into buf/ibuf
c     inbint = position of last integer written into buf/ibuf
c     map = map of the nzero elements of triangle q
c     nzero = number of nonzero elements
c     q = triangle of integrals
c
      left = lbuf - inbrel
c
c     there are left real words remaining in buffer
c
      iret = 1
      if (left.eq.0) go to 90
c
c     put value of nzero in buffer
c
 20   ibuf(lenrel(inbrel)+1) = nzero
      inbrel = inbrel + 1
      left = lbuf - inbrel
      iret = 2
      nfirst = 0
      if (left.eq.0) go to 90
 30   inbint = lenrel(inbrel)
      nspace = lenrel(lbuf) - inbint
      need = nzero - nfirst
c
c     there are need labels still to be written
c     and nspace integer words available in the buffer
c
      if (need.le.nspace) then
         do 40 i = 1 , need
            ibuf(inbint+i) = map(i+nfirst)
 40      continue
         inbrel = inbrel + lenint(need)
         left = lbuf - inbrel
c
c     now for real part
c
         iret = 3
         nfirst = 0
         if (left.eq.0) go to 90
      else
         do 50 i = 1 , nspace
            ibuf(inbint+i) = map(i+nfirst)
 50      continue
         nfirst = nfirst + nspace
         go to 90
      end if
 60   nspace = lbuf - inbrel
      need = nzero - nfirst
      if (need.le.nspace) then
         do 70 i = 1 , need
            buf(inbrel+i) = q(map(i+nfirst))
 70      continue
         inbrel = inbrel + need
         go to 100
      else
         do 80 i = 1 , nspace
            buf(inbrel+i) = q(map(i+nfirst))
 80      continue
         nfirst = nfirst + nspace
      end if
c
c
 90   call wrt3s(buf,lbuf,ifort)
      inbrel = 0
      go to (20,30,60,110) , iret
c
 100  if (last) then
         iret = 4
         go to 90
      end if
 110  return
      end
      subroutine canrd1(a,ia,maxa,iblki,ifili,ltin)
c
c     does the sorting part to get coulomb matrices
c     lower triangles only are produced
c     adaption for use in sorting to canonical order ,
c-----------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      logical ltin
INCLUDE(common/sizes)
c
c     note arrays a and ia actually overlap
c
      dimension a(maxa),ia(maxa)
INCLUDE(common/three)
INCLUDE(common/atmblk)
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
      common/blkin/gin(510),nint
      common/stak/btri,mlow,nstack,iblock
INCLUDE(common/blksiz)
      common/bufb/nkk1,mkk1,g(1)
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
INCLUDE(common/mapper)
      dimension iaad(340),ibbu(340)
      data lastb/999999/
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i)
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      call vclr(g,1,nsz340+nsz170)
c
      iblock = 0
c     ninb = no of elements in bucket(coreload)
      ninb = nadd*nij
c
      call search(iblki,ifili)
 30   call find(ifili)
      call get(gin,nw)
      if (nw.eq.0) then
c
c     empty anything remaining in buckets
c
         do 40 ibuck = 1 , nbuck
            nwb = nwbuck(ibuck)
            if (nwb.ne.0) then
               call stopbk
               mkk1 = mark(ibuck)
               nkk1 = nwb
_IFN1(iv)               call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)               call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)               call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)               call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 40      continue
c
c
         call stopbk
         return
      else
         if (ltin) then
_IFN1(iv)            call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)            call upak8v(gin(num2e+1),i205)
            do 50 int = 1 , nint
_IFN1(iv)               n4 = int + int + int + int
_IF(ibm,vax)
               ij = iky(i205(int)) + j205(int)
               kl = iky(k205(int)) + l205(int)
_ELSEIF(littleendian)
               ij = iky(labs(n4-2)) + labs(n4-3)
               kl = iky(labs(n4  )) + labs(n4-1)
_ELSE
               ij = iky(labs(n4-3)) + labs(n4-2)
               kl = iky(labs(n4-1)) + labs(n4)
_ENDIF
c
               iaddr = (ij-1)*nij + kl
c
c     iaddr is address of integral in final sequence
c
               ibuck = (iaddr-1)/ninb
               iaad(int) = iaddr - ninb*ibuck
               ibbu(int) = ibuck + 1
 50         continue
         else
_IFN1(iv)            call unpack(gin(num2e+1),16,labs,680)
_IF1(iv)            call upak4v(igin,i205)
            do 60 int = 1 , nint
_IFN1(iv)               n4 = int + int
c
_IFN1(iv)               iaddr = (labs(n4-1)-1)*nij + labs(n4)
_IF1(iv)               iaddr = (i205(int)-1)*nij + j205(int)
c
c     iaddr is address of integral in final sequence
c
               ibuck = (iaddr-1)/ninb
               iaad(int) = iaddr - ninb*ibuck
               ibbu(int) = ibuck + 1
 60         continue
         end if
c
c     element goes in bucket ibuck with modified address
c
         do 70 int = 1 , nint
            ibuck = ibbu(int)
            nwb = nwbuck(ibuck) + 1
            a(ibase(ibuck)+nwb) = gin(int)
            ia(ibasen(ibuck)+nwb) = iaad(int)
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - empty
c
               call stopbk
               mkk1 = mark(ibuck)
               nkk1 = nwb
_IFN1(iv)               call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)               call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)               call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)               call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
c
 70      continue
         go to 30
      end if
      end
      subroutine canwt1(buf,lbuf,ibuff,q,maxqq,ifort)
c----------------------------------------------------------
c     this reads back down the back-chained sort file
c     produced by canrd1 to give a final file
c     containing the coulomb matrices arranged
c     sequentially
c----------------------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical last
      dimension q(maxqq),buf(lbuf),ibuff(n2)
      common/sortpk/labin(1)
      common/stak/btri,mlow,nstack,iblock
INCLUDE(common/blksiz)
      common/bufb/nkk,mkk,g(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
INCLUDE(common/atmblk)
INCLUDE(common/mapper)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
_IF1(iv)      dimension lin(2)
_IF1(iv)      equivalence (lin(1),g(1))
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call rewedz(ifort)
c
      min = 1
      max = nadd
      inbrel = 0
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxqq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     triangles min thru max are in core - clear them out
c
            itri = 1
            do 30 n = min , max
c
c       get map of non-zero elements
c
_IFN1(c)               call dlstne(n,q(itri),1,0.0d0,nzero,ibuff)
_IF1(c)               call whenne(n,q(itri),1,0.0d0,ibuff,nzero)
               last = n.eq.nij
               call canput(buf,buf,lbuf,inbrel,ibuff,nzero,q(itri),last
     +  ,ifort)
               itri = itri + nij
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nij) max = nij
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
_IFN1(iv)            call unpack(g(nsz341),32,labin,ibl5)
_IFN1(civ)            call dsctr(nkk,g,labin,q)
_IF1(c)            call scatter(nkk,q,labin,g)
_IF1(iv)            ij = ibl5+ibl5+1
_IF1(iv)            do 4000 iword=1,nkk
_IF1(iv)            q(lin(ij)) = g(iword)
_IF1(iv) 4000       ij = ij+1
c
            go to 20
         end if
 40   continue
      return
      end
      subroutine ver_intege(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/intege.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
