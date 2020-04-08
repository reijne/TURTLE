c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util6.m,v $
c  $State: Exp $
c  
      subroutine ijconr(s,t,nx,mapnr)
      implicit REAL  (a-h,o-z)
c     condense lower triangle to n.r. members
c
      dimension s(*),t(*),mapnr(*)
c
c     nx = length of triangle
c
      do 20 i = 1 , nx
         j = mapnr(i)
         if (j.ne.0) then
c
c     s = input : t =output
c
            t(j) = s(i)
         end if
 20   continue
      return
      end
      subroutine grhfbl(scftyp)
c     set up grhf block for closed and oscf types
      implicit REAL  (a-h,o-z)
      character *8 scftyp,oscf,grhf,open
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/ghfblk)
      dimension i11(11)
      data i11/0,11,22,33,44,55,66,77,88,99,110/
      data oscf,grhf,open/'oscf','grhf','open'/
      if (scftyp.eq.grhf .or. scftyp.eq.open) return
c
c
      nact = num
      do 20 i = 1 , num
         iactiv(i) = i
 20   continue
      nbshel(1) = na
      ilfshl(1) = 0
      ilfshl(2) = na
      do 40 i = 1 , 11
         fjk(i) = 0.0d0
         fcan(i) = 0.0d0
         do 30 j = 1 , 11
            erga(i11(i)+j) = 0.0d0
            cana(i11(i)+j) = 0.0d0
            ergb(i11(i)+j) = 0.0d0
            canb(i11(i)+j) = 0.0d0
 30      continue
 40   continue
      fjk(1) = 2.0d0
      erga(i11(1)+1) = 2.0d0
      ergb(i11(1)+1) = -1.0d0
      cana(i11(1)+1) = 2.0d0
      cana(i11(2)+1) = 2.0d0
      canb(i11(1)+1) = -1.0d0
      canb(i11(2)+1) = -1.0d0
      fcan(1) = 2.0d0
      fcan(2) = 2.0d0
      if (scftyp.eq.oscf) then
         njk = 2
         njk1 = 3
         fjk(2) = 1.0d0
         fjk(3) = 0.0d0
         nbshel(2) = nb - na
         nbshel(3) = num - nb
         ilfshl(3) = nb
         erga(i11(1)+2) = 1.0d0
         erga(i11(2)+1) = 1.0d0
         erga(i11(2)+2) = 0.5d0
         erga(i11(1)+3) = 0.0d0
         erga(i11(3)+1) = 0.0d0
         erga(i11(2)+3) = 0.0d0
         erga(i11(3)+2) = 0.0d0
         erga(i11(3)+3) = 0.0d0
         ergb(i11(1)+2) = -0.5d0
         ergb(i11(2)+1) = -0.5d0
         ergb(i11(2)+2) = -0.5d0
         ergb(i11(3)+1) = 0.0d0
         ergb(i11(1)+3) = 0.0d0
         ergb(i11(2)+3) = 0.0d0
         ergb(i11(3)+2) = 0.0d0
         ergb(i11(3)+3) = 0.0d0
         return
      else
c
c     closed
c
         njk = 1
         njk1 = 2
         fjk(2) = 0.0d0
         nbshel(2) = num - na
         erga(i11(1)+2) = 0.0d0
         erga(i11(2)+1) = 0.0d0
         erga(i11(2)+2) = 0.0d0
         ergb(i11(1)+2) = 0.0d0
         ergb(i11(2)+1) = 0.0d0
         ergb(i11(2)+2) = 0.0d0
         return
      end if
      end
      subroutine dr2lag(dd,v,e,ndim)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
c     extracts the lagrangian for various types of wavefunction
c
      dimension dd(*),e(*),v(ndim,*)
c
INCLUDE(common/cigrad)
INCLUDE(common/prnprn)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/atmblk)
INCLUDE(common/infoa)
      logical exist
      character *8 open,oscf,grhf
      data zer,one,two /0.0d0,1.0d0,2.0d0/
      data open/'open'/
      data oscf /'oscf'/
      data grhf/'grhf'/
c
      if (lci .or. lmcscf .or. cigr .or. mcgr) then
         mtyp = 0
         call secloc(isecll,exist,iblok)
         if (exist) then
            call rdedx(dd,nx,iblok,ifild)
         else
            call caserr('lagrangian not found')
         end if
         return
      else if (scftyp.eq.grhf .or. scftyp.eq.oscf) then
         call secloc(isect(42),exist,iblok)
         if (exist) then
            call rdedx(dd,lds(isect(42)),iblok,ifild)
         else
            call caserr('lagrangian not found')
         end if
         return
      else
         if (mp2 .or. mp3) then
            mtyp = 0
            call secloc(isecll,exist,iblok)
            if (exist) then
               if (odebug(30). and. nprint.ne.-5)  then
                if (mp2 ) write (iwr,6010) isecll
                if (mp3 ) write (iwr,6020) isecll
               endif
               call rdedx(dd,nx,iblok,ifild)
            else
               if (mp2) call caserr('mp2 lagrangian not found')
               if (mp3) call caserr('mp3 lagrangian not found')
            end if
         else
            call vclr(dd,1,nx)
         end if
         occ = one
         if (scftyp.ne.open) occ = two
         mtyp = 0
         call secloc(isect(8),exist,iblok)
         if (exist) then
            iblvec = iblok + mvadd
            call rdedx(v,num*ncoorb,iblvec,ifild)
         else
            call caserr('vectors not found')
         end if
         mtyp = 0
         call secloc(isect(9),exist,iblok)
         if (exist) then
            call rdedx(e,lds(isect(9)),iblok,ifild)
         else
            call caserr('eigenvalues not found')
         end if
         ij = 0
         do 40 i = 1 , num
            do 30 j = 1 , i
               ij = ij + 1
               dum = zer
               do 20 k = 1 , na
                  dum = dum - e(k)*v(i,k)*v(j,k)
 20            continue
               dd(ij) = dd(ij) + dum*occ
 30         continue
 40      continue
         if (scftyp.ne.open) return
         if (nb.eq.0) return
         mtyp = 0
         call secget(isect(11),mtyp,iblok)
         iblvec = iblok + mvadd
         call rdedx(v,num*ncoorb,iblvec,ifild)
         mtyp = 0
         call secget(isect(12),mtyp,iblok)
         call rdedx(e,lds(isect(12)),iblok,ifild)
         ij = 0
         do 70 i = 1 , num
            do 60 j = 1 , i
               ij = ij + 1
               dum = zer
               do 50 k = 1 , nb
                  dum = dum - e(k)*v(i,k)*v(j,k)
 50            continue
               dd(ij) = dd(ij) + dum*occ
 60         continue
 70      continue
         return
      end if
 6010 format (/1x,'mp2 lagrangian at section ',i4)
 6020 format (/1x,'mp3 lagrangian at section ',i4)
      end
      subroutine nrmapo(mnr,nocca,noccb,nsa,iky,mapie)
      implicit REAL  (a-h,o-z)
      dimension mnr(*),mapie(*),iky(*)
c
c     nsa is total number of active m.o.'s
c     nocca is number of doubly occupied
c     noccb is total occupied
c
c     maps from lower triangle to non-redundant pairs
c
      nmax = mapie(nsa)
      lenn = iky(nmax+1)
      do 20 i = 1 , lenn
         mnr(i) = 0
 20   continue
      nsoc = noccb - nocca
      nvirta = nsa - noccb
      ntpls1 = noccb + 1
      ndpls1 = nocca + 1
      if (nocca.ne.0) then
         do 40 j = 1 , nocca
            do 30 i = ndpls1 , nsa
               it = (i-ndpls1)*nocca + j
               ij = iky(mapie(i)) + mapie(j)
               mnr(ij) = it
 30         continue
 40      continue
      end if
      if (noccb.ne.nocca) then
         do 60 j = ndpls1 , noccb
            do 50 i = ntpls1 , nsa
               it = nvirta*nocca + (i-nsoc-1)*nsoc + j - nocca
               ij = iky(mapie(i)) + mapie(j)
               mnr(ij) = it
 50         continue
 60      continue
      end if
      return
      end
      subroutine qhq1(a,q,ilifq,ncore,h,iky,nbasis)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension a(*),q(*),h(*),ilifq(*),iky(*)
      common/small/p(maxorb),y(maxorb)
c...   a=q(transpose) * h * q
c...   a and h stored in triangle form
      m1 = ilifq(1)
      nb = ilifq(2) - m1
      m1 = m1 + 1
      do 30 j = 1 , ncore
         map = ilifq(j)
         do 20 i = 1 , nbasis
            m = iky(i) + 1
            ilen = i
            p(i) = ddot(ilen,q(map+1),1,h(m),1)
            iless1 = i - 1
            call daxpy(iless1,q(map+i),h(m),1,p,1)
 20      continue
         call mxma(q(m1),nb,1,p,1,nbasis,a(iky(j)+1),1,j,j,nbasis,1)
 30   continue
      return
      end
      subroutine rdedv (v,n,nv,ibl,idev)
      implicit REAL  (a-h,o-z)
      dimension v(n,nv)
      call search(ibl,idev)
      call find(idev)
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      call get(v(i,iv),nw)
         i = i + nw
         k = k - nw
         if (k.gt.0 .or. iv.lt.nv) call find(idev)
         if (k.lt.0) go to 40
         if (k.ne.0) go to 20
 30   continue
      return
 40   call caserr('error in rdedv')
      end
      subroutine rdedvs(v,n,nv,idev)
      implicit REAL  (a-h,o-z)
      dimension v(n,nv)
      call find(idev)
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      call get(v(i,iv),nw)
         i = i + nw
         k = k - nw
         if (k.gt.0 .or. iv.lt.nv) call find(idev)
         if (k.lt.0) go to 40
         if (k.ne.0) go to 20
 30   continue
      return
 40   call caserr('error in rdedvs')
      end
      subroutine squars(t,sq,n)
      implicit REAL  (a-h,o-z)
c     triangle to square
      dimension sq(*),t(*)
      ij = 0
      ii = 0
      do 30 i = 1 , n
         jj = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 20 j = 1 , i
            ij = ij + 1
            sq(ii+j) = -t(ij)
            sq(jj+i) = t(ij)
            jj = jj + n
 20      continue
         ii = ii + n
 30   continue
      return
      end
      subroutine stvin2
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/junk/xyzsp(3,maxat),
     * xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h(28)
      common/wermit/w(28)
      dimension min(7),max(7)
      data min /1,2,4,7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero,one /0.0d0,1.0d0/
      xint = zero
      yint = zero
      zint = zero
      npts = (ni+nj)/2
      imin = min(npts)
      imax = max(npts)
      do 160 i = imin , imax
         px = one
         py = one
         pz = one
         dum = h(i)/t
         ptx = dum + x0
         pty = dum + y0
         ptz = dum + z0
         ax = ptx - xi
         ay = pty - yi
         az = ptz - zi
         bx = ptx - xj
         by = pty - yj
         bz = ptz - zj
         go to (80,70,60,50,40,30,20) , ni
 20      px = ax
         py = ay
         pz = az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      px = px*ax
         py = py*ay
         pz = pz*az
 70      px = px*ax
         py = py*ay
         pz = pz*az
 80      go to (150,140,130,120,110,100,90) , nj
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     px = px*bx
         py = py*by
         pz = pz*bz
 140     px = px*bx
         py = py*by
         pz = pz*bz
 150     dum = w(i)
         xint = xint + dum*px
         yint = yint + dum*py
         zint = zint + dum*pz
 160  continue
      return
      end
      subroutine clenms(messge)
c
c
      implicit REAL  (a-h,o-z)
      character messge*(*)
c
INCLUDE(common/iofile)
c
c
      l = len(messge)
      write (iwr,'(/1x,a)') messge(1:l)
      call clenup
      stop
      end
      subroutine wrt3z(iblk,num3,len)
      implicit REAL  (a-h,o-z)
      common/blkin/q(511)
      call search(iblk,num3)
      do 20 i = 1 , len
         call vclr(q,1,511)
         call put(q,511,num3)
 20   continue
      call search(iblk,num3)
      return
      end
      subroutine vwrt3s(v,n,nv,idev)
      implicit REAL  (a-h,o-z)
      dimension v(n,nv)
      data ilen/511/
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      nw = min(k,ilen)
         call put(v(i,iv),nw,idev)
         i = i + nw
         k = k - nw
         if (k.gt.0) go to 20
 30   continue
      return
      end
_IFN(cray,convex)
      subroutine daxpyi(n,sa,sx,indx,sy)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension sx(*),sy(*),indx(*)
      if(n.le.0.or.sa.eq.0.d0) return
      m = n - (n/4)*4
      if( m .eq. 0 ) go to  4
      do  3 i = 1,m
        sy(indx(i)) = sy(indx(i)) + sa*sx(i)
    3 continue
      if( n .lt. 4 ) return
    4 mp1 = m + 1
      do  8 i = mp1,n,4
        sy(indx(i)) = sy(indx(i)) + sa*sx(i)
        sy(indx(i+1)) = sy(indx(i+1)) + sa*sx(i + 1)
        sy(indx(i+2)) = sy(indx(i+2)) + sa*sx(i + 2)
        sy(indx(i+3)) = sy(indx(i+3)) + sa*sx(i + 3)
    8 continue
      return
      end
      subroutine dlstgt(n,x,incx,scalar,nindx,indx)
      implicit REAL  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).le.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstne(n,x,incx,scalar,nindx,indx)
      implicit REAL  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).eq.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstge(n,x,incx,scalar,nindx,indx)
      implicit REAL  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).lt.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstlt(n,x,incx,scalar,nindx,indx)
      implicit REAL  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).ge.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
_ENDIF
      subroutine ver_util6(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util6.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
