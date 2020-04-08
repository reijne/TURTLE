c 
c  $Author: jmht $
c  $Date: 2006-01-13 17:54:55 $
c  $Locker:  $
c  $Revision: 1.15 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/nag.m,v $
c  $State: Exp $
c  
      subroutine f04atf(a,ia,b,n,c,aa,iaa,wks1,ipvt,info)
      implicit real*8 (a-h,o-z)
      dimension a(ia,n),b(n),c(n),aa(iaa,n),wks1(n),ipvt(n)
c
c     this is not the nag routine f04atf, though the
c     arguments are the same
c     it is a (hopefully) less machine dependent replacement
c     based on vector algorithm of dongarra using
c     lu decomposition
c     see acm transactions of mathematical software
c     vol 10 sept 1984
c
      do 10 i=1,n
      do 10 j=1,n
10    aa(j,i)=a(j,i)
      call lu(aa,iaa,n,ipvt,info)
       if(info.eq.0) then
       call lus(aa,iaa,n,ipvt,c,b)
       endif
      return
      end
      subroutine lus(a,lda,n,ipvt,x,b)
      implicit real*8 (a-h,o-z)
      dimension a(lda,*),x(*),b(*),ipvt(*)
c
c     simultaneous equation ax=b given prior lu
c     decomposition of a
c     taken from dongarra acm trans math software vol 10
c
      do 10 k=1,n
10    x(k)=b(k)
      do 20 k=1,n
      l=ipvt(k)
      xk=x(l)
      x(l)=x(k)
20    x(k)=xk
      do 40 k=1,n
      xk=x(k)*a(k,k)
      if(xk.ne.0.0d0)then
      do 30 i=k+1,n
30    x(i)=x(i)-a(i,k)*xk
      endif
40    x(k)=xk
      do 60 k=n,1,-1
      xk=x(k)
      if(xk.ne.0.0d0)then
      do 50 i=1,k-1
50    x(i)=x(i)+a(i,k)*xk
      endif
60    continue
      return
      end
      subroutine lu(a,lda,n,ipvt,info)
      implicit real*8 (a-h,o-z)
      dimension a(lda,*),ipvt(*)
c
c     lu decomposition of a
c
      info=0
      do 40 j=1,n
      call smxpy(n-j+1,a(j,j),j-1,lda,a(1,j),a(j,1))
c
c     find pivot
c     note : this is isamax operation
c     (use cray routine or vector merge)
c
      t=dabs(a(j,j))
      k=j
      do 10 i=j+1,n
      if(dabs(a(i,j)).gt.t)then
      t=dabs(a(i,j))
      k=i
      endif
10    continue
      ipvt(j)=k
c
      if(t.eq.0.0d0)then
      write(6,*) 'zero pivot : matrix singular'
      info = 1
      return
      endif
c
c     swop rows
c     note : this is sswap operation on cray
c
      do 20 i=1,n
      t=a(j,i)
      a(j,i)=a(k,i)
      a(k,i)=t
20    continue
c
      a(j,j)=1.0d0/a(j,j)
      call sxmpy(n-j,lda,a(j,j+1),j-1,lda,a(j,1),lda,a(1,j+1))
      t=-a(j,j)
      do 30 i=j+1,n
30    a(j,i)=t*a(j,i)
40    continue
      return
      end
      subroutine smxpy(n1,y,n2,lda,x,a)
      implicit real*8 (a-h,o-z)
      dimension y(*),x(*),a(lda,*)
c
c     y = y +a x
c
      j=mod(n2,2)
      if(j.ge.1)then
      do 10 i=1,n1
10    y(i)=(y(i))+x(j)*a(i,j)
      endif
      j=mod(n2,4)
      if(j.ge.2)then
      do 20 i=1,n1
20    y(i)=((y(i))+x(j-1)*a(i,j-1))+x(j)*a(i,j)
      endif
      jmin=j+4
      do 40 j=jmin,n2,4
      do 30 i=1,n1
      y(i)=((((y(i))+x(j-3)*a(i,j-3))+x(j-2)*a(i,j-2))
     &    +x(j-1)*a(i,j-1))+x(j)*a(i,j)
30    continue
40    continue
      return
      end
      subroutine sxmpy(n1,ldy,y,n2,ldx,x,lda,a)
      implicit real*8 (a-h,o-z)
      dimension y(ldy,*),x(ldx,*),a(lda,*)
c
c
      j=mod(n2,2)
      if(j.ge.1)then
      do 10 i=1,n1
10    y(1,i)=(y(1,i))+x(1,j)*a(j,i)
      endif
      j=mod(n2,4)
      if(j.ge.2)then
      do 20 i=1,n1
20    y(1,i)=((y(1,i))+x(1,j-1)*a(j-1,i))+x(1,j)*a(j,i)
      endif
      jmin=j+4
      do 40 j=jmin,n2,4
      do 30 i=1,n1
      y(1,i)=((((y(1,i))+x(1,j-3)*a(j-3,i))+x(1,j-2)*a(j-2,i))
     &    +x(1,j-1)*a(j-1,i))+x(1,j)*a(j,i)
30    continue
40    continue
      return
      end
      subroutine f01aef(n, a, ia, b, ib, dl, ifail)
c     reduction of the general symmetric eigenvalue problem a * x =
c     lambda * b * x, with symmetric matrix a and symmetric
c     positive definite matrix b, to the equivalent standard
c     problem p * z = lambda * z. the upper triangle, including
c     diagonal elements, of a and b are given in the arrays a(n,n)
c     and b(n,n). l(b = l * lt) is formed in the
c     remaining strictly lower triangle of the array b with its
c     diagonal elements in the array dl(n) , and the lower
c     triangle of the symmetric matrix p(p = inv(l) * a * inv( lt))
c     is formed in the lower triangle of the array a, including the
c     diagonal elements. hence the diagonal elements of a are lost.
c     the subroutine will fail if b, perhaps on account of rounding
c     errors, is not positive definite - sets ifail = 0 if
c     successful else ifail = 1.
c
      implicit real*8  (a-h,o-z)
      character *8  srname
      dimension  a(ia,n), b(ib,n), dl(n)
      data srname /' f01aef '/
c     copy upper triangles of a and b into lower triangles
c     and diagonal of b into dl
      do 15 j=1,n
         do 10 i=1,j
            a(j,i) = a(i,j)
            b(j,i) = b(i,j)
   10    continue
         dl(j) = b(j,j)
   15 continue
      do 100 i=1,n
         x = b(i,i)
         if (x.le.0.0d0) go to 320
         x = dsqrt(x)
         b(i,i) = dl(i)
         dl(i) = x
         ip1 = i+1
         if (ip1.gt.n) go to 100
         do 80 j=ip1,n
            b(j,i) = b(j,i)/x
   80    continue
         call f03aez(b(i+1,1), ib, n-i, i, b(i+1,1), ib, b(i+1,i+1))
  100 continue
c     l has been formed in array b
      do 180 i=1,n
         y = dl(i)
         if (i.gt.1)
     *    call f03aez(a(i,1), ia, n-i+1, i-1, b(i,1), ib, a(i,i))
         do 160 j=i,n
            a(j,i) = a(j,i)/y
  160    continue
  180 continue
c     the transpose of the upper triangle of
c     inv(l) * a has been formed in the lower
c     triangle of array a
      do 300 j=1,n
         if (j.gt.1)
     *    call f03aez(b(j,1), ib, n-j+1, j-1, a(j,1), ia, a(j,j))
         call f01aez(b(j,j), ib, n-j+1, dl(j), a(j,j))
  300 continue
      ifail = 0
      return
  320 ifail = nagmsg(ifail,1,srname)
      do 340 j=i,n
         b(j,j) = dl(j)
  340 continue
      return
      end
      subroutine f01aez(a, ia, n, d, b)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  b = l**(-1) * b  where
c     a holds the subdiagonal elements of l and
c     d holds the diagonal elements of l.
c
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      real*8  a(ia,n), b(n), d(n)
c     ..
c
      n4 = 4*(n/4)
      do 40 j=1,n4,4
         b(j) = b(j)/d(j)
         b(j+1) = (b(j+1)-a(j+1,j)*b(j))/d(j+1)
         b(j+2) = ((b(j+2)-a(j+2,j)*b(j))-a(j+2,j+1)*b(j+1))/d(j+2)
         b(j+3) = (((b(j+3)-a(j+3,j)*b(j))-a(j+3,j+1)*b(j+1))-a(j+3,
     *    j+2)*b(j+2))/d(j+3)
         do 20 i=j+4,n
            b(i) = (((b(i)-a(i,j)*b(j))-a(i,j+1)*b(j+1))-a(i,j+2)*
     *       b(j+2)) - a(i,j+3)*b(j+3)
   20    continue
   40 continue
      if (n4+1.le.n) b(n4+1) = b(n4+1)/d(n4+1)
      if (n4+2.le.n) b(n4+2) = (b(n4+2)-a(n4+2,n4+1)*b(n4+1))
     * /d(n4+2)
      if (n4+3.le.n) b(n4+3) = ((b(n4+3)-a(n4+3,n4+1)*b(n4+1))
     * -a(n4+3,n4+2)*b(n4+2))/d(n4+3)
      return
      end
      subroutine f01aff(n, im1, im2, b, ib, dl, z, iz)
c     mark 11 revised. vectorisation (jan 1984).
c     this subroutine performs, on the matrix of eigenvectors, z,
c     stored
c     in columns im1 to im2 of the array z(n,n), a backward
c     substitution
c     lt* x = z, over- writing x on z. the diagonal elements of l
c     must be stored in the array dl(n), and the remaining
c     triangle in the strictly lower triangle of the array b(n,n).
c     the subroutines f01aef and f01bdf leave l in this
c     desired form. if x denotes any column of the resultant matrix
c     x, then x satisfies xt * b * x = zt * z, where b = l * lt.
      implicit real*8  (a-h,o-z)
      dimension b(ib,n), dl(n), z(iz,im2)
      do 80 j=im1,im2
         call f01afz(b, ib, n, dl, z(1,j))
   80 continue
      return
      end
      subroutine f01afz(a, ia, n, d, b)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  b = (l**t)**(-1) * b  where
c     a holds the subdiagonal elements of l and
c     d holds the diagonal elements of l.
c     the solution is overwritten on b.
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      dimension a(ia,n), b(n), d(n)
c     ..
c
      n4 = mod(n,4) + 1
      do 40 j=n,n4,-4
         b(j) = b(j)/d(j)
         b(j-1) = (b(j-1)-a(j,j-1)*b(j))/d(j-1)
         b(j-2) = ((b(j-2)-a(j,j-2)*b(j))-a(j-1,j-2)*b(j-1))/d(j-2)
         b(j-3) = (((b(j-3)-a(j,j-3)*b(j))-a(j-1,j-3)*b(j-1))
     *    -a(j-2,j-3)*b(j-2))/d(j-3)
         do 20 i=1,j-4
            b(i) = (((b(i)-a(j,i)*b(j))-a(j-1,i)*b(j-1))
     *       -a(j-2,i)*b(j-2)) - a(j-3,i)*b(j-3)
   20    continue
   40 continue
      if (n4.ge.2) b(n4-1) = b(n4-1)/d(n4-1)
      if (n4.ge.3) b(n4-2) = (b(n4-2)-a(n4-1,n4-2)*b(n4-1))/d(n4-2)
      if (n4.ge.4) b(n4-3) = ((b(n4-3)-a(n4-1,n4-3)*b(n4-1))
     * -a(n4-2,n4-3)*b(n4-2))/d(n4-3)
      return
      end
      subroutine f01agz(a, ia, n, b, c)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  c = a*b  where
c     a is a symmetric n-by-n matrix,
c     whose lower triangle is stored in a.
c     c must be distinct from b.
c
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      dimension a(ia,n), b(n), c(n)
c     ..
c
      n2 = mod(n,2)
      if (n2.eq.0) then
         do 20 i=1,n
            c(i) = 0.0d0
   20    continue
      else
         do 40 i=1,n
            c(i) = a(i,1)*b(1)
   40    continue
      endif
      do 100 j=n2+1,n,2
         do 60 i=1,j
            c(i) = (c(i)+a(j,i)*b(j)) + a(j+1,i)*b(j+1)
   60    continue
         do 80 i=j+1,n
            c(i) = (c(i)+a(i,j)*b(j)) + a(i,j+1)*b(j+1)
   80    continue
  100 continue
      return
      end
      subroutine f01ajf(n, atol, a, ia, d, e, z, iz)
c     mark 11 revised. vectorisation (jan 1984).
c     tred2
c     this subroutine reduces the given lower triangle of a
c     symmetric matrix, a, stored in the array a(n,n), to
c     tridiagonal form using householders reduction. the diagonal
c     of the result is stored in the array d(n) and the
c     sub-diagonal in the last n - 1 stores of the array e(n)
c     (with the additional element e(1) = 0). the transformation
c     matrices are accumulated in the array z(n,n). the array
c     a is left unaltered unless the actual parameters
c     corresponding to a and z are identical.
c
      implicit real*8  (a-h,o-z)
      dimension               a(ia,n), d(n), e(n), z(iz,n)
      do 40 i=1,n
         do 20 j=i,n
            z(j,i) = a(j,i)
   20    continue
         d(i) = a(n,i)
   40 continue
      if (n.eq.1) go to 440
      do 260 ii=2,n
         i = n - ii + 2
         l = i - 2
         f = d(i-1)
         g = 0.0d0
         if (l.eq.0) go to 80
         do 60 k=1,l
            g = g + d(k)*d(k)
   60    continue
   80    h = g + f*f
         l = l + 1
c     if g is too small for orthogonality to be
c     guaranteed the transformation is skipped
         if (g.gt.atol) go to 100
         e(i) = f
         h = 0.0d0
         do 85 j=1,l
            z(j,i) = d(j)
   85    continue
         do 90 j=1,l
            z(i,j) = 0.0d0
            d(j) = z(i-1,j)
   90    continue
         go to 240
  100    g = dsqrt(h)
         if (f.ge.0.0d0) g = -g
         e(i) = g
         h = h - f*g
         d(i-1) = f - g
c     copy u
         do 110 j=1,l
            z(j,i) = d(j)
  110    continue
c     form a*u
         call f01agz(z, iz, l, d, e)
c     form p
         do 182 j=1,l
            e(j) = e(j)/h
  182    continue
         f = 0.0d0
         do 185 j=1,l
            f = f + e(j)*d(j)
  185    continue
c     form k
         hh = f/(h+h)
c     form q
         do 190 j=1,l
            e(j) = e(j) - hh*d(j)
  190    continue
c     form reduced a
         do 220 j=1,l
            f = d(j)
            g = e(j)
            do 200 k=j,l
               z(k,j) = (z(k,j) - g*d(k)) - f*e(k)
  200       continue
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  220    continue
  240    d(i) = h
  260 continue
c     accumulation of transformation matrices
      do 400 i=2,n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h.eq.0.0d0) go to 360
         do 290 k=1,l
            d(k) = 0.0d0
  290    continue
         call f01cky(z, iz, l, l, z(1,i), 1, d)
         do 310 k=1,l
            d(k) = d(k)/h
  310    continue
         do 340 j=1,l
            do 320 k=1,l
               z(k,j) = z(k,j) - z(k,i)*d(j)
  320       continue
  340    continue
  360    do 380 j=1,l
            z(j,i) = 0.0d0
  380    continue
  400 continue
      do 420 i=1,n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  420 continue
  440 z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
      subroutine f01cky(a, ia, m, n, b, ib, c)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  c = c +  (a**t)*b  where
c     a is rectangular m by n.
c     c must be distinct from b.
c     the elements of b may be non-consecutive, with offset ib.
c
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      dimension a(ia,n), b(ib,m), c(n)
c
      m2 = 2*(m/2)
      m4 = 4*(m/4)
      do 60 j=1,m4,4
         do 40 i=1,n
            c(i) = (((c(i)+a(j,i)*b(1,j))+a(j+1,i)*b(1,j+1))+a(j+2,i)*
     *       b(1,j+2)) + a(j+3,i)*b(1,j+3)
   40    continue
   60 continue
      if (m2.gt.m4) then
         do 80 i=1,n
            c(i) = (c(i)+a(m2-1,i)*b(1,m2-1)) + a(m2,i)*b(1,m2)
   80    continue
      endif
      if (m.gt.m2) then
         do 100 i=1,n
            c(i) = c(i) + a(m,i)*b(1,m)
  100    continue
      endif
      return
      end
      subroutine f02abf(a, ia, n, r, v, iv, e, ifail)
c
c     eigenvalues and eigenvectors of a real symmetrix matrix
c     1st august 1971
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension       a(ia,n), r(n), v(iv,n), e(n)
      data srname /' f02abf '/
      isave = ifail
      ifail = 1
      tol = x02adf(xxxx)
      call f01ajf(n, tol, a, ia, r, e, v, iv)
      tol = x02aaf(xxxx)
      call f02amf(n, tol, r, e, v, iv, ifail)
      if (ifail.ne.0) ifail = nagmsg(isave,ifail,srname)
      return
      end
      subroutine f02aef(a, ia, b, ib, n, r, v, iv, dl, e, ifail)
c
c     eigenvalues and eigenvectors of a-lambda*b
c     1st december 1971
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension       a(ia,n), b(ib,n), r(n), v(iv,n), dl(n), e(n)
      data srname /' f02aef '/
      isave = ifail
      ifail = 1
      call f01aef(n, a, ia, b, ib, dl, ifail)
      if (ifail.eq.0) go to 20
      ifail = nagmsg(isave,ifail,srname)
      return
   20 tol = x02adf(xxxx)
      call f01ajf(n, tol, a, ia, r, e, v, iv)
      tol = x02aaf(xxxx)
      ifail = 1
      call f02amf(n, tol, r, e, v, iv, ifail)
      if (ifail.eq.0) go to 40
      ifail = nagmsg(isave,2,srname)
      return
   40 call f01aff(n, 1, n, b, ib, dl, v, iv)
      return
      end
      subroutine f02amf(n, acheps, d, e, z, iz, ifail)
c     mark 9 revised. ier-326 (sep 1981).
c
c     tql2
c     this subroutine finds the eigenvalues and eigenvectors of a
c     tridiagonal matrix, t, given with its diagonal elements in
c     the array d(n) and its sub-diagonal elements in the last n
c     - 1 stores of the array e(n), using ql transformations. the
c     eigenvalues are overwritten on the diagonal elements in the
c     array d in ascending order. the eigenvectors are formed in
c     the array z(n,n), overwriting the accumulated
c     transformations as supplied by the subroutine f01ajf. the
c     subroutine will fail if all eigenvalues take more than 30*n
c     iterations.
c     1st april 1972
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension                            d(n), e(n), z(iz,n)
      data srname /' f02amf ' /
      isave = ifail
      if (n.eq.1) go to 40
      do 20 i=2,n
         e(i-1) = e(i)
   20 continue
   40 e(n) = 0.0d0
      b = 0.0d0
      f = 0.0d0
      j = 30*n
      do 300 l=1,n
         h = acheps*(dabs(d(l))+dabs(e(l)))
         if (b.lt.h) b = h
c     look for small sub-diag element
         do 60 m=l,n
            if (dabs(e(m)).le.b) go to 80
   60    continue
   80    if (m.eq.l) go to 280
  100    if (j.le.0) go to 400
         j = j - 1
c     form shift
         g = d(l)
         h = d(l+1) - g
         if (dabs(h).ge.dabs(e(l))) go to 120
         p = h*0.5d0/e(l)
         r = dsqrt(p*p+1.0d0)
         h = p + r
         if (p.lt.0.0d0) h = p - r
         d(l) = e(l)/h
         go to 140
  120    p = 2.0d0*e(l)/h
         r = dsqrt(p*p+1.0d0)
         d(l) = e(l)*p/(1.0d0+r)
  140    h = g - d(l)
         i1 = l + 1
         if (i1.gt.n) go to 180
         do 160 i=i1,n
            d(i) = d(i) - h
  160    continue
  180    f = f + h
c     ql transformation
         p = d(m)
         c = 1.0d0
         s = 0.0d0
         m1 = m - 1
         do 260 ii=l,m1
            i = m1 - ii + l
            g = c*e(i)
            h = c*p
            if (dabs(p).lt.dabs(e(i))) go to 200
            c = e(i)/p
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
            go to 220
  200       c = p/e(i)
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*e(i)*r
            s = 1.0d0/r
            c = c/r
  220       p = c*d(i) - s*g
            d(i+1) = h + s*(c*g+s*d(i))
c     form vector
            do 240 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i) + c*h
               z(k,i) = c*z(k,i) - s*h
  240       continue
  260    continue
         e(l) = s*p
         d(l) = c*p
         if (dabs(e(l)).gt.b) go to 100
  280    d(l) = d(l) + f
  300 continue
c     order eigenvalues and eigenvectors
      do 380 i=1,n
         k = i
         p = d(i)
         i1 = i + 1
         if (i1.gt.n) go to 340
         do 320 j=i1,n
            if (d(j).ge.p) go to 320
            k = j
            p = d(j)
  320    continue
  340    if (k.eq.i) go to 380
         d(k) = d(i)
         d(i) = p
         do 360 j=1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  360    continue
  380 continue
      ifail = 0
      return
  400 ifail = nagmsg(isave,1,srname)
      return
      end
      subroutine f03aez(a, ia, m, n, b, ib, c)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  c = c - a*b  where
c     a is rectangular m by n.
c     the elements of b may be non-consecutive, with offset ib.
c
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      dimension a(ia,n), b(ib,n), c(m)
c
      n2 = 2*(n/2)
      n4 = 4*(n/4)
      do 60 j=1,n4,4
         do 40 i=1,m
            c(i) = (((c(i)-a(i,j)*b(1,j))-a(i,j+1)*b(1,j+1))-a(i,j+2)*
     *       b(1,j+2)) - a(i,j+3)*b(1,j+3)
   40    continue
   60 continue
      if (n2.gt.n4) then
         do 80 i=1,m
            c(i) = (c(i)-a(i,n2-1)*b(1,n2-1)) - a(i,n2)*b(1,n2)
   80    continue
      endif
      if (n.gt.n2) then
         do 100 i=1,m
            c(i) = c(i) - a(i,n)*b(1,n)
  100    continue
      endif
      return
      end
      subroutine f03aff(n, eps, a, ia, d1, id, p, ifail)
c     mark 11 revised. vectorisation (jan 1984).
c
c     unsymdet
c     the unsymmetric matrix, a, is stored in the n*n array a(i,j),
c     i=1,n, j=1,n. the decomposition a=lu, where l is a
c     lower triangular matrix and u is a unit upper triangular
c     matrix, is performed and overwritten on a, omitting the unit
c     diagonal of u. a record of any interchanges made to the rows
c     of a is kept in p(i), i=1,n, such that the i-th row and
c     the p(i)-th row were interchanged at the i-th step. the
c     determinant, d1 * 2.0**id, of a is also computed. the
c     subroutine
c     will fail if a, modified by the rounding errors, is singular
c     or almost singular. sets ifail = 0 if successful else ifail =
c     1.
c     1st december 1971
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension           a(ia,n), p(n)
      data srname /' f03aff '/
      do 10 i=1,n
         p(i) = 0.0d0
   10 continue
      do 15 j=1,n
         do 13 i=1,n
            p(i) = p(i) + a(i,j)**2
   13    continue
   15 continue
      do 20 i=1,n
         if (p(i).le.0.0d0) go to 200
         p(i) = 1.0d0/dsqrt(p(i))
   20 continue
      d1 = 1.0d0
      id = 0
      do 180 k=1,n
         l = k
         x = 0.0d0
         do 40 i=k,n
            y = dabs(a(i,k)*p(i))
            if (y.le.x) go to 40
            x = y
            l = i
   40    continue
         if (l.eq.k) go to 80
         d1 = -d1
         do 60 j=1,n
            y = a(k,j)
            a(k,j) = a(l,j)
            a(l,j) = y
   60    continue
         p(l) = p(k)
   80    p(k) = l
         d1 = d1*a(k,k)
         if (x.lt.8.0d0*eps) go to 200
  100    if (dabs(d1).lt.1.0d0) go to 120
         d1 = d1*0.0625d0
         id = id + 4
         go to 100
  120    if (dabs(d1).ge.0.0625d0) go to 140
         d1 = d1*16.0d0
         id = id - 4
         go to 120
  140    if (k.lt.n) call f03afz(a, ia, n, k, a(1,k+1))
  180 continue
      ifail = 0
      return
  200 ifail = nagmsg(ifail,1,srname)
      return
      end
      subroutine f03afz(a, ia, m, n, b)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  b = l**(-1) * b  where
c     l is a lower triangular matrix of the special form
c     ( l11  0 )
c     ( l21  i )
c     l11 is lower triangular of order n,
c     l21 is rectangular (m-n) by n.
c     a holds the elements of the first n columns of l.
c
c     .. array arguments ..
      implicit real*8  (a-h,o-z)
      dimension a(ia,n), b(m)
c
      n2 = 2*(n/2)
      n4 = 4*(n/4)
      do 40 j=1,n4,4
         b(j) = b(j)/a(j,j)
         b(j+1) = (b(j+1)-a(j+1,j)*b(j))/a(j+1,j+1)
         b(j+2) = ((b(j+2)-a(j+2,j)*b(j))-a(j+2,j+1)*b(j+1))/a(j+2,
     *    j+2)
         b(j+3) = (((b(j+3)-a(j+3,j)*b(j))-a(j+3,j+1)*b(j+1))-a(j+3,
     *    j+2)*b(j+2))/a(j+3,j+3)
         do 20 i=j+4,m
            b(i) = (((b(i)-a(i,j)*b(j))-a(i,j+1)*b(j+1))-a(i,j+2)*
     *       b(j+2)) - a(i,j+3)*b(j+3)
   20    continue
   40 continue
      if (n2.gt.n4) then
         b(n2-1) = b(n2-1)/a(n2-1,n2-1)
         b(n2) = (b(n2)-a(n2,n2-1)*b(n2-1))/a(n2,n2)
         do 60 i=n2+1,m
            b(i) = (b(i)-a(i,n2-1)*b(n2-1)) - a(i,n2)*b(n2)
   60    continue
      endif
      if (n.gt.n2) then
         b(n) = b(n)/a(n,n)
         do 80 i=n+1,m
            b(i) = b(i) - a(i,n)*b(n)
   80    continue
      endif
      return
      end
      function x02aaf(dumx)
c     * eps *
c     returns the value eps where eps is the smallest
c     positive
c     number such that 1.0 + eps > 1.0
c     the x parameter is not used
c     for fps ap-164
c     x02aaf = 2.0**(-51.0)
c     data tol/z'79c8000000000000'/
c     for ipsc  (experimental)
      implicit real*8  (a-h,o-z)
cjvl      x02aaf= 2.0d0**(-53) + 2.0d0**(-53)*1.0d-15 too accurate
c      x02aaf = 2.0d0**(-50)
      data t02aaf /0.278d-16/
c     x02aaf = 1.110223024625157E-016
      x02aaf = t02aaf
      return
      end
      function x02acf(dumx)
c     * rmax *
c     returns the value of the largest positive real  floating-
c     point number representable on the computer
c     for ipsc (experimental)
      implicit real*8  (a-h,o-z)
c     x02acf = 2.0d0**1023.0d0*(2.0d0-2.0d0**-52) really
c      x02acf = 2.0d0**1018.0d0
c     x02acf = 1.797693134862313E+308
      data t02acf/0.17d+39/
      x02acf = t02acf
      return
      end
      function x02adf(dumx)
c     * tol *
c     returns the ratio of the smallest positive real floating-
c     point number representable on the computer to eps
c     for fps ap-164
c     x02adf = 2.0**(-974.0)*(1.0+2.0**(-51))
c
c     ipsc (experimental)
c
      implicit real*8  (a-h,o-z)
c     x02adf=2.0**(-1022.0d0)/x02aaf really
c      x02adf=2.0**(-1018.0d0)/x02aaf
c     x02adf = 2.004168360009115E-292
      data t02adf/0.211758d-21/
      x02adf = t02adf
      return
      end
      subroutine nagout(i,nerr)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      if (i.eq.0) nerr = iwr
      if (i.eq.1) iwr = nerr
      return
      end
      function nagmsg(ifail, ifault, zrname)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      if (ifault.eq.0) go to 20
      call nagout (0,nout)
      if (mod(ifail,10).eq.1) go to 10
      write (nout,99999) zrname, ifault
      call caserr('error detected by nag library routine')
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nout,99999) zrname, ifault
   20 nagmsg = ifault
      return
99999 format (      ' error detected by nag library routine ', a8,
     * '   - ifail =', i5//)
      end
      subroutine f04arf(a, ia, b, n, c, wks, ifail)
c     approximate solution of a set of real linear
c     equations with one right hand side.
c     ist. april 1973
c
      implicit real*8 (a-h,o-z)
      character *8 srname
      dimension      a(ia,n), b(n), c(n), wks(n)
      data srname /' f04arf '/
      do 20 i=1,n
         c(i) = b(i)
   20 continue
      it = 1
      call f03aff(n, x02aaf(xxxx), a, ia, d1, i, wks, it)
      if (it.eq.0) go to 40
      ifail = nagmsg(ifail,1,srname)
      return
   40 ifail = 0
      call f04ajf(n, 1, a, ia, wks, c, n)
      return
      end
      subroutine f04ajf(n, ir, a, ia, p, b, ib)
c     mark 11 revised. vectorisation (jan 1984).
c
c     unsymsol
c     solves ax=b, where a is an unsymmetric matrix and b is an
c     n*ir
c     matrix of ir right-hand sides. the subroutine f04ajf must by
c     preceded by f03aff in which l and u are produced in a(i,j),
c     from a, and the record of the interchanges is produced in
c     p(i). ax=b is solved in three steps, interchange the
c     elements of b, ly=b and ux=y. the matrices y and then x are
c     overwritten on b.
c     1st august 1971
c
      implicit real*8 (a-h,o-z)
      dimension  a(ia,n), p(n), b(ib,ir)
c     interchanging of elements of b
      do 40 i=1,n
         i1 = p(i) + 0.5d0
         if (i1.eq.i) go to 40
         do 20 k=1,ir
            x = b(i,k)
            b(i,k) = b(i1,k)
            b(i1,k) = x
   20    continue
   40 continue
      do 100 k=1,ir
c     solution of ly= b
         call f03afz(a, ia, n, n, b(1,k))
c     solution of ux= y
         call f04ajz(a, ia, n, b(1,k))
  100 continue
      return
      end
      subroutine f04ajz(a, ia, n, b)
c     mark 11 release. nag copyright 1983
c     *** cray-1 version ***
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  b = u**(-1) * b  where
c     a holds the superdiagonal elements of u and
c     the diagonal elements of u are taken to be 1.0 .
c
c     .. array arguments ..
      implicit real*8 (a-h,o-z)
      dimension a(ia,n), b(n)
c
      n4 = mod(n,4) + 1
      do 40 j=n,n4,-4
         b(j-1) = b(j-1) - a(j-1,j)*b(j)
         b(j-2) = (b(j-2)-a(j-2,j)*b(j)) - a(j-2,j-1)*b(j-1)
         b(j-3) = ((b(j-3)-a(j-3,j)*b(j))-a(j-3,j-1)*b(j-1))
     *    - a(j-3,j-2)*b(j-2)
         do 20 i=1,j-4
            b(i) = (((b(i)-a(i,j)*b(j))-a(i,j-1)*b(j-1))
     *       -a(i,j-2)*b(j-2)) - a(i,j-3)*b(j-3)
   20    continue
   40 continue
      if (n4.ge.3) b(n4-2) = b(n4-2) - a(n4-2,n4-1)*b(n4-1)
      if (n4.ge.4) b(n4-3) = (b(n4-3)-a(n4-3,n4-1)*b(n4-1))
     * - a(n4-3,n4-2)*b(n4-2)
      return
      end
      subroutine a02acf(xxr, xxi, yyr, yyi, zr, zi)
      implicit real*8  (a-h,o-z)
c
c     divides one complex number by a second
c
      data one /1.0d0/
c
      if (dabs(yyr).le.dabs(yyi)) go to 20
      h = yyi/yyr
      a = one/(h*yyi+yyr)
      zr = (xxr+h*xxi)*a
      zi = (xxi-h*xxr)*a
      return
   20 h = yyr/yyi
      a = one/(h*yyr+yyi)
      zr = (h*xxr+xxi)*a
      zi = (h*xxi-xxr)*a
      return
      end
      subroutine f01akf(n, k, l, a, ia, intger)
c
c     dirhes
c     august 1st, 1971 .
c     given the unsymmetric matrix, a, stored in the array a(n,n),
c     this subroutine reduces the sub-matrix of order l - k + 1,
c     which starts at the element a(k,k) and finishes at the
c     element a(l,l), to hessenberg form, h, by the direct
c     method(an = nh). the matrix h is overwritten on a with
c     details of the transformations (n) stored in the remaining
c     triangle under h and in elements k to l of the array
c     intger(n).
c     1st august 1971
c
      implicit real*8  (a-h,o-z)
      dimension intger(n)
      dimension  a(ia,n)
      k1 = k + 1
      if (k1.gt.l) return
      do 200 j=k1,n
         m = j
         x = 0.0d0
         if (j.gt.l) go to 120
         do 20 i=j,l
            if (dabs(a(i,j-1)).le.dabs(x)) go to 20
            x = a(i,j-1)
            m = i
   20    continue
         intger(j) = m
         if (m.eq.j) go to 80
c     interchange rows and columns of a.
         do 40 i=k,n
            y = a(m,i)
            a(m,i) = a(j,i)
            a(j,i) = y
   40    continue
         do 60 i=1,l
            y = a(i,m)
            a(i,m) = a(i,j)
            a(i,j) = y
   60    continue
   80    if (x.eq.0.0d0) go to 120
         if (j.ge.l) go to 120
         j1 = j + 1
         do 100 i=j1,l
            a(i,j-1) = a(i,j-1)/x
  100    continue
         call f01ckz(a(1,j+1), ia, l, l-j, a(j+1,j-1), 1, a(1,j))
  120    call f01aky(a(k+1,k), ia, l-k, j-k, a(k+1,j))
  200 continue
      return
      end
      subroutine f01aky(a, ia, m, n, b)
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  b = l**(-1) * b  where
c     l is a lower triangular matrix of the special form
c     ( l11  0 )
c     ( l21  i )
c     l11 is lower triangular of order n,
c     l21 is rectangular (m-n) by n.
c     a holds the subdiagonal elements of the first n columns of l
c     and the diagonal elements of l are taken to be 1.0 .
c
      implicit real*8  (a-h,o-z)
c     .. scalar arguments ..
c     integer ia, m, n
c     .. array arguments ..
      dimension a(ia,n), b(m)
c     ..
c     .. local scalars ..
c     integer i, j, n2, n4
c     ..
c
      n2 = 2*(n/2)
      n4 = 4*(n/4)
      do 40 j=1,n4,4
         b(j+1) = b(j+1)-a(j+1,j)*b(j)
         b(j+2) = (b(j+2)-a(j+2,j)*b(j))-a(j+2,j+1)*b(j+1)
         b(j+3) = ((b(j+3)-a(j+3,j)*b(j))-a(j+3,j+1)*b(j+1))-a(j+3,
     *    j+2)*b(j+2)
         do 20 i=j+4,m
            b(i) = (((b(i)-a(i,j)*b(j))-a(i,j+1)*b(j+1))-a(i,j+2)*
     *       b(j+2)) - a(i,j+3)*b(j+3)
   20    continue
   40 continue
      if (n2.gt.n4) then
         b(n2) = b(n2)-a(n2,n2-1)*b(n2-1)
         do 60 i=n2+1,m
            b(i) = (b(i)-a(i,n2-1)*b(n2-1)) - a(i,n2)*b(n2)
   60    continue
      endif
      if (n.gt.n2) then
         do 80 i=n+1,m
            b(i) = b(i) - a(i,n)*b(n)
   80    continue
      endif
      return
      end
      subroutine f01apf(n, low, iupp, intger, h, ih, v, iv)
c
c     dirtrans
c     form the matrix of accumulated transformations in the array
c     v(n,n) from the information left by subroutine f01akf
c     below the upper hessenberg matrix, h, in the array h(n,n)
c     and in the integer array intger(n).
c     1st august 1971
c
      implicit real*8  (a-h,o-z)
      dimension  intger(n)
      dimension  h(ih,n), v(iv,n)
c
      do 40 i=1,n
         do 20 j=1,n
            v(i,j) = 0.0d0
   20    continue
         v(i,i) = 1.0d0
   40 continue
      low1 = low + 1
      if (low1.gt.iupp) return
      do 120 ii=low1,iupp
         i = low1 + iupp - ii
         i1 = i - 1
         if (low1.gt.i1) go to 80
         do 60 j=low1,i1
            v(i,j) = h(i,j-1)
   60    continue
   80    m = intger(i)
         if (m.eq.i) go to 120
         do 100 j=low1,iupp
            x = v(m,j)
            v(m,j) = v(i,j)
            v(i,j) = x
  100    continue
  120 continue
      return
      end
      subroutine f01bdf(n, a, ia, b, ib, dl, ifail)
c     reduc2
c     reduction of the general symmetric eigenvalue problems
c     a*b*x=lambda*x,  yt*a*b=yt*lambda,
c     b*a*y=lambda*y,  xt*b*a=xt*lambda,
c     with symmetric matrix a and symmetric positive definite
c     matrix
c     b, to the equivalent standard problem q*z=lambda*z.
c     the upper triangle, including diagonal elements, of a and b
c     are given in the arrays a(n,n) and b(n,n).
c     l (b=l*lt) is formed in the remaining strictly lower triangle
c     of the array b with its diagonal elements in the array dl(n),
c     and the lower triangle of the symmetric matrix q (q=lt*a*l)
c     is formed in the lower triangle of the array a, including the
c     diagonal elements. hence the diagonal elements of a are lost.
c     the subroutine will fail if b, perhaps on account of rounding
c     errors, is not positive definite.
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension  a(ia,n), b(ib,n), dl(n)
      data srname /' f01bdf '/
      isave = ifail
      do 100 i=1,n
         i1 = i - 1
         do 80 j=i,n
            x = b(i,j)
            if (i1.eq.0) go to 40
            do 20 kk=1,i1
               k = i1 - kk + 1
               x = x - b(i,k)*b(j,k)
   20       continue
   40       if (i.ne.j) go to 60
            if (x.lt.0.0d0) go to 320
            y = dsqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x/y
   80    continue
  100 continue
c     l has been formed in array b
      do 220 i=1,n
         i1 = i + 1
         do 200 j=1,i
            x = a(j,i)*dl(j)
            j1 = j + 1
            if (j1.gt.i) go to 140
            do 120 k=j1,i
               x = x + a(k,i)*b(k,j)
  120       continue
  140       if (i1.gt.n) go to 180
            do 160 k=i1,n
               x = x + a(i,k)*b(k,j)
  160       continue
  180       a(i,j) = x
  200    continue
  220 continue
c     the lower triangle of a*l has been formed
c     in the lower triangle of array a
      do 300 i=1,n
         y = dl(i)
         i1 = i + 1
         do 280 j=1,i
            x = y*a(i,j)
            if (i1.gt.n) go to 260
            do 240 k=i1,n
               x = x + a(k,j)*b(k,i)
  240       continue
  260       a(i,j) = x
  280    continue
  300 continue
      ifail = 0
      return
  320 ifail = nagmsg(isave,1,srname)
      return
      end
      subroutine f01ckz(a, ia, m, n, b, ib, c)
c
c     this version uses the technique of unrolling vector
c     loops, due to dongarra and eisenstat.
c     also uses fortran 77.
c
c     computes  c = c + a*b  where
c     a is rectangular m by n.
c     c must be distinct from b.
c     the elements of b may be non-consecutive, with offset ib.
c
      implicit real*8  (a-h,o-z)
c     .. scalar arguments ..
c     integer ia, ib, m, n
c     .. array arguments ..
      dimension a(ia,n), b(ib,n), c(m)
c     ..
c     .. local scalars ..
c     integer i, j, n2, n4
c     ..
c
      n2 = 2*(n/2)
      n4 = 4*(n/4)
      do 60 j=1,n4,4
         do 40 i=1,m
            c(i) = (((c(i)+a(i,j)*b(1,j))+a(i,j+1)*b(1,j+1))+a(i,j+2)*
     *       b(1,j+2)) + a(i,j+3)*b(1,j+3)
   40    continue
   60 continue
      if (n2.gt.n4) then
         do 80 i=1,m
            c(i) = (c(i)+a(i,n2-1)*b(1,n2-1))+a(i,n2)*b(1,n2)
   80    continue
      endif
      if (n.gt.n2) then
         do 100 i=1,m
            c(i) = c(i) + a(i,n)*b(1,n)
  100    continue
      endif
      return
      end
      subroutine f02agf(a, ia, n, rr, ri, vr, ivr, vi, ivi, intger,
     * ifail)
c
c     eigenvalues and eigenvectors of real unsymmetric matrix
c     1st august 1971
c
      implicit real*8  (a-h,o-z)
      dimension  intger(n)
      character *8 srname
      real*8  machep, max
      dimension  a(ia,n), vr(ivr,n),
     * vi(ivi,n), rr(n), ri(n)
      data srname /' f02agf '/
      isave = ifail
      ifail = 1
      machep = x02aaf(xxxx)
      call f01akf(n, 1, n, a, ia, intger)
      call f01apf(n, 1, n, intger, a, ia, vr, ivr)
      call f02aqf(n, 1, n, machep, a, ia, vr, ivr, rr, ri, intger,
     * ifail)
      if (ifail.eq.0) go to 20
      ifail = nagmsg(isave,ifail,srname)
      return
   20 do 140 i=1,n
         if (ri(i).eq.0.0d0) go to 60
         if (ri(i).gt.0.0d0) go to 100
         do 40 j=1,n
            vr(j,i) = vr(j,i-1)
            vi(j,i) = -vi(j,i-1)
   40    continue
         go to 140
   60    do 80 j=1,n
            vi(j,i) = 0.0d0
   80    continue
         go to 140
  100    do 120 j=1,n
            vi(j,i) = vr(j,i+1)
  120    continue
  140 continue
      do 280 i=1,n
         sum = 0.0d0
         max = 0.0d0
         do 180 j=1,n
            if (dabs(vr(j,i)).le.max) go to 160
            max = dabs(vr(j,i))
  160       if (dabs(vi(j,i)).le.max) go to 180
            max = dabs(vi(j,i))
  180    continue
         do 200 j=1,n
            vr(j,i) = vr(j,i)/max
            vi(j,i) = vi(j,i)/max
  200    continue
         max = 0.0d0
         do 240 j=1,n
            term = vr(j,i)**2 + vi(j,i)**2
            sum = sum + term
            if (term.le.max) go to 220
            max = term
            c = vr(j,i)
            d = -vi(j,i)
  220       continue
  240    continue
         sum = sum*(c**2+d**2)
         sum = dsqrt(sum)
         do 260 j=1,n
            term = vr(j,i)
            vr(j,i) = (vr(j,i)*c-vi(j,i)*d)/sum
            vi(j,i) = (d*term+c*vi(j,i))/sum
  260    continue
  280 continue
      return
      end
      subroutine f02aqf(n, low, upp, machep, h, ih, vecs, ivecs,
     * wr, wi, cnt, ifail)
c
c     hqr2
c     finds the eigenvalues and eigenvectors of a real matrix
c     which has been reduced to upper hessenberg form in the array
c     h(n,n) with the accumulated transformations stored in
c     the array vecs(n,n). the real and imaginary parts of the
c     eigenvalues are formed in the arrays wr, wi(n) and the
c     eigenvectors are formed in the array vecs(n,n) where
c     only one complex vector, corresponding to the root with
c     positive imaginary part, is formed for a complex pair. low
c     and upp are two integers produced in balancing where
c     eigenvalues are isolated in positions 1 to low-1 and upp+1
c     to n. if balancing is not used low=1, upp=n. machep is the
c     relative machine precision. the subroutine will fail if
c     all eigenvalues take more than 30*n iterations.
c     1st december 1971
c
      implicit real*8  (a-h,o-z)
c
      integer upp, en, en2, upp1, cnt
      dimension  cnt(n)
      character *8 srname
      real*8  norm, machep
      dimension  h(ih,n), vecs(ivecs,n), wr(n), wi(n)
      logical notlas
      data srname /' f02aqf '/
      isave = ifail
c     compute matrix norm
      norm = 0.0d0
      k = 1
      do 10 i = 1,n
         do 5 j = k,n
            norm = norm + dabs(h(i,j))
    5    continue
         k = i
   10 continue
      nhs = n*(n+1)/2 +n -1
c     isolated roots
      if (low.le.1) go to 40
      j = low - 1
      do 20 i=1,j
         wr(i) = h(i,i)
         wi(i) = 0.0d0
         cnt(i) = 0
   20 continue
   40 if (upp.ge.n) go to 80
      j = upp + 1
      do 60 i=j,n
         wr(i) = h(i,i)
         wi(i) = 0.0d0
         cnt(i) = 0
   60 continue
   80 en = upp
      t = 0.0d0
      itn = 30*n
  100 if (en.lt.low) go to 740
      its = 0
      na = en - 1
c     look for single small sub-diagonal element
  120 if (low+1.gt.en) go to 160
      low1 = low + 1
      do 140 ll=low1,en
         l = en + low1 - ll
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s.lt.x02adf(0.0d0)) s = norm/dble(nhs)
         if (dabs(h(l,l-1)).le.machep*s) go to 180
  140 continue
  160 l = low
  180 x = h(en,en)
      if (l.eq.en) go to 600
      y = h(na,na)
      w = h(en,na)*h(na,en)
      if (l.eq.na) go to 620
      if (itn.le.0) go to 1480
c     form shift
      if ((its.ne.10) .and. (its.ne.20)) go to 240
      t = t + x
      if (low.gt.en) go to 220
      do 200 i=low,en
         h(i,i) = h(i,i) - x
  200 continue
  220 s = dabs(h(en,na)) + dabs(h(na,en-2))
      x = 0.75d0*s
      y = x
      w = -0.4375d0*s**2
  240 its = its + 1
      itn = itn - 1
c     look for two consecutive small sub-diagonal elements
      if (l.gt.en-2) go to 280
      en2 = en - 2
      do 260 mm=l,en2
         m = l + en2 - mm
         z = h(m,m)
         r = x - z
         s = y - z
         p = (r*s-w)/h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - z - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p/s
         q = q/s
         r = r/s
       if (m.eq.l) go to 280
       if ((dabs(h(m,m-1))*(dabs(q)+dabs(r))).le.(machep*dabs(p)*
     *  (dabs(h(m-1,m-1))+dabs(z)+dabs(h(m+1,m+1))))) go to 280
  260 continue
  280 m2 = m + 2
      if (m2.gt.en) go to 320
      do 300 i=m2,en
         h(i,i-2) = 0.0d0
  300 continue
  320 m3 = m + 3
      if (m3.gt.en) go to 360
      do 340 i=m3,en
         h(i,i-3) = 0.0d0
  340 continue
  360 if (m.gt.na) go to 580
      do 560 k=m,na
         notlas = k.ne.na
         if (k.eq.m) go to 380
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x.eq.0.0d0) go to 560
         p = p/x
         q = q/x
         r = r/x
  380    s = dsqrt(p**2+q**2+r**2)
         if (p.lt.0.0d0) s = -s
         if (k.ne.m) go to 400
         if (l.ne.m) h(k,k-1) = -h(k,k-1)
         go to 420
  400    h(k,k-1) = -s*x
  420    p = p + s
         x = p/s
         y = q/s
         z = r/s
         q = q/p
         r = r/p
c     row modification
         if (notlas) go to 440
         do 430 j=k,n
            p = h(k,j) + q*h(k+1,j)
            h(k+1,j) = h(k+1,j) - p*y
            h(k,j) = h(k,j) - p*x
  430    continue
         go to 465
  440    do 460 j=k,n
            p = h(k,j) + q*h(k+1,j) + r*h(k+2,j)
            h(k+2,j) = h(k+2,j) - p*z
            h(k+1,j) = h(k+1,j) - p*y
            h(k,j) = h(k,j) - p*x
  460    continue
  465    j = en
         if (k+3.lt.en) j = k + 3
c     column modification
         if (notlas) go to 480
         do 470 i=1,j
            p = x*h(i,k) + y*h(i,k+1)
            h(i,k+1) = h(i,k+1) - p*q
            h(i,k) = h(i,k) - p
  470    continue
         go to 505
  480    do 500 i=1,j
            p = x*h(i,k) + y*h(i,k+1) + z*h(i,k+2)
            h(i,k+2) = h(i,k+2) - p*r
            h(i,k+1) = h(i,k+1) - p*q
            h(i,k) = h(i,k) - p
  500    continue
c     accumulate transformations
  505    if (low.gt.upp) go to 560
         if (notlas) go to 520
         do 510 i=low,upp
            p = x*vecs(i,k) + y*vecs(i,k+1)
            vecs(i,k+1) = vecs(i,k+1) - p*q
            vecs(i,k) = vecs(i,k) - p
  510    continue
         go to 560
  520    do 540 i=low,upp
            p = x*vecs(i,k) + y*vecs(i,k+1) + z*vecs(i,k+2)
            vecs(i,k+2) = vecs(i,k+2) - p*r
            vecs(i,k+1) = vecs(i,k+1) - p*q
            vecs(i,k) = vecs(i,k) - p
  540    continue
  560 continue
  580 go to 120
c     one root found
  600 wr(en) = x + t
      h(en,en) = wr(en)
      wi(en) = 0.0d0
      cnt(en) = its
      en = na
      go to 100
c     two roots found
  620 p = (y-x)*0.5d0
      q = p**2 + w
      z = dsqrt(dabs(q))
      x = x + t
      h(en,en) = x
      h(na,na) = y + t
      cnt(en) = -its
      cnt(na) = its
      if (q.lt.0.0d0) go to 700
c     real pair
      if (p.lt.0.0d0) z = p - z
      if (p.gt.0.0d0) z = p + z
      wr(na) = x + z
      wr(en) = wr(na)
      if (z.ne.0.0d0) wr(en) = x - w/z
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      r = a02abf(x,z)
      p = x/r
      q = z/r
c     row modification
      do 640 j=na,n
         z = h(na,j)
         h(na,j) = q*z + p*h(en,j)
         h(en,j) = q*h(en,j) - p*z
  640 continue
c     column modification
      do 660 i=1,en
         z = h(i,na)
         h(i,na) = q*z + p*h(i,en)
         h(i,en) = q*h(i,en) - p*z
  660 continue
c     accumulate transformations
      do 680 i=low,upp
         z = vecs(i,na)
         vecs(i,na) = q*z + p*vecs(i,en)
         vecs(i,en) = q*vecs(i,en) - p*z
  680 continue
      go to 720
c     complex pair
  700 wr(na) = x + p
      wr(en) = x + p
      wi(na) = z
      wi(en) = -z
  720 en = en - 2
      go to 100
c     all roots found now backsubstitute
  740 if (norm.eq.0.0d0) go to 1460
      norm = norm*machep
c     backsubstitution
      do 1160 kk=1,n
         en = n + 1 - kk
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q.ne.0.0d0) go to 980
c     real vector
         h(en,en) = 1.0d0
         if (na.lt.1) go to 1160
         do 960 ii=1,na
            i = na + 1 - ii
            i1 = i - 1
            w = h(i,i) - p
            r = h(i,en)
            if (wi(i).ge.0.0d0) go to 840
            z = w
            s = r
            go to 960
  840       if (wi(i).gt.0.0d0) go to 920
c     modification to stop overflow
            if (w.ne.0.0d0) go to 880
            if (dabs(r).lt.10.0d0*norm) go to 900
            r = -r
            do 860 j=1,en
               h(j,en) = h(j,en)*norm
  860       continue
            go to 910
  880       r = -r/w
            go to 910
  900       r = -r/norm
  910       h(i,en) = r
            if (i1.eq.0) go to 960
            do 915 j=1,i1
               h(j,en) = h(j,en) + h(j,i)*r
  915       continue
            go to 960
c     solve real equations
  920       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i)-p)**2 + wi(i)**2
            t = (x*s-z*r)/q
            h(i,en) = t
            if (dabs(x).gt.dabs(z)) go to 930
            r = (-s-y*t)/z
            go to 940
  930       r = (-r-w*t)/x
  940       h(i+1,en) = r
            if (i1.eq.0) go to 960
            do 950 j=1,i1
               h(j,en) = (h(j,en)+h(j,i+1)*r) + h(j,i)*t
  950       continue
  960    continue
c     end real vector
         go to 1160
  980    if (q.gt.0.0d0) go to 1160
c     complex vector associated with lambda=p-i*q
         if (dabs(h(en,na)).le.dabs(h(na,en))) go to 1000
         r = q/h(en,na)
         s = -(h(en,en)-p)/h(en,na)
         go to 1020
 1000    call a02acf(0.0d0, -h(na,en), h(na,na)-p, q, r, s)
 1020    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         h(na,na) = r
         h(na,en) = s
         if (na.lt.2) go to 1160
         na1 = na - 1
         do 1030 j=1,na1
            h(j,en) = h(j,en) + h(j,na)*s
            h(j,na) = h(j,na)*r
 1030    continue
         do 1140 ii=1,na1
            i = 1 + na1 - ii
            i1 = i - 1
            w = h(i,i) - p
            ra = h(i,na)
            sa = h(i,en)
            if (wi(i).ge.0.0d0) go to 1080
            z = w
            r = ra
            s = sa
            go to 1140
 1080       if (wi(i).eq.0.0d0) go to 1120
c     solve complex equations
            x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i)-p)**2 + wi(i)**2 - q**2
            vi = (wr(i)-p)*2.0d0*q
            if ((vr.eq.0.0d0) .and. (vi.eq.0.0d0)) vr =
     *       machep*norm*(dabs(w)+dabs(q)+dabs(x)+dabs(y)+dabs(z))
            call a02acf(x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi,
     *       t, u)
            if (dabs(x).le.dabs(z)+dabs(q)) go to 1100
            r = (-ra-w*t+q*u)/x
            s = (-sa-w*u-q*t)/x
            go to 1110
 1100       call a02acf(-r-y*t, -s-y*u, z, q, r, s)
 1110       h(i,na) = t
            h(i,en) = u
            h(i+1,na) = r
            h(i+1,en) = s
            if (i1.eq.0) go to 1140
            do 1115 j=1,i1
               h(j,na) = (h(j,na)+h(j,i+1)*r) + h(j,i)*t
               h(j,en) = (h(j,en)+h(j,i+1)*s) + h(j,i)*u
 1115       continue
            go to 1140
 1120       call a02acf(-ra, -sa, w, q, r, s)
            h(i,na) = r
            h(i,en) = s
            if (i1.eq.0) go to 1140
            do 1130 j=1,i1
               h(j,na) = h(j,na) + h(j,i)*r
               h(j,en) = h(j,en) + h(j,i)*s
 1130       continue
 1140    continue
c     end complex vector
 1160 continue
c     end backsubstitution
c     vectors of isolated roots
      low1 = low - 1
      upp1 = upp + 1
      do 1260 j=1,n
         m = min(j,low1)
         if (m.lt.1) go to 1220
         do 1180 i=1,m
            vecs(i,j) = h(i,j)
 1180    continue
 1220    if (upp1.gt.j) go to 1260
         do 1240 i=upp1,j
            vecs(i,j) = h(i,j)
 1240    continue
 1260 continue
c     multiply by transformation matrix to give
c     vectors of original full matrix
      do 1440 jj=low,n
         j = low + n - jj
         m = min(j,upp)
         do 1360 i=low,upp
            vecs(i,j) = vecs(i,m)*h(m,j)
 1360    continue
         m = m - 1
         if (low.le.m) call f01ckz(vecs(low,low), ivecs, upp-low+1,
     *    m-low+1, h(low,j), 1, vecs(low,j))
 1440 continue
 1460 ifail = 0
      return
 1480 ifail = nagmsg(isave,1,srname)
      return
      end
      subroutine f04aef(a, ia, b, ib, n, m, c, ic, wkspce, aa, iaa,
     * bb, ibb, ifail)
c
c     accurate solution of a set of real linear equations
c     with multiple right hand sides.
c     1st august 1971
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension a(ia,n), b(ib,m), c(ic,m), aa(iaa,n),
     * wkspce(n), bb(ibb,m)
      data srname /' f04aef '/
      isave = ifail
      ifail = 1
      eps = x02aaf(xxxx)
c     copy a to working array aa
      do 40 i=1,n
         do 20 j=1,n
            aa(i,j) = a(i,j)
   20    continue
   40 continue
      call f03aff(n, eps, aa, iaa, d1, i, wkspce, ifail)
      if (ifail.eq.0) go to 60
      ifail = nagmsg(isave,ifail,srname)
      return
   60 ifail = 1
      call f04ahf(n, m, a, ia, aa, iaa, wkspce, b, ib, eps, c, ic,
     * bb, ibb, i, ifail)
      if (ifail.eq.0) go to 80
      ifail = nagmsg(isave,2,srname)
   80 return
      end
      subroutine f04ahf(n, ir, a, ia, aa, iaa, p, b, ib, eps, x,
     * ix, bb, ibb, l, ifail)
c
c     unsymaccsolve
c     solves ax=b where a is an n*n unsymmetric matrix and b is an
c     n*ir matrix of right hand sides, using the subroutine f04ajf.
c     the subroutine must by preceded by f03aff in which l and u
c     are produced in aa(i,j) and the interchanges in p(i). the
c     residuals bb=b-ax are calculated and ad=bb is solved, over-
c     writing d on bb. the refinement is repeated, as long as the
c     maximum correction at any stage is less than half that at the
c     previous stage, until the maximum correction is less than 2
c     eps times the maximum x. sets ifail = 1 if the solution fails
c     to improve, else ifail = 0. l is the number of iterations.
c     additional precision innerproducts are absolutely necessary.
c     1st december 1971
c
      implicit real*8  (a-h,o-z)
      character *8 srname
      dimension  a(ia,n), aa(iaa,n),
     * p(n), b(ib,ir), x(ix,ir), bb(ibb,ir)
      data srname /' f04ahf '/
      isave = ifail
      ifail1 = 0
      do 40 j=1,ir
         do 20 i=1,n
            x(i,j) = 0.0d0
            bb(i,j) = b(i,j)
   20    continue
   40 continue
      l = 0
      d0 = 0.0d0
   60 call f04ajf(n, ir, aa, iaa, p, bb, ibb)
      l = l + 1
      id2 = 0
      d1 = 0.0d0
      do 100 j=1,ir
         do 80 i=1,n
            x(i,j) = x(i,j) + bb(i,j)
   80    continue
  100 continue
      do 140 j=1,ir
         xmax = 0.0d0
         bbmax = 0.0d0
         do 120 i=1,n
            if (dabs(x(i,j)).gt.xmax) xmax = dabs(x(i,j))
            if (dabs(bb(i,j)).gt.bbmax) bbmax = dabs(bb(i,j))
            bb(i,j) =  b(i,j)-ddot(n,a(i,1),ia,x(1,j),1)
  120    continue
         if (bbmax.gt.d1*xmax) d1 = bbmax/xmax
         if (bbmax.gt.2.0d0*eps*xmax) id2 = 1
  140 continue
      if ((d1.gt.0.5d0*d0) .and. (l.ne.1)) go to 160
      d0 = d1
      if (id2.eq.1) go to 60
      ifail = 0
      return
  160 ifail = nagmsg(isave,1,srname)
      return
      end
      subroutine m01ajf(a, w, ind, indw, n, nw, ifail)
c
c     sorts a vector of real numbers into ascending order and
c     provides an in
c     indicating the position of the sorted numbers in the original
c     array.
      implicit real*8  (a-h,o-z)
      dimension ind(n), indw(nw), kount(20), ist(20)
      character *8 srname
      dimension a(n), w(nw)
      data srname /' m01ajf '/
      k = 0
      if (n-1) 20, 680, 40
   20 k = 1
   40 if (nw) 640, 640, 60
   60 if (k) 660, 80, 660
   80 nn = min(n,nw)
      do 100 k=1,n
         ind(k) = k
  100 continue
      l3 = 1
      iko = 0
  120 iko = iko + 1
  140 nm = 1
      l1 = l3
  160 l2 = min(l1+nm,n+1)
      l3 = min(l2+nm,n+1)
      kount(iko) = l3 - l1
      ist(iko) = l1
  180 iw = 1
      il1 = l1
      il2 = l2
  200 ndo = min(il2-l1,l3-l2,nn-iw+1)
      if (ndo) 320, 320, 220
  220 do 300 i=1,ndo
         if (a(l1)-a(l2)) 240, 240, 260
  240    w(iw) = a(l1)
         indw(iw) = ind(l1)
         l1 = l1 + 1
         go to 280
  260    w(iw) = a(l2)
         indw(iw) = ind(l2)
         l2 = l2 + 1
  280    iw = iw + 1
  300 continue
      go to 200
  320 if (il2-l1) 340, 340, 460
  340 iw1 = iw - 1
  360 do 380 i=1,iw1
         a(il1) = w(i)
         ind(il1) = indw(i)
         il1 = il1 + 1
  380 continue
  400 if (iko-1) 420, 420, 580
  420 if (l3-n) 440, 440, 700
  440 iko = 2
      go to 140
  460 if (l3-l2) 480, 480, 520
  480 j = l3 - 1
      k = il2 - l1
      j1 = il2 - 1
      do 500 i=1,k
         a(j) = a(j1)
         ind(j) = ind(j1)
         j = j - 1
         j1 = j1 - 1
  500 continue
      iw1 = iw - 1
      if (iw1) 400, 400, 360
  520 j = il2 - 1
      k = l2 - 1
      n21 = il2 - l1
      do 540 i=1,n21
         a(k) = a(j)
         ind(k) = ind(j)
         k = k - 1
         j = j - 1
  540 continue
      do 560 i=1,nn
         a(il1) = w(i)
         ind(il1) = indw(i)
         il1 = il1 + 1
  560 continue
      l1 = k + 1
      go to 180
  580 if (l3-n) 600, 600, 620
  600 if (kount(iko-1)-nm-nm) 620, 620, 120
  620 iko = iko - 1
      nm = kount(iko)
      l1 = ist(iko)
      go to 160
  640 k = k + 2
  660 ifail = nagmsg(ifail,k,srname)
      return
  680 ind(1) = 1
  700 ifail = 0
      return
      end
      function a02abf(xxr, xxi)
c
      implicit real*8  (a-h,o-z)
c
c     returns the absolute value of a complex number via routine
c     name
c
      data zero /0.0d0/, one /1.0d0/
c
      xr = dabs(xxr)
      xi = dabs(xxi)
      if (xi.le.xr) go to 20
      h = xr
      xr = xi
      xi = h
   20 if (xi.ne.zero) go to 40
      a02abf = xr
      return
   40 h = xr*dsqrt(one+(xi/xr)**2)
      a02abf = h
      return
      end
      function x02agf(dumx)
c
c     returns the smallest positive floating-point number  r
c     exactly representable on the computer such that  -r, 1.0/r,
c     and -1.0/r can all be computed without overflow or underflow.
c     on many machines the correct value can be derived from those
c     of x02aaf, x02abf and x02acf as follows
c
c     if (x02abf(x)*x02acf(x).ge.1.0) x02agf = x02abf(x)
c     if (x02abf(x)*x02acf(x).lt.1.0)
c    *                            x02agf = (1.0+x02aaf(x))/x02acf(x)
c
c     the correct value should be defined as a constant,
c     possibly in some binary, octal or hexadecimal representation,
c     and inserted into the assignment statement below.
c
c     dumx is a dummy argument
c
c     for fps ap-164
c     x02agf = 2.0**(-1023.0)*(1.0+2.0**(-51))
c
c     ipsc (experimental)
      implicit real*8  (a-h,o-z)
c
c*f    data tol/z'0048000000000001'/
c      x02agf = 2.0d0**(-1022.0d0)
c     x02agf = 2.225073858507266E-308
      data t02agf /3.d-308/
      x02agf = t02agf
      return
      end
      subroutine ver_nag(s,r,d)
      character*80 source
      character*30 revision
      character*30 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/nag.m,v $
     +     "/
      data revision /"$Revision: 1.15 $"/
      data date /"$Date: 2006-01-13 17:54:55 $"/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
