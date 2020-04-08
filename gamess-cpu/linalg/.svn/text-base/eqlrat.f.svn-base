c***********************************************************************
c*       routine eqlrat(n,diag,e,e2in,d,ind,ierr,e2)
c*
c*    authors -
c*       this is a modification of routine eqlrat from eispack edition 3
c*       dated august 1983.
c*       tqlrat is a translation of the algol procedure tqlrat,
c*       algorithm 464, comm. acm 16, 689(1973) by reinsch.
c*       this version is by s. t. elbert (ames laboratory-usdoe)
c*
c*    purpose -
c*       this routine finds the eigenvalues of a symmetric
c*       tridiagonal matrix
c*
c*    method -
c*       rational ql
c*
c*    on entry -
c*       n      - integer*4
c*                the order of the matrix.
c*       d      - w.p. real (n)
c*                contains the diagonal elements of the input matrix.
c*       e2     - w.p. real (n)
c*                contains the squares of the subdiagonal elements of
c*                the input matrix in its last n-1 positions.
c*                e2(1) is arbitrary.
c*
c*     on exit -
c*       d      - w.p. real (n)
c*                contains the eigenvalues in ascending order.  if an
c*                error exit is made, the eigenvalues are correct and
c*                ordered for indices 1,2,...ierr-1, but may not be
c*                the smallest eigenvalues.
c*       e2     - w.p. real (n)
c*                destroyed.
c*       ierr   - integer*4
c*                set to
c*                zero       for normal return,
c*                j          if the j-th eigenvalue has not been
c*                           determined after 30 iterations.
c*
c*    differences from eispack 3 -
c*       g=g+b instead of if(g.eq.0) g=b ; b=b/4
c*       f77 backward loops instead of f66 construct
c*       generic intrinsic functions
c*       arrary  ind  added for use by einvit
c*
c*    external routines -
c*       epslon
c*       intrinsic--abs, sign, sqrt
c*
c*    note -
c*       questions and comments concerning eispack should be directed to
c*       b. s. garbow, applied math. division, argonne national lab.
c*
c***********************************************************************
      subroutine eqlrat(n,diag,e,e2in,d,ind,ierr,e2)
c
      integer i,j,l,m,n,ii,l1,ierr
      integer ind(n)
c
      double precision d(n),e(n),e2(n),diag(n),e2in(n)
      double precision b,c,f,g,h,p,r,s,t,epslon
      double precision scale,zero,one
c
      parameter (zero = 0.0d+00, scale= 1.0d+00/64.0d+00, one = 1.0d+00)
c
c-----------------------------------------------------------------------
      ierr = 0
      d(1)=diag(1)
      ind(1) = 1
      k = 0
      itag = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
         d(i)=diag(i)
  100 e2(i-1) = e2in(i)
c
      f = zero
      t = zero
      b = epslon(one)
      c = b *b
      b = b * scale
      e2(n) = zero
c
      do 290 l = 1, n
         h = abs(d(l)) + abs(e(l))
         if (t .ge. h) go to 105
            t = h
            b = epslon(t)
            c = b * b
            b = b * scale
  105    continue
c     .......... look for small squared sub-diagonal element ..........
         m = l - 1
  110    m = m + 1
         if (e2(m) .gt. c) go to 110
c     .......... e2(n) is always zero, so there is an exit
c                from the loop ..........
c
         if (m .le. k) go to 125
            if (m .ne. n) e2in(m+1) = zero
            k = m
            itag = itag + 1
  125    continue
         if (m .eq. l) go to 210
c
c           iterate
c
         do 205 j = 1, 30
c              .......... form shift ..........
            l1 = l + 1
            s = sqrt(e2(l))
            g = d(l)
            p = (d(l1) - g) / (2.0d+00 * s)
            r = sqrt(p*p+1.0d+00)
            d(l) = s / (p + sign(r,p))
            h = g - d(l)
c
            do 140 i = l1, n
  140       d(i) = d(i) - h
c
            f = f + h
c              .......... rational ql transformation ..........
            g = d(m) + b
            h = g
            s = zero
            do 200 i = m-1,l,-1
               p = g * h
               r = p + e2(i)
               e2(i+1) = s * r
               s = e2(i) / r
               d(i+1) = h + s * (h + d(i))
               g = d(i) - e2(i) / g   + b
               h = g * p / r
  200       continue
c
            e2(l) = s * g
            d(l) = h
c              .......... guard against underflow in convergence test
            if (h .eq. zero) go to 210
            if (abs(e2(l)) .le. abs(c/h)) go to 210
            e2(l) = h * e2(l)
            if (e2(l) .eq. zero) go to 210
  205    continue
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
      ierr = l
      go to 1001
c
c           converged
c
  210    p = d(l) + f
c           .......... order eigenvalues ..........
         i = 1
         if (l .eq. 1) go to 250
            if (p .lt. d(1)) go to 230
               i = l
c           .......... loop to find ordered position
  220          i = i - 1
               if (p .lt. d(i)) go to 220
c
               i = i + 1
               if (i .eq. l) go to 250
  230       continue
            do 240 ii = l, i+1, -1
               d(ii) = d(ii-1)
               ind(ii) = ind(ii-1)
  240       continue
c
  250    continue
         d(i) = p
         ind(i) = itag
  290 continue
c
 1001 return
      end
