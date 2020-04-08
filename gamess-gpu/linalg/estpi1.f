c*module eigen   *deck estpi1
c***********************************************************************
c*    funct. routine estpi1 (n,eval,d,e,x,anorm)
c*
c*    author -
c*       stephen t. elbert (ames laboratory-usdoe) date: 5 dec 1986
c*
c*    purpose -
c*       evaluate symmetric tridiagonal matrix performance index
c*       *        *         *                  *           *
c*       for 1 eigenvector
c*           *
c*
c*    method -
c*       this routine forms the 1-norm of the residual matrix a*x-x*eval
c*       where  a  is a symmetric tridiagonal matrix stored
c*       in the diagonal (d) and sub-diagonal (e) vectors, eval is the
c*       eigenvalue of an eigenvector of  a,  namely  x.
c*       this norm is scaled by machine accuracy for the problem size.
c*       all norms appearing in the comments below are 1-norms.
c*
c*    on entry -
c*       n      - integer*4
c*                the order of the matrix  a.
c*       eval   - w.p. real
c*                the eigenvalue corresponding to vector  x.
c*       d      - w.p. real (n)
c*                the diagonal vector of  a.
c*       e      - w.p. real (n)
c*                the sub-diagonal vector of  a.
c*       x      - w.p. real (n)
c*                an eigenvector of  a.
c*       anorm  - w.p. real
c*                the norm of  a  if it has been previously computed.
c*
c*    on exit -
c*       anorm  - w.p. real
c*                the norm of  a, computed if initially zero.
c*       estpi1 - w.p. real
c*          !!a*x-x*eval!! / (epslon(10*n)*!!a!!*!!x!!);
c*          where epslon(x) is the smallest number such that
c*             x + epslon(x) .ne. x
c*
c*          estpi1 .lt. 1 == satisfactory performance
c*                 .ge. 1 and .le. 100 == marginal performance
c*                 .gt. 100 == poor performance
c*          (see lect. notes in comp. sci. vol.6 pp 124-125)
c*
c***********************************************************************
      function estpi1 (n,eval,d,e,x,anorm)
c
c                    declarations
c
      double precision estpi1
      double precision anorm,eval,rnorm,size,xnorm
      double precision d(n), e(n), x(n)
      double precision epslon, one, zero
c
      parameter (zero = 0.0d+00, one = 1.0d+00)
c
c-----------------------------------------------------------------------
c
      estpi1 = zero
      if( n .le. 1 ) return
      size = 10 * n
      if (anorm .eq. zero) then
c
c              compute norm of  a
c
         anorm = max( abs(d(1)) + abs(e(2))
     *               ,abs(d(n)) + abs(e(n)))
         do 110 i = 2, n-1
            anorm = max( anorm, abs(e(i))+abs(d(i))+abs(e(i+1)))
  110    continue
         if(anorm .eq. zero) anorm = one
      end if
c
c           compute norms of residual and eigenvector
c
      xnorm = abs(x(1)) + abs(x(n))
      rnorm = abs( (d(1)-eval)*x(1) + e(2)*x(2))
     *       +abs( (d(n)-eval)*x(n) + e(n)*x(n-1))
      do 120 i = 2, n-1
         xnorm = xnorm + abs(x(i))
         rnorm = rnorm + abs(e(i)*x(i-1) + (d(i)-eval)*x(i)
     *                       + e(i+1)*x(i+1))
  120 continue
c
      estpi1 = rnorm / (epslon(int(size))*anorm*xnorm)
      return
      end
