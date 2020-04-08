c*module eigen   *deck etred3
c***********************************************************************
c*       routine etred3(n,nv,a,d,e,e2)
c*
c*    authors -
c*       this is a modification of routine tred3 from eispack edition 3
c*       dated august 1983.
c*       eispack tred3 is a translation of the algol procedure tred3,
c*       num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c*       handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c*       this version is by s. t. elbert, ames laboratory-usdoe jun 1986
c*
c*    purpose -
c*       this routine reduces a real symmetric (packed) matrix, stored
c*       as a one-dimensional array, to a symmetric tridiagonal matrix
c*       using orthogonal similarity transformations, preserving the
c*       information about the transformations in  a.
c*
c*    method -
c*       the tridiagonal reduction is performed in the following way.
c*       starting with j=n, the elements in the j-th row to the
c*       left of the diagonal are first scaled, to avoid possible
c*       underflow in the transformation that might result in severe
c*       departure from orthogonality.  the sum of squares  sigma  of
c*       these scaled elements is next formed.  then, a vector  u  and
c*       a scalar
c*                      h = u(transpose) * u / 2
c*       define a reflection operator
c*                      p = i - u * u(transpose) / h
c*       which is orthogonal and symmetric and for which the
c*       similiarity transformation  pap  eliminates the elements in
c*       the j-th row of  a  to the left of the subdiagonal and the
c*       symmetrical elements in the j-th column.
c*
c*       the non-zero components of  u  are the elements of the j-th
c*       row to the left of the diagonal with the last of them
c*       augmented by the square root of  sigma  prefixed by the sign
c*       of the subdiagonal element.  by storing the transformed sub-
c*       diagonal element in  e(j)  and not overwriting the row
c*       elements eliminated in the transformation, full information
c*       about  p  is save for later use in  etrbk3.
c*
c*       the transformation sets  e2(j)  equal to  sigma  and  e(j)
c*       equal to the square root of  sigma  prefixed by the sign
c*       of the replaced subdiagonal element.
c*
c*       the above steps are repeated on further rows of the
c*       transformed  a  in reverse order until  a  is reduced to tri-
c*       diagonal form, that is, repeated for  j = n-1,n-2,...,3.
c*
c*    complexity -
c*       2/3 n**3
c*
c*    on entry-
c*       n      - integer*4
c*                the order of the matrix.
c*       nv     - integer*4
c*                must be set to the dimension of the array parameter a
c*                as declared in the calling routine dimension statement
c*       a      - w.p. real (nv)
c*                contains the lower triangle of the real symmetric
c*                input matrix, stored row-wise as a one-dimensional
c*                array, in its first n*(n+1)/2 positions.
c*
c*    on exit-
c*       a      - w.p. real (nv)
c*                contains information about the orthogonal
c*                transformations used in the reduction.
c*       d      - w.p. real (n)
c*                contains the diagonal elements of the tridiagonal
c*                matrix.
c*       e      - w.p. real (n)
c*                contains the subdiagonal elements of the tridiagonal
c*                matrix in its last n-1 positions.  e(1) is set to zero
c*       e2     - w.p. real (n)
c*                contains the squares of the corresponding elements of
c*                e. may coincide with e if the squares are not needed.
c*
c*    differences from eispack 3 -
c*       outer loop changed from ii=1,n to i=n,3,-1
c*       parameter statement and generic intrinsic functions used
c*       scale.ne.0 test now spots tri-diagonal form
c*       values less than epslon cleared to zero
c*       use blas(1)
c*       u not copied to d, left in a
c*       e2 computed from e
c*       inner loops split into routines elau and freda
c*       inverse of h stored instead of h
c*
c*    external routines -
c*       elau, freda
c*       blas(1)--dasum, dnrm2, dscal
c*       intrinsic--abs, max, sign, sqrt
c*
c*    note -
c*       questions and comments concerning eispack should be directed to
c*       b. s. garbow, applied math. division, argonne national lab.
c*
c***********************************************************************
      subroutine etred3(n,nv,a,d,e,e2)
c
c                        declarations
c
      integer i,iia,iz0,l,n,nv
c
      double precision a(*),d(*),e(*),e2(*)
      double precision aiimax,f,g,h,hroot,scale,scalei
      real*8 dasum, dnrm2
      double precision one, zero
c
      parameter (zero = 0.0d+00, one = 1.0d+00)
c
c-----------------------------------------------------------------------
c
      if (n .le. 2) go to 310
      iz0 = (n*n+n)/2
      aiimax = abs(a(iz0))
      do 300 i = n, 3, -1
         l = i - 1
         iia = iz0
         iz0 = iz0 - i
         aiimax = max(aiimax, abs(a(iia)))
         scale = dasum(l, a(iz0+1), 1)
         if(scale .eq. abs(a(iia-1)) .or. aiimax+scale .eq. aiimax) then
c
c           this row is already in tri-diagonal form
c
            d(i) = a(iia)
            if (aiimax+d(i) .eq. aiimax) d(i) = zero
            e(i) = a(iia-1)
            if (aiimax+e(i) .eq. aiimax) e(i) = zero
            e2(i) = e(i)*e(i)
            a(iia) = zero
            go to 300
c
         end if
c
         scalei = one / scale
         call dscal(l,scalei,a(iz0+1),1)
         hroot = dnrm2(l,a(iz0+1),1)
c
         f = a(iz0+l)
         g = -sign(hroot,f)
         e(i) = scale * g
         e2(i) = e(i)*e(i)
         h = hroot*hroot - f * g
         a(iz0+l) = f - g
         d(i) = a(iia)
         a(iia) = one / sqrt(h)
c           .......... form p then q in e(1:l) ..........
         call elau(one/h,l,a(iz0+1),a,e)
c           .......... form reduced a ..........
         call freda(l,a(iz0+1),a,e)
c
  300 continue
  310 continue
      e(1) = zero
      e2(1)= zero
      d(1) = a(1)
      e(2) = a(2)
      e2(2)= a(2)*a(2)
      d(2) = a(3)
c
      return
      end
