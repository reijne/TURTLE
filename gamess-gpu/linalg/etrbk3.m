c*module eigen   *deck etrbk3
c***********************************************************************
c*       routine etrbk3(nm,n,nv,a,m,z)
c*
c*    authors-
c*       this is a modification of routine trbak3 from eispack edition 3
c*       dated august 1983.
c*       eispack trbak3 is a translation of the algol procedure trbak3,
c*       num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c*       handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c*       this version is by s. t. elbert (ames laboratory-usdoe)
c*
c*    purpose -
c*       this routine forms the eigenvectors of a real symmetric
c*       matrix by back transforming those of the corresponding
c*       symmetric tridiagonal matrix determined by  etred3.
c*
c*    method -
c*       the calculation is carried out by forming the matrix product
c*          q*z
c*       where  q  is a product of the orthogonal symmetric matrices
c*                q = prod(i)[1 - u(i)*.transpose.u(i)*h(i)]
c*       u  is the augmented sub-diagonal rows of  a  and
c*       z  is the set of eigenvectors of the tridiagonal
c*       matrix  f  which was formed from the original symmetric
c*       matrix  c  by the similarity transformation
c*                f = q(transpose) c q
c*       note that etrbk3 preserves vector euclidean norms.
c*
c*
c*    complexity -
c*       m*n**2
c*
c*    on entry-
c*       nm     - integer*4
c*                must be set to the row dimension of two-dimensional
c*                array parameters as declared in the calling routine
c*                dimension statement.
c*       n      - integer*4
c*                the order of the matrix  a.
c*       nv     - integer*4
c*                must be set to the dimension of the array  a  as
c*                declared in the calling routine dimension statement.
c*       a      - w.p. real (nv)
c*                contains information about the orthogonal
c*                transformations used in the reduction by  etred3  in
c*                its first  nv = n*(n+1)/2 positions.
c*       m      - integer*4
c*                the number of eigenvectors to be back transformed.
c*       z      - w.p real (nm,m)
c*                contains the eigenvectors to be back transformed
c*                in its first m columns.
c*
c*    on exit-
c*       z      - w.p. real (nm,m)
c*                contains the transformed eigenvectors
c*                in its first m columns.
c*
c*    differences with eispack 3 -
c*       the two inner loops are replaced by ddot and daxpy.
c*       multiplication used instead of division to find s.
c*       outer loop range changed from 2,n to 3,n.
c*       address pointers for  a  simplified.
c*
c*    external routines -
c*       blas(1)--ddot, daxpy
c*       intrinsic functions - none
c*
c*    note -
c*       questions and comments concerning eispack should be directed to
c*       b. s. garbow, applied math. division, argonne national lab.
c*
c***********************************************************************
      subroutine etrbk3(nm,n,nv,a,m,z)
c
c                        declarations
c
      integer i,ii,im1,iz,j,m,n,nm,nv
c
      double precision a(nv),z(nm,m)
      double precision h,s,zero
_IF(cray,t3e)
      REAL `sdot'
_ELSE
      REAL ddot
_ENDIF
c
      parameter (zero = 0.0d+00)
c
c-----------------------------------------------------------------------
c
      if (m .eq. 0) return
      if (n .le. 2) return
c
      ii=3
      do 140 i = 3, n
         iz=ii+1
         ii=ii+i
         h = a(ii)
         if (h .eq. zero) go to 140
            im1 = i - 1
            do 130 j = 1, m
               s = -( ddot(im1,a(iz),1,z(1,j),1) * h) * h
               call daxpy(im1,s,a(iz),1,z(1,j),1)
  130       continue
  140 continue
      return
      end
