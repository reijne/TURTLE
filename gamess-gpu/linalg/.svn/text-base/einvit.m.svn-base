c***********************************************************************
c*       routine einvit(nm,n,d,e,e2,m,w,ind,z,ierr,rv1,rv2,rv3,rv4,rv6)
c*
c*    authors-
c*       this is a modification of routine tinvit from eispack edition 3
c*       dated august 1983.
c*       tinvit is a translation of the inverse iteration technique
c*       in the algol procedure tristurm by peters and wilkinson.
c*       handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c*       this version is by s. t. elbert (ames laboratory-usdoe)
c*
c*    purpose -
c*       this routine finds those eigenvectors of a tridiagonal
c*       symmetric matrix corresponding to specified eigenvalues.
c*
c*    method -
c*       inverse iteration.
c*
c*    on entry -
c*       nm     - integer*4
c*                must be set to the row dimension of two-dimensional
c*                array parameters as declared in the calling routine
c*                dimension statement.
c*       n      - integer*4
c*       d      - w.p. real (n)
c*                contains the diagonal elements of the input matrix.
c*       e      - w.p. real (n)
c*                contains the subdiagonal elements of the input matrix
c*                in its last n-1 positions.  e(1) is arbitrary.
c*       e2     - w.p. real (n)
c*                contains the squares of corresponding elements of e,
c*                with zeros corresponding to negligible elements of e.
c*                e(i) is considered negligible if it is not larger than
c*                the product of the relative machine precision and the
c*                sum of the magnitudes of d(i) and d(i-1).  e2(1) must
c*                contain 0.0 if the eigenvalues are in ascending order,
c*                or 2.0 if the eigenvalues are in descending order.
c*                if tqlrat, bisect, tridib, or imtqlv
c*                has been used to find the eigenvalues, their
c*                output e2 array is exactly what is expected here.
c*       m      - integer*4
c*                the number of specified eigenvectors.
c*       w      - w.p. real (m)
c*                contains the m eigenvalues in ascending
c*                or descending order.
c*       ind    - integer*4 (m)
c*                contains in first m positions the submatrix indices
c*                associated with the corresponding eigenvalues in w --
c*                1 for eigenvalues belonging to the first submatrix
c*                from the top, 2 for those belonging to the second
c*                submatrix, etc.
c*       ierr   - integer*4 (logical unit number)
c*                logical unit for error messages
c*
c*    on exit -
c*       all input arrays are unaltered.
c*       z      - w.p. real (nm,m)
c*                contains the associated set of orthonormal
c*                eigenvectors. any vector which which fails to converge
c*                is left as is (but normalized) when iterating stopped.
c*       ierr   - integer*4
c*                set to
c*                zero    for normal return,
c*                -r      if the eigenvector corresponding to the r-th
c*                        eigenvalue fails to converge in 5 iterations.
c*                        (only last failure to converge is reported)
c*
c*       rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c*
c*       rv1    - w.p. real (n)
c*                diagonal elements of u from lu decomposition
c*       rv2    - w.p. real (n)
c*                super(1)-diagonal elements of u from lu decomposition
c*       rv3    - w.p. real (n)
c*                super(2)-diagonal elements of u from lu decomposition
c*       rv4    - w.p. real (n)
c*                elements defining l in lu decomposition
c*       rv6    - w.p. real (n)
c*                approximate eigenvector
c*
c*    differences from eispack 3 -
c*       eps3 is scaled by  epscal  (enhances convergence, but
c*          lowers accuracy)!
c*       one more iteration (minimum 2) is performed after convergence
c*          (enhances accuracy)!
c*       replace loop with pythagG with single call to dnrm2!
c*       if not converged, use performance index to decide on error
c*          value setting, but do not stop!
c*       l.u. for error messages passed through ierr
c*       use parameter statements and generic intrinsic functions
c*       use level 1 blas
c*       use if-then-else to clarify logic
c*       loop over subspaces made into do loop.
c*       loop over inverse iterations made into do loop
c*       zero only required portions of output vector
c*
c*    external routines -
c*       epslon
c*       blas(1)--dasum, daxpy, ddot, dnrm2, dscal
c*       intrinsic functions - abs, max, sqrt
c*
c*    note -
c*       questions and comments concerning eispack should be directed to
c*       b. s. garbow, applied math. division, argonne national lab.
c*
c***********************************************************************
      subroutine einvit(nm,n,d,e,e2,m,w,ind,z,ierr,rv1,rv2,rv3,rv4,rv6)
c
c                        declarations
c
      logical convgd,goparr,dskwrk,maswrk
c
      integer group,i,ierr,its,j,jj,m,n,nm,p,q,r,s,submat,tag
      integer ind(m)
c
      double precision d(n),e(n+1),e2(n),w(m),z(nm,m)
      double precision rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision anorm,eps2,eps3,eps4,norm,order,rho,u,uk,v
      double precision x0,x1,xu
      double precision epscal,grptol,hundrd,one,ten,zero
      double precision epslon, estpi1
_IF(t3e,cray)
      REAL `sasum', `sdot', `snrm2'
_ELSE
      REAL dasum, ddot, dnrm2
_ENDIF
c
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      parameter (zero = 0.0d+00, one = 1.0d+00, grptol = 0.001d+00)
      parameter (epscal = 0.5d+00, hundrd = 100.0d+00, ten = 10.0d+00)
c
  001 format(' eigenvector routine einvit did not converge for vector'
     *      ,i5,'.  norm =',1p,e10.2,' performance index =',e10.2/
     *      ' (an error halt will occur if the pi is greater than 100)')
c
c-----------------------------------------------------------------------
c
      luemsg = ierr
      ierr = 0
      x0 = zero
      uk = zero
      norm = zero
      eps2 = zero
      eps3 = zero
      eps4 = zero
      group = 0
      tag = 0
      order = one - e2(1)
      q = 0
      do 930 submat = 1, n
         p = q + 1
c
c        .......... establish and process next submatrix ..........
c
         do 120 q = p, n-1
            if (e2(q+1) .eq. zero) go to 140
  120    continue
         q = n
c
c        .......... find vectors by inverse iteration ..........
c
  140    continue
         tag = tag + 1
         anorm = zero
         s = 0
c
         do 920 r = 1, m
            if (ind(r) .ne. tag) go to 920
            its = 1
            x1 = w(r)
            if (s .ne. 0) go to 510
c
c        .......... check for isolated root ..........
c
            xu = one
            if (p .eq. q) then
               rv6(p) = one
               convgd = .true.
               go to 860
c
            end if
            norm = abs(d(p))
            do 500 i = p+1, q
               norm = max( norm, abs(d(i)) + abs(e(i)) )
  500       continue
c
c        .......... eps2 is the criterion for grouping,
c                   eps3 replaces zero pivots and equal
c                   roots are modified by eps3,
c                   eps4 is taken very small to avoid overflow .........
c
            eps2 = grptol * norm
            eps3 = epscal * epslon(norm)
            uk = q - p + 1
            eps4 = uk * eps3
            uk = eps4 / sqrt(uk)
            s = p
            group = 0
            go to 520
c
c        .......... look for close or coincident roots ..........
c
  510       if (abs(x1-x0) .ge. eps2) then
c
c                 roots are seperate
c
               group = 0
            else
c
c                 roots are close
c
               group = group + 1
               if (order * (x1 - x0) .le. eps3) x1 = x0 + order * eps3
            end if
c
c        .......... elimination with interchanges and
c                   initialization of vector ..........
c
  520       continue
c
            u = d(p) - x1
            v = e(p+1)
            rv6(p) = uk
            do 550 i = p+1, q
               rv6(i) = uk
               if (abs(e(i)) .gt. abs(u)) then
c
c                 exchange rows before elimination
c
c                  *** warning -- a divide check may occur here if
c                      e2 array has not been specified correctly .......
c
                  xu = u / e(i)
                  rv4(i) = xu
                  rv1(i-1) = e(i)
                  rv2(i-1) = d(i) - x1
                  rv3(i-1) = e(i+1)
                  u = v - xu * rv2(i-1)
                  v = -xu * rv3(i-1)
c
               else
c
c                    straight elimination
c
                  xu = e(i) / u
                  rv4(i) = xu
                  rv1(i-1) = u
                  rv2(i-1) = v
                  rv3(i-1) = zero
                  u = d(i) - x1 - xu * v
                  v = e(i+1)
               end if
  550       continue
c
            if (abs(u) .le. eps3) u = eps3
            rv1(q) = u
            rv2(q) = zero
            rv3(q) = zero
c
c              do inverse iterations
c
            convgd = .false.
            do 800 its = 1, 5
               if (its .eq. 1) go to 600
c
c                    .......... forward substitution ..........
c
                  if (norm .eq. zero) then
                     rv6(s) = eps4
                     s = s + 1
                     if (s .gt. q) s = p
                  else
                     xu = eps4 / norm
                     call dscal(q-p+1, xu, rv6(p), 1)
                  end if
c
c                     ... elimination operations on next vector
c
                  do 590 i = p+1, q
                     u = rv6(i)
c
c                         if rv1(i-1) .eq. e(i), a row interchange
c                         was performed earlier in the
c                         triangularization process ..........
c
                     if (rv1(i-1) .eq. e(i)) then
                        u = rv6(i-1)
                        rv6(i-1) = rv6(i)
                     else
                        u = rv6(i)
                     end if
                     rv6(i) = u - rv4(i) * rv6(i-1)
  590             continue
  600          continue
c
c           .......... back substitution
c
               rv6(q) = rv6(q) / rv1(q)
               v = u
               u = rv6(q)
               norm = abs(u)
               do 620 i = q-1, p, -1
                  rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
                  v = u
                  u = rv6(i)
                  norm = norm + abs(u)
  620          continue
               if (group .eq. 0) go to 700
c
c                 ....... orthogonalize with respect to previous
c                         members of group ..........
c
                  j = r
                  do 680 jj = 1, group
  630                j = j - 1
                     if (ind(j) .ne. tag) go to 630
                     call daxpy(q-p+1, -ddot(q-p+1
     *                         ,rv6(p),1,z(p,j),1)
     *                         ,z(p,j),1,rv6(p),1)
  680             continue
                  norm = dasum(q-p+1, rv6(p), 1)
  700          continue
c
               if (convgd) go to 840
               if (norm .ge. one) convgd = .true.
  800       continue
c
c        .......... normalize so that sum of squares is
c                   1 and expand to full order ..........
c
  840       continue
c
            xu = one / dnrm2(q-p+1,rv6(p),1)
c
  860       continue
            do 870 i = 1, p-1
               z(i,r) = zero
  870       continue
            do 890 i = p,q
               z(i,r) = rv6(i) * xu
  890       continue
            do 900 i = q+1, n
               z(i,r) = zero
  900       continue
c
            if (.not.convgd) then
               rho = estpi1(q-p+1,x1,d(p),e(p),z(p,r),anorm)
               if (rho .ge. ten .and. luemsg .gt. 0 .and. maswrk)
     *             write(luemsg,001) r,norm,rho
c
c               *** set error -- non-converged eigenvector ..........
c
               if (rho .gt. hundrd) ierr = -r
            end if
c
            x0 = x1
  920    continue
c
         if (q .eq. n) go to 940
  930 continue
  940 continue
      return
      end
