      subroutine evvrsp (msgfl,n,nvect,lena,nv,a,b,ind,root,vect,iorder
     *                 ,ierr)
c***********************************************************************
c*       routine evvrsp (msgfl,n,nvect,lena,nv,a,b,ind,root,vect,iorder
c*   *                 ,ierr)
c*
c*    author:  s. t. elbert, ames laboratory-usdoe, june 1985
c*
c*    purpose -
c*       finds   (all) eigenvalues    and    (some or all) eigenvectors
c*                     *    *                                   *
c*       of a real symmetric packed matrix.
c*            *    *         *
c*
c*    method -
c*       the method as presented in this routine consists of four steps:
c*       first, the input matrix is reduced to tridiagonal form by the
c*       householder technique (orthogonal similarity transformations).
c*       second, the roots are located using the rational ql method.
c*       third, the vectors of the tridiagonal form are evaluated by the
c*       inverse iteration technique.  vectors for degenerate or near-
c*       degenerate roots are forced to be orthogonal.
c*       fourth, the tridiagonal vectors are rotated to vectors of the
c*       original array.
c*
c*       these routines are modifications of the eispack 3
c*       routines tred3, tqlrat, tinvit and trbak3
c*
c*       for further details, see eispack users guide, b. t. smith
c*       et al, springer-verlag, lecture notes in computer science,
c*       vol. 6, 2-nd edition, 1976.  another good reference is
c*       the symmetric eigenvalue problem by b. n. parlett
c*       published by prentice-hall, inc., englewood cliffs, n.j. (1980)
c*
c*    on entry -
c*       msgfl  - integer*4 (logical unit no.)
c*                file where error messages will be printed.
c*                if msgfl is 0, error messages will be printed on lu 6.
c*                if msgfl is negative, no error messages printed.
c*       n      - integer*4
c*                order of matrix a.
c*       nvect  - integer*4
c*                number of vectors desired.  0 .le. nvect .le. n.
c*       lena   - integer*4
c*                dimension of  a  in calling routine.  must not be less
c*                than (n*n+n)/2.
c*       nv     - integer*4
c*                row dimension of vect in calling routine.   n .le. nv.
c*       a      - working precision real (lena)
c*                input matrix, rows of the lower triangle packed into
c*                linear array of dimension n*(n+1)/2.  the packed order
c*                is a(1,1), a(2,1), a(2,2), a(3,1), a(3,2), ...
c*       b      - working precision real (n,8)
c*                scratch array, 8*n elements
c*       ind    - integer*4 (n)
c*                scratch array of length n.
c*       iorder - integer*4
c*                root ordering flag.
c*                = 0, roots will be put in ascending order.
c*                = 2, roots will be put in descending order.
c*
c*    on exit -
c*       a      - destoryed.  now holds reflection operators.
c*       root   - working precision real (n)
c*                all eigenvalues in ascending or descending order.
c*                  if iorder = 0, root(1) .le. ... .le. root(n)
c*                  if iorder = 2, root(1) .ge. ... .ge. root(n)
c*       vect   - working precision real (nv,nvect)
c*                eigenvectors for root(1), ..., root(nvect).
c*       ierr   - integer*4
c*                = 0 if no error detected,
c*                = k if iteration for k-th eigenvalue failed,
c*                = -k if iteration for k-th eigenvector failed.
c*                (failures should be very rare.  contact c. moler.)
c*
c*    external routines -
c*       etred3(elau, freda, dasum, dnrm2, dscal)
c*       eqlrat(epslon)
c*       einvit(epslon, estpi1, dasum, daxpy, ddot, dnrm2, dscal)
c*       etrbk3(ddot, daxpy)
c*
c***********************************************************************
c
c                              declarations
c
      logical goparr,dskwrk,maswrk
c
      double precision a(lena)
      double precision b(n,8)
      double precision root(n)
      double precision t
      double precision vect(nv,*)
c
      integer ind(n)
c
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
  900 format(26h0*** evvrsp parameters ***/
     +       14h ***      n = ,i8,4h ***/
     +       14h ***  nvect = ,i8,4h ***/
     +       14h ***   lena = ,i8,4h ***/
     +       14h ***     nv = ,i8,4h ***/
     +       14h *** iorder = ,i8,4h ***/
     +       14h ***   ierr = ,i8,4h ***)
  901 format(37h value of lena is less than (n*n+n)/2)
  902 format(39h eqlrat has failed to converge for root,i5)
  903 format(18h nv is less than n)
  904 format(41h einvit has failed to converge for vector,i5)
  905 format(51h value of iorder must be 0 (smallest root first) or
     *      ,23h 2 (largest root first))
  906 format(' value of n is less than or equal zero')
c
c-----------------------------------------------------------------------
c
      lmsgfl=msgfl
      if (msgfl .eq. 0) lmsgfl=6
      ierr = n - 1
      if (n .le. 0) go to 800
      ierr = n + 1
      if ( (n*n+n)/2 .gt. lena) go to 810
c
c        reduce real symmetric matrix a to tridiagonal form
c
      call etred3(n,lena,a,b(1,1),b(1,2),b(1,3))
c
c        find all eigenvalues of tridiagonal matrix
c
      call eqlrat(n,b(1,1),b(1,2),b(1,3),root,ind,ierr,b(1,4))
      if (ierr .ne. 0) go to 820
c
c         check the desired order of the eigenvalues
c
      b(1,3) = float(iorder)
      if (iorder .eq. 0) go to 300
         if (iorder .ne. 2) go to 850
c
c         order roots in descending order (largest first)...
c        turn root and ind arrays end for end
c
         do 210 i = 1, n/2
            j = n+1-i
            t = root(i)
            root(i) = root(j)
            root(j) = t
            l = ind(i)
            ind(i) = ind(j)
            ind(j) = l
  210    continue
c
c           find i and j marking the start and end of a sequence
c           of degenerate roots
c
         i=0
  220    continue
            i = i+1
            if (i .gt. n) go to 300
            do 230 j=i,n
               if (root(j) .ne. root(i)) go to 240
  230       continue
            j = n+1
  240       continue
            j = j-1
            if (j .eq. i) go to 220
c
c                    turn around ind between i and j
c
            jsv = j
            klim = (j-i+1)/2
            do 250 k=1,klim
               l = ind(j)
               ind(j) = ind(i)
               ind(i) = l
               i = i+1
               j = j-1
  250       continue
            i = jsv
         go to 220
c
  300 continue
c
      if (nvect .le. 0) return
      if (nv .lt. n) go to 830
c
c        find eigenvectors of tri-diagonal matrix via inverse iteration
c
      ierr = lmsgfl
      call einvit(nv,n,b(1,1),b(1,2),b(1,3),nvect,root,ind,
     +            vect,ierr,b(1,4),b(1,5),b(1,6),b(1,7),b(1,8))
      if (ierr .ne. 0) go to 840
c
c        find eigenvectors of symmetric matrix via back transformation
c
  400 continue
      call etrbk3(nv,n,lena,a,nvect,vect)
      return
c
c        error message section
c
  800 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,906)
      go to 890
c
  810 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,901)
      go to 890
c
  820 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,902) ierr
      go to 890
c
  830 if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,903)
      go to 890
c
  840 continue
      if ((lmsgfl .gt. 0).and.(maswrk)) write(lmsgfl,904) -ierr
      go to 400
c
  850 ierr=-1
      if (lmsgfl .lt. 0) return
      if (maswrk) write(lmsgfl,905)
      go to 890
c
  890 continue
      if (maswrk) write(lmsgfl,900) n,nvect,lena,nv,iorder,ierr
      return
      end
