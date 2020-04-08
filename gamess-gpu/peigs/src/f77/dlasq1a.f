*
* $Id: dlasq1a.f,v 1.2 2000-10-26 15:38:33 psh Exp $
*
      SUBROUTINE DLASQ1A( N, D, E, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 15, 1996 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ1A computes all the eigenvalues of the symmetric positive 
*  definite tridiagonal matrix associated with the qd array Z 
*  (which is stored in WORK) to high relative accuracy.
*
*  To see the relation of Z to the tridiagonal matrix, let L be a
*  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
*  let U be an upper bidiagonal matrix with 1's above and diagonal
*  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
*  symmetric tridiagonal to which it is similar.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The number of rows and columns in the matrix. N >= 0.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, D contains the diagonal elements of the
*         tridiagonal matrix whose eigenvalues are desired.
*         On normal exit, D contains the eigenvalues in
*         decreasing order.
*
*  E      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, elements E(1:N-1) contain the off-diagonal
*         elements of the tridiagonal matrix whose eigenvalues
*         are desired.
*
*  WORK   (workspace) DOUBLE PRECISION array, dimension (4*N)
*         WORK holds the qd array. On exit, WORK(1) holds the trace
*         WORK(2) holds the sum of the eigenvalues, WORK(3) holds
*         NDIVS/NIN^2, WORK(4) holds the iteration count, and
*         WORK(5) holds the percentage of shifts that failed.
*
*  INFO   (output) INTEGER
*         = 0:  successful exit
*         < 0:  if INFO = -i, the i-th argument had an illegal value
*         > 0:  if INFO = i, the algorithm did not converge;  i
*               specifies how many superdiagonals did not converge.
*
*  INFO   (output) INTEGER
*         = 0 - successful exit
*         > 0 - indicates an error
*             =1, N<1 (prologue)
*             =2, bad data, E(I)<0 or Q(I)<=0 for some i (prologue)
*             =3, a split was marked by a positive value in E (geteigs)
*             =4, termination criterion of outer while loop not met.
*                 Program created more than N unreduced blocks
*                 (geteigs).
*             =5, current block of WORK not diagonalized after 10*N
*                 iterations (in inner while loop) (geteigs)
*
*  Further Details
*  ===============
*  Local Variables: I0:N0 defines a current unreduced segment of WORK.
*  The shifts are accumulated in SIGMA. Iteration count is in IT.
*  Ping-pong is controlled by PP (alternates between 0 and 1).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, FOUR, HNDRD
      PARAMETER          ( ZERO = 0.0D+0, TWO = 2.0D0, FOUR = 4.0D0,
     $                     HNDRD = 100.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            FAIL, I, I0, I4, IT, IWHILA, IWHILB,
     $                   J, NBIG, N0, NDIV, PP, SPLT 
      DOUBLE PRECISION   DIAGMX, DMIN, EMAX, EPS, QMIN,
     $                   SIGMA, SIGMN, SIGMX, SIGTST, TOL2, TRACE
*     ..
*     .. External Subroutines ..
      EXTERNAL           GOODSTEP
*     ..
*     .. External Functions
      DOUBLE PRECISION   DASUM, DLAMCH
      EXTERNAL           DASUM, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*      
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLASQ1', -INFO )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      ELSE IF( N.EQ.2 ) THEN
         CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      END IF
*
*     Compute EPS.
*      
      EPS = DLAMCH( 'Precision')
*
*     Convert q's and e's to Z. Interlace data for locality,
*     Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...)
*
      J = 3
      TRACE = D(1)
      WORK( 1 ) = D( 1 )
      WORK( 2 ) = ZERO
      DO 10 I = 2,N
         WORK( J   ) = E( I )
         WORK( J+1 ) = ZERO 
         WORK( J+2 ) = D( I+1 )
         WORK( J+3 ) = ZERO 
         J = J + 4
         TRACE = TRACE + E( I ) + D( I+1 )
   10 CONTINUE
      WORK( J   ) = ZERO 
      WORK( J+1 ) = ZERO 
*
*     Scale here, if desired
*
      IF( TRACE.EQ.ZERO ) WORK( 2*N+1 ) = ZERO 
*
      IF( INFO.GT.0 ) RETURN
      IF( TRACE.EQ.ZERO ) GO TO 110
*
      I0 = 1
      N0 = N
      TOL2 = EPS**2
*
*     Initialize SIGTST for desired level of accuracy.
*
*     SIGTST = ZERO   : relative accuracy
*     SIGTST = DIAGMX : absolute accuracy
*
      SIGTST = ZERO
*
*     Check for negligible off-diagonal entries, use geometric
*     mean test.
*
      SPLT = I0 - 1
      DO 20 I4 = 4*I0, 4*(N0-1), 4
         IF( WORK( I4-1 )*WORK( I4+1 ).LE.
     $      TOL2*(WORK( I4-3 )+WORK( I4-1 ))*
     $      (WORK( I4+1 )+WORK( I4+3 )) ) THEN
            WORK( I4-1 ) = ZERO 
            SPLT = I4/4
         END IF
   20 CONTINUE
*
      IT = 0
      NDIV = 0
      FAIL = 0
*
      DO 80 IWHILA = 1, N+1
         IF( N0.LT.1 ) GO TO 90
*
*        While array unfinished do 
*
*        E(N0) holds the value of SIGMA when submatrix in I0:N0
*        splits from the rest of the array, but is negated.
*      
         SIGMA = -WORK( 4*N0-1 )
         IF( SIGMA.LT.ZERO ) THEN
            INFO = 3
            RETURN
         END IF
*
*        Find last unreduced submatrix's top index I0.
*        Find Gersgorin-type bound if Q's much greater than E's.
*
         EMAX = ZERO 
         QMIN = WORK( 4*N0-3 )
         DO 30 I4 = 4*N0, 8, -4
            IF( WORK( I4-5 ).LE.ZERO ) GO TO 40
            IF( QMIN.GE.FOUR*EMAX ) THEN
               QMIN = MIN( QMIN, WORK( I4-3 ) )
               EMAX = MAX( EMAX, WORK( I4-5 ) )
            ENDIF
   30    CONTINUE
         I4 = 4 
*
   40    CONTINUE
         I0 = I4/4
*
*        Put - initial shift into dmin.
*
         DMIN = -MAX( ZERO, QMIN - TWO*SQRT(QMIN*EMAX) )
*
*        Now I0:N0 is unreduced.
*        PP = 0 for ping, PP = 1 for pong.
*
         PP = 0 
*
         NBIG = 30*(N0-I0+1)
         DO 60 IWHILB = 1, NBIG
            IF ( I0.GT.N0 ) GO TO 70
*
*           While last unreduced submatrix unfinished do 
*
            CALL GOODSTEP( I0, N0, WORK, PP, DMIN, FAIL, IT, NDIV,
     $                     SIGMA, SIGTST, EPS, INFO )
*
            IF( INFO.GT.0 ) RETURN
	    PP = 1 - PP
*
*           When Z, not ZZ, holds the latest qd-array and EMIN is
*           very small then check for splits.
*
            IF( PP.EQ.0 .AND. N0-I0.GE.2 ) THEN
               IF( WORK( 4*N0 )*WORK( 4*N0-11 ).LE.
     $            TOL2*(SIGTST+WORK( 4*N0-11 ))*
     $            (SIGTST+WORK( 4*N0-7 )) ) THEN
                  SPLT = I0 - 1
                  DO 50 I4 = 4*I0, 4*(N0-1), 4
                     IF( WORK( I4-1 )*WORK( I4+1 ).LE.
     $                  TOL2*(WORK( I4-3 )+WORK( I4-1 )+SIGTST)*
     $                  (WORK( I4+1 )+WORK( I4+3 )+SIGTST) ) THEN
                        WORK( I4-1 ) = -SIGMA 
                        SPLT = I4/4
                     END IF
   50             CONTINUE
                  I0 = SPLT + 1
               END IF
            END IF
*
   60    CONTINUE
*
         INFO = 5
         RETURN
*
*        end IWHILB
*
   70    CONTINUE
*
   80 CONTINUE
*
      INFO = 4
      RETURN
*
*     end IWHILA   
*
   90 CONTINUE
*      
      D( 1 ) = WORK( 1 )
      DO 100 I = 2, N
         D( I ) = WORK( 4*I-3 )
  100 CONTINUE
*
*     Sort and undo scaling, if needed. Compute sum of eigenvalues.
*
      CALL DLASRT( 'D', N, D, INFO )
*
      WORK( 2 ) = DASUM( N, D, 1 )
      WORK( 1 ) = TRACE 
*
*     Store information on performance.
*
      WORK( 3 ) = DBLE( NDIV )/DBLE( N**2 )
      WORK( 4 ) = DBLE( IT )
      WORK( 5 ) = HNDRD*FAIL/DBLE( IT )
  110 CONTINUE
      RETURN
*
*     End of DLASQ1A
*
      END
