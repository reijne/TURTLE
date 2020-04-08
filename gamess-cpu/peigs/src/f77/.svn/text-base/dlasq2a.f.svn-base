*
* $Id: dlasq2a.f,v 1.2 2000-10-26 15:38:33 psh Exp $
*
      SUBROUTINE DLASQ2A( N, Z, INFO )
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
      DOUBLE PRECISION   Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ2A computes all the eigenvalues of the symmetric positive 
*  definite tridiagonal matrix associated with the qd array Z 
*  to high relative accuracy.
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
*  N     (input) INTEGER
*        The number of rows and columns in the matrix. N >= 0.
*
*  Z     (workspace) DOUBLE PRECISION array, dimension ( 4*N )
*        On entry Z holds the qd array. On exit, entries 1 to N hold
*        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
*        trace, Z( 2*N+2 ) holds the sum of the eigenvalues, Z( 2*N+3 )
*        holds the iteration count, Z( 2*N+4 ) holds NDIVS/NIN^2, and
*        Z( 2*N+5 ) holds the percentage of shifts that failed.
*
*  INFO  (output) INTEGER
*        = 0: successful exit
*        < 0: if the i-th argument is a scalar and had an illegal
*             value, then INFO = -i, if the i-th argument is an
*             array and the j-entry had an illegal value, then
*             INFO = -(i*100+j)
*        > 0: the algorithm failed
*              =1, a split was marked by a positive value in E
*              =2, termination criterion of outer while loop not met 
*                  (program created more than N unreduced blocks)
*              =3, current block of Z not diagonalized after 30*N
*                  iterations (in inner while loop)
*              =4, POSITIVE TTYPE, WHAT ELSE SHOULD WE SAY?
*
*  Further Details
*  ===============
*  Local Variables: I0:N0 defines a current unreduced segment of Z.
*  The shifts are accumulated in SIGMA. Iteration count is in ITER.
*  Ping-pong is controlled by PP (alternates between 0 and 1).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR, HNDRD, TWO56
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                     FOUR = 4.0D0, HNDRD = 110.0D0,
     $                     TWO56 = 256.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I0, I4, ITER, IWHILA, IWHILB, K, N0, 
     $                   NBIG, NFAIL, NDIV, PP, SPLT, TTYPE
      DOUBLE PRECISION   D, DIAGMX, DMIN, E, EMAX, EPS, QMIN,
     $                   SCALE, SFMIN, SIGMA, SIGTST, TOL2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL, DLASQ3A, DLASRT, XERBLA
*     ..
*     .. External Functions
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*      
*     Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLASQ2A', -INFO )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         IF( Z( 1 ).LT.ZERO ) THEN
            INFO = -201
            CALL XERBLA( 'DLASQ2A', -INFO )
         ELSE
            Z( 2*N-1 ) = Z( 1 )
         END IF
         RETURN
      ELSE IF( N.EQ.2 ) THEN
         IF( Z( 2 ).LT.ZERO .OR. Z( 3 ).LT.ZERO ) THEN
            INFO = -2
            CALL XERBLA( 'DLASQ2A', -INFO )
            RETURN
         ELSE IF( Z( 2 ).GT.Z( 1 ) ) THEN
            D = Z( 2 )
            Z( 2 ) = Z( 1 )
            Z( 1 ) = D
         END IF
         IF( Z( 1 ).GT.Z( 3 ) ) THEN
            D = ( Z( 1 )-Z( 2 ) ) + Z( 3 )
         ELSE
            D = ( Z( 3 )-Z( 2 ) ) + Z( 1 )
         END IF
         D = TWO*Z( 1 )*Z( 2 ) / ( ( Z( 1 )+Z( 2 )+Z( 3 ) ) +
     $       SQRT( D**2 + FOUR*Z( 2 )*Z( 3 ) ) )
         E = Z( 2 ) / ( ( Z( 1 )-D ) + Z( 3 ) )
         Z( 2 ) = ( Z( 1 )-D )*E
         Z( 1 ) = Z( 1 ) + Z( 3 )*(ONE+E)
         RETURN
      END IF
*
*     Check for negative data and compute sums of q's and e's.
*      
      Z( 2*N ) = ZERO
      DIAGMX = ZERO
      D = ZERO
      E = ZERO
*
      DO 10 K = 1, N
         IF( Z( K ).LT.ZERO ) THEN
            INFO = -(200+K)
            CALL XERBLA( 'DLASQ2A', -INFO )
            RETURN
         ELSE IF( Z( N+K ).LT.ZERO ) THEN
            INFO = -(200+N+K)
            CALL XERBLA( 'DLASQ2A', -INFO )
            RETURN
         END IF
         D = D + Z( K )
         E = E + Z( N+K )
         DIAGMX = MAX( DIAGMX, Z( K )+Z( N+K ) )
   10 CONTINUE
*
*     Check for diagonality.
*
      IF( E.EQ.ZERO ) THEN
         CALL DLASRT( 'D', N, Z, INFO )
         Z( 2*N-1 ) = D
         RETURN
      END IF
*
      D = D + E
*
*     Check for zero data.
*
      IF( D.EQ.ZERO ) THEN
         Z( 2*N-1 ) = ZERO
         RETURN
      END IF
*         
*     Get machine parameters.
*
      EPS = DLAMCH( 'Precision' )
      SFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( ONE / ( TWO56*SFMIN ) )
      TOL2 = EPS**2
*         
*     Scale up using D (=TRACE).
*         
      CALL DLASCL( 'G', 0, 0, D, SCALE, 2*N-1, 1, Z,
     $     2*N-1, INFO )
*
*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
*
      DO 20 K = 2*N, 2, -2
         Z( 2*K ) = ZERO 
         Z( 2*K-1 ) = Z( K ) 
         Z( 2*K-2 ) = ZERO 
         Z( 2*K-3 ) = Z( K-1 ) 
   20 CONTINUE
*
*     Initialize SIGTST for desired level of accuracy.
*
*     SIGTST = DIAGMX : absolute accuracy
*     SIGTST = ZERO   : relative accuracy
*
      SIGTST = ZERO
*
      I0 = 1
      N0 = N
*
*     Check for negligible off-diagonal entries, geometric mean test.
*
      SPLT = I0 - 1
      DO 30 I4 = 4*I0, 4*(N0-1), 4
         IF( Z( I4-1 )*Z( I4+1 ).LE.
     $      TOL2*( Z( I4-3 )+Z( I4-1 ) )*( Z( I4+1 )+Z( I4+3 ) ) ) THEN
            Z( I4-1 ) = ZERO 
            SPLT = I4/4
         END IF
   30 CONTINUE
*
      ITER = 0
      NDIV = 0
      NFAIL = 0
*
      DO 90 IWHILA = 1, N+1
         IF( N0.LT.1 ) 
     $      GO TO 100
*
*        While array unfinished do 
*
*        E(N0) holds the value of SIGMA when submatrix in I0:N0
*        splits from the rest of the array, but is negated.
*      
         SIGMA = -Z( 4*N0-1 )
         IF( SIGMA.LT.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Find last unreduced submatrix's top index I0.
*        Find Gershgorin-type bound if Q's much greater than E's.
*
         EMAX = ZERO 
         QMIN = Z( 4*N0-3 )
         DO 40 I4 = 4*N0, 8, -4
            IF( Z( I4-5 ).LE.ZERO )
     $         GO TO 50
            IF( QMIN.GE.FOUR*EMAX ) THEN
               QMIN = MIN( QMIN, Z( I4-3 ) )
               EMAX = MAX( EMAX, Z( I4-5 ) )
            ENDIF
   40    CONTINUE
         I4 = 4 
*
   50    CONTINUE
         I0 = I4/4
*
*        Put -(initial shift) into dmin.
*
         DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN*EMAX ) )
*
*        Now I0:N0 is unreduced. PP = 0 for ping, PP = 1 for pong.
*
         PP = 0 
*
         NBIG = 30*(N0-I0+1)
         DO 70 IWHILB = 1, NBIG
            IF( I0.GT.N0 ) 
     $         GO TO 80
*
*           While last unreduced submatrix unfinished do 
*
            CALL DLASQ3A( I0, N0, Z, PP, DMIN, SIGMA, SIGTST,
     $                    EPS, NFAIL, ITER, NDIV, TTYPE )
*
            IF( TTYPE.GT.0 ) THEN
               INFO = 4
               RETURN
            END IF
	    PP = 1 - PP
*
*           When Z, not ZZ, holds the latest qd-array and EMIN is
*           very small then check for splits.
*
            IF( PP.EQ.0 .AND. (N0-I0).GE.2 ) THEN
               IF( Z( 4*N0 )*Z( 4*N0-11 ).LE.
     $            TOL2*( SIGTST+Z( 4*N0-11 ) )*
     $            ( SIGTST+Z( 4*N0-7 ) ) ) THEN
                  SPLT = I0 - 1
                  DO 60 I4 = 4*I0, 4*(N0-1), 4
                     IF( Z( I4-1 )*Z( I4+1 ).LE.
     $                  TOL2*( Z( I4-3 )+Z( I4-1 )+SIGTST )*
     $                  ( Z( I4+1 )+Z( I4+3 )+SIGTST ) ) THEN
                        Z( I4-1 ) = -SIGMA 
                        SPLT = I4/4
                     END IF
   60             CONTINUE
                  I0 = SPLT + 1
               END IF
            END IF
*
   70    CONTINUE
*
         INFO = 2
         RETURN
*
*        end IWHILB
*
   80    CONTINUE
*
   90 CONTINUE
*
      INFO = 3
      RETURN
*
*     end IWHILA   
*
  100 CONTINUE
*      
      DO 110 K = 2, N
         Z( K ) = Z( 4*K-3 )
  110 CONTINUE
*      
*     Undo scaling and sort. Compute sum of eigenvalues.
*
      CALL DLASCL( 'G', 0, 0, SCALE, D, N, 1, Z, N, INFO )
      CALL DLASRT( 'D', N, Z, INFO )
*
      E = ZERO
      DO 120 K = N, 1, -1
         E = E + Z( K )
  120 CONTINUE
*
*     Store trace, sum(eigenvalues) and information on performance.
*
      Z( 2*N+1 ) = D 
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = DBLE( ITER )
      Z( 2*N+4 ) = DBLE( NDIV )/DBLE( N**2 )
      Z( 2*N+5 ) = HNDRD*NFAIL/DBLE( ITER )
      RETURN
*
*     End of DLASQ2A
*
      END
