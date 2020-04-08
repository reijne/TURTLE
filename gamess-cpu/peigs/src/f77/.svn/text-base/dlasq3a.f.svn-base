*
* $Id: dlasq3a.f,v 1.2 2000-10-26 15:38:33 psh Exp $
*
      SUBROUTINE DLASQ3A( I0, N0, Z, PP, DMIN, SIGMA, SIGTST,
     $                    EPS, NFAIL, ITER, NDIV, TTYPE )
*
*  -- LAPACK auxiliary routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1996
*
*     .. Scalar Arguments ..
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP, TTYPE
      DOUBLE PRECISION   DMIN, EPS, SIGMA, SIGTST
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
*     ..
*
*  Purpose
*  =======
*  DLASQ3A checks for deflation, computes a shift (TAU) and calls dqds.
*  In case of failure it changes shifts, and tries again until output
*  is positive.
*
*  Arguments
*  =========
*
*  I0     (input) INTEGER
*         First index.
*
*  N0     (input) INTEGER
*         Last index.
*
*  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
*         Z holds the qd array. 
*
*  PP     (input) INTEGER
*         PP=0 for ping, PP=1 for pong.
*
*  DMIN   (output) DOUBLE PRECISION
*         Minimum value of d.
*
*  SIGMA  (output) DOUBLE PRECISION
*         Sum of shifts used in current segment.
*
*  SIGTST (input) DOUBLE PRECISION
*         Sum of shifts used in current segment.
*
*  EPS    (input) DOUBLE PRECISION
*         Machine epsilon
*
*  NFAIL  (output) INTEGER
*         Number of times shift was too big.
*
*  ITER   (output) INTEGER
*         Number of iterations.
*
*  NDIV   (output) INTEGER
*         Number of divisions.
*
*  TTYPE  (output) INTEGER
*         Shift type.
*
*  Further Details
*  ===============
*  See "An implementation of dqds", LAPACK report XXX.
*  CWK66 = 1 - 1/sqrt(2)) (W. Kahan, 1966)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   CWK66, CBIAS
      PARAMETER          ( CWK66 = 0.292890D0, CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IPN4, J4, NN, N0IN
      DOUBLE PRECISION   DMIN1, DMIN2, DN, DN1, DN2, GAP2,
     $                   PROD, RHO2, S, SUM, T, TAU, TEMP
      SAVE               DMIN1, DMIN2, DN, DN1, DN2, TAU
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASQ4A, DLASQ5A
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Check for deflation. 
*
      N0IN = N0
*
   10 CONTINUE
*
      IF( N0.LT.I0 ) 
     $   GO TO 50
      IF( N0.EQ.I0 ) 
     $   GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.(I0+1) ) 
     $   GO TO 40
*
*     Check whether E(N0-1) is negligible, 1-by-1 case.
*
      IF( Z( NN-3 ).GT.Z( NN-7 ) .OR.
     $   Z( NN-5 ).GT.EPS*( SIGTST+Z( NN-7 ) ) ) 
     $   GO TO 30
*
      RHO2 = CWK66*Z( NN-7 )*Z( NN-9 )
      GAP2 = QURTR*( Z( NN-7 )-Z( NN-3 ) )**2
      PROD = Z( NN-3 )*Z( NN-5 )
      SUM  = RHO2 + GAP2
      IF( PROD*( TWO*RHO2*SUM + GAP2*PROD ).GT.
     $   ( EPS*( SIGTST+Z( NN-3 ) )*SUM )**2 ) 
     $   GO TO 30
*
   20 CONTINUE
*
      Z( 4*N0-3 ) = SIGMA + Z( 4*N0+PP-3 )
      N0 = N0 - 1
      GO TO 10
*
*     Check  whether E(N0-2) is negligible, 2-by-2 case.
*
   30 CONTINUE
*
      IF( Z( NN-9 ).GT.EPS*( SIGTST+Z( NN-11 ) ) .OR.
     $   Z( NN-9 )*Z( NN-7 ).GT.( SIGTST+Z( NN-11 ) )*
     $   ( SIGTST+Z( NN-7 ) )*EPS**2 ) 
     $   GO TO 50
*
   40 CONTINUE
*
      IF( Z( NN-3 ) .GT. Z( NN-7 ) ) THEN
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      END IF
      IF( Z( NN-7 ) .GT. Z( NN-5 ) ) THEN
         S = ( Z( NN-7 )-Z( NN-3 ) ) + Z( NN-5 )
      ELSE
         S = ( Z( NN-5 )-Z( NN-3 ) ) + Z( NN-7 )
      END IF
      S = TWO*Z( NN-7 )*Z( NN-3 ) / ( (Z( NN-7 )+Z( NN-5 )+Z( NN-3 )) +
     $    SQRT( S**2 + FOUR*Z( NN-3 )*Z( NN-5 ) ) )
      T = Z( NN-3 ) / ( ( Z( NN-7 )-S )+Z( NN-5 ) )
      Z( 4*N0-3 ) = ( Z( NN-7 )-S )*T + SIGMA
      Z( 4*N0-7 ) = Z( NN-7 ) + Z( NN-5 )*( ONE+T ) + SIGMA
      N0 = N0 - 2
      GO TO 10
*
   50 CONTINUE
*
      IF( N0.LT.I0 ) RETURN
*
*     Reverse the qd-array, if warranted.
*
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
         IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
            IPN4 = 4*(I0+N0)
            DO 60 J4 = 4*I0, 2*(I0+N0-1), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            DMIN = -ZERO
         ENDIF
      ENDIF
*
*     Choose a shift.
*
      IF( 1.LE.ITER .AND. ITER.LE.5 ) THEN
         TAU = ZERO
      ELSE
         CALL DLASQ4A( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2,
     $                 DN, DN1, DN2, TAU, TTYPE )
      END IF  
*
*     Call dqds until DMIN > 0.
*
   70 CONTINUE
*
      CALL DLASQ5A( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2,
     $              DN, DN1, DN2 )
*
      ITER = ITER + 1
      NDIV = NDIV + (N0-I0)
*
*     Check for NaN: "DMIN.NE.DMIN" 
*
      IF( DMIN.NE.DMIN ) THEN
         TAU = ZERO 
         GO TO 70
      END IF
*
*     Check for convergence hidden by negative DN
*
      IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.
     $   Z( 4*(N0-1)-PP ).LT.EPS*( SIGTST+DN1 ) .AND.
     $   ABS( DN ).LT.EPS*SIGTST ) THEN
         Z( 4*(N0-1)-PP+2 ) = ZERO
         DMIN = ABS( DMIN )
      END IF
*
      IF( DMIN.LT.ZERO ) THEN
*
*        Failure. Select new TAU and try again.
*
         NFAIL = NFAIL + 1
*
*        Failed twice. Play it safe.
*
         IF( TTYPE.LT.-22 ) THEN
	    TAU = ZERO 
	    GO TO 70
         END IF
*
         IF( DMIN1.GT.ZERO ) THEN
*
*           Late failure. Gives excellent shift.
*
            TAU = ( TAU+DMIN )*( ONE-TWO*EPS ) 
            TTYPE = TTYPE - 11 
         ELSE
*
*           Early failure. Divide by 4.
*
            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
         END IF
         GO TO 70
      END IF
*
      SIGMA = SIGMA + TAU
      SIGTST = SIGTST + TAU
      RETURN
*
*     End of DLASQ3A
*
      END
