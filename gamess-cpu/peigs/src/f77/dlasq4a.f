*
* $Id: dlasq4a.f,v 1.2 2000-10-26 15:38:33 psh Exp $
*
      SUBROUTINE DLASQ4A( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2,
     $                    DN, DN1, DN2, TAU, TTYPE )
*
*  -- LAPACK auxiliary routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1996
*
*     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
*     ..
*
*  Purpose
*  =======
*  DLASQ4A computes an approximation TAU to the smallest eigenvalue 
*  using values of d from the previous transform.
*
*  I0    (input) INTEGER
*        First index.
*
*  N0    (input) INTEGER
*        Last index.
*
*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
*        Z holds the qd array.
*
*  PP    (input) INTEGER
*        PP=0 for ping, PP=1 for pong.
*
*  NOIN  (output) INTEGER
*        The value of N0 at start of EIGTEST.
*
*  DMIN  (output) DOUBLE PRECISION
*        Minimum value of d.
*
*  DMIN1 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ).
*
*  DMIN2 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
*
*  DN    (output) DOUBLE PRECISION
*        d(N)
*
*  DN1   (output) DOUBLE PRECISION
*        d(N-1)
*
*  DN2   (output) DOUBLE PRECISION
*        d(N-2)
*
*  TAU   (input) DOUBLE PRECISION
*        This is the shift.
*
*  TTYPE (output) INTEGER
*        Shift type.
*
*  Further Details
*  ===============
*  See Section 7.3 of "An implementation of dqds", LAPACK report XXX.
*  CNST1 = 9/16
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0,
     $                     CNST3 = 1.050D0 )
      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0,
     $                     HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, HNDRD = 100.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I4, NN, NP
      DOUBLE PRECISION   A2, B1, B2, G, GAM, GAP1, GAP2, S, SN1
      SAVE               G
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     A negative DMIN forces the shift to take that absolute value
*     TTYPE records the type of shift.
*
      IF( DMIN .LE. ZERO ) THEN
         TAU = -DMIN
         TTYPE = -1
         RETURN
      END IF
*       
      NN = 4*N0 + PP
      IF( N0IN .EQ. N0 ) THEN
*
*        No eigenvalues deflated.
*
         IF( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) THEN
*
            B1 = Z( NN-3 ) * Z( NN-5 )
            B2 = Z( NN-7 ) * Z( NN-9 )
            A2 = Z( NN-7 ) + Z( NN-5 )
*
*           Cases 2 and 3.
*
	    IF( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) THEN
	       GAP2 = DMIN2 - A2 - DMIN2*QURTR
	       IF( GAP2.GT.ZERO .AND. GAP2**2.GT.B2 ) THEN
	          GAP1 = A2 - DN - B2/GAP2
	       ELSE
	          GAP1 = A2 - DN - SQRT( B1+B2 )
	       END IF
               IF( GAP1.GT.ZERO .AND. GAP1**2.GT.B1 ) THEN
	          S = MAX( DN-B1/GAP1, HALF*DMIN )
	          TTYPE = -2
	       ELSE
                  S = ZERO
                  IF( DN**2.GT.B1 )
     $               S = DN-SQRT( B1 )
                  IF( A2**2.GT.(B1+B2) )
     $               S = MIN( S, A2-SQRT( B1+B2 ) )
	          S = MAX( S, THIRD*DMIN )
	          TTYPE = -3
	       END IF
            ELSE
*
*              Case 4.
*
               IF( DMIN .EQ. DN ) THEN
                  GAM = DN
                  A2 = ZERO
                  B2 = Z( NN-5 )/Z( NN-7 )
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  B2 = Z( NP-2 )-TAU
                  GAM = DN1 - TAU*Z( NP-4 )/B2
                  A2 = Z( NP-2 )*Z( NP-4 )/B2**2
                  B2 = Z( NN-9 )/Z( NN-11 )
                  NP = NN - 13
               END IF
*
*              Approximate contribution to norm squared from I < NN-1.
*
               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0-1+PP, -4
                  B1 = B2
                  B2 = B2*Z( I4 )/Z( I4-2 )
                  A2 = A2 + B2
                  IF( HNDRD*MAX( B2, B1 ).LT.A2 .OR.
     $               CNST1.LT.A2 ) 
     $               GO TO 20
   10          CONTINUE
   20          CONTINUE
               A2 = CNST3*A2
*
*              Rayleigh quotient residual bound.
*
               IF( A2 .LT. CNST1 ) THEN
                  S = GAM*(ONE-SQRT(A2))/(ONE+A2)
               ELSE
                  S = QURTR*GAM
               END IF
               TTYPE = -4
	    END IF
         ELSE IF( DMIN .EQ. DN2 ) THEN
*
*           Case 5.
*
*           Compute contribution to norm squared from I > NN-2.
*
            NP = NN - 2*PP
            B1 = Z( NP-2 ) - TAU
            SN1 = -TAU*(ONE+Z( NP-4 )/B1)
            B2 = Z( NP-6 ) + SN1
            GAM = DN2 + SN1*Z( NP-8 )/B2
            A2 = (Z( NP-8 )*Z( NP-6 )/B2**2)*
     $            (ONE+Z( NP-4 )*Z( NP-2 )/B1**2)
*
*           Approximate contribution to norm squared from I < NN-2.
*
            IF( N0-I0.GT.2 ) THEN
               B2 = Z( NN-13 )/Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN-17, 4*I0-1+PP, -4
                  B1 = B2
                  B2 = B2*Z( I4 )/Z( I4-2 )
                  A2 = A2 + B2
                  IF( HNDRD*MAX( B2, B1 ).LT.A2 .OR.
     $               CNST1.LT.A2 ) 
     $               GO TO 40
   30          CONTINUE
   40          CONTINUE
               A2 = CNST3*A2
            END IF
*
            IF( A2.LT.CNST1 ) THEN
               S = GAM*(ONE-SQRT(A2))/(ONE+A2)
            ELSE
               S = QURTR*GAM/(ONE+A2)
            END IF
            TTYPE = -5
         ELSE
*
*           Case 6, no information to guide us.
*
            IF( TTYPE .EQ. -6 ) THEN
               G = G + THIRD*(ONE-G)
            ELSE IF ( TTYPE .EQ. -18 ) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            ENDIF
            S = G*DMIN
            TTYPE = -6
	 END IF
*
      ELSE IF( N0IN .EQ. (N0+1) ) THEN
*
*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
*
         IF( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) THEN 
*
*           Cases 7 and 8.
*
            B1 = Z( NN-5 )/Z( NN-7 )
	    B2 = B1
	    DO 50 I4 = 4*N0-9+PP, 4*I0-1+PP, -4
               A2 = B1
	       B1 = B1*Z( I4 )/Z( I4-2 )
	       B2 = B2 + B1
               IF( HNDRD*MAX( B1, A2 ).LT.B2 ) 
     $            GO TO 60
   50       CONTINUE
   60       CONTINUE
   	    B2 = CNST3*B2
	    A2 = DMIN1/(ONE + B2)
	    GAP2 = HALF*DMIN2 - A2
	    IF( GAP2.GT.ZERO .AND. GAP2**2.GT.B2*(A2**2)) THEN
               S = MAX( A2*(ONE-CNST2*A2*B2/GAP2), THIRD*DMIN1 )
               TTYPE = -7
	    ELSE 
               S = MAX( A2*(ONE-CNST2*SQRT(B2)), THIRD*DMIN1 )
               TTYPE = -8
            END IF
         ELSE
*
*           Case 9.
*
            S = QURTR*DMIN1
            IF( DMIN1 .EQ. DN1 )
     $         S = HALF*DMIN1
            TTYPE = -9
         END IF
*
      ELSE IF( N0IN .EQ. (N0+2) ) THEN
*
*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
*
*        Cases 10 and 11.
*
         IF( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) THEN 
            B1 = Z( NN-5 )/Z( NN-7 )
            B2 = B1
            DO 70 I4 = 4*N0-9+PP, 4*I0-1+PP, -4
               B1 = B1*Z( I4 )/Z( I4-2 )
               B2 = B2 + B1
               IF( HNDRD*B1.LT.B2 )
     $            GO TO 80
   70       CONTINUE
   80       CONTINUE
            B2 = CNST3*B2
            A2 = DMIN2/(ONE + B2)
            GAP2 = Z( NN-7 )+Z( NN-9 )-SQRT( Z( NN-11 )*Z( NN-9 ) )-A2
            IF( GAP2.GT.ZERO .AND. GAP2**2.GT.B2*(A2**2) ) THEN
	       S = MAX( A2*(ONE - CNST2*A2*B2/GAP2), THIRD*DMIN2 )
            ELSE 
               S = MAX( A2*(ONE - CNST2*SQRT(B2)), THIRD*DMIN2 )
            END IF
            TTYPE = -10
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         END IF
      ELSE IF( N0IN.GT.(N0+2) ) THEN
*
*        Case 12, more than two eigenvalues deflated. No information.
*
         S = ZERO 
         TTYPE = -12
      END IF
*
      TAU = S
      RETURN
*
*     End of DLASQ4A
*
      END
