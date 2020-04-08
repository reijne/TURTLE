*
* $Id: dlasq5a.f,v 1.2 2000-10-26 15:38:33 psh Exp $
*
      SUBROUTINE DLASQ5A( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, 
     $                    DN, DNM1, DNM2 )
*
*  -- LAPACK auxiliary routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1996
*
*     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
*     ..
*
*  Purpose
*  =======
*  DLASQ5A computes one dqds transform in ping-pong form.
*
*  Arguments
*  =========
*
*  I0    (input) INTEGER
*        First index.
*
*  N0    (input) INTEGER
*        Last index.
*
*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
*        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
*        an extra argument.
*
*  PP    (input) INTEGER
*        PP=0 for ping, PP=1 for pong.
*
*  TAU   (input) DOUBLE PRECISION
*        This is the shift.
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
*        d(N0), the last value of d.
*
*  DNM1  (output) DOUBLE PRECISION
*        d(N0-1).
*
*  DNM2  (output) DOUBLE PRECISION
*        d(N0-2).
*
*  Further Details
*  ===============
*  See "An implementation of dqds", LAPACK report XXX.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( (N0-I0-1) .LE. 0 ) RETURN
*
      J4 = 4*I0 + PP - 3
*
      EMIN = Z( J4+4 ) 
      D = Z( J4 ) - TAU
      DMIN = D
*
      IF( PP.EQ.0 ) THEN
         DO 10 J4 = 4*I0, 4*(N0-3), 4
            Z( J4-2 ) = D + Z( J4-1 ) 
            TEMP = Z( J4+1 )/Z( J4-2 )
            D = D * TEMP  - TAU
            DMIN = MIN( DMIN, D )
            Z( J4 ) = Z( J4-1 ) * TEMP
            EMIN = MIN( Z( J4 ), EMIN )
   10    CONTINUE
      ELSE
         DO 20 J4 = 4*I0, 4*(N0-3), 4
            Z( J4-3 ) = D + Z( J4 ) 
            TEMP = Z( J4+2 )/Z( J4-3 )
            D = D * TEMP  - TAU
            DMIN = MIN( DMIN, D )
            Z( J4-1 ) = Z( J4 ) * TEMP
            EMIN = MIN( Z( J4-1 ), EMIN )
   20    CONTINUE
      END IF
*
*     Unroll last two steps. 
*
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*(N0-2) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      TEMP = Z( J4P2+2 )/Z( J4-2 )
      Z( J4 ) = Z( J4P2 ) * TEMP
      DNM1 = DNM2 * TEMP  - TAU
      DMIN = MIN( DMIN, DNM1 )
*
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      TEMP = Z( J4P2+2 )/Z( J4-2 )
      Z( J4 ) = Z( J4P2 ) * TEMP
      DN = DNM1 * TEMP  - TAU
      DMIN = MIN( DMIN, DN )
*
      Z( J4+2 ) = DN
      Z( J4+4 ) = EMIN
      RETURN
*
*     End of DLASQ5A
*
      END
