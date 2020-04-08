*
* $Id: ddstqds.f,v 1.2 2000-10-26 15:38:31 psh Exp $
*
c
      SUBROUTINE DDSTQDS( LAMBDA, DELTA, N, B1, BN, L, D,
     $     LD, LPLUS, DPLUS, T )
*
*     -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*     
*     .. Scalar Arguments ..
	implicit none
      INTEGER            N, B1, BN
      DOUBLE PRECISION   DELTA, LAMBDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), LPLUS( * ),
     $                   LD( * ), T(*), DPLUS( * )
*     ..
*
*  Purpose
*  =======
*
*     DDSTQDS computes the L+ D+ L+^T decomposition of the submatrix
*     indexed from B1 to BN of (L*D*L^T - (LAMBDA+DELTA)*I)
*     by the differential form of the stationary qd algorithm.
*     
*     Arguments
*     =========
*     
*     LAMBDA  (input) DOUBLE PRECISION
*     The shift.
*     
*     DELTA   (input) DOUBLE PRECISION
*     Lower order bits of the shift.
*     
*     N       (input) INTEGER
*     The order of the matrix L * D * L^T.
*     
*     B1      (input) INTEGER
*     Starting index of the submatrix (of L * D * L^T).
*     
*     BN      (input) INTEGER
*     Last index of the submatrix (of L * D * L^T).
*     
*     L       (input) DOUBLE PRECISION array, dimension (N-1)
*     The (n-1) subdiagonal elements of the bidiagonal matrix
*     L, in elements 1 to N-1.  L(N) need not be set.
*     
*     D       (input) DOUBLE PRECISION array, dimension (N)
*     The n diagonal elements of the diagonal matrix D.
*     
*     LD      (input) DOUBLE PRECISION array, dimension (N-1)
*     The n-1 elements L(i)*D(i).
*     
*     LPLUS   (output) DOUBLE PRECISION array, dimension (N-1)
*     The (n-1) subdiagonal elements of L+.
*     
*     DPLUS   (output) DOUBLE PRECISION array, dimension (N)
*     The n diagonal elements of D+.
*     
*     T       (output) DOUBLE PRECISION array, dimension (N)
*     Intermediate results of the algorithm.
*     
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, ONE = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   TMP, EPS
*     ..
*     .. External Functions ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
      integer            doprnt1,doprnt2
      common             doprnt1,doprnt2

*     ..
*     .. Executable Statements ..
*
 310  CONTINUE
      IF( B1.EQ.1 ) THEN
         T( B1 ) = -LAMBDA
      ELSE
         T( B1 ) = LD( B1-1 ) * L( B1-1 ) - LAMBDA
      END IF

      DO I = B1, BN-1
         DPLUS( I ) = ( D( I ) + T( I ) ) - DELTA
         LPLUS( I ) = LD( I ) / DPLUS( I )
         T ( I+1 ) = ( T( I )-DELTA ) * LPLUS( I ) * L( I ) - LAMBDA
      END DO
      DPLUS( BN ) = ( D( BN ) + T( BN ) ) - DELTA

      IF(.NOT.( ( DPLUS( BN ).GT.ZERO) .OR.
     $          ( DPLUS( BN ).LT.ONE ) ) ) THEN
         IF( B1.EQ.1 ) THEN
            T( B1 ) = -LAMBDA
         ELSE
            T( B1 ) = LD( B1-1 ) * L( B1-1 ) - LAMBDA
         END IF

         DO I = B1, BN-1
            DPLUS( I ) = ( D( I ) + T( I ) ) - DELTA
            LPLUS( I ) = LD( I ) / DPLUS( I )
*
*  Need to check the next few lines
*
            IF(LPLUS( I ).EQ.0 .AND. ONE/(T( I )-DELTA).EQ.ZERO ) THEN
               T ( I+1 ) = L( I ) * LD( I ) - LAMBDA
            ELSE IF((T(I)-DELTA).EQ.ZERO .AND. ONE/LPLUS(I).EQ.0 ) THEN
               T ( I+1 ) = -ONE / DPLUS( I )
            ELSE
              T ( I+1 ) = ( T( I )-DELTA )*LPLUS( I )*L( I ) - LAMBDA
            END IF
         END DO
         DPLUS( BN ) = ( D( BN ) + T( BN ) ) - DELTA
      END IF

      IF(.NOT.( ( DPLUS( BN ).GT.ZERO) .OR.
     $          ( DPLUS( BN ).LT.ONE ) ) ) THEN
         print *,"DDSTQDS : NaN detected!"
         print *,"dplus = [ "
         do i = b1,bn
            write(*,101) dplus(i)
         end do
         print *,"];"
         print *,"lplus = [ "
         do i = b1,bn-1
           write(*,101) lplus(i)
         end do
         print *,"];"
         print *,"t  = [ "
         do i = b1,bn
           write(*,101) t(i)
         end do
         print *,"];"
         stop
         EPS = 0.111022302462515654E-15
         LAMBDA = LAMBDA*(ONE+EPS)
         print *,"New Lambda = "
         write(*,101) LAMBDA
         GO TO 310
      END IF


      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*
*     End of DDSTQDS
*
      END
