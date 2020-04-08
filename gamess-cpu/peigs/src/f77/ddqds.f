*
* $Id: ddqds.f,v 1.2 2000-10-26 15:38:31 psh Exp $
*
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 11/ 3/97 at 12:14:42
C****************************************************************
C
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2
C   Version 6.1C1 on 11/ 3/97 at 12:14:02
C****************************************************************
C      
       
      SUBROUTINE DDQDS( LAMBDA, DELTA, N, B1, BN, L, D,
     $                  LLD, UMINUS, DMINUS, P )
*      
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*      
*     .. Scalar Arguments ..
      implicit none
      INTEGER            N, B1, BN
      DOUBLE PRECISION   DELTA, LAMBDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), P( * ), DMINUS( * ),
     $                   LLD( * ), UMINUS( * )
*     ..
*      
*  Purpose
*  =======
*      
*  DDQDS computes the U- D- U-^T decomposition of the submatrix
*  indexed from B1 to BN of (L*D*L^T - (LAMBDA+DELTA)*I)
*  by the differential form of the progressive qd algorithm.
*      
*  Arguments
*  =========
*      
*  LAMBDA  (input) DOUBLE PRECISION
*          The shift.
*      
*  DELTA   (input) DOUBLE PRECISION
*          Lower order bits of the shift.
*      
*  N       (input) INTEGER
*          The order of the matrix L * D * L^T.
*      
*  B1      (input) INTEGER
*          Starting index of the submatrix (of L * D * L^T).
*      
*  BN      (input) INTEGER
*          Last index of the submatrix (of L * D * L^T).
*      
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of the bidiagonal matrix
*          L, in elements 1 to N-1.  L(N) need not be set.
*      
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D.
*      
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*L(i)*D(i).
*      
*  UMINUS  (output) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) superdiagonal elements of U-.
*      
*  DMINUS  (output) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D-.
*      
*  P       (output) DOUBLE PRECISION array, dimension (N)
*          Intermediate results of the algorithm.
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
 210  CONTINUE
      P( BN ) = D( BN ) - LAMBDA
      DO I = BN-1, B1, -1
         DMINUS( I+1 ) = ( LLD( I ) + P( I+1 ) ) - DELTA
         TMP = D( I ) / DMINUS( I+1 )
         UMINUS( I ) = L( I ) * TMP
         P ( I ) = ( P( I+1 ) - DELTA ) * TMP - LAMBDA
      END DO
      IF( B1.EQ.1 ) THEN
         DMINUS( B1 ) = P( B1 ) - DELTA
      ELSE
         DMINUS( B1 ) = ( LLD( B1-1 ) + P( B1 ) ) - DELTA
      END IF
       
      IF(.NOT.( ( DMINUS( B1 ).GT.ZERO) .OR.
     $          ( DMINUS( B1 ).LT.ONE ) ) ) THEN
         P( BN ) = D( BN ) - LAMBDA
         DO I = BN-1, B1, -1
            DMINUS( I+1 ) = ( LLD( I ) + P( I+1 ) ) - DELTA
            TMP = D( I ) / DMINUS( I+1 )
            UMINUS( I ) = L( I ) * TMP
*      
*  Need to check the next few lines
*      
            IF( TMP.EQ.0 .AND. ONE/(P( I+1) - DELTA).EQ.ZERO ) THEN
               P ( I ) = ( L( I ) ) - LAMBDA
            ELSE IF((P( I+1)-DELTA).EQ.ZERO .AND. ONE/TMP.EQ.0) THEN
               P ( I ) = - ONE / DMINUS( I+1 )
            ELSE
              P ( I ) = ( P( I+1 ) - DELTA ) * TMP - LAMBDA
            END IF
         END DO
         IF( B1.EQ.1 ) THEN
            DMINUS( B1 ) = P( B1 ) - DELTA
         ELSE
            DMINUS( B1 ) = ( LLD( B1-1 ) + P( B1 ) ) - DELTA
         END IF
      END IF
       
      IF(.NOT.( ( DMINUS( B1 ).GT.ZERO) .OR.
     $          ( DMINUS( B1 ).LT.ONE ) ) ) THEN
         print *,"DDQDS : NaN detected!"
         print *,"dminus : "
         do i = b1,bn
           write(*,101) dminus(i)
         end do
         print *,"];"
         print *,"uminus = [ "
         do i = b1,bn-1
           write(*,101) uminus(i)
         end do
         print *,"];"
         print *,"p  = [ "
         do i = b1,bn
           write(*,101) p(i)
         end do
         print *,"];"
         stop
         EPS = 0.111022302462515654E-15
         LAMBDA = LAMBDA*(ONE+EPS)
         print *,"New Lambda = "
         write(*,101) LAMBDA
         GO TO 210
      END IF
       
       
      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*      
*     End of DDQDS
*      
      END
