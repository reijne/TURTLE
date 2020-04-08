      SUBROUTINE DHQRII1(N,M,A,E,V,MX) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C
      DIMENSION A(*), E(MX), V(MX,MX)
************************************************************* 
* 
* HQRII IS A DIAGONALISATION ROUTINE, WRITTEN BY YOSHITAKA BEPPU OF 
*       NAGOYA UNIVERSITY, JAPAN. 
*       FOR DETAILS SEE 'COMPUTERS & CHEMISTRY' VOL.6 1982. PAGE 000. 
* 
* ON INPUT    A       = MATRIX TO BE DIAGONALIZED 
*             N       = SIZE OF MATRIX TO BE DIAGONALIZED. 
*             M       = NUMBER OF EIGENVECTORS NEEDED. 
*             E       = ARRAY OF SIZE AT LEAST N 
*             V       = ARRAY OF SIZE AT LEAST NMAX*M 
* 
* ON OUTPUT   E       = EIGENVALUES 
*             V       = EIGENVECTORS IN ARRAY OF SIZE NMAX*M 
* 
* ORIGIN : MOPAC 5.0 / IBM 3090 
************************************************************************ 

      DIMENSION W(5,10000)        
 
      IF(N.EQ.1 .AND. M.EQ.1) THEN 
            E(1)=A(1) 
            V(1,1)=1.D0 
      RETURN 
      ENDIF 
* 
* EPS3 AND EPS ARE MACHINE-PRECISION DEPENDENT 
* 
      EPS3=1.D-16 
      ZERO=0.D0 
      LL=(N*(N+1))/2+1 
      EPS=1.D-10
      IORD=-1 
      NM1=N-1 
      IF(N.EQ.2) GOTO 80 
      NM2=N-2 
C     HOUSEHOLDER TRANSFORMATION 
      DO 70 K=1,NM2 
         KP1=K+1 
         W(2,K)=A((K*(K+1))/2) 
         SUM=0. 
         DO 10 J=KP1,N 
            W(2,J)=A((J*(J-1))/2+K) 
   10    SUM=W(2,J)**2+SUM 
         S=SIGN(DSQRT(SUM),W(2,KP1)) 
         W(1,K)=-S 
         W(2,KP1)=W(2,KP1)+S 
         A(K+(KP1*(KP1-1))/2)=W(2,KP1) 
         H=W(2,KP1)*S 
         IF(ABS(H).LT.1.D-35) GOTO 70 
C#      IF(H.EQ.0.D0) GOTO 70 
         SUMM=0.D0 
         DO 50 I=KP1,N 
            SUM=0.D0 
            DO 20 J=KP1,I 
   20       SUM=SUM+A(J+(I*(I-1))/2)*W(2,J) 
            IF(I.GE.N) GOTO 40 
            IP1=I+1 
            DO 30 J=IP1,N 
   30       SUM=SUM+A(I+(J*(J-1))/2)*W(2,J) 
   40       W(1,I)=SUM/H 
   50    SUMM=W(1,I)*W(2,I)+SUMM 
         U=SUMM*0.5D0/H 
         DO 60 J=KP1,N 
            W(1,J)=W(2,J)*U-W(1,J) 
            DO 60 I=KP1,J 
   60    A(I+(J*(J-1))/2)=W(1,I)*W(2,J)+W(1,J)*W(2,I)+A(I+(J*(J-1))/2) 
   70 A((K*(K+1))/2)=H 
   80 W(2,NM1)=A((NM1*(NM1+1))/2) 
      W(2,N)=A((N*(N+1))/2) 
      W(1,NM1)=A(NM1+(N*(N-1))/2) 
      W(1,N)=0.D0 
      GERSCH=ABS(W(2,1))+ABS(W(1,1)) 
      DO 90 I=1,NM1 
   90 GERSCH=MAX(ABS(W(2,I+1))+ABS(W(1,I))+ABS(W(1,I+1)),GERSCH) 
      DEL=EPS*GERSCH 
      DO 100 I=1,N 
         W(3,I)=W(1,I) 
         E(I)=W(2,I) 
  100 V(I,M)=E(I) 
      IF(DEL.EQ.ZERO)  GOTO  210 
C     QR-METHOD WITH ORIGIN SHIFT 
      K=N 
  110 L=K 
  120 IF(ABS(W(3,L-1)).LT.DEL) GOTO 130 
      L=L-1 
      IF(L.GT.1)  GOTO 120 
  130 IF(L.EQ.K)  GOTO 160 
      WW=(E(K-1)+E(K))*0.5D0 
      R=E(K)-WW 
      Z=SIGN(DSQRT(W(3,K-1)**2+R*R),R)+WW 
      EE=E(L)-Z 
      E(L)=EE 
      FF=W(3,L) 
      R=DSQRT(EE*EE+FF*FF) 
      J=L 
      GOTO 150 
  140 R=DSQRT(E(J)**2+W(3,J)**2) 
      W(3,J-1)=S*R 
      EE=E(J)*C 
      FF=W(3,J)*C 
  150 R=R+1.D-15 
      C=E(J)/R 
      S=W(3,J)/R 
      WW=E(J+1)-Z 
      E(J)=(FF*C+WW*S)*S+EE+Z 
      E(J+1)=C*WW-S*FF 
      J=J+1 
      IF(J.LT.K) GOTO 140 
      W(3,K-1)=E(K)*S 
      E(K)=E(K)*C+Z 
      GOTO 110 
  160 K=K-1 
      IF(K.GT.1) GOTO 110 
*    *    *    *    *    *    *    *    *    *    *    *    * 
* 
*   AT THIS POINT THE ARRAY 'E' CONTAINS THE UN-ORDERED EIGENVALUES 
* 
*    *    *    *    *    *    *    *    *    *    *    *    * 
C     STRAIGHT SELECTION SORT OF EIGENVALUES 
      SORTER=1.D0 
      IF(IORD.LT.0) SORTER=-1.D0 
      J=N 
  170 L=1 
      II=1 
      LL=1 
      DO 190 I=2,J 
         IF((E(I)-E(L))*SORTER .GT. 0.D0) GOTO 180 
         L=I 
         GOTO 190 
  180    II=I 
         LL=L 
  190 CONTINUE 
      IF(II.EQ.LL) GOTO 200 
      WW=E(LL) 
      E(LL)=E(II) 
      E(II)=WW 
  200 J=II-1 
      IF(J.GE.2) GOTO 170 
  210 IF(M.EQ.0) RETURN 
*************** 
*  ORDERING OF EIGENVALUES COMPLETE. 
*************** 
C      INVERSE-ITERATION FOR EIGENVECTORS 
      FN=FLOAT(N) 
      EPS1=1.D-5 
      SEPS=DSQRT(EPS) 
      EPS2=0.05D0 
      RN=0.D0 
      RA=EPS*0.6180339887485D0 
C    0.618... IS THE FIBONACCI NUMBER (-1+DSQRT(5))/2. 
      IG=1 
      DO 430 I=1,M 
         IM1=I-1 
         DO 220 J=1,N 
            W(3,J)=0.D0 
            W(4,J)=W(1,J) 
            W(5,J)=V(J,M)-E(I) 
            RN=RN+RA 
            IF(RN.GE.EPS) RN=RN-EPS 
  220    V(J,I)=RN 
         DO 250 J=1,NM1 
            IF(ABS(W(5,J)).GE.ABS(W(1,J))) GOTO 230 
            W(2,J)=-W(5,J)/W(1,J) 
            W(5,J)=W(1,J) 
            T=W(5,J+1) 
            W(5,J+1)=W(4,J) 
            W(4,J)=T 
            W(3,J)=W(4,J+1) 
            IF(W(3,J).EQ.ZERO) W(3,J)=DEL 
            W(4,J+1)=0.D0 
            GOTO 240 
  230       IF(W(5,J).EQ.ZERO) W(5,J)=DEL 
            W(2,J)=-W(1,J)/W(5,J) 
  240       W(4,J+1)=W(3,J)*W(2,J)+W(4,J+1) 
  250    W(5,J+1)=W(4,J)*W(2,J)+W(5,J+1) 
         IF(ABS(W(5,N)) .LT. EPS3) W(5,N)=DEL 
         DO 310 ITERE=1,5 
            IF(ITERE.EQ.1) GOTO 270 
            DO 260 J=1,NM1 
               IF(W(3,J).EQ.ZERO) GOTO 260 
               T=V(J,I) 
               V(J,I)=V(J+1,I) 
               V(J+1,I)=T 
  260       V(J+1,I)=V(J,I)*W(2,J)+V(J+1,I) 
  270       V(N,I)=V(N,I)/W(5,N) 
            V(NM1,I)=(V(NM1,I)-V(N,I)*W(4,NM1))/W(5,NM1) 
            VN=MAX(ABS(V(N,I)),ABS(V(NM1,I)),1.D-20) 
            IF(N.EQ.2) GOTO 290 
            K=NM2 
  280       V(K,I)=(V(K,I)-V(K+1,I)*W(4,K)-V(K+2,I)*W(3,K))/W(5,K) 
            VN=MAX(ABS(V(K,I)),VN,1.D-20) 
            K=K-1 
            IF(K.GE.1) GOTO 280 
  290       S=EPS1/VN 
            DO 300 J=1,N 
  300       V(J,I)=V(J,I)*S 
            IF(ITERE.GT.1 .AND. VN.GT.1) GOTO 320 
  310    CONTINUE 
C     TRANSFORMATION OF EIGENVECTORS 
  320    IF(N.EQ.2) GOTO 360 
         DO 350 J=1,NM2 
            K=N-J-1 
            IF(A((K*(K+1))/2).EQ.ZERO) GOTO 350 
            KP1=K+1 
            SUM=0.D0 
            DO 330 KK=KP1,N 
  330       SUM=SUM+A(K+(KK*(KK-1))/2)*V(KK,I) 
            S=-SUM/A((K*(K+1))/2) 
            DO 340 KK=KP1,N 
  340       V(KK,I)=A(K+(KK*(KK-1))/2)*S+V(KK,I) 
  350    CONTINUE 
  360    DO 370 J=IG,I 
            IF(ABS(E(J)-E(I)) .LT. EPS2) GOTO 380 
  370    CONTINUE 
         J=I 
  380    IG=J 
         IF(IG .EQ. I) GOTO 410 
C     RE-ORTHOGONALISATION 
         DO 400 K=IG,IM1 
            SUM=0.D0 
            DO 390 J=1,N 
  390       SUM=V(J,K)*V(J,I)+SUM 
            S=-SUM 
            DO 400 J=1,N 
  400    V(J,I)=V(J,K)*S+V(J,I) 
C     NORMALISATION 
  410    SUM=1.D-24 
         DO 420 J=1,N 
  420    SUM=SUM+V(J,I)**2 
         SINV=1.D0/DSQRT(SUM) 
         DO 430 J=1,N 
  430 V(J,I)=V(J,I)*SINV 
      RETURN 
      END 
