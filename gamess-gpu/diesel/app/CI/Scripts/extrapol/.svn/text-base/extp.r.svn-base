      SUBROUTINE EXTP
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'prolog.h'
      INCLUDE 'param.h'
      INCLUDE 'common.h'
      REAL*8 s(7),del(7),p1(7),ph(7),pg(7),txsum(7,mxroot),ex(mxroot,3)
      CHARACTER*4 itext(4)
      REAL*8 EL1(0:7),EL2(0:7),BETA,ESTLIM
      LOGICAL ZBW,ZCONT,ZPT
      LOGICAL lhead
      CHARACTER*4 ITYES,ITNO
      DATA ITYES /' YES'/,ITNO /'  NO'/
      IF (inflg(lndeci).NE.0) RETURN
      CALL tenter('extp')
      ZPT=MOD(JXT,2).EQ.1
      ZCONT=MOD(JXT,8).GE.4
      ZBW=MOD(JXT,16).GE.8
      IF (ZBW    .AND. .NOT. ZPT) ZPT=.TRUE.
      IF (ZCONT  .AND. .NOT. ZPT) ZPT=.TRUE.
      DO 540 I=1,4
540   ITEXT(I)=ITNO
      IF (ZPT) ITEXT(1)=ITYES
      IF (inflg(lmerge).NE.0) ITEXT(2)=ITYES
      IF (ZCONT)  ITEXT(2)=ITYES
      IF (ZCONT)  ITEXT(3)=ITYES
      IF (ZBW)    ITEXT(4)=ITYES
      WRITE (JTAPE,541) ITEXT
541   FORMAT(/' ** EXTRAPOLATION TYPE **'//,
     *        ' SPECIAL PERTURBATION',T40,A4/,
     *        ' SPECIAL MERGE',T40,A4/,
     *        ' SPECIAL MERGE CONTINUATION',T40,A4/,
     *        ' BW PERTURBATION',T40,A4/)
      NXT=3
      MTAPE=11
      IF (ZPT) MTAPE=IFSAVE
      IFSAVE=MTAPE
      CALL ICREWI(MTAPE,IRC)
      IF (.NOT. ZCONT) GOTO 81
C
C     GET RESTART DATA
C
      READ (MTAPE) NTCH,NROOT,NROTR,IBG,IPERM
      READ (MTAPE) EW,UEL,TSSUM,TXSUM,PROJ
      GOTO 80
C
C     NORMAL RUN
C
81    NTCH=NTCH+1
      IBG=IBG-1
      WRITE (JTAPE,1) (IPERM(I),I=1,nroot)
  1   FORMAT(///5X,'PERMUTED ORDER OF GEYSER CASES IS:',9I5//)
C
C     SAVE DATA FOR RESTART
C
      WRITE (MTAPE) NTCH,NROOT,NROTR,IBG,IPERM
      WRITE (MTAPE) EW,UEL,TSSUM,TXSUM,PROJ
80    IF (ZBW) IBG=-IBG
      IF (.NOT. ZPT) GOTO 96
      STOP 'PTC6'
C6    CALL PTC6(EW,UEL,TSSUM,TRSUM,TXSUM,PTSUM,NTCH,NROOT,NROTR,IBG)
C6    IBG=0
C6    DO 95 I=1,nroot
C6 95 IPERM(I)=I
C6    GOTO 97
C
C     USUAL WEIGHTED EXTRAPOLATION (AS IN CONDOXC4)
C
96    DO 11 I=1,NTCH
      DO 12 J=1,NROOT
      TXSUM(I,J)=0.D0
      X=0D0
      DO 13 K=1,NROTR
C#    WRITE (JTAPE,106) I,J,K,UEL(I,J,K),RSUM(IBG+I,K),TRSUM(IBG+I,K)
C#106 FORMAT(' I,J,K,UEL,RSUM,TRSUM ',3I5,3D14.6)
      Y=UEL(I,J,K)
      X=X+Y*Y
13    TXSUM(I,J)=TXSUM(I,J)+Y*Y*RSUM(IBG+I,K)
12    TXSUM(I,J)=TXSUM(I,J)/X
11    CONTINUE
C
C     LOOP OVER EXTRAPOLATION TYPES
C
97    DO 33 K=1,NXT
      lhead=.TRUE.
      DO 2 I=1,nroot
      JBG=IBG
      DO 3 J=1,NTCH
      JBG=JBG+1
      S(J)=TSSUM(J, I)+EZERO
      JGPE=IPERM(I)
      GOTO (41,42,43),K
41    DEL(J)=TRSUM(JBG,JGPE)
      GOTO 3
42    DEL(J)=TXSUM(J,I)
      GOTO 3
CBAH  DEL(J)=PTSUM(J,I)
43    DEL(J)=EMP10(JBG,JGPE)
  3   CONTINUE
      IF (DABS(del(1)).LT.1.D-8) THEN
       ex(i,k)=s(1)
       GOTO 2
      ENDIF
      IF (lhead) THEN
       WRITE (JTAPE,31) K
31     FORMAT(//' --- EXTRAPOLATION TYPE ',I1,' ---'//)
       lhead=.FALSE.
      ENDIF
      WRITE (JTAPE,4) I
  4   FORMAT(//5X,'ROOT NO.',I3/)
      WRITE (JTAPE,5) (S(J),J=1,NTCH)
      WRITE (JTAPE,5) (DEL(J),J=1,NTCH)
    5 FORMAT(/10X,7F14.6/)
C      LEVIN EXTRAPOLATION
      BETA=1.1D0
      DO J=1,NTCH
       MU=NTCH-J+1
       CALL GLEVIN(S(MU),DEL(MU),BETA,J-1,EL1,EL2,NTCH,ESTLIM)
      ENDDO
      ex(i,k)=estlim
      WRITE (JTAPE,122) ntch,ex(i,k)
  122  FORMAT(//10X,
     * ' ENERGY (',I1,'-POINT LEVINE EXTRAPOLATION): ',G20.8)
      IF (inflg(l2point).NE.0) THEN
C      2 point extrapolation
       WRITE (JTAPE,5) ESTLIM
       CALL GLEVIN(S(2),DEL(2),BETA,0,EL1,EL2,NTCH,ESTLIM)
       CALL GLEVIN(S(1),DEL(1),BETA,1,EL1,EL2,NTCH,ESTLIM)
       WRITE (JTAPE,5) ESTLIM
       IF (NTCH.GT.1) GOTO 90
       X=0.D0
       GOTO 91
90     X=S(2)-S(1)
       Y=DEL(2)-DEL(1)
       X=X/Y
91     DO 6 J=1,NTCH
       XY=S(J)
       XZ=-DEL(J)
       P1(J)=XY+XZ
       PH(J)=XY+0.50D0*XZ
    6  PG(J)=XY+X*XZ
       WRITE (JTAPE,7) (P1(J),J=1,NTCH)
    7  FORMAT(/5X,'AT LAMBDA=1.0',7F14.6)
       WRITE (JTAPE,8) (PH(J),J=1,NTCH)
    8  FORMAT(/5X,'AT LAMBDA=0.5',7F14.6)
       WRITE (JTAPE,10) X
   10  FORMAT(//20X,'LAMBDA OPTIMUM= ',F11.8)
       XZ=X+0.05D0
       Y=X-0.05D0
       DO 20 J=1,NTCH
       P1(J)=S(J)-XZ*DEL(J)
   20  PH(J)=S(J)-Y*DEL(J)
       WRITE (JTAPE,21) (P1(J),J=1,NTCH)
   21  FORMAT(/5X,'AT L OPT+0.05',7F14.6)
       WRITE (JTAPE,22) (PH(J),J=1,NTCH)
   22  FORMAT(/5X,'AT L OPT-0.05',7F14.6)
       WRITE (JTAPE,9) (PG(J),J=1,NTCH)
    9  FORMAT(/5X,'AT LAMBDA=OPT',7F14.6)
       XZ=PG(1)
       Y=DEL(1)*0.05D0*X
       EX(I,K)=XZ
       WRITE (JTAPE,121) XZ,Y
  121  FORMAT(//20X,' EXTRAPOLATED ENERGY =',F14.6,'+/-',F8.6)
      ENDIF
2     CONTINUE
33    CONTINUE
      WRITE(MTAPE) EX
      CALL texit('extp')
      RETURN
      END
      SUBROUTINE GLEVIN(SOFN,ROFN,BETA,N,ARUP,ARLO,LARRAY,ESTLIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARUP(0:LARRAY),ARLO(0:LARRAY)
      PARAMETER ( HUGE = 1.E+60 , TINY = 1.E-60 , ONE = 1.E0 )
      ARUP(N)   = SOFN / ROFN
      ARLO(N)   = ONE  / ROFN
      IF ( N .GT. 0 ) THEN
        ARUP(N-1) = ARUP(N) - ARUP(N-1)
        ARLO(N-1) = ARLO(N) - ARLO(N-1)
        IF (N.GT.1) THEN
          BN1     = BETA + FLOAT(N-1)
          BN2     = BETA + FLOAT(N)
          COEF    = BN1 / BN2
          DO 10 J = 2,N
            FACT  = (BETA+FLOAT(N-J)) * COEF**(J-2) / BN2
            ARUP(N-J)  = ARUP(N-J+1) - FACT*ARUP(N-J)
            ARLO(N-J)  = ARLO(N-J+1) - FACT*ARLO(N-J)
10        CONTINUE
        END IF
      END IF
        IF (ABS(ARLO(0)).LT.TINY) THEN
          ESTLIM = HUGE
        ELSE
          ESTLIM = ARUP(0)/ARLO(0)
        END IF
      RETURN
      END
