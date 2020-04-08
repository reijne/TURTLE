      SUBROUTINE DUMPER(NIT,OIIR,IDK,IBIAS,LJ)
CSUT                    LG=NIT(667)
CVP   
CVP  NEUE VERSION VON HUS OHNE AUFRUF VON GETMEM
CVP
C
C --- SORTIERUNG UND AUSCHRIFT DER 4INDEXINTEGRALE
C
C   : LG ANZAHL DER INTEGRALE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'wrkspc.inc'
      REAL*8 ONE0
      INTEGER*4 OIIR,NIT,IDK,IBIAS,LJ
      INTEGER*4 MOMAX, MXO
      PARAMETER (MOMAX=400)
      PARAMETER (MXO = MOMAX*(MOMAX+1)/2+1)
      PARAMETER (MDONE = MXO)
      DIMENSION NIT(*),OIIR(*),IDK(*),IBIAS(*),LJ(*)
      DIMENSION IOLIST(20)
      DIMENSION IADR(106)
      DIMENSION ONE0(2*MXO)
      DIMENSION NBAS(8),NORB(8),NFRO(8),NDEL(8)
CSUT
C BE SURE, THAT MAXBUF IS LARGER THAN THE ACTUAL NUMBER OF INTEGRALS.
C
      PARAMETER (NDIMH=2000)
CVP   PARAMETER (MAXREC=5000)
C      PARAMETER (MAXREC=3000)
cbe      PARAMETER (MAXREC=10000)
C      PARAMETER (MAXREC=20000)
C      PARAMETER (MAXREC=60000)
      PARAMETER (MAXREC=20000)
CVP   PARAMETER (MAXREC=20)
      PARAMETER (MAXBUF=MAXREC*NDIMH)
C
      DOUBLE PRECISION BOB,CORE
      REAL*4 COUL,EXC
      COMMON /CONEI/ COUL(MDONE),EXC(MDONE),BOB(MDONE),CORE
C
      REAL*4 STONN
C      COMMON /BUF31/ STONN(MAXBUF)
      DIMENSION STONN(MAXBUF)
CSUT NEUER ZWEIELEKTRONBUFFER FUER TRAINT-FILE
C   NFILM : FORMAT DES MOLCAS-2 FILES
C   RINTMO : ARRAY DER MOLCAS-2 INTEGRALE
C
       INTEGER*4 NFILM
       PARAMETER (NFILM=96000)
       REAL*8 RINTMO
       DIMENSION RINTMO(NFILM)
C  BUFFER FUER ZWEIELEKTRONENINTEGRALE
*
*     DUMP THE CONTENT OF THE ONE-ELECTRON FILE
*
      WRITE(6,*) '<><><><><><><><><><><><><><><><>'
      WRITE(6,*) '<<<<<< D U M P E R   V1.0 >>>>>>'
      WRITE(6,*) '<><><><><><><><><><><><><><><><>'
      WRITE(6,*)
      CALL DANAME(33,'TRAONE')
      CALL GSLIST(IOLIST,7,IADR,64,ECOR,2,NSYM,1,
     *            NBAS,8,NORB,8,NFRO,8,NDEL,8)
      IDISK=0
      CALL DAFILE(33,4,IOLIST,IDUM,IDISK)
      PRINT *,' NSYM',NSYM
      PRINT *,' NBAS',NBAS
      PRINT *,' NORB',NORB
      PRINT *,' NFRO',NFRO
      PRINT *,' NDEL',NDEL
      PRINT *,' ECOR',ECOR
      CALL DACLOS(33,'TRAONE')
      LG = NIT(667)
      WRITE(6,*) 'DUMPER ARBEITET FUER ',LG,' INTEGRALE.'
CVP FEHLER, FALLS LG=MAXBUF ?
CVP   NPASS = LG/MAXBUF + 1
      NPASS = (LG-1)/MAXBUF + 1
      WRITE(6,*) 'DUMPER BENOETIGT ',NPASS,' DURCHGAENGE.'
*
*     DUMP THE CONTENT OF THE TWO-ELECTRON FILE
*
C ------------------------------------------
C ------------------------------------------
C ------------------------------------------
C ------------------------------------------
C ------------------------------------------
* DUMPER PASSES THROUGH MOTRA FILE
      DO 8000 IIPASS=1,NPASS
C ------------------------------------------
C BERECHNE DIE ANZAHL INTEGRALE IM PASS
       IF (IIPASS.EQ.NPASS) THEN
CVP     FEHLER ?
CVP      NIT667 = LG - (IIPASS-1)*MAXBUF
         NIT667 = LG 
       ELSE
         NIT667 = IIPASS*MAXBUF
       END IF
C - BEGINN DER INTEGRALE
       NIT000 = (IIPASS-1)*MAXBUF + 1
C - ENDE DER INTEGRALE
       NITEE = IIPASS*MAXBUF
C ------------------------------------------
       WRITE(6,*) 'ARBEITE AM ',IIPASS,'. DURCHGANG'
       WRITE(6,*) 'INTEGRAL-NR.:',NIT000,' BIS',NIT667
      IDISK=0
      CALL DANAME(34,'TRAINT')
      CALL DAFILE(34,2,IADR,106,IDISK)
      PRINT *,' IADR',IADR
      PRINT *,'  ----------------------------------------------',
     *                       '---------------------------------'
      IBATCH=0
      DO 1100 ISYM=1,NSYM
        IB=NBAS(ISYM)
        IO=NORB(ISYM)
        DO 1200 JSYM=1,ISYM
          JB=NBAS(JSYM)
          JO=NORB(JSYM)
          IF( ISYM.EQ.JSYM ) THEN
            IJB=IB*(IB+1)/2
            IJO=IO*(IO+1)/2
          ELSE
            IJB=IB*JB
            IJO=IO*JO
          END IF
          IJSYM=1+IEOR(ISYM-1,JSYM-1)
          DO 1300 KSYM=1,ISYM
            KB=NBAS(KSYM)
            KO=NORB(KSYM)
            LSYMMX=JSYM
            IF( KSYM.NE.ISYM ) LSYMMX=KSYM
            DO 1400 LSYM=1,LSYMMX
              LB=NBAS(LSYM)
              LO=NORB(LSYM)
              IF( KSYM.EQ.LSYM ) THEN
                KLB=KB*(KB+1)/2
                KLO=KO*(KO+1)/2
              ELSE
                KLB=KB*LB
                KLO=KO*LO
              END IF
              KLSYM=1+IEOR(KSYM-1,LSYM-1)
              IF ( IEOR(IJSYM-1,KLSYM-1).EQ.0 ) THEN
                IF ( IB*JB*KB*LB.NE.0 ) THEN
                  IF ( IO*JO*KO*LO.NE.0 ) THEN
                    IBATCH=IBATCH+1
                    NINTS=IJO*KLO
                    IF ( ISYM.EQ.KSYM .AND. JSYM.EQ.LSYM )
     *              NINTS=IJO*(IJO+1)/2
CSUT                    PRINT
CSUT     *              '(2X,I3,2X,4I2,2X,4I3,2I4,I7,2X,4I3,2I4,I7,I9)',
CSUT     *                      IBATCH,
CSUT     *                      ISYM,JSYM,KSYM,LSYM,
CSUT     *                      IB,JB,KB,LB,IJB,KLB,IJB*KLB,
CSUT     *                      IO,JO,KO,LO,IJO,KLO,IJO*KLO,
CSUT     *                      NINTS
C                    WRITE(6,*) 'ISYM,JSYM,KSYM,LSYM',ISYM,JSYM,KSYM,LSYM
C                    WRITE(6,*) 'IO,JO,KO,LO',IO,JO,KO,LO
C                    WRITE(6,*) 'IB,JB,KB,LB',IB,JB,KB,LB
C                    WRITE(6,*) 'NINTS = ',NINTS,' IBATCH=',IBATCH
CSUT94                    CALL GETMEM('TUVX','ALLO','REAL',LW1,NINTS)
*                   NBUF=(NINTS+9600-1)/9600
                    IDISK=IADR(IBATCH)
                    IST = 1
CERROR              CALL DAFILE(34,2,WORK(IST),9600,IDISK)
CVP             19200, DA DER NEARLY I/O UEBER INTEGER*4 WORTE GEHT
                    CALL DAFILE(34,2,RINTMO(IST),19200,IDISK)
C		write(*,*) "OK"
CSUT
*                   IST=LW1
*                   DO IBUF=1,NBUF
*                     CALL DAFILE(34,2,WORK(IST),9600,IDISK)
*                     IST=IST+9600
*                   END DO
C                   PRINT '(10F12.8)',(WORK(IXXX),IXXX=LW1,LW1+NINTS-1)
CSUT
                    IOUT=0
                    NOR = KO
                    NOS = LO
                    NOP = IO
                    NOQ = JO
                    NSP = ISYM
                    NSQ = JSYM
                    NSR = KSYM
                    NSS = LSYM
CSUT
C              WRITE(6,*) 'NOP,NOQ,NOR,NOS',NOP,NOQ,NOR,NOS
C   NOP : NUMBER OF ORBITALS IN SYMMETRY I
C   NOQ : NUMBER OF ORBITALS IN SYMMETRY J
C   NOR : NUMBER OF ORBITALS IN SYMMETRY K
C   NOS : NUMBER OF ORBITALS IN SYMMETRY L
                    IOUT=0
                    DO 2300 NV=1,NOR
                      NXM=NOS
                      IF(NSR.EQ.NSS)NXM=NV
                      DO 2400 NX=1,NXM
                        NTM=1
                        IF(NSP.EQ.NSR)NTM=NV
                        DO 2100 NT=NTM,NOP
                         NUMIN=1
                         IF(NSP.EQ.NSR.AND.NT.EQ.NV)NUMIN=NX
                         NUMAX=NOQ
                         IF(NSP.EQ.NSQ)NUMAX=NT
                         DO 2200 NU=NUMIN,NUMAX
                            IOUT=IOUT+1
CVP                 DAS HOCHZAEHLEN VON IDISK PASSIERT WOHL IM NEARLY I/O
                            IF (IOUT.GT.9600) THEN
C                            write(*,*) "iout=",iout
                              CALL DAFILE(34,2,RINTMO(IST),19200,IDISK)
C                             WRITE(6,*) 'READDA : IDISK=',IDISK
                              IOUT = 1
                            END IF
CSUT : DURCH FELD --- MJ??? ERSETZEN !
                            NADD = 0
                            DO II0=1,(ISYM-1)
                               NADD = NADD + NORB(II0)
                            END DO
                            II = NT + NADD
                            NADD = 0
                            DO II0=1,(JSYM-1)
                               NADD = NADD + NORB(II0)
                            END DO
                            JJ = NU + NADD
                            NADD = 0
                            DO II0=1,(KSYM-1)
                               NADD = NADD + NORB(II0)
                            END DO
                            KK = NV + NADD
                            NADD = 0
                            DO II0=1,(LSYM-1)
                               NADD = NADD + NORB(II0)
                            END DO
                            LL = NX + NADD
                            SINT = RINTMO(IOUT)
CDEBUG                IF (DABS(SINT).LT.0.0001D0) THEN
CDEBUG                 WRITE(6,*)IOUT,RINTMO(IOUT)
CDEBUG                END IF
*  --- SORTIERE INDEXE IN CANONISCHE ORDNUNG
                            CALL CSORT(II,JJ,KK,LL)
C                        IF(IBATCH.EQ.46) THEN
C                        write(*,*) "vorher"
C                            WRITE(6,*) 'IJKL,SINT,NEWPOS',
C     &                      II,JJ,KK,LL,SINT,NEWPOS
C    					end if
                            CALL RS4IND (II,JJ,KK,LL,SINT,NEWPOS,
     &                                   NIT,OIIR,IDK,IBIAS,LJ)
C                        IF(IBATCH.EQ.46) write(*,*) "nachher"
CDEBUG                            WRITE(6,*) 'IJKL,SINT,NEWPOS',
CDEBUG     &                      II,JJ,KK,LL,SINT,NEWPOS
** <---- ABBILDUNG AUF BUFFER
CSUT : HIER BUFFER - ABFRAGE !!!A
CSUT : DIE GUTEN INS KOERBCHEN, DIE SCHLECHTEN ...
CVP     .GT. WAR FALSCH: CHRISTOF HAETTIG, PRIVATE MITTEILUNG
CVP                         IF (NEWPOS.GT.NIT000) THEN
                            IF (NEWPOS.GE.NIT000) THEN
                             IF (NEWPOS.LE.NIT667) THEN
                               NEWPOS = NEWPOS - NIT000 + 1
                               STONN(NEWPOS) = (SINT)
                             END IF
                            END IF
 2200                     CONTINUE
 2100                   CONTINUE
 2400                 CONTINUE
 2300               CONTINUE
CSUT
CSUT                WRITE(6,*) 'IOUT=',IOUT
CSUT94                    CALL GETMEM('TUVX','FREE','REAL',LW1,NINTS)
                  END IF
                END IF
              END IF
 1400       CONTINUE
 1300     CONTINUE
 1200   CONTINUE
 1100 CONTINUE
      PRINT *,'  ----------------------------------------------',
     *                       '---------------------------------'
      CALL DACLOS(34,'TRAINT')
CDEBUG
*     WRITE(6,*) 'ENDE DER ZWEIELEKTRONENINTEGRALE !'
      WRITE(6,*)
      WRITE(6,*) 'BIN AM DUMPEN'
      WRITE(6,*)
*     WRITE(6,*) 'LG=,IOUT=',NIT(667),IOUT
C NUMBER OF  RECORDS
CVP   N2R = NIT667/NDIMH
      N2R = (NIT667-NIT000+1)/NDIMH
C REST
CVP   N2RR = NIT667- N2R*NDIMH
      N2RR = NIT667-NIT000+1 - N2R*NDIMH
CVP
      WRITE(6,*) ' ANZAHL ZU SCREIBENDER RECORDS ',N2R
      WRITE(6,*) ' ANZAHL DER VERBLEIBENDEN INTEGRALE FUER RESTRECORD '
     &            ,N2RR
C
      NSTART=1
      NEND  = 2000
      DO I=1,N2R
        WRITE(31) (STONN(II),II=NSTART,NEND)
        NSTART=NSTART+ NDIMH
        NEND  = NEND + NDIMH
      END DO
C --- END PASS IIPASS
 8000 CONTINUE
C
      IF (N2RR.NE.0) THEN
        WRITE(31) (STONN(II),II=NSTART,NEND)
        WRITE(6,*) ' RESTREKORD GESCHRIEBEN'
      END IF
C
*******************************************************************
*
*     ONE-ELECTRON INTEGRAL SECTION
*
*******************************************************************
      CALL GETH0(ONE0,ECORE)
      WRITE(6,*) 'ONE-ELECTRON PART'
      WRITE(6,'(4F18.9,2X)') (ONE0(I),I=1,200)
* <--- ABBILDUNG AUF BUFFER
      DO I=1,MXO
        BOB(I) = (ONE0(I))
      END DO
* <--- SCHREIBE LETZTE RECORDS
      CALL KOTZ
      RETURN
*******************************************************************
      END
