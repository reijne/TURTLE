      SUBROUTINE RDVEC(LU,NSYM,NBAS,NORB,CMO,OCC,LOCC,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*)
      CHARACTER*(*) TITLE
      CHARACTER LINE*80,FMT*40
      LOGICAL SET
      FMT='(4E18.12)'
      KCMO  = 0
      NDIV  = 4
      TITLE = ' '
      SET   = .FALSE.
      DO 100 ISYM=1,NSYM
         DO 110 IORB=1,NORB(ISYM)
            DO 111 IBAS=1,NBAS(ISYM),NDIV
112            READ(LU,'(A80)',END=999,ERR=999) LINE
               IF(LINE(1:1).EQ.'*') THEN
                  IF(.NOT. SET) THEN
                     TITLE=LINE
                     SET=.TRUE.
                  END IF
                  GOTO 112
               END IF
               READ(LINE,FMT)
     &             (CMO(I+KCMO),I=IBAS,MIN(IBAS+3,NBAS(ISYM)))
111         CONTINUE
            KCMO=KCMO+NBAS(ISYM)
110      CONTINUE
100   CONTINUE
      IF(LOCC.EQ.0) RETURN
      KOCC=0
      DO 200 ISYM=1,NSYM
         DO 210 IORB=1,NORB(ISYM),NDIV
212         READ(LU,'(A80)',END=999,ERR=999) LINE
            IF(LINE(1:1).EQ.'*') THEN
               IF(.NOT. SET) THEN
                  TITLE=LINE
                  SET=.TRUE.
               END IF
               GOTO 212
            END IF
            READ(LINE,FMT) (OCC(I+KOCC),I=IORB,MIN(IORB+3,NORB(ISYM)))
210      CONTINUE
         KOCC=KOCC+NORB(ISYM)
200   CONTINUE
      RETURN
999   CONTINUE
      WRITE(*,*) '* ERROR IN RDVEC WHILE READING VECTOR SOURCE FILE *'
      STOP 20
      END
