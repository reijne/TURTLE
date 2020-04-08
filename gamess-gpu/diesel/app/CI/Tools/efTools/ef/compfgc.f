      SUBROUTINE COMPFGC(GRAD,rfcrd)
C     this routine seems to control calculation
C     of cartesian gradients and hessians @ various levels of theory.
C     from here, rimp2 programs etc are called to determine energies
C     and dipole momenta.
C     ghost
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      DIMENSION rfcrd(*),GRAD(*)
      COMMON /KEYWRD/KEYWRD
      CHARACTER*241 KEYWRD
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /CYCLE / ICYC
      COMMON /ENERG / ELAST, CGRAD
      COMMON /PRTFLG/ IPRT
      CHARACTER*2 ELEMNT              
      COMMON /GEOM  / GEO(3,NUMATM)
      DIMENSION GXYZ(3,NUMATM)
      DIMENSION GEOG(3,NUMATM)
      DIMENSION COORD(3,NUMATM)
      f=0.529177

      ICYC=ICYC+1

      write(*,*)
      write(*,1)
  1   format(70('-'))
      write(*,'(30x,'' CYCLE '',I3)') ICYC
      write(*,1)

      IF(INDEX(KEYWRD,'DSCF').NE.0) MODE=0                                      
      IF(INDEX(KEYWRD,'TCSCF').NE.0) MODE=1   
      IF(INDEX(KEYWRD,'MOLPRO').NE.0) MODE=2                                    
      IF(INDEX(KEYWRD,'MOPAC').NE.0) MODE=3                                     
      IF(INDEX(KEYWRD,'ZINDO').NE.0) MODE=4                                     
      IF(INDEX(KEYWRD,'MP2').NE.0) MODE=5
      IF(INDEX(KEYWRD,'RIMP2').NE.0) MODE=6
      IF(INDEX(KEYWRD,'RIDFT').NE.0) MODE=7
      IF(INDEX(KEYWRD,'DIESEL').NE.0) MODE=8
	  
c	  if(MODE.EQ.8) write(*,*) 'switch to mode for diesel'

      kl=0
      do 210 i=1,natoms
      do 210 j=1,3
      kl=kl+1
      coord(j,i)=f*rfcrd(kl)
c     write(*,*) 'cfc coord(j,i)=',coord(j,i)
 210  continue

c     write(*,*)'natoms=',natoms
      WRITE(6,'(//10X,''CARTESIAN COORDINATES'',/)')
      WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''X'',
     .9X,''Y'',9X,''Z'',/)')
      DO 240 I=1,NATOMS
         WRITE(6,'(I6,8X,A2,4X,3F12.6)')
     .  I,ELEMNT(LABELS(I)),(COORD(J,I),J=1,3)
  240 CONTINUE


C output the coordinates in the proper format

      CALL OUT(1,MODE,COORD)

C calc energy and gradients

      ESCF=0.0D0


      CALL EGRAD(MODE,NUMAT,ESCF,GXYZ)
      WRITE(6,'(/10X,''TOTAL ENERGY (AU) = '',F15.7)') ESCF
      WRITE(6,'(/10X,''ENERGY CHANGE     = '',F15.7/)') ESCF-ELAST
      ELAST=ESCF

      IF(DABS(ESCF).LT.1.0D-10) STOP 'EXTERNAL PROGRAM ERROR'

      IF(MODE.EQ.2) CALL perm(natoms,labels,coord,gxyz)

      SUM=0.0D0
      DO I=1,NUMAT
         DO J=1,3
            SUM=SUM+GXYZ(J,I)**2
         ENDDO
      ENDDO

      CGRAD=DSQRT(SUM)

      WRITE(6,250) DSQRT(SUM)
  250 FORMAT(/,' NORM OF CARTESIAN GRADIENT = ',F9.2,/) 

      WRITE(6,'(/10X,''CARTESIAN GRADIENTS (KCAL/ANG)'',/)')
      WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''FX'',
     .9X,''FY'',9X,''FZ'',/)')

      K=0
      DO 241 I=1,NATOMS 
c     do 241 i=1,numat
         IF(LABELS(I).LT.99) THEN
         K=K+1
         WRITE(6,'(I6,8X,A2,4X,3F10.4)')
     .  I,ELEMNT(LABELS(I)),(GXYZ(J,k),J=1,3)
        ENDIF
  241 CONTINUE

      ij=0
      do 300 i=1,numat
	 do 300 j=1,3
	 ij=ij+1
	 grad(ij)=gxyz(j,i)
c        write(*,*)'grad(ij)=',grad(ij)
 300  continue

      escf=escf*627.5095410d0

      RETURN
      END

