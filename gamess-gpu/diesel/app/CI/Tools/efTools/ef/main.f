      PROGRAM MOPAC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      COMMON /KEYWRD/ KEYWRD
      CHARACTER*241 KEYWRD
      COMMON /OKMANY/ ISOK
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XPARAM(MAXPAR)
      COMMON /MESAGE/ IFLEPO,ISCF
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /CGOORD/ COORD(3,NUMATM)
      COMMON /ELEMTS/ ELEMNT(107)
      CHARACTER*2 ELEMNT              
      COMMON /LAST  / LAST
      COMMON /ENERG / ELAST
      COMMON /NUMCAL/ NUMCAL
      COMMON /TIMEXX  / TIME0
      COMMON /PATH  / LATOM,LPARAM,REACT(200)
      COMMON /PRTFLG/ IPRT
      LOGICAL ISOK, IPRT
      DIMENSION NAME(3)
      CHARACTER*10 NAME
      CHARACTER*80 A80 
      DATA NAME /'LENGTH','BEND','DIHEDRAL'/



      call system('date>startandstop')

	
      write(*,1)
  1   format(/,70('#'))
      print*
      print*, '                OPTIMIZATION BY EIGENVECTOR FOLLOWING'
      print*
      print*, '      REFERENCE: J. BAKER, J. COMP. CHEM. 7 (1985),385'
      print*, '      IMPLEMT. : S. GRIMME (NOV.-DEC. 1994)'
      print* 
      print*, '                       Add-Ons and modifications by'
      print*, '                       CM Marian and M Gastreich  '
      print*,' '
      print*, '                       cm & ghost@silly.thch.uni-bonn.de'
      print*,' '
      print*,' '
      print*,'  Tools:'
      print*,'  isovib: re-calculates numerically determined IR spectra'
      print*,'          (just call it without arguments)'
      print*,'  pse:    prints out the PSE and corresponding atomic '
      print*,'          masses we are working with'
      print*,' '
      write(*,1)

      write(*,*)
      write(*,*) 'IMPORTANT KEYWORDS'
      write(*,*)
      write(*,*) '0SCF    : ECHO INPUT ONLY'     
      write(*,*) '          COORDINATES FOR EXTERNAL PROGRAM ON <coord>'
      write(*,*) 'EF      : SEARCH MINIMUM '
      write(*,*) 'TS      : SEARCH TRANSITION STATE'
      write(*,*)
      write(*,*) 'HESSI   : (INITIAL) HESSIAN CALCULATED NUMERICALLY'
      write(*,*) '          IN INTERNAL COORDINATES'
      write(*,*) '          WRITTEN ON FILE <hessian>'
      write(*,*) 'RHESS   : INITIAL HESSIAN READ FROM <hessian>'   
      write(*,*)
      write(*,*) 'HESSC   : HESSIAN CALCULATED NUMERICALLY'
      write(*,*) '          IN CARTESIAN COORDINATES'
      write(*,*) '          WRITTEN ON FILE <hessianc>'
      write(*,*) '          WARNING: <hessianc> CANNOT BE USED AS START'
      write(*,*) '          HESSIAN FOR GEOMETRY OPTIMIZATION'
      write(*,*) 'INT     : DETERMINE NUMERICAL INFRA-RED INTENSITIES'
      write(*,*) '          - MP2 AND RIMP2 METHODS ONLY'
      write(*,*) '          - NON-LINEAR MOLECULES  ONLY'
      write(*,*) '          - VALUE OF DOUBLE (see below) WILL BE SET '
      write(*,*) '            EQUAL TO [1] TO MINIMIZE NUMERICAL ERRORS'
      write(*,*) '          - RESTARTS WILL FAIL'
c     ghost half the step size (cp. hessc.f)
      write(*,*) 'CSTEP=  : STEP FOR HESSC IN ANGSTROMS'
      write(*,*) '          (DEFAULT=0.010)'
      write(*,*) 'DIR     : START WITH ELONGATION OF DIRECTION XX'
      write(*,*) '          (DEF.=0) THIS OPTION MAKES ONLY SENSE IF'
      write(*,*) '          IT IS USED AS A RESTART OPTION FOR HESSC'
      write(*,*) '          XX IS THE LAST INTEGER OF FILE'
      write(*,*) '          <hesscinfo>'
      write(*,*) 'DOUBLE= : VALUES OF 0 (DEFAULT) AND 1 ARE ALLOWED'
      WRITE(*,*) '           0 MEANS ELONGATIONS IN JUST ONE DIRECTION;'
      WRITE(*,*) '           1 MEANS ELONGATIONS IN TWO  DIRECTIONS'
      write(*,*) '          - NOTE THAT IN THE CASE OF A RESTART, '
      write(*,*) '            DOUBLE MUST HAVE THE SAME VALUE AS'
      write(*,*) '            IN THE ORIGINAL RUN'
      write(*,*) 'ISOAV   : USE AVERAGE ISOTOPIC MASSES IN FREQ'
      write(*,*)
      write(*,*) 'DSCF    : CALL DSCF PROGRAM (DEFAULT)'                     
      write(*,*) 'TCSCF   : CALL DSCF-TC PROGRAM'                     
      write(*,*) 'MOLPRO  : CALL MOLPRO PROGRAM'                     
      write(*,*) 'MOPAC   : CALL MOPAC  PROGRAM'                     
      write(*,*) 'ZINDO   : CALL ZINDO PROGRAM'                     
      write(*,*) 'MP2     : CALL MP2 PROGRAM'                     
      write(*,*) 'RIMP2   : CALL RIMP2 PROGRAM'                     
      write(*,*) 'RIDFT   : CALL RIDFT PROGRAM'                     
      write(*,*)
      write(*,*) 'MODE=N  : FOLLOW MODE N IN TS SEARCH'              
      write(*,*) 'GNORM=X : STOP IF GRADIENT NORM BELOW X (DEF=1.0)'           
      write(*,*) 'DMAX=X  : MAXIMUM INT COORD CHANGE (ANG/RAD, DEF=0.2)' 
      write(*,*) 'XZ      : PERMUTE X AND Z COORDIANTES'           
      write(*,*) 'YZ      : PERMUTE Y AND Z COORDIANTES'           
      write(*,*) 'XY      : PERMUTE X AND Y COORDIANTES'           
      write(*,*) 'NOCEMA  : DO NOT TRANSFORM TO CENTER OF MASS'    
      write(*,*) 'MAXCYCLE: MAXIMUM NUMBER OF OPT. CYCLES (>1;DEF=25)'
      write(*,*)
      write(*,*) 'IROOT=N : OPTIMIZE CI ROOT N (MRDCI, MOLPRO, DEF=1)'           
      write(*,*)
      
      write(*,*)
      write(*,*) 'IMPORTANT FILES'
      write(*,*)
      write(*,*) 'INPUT, OUTPUT : STANDARD IO'
      write(*,*) '<intcoord>    : ACTUAL (OPTIMIZED) GEOMETRY'
      write(*,*) '<hessian>     : ACTUAL HESSIAN IN INTERNAL COORDINATE'
      write(*,*) '<hessianc>    : HESSIAN IN CARTESIAN COORDINATES'
      write(*,*) '<hesscinfo>   : NEEDED IN THE CASE OF A RESTART'
      write(*,*) '                CALCULATION OF HESSC'
      write(*,*) '<coord>       : COORDINTES FOR EXTERNAL PROGRAM'
      write(*,*) '<energy>      : OPTIMIZATION STATISTICS'          
      write(*,*)
      write(*,*) 'CARTESIAN COORDINATE SYSTEM CONVENTION'
      write(*,*)
      write(*,*) 'ATOM 1 AND 2 DEFINE THE Z-AXES'
      write(*,*) 'ATOM 1,2 AND 3 DEFINE THE YZ-PLANE'
      write(*,*)
      write(*,*) 'CARTESIAN COORDINATE SYSTEM CONVENTION (WITH NOCEMA)'
      write(*,*)
      write(*,*) 'ATOM 1 AND 2 DEFINE THE X-AXES'
      write(*,*) 'ATOM 1,2 AND 3 DEFINE THE XY-PLANE'

      ELAST=0.0D0
      NUMCAL=0
      ISOK=.TRUE.
      IPRT=.TRUE.
   10 NUMCAL=NUMCAL+1
C
      TIME0=SECOND()
C
C READ AND CHECK INPUT FILE, EXIT IF NECESSARY.
C     WRITE INPUT FILE TO UNIT 6 AS FEEDBACK TO USER
C
      CALL READMO

      MODE=0
      IF(INDEX(KEYWRD,'TCSCF').NE.0) MODE=1
      IF(INDEX(KEYWRD,'MOLPRO').NE.0) MODE=2
      IF(INDEX(KEYWRD,'MOPAC').NE.0) MODE=3
      IF(INDEX(KEYWRD,'ZINDO').NE.0) MODE=4
      IF(INDEX(KEYWRD,'MP2').NE.0) MODE=5

c     ghost IF(INDEX(KEYWRD,'RIMP2').NE.0) MODE=6
      IF(INDEX(KEYWRD,'RIMP2').NE.0) then
         mode=6
C        write(*,*) 'starting at MODE:', mode
      endif
      
      IF(INDEX(KEYWRD,'RIDFT').NE.0) MODE=7
      IF(INDEX(KEYWRD,'DIESEL').NE.0) MODE=8

      CALL OUT(1,MODE,COORD)
c     ghost
      write(*,*) '  MODE control:', mode

      write(*,'(//,5X,''CARTESIAN COORDINATES FOR EXTERNAL CODE'',/)')
      CALL OUT(6,MODE,COORD)
c     ghost
      write(*,*) '  MODE control', mode

      L=0
      DO I=1,NATOMS 
         IF(LABELS(I).NE.99) L=L+1
      ENDDO
      NUMAT=L

      write(*,'(//,'' NUMBER OF ATOMS (INCL. DUMMY) = '',I3)') NATOMS
      write(*,'(   '' NUMBER OF ATOMS               = '',I3)') NUMAT  
      write(*,'(   '' NUMBER OF GEOMETRIC VARIABLES = '',I3)') NVAR  
      write(*,'(/  '' GEOMETRIC VARIABLES (ANG/RAD)'',/)')
 
      do i=1,nvar
         write(*,'(I3,5X,A10,F12.6)')I,NAME(LOC(2,I)),xparam(i)
      enddo
      write(*,*)

      IF(INDEX(KEYWRD,'0SCF').NE.0) GOTO 99                       

      IF(INDEX(KEYWRD,'GRAD').NE.0) THEN
         CALL COMPFG(XPARAM,GRAD)
         GOTO 99
      ENDIF

      IF(INDEX(KEYWRD,'HESSI').NE.0) CALL HESSIAN(XPARAM)        

      IF(INDEX(KEYWRD,'DYNAM').NE.0) THEN
         CALL IRC(MODE,NUMAT,XPARAM)        
         GOTO 99
      ENDIF
 
      IF(INDEX(KEYWRD,'EF').NE.0.OR.
     .   INDEX(KEYWRD,'TS').NE.0) THEN
         NL=0
         OPEN(UNIT=99,FILE='energie')
 990     READ(99,'(A)',END=991) A80
         NL=NL+1
         GOTO 990
 991     CONTINUE
         IF(NL.EQ.0) WRITE(99,'(7X,''CYCLE'',6X,''ENERGY'',
     .                      6X,''CART GRAD'',6X,''INT GRAD'')')
         CALL EF(XPARAM,COORD,NVAR)        
         CLOSE (99)
         GOTO 99
      ENDIF

      IF(INDEX(KEYWRD,'HESSC').NE.0) then
         CALL HESSIANC(coord,numat)
c     ghost
         write(*,*) '  MODE final value:', mode
         write(*,*) '  good accuracy can be obtained with'
         write(*,*) '  scfconv=8 or denconv=.1e-6 only ! '
      ENDIF

99    WRITE(*,*) 
      call system('date>>startandstop')
      WRITE(*,*) '***** ALL DONE *****'
      STOP 
      END
