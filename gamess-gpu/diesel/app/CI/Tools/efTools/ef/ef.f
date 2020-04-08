      SUBROUTINE EF(XPARAM, COORD, NVAR)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      DIMENSION XPARAM(*), COORD(3,*)                                      
**********************************************************************          
*                                                                               
*   BAKER IS A QUASI NEWTON RAPHSON OPTIMIZATION ROUTINE BASED ON               
*         JON BAKERS ALGORITHM (J.COMP.CHEM. 7, 385).                           
*                                                                               
*  ON ENTRY XPARAM = VALUES OF PARAMETERS TO BE OPTIMISED.                      
*           NVAR   = NUMBER OF PARAMETERS TO BE OPTIMISED.                      
*                                                                               
*  ON EXIT  XPARAM = OPTIMISED PARAMETERS.                                      
*           FUNCT  = HEAT OF FORMATION IN KCAL/MOL.                             
*                                                                               
*                    IMPLEMENTED IN MOPAC AUGUST 1990 BY JJPS                   
*                                                                               
**********************************************************************          
C                                                                               
C   THE BAKER ROUTINE IS CONTROLED BY THE FOLLOWING KEYWORDS :                  
C   EF,TS           : SEARCH FOR MINIMUM OR TRANSITION STRUCTURE                
C   CYCLES          : MAXIMUM NUMBER OF ITERATIONS IN OPTIMIZATION              
C   RECALC=N        : NUMBER OF CYCLES BETWEEN HESSIAN CALCULATION              
C   DMAX=N.N        : MAXIMUM STEPSIZE IN ANG OR RAD                            
C   IUPD=0,1,2      : HESSIAN UPDATING ; NONE, POWELL, BFGS                     
C   HESS=0,1,2,3    : START HESSIAN ; DIAGONAL, CALCULATE,                      
C                     READ FROM DISK, READ FROM INPUT.                          
C                     NOTE THAT CALCULATION OF SELECTIVE ROWS (COLOUMNS)        
C                     OF THE HESSIAN COMBINED WITH DIAGONAL GUESS IS            
C                     CONTROLED BY SETTING OPTIMIZATION PARAMETERS = 2.         
C                     IF ANY OPT PARAMETERS = 2, THEN ABOVE HESS OPTION         
C                     IS OVERRIDDEN.                                            
C   MODE=N          : FOLLOW THE N'TH LOWEST HESSIAN MODE TOWARD TS             
C   NONR            : DO NOT USE NEWTON-RAPHSON STEP IF FORMD FAILS AND         
C                     MAX GRADIENT IS LESS THAN RCUT.                           
C   EIGINV          : USE HESSIAN EIGENVALUE REVERSION RATHER THAN              
C                     DIAGONAL SHIFT IF FORMD FAILS AND NONR OR GNORM           
C                     IS LARGER THAN RCUT.                                      
C   ORIDE           : USE WHATEVER LAMDA(S) FORMD PRODUCES, EVEN IF THEY        
C                     ARE 'UNACCEPTABLE'.                                       
C   GNORM=N.N       : EXIT WHEN LARGEST GRADIENT IS LESS THAN N.N               
C   PRNT=0-5        : 0 : PRINT SUMMERY EVERY CYCLE                             
C                     1 : + GEOMETRY AND TS MODE                                
C                     2 : + GRADIENTS                                           
C                     3 : + HESSIAN EIGENVALUES AND -VECTORS                    
C                     4 : + HESSIAN MATRIX                                      
C                     5 : MAXIMUM PRINTING = DEBUG PRINTING                     
C                                                                               
C   DEFAULTS :                                                                  
C   PRNT     : 0                                                                
C   CYCLES   : 100                                                              
C   DMAX     : 0.20                                                             
C   RECALC   : 999999                                                           
C   IUPD     : 2 FOR MINIMUM SEARCH, 1 FOR TS SEARCH                            
C   HESS     : 0 FOR MINIMUM SEARCH, 1 FOR TS SEARCH                            
C   MODE     : 1, ONLY USED IN TS SEARCH                                        
C   GNORM    : 0.4 OR 0.01 IF "PRECISE" IS SPECIFIED                            
C   NONR     : FALSE                                                            
C   EIGINV   : FALSE                                                            
C   ORIDE    : FALSE                                                            
C                                                                               
C   HINTS AND COMMENTS :                                                        
C   NORMALLY THE DEFAULT OPTIONS ARE SUFFICIENT. NONR, EIGINV AND               
C      ORIDE OPTIONS ARE SPECIAL AND SHOULD RARELY BE TOUCHED.                  
C   CONVERGENCE CRITERIA IS THE SAME AS FOR "SIGMA" OPTIMIZATION,               
C      NAMELY THAT ALL GRADIENTS MUST BE BELOW 0.4 OR 0.01 IF                   
C      "PRECISE" IS SPECIFIED. RESULTS OBTAINED WITHOUT THE "PRECISE"           
C      OPTION TURNED ON MAY BE UNRELIABLE DUE TO ACCUMULATED ERRORS             
C      FROM VARIOUS PARTS OF THE PROGRAM (SEE E.G. J.COMP.CHEM. 9,X)            
C      BUT CAN BE USED IN PRELIMINARY STAGES. ONLY RESULTS FROM                 
C      PRECISE RUNS SHOULD EVER BE PUBLISHED!                                   
C      WITH THE PRECISE KEYWORD GEOMETRIES ARE USUALLY CONVERGED TO             
C      0.0001 IN BOND LENGTH, 0.01 IN BOND ANGLES AND 0.1 IN TORSIONAL          
C      ANGLES. ENERGIES ARE CONVERGED TO AT LEAST 0.001 KCAL/MOL.               
C   CONVERGENCE CAN BE MONITORED IN TWO WAYS : THE GRADIENT NORM SHOULD         
C      GO TOWARD ZERO, AND THE LAMDA VALUES SHOULD ALSO GO TO ZERO.             
C   HESSIAN AND GEOMETRY WILL BE SAVED ON DISK (AND CAN BE RESTARTED) :         
C      1) AFTER HESSIAN CALCULATION (PARTLY OR WHOLE)                           
C      2) IN AN OPTIMIZATION EVERY X MIN. CPU USED (MACHINE DEPENDENT)          
C      3) BEFORE TERMINATION, EITHER DUE TO CONVERGENCE, TIME OR CYCLES         
C   IF HESS=2 THEN THE GEOMETRY IS TAKEN FROM THE INPUT WHILE THE HESSIA        
C      IS READ FROM DISK. UPON RESTART BOTH GEOMETRY AND HESSIAN IS             
C      AUTOMATICALLY PICKED UP FROM DISK.                                       
C   TS SEARCH REQUIRES HESS=1, 2 OR 3, IT CANNOT WORK WITH HESS=0               
C   HESS=1 FOR BAKERMIN IS RECOMMENDED FOR DIFFICULT CASES (E.G. CYCLIC)        
C   HESS=3 WILL READ THE HESSIAN FROM INPUT, MOSTLY USED FOR TESTING.           
C      TO READ IT IN, PLACE THE WORD 'HESSIAN' SOMEWHERE IN THE INPUT           
C      FILE IN FORMAT (A7), FOLLOWED BY THE HESSIAN ELEMENTS IN                 
C      FORMAT (8F10.2).                                                         
C   THE ROUTINE (ESPECIALLY TS SEARCH) WORKS BEST IN INTERNAL COORDINATE        
C   MOST COMMON ERROR FROM BAKER IS "UNABLE TO DETERMINE LAMDA IN               
C      FORMD". THIS WILL OFTEN HAPPENS IF :                                     
C      1) SYMMETRY CONSTRAINED VARIABLES ARE MARKED FOR OPTIMIZATION            
C      2) VARIABLE DETERMINED BY SYMMETRY IS MARKED FOR OPTIMIZATION            
C      3) MORE VARIABLES THAN INTERNAL DEGREES OF FREEDOM SPECIFIED             
C      4) TOO LARGE STEP SIZE HAS BEEN CHOSEN, UPDATING GOES HAYWIRE            
C      5) BAD HESSIAN IS USED AS START HESSIAN                                  
C   EVEN WHEN THESE REASONS HAVE BEEN ELIMINATED, THE FAIL OCCURS QUITE         
C      FREQUENTLY IN THE ***** ENVIROMENT, AND USUALLY IF SOFT MODES            
C      SUCH AS METHYL ROTATIONS ARE PRESENT. SOME PRELIMINARY PATCHES           
C      FOR CIRCUMVENTING THESE FAILURES HAVE BEEN IMPLEMENTED IN THIS           
C      VERSION. IF THE STRUCTURE IS CLOSE TO OPTIMIZED (MAX GRAD.LT.RCUT        
C      = 10), THEN SIMPLY TAKE THE NEWTON-RAPHSON STEP, UNLESS NONR IS          
C      TURNED ON. IF MAX GRAD.GT.RCUT OR NONR IS ON, THEN USE EITHER            
C      A HESSIAN SHIFT OR A HESSIAN EIGENVALUE INVERSION.                       
C      THE DEFAULT, AND ONLY OPTION FOR MINIMUM SEARCHES, IS TO                 
C      INCREASE THE HESSIAN DIAGONAL WITH LOWEST HESSIAN EIGENVALUE             
C      (SECOND LOWEST EIGENVALUE IF TS SEARCH) + 2*EIGMIN.                      
C      THIS SHOULD ALWAYS WORK FOR A MIN SEARCH, BUT MAY GO WRONG IN A          
C      TS SEARCH. IF THIS HAPPENS ONE CAN TRY TO TURN ON THE EIGINV             
C      OPTION. THIS WILL REVERSE THE SIGN OF ALL NEGATIVE EIGENVALUES           
C      OTHER THAN THE TS MODE BEING FOLLOWED. AGAIN THERE IS NO GARANTI         
C      THAT THIS IS WILL LEAD TO A 'BETTER' GEOMETRY.                           
C      FINALLY AN OVERRIDE SWITCH IS IMPLEMENTED. IF KEYWORD ORIDE IS           
C      USED THEN WHATEVER LAMDA'S FORMD PRODUCES WILL BE USED,                  
C      EVEN IF THEY ARE 'UNACCEPTABLE', AS DISCUSSED BY SIMONS.                 
C   THESE PATCHES ARE AS MENTIONED STILL EXPERIMENTAL, BUT HAVE                 
C      BEEN FOUND TO WORK IN MOST CASES.                                        
C   INCREASING DMAX CAN LEAD TO FASTER CONVERGENCE BUT CAN ALSO MAKE            
C      THE OPTIMIZATION GO BAD VERY FAST. FURTHERMORE THE HESSIAN               
C      UPDATING DETERIORATE WHEN USING LARGE STEPSIZES. REDUCING THE            
C      STEPSIZE TO 0.10 OR 0.05 IS RECOMMENDED WHEN ENCOUNTERING                
C      CONVERGENCE PROBLEMS.                                                    
C   RECALC=N CALCULATES THE HESSIAN EVERY N STEPS IN THE OPTIMIZATION.          
C      FOR SMALL N THIS IS COSTLY BUT IS ALSO VERY EFFECTIVE IN                 
C      TERMS OF CONVERGENCE. RECALC=10 AND DMAX=0.10 CAN BE USEFUL              
C      FOR DIFFICULT CASES. IN EXTREME CASES RECALC=1 AND DMAX=0.05             
C      WILL ALWAYS FIND A STATIONARY POINT, IF IT EXIST.                        
C   IF READING IN GNORM=N.N NOTICE THAT IT IS MAX GRADIENT THRESHOLD            
C      THAT IS BEING MODIFIED, NOT GRADIENT NORM. IGNORE THE MESSAGE            
C      ABOUT "FLEPO EXIT WHEN...".                                              
C   IUPD SHOULD VERY RARELY BE TOUCHED. IUPD=1 CAN BE USED IN MINIMUM           
C      SEARCHES IF THE THE MESSAGE "HEREDITARY POSITIVE DEFINITENESS            
C      ENDANGERED. UPDATE SKIPPED THIS CYCLE" OCCURS EVERY CYCLE FOR            
C      10-20 ITERATIONS. NEVER USE IUPD=2 FOR TS SEARCH!                        
C                                                                               
C   **** NOTICE THAT VERY LITTLE ERROR CHECKING IS DONE.                        
C   **** ONLY RUNS REQUESTING INCOMPATIBLE OPTIONS ARE TERMINATED.              
C   **** A MINIMUM KNOWLEDGE OF THE THEORY BEHIND THE ROUTINE                   
C   **** IS RECOMMENDED TO GET THE OPTIMUM USAGE OUT OF IT.                     
C                                                                               
      COMMON /ENERG / ELAST, CGRAD, XGRAD                                       
      COMMON /MESAGE/ IFLEPO,ISCF                                               
      COMMON /GEOVAR/ NDUM,LOC(2,MAXPAR), IDUMY, XARAM(MAXPAR)                  
      COMMON /GEOM  / GEO(3,NUMATM)                                             
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)         
      COMMON /ISTOPE/ AMS(107)                                                  
      COMMON /LAST  / LAST                                                      
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /TIMEXX  / TIME0                                                     
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA                                       
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),          
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,                
     2                NCLOSE,NOPEN,NDUMY,FRACT                                  
      COMMON /NUMCAL/ NUMCAL                                                    
      COMMON /TIMDMP/ TLEFT, TDUMP                                              
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GMIN1(MAXPAR)                             
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),                  
     1PMAT(MAXPAR**2)                                                           
      COMMON /SCRACH/ PVEC                                                      
      COMMON /SCFTYP/ EMIN, LIMSCF                                              
      COMMON /FMATRX/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),                      
     1U(MAXPAR,MAXPAR),DD,MODE,NSTEP,NEGREQ,IPRNT                               
      DIMENSION IPOW(9), EIGVAL(MAXPAR),TVEC(MAXPAR),SVEC(MAXPAR),              
     1FX(MAXPAR),HESSC(MAXPAR**2),UC(MAXPAR**2), PVEC(MAXPAR*MAXPAR)            
      LOGICAL RESTRT,SCF1,FAIL,HSSHFT,EIGINV,NR,ORIDE                           
      LOGICAL  LIMSCF, LOG                                                      
      DIMENSION NAME(3), XX(20)                                          
      CHARACTER*10 NAME                                                         
      CHARACTER KEYWRD*241                                                      
      CHARACTER*80 A80                                                          
      CHARACTER HESWRD*7,WRD*7                                                  
      EQUIVALENCE(IPOW(1),IHESS)
C, (UC,U)                                        
      DATA HESWRD   /'HESSIAN'/                                                 
      DATA  ICALCN,ZZERO,ONE,TWO    /0,0.D0,1.D0,2.D0/                          
      DATA EIGMIN,RCUT /1.D-1,1.D+01/                                           
      DATA NAME /'LENGTH','BEND','DIHEDRAL'/                                    
                                                                                
C     DGHSX IS HESSIAN DIAGONAL FOR HESS=0. FOR STRETCHING, ANGLE,              
C     DIHEDRAL.  THE VALUES GIVEN BELOW ARE 'OPTIMUM' FOR CYCLOHEXANONE         
                                                                                
      DGHSS=2500.D0                                                             
      DGHSA= 900.D0                                                             
      DGHSD=1000.D0                                                             
      TLEFT=18000.D0                                                            
                                                                                
C     GET ALL INITIALIZATION DATA                                               
      IF(ICALCN.NE.NUMCAL) THEN                                                 
         NVAR=ABS(NVAR)                                                         
         LDUMP=0                                                                
         ICALCN=NUMCAL                                                          
         NSTEP=0                                                                
         IHESS=0                                                                
         MODE =0                                                                
         LAST=0                                                                 
         NTIME=0                                                                
         ILOOP=1                                                                
         IMIN=INDEX(KEYWRD,' EF')                                               
         IF(IMIN.NE.0) THEN                                                     
            WRITE(6,'(//,5X,''SEARCHING FOR A MINIMUM'',//)')                   
            IGTHES=0                                                            
            IUPD  =2                                                            
            WRITE(6,'(5X,''BFGS UPDATE USED '',/)')                             
            NEGREQ=0                                                            
         ENDIF                                                                  
         LIMSCF=.TRUE.                                                          
         ITS=INDEX(KEYWRD,' TS')                                                
         IF(ITS.NE.0) THEN                                                      
            WRITE(6,'(//,5X,''SEARCHING FOR A TRANSITION STATE'',//)')          
            IGTHES=1                                                            
            IUPD  =1                                                            
            WRITE(6,'(5X,''POWELL UPDATE USED '',/)')                           
            NEGREQ=1                                                            
         ENDIF                                                                  
         I=INDEX(KEYWRD,' DEBUG')*INDEX(KEYWRD,' BAKER')                        
         IPRNT=0                                                                
         IF(I.NE.0) IPRNT=2                                                     
         IRECLC=999999                                                          
         I=INDEX(KEYWRD,' RECALC=')                                             
         IF(I.NE.0) IRECLC=READA(KEYWRD,I)                                      
         I=INDEX(KEYWRD,' IUPD=')                                               
         IF(I.NE.0) IUPD=READA(KEYWRD,I)                                        
         I=INDEX(KEYWRD,' MODE=')                                               
         IF(I.NE.0) MODE=READA(KEYWRD,I)                                        
         DMAX=0.2                                                               
         I=INDEX(KEYWRD,' DMAX=')                                               
         IF(I.NE.0) DMAX=READA(KEYWRD,I)                                        
cmmg1                                                                           
         maxcycle=35                                                            
         I=INDEX(KEYWRD,' MAXCYCLE=')                                           
         IF(I.NE.0) MAXCYCLE=READA(KEYWRD,I)                                    
cmmg2                                                                           
         TOL2=1.0D0                                                             
         IF(INDEX(KEYWRD,' PREC') .NE. 0) TOL2=1.D-2                            
         I=INDEX(KEYWRD,' GNORM=')                                              
         IF(I.NE.0) TOL2=READA(KEYWRD,I)                                        
         IF(INDEX(KEYWRD,' LET').EQ.0.AND.TOL2.LT.0.01D0)THEN                   
            WRITE(6,'(/,A)')'  GNORM HAS BEEN SET TOO LOW, RESET TO 0           
     1.01'                                                                      
            TOL2=0.01D0                                                         
         ENDIF                                                                  
         HSSHFT=.TRUE.                                                          
         NR=.TRUE.                                                              
         IF(INDEX(KEYWRD,'NONR') .NE. 0)NR=.FALSE.                              
         EIGINV=.FALSE.                                                         
         EIGINV=(INDEX(KEYWRD,'EIGINV') .NE. 0)                                 
         IF(EIGINV) THEN                                                        
            HSSHFT=.FALSE.                                                      
            IF(IPRNT.EQ.2)THEN                                                  
               WRITE(6,*)' '                                                    
               WRITE(6,*)'    HESSIAN EIGENVALUE INVERSION TURNED ON'           
            ENDIF                                                               
         ENDIF                                                                  
         I=INDEX(KEYWRD,' HESS=')                                               
         IF(I.NE.0) IGTHES=READA(KEYWRD,I)                                      
         IF(INDEX(KEYWRD,'RHESS') .NE. 0) IGTHES=2                              
                                                                                
         TIME1=TIME0                                                            
         TIME2=TIME1                                                            
                                                                                
C   DONE WITH ALL INITIALIZING STUFF.                                           
C   CHECK THAT OPTIONS REQUESTED ARE RESONABLE                                  
                                                                                
         IF(NVAR.GT.(3*NUMATM-6))WRITE(6,10)                                    
   10    FORMAT(/,'*** WARNING! MORE VARIABLES THAN DEGREES OF FREEDOM',        
     1/)                                                                        
         IF((ITS.NE.0).AND.(IUPD.EQ.2))THEN                                     
            WRITE(6,*)' TS SEARCH AND BFGS UPDATE WILL NOT WORK'                
            STOP                                                                
         ENDIF                                                                  
         IF((ITS.NE.0).AND.(IGTHES.EQ.0))THEN                                   
            WRITE(6,*)' TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN'         
            STOP                                                                
         ENDIF                                                                  
         IF((IMIN.NE.0).AND.(EIGINV))THEN                                       
            WRITE(6,*)' MIN SEARCH AND EIGINV WILL NOT WORK'                    
            STOP                                                                
         ENDIF                                                                  
         IF((IGTHES.LT.-1).OR.(IGTHES.GT.3))THEN                                
            WRITE(6,*)' UNRECOGNIZED HESS OPTION',IGTHES                        
            STOP                                                                
         ENDIF                                                                  
      ENDIF                                                                     
                                                                                
C     GET THE HESSIAN. DEPENDING ON IGTHES WE GET IT FROM :                     
C                                                                               
C     -1 : DIAGONAL MATRIX, DGHSX*I, BUT ROWS AND COLOUMNS OF                   
C          HESSIAN CALCULATED ACCORDING TO VARIABLES MARKED '2'                 
C      0 : DIAGONAL MATRIX, DGHSX*I (DEFAULT FOR MIN-SEARCH)                    
C      1 : CALCULATE IT NUMERICALLY (DEFAULT FOR TS-SEARCH)                     
C      2 : READ IN FROM FTN009                                                  
C      3 : READ IT IN FROM INPUT FILE (FTN005)                                  
C      4 : READ IN FROM FTN009 (DURING RESTART, PARTLY OR WHOLE,                
C          ALREADY DONE AT THIS POINT)                                          
C                                                                               
      IF (IGTHES.EQ.0) THEN                                                     
C        write(6,60)                                                            
C  60    FORMAT(/,10X,'DIAGONAL MATRIX USED AS START HESSIAN',/)                
         call esthess(nvar,xparam,coord,hess)                                   
c        DO 70 I=1,NVAR                                                         
c           DO 70 J=1,NVAR                                                      
c              HESS(I,J)=ZZERO                                                  
c  70    CONTINUE                                                               
c        DO 80 IJ=1,NVAR                                                        
c              I=LOC(2,IJ)                                                      
c              IF (I.EQ.1)HESS(IJ,IJ)=DGHSS                                     
c              IF (I.EQ.2)HESS(IJ,IJ)=DGHSA                                     
c              IF (I.EQ.3)HESS(IJ,IJ)=DGHSD                                     
c  80    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF (IGTHES.EQ.1.OR.IGTHES.EQ.2) THEN                                      
         WRITE(6,100)                                                           
  100    FORMAT(/,10X,'HESSIAN READ FROM <hessian>',/)                          
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      read (4,'(A)') a80                                                        
      call readl(80,a80,xx,nn)                                                  
      if(idint(xx(1)).ne.nvar) stop 'WRONG HESSIAN'                             
      read (4,'(A)') A80                                                        
      write(6,*) A80                                                            
      do i=1,nvar                                                               
         read (4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 
      DO 220 I=1,NVAR                                                           
         DO 220 J=1,I                                                           
            HESS(J,I)=HESS(I,J)                                                 
  220 CONTINUE                                                                  
      ENDIF                                                                     
                                                                                
C   NOT A RESTART, WE NEED TO GET THE GRADIENTS                                 
                                                                                
      DO 20 I=1,NVAR                                                            
   20       GRAD(I)=ZZERO                                                       
      CALL COMPFG(XPARAM, GRAD)                  
                                                                                
      IF(ITS.EQ.0)THEN                                                          
C     CHECK THAT GEOMETRY IS NOT ALREADY OPTIMIZED                              
         RMX=ZZERO                                                              
         DO 40 I=1,NVAR                                                         
            IF (ABS(GRAD(I)).GT.RMX) RMX=ABS(GRAD(I))                           
   40    CONTINUE                                                               
         IF (RMX.LT.TOL2) THEN                                                  
            IFLEPO=2                                                            
            LAST=1                                                              
            write(*,*)' '
            write(*,*)'This is an already optimized geometry'
            write(*,*)' '
            RETURN                                                              
         ENDIF                                                                  
      ENDIF                                                                     
                                                                                
      ICYCLE=0                                                                  
cmmg  MAXCYCLE=25                                                               
                                                                                
                                                                                
                                                                                
  230 CONTINUE                                                                  
C     START OF MAIN LOOP                                                        
C     WE NOW HAVE GRADIENTS AND A HESSIAN. IF THIS IS THE FIRST                 
C     TIME THROUGH DON'T UPDATE THE HESSIAN. FOR LATER LOOPS ALSO               
C     CHECK IF WE NEED TO RECALCULATE THE HESSIAN                               
C     TIME2=SECOND()                                                            
C     TSTEP=TIME2-TIME1                                                         
C     TLEFT=TLEFT-TSTEP                                                         
C     TIME1=TIME2                                                               
c     IFLEPO=0                                                                  
                                                                                
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
                                                                                
                                                                                
  240 CONTINUE                                                                  
      ICYCLE=ICYCLE+1                                                           
      IF (IHESS.GE.IRECLC.AND.IFLEPO.NE.15) THEN                                
         ILOOP=1                                                                
         IHESS=0                                                                
         IGTHES=1                                                               
      ELSE                                                                      
        IF(IFLEPO.EQ.15)GOTO 270                                                
C                                                                               
C        PRINT RESULTS IN CYCLE                                                 
         GNFINA=SQRT(DOT(GRAD,GRAD,NVAR))                                       
C        WRITE(6,250) GNFINA,TOL2                                               
C 250    FORMAT(/,' NORM OF INTERNAL GRADIENT = ',F9.2,                         
C    .   5X,' (THRESHOLD =',F6.2,')')                                           
                                                                                
         WRITE(99,'(5X,I4,5X,F12.6,5X,F8.3,5X,F8.3)')                           
     .   ICYCLE,ELAST,CGRAD,GNFINA                                              
C                                                                               
C        TEST FOR CONVERGENCE                                                   
C                                                                               
         RMX=SQRT(DOT(GRAD,GRAD,NVAR))                                          
         IF (RMX.LT.TOL2)GOTO 580                                               
         IF (ICYCLE.GE.MAXCYCLE) THEN                                           
             WRITE(6,'(//5X,'' MAXCYCLE EXCEEDED '',I4)') MAXCYCLE              
             GOTO 580                                                           
         ENDIF                                                                  
  270    CONTINUE                                                               
         IHESS=IHESS+1                                                          
         NSTEP=NSTEP+1                                                          
         IF (IHESS.GT.1) CALL UPDHES(SVEC,TVEC,NVAR,IUPD)                       
      ENDIF                                                                     
C                                                                               
      NHSINC=0                                                                  
  280 IJ=0                                                                      
      DO 290 I=1,NVAR                                                           
         DO 290 J=1,I                                                           
            IJ=IJ+1                                                             
            HESSC(IJ)=HESS(J,I)                                                 
  290 CONTINUE                                                                  
      CALL HQRII(HESSC,NVAR,NVAR,EIGVAL,UC)                                     
      DO 300 I=NVAR,1,-1                                                        
      DO 300 J=NVAR,1,-1                                                        
  300 U(J,I)=UC(J+(I-1)*NVAR)                                                   
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      write(4,*) 'hessian dim=',nvar                                            
      write(4,*) 'lowest eigenvalue=',eigval(1)                                 
      do i=1,nvar                                                               
         write(4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 
      IF (IPRNT.GE.3) THEN                                                      
         IF (IPRNT.GE.4) THEN                                                   
            WRITE(6,*)' '                                                       
            WRITE(6,*)'              HESSIAN MATRIX'                            
            LOW=1                                                               
            NUP=8                                                               
  310       NUP=MIN(NUP,NVAR)                                                   
            WRITE(6,330) (I,I=LOW,NUP)                                          
            DO 320 I=1,NVAR                                                     
               WRITE(6,340) I,(HESS(I,J),J=LOW,NUP)                             
  320       CONTINUE                                                            
            NUP=NUP+8                                                           
            LOW=LOW+8                                                           
            IF(LOW.LT.NVAR) GOTO 310                                            
  330       FORMAT(/,3X,8I9)                                                    
  340       FORMAT(1X,I3,8F9.1)                                                 
  350       FORMAT(/,4X,8F9.1,/)                                                
  360       FORMAT(1X,I3,8F9.4)                                                 
         ENDIF                                                                  
         WRITE(6,*)' '                                                          
         WRITE(6,*)'              HESSIAN EIGENVALUES AND -VECTORS'             
         LOW=1                                                                  
         NUP=8                                                                  
  370    NUP=MIN(NUP,NVAR)                                                      
         WRITE(6,330) (I,I=LOW,NUP)                                             
         WRITE(6,350) (EIGVAL(I),I=LOW,NUP)                                     
         DO 380 I=1,NVAR                                                        
            WRITE(6,360) I,(U(I,J),J=LOW,NUP)                                   
  380    CONTINUE                                                               
         NUP=NUP+8                                                              
         LOW=LOW+8                                                              
         IF(LOW.LT.NVAR) GOTO 370                                               
      ENDIF                                                                     
      DO 390 I=1,NVAR                                                           
         IF (ABS(EIGVAL(I)).LT.EIGMIN)THEN                                      
            IF (EIGVAL(I).LT.ZZERO)EIGVAL(I)=-EIGMIN                            
            IF (EIGVAL(I).GT.ZZERO)EIGVAL(I)=EIGMIN                             
            WRITE(6,400)I,EIGVAL(I)                                             
         ENDIF                                                                  
  390 CONTINUE                                                                  
  400 FORMAT(5X,'WARNING! EIGENVALUE',I3,' TOO SMALL, REPLACED BY',             
     1F7.3)                                                                     
      NEG=0                                                                     
      DO 410 I=1,NVAR                                                           
         IF (EIGVAL(I) .LT. ZZERO)NEG=NEG+1                                     
  410 CONTINUE                                                                  
      IF (RMX.LT.TOL2)GOTO 580                                                  
      IF (IHESS.EQ.1.AND.IPRNT.EQ.2)WRITE(6,420)NEG                             
  420 FORMAT(/,10X,'HESSIAN HAS',I3,' NEGATIVE EIGENVALUE(S)',/)                
      DO 430 I=1,NVAR                                                           
         FX(I)=DOT(U(1,I),GRAD,NVAR)                                            
  430 CONTINUE                                                                  
C  TAKE THE P-RFO STEP FOR A TS SEARCH AND                                      
C  THE SIMPLE RFO STEP FOR A MINIMUM SEARCH                                     
C                                                                               
      IF(IPRNT.GE.5.AND.NEG.EQ.NEGREQ.AND.NEGREQ.EQ.1) WRITE(6,690)             
      IF(IPRNT.GE.5.AND.NEG.EQ.NEGREQ.AND.NEGREQ.EQ.0) WRITE(6,700)             
      IF(IPRNT.GE.5.AND.NEG.NE.NEGREQ.AND.NEGREQ.EQ.1) WRITE(6,710)             
      IF(IPRNT.GE.5.AND.NEG.NE.NEGREQ.AND.NEGREQ.EQ.0) WRITE(6,720)             
  440 CALL FORMD(EIGVAL,FX,NVAR,FAIL,ORIDE)                                     
C     IF FAIL IS TRUE THEN NO ACCEPTABLE LAMDA VALUE COULD BE FOUND.            
C     IF GNORM IS .LT. RCUT AND NR .EQ. TRUE, THEN TAKE SIMPEL NR-STEP.         
C     IF GNORM IS .GT. RCUT OR NR .EQ. FALSE, THEN MODIFY HESSIAN.              
C     IF HSSHFT IS ON, THEN INCREASE HESSIAN DIAGONAL ELEMENTS WITH             
C     EIGVAL(1) + 2*EIGMIN FOR MIN SEARCH, OR EIGVAL(1+MODE) + 2*EIGMIN         
C     FOR TS SEARCH. THEN CALL FORMD AGAIN.                                     
C     IF EIGINV IS ON, THEN REVERSE THE SIGN OF ALL NEGATIVE EIGENVALUES        
C     IN THE HESSIAN OTHER THAN THE TS MODE BEING FOLLOWED, AND CALL            
C     FORMD AGAIN. IF NO ACCEPTABLE LAMDA CAN BE FOUND AFTER HESSIAN            
C     HAS BEEN MODIFIED, WE ARE OUT OF LUCK...                                  
C                                                                               
      IF (FAIL) THEN                                                            
         IF (IFLEPO.EQ.7) RETURN                                                
         IF (NR.AND.(RMX.LT.RCUT)) THEN                                         
            WRITE(6,450)RCUT                                                    
  450       FORMAT(5X,'MAX GRADIENT IS LESS THAN RCUT =',F5.1,                  
     1             ' TAKING SIMPLE NR-STEP')                                    
            DO 460 I=1,NVAR                                                     
               D(I)=ZZERO                                                       
  460       CONTINUE                                                            
            DO 480 I=1,NVAR                                                     
               TEMP=-FX(I)/EIGVAL(I)                                            
               DO 470 J=1,NVAR                                                  
                  D(J)=D(J)+TEMP*U(J,I)                                         
  470          CONTINUE                                                         
  480       CONTINUE                                                            
            FAIL=.FALSE.                                                        
            GOTO 530                                                            
         ENDIF                                                                  
         IF (HSSHFT) THEN                                                       
            NHSINC=NHSINC+1                                                     
            IF (NHSINC.GT.1) THEN                                               
               WRITE(6,640)                                                     
               WRITE(6,660)                                                     
               GOTO 620                                                         
            ENDIF                                                               
            IF (IMIN.NE.0)TEMP=TWO*EIGMIN+ABS(EIGVAL(1))                        
            IF (ITS .NE.0) THEN                                                 
               IF(MODE.EQ.0)TEMP=TWO*EIGMIN+ABS(EIGVAL(2))                      
               IF(MODE.NE.0)TEMP=TWO*EIGMIN+ABS(EIGVAL(MODE+1))                 
            ENDIF                                                               
            DO 490 I=1,NVAR                                                     
               HESS(I,I)=HESS(I,I)+TEMP                                         
               EIGVAL(I)=EIGVAL(I)+TEMP                                         
  490       CONTINUE                                                            
            WRITE(6,500)TEMP                                                    
  500       FORMAT(5X,'HESSIAN DIAGONAL INCREASED BY',F8.1)                     
            FAIL=.FALSE.                                                        
            GOTO 440                                                            
         ENDIF                                                                  
         IF (EIGINV) THEN                                                       
            NNEG=0                                                              
            DO 520 I=1,NVAR                                                     
               IF ((MODE.EQ.0).AND.(I.EQ.1)) GOTO 520                           
               IF ((MODE.NE.0).AND.(I.EQ.MODE)) GOTO 520                        
               IF (EIGVAL(I).LT.ZZERO) THEN                                     
                  NNEG=NNEG+1                                                   
                  EIGVAL(I)=-EIGVAL(I)                                          
                  WRITE(6,510)EIGVAL(I)                                         
  510             FORMAT(5X,'HESSIAN EIGENVALUE REVERSED, NEW VALUE =',F        
     18.1)                                                                      
               ENDIF                                                            
  520       CONTINUE                                                            
            IF (NNEG.GT.0) THEN                                                 
               FAIL=.FALSE.                                                     
               GOTO 440                                                         
            ELSE                                                                
               WRITE(6,650)                                                     
               WRITE(6,660)                                                     
               GOTO 620                                                         
            ENDIF                                                               
         ENDIF                                                                  
      ENDIF                                                                     
C                                                                               
C  WE NOW HAVE A NEW STEP IN D                                                  
C  CHECK THAT THE STEPSIZE DOES NOT EXCEED DMAX                                 
C  IF SO, SCALE                                                                 
C                                                                               
  530 SKAL=ONE                                                                  
      DD=DOT(D,D,NVAR)                                                          
      DD=SQRT(DD)                                                               
C                                                                               
      IF(DD.GT.DMAX) THEN                                                       
         SKAL=DMAX/DD                                                           
         DO 540 I=1,NVAR                                                        
            D(I)=D(I)*SKAL                                                      
  540    CONTINUE                                                               
         DD=DMAX                                                                
         IF(IPRNT.GE.2)WRITE(6,550)SKAL                                         
  550    FORMAT(5X,'CALCULATED STEP SIZE TOO LARGE, SCALED WITH',F9.5)          
      ENDIF                                                                     
C                                                                               
      IF(IPRNT.GE.2)WRITE(6,560)DD                                              
  560 FORMAT(5X,'STEPSIZE USED IS',F9.5)                                        
      IF(IPRNT.GE.2) THEN                                                       
         WRITE(6,'('' CALCULATED STEP'')')                                      
         WRITE(6,'(3X,8F9.5)')(D(I),I=1,NVAR)                                   
      ENDIF                                                                     
C                                                                               
C  FORM NEW XPARAM AND GEO, SAVE CURRENT GRADIENTS                              
C                                                                               
      DO 570 I=1,NVAR                                                           
         XPARAM(I)=XPARAM(I)+D(I)                                               
         OLDF(I)=GRAD(I)                                                        
         K=LOC(1,I)                                                             
         L=LOC(2,I)                                                             
         GEO(L,K)=XPARAM(I)                                                     
  570 CONTINUE                                                                  
      IF(NDEP.NE.0) CALL SYMTRY                                                 
C     CHECK STEPS AND ENOUGH TIME FOR ANOTHER PASS                              
C     TIME2=SECOND()                                                            
c     TSTEP=TIME2-TIME1                                                         
c     TLEFT=TLEFT-TSTEP                                                         
c     TIME1=TIME2                                                               
c     IF (TSTEP.LT.ZZERO)TSTEP=ZZERO                                            
c     IF (TLEFT .LT. TSTEP*TWO) GOTO 620                                        
C     IN USER UNFRIENDLY ENVIROMENT, SAVE RESULTS EVERY 1 CPU HRS               
C     ITTEST=AINT((TIME2-TIME0)/TDUMP)                                          
C     IF (ITTEST.GT.NTIME) THEN                                                 
C        LDUMP=1                                                                
C        NTIME=MAX(ITTEST,(NTIME+1))                                            
C        IPOW(9)=2                                                              
C        TT0=SECOND()-TIME0                                                     
C        CALL BKRSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,-NSTEP,NSTEP,BMAT,I        
C    1POW)                                                                      
C     ELSE                                                                      
C        LDUMP=0                                                                
C     ENDIF                                                                     
C     GET GRADIENT FOR NEW GEOMETRY AND RETURN FOR ANOTHER CYCLE                
      CALL COMPFG(XPARAM, GRAD)                   
      GOTO 240                                                                  
  580 CONTINUE                                                                  
      IFLEPO=15                                                                 
C                                                                               
C  CHECK: HAS EIGVAL BEEN SET, IF NOT THEN GO BACK AND SET IT.                  
C                                                                               
      IF(DOT(EIGVAL,EIGVAL,NVAR).LT.1.D-15)GOTO 240                             
      IF (NEG.NE.NEGREQ) THEN                                                   
         IFLEPO=16                                                              
         WRITE(6,590)                                                           
         WRITE(6,600)(EIGVAL(I),I=1,5)                                          
  590    FORMAT(/,5X,'WARNING! HESSIAN DOES NOT HAVE THE REQUIRED ',            
     1'STRUCTURE',/)                                                            
  600    FORMAT(5X,'LOWEST 5 EIGENVALUES OF HESSIAN',3X,5F8.2)                  
      ENDIF                                                                     
      LAST=1                                                                    
      IF(IPRNT.GE.2)WRITE(6,610)RMX,TOL2                                        
  610 FORMAT(/,5X,'MAXIMUM GRADIENT =',F9.5,'  IS LESS THAN CUTOFF =',          
     1F9.5,//)                                                                  
      WRITE(6,'(//,70(''-''),//)')                                              
      CALL HESSOUT(HESS,0,numat)                                                        
      RETURN                                                                    
  620 CONTINUE                                                                  
C     WE RAN OUT OF TIME, OR FORMD FAIL. DUMP RESULTS                           
      IF (TLEFT .LT. TSTEP*TWO) THEN                                            
         WRITE(6,680)                                                           
      ENDIF                                                                     
      IPOW(9)=1                                                                 
C     TT0=SECOND()-TIME0                                                        
C     CALL BKRSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,-NSTEP,NSTEP,BMAT,I           
C    1POW)                                                                      
      STOP                                                                      
  630 CONTINUE                                                                  
      WRITE(6,*)'     ERROR DURING READ OF HESSIAN FROM INPUT'                  
      STOP                                                                      
  640 FORMAT(5X,'HESSIAN DIAGONAL INCREASED TWICE WITHOUT STEP',/,5X,           
     1'STILL NO ACCEPTABLE LAMDA VALUE. TERMINATING')                           
  650 FORMAT(5X,'ALL NON-TS MODES HAVE POSITIVE EIGENVALUES',/,5X,              
     1'YET FORMD CANNOT FIND ACCEPTABLE LAMDA. TERMINATING')                    
  660 FORMAT(5X,'EF OPTIMIZATION FAILED.... SORRY....')                         
  670 FORMAT(/,5X,'EXCESS NUMBER OF OPTIMIZATION CYCLES')                       
  680 FORMAT(/,5X,'NOT ENOUGH TIME FOR ANOTHER CYCLE')                          
  690 FORMAT(5X,'TS SEARCH. TAKING P-RFO STEP')                                 
  700 FORMAT(5X,'MINIMUM SEARCH. TAKING SIMPLE RFO STEP')                       
  710 FORMAT(5X,'HESSIAN DOES NOT HAVE THE DESIRED LOCAL STRUCTURE'/            
     1       5X,'TAKING P-RFO STEP')                                            
  720 FORMAT(5X,'HESSIAN DOES NOT HAVE THE DESIRED LOCAL STRUCTURE'/            
     1       5X,'TAKING SIMPLE RFO STEP')                                       
      END                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
      SUBROUTINE FORMD(EIGVAL,FX,NVAR,FAIL,ORIDE)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      DIMENSION EIGVAL(MAXPAR),FX(MAXPAR)                                       
      DOUBLE PRECISION LAMDA,LAMDA0,LAMDA1,LAMDA2                               
      LOGICAL FAIL,ORIDE                                                        
      COMMON /MESAGE/ IFLEPO,ISCF                                               
      COMMON /FMATRX/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),                      
     1U(MAXPAR,MAXPAR),DD,MODE,NSTEP,NEGREQ,IPRNT                               
C                                                                               
      DATA ZERO/0.0D0/, HALF/0.5D0/, TOLL/1.0D-8/                               
      DATA STEP/0.05D0/, TEN/1.0D+0/, ONE/1.0D+0/, BIG/1.0D+3/                  
C                                                                               
C                                                                               
C  TS SEARCH    FORMS A STEP BY P-RFO THAT TRIES TO MAXIMIZE                    
C               ALONG THE DIRECTION OF A CHOOSEN HESSIAN MODE                   
C               AND MINIMIZE ALONG ALL OTHER MODES                              
C  MIN SEARCH   FORMS A STEP BY SIMPLE RFO THAT ATTEMPTS TO                     
C               MINIMIZE ALONG ALL HESSIAN MODES                                
C                                                                               
      MAXIT=999                                                                 
      IPRT=0                                                                    
      NUMIT=0                                                                   
      IT=0                                                                      
      IF(NEGREQ.EQ.0) GO TO 30                                                  
C                                                                               
C  (A) MAXIMIZATION ALONG ONE OF THE HESSIAN MODES                              
C                                                                               
      IF(MODE.NE.0) THEN                                                        
         CALL OVERLP(NEWMOD,NVAR)                                               
C                                                                               
C  ON RETURN FROM OVERLP, NEWMOD IS THE MODE ALONG WHICH                        
C  THE ENERGY IS TO BE MAXIMIZED                                                
C                                                                               
         IF(NEWMOD.NE.MODE.AND.IPRT.GE.2) WRITE(6,210) MODE,NEWMOD              
         MODE=NEWMOD                                                            
         IF (IPRNT.GE.5) WRITE(6,220) MODE                                      
         IT=MODE                                                                
C                                                                               
C  IF THE MODE BEING FOLLOWED IS NOW THE LOWEST MODE,                           
C  THEN SWITCH OFF MODE FOLLOWING                                               
C                                                                               
         IF(MODE.EQ.1) THEN                                                     
            MODE=0                                                              
            IF(IPRNT.GE.2)WRITE(6,230)                                          
         ENDIF                                                                  
C                                                                               
      ELSE                                                                      
C                                                                               
         IF (IPRNT.GE.5) WRITE(6,240)                                           
         IT=1                                                                   
      ENDIF                                                                     
C                                                                               
      IF (IPRNT.GE.1) THEN                                                      
         WRITE(6,10)IT,EIGVAL(IT)                                               
         WRITE(6,20)(U(I,IT),I=1,NVAR)                                          
   10    FORMAT(/,5X,'TS MODE IS NUMBER',I3,' WITH EIGENVALUE',F9.1,/,          
     15X,'AND COMPONENTS',/)                                                    
   20    FORMAT(5X,8F9.4)                                                       
      ENDIF                                                                     
C                                                                               
      LAMDA0=EIGVAL(IT)+SQRT(EIGVAL(IT)*EIGVAL(IT)+4.0D0*FX(IT)*FX(IT))         
      LAMDA0=HALF*LAMDA0                                                        
      IF(IPRNT.GE.2)WRITE(6,250) LAMDA0                                         
      IF(NVAR.EQ.1) GO TO 60                                                    
C                                                                               
C  (B) MINIMIZATION ALONG ALL OTHER MODES                                       
C                                                                               
   30 CONTINUE                                                                  
      JT=1+IT                                                                   
      IF(JT.GT.2) JT=1                                                          
C                                                                               
      IF ((IPRNT.GE.5).AND.(NEGREQ.EQ.1)) WRITE(6,270)                          
      IF ((IPRNT.GE.5).AND.(NEGREQ.EQ.0)) WRITE(6,280)                          
C                                                                               
C  SOLVE ITERATIVELY FOR LAMDA                                                  
C  INITIAL GUESS FOR LAMDA IS ZERO EXCEPT NOTE THAT                             
C  LAMDA SHOULD BE LESS THAN EIGVAL(JT)                                         
C                                                                               
      LAMDA=ZERO                                                                
      IF(EIGVAL(JT).LT.ZERO) LAMDA=EIGVAL(JT)-STEP                              
C     LETS FIRST TRY A NEWTON-RAPHSON WITH ABOVE LAMDA AS GUESS.                
C     IT WILL ALMOST ALWAYS WORK AND IT IS FAST. FFF IS FUNCTION,               
C     FF1 THE GRADIENT OF THE LAMDA FUNCTION.                                   
   40 NUMIT=NUMIT+1                                                             
      ITRY=0                                                                    
      FFF = ZERO                                                                
      FF1 = ZERO                                                                
      DO 50 I=1,NVAR                                                            
         IF (I.EQ.IT) GOTO 50                                                   
         FFF = FFF + (FX(I)*FX(I))/(LAMDA-EIGVAL(I))                            
         FF1 = FF1 - (FX(I)*FX(I))/((LAMDA-EIGVAL(I))**2)                       
   50 CONTINUE                                                                  
      FFF = FFF - LAMDA                                                         
      FF1 = FF1 - ONE                                                           
      TEMP = FFF/FF1                                                            
      LAMDA = LAMDA -TEMP                                                       
      IF (IPRNT.GE.5) WRITE(6,310)LAMDA,TEMP                                    
      IF (ABS(TEMP) .LT. TOLL) GOTO 130                                         
C     CHECK MAX NUMBER OF ITERATIONS                                            
      IF (NUMIT .GT. MAXIT) GOTO 60                                             
      GOTO 40                                                                   
C     NEWTON-RAPHSON FAILED OR PRODUCED AN UNACCEPTABLE ROOT.                   
C     FIRST TRY TO DETERMINE IF IT IS LIKELY THAT THERE IS A                    
C     A ROOT IN THE INTERVAL EIGVAL(JT) TO "-INFINITY" (=-10**6).               
   60 MAXIT=9999                                                                
      ITRY=ITRY+1                                                               
      NUMIT=0                                                                   
      LAMDA1=ZERO                                                               
      IF(EIGVAL(JT).LT.ZERO) LAMDA1=EIGVAL(JT)-STEP                             
      FF1 = ZERO                                                                
      DO 70 I=1,NVAR                                                            
         IF (I.EQ.IT) GOTO 70                                                   
         FF1 = FF1 + (FX(I)*FX(I))/(LAMDA1 - EIGVAL(I))                         
   70 CONTINUE                                                                  
      FF1 = FF1 - LAMDA1                                                        
      XSTEP = ABS(LAMDA1)/BIG                                                   
      IF (XSTEP.LT.TEN)XSTEP=TEN                                                
   80 NUMIT=NUMIT+1                                                             
      IF (NUMIT.GT.MAXIT) GOTO 190                                              
      FF2 = ZERO                                                                
      LAMDA2 = LAMDA1 - NUMIT*XSTEP                                             
      DO 90 J=1,NVAR                                                            
         IF (J.EQ.IT) GOTO 90                                                   
         FF2 = FF2 + (FX(J)*FX(J))/(LAMDA2 - EIGVAL(J))                         
   90 CONTINUE                                                                  
      FF2 = FF2 - LAMDA2                                                        
      IF (FF2*FF1 .LT. ZERO) GOTO 100                                           
      GOTO 80                                                                   
C     SINCE FF1(LAMDA1) HAS OPPOSITE SIGN AS FF2(LAMDA2) THERE                  
C     MUST BE A ROOT IN THE INTERVAL. NOW HUNT IT DOWN WITH                     
C     BRUTE FORCE BISECT METHOD.                                                
  100 NUMIT=0                                                                   
      IF(IPRNT.GE.1)WRITE(6,*)'BISECT FOUND BOUNDARIES',LAMDA1,LAMDA2           
  110 NUMIT=NUMIT+1                                                             
      IF (NUMIT .GT. MAXIT) GOTO 190                                            
      TEMP = (LAMDA1 + LAMDA2)/2                                                
      FF3 = ZERO                                                                
      DO 120 I=1,NVAR                                                           
         IF (I.EQ.IT) GOTO 120                                                  
         FF3 = FF3 + (FX(I)*FX(I))/(TEMP - EIGVAL(I))                           
  120 CONTINUE                                                                  
      FF3 = FF3 - TEMP                                                          
      IF (ABS(TEMP-LAMDA2).LT.TOLL) THEN                                        
         LAMDA=TEMP                                                             
         GOTO 130                                                               
      ENDIF                                                                     
      IF (FF3*FF1 .LT. ZERO) THEN                                               
         LAMDA2=TEMP                                                            
      ELSE                                                                      
         LAMDA1=TEMP                                                            
      ENDIF                                                                     
      GOTO 110                                                                  
C                                                                               
C  AT THIS POINT WE HAVE A SOLUTION. IF ITRY = 0 THEN IT IS FROM                
C  NEWTON-RAPHSON. CHECK TO SEE IF IT IS ACCEPTABLE. IT NOT, CHECK              
C  TO SEE IF WE MISSED ANY ROOTS. IF ITRY .GT. 0 THEN BISECT HAS                
C  FOUND ANOTHER ROOT. IT SHOULD BE OK, BUT CHECK IT ANYWAY.                    
C                                                                               
  130 IF(IPRNT.GE.2)WRITE(6,260) LAMDA                                          
      IF((IPRNT.GE.1).AND.(ITRY.GT.0))WRITE(6,*)'LAMDA FOUND BY BISECT'         
      IF ((LAMDA.GT.EIGVAL(JT)) .OR.                                            
     1   (EIGVAL(JT).GT.ZERO.AND.LAMDA.GT.ZERO)) THEN                           
         IF (ITRY.GT.0) THEN                                                    
            GOTO 180                                                            
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
      ENDIF                                                                     
C                                                                               
C  CALCULATE THE STEP                                                           
C                                                                               
  140 DO 150 I=1,NVAR                                                           
         D(I)=ZERO                                                              
  150 CONTINUE                                                                  
      DO 170 I=1,NVAR                                                           
         TEMP=FX(I)/(LAMDA-EIGVAL(I))                                           
         IF(I.EQ.IT) TEMP=FX(I)/(LAMDA0-EIGVAL(I))                              
         DO 160 J=1,NVAR                                                        
            D(J)=D(J)+TEMP*U(J,I)                                               
  160    CONTINUE                                                               
  170 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C                                                                               
C                                                                               
C    ERROR SECTION                                                              
C                                                                               
  180 CONTINUE                                                                  
C     IF WE GET TO HERE BISECT HAS FOUND AN UNACCEPTABLE ROOT.                  
C     THIS SHOULD NOT COULD HAPPEN..... SET IFLEPO=7                            
      WRITE(6,290)                                                              
      IF((IPRNT.GE.5).AND.(ITRY.GT.0))WRITE(6,*)'BISECT METHOD FAILED,          
     1NUMIT,FF1,FF2,FF3,LAMD1,LAMD2',NUMIT,FF1,FF2,FF3,LAMDA1,LAMDA2            
      FAIL=.TRUE.                                                               
      IFLEPO=7                                                                  
      RETURN                                                                    
  190 CONTINUE                                                                  
      WRITE(6,300)                                                              
      IF((IPRNT.GE.5).AND.(ITRY.GT.0))WRITE(6,*)'BISECT METHOD FAILED,          
     1NUMIT,FF1,FF2,FF3,LAMD1,LAMD2',NUMIT,FF1,FF2,FF3,LAMDA1,LAMDA2            
C     IF OVERRIDE IS TURNED ON, USE LAMDA ANYWAY.....                           
      IF (ORIDE) WRITE(*,200)                                                   
  200 FORMAT(5X,'** OVERIDE IS TURNED ON, USING LAMDA(S) ANYWAY **')            
      IF (ORIDE) GOTO 140                                                       
      FAIL=.TRUE.                                                               
      RETURN                                                                    
C                                                                               
  210 FORMAT(5X,'WARNING!! MODE SWITCHING. WAS FOLLOWING MODE ',I3,             
     1       ' NOW FOLLOWING MODE ',I3)                                         
  220 FORMAT(5X,'SEARCHING FOR LAMDA THAT MAXIMIZES ALONG MODE ',I3)            
  230 FORMAT(5X,'MODE FOLLOWING SWITCHED OFF')                                  
  240 FORMAT(' SEARCHING FOR LAMDA THAT MAXIMIZES ALONG THE',                   
     1       ' LOWEST MODE')                                                    
  250 FORMAT(5X,'LAMDA THAT MAXIMIZES ALONG TS MODE =           ',F13.7)        
  260 FORMAT(5X,'LAMDA THAT MINIMIZES ALONG ALL (OTHER) MODES = ',F13.7)        
  270 FORMAT(' SEARCHING FOR LAMDA THAT MINIMIZES ALONG ALL',                   
     1       ' OTHER MODES')                                                    
  280 FORMAT(' SEARCHING FOR LAMDA THAT MINIMIZES ALONG ALL MODES')             
  290 FORMAT(//' *****************************************'/                    
     1         ' ** ERROR IN DETERMINING LAMDA IN FORMD **'/                    
     2         ' *****************************************'//)                  
  300 FORMAT(  ' *** UNABLE TO DETERMINE LAMDA IN FORMD ***')                   
C070  FORMAT( /' *****************************************'/                    
C    $         ' *** UNABLE TO DETERMINE LAMDA IN FORMD **'/                    
C    $         ' *****************************************'/ )                  
  310 FORMAT(' IN ITERATIVE CYCLE:  LAMDA= ',F20.5,' TEMP= ',F25.10)            
      END                                                                       
      SUBROUTINE OVERLP(NEWMOD,NVAR)                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      COMMON /FMATRX/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),                      
     1U(MAXPAR,MAXPAR),DD,MODE,NSTEP,NEGREQ,IPRNT                               
C                                                                               
C  ON THE FIRST STEP SIMPLY DETERMINE WHICH MODE TO FOLLOW                      
C                                                                               
      IF(NSTEP.EQ.1) THEN                                                       
C                                                                               
         IF(MODE.GT.NVAR)THEN                                                   
            WRITE(6,*)'ERROR!! MODE IS LARGER THAN NVAR',MODE                   
            STOP                                                                
         ENDIF                                                                  
C                                                                               
         IT=MODE                                                                
         WRITE(6,40) MODE                                                       
C                                                                               
      ELSE                                                                      
C                                                                               
C  ON SUBSEQUENT STEPS DETERMINE WHICH HESSIAN EIGENVECTOR HAS                  
C  THE GREATEST OVERLAP WITH THE MODE WE ARE FOLLOWING                          
C                                                                               
         IT=1                                                                   
         TOVLP=DOT(U(1,1),VMODE,NVAR)                                           
         TOVLP=ABS(TOVLP)                                                       
         DO 10 I=2,NVAR                                                         
            OVLP=DOT(U(1,I),VMODE,NVAR)                                         
            OVLP=ABS(OVLP)                                                      
            IF(OVLP.GT.TOVLP) THEN                                              
               TOVLP=OVLP                                                       
               IT=I                                                             
            ENDIF                                                               
   10    CONTINUE                                                               
C                                                                               
         WRITE(6,30) TOVLP                                                      
      ENDIF                                                                     
C                                                                               
C  SAVE THE EIGENVECTOR IN VMODE                                                
C                                                                               
      DO 20 I=1,NVAR                                                            
         VMODE(I)=U(I,IT)                                                       
   20 CONTINUE                                                                  
C                                                                               
      NEWMOD=IT                                                                 
      RETURN                                                                    
C                                                                               
   30 FORMAT(5X,'OVERLAP OF CURRENT MODE WITH PREVIOUS MODE IS ',F12.6)         
C                                                                               
   40 FORMAT(5X,'HESSIAN MODE FOLLOWING SWITCHED ON'/                           
     1     '     FOLLOWING MODE ',I3)                                           
C                                                                               
      END                                                                       
      SUBROUTINE UPDHES(SVEC,TVEC,NVAR,IUPD)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      DIMENSION TVEC(*),SVEC(*)                                                 
      LOGICAL FIRST                                                             
      COMMON /FMATRX/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),                      
     1U(MAXPAR,MAXPAR),DD,MODE,NSTEP,NEGREQ,IPRNT                               
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR*3)                                     
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA                                       
C                                                                               
      DATA ZERO/0.0D0/                                                          
C                                                                               
C  UPDATING OF THE HESSIAN                                                      
C  DEPENDS ON CURRENT GRADIENTS, OLD GRADIENTS AND THE                          
C  CORRECTION VECTOR USED ON THE LAST CYCLE                                     
C  SVEC & TVEC ARE FOR TEMPORARY STORAGE                                        
C                                                                               
C  2 UPDATING PROCEDURES ARE POSSIBLE                                           
C  (I)   THE POWELL UPDATE                                                      
C        THIS PRESERVES THE SYMMETRIC CHARACTER OF THE HESSIAN                  
C        WHILST ALLOWING ITS EIGENVALUE STRUCTURE TO CHANGE.                    
C        IT IS THE DEFAULT UPDATE FOR A TRANSITION STATE SEARCH                 
C  (II)  THE BFGS UPDATE                                                        
C        THIS UPDATE HAS THE IMPORTANT CHARACTERISTIC OF RETAINING              
C        POSITIVE DEFINITENESS (NOTE: THIS IS NOT RIGOROUSLY                    
C        GUARANTEED, BUT CAN BE CHECKED FOR BY THE PROGRAM).                    
C        IT IS THE DEFAULT UPDATE FOR A MINIMUM SEARCH                          
C                                                                               
C     SWITCH : IUPD                                                             
C       IUPD = 0  :  SKIP UPDATE                                                
C       IUPD = 1  :  POWELL                                                     
C       IUPD = 2  :  BFGS                                                       
C                                                                               
      IF (.NOT. FIRST) THEN                                                     
         FIRST=.TRUE.                                                           
         IF(IPRNT.GE.2) THEN                                                    
            IF (IUPD.EQ.0)WRITE(6,90)                                           
            IF (IUPD.EQ.1)WRITE(6,80)                                           
            IF (IUPD.EQ.2)WRITE(6,120)                                          
         ENDIF                                                                  
      ENDIF                                                                     
      IF(IUPD.EQ.0) RETURN                                                      
      DO 10 I=1,NVAR                                                            
         TVEC(I)=ZERO                                                           
         DO 10 J=1,NVAR                                                         
            TVEC(I)=TVEC(I) + HESS(I,J)*D(J)                                    
   10 CONTINUE                                                                  
C                                                                               
      IF(IUPD.EQ.1) THEN                                                        
C                                                                               
C   (I) POWELL UPDATE                                                           
C                                                                               
         DO 20 I=1,NVAR                                                         
            TVEC(I)=GRAD(I)-OLDF(I)-TVEC(I)                                     
   20    CONTINUE                                                               
         DDS=DD*DD                                                              
         DDTD=DOT(TVEC,D,NVAR)                                                  
         DDTD=DDTD/DDS                                                          
C                                                                               
         DO 40 I=1,NVAR                                                         
            DO 30 J=1,I                                                         
               TEMP=TVEC(I)*D(J) + D(I)*TVEC(J) - D(I)*DDTD*D(J)                
               HESS(I,J)=HESS(I,J)+TEMP/DDS                                     
               HESS(J,I)=HESS(I,J)                                              
   30       CONTINUE                                                            
   40    CONTINUE                                                               
C                                                                               
      ENDIF                                                                     
      IF (IUPD.EQ.2) THEN                                                       
C                                                                               
C  (II) BFGS UPDATE                                                             
C                                                                               
         DO 50 I=1,NVAR                                                         
            SVEC(I)=GRAD(I)-OLDF(I)                                             
   50    CONTINUE                                                               
         DDS=DOT(SVEC,D,NVAR)                                                   
C                                                                               
C  IF DDS IS NEGATIVE, RETENTION OF POSITIVE DEFINITENESS IS NOT                
C  GUARANTEED. PRINT A WARNING AND SKIP UPDATE THIS CYCLE.                      
C                                                                               
         IF(DDS.LT.ZERO) THEN                                                   
            WRITE(6,100)                                                        
            WRITE(6,110)                                                        
            RETURN                                                              
         ENDIF                                                                  
C                                                                               
         DDTD=DOT(D,TVEC,NVAR)                                                  
C                                                                               
         DO 70 I=1,NVAR                                                         
            DO 60 J=1,I                                                         
               TEMP= (SVEC(I)*SVEC(J))/DDS - (TVEC(I)*TVEC(J))/DDTD             
               HESS(I,J)=HESS(I,J)+TEMP                                         
               HESS(J,I)=HESS(I,J)                                              
   60       CONTINUE                                                            
   70    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
C                                                                               
   80 FORMAT(/,5X,'HESSIAN IS BEING UPDATED USING THE POWELL UPDATE',/)         
   90 FORMAT(/,5X,'HESSIAN IS NOT BEING UPDATED',/)                             
  100 FORMAT(5X,'WARNING! HEREDITARY POSITIVE DEFINITENESS ENDANGERED')         
  110 FORMAT(5X,'UPDATE SKIPPED THIS CYCLE')                                    
  120 FORMAT(/,5X,'HESSIAN IS BEING UPDATED USING THE BFGS UPDATE',/)           
      END                                                                       
