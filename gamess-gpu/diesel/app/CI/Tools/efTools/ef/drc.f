      SUBROUTINE DRC(MODE,COORD)                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
                                                                                
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /MASS  / AMSS(107)                                                 
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)                          
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),          
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,                
     2                NCLOSE,NOPEN,NDUMY,XRACT                                  
      CHARACTER KEYWRD*241                                                      
                                                                                
      DIMENSION  VELO0(MAXPAR), VELO1(MAXPAR), STARTV(MAXPAR),                  
     1VELO2(MAXPAR), VELO3(MAXPAR), GERROR(MAXPAR), XPARAM(MAXPAR),             
     2COORD(3,NUMATM), GROLD2(MAXPAR), PAST10(10),                              
     3GROLD(MAXPAR), GXYZ(3,NUMATM), GRAD(MAXPAR), ATMASS(MAXPAR)               
                                                                                
      LOGICAL  ADDK, LETOT, LET, VELRED                                         
                                                                                
      DATA VELO0 /MAXPAR*0.D0/                                                  
      DATA STARTV/MAXPAR*0.D0/                                                  
      DATA ADDK/.TRUE./                                                         
                                                                                
      DELOLD=10.D0                                                              
      GTOT=0.D0                                                                 
      ACCU=0.25D0                                                               
      GNLIM=0.5D0                                                               
      PAST10(5)=100.D0                                                          
      VELRED=(INDEX(KEYWRD,'VELO').NE.0)                                        
      MAXLOOP=10                                                                
      I=INDEX(KEYWRD,'MAXLOOP=')                                                
      IF(I.NE.0) MAXLOOP=IDINT(READA(KEYWRD,I))                                 
      IF(VELRED) THEN                                                           
         OPEN(UNIT=10,FILE='startv')                                            
         READ(10,'(3F13.5)')(STARTV(I),I=1,NUMAT*3)                             
         CLOSE (10)                                                             
      ENDIF                                                                     
      IF(DOT(STARTV,STARTV,3*NUMAT).GT.0.001D0)THEN                             
C                                                                               
C     PRINT OUT INITIAL VELOCITIES                                              
C                                                                               
         WRITE(6,'(A)')' INITIAL VELOCITY IN DRC'                               
         WRITE(6,'(3F13.5)')(STARTV(I),I=1,NUMAT*3)                             
      ENDIF                                                                     
      LET=(INDEX(KEYWRD,' GEO-OK').NE.0.OR.VELRED)                              
C                                                                               
C  TRANSFER COORDINATES TO XPARAM AND LOC                                       
C                                                                               
      L=0                                                                       
      DO 50 I=1,NATOMS                                                          
         IF(LABELS(I).GT.98) GOTO 50                                            
         XPARAM(L+1)=COORD(1,I)                                                 
         XPARAM(L+2)=COORD(2,I)                                                 
         XPARAM(L+3)=COORD(3,I)                                                 
         ATMASS(L+1)=AMSS(LABELS(I))                                            
         ATMASS(L+2)=AMSS(LABELS(I))                                            
         ATMASS(L+3)=AMSS(LABELS(I))                                            
c        write(*,'(i3,4F8.4)') k,atmass(k),                                     
c    .xparam(l+1),xparam(l+2),xparam(l+3)                                       
         L=L+3                                                                  
   50 CONTINUE                                                                  
      NVAR=NUMAT*3                                                              
C                                                                               
C DETERMINE DAMPING FACTOR                                                      
C                                                                               
      HALF=0.0D0                                                                
      IF(INDEX(KEYWRD,'DRC=').NE.0) THEN                                        
         HALF=READA(KEYWRD,INDEX(KEYWRD,'DRC='))                                
         WRITE(6,'(//10X,'' DAMPING FACTOR FOR KINETIC ENERGY ='',F12.6)        
     1')HALF                                                                    
      ELSEIF (INDEX(KEYWRD,'DRC').EQ.0) THEN                                    
         HALF=0.D0                                                              
      ELSE                                                                      
         HALF=1.D6                                                              
      ENDIF                                                                     
C                                                                               
C  LETOT IS TRUE IF CORRECTIONS ARE NOT TO BE MADE PART WAY INTO                
C        THE CALCULATION                                                        
C                                                                               
C  USAGE OF LETOT:                                                              
C (1) WHILE LETOT IS FALSE, NO DAMPING WILL BE DONE                             
C (2) WHEN LETOT IS TURNED TRUE,                                                
C     IF AN IRC, THEN ETOT IS RESET SO THE ERROR IS ZERO.                       
C     IF A  DRC, EXCESS KINETIC ENERGY USED TO START THE RUN IS REMOVED.        
C                                                                               
      LETOT=(INDEX(KEYWRD,'IRC=').EQ.0 .AND. .NOT. LET)                         
      HALF=SIGN(MAX(0.000001D0,ABS(HALF)),HALF)                                 
C                                                                               
C DETERMINE EXCESS KINETIC ENERGY                                               
C                                                                               
      ISKIN=0                                                                   
      ADDONK=0.D0                                                               
C                                                                               
C   LOOP OVER TIME-INTERVALS OF DELTAT SECOND                                   
C                                                                               
      DELTAT=1.D-16                                                             
      QUADR=1.D0                                                                
      ETOT=0.D0                                                                 
      ESCF=0.D0                                                                 
      CONST=1.D0                                                                
C                         NOT A RESTART                                         
         ILOOP=1                                                                
         IF(INDEX(KEYWRD,'IRC=').NE.0.OR.VELRED)THEN                            
C                                                                               
C  GET HOLD OF VELOCITY VECTOR                                                  
C                                                                               
            IF(INDEX(KEYWRD,'IRC=').NE.0)THEN                                   
               K=READA(KEYWRD,INDEX(KEYWRD,'IRC='))                             
            ELSE                                                                
               K=1                                                              
            ENDIF                                                               
            IF(K.LT.0)THEN                                                      
               K=-K                                                             
               ONE=-1.D0                                                        
            ELSE                                                                
               ONE=1.D0                                                         
            ENDIF                                                               
            KL=(K-1)*NVAR                                                       
            SUMM=0.D0                                                           
            VELO1(1)=0                                                          
            VELO1(2)=0                                                          
            VELO1(3)=0                                                          
            SUMMAS=0.D0                                                         
            I=0                                                                 
            DO 60 II=1,NUMAT                                                    
               AMS=ATMASS(II)                                                   
               SUMMAS=SUMMAS+AMS                                                
               VELO0(I+1)=STARTV(KL+I+1)*ONE                                    
               VELO1(1)=VELO1(1)+VELO0(I+1)*AMS                                 
C                                                                               
               VELO0(I+2)=STARTV(KL+I+2)*ONE                                    
               VELO1(2)=VELO1(2)+VELO0(I+2)*AMS                                 
C                                                                               
               VELO0(I+3)=STARTV(KL+I+3)*ONE                                    
               VELO1(3)=VELO1(3)+VELO0(I+3)*AMS                                 
C                                                                               
               I=I+3                                                            
   60       CONTINUE                                                            
C$DOIT ASIS                                                                     
            DO 70 I=1,3                                                         
   70       VELO1(I)=-VELO1(I)/SUMMAS                                           
            I=0                                                                 
C$DOIT ASIS                                                                     
            DO 80 II=1,NUMAT                                                    
               AMS=ATMASS(II)                                                   
C$DOIT ASIS                                                                     
               DO 80 I1=1,3                                                     
                  I=I+1                                                         
                  IF(ADDONK.GT.1.D-5.OR..NOT.VELRED)VELO0(I)=VELO0(I)+VE        
     1LO1(I1)                                                                   
   80       SUMM=SUMM+VELO0(I)**2*AMS                                           
            IF(ADDONK.LT.1.D-5.AND.VELRED)ADDONK=0.5D0*SUMM/4.184D10            
            ENDIF                                                               
C                                                                               
C   AT THIS POINT ADDONK IS IN KCAL/MOLE                                        
C   NORMALIZE SO THAT TOTAL K.E. = ONE QUANTUM (DEFAULT) (DRC ONLY)             
C                              OR 0.3KCAL/MOLE (IRC ONLY)                       
C                              OR ADDONK IF KINETIC=NN SUPPLIED                 
C                                                                               
C           IF(SUMM.LT.1.D-4) THEN                                              
C              WRITE(6,'(A)')' SYSTEM IS APPARENTLY NOT MOVING!'                
C              RETURN                                                           
C           ENDIF                                                               
C                                                                               
C  ADDONK IS EXCESS KINETIC ENERGY.  IF THE CALCULATION IS AN IRC,              
C  THIS ENERGY MUST BE REMOVED AFTER A SHORT 'TIME'.                            
C                                                                               
C  MAKE AN AD-HOC CORRECTION: IF ADDONK IS NON-ZERO AND HALF IS LARGER          
C  THAN 0.1, MODIFY ADDONK TO REFLECT ERRORS DUE TO START-UP.                   
C                                                                               
            IF(HALF.GT.0.1D0.AND.HALF.LT.10000.D0)                              
     1ADDONK=ADDONK*(1.D0+0.06972D0/HALF)                                       
C                                                                               
C  MAKE AN AD-HOC CORRECTION: IF ADDONK IS NON-ZERO AND HALF IS LESS            
C  THAN -0.1, MODIFY ADDONK TO REFLECT ERRORS DUE TO START-UP.                  
C                                                                               
            IF(HALF.LT.-0.1D0.AND.HALF.GT.-10000.D0)                            
     1ADDONK=ADDONK*(1.D0+0.06886D0/HALF)                                       
            SUMM=SQRT(ADDONK/(0.5D0*SUMM/4.184D10))                             
            ADDK=.FALSE.                                                        
            IF(SUMM.GT.1.D-10)THEN                                              
               DO  90 I=1,NVAR                                                  
   90          VELO0(I)=VELO0(I)*SUMM                                           
C                                                                               
C  IF IT IS A DRC, DESTROY ADDONK.  THE KINETIC ENERGY USED WILL COME           
C  FROM THE VELOCITY ONLY.                                                      
C                                                                               
               IF(HALF.GT.1.D-3)ADDONK=0.D0                                     
            ENDIF                                                               
  100 CONTINUE                                                                  
                                                                                
                                                                                
      WRITE(6,101) ADDONK/627.51                                                
  101 FORMAT('INITIAL KINETIC ENERGY : ',F12.8,' AU')                           
      WRITE(6,102)                                                              
  102 FORMAT(//,6x,'STEP',5X,'DELTA T',5X,'EPOT',5X,                            
     .             'EKIN',5X,'ETOT',5X,'ERR (%)',//)                            
                                                                                
         K=0                                                                    
         DO I=1,NATOMS                                                          
            IF(LABELS(I).LT.99) THEN                                            
               COORD(1,I)=XPARAM(K+1)                                           
               COORD(2,I)=XPARAM(K+2)                                           
               COORD(3,I)=XPARAM(K+3)                                           
               K=K+3                                                            
            ENDIF                                                               
         ENDDO                                                                  
         CALL OUT(1,MODE,COORD)                                                 
         CALL EGRAD(MODE,NUMAT,ESCF,GXYZ)                                       
         IF(DABS(ESCF).LT.1.0D-8) STOP 'EXTERNAL PROGRAM ERROR'                 
         K=0                                                                    
         DO I=1,NUMAT                                                           
            GRAD(K+1)=GXYZ(1,I)                                                 
            GRAD(K+2)=GXYZ(2,I)                                                 
            GRAD(K+3)=GXYZ(3,I)                                                 
            K=K+3                                                               
         ENDDO                                                                  
         ESCF=ESCF*627.51                                                       
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                      BEGIN OF TIME LOOP                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
                                                                                
      ONE=1.D0                                                                  
      EKIN=0.0D0                                                                
      ETOT=ESCF                                                                 
      DO 190 ILOOP=1,MAXLOOP                                                    
C                                                                               
C  MOVEMENT OF ATOMS WILL BE PROPORTIONAL TO THE AVERAGE VELOCITIES             
C  OF THE ATOMS BEFORE AND AFTER TIME INTERVAL                                  
C                                                                               
C                                                                               
C  RAPID CHANGE IN GRADIENT IMPLIES SMALL STEP SIZE FOR DELTAT                  
C                                                                               
C   KINETIC ENERGY = 1/2 * M * V * V                                            
C                  = 0.5 / (4.184D10) * M * V * V                               
C   NEW VELOCITY = OLD VELOCITY + GRADIENT * TIME / MASS                        
C                = KCAL/ANGSTROM*SECOND/(ATOMIC WEIGHT)                         
C                =4.184*10**10(ERGS)*10**8(PER CM)*DELTAT(SECONDS)              
C   NEW POSITION = OLD POSITION - AVERAGE VELOCITY * TIME INTERVAL              
C                                                                               
C                                                                               
C   ESTABLISH REFERENCE TOTAL ENERGY                                            
C                                                                               
         ERROR=(ETOT-(EKIN+ESCF))                                               
         IF(ILOOP.GT.2)THEN                                                     
            QUADR = 1.D0+ERROR/(EKIN*CONST+0.001D0)*0.5D0                       
            QUADR = MIN(1.3D0,MAX(0.8D0,QUADR))                                 
         ELSE                                                                   
            QUADR=1.D0                                                          
         ENDIF                                                                  
         IF((LET.OR.EKIN.GT.0.2).AND.ADDK)THEN                                  
C                                                                               
C   DUMP IN EXCESS KINETIC ENERGY                                               
C                                                                               
            ETOT=ETOT+ADDONK                                                    
            ADDK=.FALSE.                                                        
            ADDONK=0.D0                                                         
         ENDIF                                                                  
C                                                                               
C  CALCULATE THE DURATION OF THE NEXT STEP.                                     
C  STEP SIZE IS THAT REQUIRED TO PRODUCE A CONSTANT CHANGE IN GEOMETRY          
C                                                                               
C                                                                               
C  IF DAMPING IS USED, CALCULATE THE NEW TOTAL ENERGY AND                       
C  THE RATIO FOR REDUCING THE KINETIC ENERGY                                    
C                                                                               
         CONST=MAX(1.D-36,0.5D0**(DELTAT*1.D15/HALF))                           
         CONST=SQRT(CONST)                                                      
         VELVEC=0.D0                                                            
         EKIN=0.D0                                                              
         DELTA1=DELOLD+DLOLD2                                                   
         ELOST=0.D0                                                             
         DO 110 I=1,NVAR                                                        
C                                                                               
C   CALCULATE COMPONENTS OF VELOCITY AS                                         
C   V = V(0) + V'*T + V"*T*T                                                    
C   WE NEED ALL THREE TERMS, V(0), V' AND V"                                    
C                                                                               
            VELO1(I) = 1.D0/ATMASS(I)*GRAD(I)                                   
            IF(ILOOP.GT.3) THEN                                                 
               VELO3(I) = 2.D0/ATMASS(I)*                                       
     1(DELTA1*(GROLD(I)-GRAD(I))-DELOLD*(GROLD2(I)-GRAD(I)))/                   
     2(DELTA1*(DELOLD**2*1.D30)-DELOLD*(DELTA1**2*1.D30))                       
               VELO2(I)=1.D0/ATMASS(I)*                                         
     1(GRAD(I)-GROLD(I)-0.5D0*VELO3(I)*(1.D30*DELOLD**2))/(DELOLD*1.D15)        
            ELSE                                                                
               VELO2(I) = 1.D0/ATMASS(I)*                                       
     1                 (GRAD(I)-GROLD(I))/(1.D15*DELOLD)                        
               VELO3(I)=0.D0                                                    
            ENDIF                                                               
C                                                                               
C  MOVE ATOMS THROUGH DISTANCE EQUAL TO VELOCITY * DELTA-TIME, NOTE             
C  VELOCITY CHANGES FROM START TO FINISH, THEREFORE AVERAGE.                    
C                                                                               
            XPARAM(I)=XPARAM(I)                                                 
     1             -1.D8*(DELTAT*VELO0(I)*ONE                                   
     2             +0.5D0*DELTAT**2*VELO1(I)                                    
     3             +0.16666D0*(DELTAT**2*1.D15)*DELTAT*VELO2(I)                 
     4             +0.0416666D0*DELTAT**2*(1.D30*DELTAT**2)*VELO3(I))           
C                                                                               
C   CORRECT ERRORS DUE TO CUBIC COMPONENTS IN ENERGY GRADIENT,                  
C   ALSO TO ADD ON EXCESS ENERGY, IF NECESSARY.                                 
C                                                                               
            VELVEC=VELVEC+VELO0(I)**2                                           
C                                                                               
C   MODIFY VELOCITY IN LIGHT OF CURRENT ENERGY GRADIENTS.                       
C                                                                               
C   VELOCITY = OLD VELOCITY + (DELTA-T / ATOMIC MASS) * CURRENT GRADIENT        
C                           + 1/2 *(DELTA-T * DELTA-T /ATOMIC MASS) *           
C                             (SLOPE OF GRADIENT)                               
C              SLOPE OF GRADIENT = (GRAD(I)-GROLD(I))/DELOLD                    
C                                                                               
C                                                                               
C   THIS EXPRESSION IS ACCURATE TO SECOND ORDER IN TIME.                        
C                                                                               
            VELO0(I) = VELO0(I) + DELTAT*VELO1(I) + 0.5D0*DELTAT**2*VELO        
     12(I)*1.D15           + 0.166666D0*DELTAT*(1.D30*DELTAT**2)*VELO3(         
     2I)                                                                        
            IF(LET.OR.GNORM.GT.3.D0)THEN                                        
               LET=.TRUE.                                                       
               ELOST=ELOST+VELO0(I)**2*ATMASS(I)*(1-CONST**2)                   
               VELO0(I)=VELO0(I)*CONST*QUADR                                    
            ENDIF                                                               
C                                                                               
C  CALCULATE KINETIC ENERGY (IN 2*ERGS AT THIS POINT)                           
C                                                                               
            EKIN=EKIN+VELO0(I)**2*ATMASS(I)                                     
  110    CONTINUE                                                               
         ONE=1.D0                                                               
         IF(LET.OR.GNORM.GT.3.D0)THEN                                           
            IF(.NOT.LETOT) THEN                                                 
               IF(ABS(HALF).LT.1.D-3)THEN                                       
C                                                                               
C  IT IS AN IRC, SO RESET THE TOTAL ENERGY                                      
C                                                                               
                  ETOT=ESCF+ELOST1                                              
                  ADDONK=0.D0                                                   
                  ELOST1=0.D0                                                   
                  ELOST=0.D0                                                    
               ELSEIF(ISKIN.EQ.0)THEN                                           
C                                                                               
C  IT IS A DRC AND KINETIC NOT USED, SO REMOVE EXTRA KINETIC ENERGY             
C                                                                               
                  ETOT=ETOT-ADDONK                                              
               ENDIF                                                            
            ENDIF                                                               
            LETOT=.TRUE.                                                        
         ENDIF                                                                  
C                                                                               
C  CONVERT ENERGY INTO KCAL/MOLE                                                
C                                                                               
         EKIN=0.5*EKIN/4.184D10                                                 
C                                                                               
C  IF IT IS A DAMPED DRC, MODIFY ETOT TO REFLECT LOSS OF KINETIC ENERGY         
C                                                                               
         IF(LETOT.AND.ABS(HALF).GT.0.00001D0)                                   
     1ETOT=ETOT-EKIN/CONST**2+EKIN                                              
         ELOST1=ELOST1+0.5D0*ELOST/4.184D10                                     
C                                                                               
C STORE OLD GRADIENTS FOR DELTA - VELOCITY CALCULATION                          
C                                                                               
         DO 120 I=1,NVAR                                                        
            GROLD2(I)=GROLD(I)                                                  
            GROLD(I)=GRAD(I)                                                    
  120    GRAD(I)=0.D0                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
C   CALCULATE ENERGY AND GRADIENTS                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
         SCFOLD=ESCF                                                            
         K=0                                                                    
         DO I=1,NATOMS                                                          
            IF(LABELS(I).LT.99) THEN                                            
               COORD(1,I)=XPARAM(K+1)                                           
               COORD(2,I)=XPARAM(K+2)                                           
               COORD(3,I)=XPARAM(K+3)                                           
               K=K+3                                                            
            ENDIF                                                               
         ENDDO                                                                  
c        do i=1,natoms                                                          
c        write(*,'(i3,4F8.4)') i,coord(1,i),coord(2,i),coord(3,i)               
c        enddo                                                                  
c        stop                                                                   
         CALL OUT(1,MODE,COORD)                                                 
         CALL EGRAD(MODE,NUMAT,ESCF,GXYZ)                                       
         IF(DABS(ESCF).LT.1.0D-8) STOP 'EXTERNAL PROGRAM ERROR'                 
         K=0                                                                    
         DO I=1,NUMAT                                                           
            GRAD(K+1)=GXYZ(1,I)                                                 
            GRAD(K+2)=GXYZ(2,I)                                                 
            GRAD(K+3)=GXYZ(3,I)                                                 
            K=K+3                                                               
         ENDDO                                                                  
         ESCF=ESCF*627.510D0                                                    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
         IF(ILOOP.GT.2)THEN                                                     
            GNORM=0.D0                                                          
            DO 140 I=1,NVAR,3                                                   
               SUM=SQRT(DOT(GRAD(I),GRAD(I),3)/                                 
     1(DOT(VELO0(I),VELO0(I),3)+1.D-20))                                        
               DO 130 J=I,I+2                                                   
  130          GERROR(J)=GERROR(J)+GRAD(J)+VELO0(J)*SUM                         
  140       CONTINUE                                                            
            GNORM=SQRT(DOT(GERROR,GERROR,NVAR))                                 
            GTOT=GNORM                                                          
         ENDIF                                                                  
         GNORM=SQRT(DOT(GRAD,GRAD,NVAR))                                        
C                                                                               
C   CONVERT GRADIENTS INTO ERGS/CM                                              
C                                                                               
         DO 150 I=1,NVAR                                                        
  150    GRAD(I)=GRAD(I)*4.184D18                                               
C                                                                               
C   SPECIAL TREATMENT FOR FIRST POINT - SET "OLD" GRADIENTS EQUAL TO            
C   CURRENT GRADIENTS.                                                          
C                                                                               
         IF(ILOOP.EQ.1) THEN                                                    
            DO 160 I=1,NVAR                                                     
  160       GROLD(I)=GRAD(I)                                                    
         ENDIF                                                                  
         DLOLD2=DELOLD                                                          
         DELOLD=DELTAT                                                          
         SUM=0.D0                                                               
         DO 170 I=1,NVAR                                                        
  170    SUM=SUM + ((GRAD(I)-GROLD(I))/4.184D18)**2                             
         IF(ABS(HALF).LT.0.001D0)THEN                                           
            DELTAT= DELTAT*                                                     
     1MIN(2.D0, (5.D-5*ACCU/(ABS(ESCF+ELOST1-ETOLD)+1.D-20)))**0.25D0           
            ETOLD=ESCF+ELOST1                                                   
            IF(ILOOP.GT.5.AND.SCFOLD-ESCF.LT.-1.D-3 .OR.                        
     1      ILOOP.GT.30.AND.SCFOLD-ESCF.LT.0.D0)  THEN                          
               WRITE(6,'(//,'' IRC CALCULATION COMPLETE '')')                   
               RETURN                                                           
            ENDIF                                                               
         ELSE                                                                   
            DELTAT= DELTAT*MIN(1.05D0, 10.D0*ACCU/(SUM+1.D-4))                  
            DELTAT=MIN(DELTAT,3.D-15*ACCU)                                      
            PAST10(10)=GNORM                                                    
            SUM=0.D0                                                            
            DO 180 I=1,9                                                        
               SUM=SUM+ABS(PAST10(I)-PAST10(I+1))                               
  180       PAST10(I)=PAST10(I+1)                                               
            IF(SUM.LT.GNLIM)THEN                                                
               WRITE(6,'(//,A)')' GRADIENT CONSTANT AND SMALL -- ASSUME'        
     1//' ALL MOTION STOPPED'                                                   
               RETURN                                                           
            ENDIF                                                               
            DELTAT=MIN(DELTAT,2.D-15)                                           
         ENDIF                                                                  
         DELTAT=MAX(1.D-16,DELTAT)                                              
C                                                                               
C ESCF   = POTENTIAL ENERGY                                                     
C EKIN   = CURRENT KINETIC ENERGY                                               
C ETOT   = COMPUTED TOTAL ENERGY = STARTING POTENTIAL ENERGY -                  
C          KINETIC ENERGY LOST THROUGH DAMPING, IF PRESENT.                     
C                                                                               
C   IN DRCOUT  'TOTAL' = ESCF + EKIN                                            
C              'ERROR' = ESCF + EKIN - ETOT                                     
C                                                                               
                                                                                
            WRITE(6,300) ILOOP,DELTAT*1.0d15,ESCF/627.51,EKIN/627.51,           
     .                   ETOT/627.51,100.*(ESCF+EKIN-ETOT)/(ETOT+1.0D-5)        
  300       FORMAT(3X,I4,F6.2,5X,3F12.8,5X,F6.2)                                
                                                                                
            IF(MOD(ILOOP,5).EQ.0) THEN                                          
               WRITE(6,'(//,'' ATOMIC VELOCITIES (M/S)'',//)')                  
               WRITE(6,'(3F10.4)')(0.01*VELO0(I),I=1,NUMAT*3)                   
            ENDIF                                                               
                                                                                
  190 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
