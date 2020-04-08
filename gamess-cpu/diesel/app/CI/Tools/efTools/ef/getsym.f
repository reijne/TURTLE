      SUBROUTINE GETSYM                                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),                     
     1                LOCDEP(MAXPAR)                                            
C***********************************************************************        
C                                                                               
C   GETSYM READS IN THE SYMMETRY DEPENDENCE RELATIONSHIPS.                      
C                                                                               
C   ON EXIT     NDEP    = NUMBER OF SYMMETRY RELATIONS.                         
C               LOCPAR  = ARRAY OF REFERENCE FUNCTION INDICES.                  
C               IDEPFN  = ARRAY OF REFERENCE ATOM LOCATIONS.                    
C               LOCDEP  = ARRAY OF DEPENDENT ATOM LOCATIONS.                    
C                                                                               
C***********************************************************************        
C                                                                               
C     LOCDEP IS THE ATOM WHOSE COORDINATES DEPEND ON THE COORDINATES OF         
C     LOCPAR.                                                                   
C     LOCPAR IS THE ATOM WHOSE COORDINATES ARE USED TO CALCULATE THOSE          
C     OF LOCDEP                                                                 
C     IDEPFN POINTS TO THE PARTICULAR FUNCTION TO BE USED (SEE NDDO)            
C                                                                               
C***********************************************************************        
      DIMENSION IVALUE(40),VALUE(40)                                            
      CHARACTER  TEXT(18)*60, LINE*80                                           
      SAVE TEXT                                                                 
      DATA TEXT/                                                                
     1' BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ',           
     2' BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ',           
     3' DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE',           
     4' DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ',           
     5' DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ',           
     6' DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ',           
     7' DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ',           
     8' DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ',           
     9' DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ',           
     1' DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ',           
     2' DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ',           
     3' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',           
     4' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',           
     5' DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ',           
     6' BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ',           
     7' BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ',           
     8' BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ',           
     9' THE USER HAS TO SUPPLY THIS FUNCTION IN DEPVAR             '/           
C                                                                               
C TITLE OUTPUT                                                                  
      WRITE (6,10)                                                              
c     ghost prepare for F90:
c   10 FORMAT (///5X,25HPARAMETER DEPENDENCE DATA//                              
   10 FORMAT (///5X,'PARAMETER DEPENDENCE DATA'//
     1'        REFERENCE ATOM      FUNCTION NO.    DEPENDENT ATOM(S)')          
C                                                                               
C INPUT SYMMETRY : FUNCTION, REFERANCE PARAMETER, AND DEPENDENT ATOMS           
C                                                                               
      NDEP=0                                                                    
   20 READ(5,'(A)',END=70) LINE                                                 
      CALL NUCHAR(LINE,VALUE,NVALUE)                                            
C   INTEGER VALUES                                                              
      DO 30 I=1,NVALUE                                                          
   30 IVALUE(I)=VALUE(I)                                                        
C   FILL THE LOCDEP ARRAY                                                       
      IF(NVALUE.EQ.0.OR.IVALUE(3).EQ.0) GO TO 70                                
      DO 40 I=3,NVALUE                                                          
         IF(IVALUE(I).EQ.0) GOTO 50                                             
         NDEP=NDEP+1                                                            
         LOCDEP(NDEP)=IVALUE(I)                                                 
         LOCPAR(NDEP)=IVALUE(1)                                                 
         IDEPFN(NDEP)=IVALUE(2)                                                 
   40 CONTINUE                                                                  
   50 LL=I-1                                                                    
      WRITE(6,60)IVALUE(1),IVALUE(2),(IVALUE(J),J=3,LL)                         
   60 FORMAT(I13,I19,I14,11I3,10(/,43X,12I3))                                   
      GO TO 20                                                                  
C                                                                               
C CLEAN UP                                                                      
   70 CONTINUE                                                                  
      WRITE(6,80)                                                               
   80 FORMAT(/10X,'   DESCRIPTIONS OF THE FUNCTIONS USED',/)                    
      DO 120 J=1,18                                                             
         DO 90 I=1,NDEP                                                         
            IF(IDEPFN(I).EQ.J) GOTO 100                                         
   90    CONTINUE                                                               
         GOTO 120                                                               
  100    WRITE(6,110)J,TEXT(J)                                                  
  110    FORMAT(I4,5X,A)                                                        
  120 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
