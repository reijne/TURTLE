C     *****************************************************************         
                                                                                
      SUBROUTINE AXIS(NUMAT,NAT,COORD)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      COMMON /KEYWRD/ KEYWRD                                                    
      CHARACTER*241 KEYWRD                                                      
      COMMON /MASS/ AMS(107)                                                    
                                                                                
      DIMENSION T(6), ROT(3), XYZMOM(3), EIG(3), EVEC(3,3)                      
      DIMENSION COORD(3,NUMATM), COORD1(3,NUMATM), NAT(*)                       
      DIMENSION X(NUMATM), Y(NUMATM), Z(NUMATM)                                 
      DATA T /6*0.D0/                                                           
************************************************************************        
*     CONST1 =  10**40/(N*A*A)                                                  
*               N = AVERGADRO'S NUMBER                                          
*               A = CM IN AN ANGSTROM                                           
*               10**40 IS TO ALLOW UNITS TO BE 10**(-40)GRAM-CM**2              
*                                                                               
************************************************************************        
      CONST1 = 1.66053D0                                                        
************************************************************************        
*                                                                               
*     CONST2 = CONVERSION FACTOR FROM ANGSTROM-AMU TO CM**(-1)                  
*                                                                               
*            = (PLANCK'S CONSTANT*N*10**16)/(8*PI*PI*C)                         
*            = 6.62618*10**(-27)[ERG-SEC]*6.02205*10**23*10**16/                
*              (8*(3.1415926535)**2*2.997925*10**10[CM/SEC])                    
*                                                                               
************************************************************************        
      CONST2=16.8576522D0                                                       
C    FIRST WE CENTRE THE MOLECULE ABOUT THE CENTRE OF GRAVITY,                  
C    THIS DEPENDS ON THE ISOTOPIC MASSES, AND THE CARTESIAN GEOMETRY.           
C                                                                               
      SUMW=1.D-20                                                               
      SUMWX=0.D0                                                                
      SUMWY=0.D0                                                                
      SUMWZ=0.D0                                                                
      AMS(99)=10.0D0                                                            
                                                                                
      DO I=1,NUMAT                                                              
         X(I)=COORD(1,I)                                                        
         Y(I)=COORD(2,I)                                                        
         Z(I)=COORD(3,I)                                                        
      ENDDO                                                                     
                                                                                
C                                                                               
         DO 10 I=1,NUMAT                                                        
            ATMASS=AMS(NAT(I))                                                  
C           write(*,*) i,atmass,nat(i)                                          
            SUMW=SUMW+ATMASS                                                    
            SUMWX=SUMWX+ATMASS*COORD(1,I)                                       
            SUMWY=SUMWY+ATMASS*COORD(2,I)                                       
            SUMWZ=SUMWZ+ATMASS*COORD(3,I)                                       
   10    CONTINUE                                                               
C                                                                               
C     WRITE(6,'(/10X,''MOLECULAR WEIGHT ='',F8.2,/)')                           
C    2MIN(99999.99D0,SUMW)                                                      
      SUMWX=SUMWX/SUMW                                                          
      SUMWY=SUMWY/SUMW                                                          
      SUMWZ=SUMWZ/SUMW                                                          
      DO 30 I=1,NUMAT                                                           
         X(I)=COORD(1,I)-SUMWX                                                  
         Y(I)=COORD(2,I)-SUMWY                                                  
   30 Z(I)=COORD(3,I)-SUMWZ                                                     
************************************************************************        
*                                                                               
*    MATRIX FOR MOMENTS OF INERTIA IS OF FORM                                   
*                                                                               
*           |   Y**2+Z**2                         |                             
*           |    -Y*X       Z**2+X**2             | -I =0                       
*           |    -Z*X        -Z*Y       X**2+Y**2 |                             
*                                                                               
************************************************************************        
C                                                                               
C$DOIT ASIS                                                                     
      DO 40 I=1,6                                                               
   40 T(I)=DBLE(I)*1.0D-10                                                      
C                                                                               
         DO 50 I=1,NUMAT                                                        
            ATMASS=AMS(NAT(I))                                                  
            T(1)=T(1)+ATMASS*(Y(I)**2+Z(I)**2)                                  
            T(2)=T(2)-ATMASS*X(I)*Y(I)                                          
            T(3)=T(3)+ATMASS*(Z(I)**2+X(I)**2)                                  
            T(4)=T(4)-ATMASS*Z(I)*X(I)                                          
            T(5)=T(5)-ATMASS*Y(I)*Z(I)                                          
            T(6)=T(6)+ATMASS*(X(I)**2+Y(I)**2)                                  
   50    CONTINUE                                                               
C                                                                               
      CALL RSP(T,3,3,EIG,EVEC)                                                  
C$DOIT ASIS                                                                     
         DO 70 I=1,3                                                            
            IF(EIG(I).LT.3.D-4) THEN                                            
               EIG(I)=0.D0                                                      
               ROT(I)=0.D0                                                      
            ELSE                                                                
               ROT(I)=CONST2/EIG(I)                                             
            ENDIF                                                               
   70    XYZMOM(I)=EIG(I)*CONST1                                                
C        WRITE(6,'(10X,''A ='',F12.6,''   B ='',F12.6,                          
C    1''   C ='',F12.6,/)')(ROT(I),I=1,3)                                       
C                                                                               
C   NOW TO ORIENT THE MOLECULE SO THE CHIRALITY IS PRESERVED                    
C                                                                               
      SUM=EVEC(1,1)*(EVEC(2,2)*EVEC(3,3)-EVEC(3,2)*EVEC(2,3)) +                 
     1    EVEC(1,2)*(EVEC(2,3)*EVEC(3,1)-EVEC(2,1)*EVEC(3,3)) +                 
     2    EVEC(1,3)*(EVEC(2,1)*EVEC(3,2)-EVEC(2,2)*EVEC(3,1))                   
      IF( SUM .LT. 0) THEN                                                      
C$DOIT ASIS                                                                     
         DO 80 J=1,3                                                            
   80    EVEC(J,1)=-EVEC(J,1)                                                   
      ENDIF                                                                     
      DO 90 I=1,NUMAT                                                           
         COORD(1,I)=X(I)                                                        
         COORD(2,I)=Y(I)                                                        
         COORD(3,I)=Z(I)                                                        
   90 CONTINUE                                                                  
                                                                                
      DO 120 I=1,NUMAT                                                          
         DO 120 J=1,3                                                           
            SUM=0.D0                                                            
            DO 110 K=1,3                                                        
  110       SUM=SUM+COORD(K,I)*EVEC(K,J)                                        
  120 COORD1(J,I)=SUM                                                           
                                                                                
C default: interchange x and z                                                  
                                                                                
      do i=1,numat                                                              
             XXX=COORD1(3,I)                                                    
             ZZZ=COORD1(1,I)                                                    
             COORD1(1,I)=XXX                                                    
             COORD1(3,I)=ZZZ                                                    
      ENDDO                                                                     
                                                                                
      i1=1                                                                      
      i2=2                                                                      
      i3=3                                                                      
                                                                                
      if(index(keywrd,'XZ').ne.0) then                                          
         I1=3                                                                   
         I2=2                                                                   
         I3=1                                                                   
      endif                                                                     
      if(index(keywrd,'XY').ne.0) then                                          
         I1=2                                                                   
         I2=1                                                                   
         I3=3                                                                   
      endif                                                                     
      if(index(keywrd,'YZ').ne.0) then                                          
         I1=1                                                                   
         I2=3                                                                   
         I3=2                                                                   
      endif                                                                     
      do i=1,numat                                                              
             COORD(1,I)=COORD1(i1,I)                                            
             COORD(2,I)=COORD1(i2,I)                                            
             COORD(3,I)=COORD1(i3,I)                                            
      ENDDO                                                                     
                                                                                
                                                                                
      RETURN                                                                    
      END                                                                       
