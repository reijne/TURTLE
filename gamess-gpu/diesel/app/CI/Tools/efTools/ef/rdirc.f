C     *****************************************************************         
                                                                                
      PROGRAM MAIN                                                              
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      character*80 a,b                                                          
      dimension xx(10),c(3,100),nr(4),t(10000)                                  
                                                                                
      f=0.529177260d0                                                           
      open(unit=1,file='ircgeo')                                                
      open(unit=2,file='ircvelo')                                               
                                                                                
      read(1,*) ngeo, numat                                                     
      read(1,'(A)') b                                                           
      call readl(80,b,xx,nn)                                                    
      do i=1,nn                                                                 
         nr(i)=idint(xx(i))                                                     
      enddo                                                                     
                                                                                
      k=0                                                                       
  10  read(2,'(A)',end=100) a                                                   
      if(index(a,'IRC').ne.0) then                                              
         call readl(80,a,xx,n)                                                  
         k=k+1                                                                  
         t(k)=xx(2)                                                             
      endif                                                                     
      goto 10                                                                   
 100  continue                                                                  
                                                                                
      do i=1,ngeo                                                               
         do j=1,numat                                                           
            read(1,'(A)') a                                                     
            call readl(80,a,xx,n)                                               
            c(3,j)=xx(n)*f                                                      
            c(2,j)=xx(n-1)*f                                                    
            c(1,j)=xx(n-2)*f                                                    
         enddo                                                                  
         if(nn.eq.2) then                                                       
            rr=dsqrt( (c(1,nr(1))-c(1,nr(2)))**2                                
     .               +(c(2,nr(1))-c(2,nr(2)))**2                                
     .               +(c(3,nr(1))-c(3,nr(2)))**2)                               
            write(*,*) t(i),rr                                                  
         endif                                                                  
         if(nn.eq.3) then                                                       
            n1=nr(1)                                                            
            n2=nr(2)                                                            
            n3=nr(3)                                                            
            call bangle(c,n1,n2,n3,ang)                                         
            write(*,*) t(i),ang*180/3.14159265358979D0                          
         endif                                                                  
         if(nn.eq.4) then                                                       
            n1=nr(1)                                                            
            n2=nr(2)                                                            
            n3=nr(3)                                                            
            n4=nr(4)                                                            
            call dihed (c,n1,n2,n3,n4,ang)                                      
            ang=ang*180/3.14159265358979D0                                      
            if(ang.gt.180) ang=360-ang                                          
            write(*,*) t(i),ang                                                 
         endif                                                                  
      enddo                                                                     
                                                                                
                                                                                
                                                                                
      END                                                                       
                                                                                
C     *****************************************************************         
                                                                                
      SUBROUTINE READL(L,A1,X,N)                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*1 A1(*)                                                         
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READA2(L,A1,IS,IB,IE)                                                
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READA2(L,A,ISTART,IEND,IEND2)                                    
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READA2                                                             
      CHARACTER*1 A(*)                                                          
      NINE=ICHAR('9')                                                           
      IZERO=ICHAR('0')                                                          
      MINUS=ICHAR('-')                                                          
      IDOT=ICHAR('.')                                                           
      ND=ICHAR('D')                                                             
      NE=ICHAR('E')                                                             
      IBL=ICHAR(' ')                                                            
      IEND=0                                                                    
      IEND2=0                                                                   
      IDIG=0                                                                    
      C1=0                                                                      
      C2=0                                                                      
      ONE=1.D0                                                                  
      X = 1.D0                                                                  
      DO 10 J=ISTART,L-1                                                        
         N=ICHAR(A(J))                                                          
         M=ICHAR(A(J+1))                                                        
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO                            
     1 .OR. M.EQ.IDOT)) GOTO 20                                                 
   10 CONTINUE                                                                  
      READA2=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,L                                                               
         N=ICHAR(A(I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IDIG=0                                                                    
      DO 50 II=I+1,L                                                            
         N=ICHAR(A(II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   50 CONTINUE                                                                  
C                                                                               
C PUT THE PIECES TOGETHER                                                       
C                                                                               
   60 CONTINUE                                                                  
      READA2= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,L                                                            
         N=ICHAR(A(J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,L                                                             
         N=ICHAR(A(I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READA2=READA2*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       
                                                                                
      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION XYZ(3,*)                                                        
*********************************************************************           
*                                                                               
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE                     
*        CARTESIAN COORDINATES ARE IN XYZ.                                      
*                                                                               
*********************************************************************           
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+                                            
     1       (XYZ(2,I)-XYZ(2,J))**2+                                            
     2       (XYZ(3,I)-XYZ(3,J))**2                                             
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+                                            
     1       (XYZ(2,J)-XYZ(2,K))**2+                                            
     2       (XYZ(3,J)-XYZ(3,K))**2                                             
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+                                            
     1       (XYZ(2,I)-XYZ(2,K))**2+                                            
     2       (XYZ(3,I)-XYZ(3,K))**2                                             
      XY = SQRT(D2IJ*D2JK)                                                      
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY                                      
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0                                           
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0                                         
      ANGLE = ACOS( TEMP )                                                      
      RETURN                                                                    
      END                                                                       
      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION XYZ(3,*)                                                        
*********************************************************************           
*                                                                               
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,               
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS                   
*            ARE IN ARRAY XYZ.                                                  
*                                                                               
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME              
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.                          
*                                                                               
*********************************************************************           
      DATA PI/3.14159265358979D0/                                               
      XI1=XYZ(1,I)-XYZ(1,K)                                                     
      XJ1=XYZ(1,J)-XYZ(1,K)                                                     
      XL1=XYZ(1,L)-XYZ(1,K)                                                     
      YI1=XYZ(2,I)-XYZ(2,K)                                                     
      YJ1=XYZ(2,J)-XYZ(2,K)                                                     
      YL1=XYZ(2,L)-XYZ(2,K)                                                     
      ZI1=XYZ(3,I)-XYZ(3,K)                                                     
      ZJ1=XYZ(3,J)-XYZ(3,K)                                                     
      ZL1=XYZ(3,L)-XYZ(3,K)                                                     
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS                              
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)                                          
      COSA=ZJ1/DIST                                                             
      IF(COSA.GT.1.0D0) COSA=1.0D0                                              
      IF(COSA.LT.-1.0D0) COSA=-1.0D0                                            
      DDD=1.0D0-COSA**2                                                         
      IF(DDD.LE.0.0) GO TO 10                                                   
      YXDIST=DIST* SQRT(DDD)                                                    
      IF(YXDIST.GT.1.0D-9) GO TO 20                                             
   10 CONTINUE                                                                  
      XI2=XI1                                                                   
      XL2=XL1                                                                   
      YI2=YI1                                                                   
      YL2=YL1                                                                   
      COSTH=COSA                                                                
      SINTH=0.D0                                                                
      GO TO 30                                                                  
   20 COSPH=YJ1/YXDIST                                                          
      SINPH=XJ1/YXDIST                                                          
      XI2=XI1*COSPH-YI1*SINPH                                                   
      XJ2=XJ1*COSPH-YJ1*SINPH                                                   
      XL2=XL1*COSPH-YL1*SINPH                                                   
      YI2=XI1*SINPH+YI1*COSPH                                                   
      YJ2=XJ1*SINPH+YJ1*COSPH                                                   
      YL2=XL1*SINPH+YL1*COSPH                                                   
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS                  
      COSTH=COSA                                                                
      SINTH=YJ2/DIST                                                            
   30 CONTINUE                                                                  
      YI3=YI2*COSTH-ZI1*SINTH                                                   
      YL3=YL2*COSTH-ZL1*SINTH                                                   
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)                                          
      IF (ANGLE .LT. 0.) ANGLE=2.0D0*PI+ANGLE                                   
      IF (ANGLE .GE. 2.0d0*PI    ) ANGLE=0.D0                                   
      RETURN                                                                    
      END                                                                       
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
**********************************************************************          
*                                                                               
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),              
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.                             
*                                                                               
**********************************************************************          
      DATA PI/3.14159265358979D0/                                               
      ZERO=1.0D-6                                                               
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10                        
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10                        
      ANORM=1.0D0/ SQRT(A1**2+A2**2)                                            
      BNORM=1.0D0/ SQRT(B1**2+B2**2)                                            
      A1=A1*ANORM                                                               
      A2=A2*ANORM                                                               
      B1=B1*BNORM                                                               
      B2=B2*BNORM                                                               
      SINTH=(A1*B2)-(A2*B1)                                                     
      COSTH=A1*B1+A2*B2                                                         
      IF(COSTH.GT.1.0D0) COSTH=1.0D0                                            
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0                                          
      RCOS= ACOS(COSTH)                                                         
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10                                         
      IF(SINTH.GT.0.D0) RCOS=2.0D0*PI-RCOS                                      
      RCOS=-RCOS                                                                
      RETURN                                                                    
   10 RCOS=0.0D0                                                                
      RETURN                                                                    
      END                                                                       
