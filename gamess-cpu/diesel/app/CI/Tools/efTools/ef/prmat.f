                                                                                
      SUBROUTINE PREIG(IO,E,NORBS)                                              
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      include 'param'                                                           
      DIMENSION E(*)                                                            
      CHARACTER*3 MOL                                                           
                                                                                
      WRITE(io,'(/,10x,''EIGENVALUES'')')                                       
      N=6                                                                       
      NTIMES=NORBS/N                                                            
      NREST=MOD(NORBS,N)                                                        
      IF(NTIMES.EQ.0) NREST=NORBS                                               
      IE=J+NREST-1                                                              
                                                                                
      J=1                                                                       
      N2=N                                                                      
                                                                                
      DO K=1,NTIMES                                                             
      WRITE(IO,*)                                                               
      WRITE(IO,100)(I,I=J,N2)                                                   
      WRITE(IO,300)(E(I),I=J,N2)                                                
      J =J +N                                                                   
      N2=N2+N                                                                   
      ENDDO                                                                     
                                                                                
      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN                                        
      WRITE(IO,*)                                                               
      IE=J+NREST-1                                                              
      WRITE(IO,100)(I,I=J,IE)                                                   
      WRITE(IO,300)(E(I),I=J,IE)                                                
      ENDIF                              

         
                                                                                
      RETURN                                                                    
                                                                                
 100  FORMAT(' ROOT NR.: ',2X,6(4X,I4,2X))                                      
 300  FORMAT(' EIGVAL  : ',2X,6F10.3)                                           
      END                                                                       
                                                                                






      SUBROUTINE PRMAT(IUOUT,AOL,R,N,M,HEAD)                                    
      include 'param'                                                           
      REAL*8 R                                                                  
      CHARACTER*(*) HEAD                                                        
      DIMENSION R(*), AOL(*)                                                    
      CHARACTER*3 AOL                                                           
                                                                                
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED                              
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND                              
C     ((N+1)*N)/2 WHEN M IS ZERO                                                
                                                                                
      WRITE(IUOUT,1001) HEAD                                                    
      NKPB=6                                                                    
c     ghost prepare for F90:
c     IF(M)10,10,80                                                             
      IF(M.lt.0) goto 10
      if(m.eq.0) goto 10
      if(m.gt.0) goto 80                                                             
C                                                                               
   10 CONTINUE                                                                  
      IBL=N/NKPB                                                                
      IR=N-IBL*NKPB                                                             
      J1=1                                                                      
      K1S=1                                                                     
      KD=0                                                                      
      IF(IBL.EQ.0) GO TO 50                                                     
      J2=NKPB                                                                   
      DO 40 I=1,IBL                                                             
      IF(M.EQ.0) WRITE(IUOUT,1002)(J,AOL(J),J=J1,J2)                            
      IF(M.NE.0) WRITE(IUOUT,1004)(J,J=J1,J2)                                   
      K1=K1S                                                                    
      K2=K1                                                                     
      KK=0                                                                      
      DO 20 J=J1,J2                                                             
      WRITE(IUOUT,1003)J,AOL(J),(R(K),K=K1,K2)                                  
      KK=KK+1                                                                   
      K1=K1+KD+KK                                                               
   20 K2=K1+KK                                                                  
      J1=J1+NKPB                                                                
      IF(J1.GT.N) RETURN                                                        
      J2=J2+NKPB                                                                
      K2=K1-1                                                                   
      K1=K2+1                                                                   
      K2=K1+(NKPB-1)                                                            
      K1S=K2+1                                                                  
      KK=KD+NKPB                                                                
      DO 30 J=J1,N                                                              
      WRITE(IUOUT,1003)J,AOL(J),(R(K),K=K1,K2)                                  
      KK=KK+1                                                                   
      K1=K1+KK                                                                  
   30 K2=K2+KK                                                                  
   40 KD=KD+NKPB                                                                
   50 IF(IR.EQ.0) GO TO 70                                                      
      K1=K1S                                                                    
      J2=J1+IR-1                                                                
      KK=0                                                                      
      K2=K1                                                                     
      IF(M.EQ.0)WRITE(IUOUT,1002)(J,AOL(J),J=J1,J2)                             
      IF(M.NE.0)WRITE(IUOUT,1004)(J,J=J1,J2)                                    
      WRITE(IUOUT,1003)                                                         
      DO 60 J=J1,J2                                                             
      WRITE(IUOUT,1003)J,AOL(J),(R(K),K=K1,K2)                                  
      KK=KK+1                                                                   
      K1=K1+KD+KK                                                               
   60 K2=K1+KK                                                                  
   70 RETURN                                                                    
   80 IBL=M/NKPB                                                                
      IR=M-IBL*NKPB                                                             
      I2=0                                                                      
      K2=0                                                                      
      IF(IBL.EQ.0) GO TO 100                                                    
      DO 90 I=1,IBL                                                             
      I1=(I-1)*N*NKPB+1                                                         
      I2=I1+(NKPB-1)*N                                                          
      K1=K2+1                                                                   
      K2=K1+(NKPB-1)                                                            
      IF(M.EQ.0)WRITE(IUOUT,1002)(K,AOL(K),K=K1,K2)                             
      IF(M.NE.0)WRITE(IUOUT,1004)(K,K=K1,K2)                                    
      DO 90 J=1,N                                                               
      WRITE(IUOUT,1003)J,AOL(J),(R(IJ),IJ=I1,I2,N)                              
      I1=I1+1                                                                   
   90 I2=I1+(NKPB-1)*N                                                          
  100 IF(IR.EQ.0) GO TO 120                                                     
      I1=IBL*N*NKPB+1                                                           
      I2=I1+(IR-1)*N                                                            
      K1=K2+1                                                                   
      K2=M                                                                      
      IF(M.EQ.0)WRITE(IUOUT,1002)(K,AOL(K),K=K1,K2)                             
      IF(M.NE.0)WRITE(IUOUT,1004)(K,K=K1,K2)                                    
      WRITE(IUOUT,1003)                                                         
      DO 110 J=1,N                                                              
      WRITE(IUOUT,1003)J,AOL(J),(R(IJ),IJ=I1,I2,N)                              
      I1=I1+1                                                                   
      I2=I1+(IR-1)*N                                                            
  110 CONTINUE                                                                  
  120 WRITE(IUOUT,1003)                                                         
      RETURN                                                                    
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)                                          
 1002 FORMAT(/,' ',11X,6(2X,I4,1X,A3),/)                                        
 1003 FORMAT(' ',I4,1X,A3,3X,6F11.3)                                            
 1004 FORMAT(/,' ',11X,6(2X,I4,7X),/)                                           
      END                                                                       
                                                                                
