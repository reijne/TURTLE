CCCCCCCCCCCCCCCCCC     MATHEMATICAL PACKAGE     CCCCCCCCCCCCCCCCCCCCCCC         
CCCCCCC    MOST OF THESE ROUTINE ARE FULLY VECTORIZED ON CRAY-1  CCCCCC         
CCCCCCC    THEY ARE ROUGHLY RESPONSIBLE OF 70% OF THE CPU TIME   CCCCCC         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
      DOUBLE PRECISION FUNCTION SDOT (N,X,IX,Y,IY)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION X(*),Y(*)                                                       
C     SDOT=DOT PRODUCT OF VECTOR X, STEP IX, BY VECTOR Y, STEP IY,              
C     N ELEMENTS.                                                               
C     SIMULATE ROUTINE ON CRAY (SAME NAME AND CALLING SEQUENCE).                
      J=1                                                                       
      SDOT=0.D0                                                                 
      DO 10 I=1,(N-1)*IX+1,IX                                                   
         SDOT=SDOT+X(I)*Y(J)                                                    
   10 J=J+IY                                                                    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SCOPY (N,X,IX,Y,IY)                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     COPY VECTOR X, STEP IX ONTO VECTOR Y, STEP IY, N ELEMENTS.                
C     SIMULATE ROUTINE ON CRAY (SAME NAME AND CALLING SEQUENCE).                
      DIMENSION X(*),Y(*)                                                       
      I=1                                                                       
      DO 10 J=1,IY*(N-1)+1,IY                                                   
         Y(J)=X(I)                                                              
   10 I=I+IX                                                                    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SAXPY(N,A,X,IX,Y,IY)                                           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     VECTOR INCREMENT Y=Y+A*X WITH X & Y VECTORS OF LENGTH N, A SCALAR.        
C     IX STEP OF X, IY STEP OF Y.                                               
C     SIMULATE ROUTINE ON CRAY (SAME NAME AND CALLING SEQUENCE).                
      DIMENSION X(*),Y(*)                                                       
      I=1                                                                       
      DO 10 J=1,IY*(N-1)+1,IY                                                   
         Y(J)=Y(J)+A*X(I)                                                       
   10 I=I+IX                                                                    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE MXM(A,NAR,B,NBR,C,NCC)                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     RECTANGULAR MATRIX PRODUCT C=A*B.                                         
C     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.                            
C     SIMULATE ROUTINE ON CRAY (SAME NAME AND CALLING SEQUENCE).                
      DIMENSION A(NAR,NBR),B(NBR,NCC),C(NAR,NCC)                                
      DO 20 J=1,NCC                                                             
         DO 10 I=1,NAR                                                          
   10    C(I,J)=0.D0                                                            
         DO 20 K=1,NBR                                                          
            DO 20 I=1,NAR                                                       
   20 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE MXMT (A,NAR,B,NBR,C,NCC)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     MATRIX PRODUCT C(NAR,NCC) = A(NAR,NBR) * (B(NCC,NBR))'                    
C     ALL MATRICES RECTANGULAR , PACKED.                                        
      DIMENSION A(NAR,NBR),B(NCC,NBR),C(NAR,NCC)                                
      DO 20 J=1,NCC                                                             
         DO 10 I=1,NAR                                                          
   10    C(I,J)=0.D0                                                            
         DO 20 K=1,NBR                                                          
            DO 20 I=1,NAR                                                       
   20 C(I,J)=C(I,J)+A(I,K)*B(J,K)                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE MTXM (A,NAR,B,NBR,C,NCC)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     MATRIX PRODUCT C(NAR,NCC) = (A(NBR,NAR))' * B(NBR,NCC)                    
C     ALL MATRICES RECTANGULAR , PACKED.                                        
      DIMENSION A(NBR,NAR),B(NBR,NCC),C(NAR,NCC)                                
      DO 20 J=1,NCC                                                             
         DO 10 I=1,NAR                                                          
   10    C(I,J)=0.D0                                                            
         DO 20 K=1,NBR                                                          
            DO 20 I=1,NAR                                                       
   20 C(I,J)=C(I,J)+A(K,I)*B(K,J)                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE MTXMC (A,NAR,B,NBR,C)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     MATRIX PRODUCT C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR)                    
C     A AND B RECTANGULAR , PACKED,                                             
C     C LOWER LEFT TRIANGLE ONLY, PACKED IN CANONICAL ORDER.                    
      DIMENSION A(NBR,NAR),B(NBR,NAR),C(*)                                      
C  NOTE ... THIS IS THE BEST VERSION ON CRAY 1.                                 
      L=1                                                                       
      DO 10 I=1,NAR                                                             
         CALL MXM (A(1,I),1,B,NBR,C(L),I)                                       
   10 L=L+I                                                                     
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SUPDOT(S,H,G,N,IG)                                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     (S)=(H)*(G) WITH  H  IN PACKED FORM (CANONICAL ORDER).                    
C     IG IS THE INCREMENT FOR THE VECTOR G.                                     
      DIMENSION S(*),H(*),G(*)                                                  
C     CRAY-1 VERSION                                                            
CCC      K=1                                                                    
CCC      L=1                                                                    
CCC      DO 10 I=1,N                                                            
CCC      S(I)=SDOT(I,H(K),1,G,IG,I)                                             
CCC      IF(I.GT.1) THEN                                                        
CCC         L=L+IG                                                              
CCC         CALL SAXPY(I-1,G(L),H(K),1,S,1)                                     
CCC      ENDIF                                                                  
CCC   10 K=K+I                                                                  
CCC      RETURN                                                                 
CCC      END                                                                    
C     SCALAR VERSION OK WITH IG=1 ONLY.                                         
      K=0                                                                       
      DO 20 I=1,N                                                               
         SUM=0.D0                                                               
         DO 10 J=1,I                                                            
   10    SUM=SUM+G(J)*H(K+J)                                                    
         S(I)=SUM                                                               
   20 K=K+I                                                                     
      IF (N.EQ.1) RETURN                                                        
      K=1                                                                       
      DO 40 I=2,N                                                               
         GI=G(I)                                                                
         DO 30 J=1,I-1                                                          
   30    S(J)=S(J)+H(K+J)*GI                                                    
   40 K=K+I                                                                     
      RETURN                                                                    
      END                                                                       
