      SUBROUTINE RS4IND(i,j,k,l,sint,npos,                                                         
     &                  nit,oiir,idk,ibias,oij)                                                              
C                                                                               
C DIESE SUBROUTINE BEARBEITET DIE 4-INDEX - INTEGRALE                           
C
C     i,j,k,l : Indexe des Integrals
C     sint : Ist das Integral selber
C     npos : Ist die Position im FT31  
C                                                                               
      IMPLICIT REAL*8 (A-H,P-Z),INTEGER*2 (O)                                   
C 
      real*8 sint
      integer*4 i,j,k,l,npos,mdone
      integer*4 mxo,mxorb             
      integer*4 nit,oiir,idk,oij
      INTEGER*4 Ibias                        
      dimension nit(*),oiir(*),idk(*),ibias(*),oij(*)
CSUT  PARAMETER (MXORB=  99)                                                    
      PARAMETER (MXORB= 400)                                                    
      PARAMETER (MXO  = 401)                                                    
      PARAMETER (MDONE = MXORB*(MXORB+1)/2+1)                                   
      PARAMETER (NDIMH=2000)                                                    
CC                                                                              
*     COMMON /I2/ OIJ(8),OXX(8),OISYM(8),OJSYM(36),OIIR(MXORB),                 
*    *OMJ(8),OKJ(8),ONJ(8),OLJ(8),OIBAL(8),                                     
*    *ONTIL(8),ONBAL(9),ODUM1(3),ONCOMP(100),                                   
*    *OITIL(8),OMCOMP(100)                                                      
*                                                                               
*     INTEGER*4 IDK,ICODE,IYXE,JOD,NIT,II,JG                                    
*     COMMON /I4/ IDK(MXO),NIT(667)                                       
      DOUBLE PRECISION BOB,CORE                                                 
      REAL*4 COUL,EXC                                                           
      COMMON /CONEI/ COUL(MDONE),EXC(MDONE),BOB(MDONE),CORE                     
C                                                                               
C --- AUSDRUCK DER NIT-FELD-INFORMATION                                         
*     WRITE(6,*)'NIT(1--20)'                                                    
*     WRITE(6,*) (NIT(II),II=1,20)                                              
*     WRITE(6,*)'NIT(666),NIT(667)'                                             
*     WRITE(6,*) NIT(666),NIT(667)                                              
C     IR BIAS FOR ORBITAL NUMBERS AND IRREP OF ORBITALS                         
*     IBIAS(1)=0  
*     do i=1,8
*       ibias(i+1) = ibias(i) + i
*     end do                                                              
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=1,10)                                   
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=11,20)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=21,30)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=31,40)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=41,50)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=51,60)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=61,70)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=71,80)                                  
*     WRITE(6,*) 'OIIR=',(OIIR(III),III=81,90)                                  
C --- ABFANGEN FALLS ALLE INTEGRALE IN DEN CORE PASSEN!                         
      IF (I.EQ.J .AND. K.EQ.L) THEN                                          
          COUL(IDK(I)+K) = (SINT)                                              
C          IF (I.EQ.K) THEN                                                   
C             WRITE(6,*) 'J=K  (',I,')=J(',idk(i)+k,')'
C     .          ,SINT ,COUL(IDK(I)+K)                                
C          END IF                                                             
      END IF                                                                 
      IF (I.EQ.K .AND. J.EQ.L) THEN                                          
          EXC(IDK(I)+J) = (SINT)                                             
      END IF                                                                 
C      CALCULATE SUBLIST AND INDEX OF THE CURRENT INTEGRAL                      
      IR=OIIR(I)                                                               
      JR=OIIR(J)                                                               
      KR=OIIR(K)                                                               
      LR=OIIR(L)                                                               
C --- IP : WELCHER INTEGRALBLOCK (1--667)                                       
CSUT   IP=IDK(IR)*(IDK(IR)+1)/2 + IDK(IR)*(JR-1)+IDK(JR)+IDK(KR)+LR             
       IP=IDK(IR)+JR                                                            
       IP = IP*(IP-1)/2 + IDK(KR) + LR                                          
       IX=I-IBIAS(IR)                                                           
       JX=J-IBIAS(JR)                                                           
       KX=K-IBIAS(KR)                                                           
       LX=L-IBIAS(LR)                                                           
C                                                                               
C      CALCULATE LOCATION OF INTEGRAL IN THE QUADRUPLE BLOCK                    
C                                                                               
C                                                                               
CSUT   IJ=QSUBLI(IP)                                                            
C                                                                               
CSUT : IM MOMENT KEINE TESTS : MUSS NOCH HINZUGEFUEGT WERDEN                    
       IJ = 1                                                                   
       IF (IJ.EQ.0) THEN                                                        
        WRITE (*,*) 'INVALID INDEX COMBINATION'                                 
        WRITE (*,*) 'INDICES: ',I,J,K,L                                         
        WRITE (*,*) 'IRREPS : ',IR,JR,KR,LR                                     
        STOP 'INVALID INDEX COMBINATION'                                        
       ENDIF                                                                    
       IF (IR.EQ.JR) GOTO 60                                                    
       IF (IR.EQ.KR) GOTO 61                                                    
C                                                                               
C      IJKL                                                                     
C                                                                               
CSUT   INDEX=(IX-1)*OIJ(JR)*OIJ(KR)*OIJ(LR)+(JX-1)*OIJ(KR)*OIJ(LR)+             
CSUT * (KX-1)*OIJ(LR)+LX                                                        
       INDEX=OIJ(LR)*(KX-1+OIJ(KR)*(JX-1+OIJ(JR)*(IX-1)))+ LX                   
       GOTO 70                                                                  
60     IF (IR.EQ.KR) GOTO 62                                                    
C                                                                               
C      IIJJ                                                                     
CSUT   INDEX=(IDK(IX)+JX-1)*IDK(OIJ(KR)+1)+IDK(KX)+LX                           
       INDEX=(IDK(IX)+JX-1)*IDK(OIJ(KR)+1)+IDK(KX)+LX                           
       GOTO 70                                                                  
62     CONTINUE                                                                 
C                                                                               
C      IIII                                                                     
C                                                                               
       IO=IDK(IX)                                                               
       INDEX=IO*(IO+1)/2                                                        
       INDEX=INDEX+IO*(JX-1)+IDK(JX)+IDK(KX)+LX                                 
       GOTO 70                                                                  
61     CONTINUE                                                                 
C                                                                               
C      IJIJ                                                                     
C  
CDEBUG write(6,*) 'Fall ijij',ix,jx,lx,kx                                                                        
cold       INDEX = IDK(IX) + JX + 1                                                 
cold       INDEX = INDEX*(INDEX-1)/2 + LX  
       nj = oij(jr)
       index = nj*(ix-1)+jx 
       INDEX = INDEX*(INDEX-1)/2 + nj*(kx-1)+lx
CDEBUG write(6,*) 'Index=',index                                         
CSUT   INDEX=0                                                                  
CSUT   IF(IX.NE.1) INDEX=IDK(IX-1)*OIJ(JR)*OIJ(JR)+(IX-1)*IDK(OIJ(JR)+1)        
CSUT   INDEX=INDEX + (IX-1)*(JX-1)*OIJ(JR) + IDK(JX)                            
CSU    INDEX=INDEX+(KX-1)*OIJ(JR)+LX                                            
C                                                                               
C -- GELESENES INTEGRAL IM BEREICH                                              
C -- NPOS : POSITION IM GESAMT-INTEGRALFILE                                     
   70  CONTINUE                                                                 
C      IF (INDEX.LT.NINTS) GOTO 71                                              
C       INDEX=INDEX-NINTS                                                       
C       IJ=IJ+1                                                                 
C       GOTO 70                                                                 
C 71   CONTINUE                                                                 
       NPOS = INDEX + NIT(IP)                                                   
CDEBUGWRITE(6,*) 'IRECC=',IRECC,' ABGESCHLOSSEN'                                
      RETURN                                                                    
      END                                                                       
