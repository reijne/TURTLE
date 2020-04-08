      SUBROUTINE KOTZ                                                           
C                                                                               
C  DIESE SUBROUTINE SCHREIBT DIE EINELEKTRONEN INTEGRAL RAUS                    
C     
      implicit real*8 (a-h,o-z) 
      include 'common.inc'
C
      integer*4 mdon,mxor
      parameter (mxor=400)
      parameter (mdon = mxor*(mxor+1)/2+1)                                                                  
      DOUBLE PRECISION BOB                                                 
      REAL*4 COUL,EXC                                                           
      COMMON /CONEI/ COUL(MDON),EXC(MDON),BOB(MDON)                     
      REAL*4 STON(2000) 
C Korrektur des Kernpotentials
      CORE = CORE + POTNUC
C                                                        
      WRITE(31) BOB                                                             
      WRITE(31) COUL                                                            
      WRITE(31) EXC                                                             
      WRITE(6,*)                                                                
      WRITE(6,*) 'EINELEKTRONENINTEGRALE GESCHRIEBEN'                           
      WRITE(6,*)                                                                
      WRITE(31) PotNuc                                                          
      WRITE(6,*)                                                                
      WRITE(6,*) 'potnuc =',potnuc
      WRITE(6,*)                                                                
C
CDEBUGWRITE(6,'(5(1x,f12.8))') 'BOB :',(BOB(I),I=1,100)                                        
CDEBUGWRITE(6,'(5(1x,f12.8))') 'COUL:',(COUL(I),I=1,100)                                       
CDEBUGWRITE(6,'(5(1x,f12.8))') 'EXC :',(EXC(I),I=1,100)                                        
      REWIND(31)                                                                
      READ(31)                                                                  
      READ(31)                                                                  
      READ(31) STON                                                             
CDEBUGWRITE(6,*) 'STON:',(STON(I),I=1,10)                                       
      RETURN                                                                    
      END                                                                       
