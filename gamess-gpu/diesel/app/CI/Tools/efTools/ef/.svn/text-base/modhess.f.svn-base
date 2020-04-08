      implicit real*8 (a-h,o-z)
      dimension hess(1000,1000),xx(10)
      character*80 a8

      open(unit=4,file='hessian')                                               
      rewind 4                                                                  

      read (4,'(A)') a8                                                        
      call readl(80,a8,xx,nn)
      nvar=idint(xx(1))
      read (4,'(A)') a8                                                        
      call readl(80,a8,xx,nn)
      eigval=xx(1)

      do i=1,nvar                                                               
         read (4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 

      write(*,*) 'modify element i,j'       

 10   read (*,*) mi,mj                      
      if(mj.gt.mi) then
         write(*,*) 'mi should be larger than mj'
         goto 10
      endif

      write(*,*) 'element i,j',mi,mj,'=',hess(mi,mj)       
      write(*,*) 'new value        '       
      read(*,*) hessnew

      hess(mi,mj)=hessnew

      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      write(4,*) 'hessian dim=',nvar                                            
      write(4,*) 'lowest eigenvalue=',eigval                                 
      do i=1,nvar                                                               
         write(4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 

      end
C     *****************************************************************         
                                                                                
      SUBROUTINE READL(NL,A1,X,N)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*1 A1(*)                                                         
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READAA(NL,A1,IS,IB,IE)                                               
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READAA(NL,A,ISTART,IEND,IEND2)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READAA                                                             
      CHARACTER*1 A(NL)                                                         
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
      DO 10 J=ISTART,NL-1                                                       
         N=ICHAR(A(J))                                                          
         M=ICHAR(A(J+1))                                                        
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO                            
     1 .OR. M.EQ.IDOT)) GOTO 20                                                 
   10 CONTINUE                                                                  
      READAA=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,NL                                                              
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
      DO 50 II=I+1,NL                                                           
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
      READAA= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READAA=READAA*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       
