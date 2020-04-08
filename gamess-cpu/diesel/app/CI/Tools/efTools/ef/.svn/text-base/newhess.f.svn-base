      PROGRAM MAIN                                                              
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      character*80 A8                                                           
      dimension xx(20),h1(100,100),h2(100,100),ip(100),ipos(100),it(100)        
      dimension xinit(3)                                                        
      data xinit /2500.0d0,1200.0d0,1000.0d0/                                   
cmmg ip=internal parameter (0,1,2,3);it=internal type (stre,bend,tors)          
cmmg ipos=Position in der Hessian                                               
cmmg h1: Hessian rein (alt); h2: hessian raus (neu)                             
cmmg nvar: dimension of h1; ntot: dimension of h2                               
                                                                                
      write(*,*) 'transform a hessian with frozen variables'                    
      write(*,*) 'to a new hessian including those variables'                   
      write(*,*) 'which are marked with 2 in the input file'                    
      write(*,*) 'additional feature: transform a hessian with'                 
      write(*,*) 'redundancy'                                                   
      write(*,*) ' - which are marked by 3 in the input file - '                
      write(*,*) 'to a non-redundant hessian'                                   
      write(*,*) 'input : standard io (z-matrix)'                               
      write(*,*) '        old hessian from <hessian>'                           
      write(*,*) '        new hessian  on  <hessian>'                           
      write(*,*)                                                                
      write(*,*) 's. grimme, nov 1994'                                          
                                                                                
                                                                                
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      read (4,'(A)') a8                                                         
      call readl(80,a8,xx,nn)                                                   
      nvar=idint(xx(1))                                                         
      read (4,'(A)') A8                                                         
      do i=1,nvar                                                               
         read (4,'(6F12.3)')(h1(i,j),j=1,i)                                     
      enddo                                                                     
      close (4)                                                                 
      DO 220 I=1,NVAR                                                           
         DO 220 J=1,I                                                           
            H1(J,I)=H1(I,J)                                                     
  220 CONTINUE                                                                  
                                                                                
      k=0                                                                       
      l=0                                                                       
cmmg l is supposed to count the redundant parameters                            
cmmg which are marked by 3                                                      
   10 read(*,'(A)',end=100) A8                                                  
      call readl(80,a8,xx,nn)                                                   
      if(nn.gt.8) then                                                          
         itest=idint(xx(2))                                                     
         if(itest.gt.0) then                                                    
            k=k+1                                                               
            ip(k)=itest                                                         
            it(k)=1                                                             
         endif                                                                  
cmmg1                                                                           
         if(itest.eq.3) then                                                    
         l=l+1                                                                  
         endif                                                                  
cmmg2                                                                           
         itest=idint(xx(4))                                                     
         if(itest.gt.0) then                                                    
            k=k+1                                                               
            ip(k)=itest                                                         
            it(k)=2                                                             
         endif                                                                  
cmmg1                                                                           
         if(itest.eq.3) then                                                    
         l=l+1                                                                  
         endif                                                                  
cmmg2                                                                           
         itest=idint(xx(6))                                                     
         if(itest.gt.0) then                                                    
            k=k+1                                                               
            ip(k)=itest                                                         
            it(k)=3                                                             
         endif                                                                  
cmmg1                                                                           
         if(itest.eq.3) then                                                    
         l=l+1                                                                  
         endif                                                                  
cmmg2                                                                           
      endif                                                                     
      goto 10                                                                   
 100  continue                                                                  
      ntot=k-l                                                                  
                                                                                
      k=0                                                                       
      do i=1,ntot                                                               
         ipos(i)=0                                                              
         if(ip(i).eq.1) then                                                    
            k=k+1                                                               
            ipos(i)=k                                                           
         endif                                                                  
cmmg1                                                                           
         if(ip(i).eq.3) then                                                    
            k=k+2                                                               
            ipos(i)=k                                                           
         endif                                                                  
cmmg2                                                                           
      enddo                                                                     
      nvar=k                                                                    
                                                                                
      do i=1,ntot                                                               
           do j=1,i                                                             
              if(ipos(i).gt.0.and.ipos(j).gt.0)then                             
                 h2(i,j)=h1(ipos(i),ipos(j))                                    
              else                                                              
                 if(i.eq.j) h2(i,j)=xinit(it(i))                                
              endif                                                             
           enddo                                                                
      enddo                                                                     
                                                                                
                                                                                
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      write(4,*) 'hessian dim=',ntot                                            
      write(4,*) 'lowest eigenvalue=',eigval                                    
      do i=1,ntot                                                               
         write(4,'(6F12.3)')(h2(i,j),j=1,i)                                     
      enddo                                                                     
      close (4)                                                                 
                                                                                
      END                                                                       
                                                                                
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
