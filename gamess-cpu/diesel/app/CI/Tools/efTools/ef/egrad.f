      SUBROUTINE EGRAD(MODE,NUMAT,ESCF,GRAD)                                    
C     this routine reads from other codes energies,
C     gradient components etc. and memorizes them.
C     ghost
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /EXTCO / COMP(3,NUMATM)                                            
      CHARACTER*241 KEYWRD                                                      
      DIMENSION GRAD(3,NUMATM),XX(10),ID(10)                                    
c     dimension ddipx(3*numat),ddipy(3*numat),ddipz(3*numat)
      CHARACTER*128 A                                                           
      
      do i=1,numat                                                              
         do j=1,3                                                               
            grad(j,i)=0.0d0                                                     
         enddo                                                                  
      enddo                                                                     
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
C mode = 0 dscf                                                                 
C        1 tc dscf                                                              
C        2 molpro                                                               
C        3 mopac                                                                
C        4 zindo                                                                
C        5 mp2
C        6 rimp2
C        7 ridft
C        8 diesel
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
                                                                                
C au/bohr -> kcal/angst                                                         
                                                                                
      f=627.509541d0/0.529177260d0                                              
      f2=0.529177260d0                                                          

C----------------- external program ---------------                             
                                                                                
      if(mode.eq.1) then                                                        
         call system('tcscf > dscf.out')                                      
         call system('tcgrad > grad.out')                                      
      endif                                                                     
                                                                                
      if(mode.eq.0) then                                                        
         call system('dscf > dscf.out')                                         
         call system('grad > grad.out')                                         
      endif                                                                     
                                                                                
      if(mode.eq.2) then                                                        
         call system('molpro -omp.out -W. -I. -d. mp.inp')                          
      endif                                                                     
                                                                                
      if(mode.eq.3) then                                                        
         call system('runmopac')                                                  
      endif                                                                     
                                                                                
      if(mode.eq.4) then                                                        
         call system('zindo1 > zindo.out')                                      
      endif                                                                     

      if(mode.eq.5) then                                                        
         call system('dscf > dscf.out')                                      
         call system('mpgrad > grad.out')                                      
      endif                                                                     

      if(mode.eq.6) then                                                        
         call system('dscf > dscf.out')                                      
         call system('rimp2 > grad.out')                                      
      endif                                                                     

      if(mode.eq.7) then                                                        
         call system('ridft > dscf.out')                                      
         call system('rdgrad > grad.out')                                      
      endif
	  
      if(mode.eq.8) then 
	  	write(*,*) 'starting diesel'                                                       
         call system('executeNG < numgrad.in > numgrad.out')                             
      endif                                                                     
                                                                                
C---------------- turbomole 2.1 code  --------------                            
                                                                                
      if(mode.le.1.or.mode.eq.7) then                                                        
                                                                                
      open(unit=2,file='dscf.out')                                              
      open(unit=3,file='grad.out')                                              
                                                                                
   10 read(2,'(A)',end=100) A                                                   
      if(index(A,' total energy = ').ne.0 .or.
     .   index(A,'|  total energy ').ne.0) then                                  
         call readl(128,a,xx,n)                                                 
         escf=xx(1)                                                             
         goto 100                                                               
      endif                                                                     
      goto 10                                                                   
  100 continue                                                                  
      close (2)                                                                 
                                                                                
   20 read(3,'(A)',end=200) A                                                   
      if(index(A,' ATOM').ne.0) then                                            
         call readl(128,a,xx,n)                                                 
         do i=1,n                                                               
            id(i)=idint(xx(i))                                                  
         enddo                                                                  
C x                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(1,id(j))=f*xx(j)                                               
         enddo                                                                  
C y                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(2,id(j))=f*xx(j)                                               
         enddo                                                                  
C z                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(3,id(j))=f*xx(j)                                               
         enddo                                                                  
      endif                                                                     
      goto 20                                                                   
  200 continue                                                                  
      close (3)                                                                 
      return                                                                    
      endif                                                                     
                                                                                
C---------------- end of turbomole 2.1  code ----------                         
                                                                                
C---------------- molpro94 code  --------------                                 
                                                                                
      if(mode.eq.2) then                                                        
                                                                                
         nroot=1                                                                
         idum=index(keywrd,'IROOT=')                                            
         if(idum.gt.0) nroot=reada(keywrd,idum)                                 
         if(nroot.eq.0) nroot=1                                                 
         nr=0                                                                   
                                                                                
      open(unit=3,file='mp.out')                                                
  5   read(3,'(A)',end=50) A                                                    
      if(index(a,' NR  ATOM  CHARGE   ').ne.0) then                             
         read(3,'(A)') A                                                        
         do i=1,numat                                                           
            read(3,'(A)',end=50) A                                              
            call readl(128,A,XX,NN)                                             
            comp(1,i)=xx(3)                                                     
            comp(2,i)=xx(4)                                                     
            comp(3,i)=xx(5)                                                     
         enddo                                                                  
      endif                                                                     
                                                                                
      if(index(a,'calculation of the energy grad').ne.0) then                   
         read(3,'(A)') A                                                        
         read(3,'(A)') A                                                        
         read(3,'(A)') A                                                        
         do i=1,numat                                                           
            read(3,*)idum,grad(1,i),grad(2,i),grad(3,i)                         
            do j=1,3                                                            
               grad(j,i)=f*grad(j,i)                                            
            enddo                                                               
         enddo                                                                  
      endif                                                                     
      if(index(a,'MC STATE').ne.0.and.index(a,'ENERGY').ne.0) then              
         nr=nr+1                                                                
         call readl(128,a,xx,nn)                                              
         if(nr.eq.nroot) escf=xx(nn)                                            
      endif                                                                     
      goto 5                                                                    
 50   close (3)                                                                 
      return                                                                    
      endif                                                                     
                                                                                
C---------------- end of molpro94 code  --------------                          

C---------------- begin of mopac  -------------                                 

      if(mode.eq.3) then                                                        
	 open(unit=3,file='grad.out')                                              
         read(3,*) escf
         escf=escf/627.509541d0
         do j=1,numat
            read(3,*) (grad(i,j),i=1,3)
         enddo
         close (3)
         return
      endif

C---------------- end of mopac  -------------                                 
                                                                                
C---------------- begin of zindo  -------------                                 
                                                                                
      if(mode.eq.4) then                                                        
                                                                                
	    open(unit=3,file='zgrad.out')                                              
C x                                                                             
	 read(3,*) n                                                                   
	   do j=1,n                                                                    
	    read(3,*) (grad(i,j),i=1,3)                                                
           enddo                                                                
	   do j=1,n                                                                    
	    do i=1,3                                                                   
cmmg	    write(*,*) grad(i,j)                                                   
	    grad(i,j)=f*grad(i,j)                                                      
cmmg	    write(*,*) grad(i,j)                                                   
	    enddo                                                                      
           enddo                                                                
	   read(3,*) escf                                                              
      close (3)                                                                 
      return                                                                    
      endif                                                                     
C---------------- end of zindo  -------------                                   
                                                                                
C---------------- turbomole 2.1 code  --------------                            
C     ghost MP2 or RIMP2

      if(mode.eq.5.or.mode.eq.6) then                                              
                                                                                
      open(unit=3,file='grad.out')                                              
                                                                                
   11 read(3,'(A)',end=101) A                                                   
      if(index(A,'*      total        ').ne.0) then
         call readl(128,a,xx,n)                                                 
         escf=xx(1)                                                             
         goto 101                                                               
      endif                                                                     
      goto 11                                                                   
  101 continue                                                                  
                                                                                
   21 read(3,'(A)',end=201) A                                                   
      if(index(A,' ATOM').ne.0) then                                            
         call readl(128,a,xx,n)                                                 
         do i=1,n                                                               
            id(i)=idint(xx(i))                                                  
         enddo                                                                  
C x                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(1,id(j))=f*xx(j)                                               
         enddo                                                                  
C y                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(2,id(j))=f*xx(j)                                               
         enddo                                                                  
C z                                                                             
         read(3,'(A)') A                                                        
         call readl(128,a,xx,n)                                                 
         do j=1,n                                                               
            grad(3,id(j))=f*xx(j)                                               
         enddo                                                                  
      endif                                                                     
      goto 21                                                                   
  201 continue                                                                  
      close (3)                                                                 
      return                                                                    
      endif                                                                     
                                                                                
                                                                                
C---------------- diesel code  --------------                            
C     jan 

      if(mode.eq.8) then                                              
                                                                                
      open(unit=3,file='energy.diesel')
	  read(3,*) escf                                              
      close(3)                                                                   
c      CHARACTER*80 buffer
      open(unit=3,file='gradient.diesel')
	  read(3,*) 
	  read(3,*)
	  do i=1,NUMAT
	    read(3,*)
	  enddo
	  do i=1,NUMAT
	    read(3,*) xx(1),xx(2),xx(3)
		do j=1,3
		  grad(j,i)=f*xx(j)
		enddo
	  enddo
      close(3)                                                                 


      return                                                                    
      endif                                                                     
                                                                                
C---------------- error section -------------                                   
                                                                                
C     ghost 
C     write(*,*)'MODE =',MODE                                                   
      stop 'UNKNOWN EXTERNAL PROGRAM'                                           
      end                                                                       
                                                                                
