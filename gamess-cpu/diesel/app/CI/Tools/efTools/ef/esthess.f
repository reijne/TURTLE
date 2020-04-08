      SUBROUTINE ESTHESS(NVAR,XPARAM,C,HESS)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      DIMENSION XPARAM(*)                                                       
      COMMON /KEYWRD/KEYWRD                                                     
      CHARACTER*241 KEYWRD                                                      
      COMMON /GEOVAR/ NDUM, LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)                  
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1NAA(NUMATM),NBB(NUMATM),NCC(NUMATM)                                       
      DIMENSION HESS(MAXPAR,MAXPAR), C(3,NUMATM)                                
      DIMENSION RAD(85)                                                         
      COMMON /ELEMTS/ EL(107)                                                   
      CHARACTER*2 EL                                                            
                                                                                
      Data Rad/                                                                 
     $  0.643d0,0.643d0,2.457d0,1.909d0,1.587d0,1.436d0,1.209d0,                
     $  1.096d0,1.020d0,0.945d0,2.986d0,2.646d0,2.400d0,2.192d0,                
     $  2.060d0,1.890d0,1.795d0,1.701d0,3.836d0,3.288d0,2.721d0,                
     $  2.494d0,2.305d0,2.230d0,2.211d0,2.211d0,2.192d0,2.173d0,                
     $  2.211d0,2.362d0,2.381d0,2.305d0,2.268d0,2.192d0,2.154d0,                
     $  2.116d0,4.082d0,3.609d0,3.061d0,2.740d0,2.532d0,2.457d0,                
     $  2.400d0,2.362d0,2.362d0,2.419d0,2.532d0,2.797d0,2.721d0,                
     $  2.665d0,2.646d0,2.570d0,2.513d0,2.476d0,4.441d0,3.742d0,                
     $  3.194d0,3.118d0,3.118d0,3.099d0,3.080d0,3.061d0,3.496d0,                
     $  3.042d0,3.005d0,3.005d0,2.986d0,2.967d0,2.948d0,2.948d0,                
     $  2.948d0,2.721d0,2.532d0,2.457d0,2.419d0,2.381d0,2.400d0,                
     $  2.457d0,2.532d0,2.816d0,2.797d0,2.778d0,2.759d0,2.759d0,                
     $  2.740d0/                                                                
                                                                                
      WRITE(6,'(/''ESTIMATING DIAGONAL HESSIAN (JPC 96, 1992, 9768)'')')        
      WRITE(6,'(/''FORCE CONSTANT (KCAL/ANG^2 OR KCAL/RAD^2)''/)')              
                                                                                
      DO I=1,NVAR                                                               
         DO J=1,NVAR                                                            
            HESS(I,J)=0.0D0                                                     
         ENDDO                                                                  
      ENDDO                                                                     
                                                                                
      xmin=1000000.0d0                                                          
                                                                                
      DO 210 I=1,NVAR                                                           
            L2=LOC(2,I)                                                         
            L1=LOC(1,I)                                                         
            ia=0                                                                
            ib=0                                                                
            ic=0                                                                
            id=0                                                                
            xh=0.0d0                                                            
            if(l2.eq.1) then                                                    
               ia=l1                                                            
               ib=naa(l1)                                                       
               ita=labels(ia)                                                   
               itb=labels(ib)                                                   
               rab=dsqrt((c(1,ia)-c(1,ib))**2+                                  
     .                   (c(2,ia)-c(2,ib))**2+                                  
     .                   (c(3,ia)-c(3,ib))**2)/0.529177260d0                    
               if(ita.eq.99.or.itb.eq.99) then                                  
                  xh=6000.0d0/rab                                               
               else                                                             
                  rcab=rad(ita)+rad(itb)                                        
                  xh=0.361*exp(-1.944*(rab-rcab))                               
                  xh=2.0d0*xh*627.51/0.529177260d0**2                           
               endif                                                            
               write(*,'(F12.1,9X,4(2X,A2))') xh,el(ita),el(itb)                
            endif                                                               
            if(l2.eq.2) then                                                    
               ib=l1                                                            
               ia=naa(l1)                                                       
               ic=nbb(l1)                                                       
               ita=labels(ia)                                                   
               itb=labels(ib)                                                   
               itc=labels(ic)                                                   
               if(ita.eq.99.or.itb.eq.99.or.itc.eq.99) then                     
                  xh=1000.0d0                                                   
               else                                                             
                  rabc=rad(ita)+rad(itb)                                        
                  racc=rad(ita)+rad(itc)                                        
                  rab=dsqrt((c(1,ia)-c(1,ib))**2+                               
     .                      (c(2,ia)-c(2,ib))**2+                               
     .                      (c(3,ia)-c(3,ib))**2)/0.529177260d0                 
                  rac=dsqrt((c(1,ia)-c(1,ic))**2+                               
     .                      (c(2,ia)-c(2,ic))**2+                               
     .                      (c(3,ia)-c(3,ic))**2)/0.529177260d0                 
                  xh=0.089+0.11/(rabc*racc)**(-0.42)                            
     .               *exp(-0.44*(rab+rac-rabc-racc))                            
                  xh=xh*627.51                                                  
               endif                                                            
               write(*,'(F12.1,9X,4(2X,A2))') xh,el(ita),el(itb),el(itc)        
            endif                                                               
            if(l2.eq.3) then                                                    
               ia=l1                                                            
               ic=naa(l1)                                                       
               id=nbb(l1)                                                       
               ib=ncc(l1)                                                       
               ita=labels(ia)                                                   
               itb=labels(ib)                                                   
               itc=labels(ic)                                                   
               itd=labels(id)                                                   
               xh=500.0d0                                                       
               write(*,'(F12.1,9X,4(2X,A2))')xh,el(ita),el(itb),                
     .                                          el(itc),el(itd)                 
            endif                                                               
            hess(i,i)=xh                                                        
            if(xh.lt.xmin) xmin=xh                                              
  210 CONTINUE                                                                  
                                                                                
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      write(4,*) 'hessian dim=',nvar                                            
      write(4,*) 'lowest eigenvalue=',xmin                                      
      do i=1,nvar                                                               
         write(4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 
                                                                                
      return                                                                    
      end                                                                       
