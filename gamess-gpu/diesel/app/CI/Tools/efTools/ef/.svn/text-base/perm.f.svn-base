      subroutine perm(numat,lab,coord,gxyz)                                     
      implicit real*8 (a-h,o-z)                                                 
      include 'param'                                                           
                                                                                
      COMMON /EXTCO / COMP(3,NUMATM)                                            
                                                                                
      dimension coord(3,*), gxyz(3,*), lab(*)                                   
      dimension coord2(3,numatm), iperm(numatm)                                 
      dimension g(3,numatm)                                                     
                                                                                
                                                                                
C remove dummys                                                                 
                                                                                
                                                                                
      nat=0                                                                     
      do i=1,numat                                                              
         if(lab(i).lt.99) then                                                  
            nat=nat+1                                                           
            do l=1,3                                                            
               coord2(l,nat)=coord(l,i)/0.529177260d0                           
            enddo                                                               
         endif                                                                  
      enddo                                                                     
                                                                                
                                                                                
      do i=1,nat                                                                
         do l=1,3                                                               
            g(l,i)=gxyz(l,i)                                                    
         enddo                                                                  
      enddo                                                                     
                                                                                
                                                                                
c     write(*,*) 'comp'                                                         
c     do i=1,nat                                                                
c        write(*,*)(comp(l,i),l=1,3)                                            
c     enddo                                                                     
c     write(*,*) 'coord2'                                                       
c     do i=1,nat                                                                
c        write(*,*)(coord2(l,i),l=1,3)                                          
c     enddo                                                                     
                                                                                
                                                                                
C where is atom i in molpro arrangement                                         
                                                                                
      do i=1,nat                                                                
         iperm(i)=0                                                             
         do j=1,nat                                                             
            id=0                                                                
            do l=1,3                                                            
               sum=coord2(l,i)-comp(l,j)                                        
               if(dabs(sum).lt.1.0d-4) id=id+1                                  
            enddo                                                               
            if(id.eq.3) then                                                    
               iperm(i)=j                                                       
               goto 100                                                         
            endif                                                               
         enddo                                                                  
 100     continue                                                               
      enddo                                                                     
                                                                                
c     do i=1,nat                                                                
c        write(*,*)i, iperm(i)                                                  
c     enddo                                                                     
                                                                                
      do i=1,nat                                                                
         if(iperm(i).eq.0) stop 'error inside perm'                             
         do l=1,3                                                               
            gxyz(l,i)=g(l,iperm(i))                                             
         enddo                                                                  
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
                                                                                
                                                                                
