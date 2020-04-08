                                                                                
      subroutine cekin(v,ekin)                                                  
      implicit real*8 (a-h,o-z)                                                 
      include 'param'                                                           
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)                          
      COMMON /MASS/ AMS(107)                                                    
      dimension v(*)                                                            
      m=0                                                                       
      ekin=0.0d0                                                                
      do i=1,natoms                                                             
         if(labels(i).lt.99) then                                               
            amass=ams(labels(i))*1.6605402d-27*0.50D0                           
            tx=amass*v(m+1)**2                                                  
            ty=amass*v(m+2)**2                                                  
            tz=amass*v(m+3)**2                                                  
            m=m+3                                                               
            ekin=ekin+dsqrt(tx*tx+ty*ty+tz*tz)                                  
         endif                                                                  
      enddo                                                                     
      ekin=ekin*2.29371049d+17                                                  
      return                                                                    
      end                                                                       
                                                                                
      subroutine deriv(mode,numat,x,f)                                          
      implicit real*8 (a-h,o-z)                                                 
      include 'param'                                                           
      dimension x(*), f(*), g(3,numatm), coord(3,numatm)                        
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)                          
      COMMON /MASS/ AMS(107)                                                    
      common /rk4   / escf                                                      
                                                                                
      m=0                                                                       
      do i=1,natoms                                                             
         if(labels(i).lt.99) then                                               
            m=m+1                                                               
            coord(1,i)=x(m)*1.0d+10                                             
            m=m+1                                                               
            coord(2,i)=x(m)*1.0d+10                                             
            m=m+1                                                               
            coord(3,i)=x(m)*1.0d+10                                             
         endif                                                                  
      enddo                                                                     
                                                                                
      CALL OUTI(1,MODE,COORD)                                                   
                                                                                
C calc energy and gradients                                                     
                                                                                
      ESCF=0.0D0                                            
C     ghost
      write(*,*) "I am subroutine DERIV. Calling egrad now!"
      CALL EGRAD(MODE,NUMAT,ESCF,G)                                             
                                                                                
      IF(DABS(ESCF).LT.1.0D-10) STOP 'EXTERNAL PROGRAM ERROR'                   
                                                                                
      IF(MODE.EQ.2) CALL perm(natoms,labels,coord,g)                            
                                                                                
      SUM=0.0D0                                                                 
      DO I=1,NUMAT                                                              
         DO J=1,3                                                               
            SUM=SUM+G(J,I)**2                                                   
         ENDDO                                                                  
      ENDDO                                                                     
                                                                                
      CGRAD=DSQRT(SUM)                                                          
                                                                                
      IF(DABS(CGRAD).LT.1.0D-3) STOP 'EXTERNAL PROGRAM ERROR (GRAD)'            
                                                                                
      i=0                                                                       
      m=0                                                                       
      do 30 k=1,natoms                                                          
         if(labels(k).gt.98) goto 30                                            
         i=i+1                                                                  
         amass=ams(labels(k))                                                   
         do j=1,3                                                               
            m=m+1                                                               
            force=g(j,i)*4.1840d13/6.0221367d23                                 
            xmass=amass*1.6605402d-27                                           
            f(m) =-force/xmass                                                  
         enddo                                                                  
   30  continue                                                                 
                                                                                
       return                                                                   
       end                                                                      
                                                                                
       subroutine irc(mode,numat,xparam)                                        
       implicit real*8 (a-h,o-z)                                                
      include 'param'                                                           
       dimension xparam(*), f(3*numatm), x(3*numatm), xm(3*numatm)              
       dimension g(3,numatm),v(3*numatm),coord(3,numatm), pm(3*numatm)          
       dimension xmem(maxpar)                                                   
       common /rk4   / escf                                                     
                                                                                
       COMMON /KEYWRD/ KEYWRD                                                   
       CHARACTER*241 KEYWRD                                                     
       COMMON /ELEMTS/ ELEMNT(107)                                              
       CHARACTER*2 ELEMNT                                                       
       COMMON /MASS/ AMS(107)                                                   
       COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                   
     1                 NA(NUMATM),NB(NUMATM),NC(NUMATM)                         
       COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)                  
       COMMON /GEOM  / GEO(3,NUMATM)                                            
                                                                                
C                                                                               
C SET UP COORDINATES FOR CURRENT CALCULATION                                    
C                                                                               
      maxloop=1000                                                              
      I=index(keywrd,'MAXLOOP=')                                                
      if(i.ne.0) maxloop=idint(reada(keywrd,i))                                 
                                                                                
      ttot=100.d-15                                                             
      I=index(keywrd,'TTOT=')                                                   
      if(i.ne.0) ttot=1.0D-15*reada(keywrd,i)                                   
                                                                                
      tstep=0.5d-15                                                             
      I=index(keywrd,'TINC=')                                                   
      if(i.ne.0) tstep=1.0D-15*reada(keywrd,i)                                  
                                                                                
      write(*,*)                                                                
      write(*,1)                                                                
  1   format(70('-'))                                                           
      write(*,'(30x,'' IRC PROGRAM '')')                                        
      write(*,1)                                                                
                                                                                
      nat3=3*numat                                                              
                                                                                
      do i=1,nat3                                                               
         v(i)=0.0d0                                                             
      enddo                                                                     
                                                                                
      open(unit=33,file='ircgeo')                                               
      open(unit=34,file='ircvelo')                                              
                                                                                
      if(index(keywrd,'RVELO').ne.0) then                                       
         write(*,*)                                                             
         write(*,*) 'RESTART FROM <velo>'                                       
         write(*,*)                                                             
         open(unit=35,file='velo')                                              
         read(35,'(3d15.6)')(x(i),i=1,nat3)                                     
         read(35,'(3d15.6)')(v(i),i=1,nat3)                                     
         close (35)                                                             
      else                                                                      
                                                                                
      DO 50 I=1,NVAR                                                            
         K=LOC(1,I)                                                             
         L=LOC(2,I)                                                             
   50 GEO(L,K)=XPARAM(I)                                                        
                                                                                
                                                                                
C      IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS        
                                                                                
      IF(NDEP.NE.0) CALL SYMTRY                                                 
                                                                                
C      NOW COMPUTE THE ATOMIC COORDINATES.                                      
                                                                                
      CALL GMETRY(GEO,COORD)                                                    
                                                                                
      WRITE(6,'(//10X,''INTERNAL COORDINATES'',/)')                             
      CALL GEOUT(1)                                                             
                                                                                
      WRITE(6,'(//10X,''CARTESIAN COORDINATES'',/)')                            
      WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''X'',                                
     .9X,''Y'',9X,''Z'',/)')                                                    
                                                                                
      M=0                                                                       
      DO 240 I=1,NATOMS                                                         
         WRITE(6,'(I6,8X,A2,4X,3F12.6)')                                        
     .  I,ELEMNT(LABELS(I)),(COORD(J,I),J=1,3)                                  
         IF(LABELS(I).LT.99) THEN                                               
            M=M+1                                                               
            X(M)=COORD(1,I)                                                     
            M=M+1                                                               
            X(M)=COORD(2,I)                                                     
            M=M+1                                                               
            X(M)=COORD(3,I)                                                     
         ENDIF                                                                  
  240 CONTINUE                                                                  
                                                                                
      endif                                                                     
                                                                                
      do i=1,nat3                                                               
         x(i)=x(i)*1.0d-10                                                      
         xmem(i)=x(i)                                                           
      enddo                                                                     
                                                                                
      icount=0                                                                  
                                                                                
C EVALUATE EPOT AT FIRST POINT                                                  
                                                                                
      call deriv(mode,numat,x,f)                                                
                                                                                
C EVALUATE EKIN AT FIRST POINT                                                  
                                                                                
      call cekin(v,ekin)                                                        
                                                                                
      efirst=escf+ekin                                                          
                                                                                
      write(*,*)                                                                
      write(*,'('' TIMESTEP (fs)   ='',f12.2)')tstep*1.0d15                     
      write(*,'('' TOTAL T  (fs)   ='',f12.2)')ttot*1.0d15                      
      write(*,'('' TOTAL E POT (au)='',f12.6)')escf                             
      write(*,'('' TOTAL E KIN (au)='',f12.6)')ekin                             
      write(*,'('' TOTAL E     (au)='',f12.6)')efirst                           
      write(*,*)                                                                
      write(*,22)                                                               
                                                                                
 20   format(i4,3x,F6.2,2x,3(2x,f12.6),2X,f12.4)                                
 21   format(   3x,F6.2,2x,3(5x,f12.6))                                         
 22   format(9x,'TIME (fs)',5X,'ETOT',10X,'EKIN',10X,'ERR (%)',5X,              
     .'RMS DISPLACEMENT')                                                       
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                  LOOP BEGINS HERE                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
                                                                                
 999  CONTINUE                                                                  
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                  implicit midpoint integration                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
                                                                                
       if(icount.lt.5) then                                                     
          tinc=tstep*0.5                                                        
       else                                                                     
          tinc=tstep                                                            
       endif                                                                    
                                                                                
       tsum=tsum+tinc                                                           
                                                                                
       h=tinc                                                                   
       h5=0.50d0*h                                                              
                                                                                
       do i=1,nat3                                                              
          xm(i)=x(i)+h5*v(i)                                                    
       enddo                                                                    
                                                                                
       call deriv(mode,numat,xm,f)                                              
                                                                                
       do i=1,nat3                                                              
          pm(i)=v(i)+h5*f(i)                                                    
       enddo                                                                    
                                                                                
       call cekin(v,ekin1)                                                      
                                                                                
       do i=1,nat3                                                              
          v(i)=v(i)+h*f(i)                                                      
          x(i)=x(i)+h*pm(i)                                                     
       enddo                                                                    
                                                                                
       call cekin(v,ekin2)                                                      
                                                                                
       ekin=(ekin1+ekin2)*0.50d0                                                
                                                                                
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C         end of midpoint integration                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
                                                                                
      icount=icount+1                                                           
                                                                                
      etot =escf+ekin                                                           
      err  =100.*(efirst-etot)/ekin                                             
                                                                                
      m=0                                                                       
      rms=0.0d0                                                                 
      do i=1,nat3                                                               
         dx=(xmem(m+1)-x(m+1))**2                                               
         dy=(xmem(m+2)-x(m+2))**2                                               
         dz=(xmem(m+3)-x(m+3))**2                                               
         rms=rms+dsqrt(dx+dy+dz)                                                
         m=m+3                                                                  
      enddo                                                                     
                                                                                
      write(*,20) icount,tsum*1.d15,escf,ekin,err,1.0d10*rms/nat3               
      write(7,21)        tsum*1.d15,escf,ekin                                   
                                                                                
C OUTPUT                                                                        
                                                                                
         open(unit=35,file='velo')                                              
         rewind 35                                                              
         write(35,'(3d15.6)')(1.0D10*x(i),i=1,nat3)                             
         write(35,'(3d15.6)')(v(i),i=1,nat3)                                    
         close (35)                                                             
         write(34,*)'IRC CYCLE',icount,tsum*1.d15                               
         write(34,'(3d15.6)')(1.0D10*x(i),i=1,nat3)                             
         write(34,'(3d15.6)')(v(i),i=1,nat3)                                    
                                                                                
         m=0                                                                    
         do i=1,natoms                                                          
            if(labels(i).lt.99) then                                            
            m=m+1                                                               
            coord(1,i)=x(m)*1.0d+10                                             
            m=m+1                                                               
            coord(2,i)=x(m)*1.0d+10                                             
            m=m+1                                                               
            coord(3,i)=x(m)*1.0d+10                                             
            endif                                                               
        enddo                                                                   
                                                                                
        call out(33,mode,coord)                                                 
                                                                                
C                                                                               
                                                                                
      if(tsum.le.ttot.and.icount.le.maxloop) goto 999                           
                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                  LOOP ENDS HERE                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
                                                                                
      return                                                                    
      end                                                                       
