      SUBROUTINE HESSIAN(XPARAM)                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      integer dir,double
      DIMENSION XPARAM(*)                                                       
      COMMON /KEYWRD/KEYWRD                                                     
      CHARACTER*241 KEYWRD                                                      
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)                   
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GRAD(MAXPAR)                              
      DIMENSION HESS(MAXPAR,MAXPAR)                                             
      DIMENSION XINC(3),UC(MAXPAR,MAXPAR),EIGVAL(MAXPAR)                        
      DIMENSION HESSC(MAXPAR*(MAXPAR+1)/2)                                      
      logical hessinfo
                                                                                
      DATA XINC/0.010,0.020,0.020/                                              
                                                                                
      WRITE(6,'(/'' CALCULATING HESSIAN NUMERICALLY IN INT. C.''/)')            
                                                                                
c get the initialization data

      double=0
      i=index(keywrd,' DOUBLE=')
      if(i.ne.0) double=READA(KEYWRD,I)
      if(double.ne.0.and.double.ne.1)then
      write(*,*)
      write(*,*)'*******WARNING********'
      write(*,*)'RECONSIDER THE VALUE OF DOUBLE'
      write(*,*)'internal hessian cannot be calculated'
      write(*,*)'program stops'
      write(*,*)'*    END OF WARNING      * '
      write(*,*)
      goto 999
      endif

      dir=0
      i=index(keywrd,' DIR=')
      if(i.ne.0) dir=reada(keywrd,i)
      if(dir.gt.nvar) then
      write(*,*)
      write(*,*)'*******WARNING********'
      write(*,*)'RECONSIDER THE VALUE OF DIR'
      write(*,*)'internal hessian cannot be calculated'
      write(*,*)'program stops'
      write(*,*)'*    END OF WARNING      * '
      write(*,*)
      goto 999
      endif

c enough initialization data

c begin of "no restart" case

      if(dir.eq.0)then

      CALL COMPFG(XPARAM, GRAD)                  
      call iout(grad,nvar,0,0)
                                                                                
      kk=0
      DO 210 I=1,NVAR                                                           
            kk=kk+1
            L=LOC(2,I)                                                          
            WRITE(6,'(/'' STEP '',I3)') I                                       
            XPARAM(I)=XPARAM(I) + XINC(L)                                       
            CALL COMPFG(XPARAM, GNEXT1)         
            call iout(gnext1,nvar,i,kk)
            XPARAM(I)=XPARAM(I) - XINC(L)                                       
            DO 204 J=1,NVAR                                                     
  204       HESS(I,J)= (GNEXT1(J)-GRAD(J))/XINC(L)                              
        if(double.eq.1) then
            kk=kk+1
            XPARAM(I)=XPARAM(I) - XINC(L)                                       
            CALL COMPFG(XPARAM, GNEXT1)         
            call iout(gnext1,nvar,i,kk)
            XPARAM(I)=XPARAM(I) + XINC(L)                                       
            do 205 j=1,nvar
 205        hess(i,j)=0.5*(hess(i,j)+(-gnext1(j)+grad(j))/xinc(l))
        endif
  210 CONTINUE                                                                  
      endif

c end of "no restart" case

c in the case of a restart, file <hessinfo> will be read
c and the remaining elements for file <hessian> will
c be calculated

      if(dir.ne.0) then
        inquire(file='hessinfo',exist=hessinfo)
        if(.not.hessinfo) then
          write(*,*)'*******WARNING********'
          write(*,*)'FILE <hessinfo> CANNOT BE READ'
          write(*,*)'program stops'
          write(*,*)'*   END OF WARNING   *'
          goto 999
        endif
       open(unit=13,file='hessinfo',status='old')
       rewind 13
       write(*,*) ' '
       write(*,*)'read gradient from reference geometry'
       write(*,*) ' '
       read(13,*) i
       read(13,'(5f14.8)') (grad(j),j=1,nvar)
c      write(*,'(5f14.8)') (grad(j),j=1,nvar)
C kk=index for number of calculated
C gradients apart from reference geometry
       kk=0
       do 110 i=1,dir
         L=LOC(2,I)                                                          
         kk=kk+1
         read(13,*) j
         read(13,'(5f14.8)') (gnext1(j),j=1,nvar)
c        write(*,'(5f14.8)') (gnext1(j),j=1,nvar)
         write(*,*) ' '
         write(*,*)'read gradient from elongation geometry',i
         do 101 j=1,nvar
 101     hess(i,j)=(gnext1(j)-grad(j))/xinc(l)
         if(double.eq.1) then
           kk=kk+1
           read(13,*) j
           read(13,'(5f14.8)') (gnext1(j),j=1,nvar)
         write(*,'(5f14.8)') (gnext1(j),j=1,nvar)
           write(*,*)'read gradient from second elongation geometry',i
           do 102 j=1,nvar
 102       hess(i,j)=0.5*(hess(i,j)+(-gnext1(j)+grad(j))/xinc(l))
         endif
 110  continue
      write(*,*) ' '
      close(13)

      do 304 i=dir+1,nvar
         kk=kk+1
         write(*,*)' '
         write(*,*)'******* '
         write(*,*) 'START WITH ELONGATION OF DIRECTION.',i
         write(*,*)'******* '
         write(*,*)' '
            L=LOC(2,I)
            WRITE(6,'(/'' STEP '',I3)') I
            XPARAM(I)=XPARAM(I) + XINC(L)
            CALL COMPFG(XPARAM, GNEXT1)
            call iout(gnext1,nvar,i,kk)
            XPARAM(I)=XPARAM(I) - XINC(L)
            DO 200 J=1,NVAR
  200       HESS(I,J)= (GNEXT1(J)-GRAD(J))/XINC(L)
        if(double.eq.1) then
            kk=kk+1
            XPARAM(I)=XPARAM(I) - XINC(L)
            CALL COMPFG(XPARAM, GNEXT1)
            call iout(gnext1,nvar,i,kk)
            XPARAM(I)=XPARAM(I) + XINC(L)
            do 201 j=1,nvar
 201        hess(i,j)=0.5*(hess(i,j)+(-gnext1(j)+grad(j))/xinc(l))
        endif
  304 CONTINUE
      endif

C end of restart case

C     SYMMETRIZE HESSIAN                                                        
                                                                                
      DO 230 I=1,NVAR                                                           
         DO 220 J=1,I-1                                                         
            HESS(I,J)=(HESS(I,J)+HESS(J,I))*0.50D0                              
            HESS(J,I)=HESS(I,J)                                                 
  220 CONTINUE                                                                  
  230 CONTINUE                                                                  
                                                                                
      IJ=0                                                                      
      DO 290 I=1,NVAR                                                           
         DO 290 J=1,I                                                           
            IJ=IJ+1                                                             
            HESSC(IJ)=HESS(J,I)                                                 
  290 CONTINUE                                                                  
                                                                                
      CALL HQRII(HESSC,NVAR,NVAR,EIGVAL,UC)                                     
                                                                                
      open(unit=4,file='hessian')                                               
      rewind 4                                                                  
      write(4,*) 'hessian dim=',nvar                                            
      write(4,*) 'lowest eigenvalue=',eigval(1)                                 
      do i=1,nvar                                                               
         write(4,'(6F12.3)')(hess(i,j),j=1,i)                                   
      enddo                                                                     
      close (4)                                                                 
                                                                                
      call hessout(hess,0,0)                                                      
                                                                                
 999  continue

      return                                                                    
      end                                                                       
