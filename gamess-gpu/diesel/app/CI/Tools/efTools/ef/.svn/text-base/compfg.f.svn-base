                                                                                
      SUBROUTINE COMPFG(XPARAM,GRAD)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
      DIMENSION XPARAM(*),GRAD(*)                                               
      COMMON /KEYWRD/KEYWRD                                                     
      CHARACTER*241 KEYWRD                                                      
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)                   
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),          
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,                
     2                NCLOSE,NOPEN,NDUMY,FRACT                                  
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)                                          
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)         
      COMMON /ELEMTS/ ELEMNT(107)                                               
      COMMON /CYCLE / ICYC                                                      
      COMMON /ENERG / ELAST, CGRAD                                              
      COMMON /PRTFLG/ IPRT                                                      
      CHARACTER*2 ELEMNT                                                        
      COMMON /GEOM  / GEO(3,NUMATM)                                             
      DIMENSION NAME(3),NAM(3)
      CHARACTER*10 NAME                                                         
      CHARACTER*3  NAM, AOL                                                    
      DIMENSION GXYZ(3,NUMATM)                                                  
      DIMENSION GXYZL(3*NUMATM)                                                  
      DIMENSION GEOG(3,NUMATM)                                                  
      DIMENSION COORD(3,NUMATM)                                                 
      DIMENSION B(MAXPAR,MAXPAR) 
      DIMENSION AOL(MAXPAR) 
      DIMENSION U(MAXPAR*MAXPAR), G(MAXPAR*(MAXPAR+1)/2), E(MAXPAR)
      DATA NAME /'LENGTH','BEND','DIHEDRAL'/                                    
      DATA NAM  /'LEN','BEN','DIH'/                                    


ccccccccc intern
      logical intern
      intern = .false.
ccccccccc intern

                                                                                
      ICYC=ICYC+1                                                               
                                                                                
      write(*,*)                                                                
      write(*,1)                                                                
  1   format(70('-'))                                                           
      write(*,'(30x,'' CYCLE '',I3)') ICYC                                      
      write(*,1)                                                                
                                                                                
      CALL OPTOUT(XPARAM)                                                       

      IF(INDEX(KEYWRD,'DSCF').NE.0) MODE=0                                      
      IF(INDEX(KEYWRD,'TCSCF').NE.0) MODE=1   
      IF(INDEX(KEYWRD,'MOLPRO').NE.0) MODE=2                                    
      IF(INDEX(KEYWRD,'MOPAC').NE.0) MODE=3                                     
      IF(INDEX(KEYWRD,'ZINDO').NE.0) MODE=4                                     
      IF(INDEX(KEYWRD,'MP2').NE.0) MODE=5
      IF(INDEX(KEYWRD,'RIMP2').NE.0) MODE=6
      IF(INDEX(KEYWRD,'RIDFT').NE.0) MODE=7
      IF(INDEX(KEYWRD,'DIESEL').NE.0) MODE=8
	  
c	  if(MODE.EQ.8) write(*,*) 'switch to mode for diesel'
c      logical   intern = .false.
C       check, if gradient should be calculated in internal coordinates
      if(index(keywrd,'INTERN').ne.0) then
        write(*,*) 'calculating gradient in internal coordinates'
        intern = .true.
      endif
      
      
              
                                                                                     
C                                                                               
C SET UP COORDINATES FOR CURRENT CALCULATION                                    
C                                                                               
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
      DO 240 I=1,NATOMS                                                         
         WRITE(6,'(I6,8X,A2,4X,3F12.6)')                                        
     .  I,ELEMNT(LABELS(I)),(COORD(J,I),J=1,3)                                  
  240 CONTINUE                                                                  
                                                                                
                                                                                
C output the coordinates in the proper format                                   
                                                                                
      CALL OUT(1,MODE,COORD)                                                    

C make the B-matrix

      NCOL=3*NUMAT                                                              
      STEP=1.0D-3                                                               
      CALL JCARIN (XPARAM,STEP,.true.,B,NCOL)                               

C check fo linear dependencies

      k=0 
      do i=1,nvar
         IT=LOC(2,I)                                                             
         aol(i)=NAM(IT)
         do j=1,i   
            k=k+1
            g(k)=0.0d0
            do m=1,ncol
               g(k)=g(k)+b(i,m)*b(j,m)
            enddo
            g(k)=g(k)*0.5/step
         enddo
      enddo

      CALL HQRII(G,NVAR,NVAR,E,U)

      WRITE(*,'(/,''LOWEST EIGENVALUE OF THE B*B(T) MATRIX ='',E12.6)')
     .E(1)

      IF(ABS(E(1)).LT.1.0D-9) THEN
         WRITE(*,*)
         WRITE(*,*) 'LINEAR DEPENDENCIES IN INTERNAL COORDINATES'
         WRITE(*,*) 'DETECTED.'
         WRITE(*,'(5E12.6)')(E(I),I=1,NVAR)       
         call prmat(6,aol,u,nvar,nvar,'EIGENVECTORS OF B*B(T) MATRIX')
         STOP 'LINEAR DEPENDENCIES IN INTERNAL COORDINATES'
      ENDIF
                                                                                
C calc energy and gradients                                                     
                                                                                
      ESCF=0.0D0
      
cccccccccccccccc       
c       use different routines for gradient-evaluation in internal coordinates
      if (intern) goto   111
cccccccccccccccc      
      
c     write (*,*) "I am compfg, i will call egrad now!"
      CALL EGRAD(MODE,NUMAT,ESCF,GXYZ)                                          
C     if(iprt) then                                                             
C     WRITE(6,'(/10X,''TOTAL ENERGY (AU) = '',F15.7)') ESCF                     
C     WRITE(6,'(/10X,''ENERGY CHANGE     = '',F15.7/)') ESCF-ELAST              
C     endif                                                                     
      EDIFF=ESCF-ELAST
      ELAST=ESCF                                                                
                                                                                
      IF(DABS(ESCF).LT.1.0D-10) STOP 'EXTERNAL PROGRAM ERROR'                   
                                                                                
      IF(MODE.EQ.2) CALL perm(natoms,labels,coord,gxyz)                         
                                                                                
      SUM=0.0D0                                                                 
      K=0
      DO I=1,NUMAT                                                              
         DO J=1,3                                                               
            K=K+1
            SUM=SUM+GXYZ(J,I)**2                                                
            GXYZL(K)=GXYZ(J,I)
            GRAD(K)=0.0D0
         ENDDO                                                                  
      ENDDO                                                                     
                                                                                
      CGRAD=DSQRT(SUM)                                                          
                                                                                
      WRITE(6,'(/10X,''CARTESIAN GRADIENTS (KCAL/ANG)'',/)')                    
      WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''FX'',                               
     .9X,''FY'',9X,''FZ'',/)')                                                  
                                                                                
      K=0                                                                       
      DO 241 I=1,NATOMS                                                         
         IF(LABELS(I).LT.99) THEN                                               
         K=K+1                                                                  
         WRITE(6,'(I6,8X,A2,4X,3F10.4)')                                        
     .  I,ELEMNT(LABELS(I)),(GXYZ(J,K),J=1,3)                                   
        ENDIF                                                                   
  241 CONTINUE                                                                  
                                                                                

C transform cart grad -> internal grad                                          
  
      DO 110 K=1,NCOL                                                         
         DO 110 I=1,NVAR                                                      
  110 grad(I)=grad(I)+B(I,K)*GXYZL(K)                                               

      step=0.5/step                                                             
      do i=1,nvar                                                               
         grad(i)=grad(i)*step                                                   
      enddo                                                                     


c      evaluation of gradient in internal coordinates
  111 continue
      if (intern) then
        call egradint(MODE,xparam,ESCF,grad)
        write(*,*) 'approximating cart. grad. norm with int. grad. norm'
        gnfina = SQRT(DOT(GRAD,GRAD,NVAR)) 
        cgrad = gnfina
        EDIFF=ESCF-ELAST
        ELAST=ESCF                                                                
      endif
    
                                                                                
      WRITE(6,'(/5X,''INTERNAL GEOMETRIC VARIABLES AND GRADIENTS''/)')        
      WRITE(6,'(20X,'' ANG OR RAD      KCAL/ANG OR KCAL/RAD''/)')             
      DO I=1,NVAR                                                             
         IT=LOC(2,I)                                                             
         WRITE(6,'(5X,I3,3X,A10,2X,2F12.5)')I,NAME(IT),XPARAM(I),GRAD(I)         
      ENDDO                                                                   
C                                                                               
      GNFINA=SQRT(DOT(GRAD,GRAD,NVAR))                                       
C     WRITE(6,251) GNFINA                                               
C 251 FORMAT(/,' NORM OF INTERNAL GRADIENT = ',F9.2)                         

      WRITE(*,350) ESCF, EDIFF, GNFINA, CGRAD
  350 FORMAT(/,'# TOTAL E = ',F14.6,  
     .'   EDIFF = ',F10.6,'   GRAD I = ',F7.2, 
     .'   GRAD C = ',F7.2)
                                                                                
      escf=escf*627.5095410d0                                                   
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
