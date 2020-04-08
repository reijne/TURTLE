      SUBROUTINE EGRADINT(MODE,xparam,ESCF,GRAD)                                    
C     this routine writes info-file to perform 
C     gradient evaluation in internal coordinates,
C     starts the interface programs
C     and reads back the gradient-components
cccccccccccc 
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

      dimension xinc(3)
      data xinc/0.010,0.020,0.020/

      k=0
      do i=1,nvar
         k=k+1
         grad(i)=0.0d0
         xparam(i)=xparam(i)
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
      auToKcal =  627.509541d0                                                         

C----------------- external program ---------------                             
      if (mode.ne.8) then
        write(*,*) 'ERROR: internal gradient evaluation is only'
        write(*,*) ' possible with diesel'
      endif
      
      open(14,file='gradinfo')
c      write(*,*) 'starte mit dem rausschreiben' 
c      write(14,*) 'base' 
c        CALL OPTOUTint(XPARAM)
      
      do 100 i=1,nvar
        l = loc(2,i)
        xparam(i) = xparam(i) + xinc(l)
c        write(14,*)
c        write(14,*) 'var' , i,'p ',xinc(l)
        write(14,'(A,I8,A,F9.5)') 'var' , i,'+ ',xinc(l)
        CALL OPTOUTint(XPARAM)
        write(14,*)
        xparam(i) = xparam(i) - xinc(l)
c  weiterer Durchlauf fuer doppelte Genauigkeit    
c  eventuell abfrage einbauen    
        xparam(i) = xparam(i) - xinc(l)
c        write(14,*)
c        write(14,*) 'var' , i,'m ',xinc(l)
        write(14,'(A,I8,A,F9.5)') 'var' , i,'- ',-xinc(l)
c        write(14,'(A,I2,A,F9.5)') 'var' , 20,'- ',xinc(l)
        CALL OPTOUTint(XPARAM)
        write(14,*)
        xparam(i) = xparam(i) + xinc(l)
        
 100  enddo
      close(14)                                                                          
      
c       diesel starten
      write(*,*) 'starting diesel'
         call system('prepareDir  >> numgrad.out')                             
      
c        alles wieder einlesen
      escf = 0
      open(unit=3,file='energy.diesel')
	  read(3,*) escf                                              
      close(3)                                                                   
c        nicht vergessen den Gradienten mit f zu multiplizieren, weil alle weiteren
c        Rechnungen in (kcal/mol)/A sind, grad aber in a.u./bohr
      open(unit=3,file='gradient.diesel')
      do 110 i=1,nvar
        read(3,*) xx
        grad(i) = auToKcal*xx
 110  enddo
      close(3)                                                                 
      
      return
      end
