      SUBROUTINE JCARIN (XPARAM,STEP,PRECI,B,NCOL)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'param'                                                           
C     JACOBIAN dCARTESIAN/dINTERNAL, WORKED OUT BY FINITE DIFFERENCE.           
C  INPUT                                                                        
C     XPARAM(*) : INTERNAL COORDINATES                                          
C     STEP      : STEP SIZE FOR FINITE DIFFERENCE METHOD                        
C     PRECI     : .TRUE. IF 2-POINTS FINITE DIFFERENCES TO BE USED,             
C                 .FALSE. OTHERWISE.                                            
C  OUTPUT                                                                       
C     B(NVAR,NCOL) : JACOBIAN, STEP TIME TOO LARGE.                             
C                                                                               
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),                     
     1                      LOCDEP(MAXPAR)                                      
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),          
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,                
     2                NCLOSE,NOPEN,NDUMY,FRACT                                  
     3       /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, YPARAM(MAXPAR)                 
     4       /EULER / TVEC(3,3),ID                                              
     5       /UCELL / L1L,L2L,L3L,L1U,L2U,L3U                                   
     6       /GEOM  / GEO(3,NUMATM)                                             
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),                                    
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)                                          
      DIMENSION COORD(3,NUMATM),XPARAM(*),COORD1(3,NUMATM)            
      DIMENSION B(MAXPAR,MAXPAR)
      LOGICAL PRECI                                                             
C                                                                               
C     INTERNAL OF CENTRAL POINT                                                 
C     DO 10 IVAR=1,NVAR                                                         
C  10 GEO(LOC(2,IVAR),LOC(1,IVAR))=XPARAM(IVAR)                                 
C     IF(NDEP.NE.0) CALL SYMTRY                                                 
C                                                                               
C                                                                               
C        MOLECULAR SYSTEM                                                       
C        ----------------                                                       
         DO 30 IVAR=1,NVAR                                                      
C        STEP FORWARD                                                           
            GEO(LOC(2,IVAR),LOC(1,IVAR))=XPARAM(IVAR)+STEP                      
            IF(NDEP.NE.0) CALL SYMTRY                                           
            CALL GMETRY (GEO,COORD)                                             
            K=0                                                                 
            DO I=1,NATOMS                                                       
               IF(LABELS(I).LT.99) THEN                                         
                  K=K+1                                                         
                  DO L=1,3                                                      
                     COORD1(L,K)=COORD(L,I)                                     
                  ENDDO                                                         
               ENDIF                                                            
            ENDDO                                                               
C HP
            kk=0
            do 20 j=1,numat 
            do 20 jj=1,3
               kk=kk+1
   20       b(ivar,kk)=coord1(jj,j)

   30    GEO(LOC(2,IVAR),LOC(1,IVAR))=XPARAM(IVAR)                              
            DO 50 IVAR=1,NVAR                                                   
C        STEP BACKWARD                                                          
            GEO(LOC(2,IVAR),LOC(1,IVAR))=XPARAM(IVAR)-STEP                      
            IF(NDEP.NE.0) CALL SYMTRY                                           
            CALL GMETRY (GEO,COORD)                                             
            K=0                                                                 
            DO I=1,NATOMS                                                       
               IF(LABELS(I).LT.99) THEN                                         
                  K=K+1                                                         
                  DO L=1,3                                                      
                     COORD1(L,K)=COORD(L,I)                                     
                  ENDDO                                                         
               ENDIF                                                            
            ENDDO                                                               
C HP 
            kk=0
            do 40 j=1,numat 
            do 40 jj=1,3
               kk=kk+1
   40       b(ivar,kk)=b(ivar,kk)-coord1(jj,j)

   50    GEO(LOC(2,IVAR),LOC(1,IVAR))=XPARAM(IVAR)                              
                                                                                
      RETURN                                                                    
      END                                                                       
