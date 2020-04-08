      SUBROUTINE HESSOUT(hess,modehess,numat)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      COMMON /KEYWRD/KEYWRD
      CHARACTER*241 KEYWRD
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)
      DIMENSION HESS(MAXPAR,MAXPAR)
      DIMENSION UC(MAXPAR,MAXPAR),EIGVAL(MAXPAR)
      DIMENSION HESSC(MAXPAR*(MAXPAR+1)/2)
      DIMENSION NAME(3),AOL(MAXPAR)
      CHARACTER*3 NAME,AOL
      DATA NAME /'LEN','BEN','DIH'/
      
C     write(6,*)'modehess=',modehess
      pi = 4.0d0*atan(1.d0)
      speedoflight=2.99792458d10
      planck = 6.62607550D-27

      if(modehess.eq.1) nvar=3*numat
      IJ=0
      DO 280 I=1,NVAR
         DO 290 J=1,I
            IJ=IJ+1
            HESSC(IJ)=HESS(J,I)
 290     CONTINUE
 280  CONTINUE
      
      
c     call prmat(6,aol,hessc,nvar,0,'FINAL HESSIAN')    
      CALL HQRII(HESSC,NVAR,NVAR,EIGVAL,UC)
      

c     ghost
c     an dieser stelle sind die eigenwerte in einheiten
c     von kcal/({\AA}**2)



c     UNITS
      if(modehess.eq.1)then
c     kcal to kJ
         xkcaltokJ=4.1840d0
c     \AA to m
         aatom=1.0d-10
         f=xkcaltokJ/(aatom*aatom)
c     first check if at all we have neg ev's
         nNegflag = 0
         do 294 i=1,nvar
            if(eigval(i).lt.0) then 
               nNegCount=nNegCount+1
            endif
 294     continue


c     now print them out!
c         if (nNegCount.ne.0) then
c            write(*,*) 'No. of negative EVs:',nNegCount
c         endif

         do 295 i=1,nvar
            if (eigval(i).lt.0) then
               eigval(i)=f*eigval(i)
               eigval(i)=-( 1000.d0/(2.0d0*pi) *sqrt( -eigval(i) ) /
     $              speedoflight)
c     write(*,*) eigval(i), 'NEGATIVE!'
            else
               eigval(i)=f*eigval(i)
               eigval(i)=( 1000.d0/(2.0d0*pi) *sqrt( eigval(i) ) /
     $              speedoflight)
            endif
 295     continue



Cmmg    Berechnung der Nullpunktsschwingungsenergie
C       Werte die 'Null'sein sollen bleiben unberuecksichtigt
         
C     ghost:
c     quatsch. wir schmeissen einfach die 1. sechs weg,
c     wir sind ja noch in der modehess-abfrage!
         schw=0.0d0
         do 296 i=1,nvar
            if (eigval(i).gt.0.0d0) then
               schw=schw+.5d0*eigval(i)
            endif
 296     continue
C     Ausgabe der Nullpunktsschwingungsenergie in den 
C     richtigen Einheiten
         wntokcal =  2.85914392D-03
         wntokjmol =  1.19626582D-02
         
         write(*,*)' '
         write(*,*)'Nullpunktsschwingungsenergie:'
         write(*,*) wntokcal*schw,'kcal/mol ',wntokjmol*schw,'kJ/mol'
         write(*,*)' '
      endif
      
      if(modehess.ne.1) then
         do i=1,nvar
            aol(i)=name(loc(2,i))
         enddo
      endif
      if(modehess.eq.1) then
         do i=1,nvar
            aol(i)='CCO'
c     write(6,*) aol(i)
         enddo
      endif
      
      write(*,*) ' '
      write(*,*) 'NO. of negative EVs:',nNegCount
      
      call preig(6,eigval,nvar)                     
      call prmat(6,aol,uc,nvar,nvar,'FINAL HESSIAN EIGENVECTORS')
      
      return
      end
      
