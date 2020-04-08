      subroutine calcint(numat,nvarn,dipxl,dipyl,dipzl,dipxr,dipyr,dipzr,uc,DipRefx,DipRefy,DipRefz)
      parameter (maxvar=12)

c     (dipxl,dipyl,dipzl,dipxr,dipyr,dipzr,uc
c     $     ,numat,nvarn,ddipQx,ddipQy,ddipQz,cstep,DipRefx,DipRefy
c     $     ,DipRefz)
C     --------------------
C     transforms and derivates \mu as \mu(X) to d\mu/dQ
C     --------------------
C     . nvar=nvarn in hessout, wenn in cart coord gearbeitet wird,
C       sonst ist nvar=n(interne coords.)
C     . in hessc existiert nvar sowieso nicht.


      implicit  none
      
      common /mass/ ams(107)
      common /elemts/ elemnt(107)
      common /geokst/ natoms,labels(numatm)

      real*8 ams
      character*2 elemnt

      real*8 uc(maxvar,maxvar)
      real*8 ddipQx(7:maxvar)
      real*8 ddipQy(7:maxvar)
      real*8 ddipQz(7:maxvar)
      real*8 dipxl(maxvar),dipyl(maxvar),dipzl(maxvar)
      real*8 dipxr(maxvar),dipyr(maxvar),dipzr(maxvar)
      real*8 cstep
      real*8 atobohr
      real*8 aIntensity
      real*8 a,b,c,d,e,f,diprefx,diprefy,diprefz
      integer*4 mode
      integer*4 i,j,k,m
      integer*4 maxvar
      integer*4 numat
      integer*4 nvarn
      dimension aIntensity(7:maxvar)
c      real*8 amass(4)
c      data amass /14.d0,3*1.d0/

      Atobohr = .529177249
     

      open(file='dipoles',form='formatted',unit=70)
      read(70,*) numat,nvarn
      read(70,*) dipxl,dipyl,dipzl,dipxr,dipyr,dipzr
      read(70,*) uc
      read(70,*) DipRefx,DipRefy,DipRefz
      close(70)


C     DRAN DENKEN!!!
      cstep = 0.0025d0
C     ===============
      write(*,*) 'Distortion length was:',cstep,'Angstroms'





      i = 0
      do j = 1,numat
         do k = 1,3
            i=i+1
      a=  DipRefx - dipxl(i)
      b= -DipRefx + dipxr(i)
      c=  DipRefy - dipyl(i)
      d= -DipRefy + dipyr(i)
      e=  DipRefz - dipzl(i)
      f= -DipRefz + dipzr(i)
c     write(*,*)'Atom',j,'x li:',a,'re:',b
c     write(*,*)'y li:',c,'re',d
c     write(*,*)'z li:',e,'re',f
c     write(*,*) ' '
      write(*,*) 'ddipXx', (-(a+b)/(2.d0*cstep/atobohr))*4.803242d0
      write(*,*) 'ddipXy', (-(c+d)/(2.d0*cstep/atobohr))*4.803242d0
      write(*,*) 'ddipXz', (-(e+f)/(2.d0*cstep/atobohr))*4.803242d0
      write(*,*) ' '
      enddo
      enddo




C     laengster nach innen!
C     und da bilden wir ein vor-symmetrisiertes dipolmoment
C     mit intrinsischer ableitung...
      do mode = 7,3*numat
         ddipQx(mode)=0.d0
         ddipQy(mode)=0.d0
         ddipQz(mode)=0.d0
         i=0
         do m=1,numat
            do j=1,3
               i=i+1
            ddipQx(mode)=ddipQx(mode)+(dipxl(i)-dipxr(i))*uc(i,mode)/(2
     $           .d0*cstep*sqrt(amass(m))/Atobohr)
            ddipQy(mode)=ddipQy(mode)+(dipyl(i)-dipyr(i))*uc(i,mode)/(2
     $           .d0*cstep*sqrt(amass(m))/Atobohr)
            ddipQz(mode)=ddipQz(mode)+(dipzl(i)-dipzr(i))*uc(i,mode)/(2
     $           .d0*cstep*sqrt(amass(m))/Atobohr)
         enddo 
         enddo
         write(*,*)'ddipQx von',mode,':',ddipQx(mode)*4.803242d0
         write(*,*)'ddipQy von',mode,':',ddipQy(mode)*4.803242d0
         write(*,*)'ddipQz von',mode,':',ddipQz(mode)*4.803242d0
      enddo


      call calcint(ddipQx,ddipQy,ddipQz,numat,cstep,aIntensity)


      stop
      end















C     ---------------------------------------------------------------------------
C     ---------------------------------------------------------------------------
      subroutine calcint(ddipQx,ddipQy,ddipQz,numat,cstep,aIntensity)
C     --------------------
C     ghost
C     . calculate IR intensities from derivated cmponents 
c       of dipole momenta from subroutine transformmu
C     . throw away 1st 6 eigvals.
C     . print out
C     --------------------
      implicit real*8 (a-h,o-z)
      
      dimension ddipQx(7:3*numat)
      dimension ddipQy(7:3*numat)
      dimension ddipQz(7:3*numat)
      dimension aIntensity(7:3*numat)

C     data taken from gaussian94, files utilnz.f and phyfil.f
C     ------------------
c     pi=4.d0*datan(1.0d0)
c     Avo=6.02213670d23
c     c=2.99792458d10
c     toe=4.803242D-10
c     FIR1 = toe * 10**10
C     ------------------
      pi=4.d0*datan(1.0d0)
      Avo=6.02213670d0
      c=2.99792458d0
      toe=4.803242D0

      umrechnung = toe * dsqrt((avo*avo*10.0d0*pi)/(3.0d0*c*c))

c     auf Debye/(Angstrom*amu)
      iEVal = 3*numat
      
      write (*,100)
      do i = 7, iEVal
         aIntensity(i) = Umrechnung*Umrechnung*(ddipQx(i)*ddipQx(i)
     $        +ddipQy(i)*ddipQy(i)+ddipQz(i)*ddipQz(i))
         write (*,200) i-6,aIntensity(i)
      enddo
      write(*,300)
 100  format (// 'Numerical Infra-Red Intensities:')
 200  format (' Mode',i3,':',f14.6)
 300  format ('Units are km','/'
     $     ,'mol.', /,'Valid for non-linear molecules only' //)
      return
      end

      

                                                                                
      BLOCK DATA                                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      COMMON /MASS/ AMS(107)                                                    
      COMMON /ELEMTS/ ELEMNT(107)                                               
      CHARACTER*2 ELEMNT                                                        
                                                                                
      DATA  AMS /  1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,               
     .10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,               
     .20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,               
     .30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,               
     .40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,               
     .54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,               
     .65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,               
     .79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,               
     .91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,               
     .102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,               
     .118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,               
     .132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,               
     .183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,               
     .196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0,               
     .18*0.000D0,   1.0079D0,  5*0.000D0/                                       
                                                                                
      DATA ELEMNT/' H','He',                                                    
     1 'Li','Be',' B',' C',' N',' O',' F','Ne',                                 
     2 'Na','Mg','Al','Si',' P',' S','Cl','Ar',                                 
     3 ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',                  
     4 'Zn','Ga','Ge','As','Se','Br','Kr',                                      
     5 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',                  
     6 'Cd','In','Sn','Sb','Te',' I','Xe',                                      
     7 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',             
     8 'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt',             
     9 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',                                 
     1 'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','XX',        
     2 'Fm','Md','Cb','++',' +','--',' -','Tv'/                                 
      END                                                                       
                                                                                
