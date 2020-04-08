      subroutine transformmu(dipxl,dipyl,dipzl,dipxr,dipyr,dipzr,uc
     $     ,numat,nvarn,ddipQx,ddipQy,ddipQz,cstep,amass)

C     --------------------
C     transforms and derivates \mu as \mu(X) to d\mu/dQ
C     --------------------
C     . nvar=nvarn in hessout, wenn in cart coord gearbeitet wird,
C       sonst ist nvar=n(interne coords.)
C     . in hessc existiert nvar sowieso nicht.
C     . AMASS ist ein vektor, auf dem die massen stehen wie:
C       masse (atom1) masse (atom1) masse(atom 1)
C       masse (atom2) masse (atom2) masse(atom 2)
C       masse (atom3) masse (atom3) masse(atom 3)
C       etc.
C     . gridsize kommt als cstep und in {\AA} an!

      implicit  none
      integer*4 numat
      real*8 uc(3*numat,3*numat)
      real*8 ddipQx(7:3*numat)
      real*8 ddipQy(7:3*numat)
      real*8 ddipQz(7:3*numat)
      real*8 dipxl(*),dipyl(*),dipzl(*)
      real*8 dipxr(*),dipyr(*),dipzr(*)
      real*8 cstep
      real*8 atoBohr
c     real*8 aIntensity
      real*8 amass(3*numat)
      integer*4 mode
      integer*4 i
      integer*4 nvarn
c     dimension aIntensity(7:3*numat)

      AtoBohr = .52917726d0

C     laengster nach innen!
C     und da bilden wir ein vor-symmetrisiertes dipolmoment
C     mit intrinsischer ableitung...
      do mode = 7,3*numat
         ddipQx(mode)=0.d0
         ddipQy(mode)=0.d0
         ddipQz(mode)=0.d0
         i=0
         do i=1,3*numat
            ddipQx(mode)=ddipQx(mode)+(dipxl(i)-dipxr(i))*uc(i,mode)
     $           /(2.d0*cstep*sqrt(amass(i))/AtoBohr)
            ddipQy(mode)=ddipQy(mode)+(dipyl(i)-dipyr(i))*uc(i,mode)
     $           /(2.d0*cstep*sqrt(amass(i))/AtoBohr)
            ddipQz(mode)=ddipQz(mode)+(dipzl(i)-dipzr(i))*uc(i,mode)
     $           /(2.d0*cstep*sqrt(amass(i))/AtoBohr)
         enddo 
      enddo
      return
      end
