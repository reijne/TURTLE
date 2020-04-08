      subroutine calcint(ddipQx,ddipQy,ddipQz,numat,cstep,aIntensity)
C     --------------------
C     ghost
C     . calculates IR intensities from derivated cmponents 
c       of dipole momenta from subroutine transformmu
C     . throws away 1st 6 eigvals.
C     . prints out
C     --------------------
      implicit real*8 (a-h,o-z)

      dimension ddipQx(7:3*numat)
      dimension ddipQy(7:3*numat)
      dimension ddipQz(7:3*numat)
      dimension aIntensity(7:3*numat)

      iEVal = 3*numat
C     data taken from gaussian94, files utilnz.f and phyfil.f:
      pi=4.d0*datan(1.0d0)
      Avo=6.02213670d0
      c=2.99792458d0
      toe=4.803242D0
      umrechnung = toe * dsqrt((avo*avo*10.0d0*pi)/(3.0d0*c*c))
c     auf Debye/(Angstrom*amu)

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
     $     ,'mol.', //)
      return
      end
