_IF(charmmtest)
c
c These routines are replaced by code in the 
c charmm interface.
c
c Using the keyword charmmtest
c will get a version of gamess containing all
c the internal charmm interface code but linked to 
c run standalone (of testing interest only)

      subroutine chmdat(ztag10,czan1,cnew1,nat, nel,expo, wght)

      implicit none

INCLUDE(common/sizes)

c!!!!!!!!!!!
INCLUDE(common/blur)
c!!!!!!!!!!!

cINCLUDE(common/infoa)
cINCLUDE(common/infob)
INCLUDE(common/phycon)
cINCLUDE(common/runlab)

      character*10 ztag10(*)
      REAL czan1(*), cnew1(3,*)
      integer nat, nel
      integer i, j
      REAL expo(maxat)
      REAL wght(maxat)
c
c  Definition of the quantum region
c
      nat = 1

      ztag10(1) = 'he'
      cnew1(1,nat) = 0.0d0
      cnew1(2,nat) = 0.0d0
      cnew1(3,nat) = 0.0d0
      czan1(nat) =  2.0d0

      nel = 2
c
c Classical region (point charges)
c
      if ( .false. ) then
      ztag10(4) = 'bq'
      cnew1(1,4) = 3.0d0
      cnew1(2,4) = 0.0d0
      cnew1(3,4) = 0.0d0
      czan1(4) = -0.6d0

      ztag10(5) = 'bq'
      cnew1(1,5) = 3.0d0
      cnew1(2,5) = 0.0d0
      cnew1(3,5) = 1.0d0
      czan1(5) = 0.3d0

      ztag10(6) = 'bq'
      cnew1(1,6) = 3.0d0
      cnew1(2,6) = 0.5d0
      cnew1(3,6) = -0.6d0
      czan1(6) = 0.3d0
      nat = nat + 3
      endif
c
c Classical region (gaussian blur)
c
      if(.true.)then

         oblur = .true.

         do i=1,nat
            expo(i) = -1.0d0;
         enddo

c
c add 1 blurred charges
c
         nat = nat + 1

         ztag10(nat) = 'bq'
         cnew1(1,nat) = 0.0d0
         cnew1(2,nat) = 0.0d0
         cnew1(3,nat) = 1.0d0

         czan1(nat) =  0.0d0
         expo(nat) =   10.0d0
         wght(nat) =   1.0d0

      endif

      end

      program gamess_top
INCLUDE(common/chmgms)
      qinigm = .true.
      call gamess
      qinigm = .false.
      call gamess
      end
      subroutine bits

INCLUDE(common/phycon)

      character*10 ztag10(10)
      REAL cnew1(5,3), czan1(5), expo(5), wght(5)

      ztag10(2) = 'h'
      cnew1(2,1) = 0.0d0
      cnew1(2,2) = 0.0d0
      cnew1(2,3) = 0.01d0
      czan1(2) = 1.0d0

         nat = nat + 1
         ztag10(nat) = 'bq'
         cnew1(nat,1) = 0.0d0
         cnew1(nat,2) = 0.0d0
         cnew1(nat,3) = 30.0d0

         czan1(nat) =  0.0d0
         expo(nat) =   10.0d0
         wght(nat) =   1.0d0

      end

      subroutine nbndgcm
      call caserr ('nbndgcm')
      end

c
c Pass gamess-uk results back to charmm
c
      subroutine GMS2CHM(GTOT,DX,DY,DZ,NATOM)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/funct)

      REAL gtot, DX(*), DY(*), DZ(*)
      integer natom

      integer i
      REAL fac
c
c These are taken from charmm
c
      REAL BOHRR
      PARAMETER (BOHRR = 0.529177249D0 )

      REAL TOKCAL
      PARAMETER (TOKCAL = 627.5095D0 )

      fac = TOKCAL / BOHRR

      write(6,*)'func:',enrgy
      do i = 1,natom
         write(6,*)'grad:',egrad(3*(i-1)+1),egrad(3*(i-1)+2),
     &        egrad(3*(i-1)+3)
         dx(i) = dx(i) + fac*egrad(3*(i-1)+1)
         dy(i) = dy(i) + fac*egrad(3*(i-1)+2)
         dz(i) = dz(i) + fac*egrad(3*(i-1)+3)
      enddo

      fac = TOKCAL

      gtot = fac*enrgy

      end
_ENDIF
