c
c  Modified handling of weights directive (original by Fritz Daalmans,
c  1994) modified P.Sherwood 1997, Principal changes:
c
c  1) input may occur either in the geometry block, or using weights 
c     directive placed either before or after the basis section, always
c     provided in the input atom ordering.
c
c  2) modifications and additions to the mass table are possible from 
c     the input stream
c
c  **these are reset for each component of a multi-step calculation**
c
c  3) masses up to z=86 are tabulated
c
c  4) masses averaged over natural abundances are available
c    
c  5) tables also contain abundant elements and so defaults
c     can replace separate tables used by individual modules.
c  
c**********************************************************************
c wnucld(nz,isowt)
c
c   nz     - atomic number
c   isowt - nint(mass required)   or 
c            0 (most abundant)    or
c           -1 (isotopic average)
c
c Atomic isotopes lookup table function.
c Author: Frits Daalmans, 1994
c modifed PS 1997 to take integer atomic number, reorganised data to
c to allow users to augment the table from the input deck
c**********************************************************************

      function wnucld(nz,isowt)
c
      implicit none

INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
      REAL wnucld
      integer nz,isowt

      integer iso,itest, indx
      REAL weight
      character*72 zerror
c
      if(isowt .eq. 0)then
c
c return the most abundant mass for the element
c
         iso = abundant(nz)

         indx=0
         itest=-1
         do while (itest .ne. 999)
	    indx = indx + 1
            itest = isotpz(indx)
            if(itest .eq. nz .and. nint(wisotp(indx)) .eq. iso)then
               wnucld = wisotp(indx) 
               return
            endif
         enddo

         write (zerror,'(1x,a2,i3,a18,i5,a20)')
     *        'z=',nz,', isotope approx. ',iso,
     *        ' not found in table.'

         call caserr(zerror)

      else if(isowt .eq. -1) then
c
c  return the average weight for the element
c
         wnucld = wave(nz)
         return

      else
c
c Now look up the isotopic weight that is closest
c
         weight = dfloat(isowt)
         indx=1
 10      continue
         if (isotpz(indx).eq.nz) then
            if (dabs(wisotp(indx)-weight).lt.0.5d0) goto 20
c     Trick to identify isotopes that are known, but for which the mass
c     is unknown... In the wisotp() table these have negative weights.
            if (dabs(wisotp(indx)+weight).lt.0.5d0) then
               call caserr(' isotope was found, but mass is unknown')
            end if
         end if
         indx=indx+1
         if (isotpz(indx).ne.999) goto 10
c     Give an informative error message
         write (zerror,'(1x,a6,i3,a18,i5,a20)')
     *        'atom z=',nz,', isotope approx. ',isowt,
     *        ' not found in table.'
         call caserr(zerror)
 20      continue
         wnucld = wisotp(indx)
      endif

      return
      end
c***********************************************************************
c
c Average masses (used to be contents of atwt before most abundant
c isotope masses were substituted
c
c See mass_init for most data initialisation
c
c***********************************************************************
      block data mass_data
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
      data nmassv/0/
c
c     data wave/1.00797d0,4.0026d0,6.939d0,9.0122d0,10.811d0,
c    *     12.01115d0,14.0067d0,15.9994d0,18.9984d0,20.183d0,
c    *     22.9898d0,24.312d0,26.9815d0,28.086d0,30.9738d0,
c    *     32.064d0,35.453d0,39.948d0,39.102d0,40.08d0,44.956d0,
c    *     47.90d0,50.942d0,51.996d0,54.9380d0,55.85d0,58.9332d0,
c    *     58.71d0,63.54d0,65.37d0,69.72d0,72.59d0,74.9216d0,
c    *     78.96d0,79.909d0,83.80d0,85.47d0,87.62d0,88.905d0,
c    *     91.22d0,92.906d0,95.94d0,99.0d0,101.07d0,102.905d0,
c    *     106.4d0,107.87d0,112.40d0,114.82d0,118.69d0,121.75d0,
c    *     127.60d0,126.9044d0,131.30d0/
c
c modified 26 Oct 1997 by C. Bayse
c
      data wave/1.00797d0,4.0026d0,6.939d0,9.0122d0,10.811d0,
     *     12.01115d0,14.0067d0,15.9994d0,18.9984d0,20.183d0,
     *     22.9898d0,24.312d0,26.9815d0,28.086d0,30.9738d0,
     *     32.064d0,35.453d0,39.948d0,39.102d0,40.08d0,44.956d0,
     *     47.90d0,50.942d0,51.996d0,54.9380d0,55.85d0,58.9332d0,
     *     58.71d0,63.54d0,65.37d0,69.72d0,72.59d0,74.9216d0,
     *     78.96d0,79.909d0,83.80d0,85.47d0,87.62d0,88.905d0,
     *     91.22d0,92.906d0,95.94d0,99.0d0,101.07d0,102.905d0,
     *     106.4d0,107.87d0,112.40d0,114.82d0,118.69d0,121.75d0,
     *     127.60d0,126.9044d0,131.30d0,132.9054d0,137.33d0,
     *     138.9055d0,140.12d0,140.9077d0,144.24d0,145.d0,
     *     150.4d0,151.96d0,157.25d0,158.9254d0,162.50d0,
     *     164.9304d0,167.26d0,168.9342d0,173.04d0,174.967d0,
     *     178.49d0,180.9479d0,183.85d0,186.207d0,190.2d0,
     *     192.22d0,195.09d0,196.9665d0,200.59d0,204.37d0,
     *     207.2d0,208.9804d0,209.d0,210.d0,222.d0/
c
      end

c***********************************************************************
c
c mass_newiso(nz,w) :  define a new isotopic mass
c
c if nint(w) distinguishes the new mass from previously stored masses
c of atomic number z, a new mass is added, otherwise the existing entry is
c overwritten.
c
c***********************************************************************

      subroutine mass_newiso(nz,w)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/iofile)
      integer nz
      REAL w

      integer indx, itest, iso
      logical onew
   
      iso = nint(w)
      onew = .true.

      do indx = 1,niso
         itest = isotpz(indx)
         if(itest .eq. nz .and. nint(wisotp(indx)) .eq. iso)then
            write(iwr,*)'Redefining mass for isotope: z=',
     &           nz,' iw=',iso,' w=',w
            wisotp(indx) = w
            onew = .false.
         endif
      enddo
      
      if(onew)then
         write(iwr,*)'Defining new isotope: z=',nz,'w=',w
         call mass_addiso(nz,w)
      endif

      end
c***********************************************************************
c
c  mass_init: initialise mass tables
c
c  imassv(ivec, iat) - set of vectors of nuclide specifications
c
c  initialises to a single table of 0s signifying the most abundant
c  nuclei.  This routine is run before input processing starts, so the
c  number of atoms is not known.
c
c isotopic masses (on the scale 12-C = 12.0d0) and most abundant isotopes 
c are stored, taken from:
c
c Z=1,18:
c   table 11-28 (table of the isotopes) in 
c   CRC handbook of Chem. & Phys., 79th edition.
c
c Z=19,54
c   IUPAC Quantities, Units and Symbols in Physical Chemistry
c   2nd edition (1993)  (72-Ge has been corrected)
c
c Z=55,86
c   Modern Physics, F.J. Blatt
c
c  Commented entries denoted (ams   x.xxx) were the default masses for
c  previous versions. They are retained for cross-checking if required.
c
c***********************************************************************
      subroutine mass_init
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
      integer iat

      nmassv=1
c
c initialise default mass tanle
c
      do iat = 1,maxat
        imassv(1, iat) = 0
      enddo     

      niso = 0
c
c store atomic masses
c
c
c H
c   (ams 1.0078252d0)
      call mass_addiso(1, 1.007825d0)
      call mass_addiso(1, 2.0140d0)
      call mass_addiso(1, 3.01605d0)
      abundant(1) = 1

c He  (ams 4.0026036d0)
      call mass_addiso(2,3.01603d0) 
      call mass_addiso(2,4.00260d0) 
      call mass_addiso(2,5.01222d0) 
      call mass_addiso(2,6.018886d0) 
      call mass_addiso(2,7.02803d0) 
      call mass_addiso(2,8.0375d0) 
      call mass_addiso(2,9.0438d0)
      abundant(2) = 4

c Li  (ams 7.016005d0)
      call mass_addiso(3,5.01254d0) 
      call mass_addiso(3,6.015121d0) 
      call mass_addiso(3,7.016003d0) 
      call mass_addiso(3,8.022485d0) 
      call mass_addiso(3,9.026789d0)
      call mass_addiso(3,11.0439d0)
      abundant(3) = 7

c Be (ams 9.012186d0)
      call mass_addiso(4,6.019725d0) 
      call mass_addiso(4,7.016928d0) 
      call mass_addiso(4,8.005305d0) 
      call mass_addiso(4,9.012182d0)
      call mass_addiso(4,10.013534d0)
      call mass_addiso(4,11.021658d0)
      call mass_addiso(4,12.02692d0)
      call mass_addiso(4,13.0375d0)
      abundant(4) = 9

c B (ams 11.0093051d0)
      call mass_addiso(5,7.0299d0) 
      call mass_addiso(5,8.024605d0) 
      call mass_addiso(5,9.013328d0)
      call mass_addiso(5,10.012937d0)
      call mass_addiso(5,11.00931d0)
      call mass_addiso(5,12.0143d0)
      call mass_addiso(5,13.0178d0)
      abundant(5) = 11

c C (ams 12.000000d0)
      call mass_addiso(6,8.03767d0) 
      call mass_addiso(6,9.031039d0)
      call mass_addiso(6,10.01686d0)
      call mass_addiso(6,11.01143d0)
      call mass_addiso(6,12.000000d0)
      call mass_addiso(6,13.003355d0)
      call mass_addiso(6,14.003241d0)
      call mass_addiso(6,15.010599d0)
      call mass_addiso(6,16.014701d0)
      call mass_addiso(6,17.02257d0)
      call mass_addiso(6,18.0267d0)
      call mass_addiso(6,19.0350d0)
      call mass_addiso(6,20.0398d0) 
      abundant(6) = 12

c N (ams 14.0030744d0)
      call mass_addiso(7,12.018613d0)
      call mass_addiso(7,13.005738d0)
      call mass_addiso(7,14.003074d0)
      call mass_addiso(7,15.000108d0)
      call mass_addiso(7,16.006099d0)
      call mass_addiso(7,17.008450d0)
      call mass_addiso(7,18.014081d0)
      call mass_addiso(7,19.017040d0)
      call mass_addiso(7,20.0237d0)
      call mass_addiso(7,21.0266d0)
      call mass_addiso(7,22.0343d0)
      call mass_addiso(7,-23.0d0)   
      abundant(7) = 14

c O (ams 15.9949149d0)
      call mass_addiso(8,12.03442d0)
      call mass_addiso(8,13.02810d0)
      call mass_addiso(8,14.008595d0)
      call mass_addiso(8,15.003065d0)
      call mass_addiso(8,15.994915d0)
      call mass_addiso(8,16.999131d0)
      call mass_addiso(8,17.999160d0)
      call mass_addiso(8,19.003577d0)
      call mass_addiso(8,20.004075d0)
      call mass_addiso(8,21.008730d0)
      call mass_addiso(8,22.0101d0)
      call mass_addiso(8,23.0157d0)
      call mass_addiso(8,24.0202d0)
      abundant(8) = 16

c F (ams 18.998405d0)
      call mass_addiso(9,15.0180d0)
      call mass_addiso(9,16.01147d0)
      call mass_addiso(9,17.002095d0)
      call mass_addiso(9,18.000937d0)
      call mass_addiso(9,18.998403d0)
      call mass_addiso(9,19.999981d0)
      call mass_addiso(9,20.999948d0)
      call mass_addiso(9,22.003030d0)
      call mass_addiso(9,23.003600d0)
      call mass_addiso(9,24.0082d0)
      call mass_addiso(9,25.0122d0)
      call mass_addiso(9,26.0198d0)  
      call mass_addiso(9,27.0275d0) 
      abundant(9) = 19

c Ne  (ams 19.9924404d0)
      call mass_addiso(10,16.02575d0)
      call mass_addiso(10,17.017690d0)
      call mass_addiso(10,18.005710d0)
      call mass_addiso(10,19.001879d0)
      call mass_addiso(10,19.992435d0)
      call mass_addiso(10,20.993843d0)
      call mass_addiso(10,21.991383d0)
      call mass_addiso(10,22.994465d0)
      call mass_addiso(10,23.993613d0)
      call mass_addiso(10,24.997690d0)
      call mass_addiso(10,26.00048d0)
      call mass_addiso(10,27.0075d0)
      call mass_addiso(10,28.0115d0)
      abundant(10) = 20

c Na (ams 22.989773d0)
      call mass_addiso(11,19.013879d0)
      call mass_addiso(11,20.007344d0)
      call mass_addiso(11,20.997650d0)
      call mass_addiso(11,21.994434d0)
      call mass_addiso(11,22.989767d0)
      call mass_addiso(11,23.990961d0)  
      call mass_addiso(11,24.989953d0)
      call mass_addiso(11,25.992586d0)
      call mass_addiso(11,26.993940d0)
      call mass_addiso(11,27.978780d0)
      call mass_addiso(11,29.002830d0)
      call mass_addiso(11,30.008800d0) 
      call mass_addiso(11,31.012680d0)
      call mass_addiso(11,32.0177d0)
      call mass_addiso(11,33.0230d0)
      call mass_addiso(11,34.0286d0)
      abundant(11) = 23

c Mg (ams 23.985045d0)
      call mass_addiso(12,20.018864d0)
      call mass_addiso(12,21.011716d0)
      call mass_addiso(12,21.999574d0)
      call mass_addiso(12,22.994124d0)
      call mass_addiso(12,23.985042d0)
      call mass_addiso(12,24.985837d0)
      call mass_addiso(12,25.982593d0)
      call mass_addiso(12,26.984341d0)
      call mass_addiso(12,27.983876d0)
      call mass_addiso(12,28.98848 d0)
      call mass_addiso(12,29.990230d0)
      call mass_addiso(12,30.995930d0)
      call mass_addiso(12,31.9992d0)
      call mass_addiso(12,33.0055d0)
      call mass_addiso(12,34.0091d0)
      abundant(12) = 24

c Al (ams 26.981535d0)
      call mass_addiso(13,22.079370d0)
      call mass_addiso(13,23.007265d0)
      call mass_addiso(13,23.999941d0)
      call mass_addiso(13,24.990429d0)
      call mass_addiso(13,25.986892d0)
      call mass_addiso(13,26.981539d0)
      call mass_addiso(13,27.981910d0)
      call mass_addiso(13,28.980446d0)
      call mass_addiso(13,29.982940d0)
      call mass_addiso(13,30.983800d0)
      call mass_addiso(13,31.9880d0)
      call mass_addiso(13,32.9905d0)
      call mass_addiso(13,33.9965d0)
      call mass_addiso(13,34.9997d0)  
      call mass_addiso(13,36.0054d0)
      abundant(13) = 27

c Si (ams 27.976927d0)
      call mass_addiso(14,-22.0d0)
      call mass_addiso(14,24.011546d0)
      call mass_addiso(14,25.004109d0)
      call mass_addiso(14,25.992330d0)
      call mass_addiso(14,26.986704d0)
      call mass_addiso(14,27.976927d0)
      call mass_addiso(14,28.976495d0)
      call mass_addiso(14,29.973770d0)
      call mass_addiso(14,30.975362d0)
      call mass_addiso(14,31.974148d0)
      call mass_addiso(14,32.977920d0)
      call mass_addiso(14,33.97636d0)
      call mass_addiso(14,34.9844d0)
      call mass_addiso(14,35.9863d0)
      call mass_addiso(14,36.9925d0)
      call mass_addiso(14,37.9942d0)
      call mass_addiso(14,-39.0d0)
      abundant(14) = 28

c P (ams 30.973763d0)
      call mass_addiso(15,26.012080d0)
      call mass_addiso(15,26.99919 d0)
      call mass_addiso(15,27.992313d0)
      call mass_addiso(15,28.981803d0)
      call mass_addiso(15,29.978307d0)
      call mass_addiso(15,30.973762d0)
      call mass_addiso(15,31.973907d0)
      call mass_addiso(15,32.971725d0)
      call mass_addiso(15,33.973636d0)
      call mass_addiso(15,34.973232d0)
      call mass_addiso(15,35.977570d0)
      call mass_addiso(15,36.9795d0)
      call mass_addiso(15,37.9843d0)
      call mass_addiso(15,38.9866d0)
      call mass_addiso(15,39.9925d0)
      call mass_addiso(15,-41.0d0)
      call mass_addiso(15,-42.0d0)
      abundant(15) = 31

c S (ams 31.972074d0)
      call mass_addiso(16,28.996610d0)
      call mass_addiso(16,29.984903d0)
      call mass_addiso(16,30.979554d0)
      call mass_addiso(16,31.972070d0)
      call mass_addiso(16,32.971456d0)
      call mass_addiso(16,33.967866d0)
      call mass_addiso(16,34.969031d0)
      call mass_addiso(16,35.967080d0)
      call mass_addiso(16,36.971125d0)
      call mass_addiso(16,37.971162d0)
      call mass_addiso(16,38.975310d0)
      call mass_addiso(16,39.97582d0)
      call mass_addiso(16,40.9808d0)
      call mass_addiso(16,41.9824d0)
      call mass_addiso(16,-43.0d0)
      call mass_addiso(16,-44.0d0)
      abundant(16) = 32

c Cl (ams 34.968854d0)
      call mass_addiso(17,30.992410d0)
      call mass_addiso(17,31.985690d0)
      call mass_addiso(17,32.977451d0)
      call mass_addiso(17,33.973763d0)
      call mass_addiso(17,34.968852d0)
      call mass_addiso(17,35.968306d0)
      call mass_addiso(17,36.965903d0)
      call mass_addiso(17,37.968010d0)
      call mass_addiso(17,38.968005d0)
      call mass_addiso(17,39.970440d0)
      call mass_addiso(17,40.970590d0)
      call mass_addiso(17,41.9735d0)
      call mass_addiso(17,42.97516d0)
      call mass_addiso(17,43.9785d0)
      call mass_addiso(17,-45.0d0)
      abundant(17) = 35

c Ar (ams 39.962384d0)
      call mass_addiso(18,-31.0d0)
      call mass_addiso(18,31.997660d0)
      call mass_addiso(18,32.989930d0)
      call mass_addiso(18,33.980269d0)
      call mass_addiso(18,34.975256d0)
      call mass_addiso(18,35.967545d0)
      call mass_addiso(18,36.966776d0)
      call mass_addiso(18,37.962732d0)
      call mass_addiso(18,38.964314d0)
      call mass_addiso(18,39.962384d0)
      call mass_addiso(18,40.964501d0)
      call mass_addiso(18,41.963050d0)
      call mass_addiso(18,42.965670d0)
      call mass_addiso(18,43.96365 d0)
      call mass_addiso(18,44.968090d0)
      call mass_addiso(18,45.968090d0)
      call mass_addiso(18,46.9722d0)  
      abundant(18) = 40

c K (ams 38.963714d0)
      call mass_addiso(19,38.963707d0)
      call mass_addiso(19,40.961825d0)
      abundant(19) = 39

c Ca (ams 39.962589d0)
      call mass_addiso(20,39.962591d0)
      call mass_addiso(20,43.955481d0)
      abundant(20) = 40

c Sc (ams 44.955919d0)
      call mass_addiso(21,44.955910d0)
      abundant(21) = 45

c Ti (ams 47.947948d0)
      call mass_addiso(22,45.952629d0) 
      call mass_addiso(22,46.951764d0) 
      call mass_addiso(22,47.947947d0)
      call mass_addiso(22,48.947871d0) 
      call mass_addiso(22,49.944792d0) 
      abundant(22) = 48

c V (ams 50.943978d0)
      call mass_addiso(23,50.943962d0)
      abundant(23) = 51

c Cr (ams 51.940514d0)
      call mass_addiso(24,49.946046d0) 
      call mass_addiso(24,51.940510d0)
      call mass_addiso(24,52.940651d0) 
      call mass_addiso(24,53.938883d0) 
      abundant(24) = 52

c Mn  (ams 54.938054d0)
      call mass_addiso(25,54.938047d0)
      abundant(25) = 55

c Fe (ams 55.93493d0)
      call mass_addiso(26,53.939613d0) 
      call mass_addiso(26,55.934939d0)
      call mass_addiso(26,56.935396d0) 
      abundant(26) = 56

c Co (ams 58.933189d0)
      call mass_addiso(27,58.933198d0)
      abundant(27) = 59

c Ni (ams 57.93534d0)
      call mass_addiso(28,57.935346d0)
      call mass_addiso(28,59.930788d0) 
      call mass_addiso(28,60.931058d0) 
      call mass_addiso(28,61.928346d0) 
      abundant(28) = 58

c Cu (ams 62.92959d0)
      call mass_addiso(29,62.929599d0)
      call mass_addiso(29,64.927793d0) 
      abundant(29) = 63

c Zn (ams 63.929145d0)
      call mass_addiso(30,63.929145d0)
      call mass_addiso(30,65.926035d0) 
      call mass_addiso(30,66.927129d0) 
      call mass_addiso(30,67.924846d0) 
      abundant(30) = 64

c Ga (ams 68.9257d0)
      call mass_addiso(31,68.925580d0)
      call mass_addiso(31,70.924701d0) 
      abundant(31) = 69

c Ge (ams 73.921d0)
      call mass_addiso(32,69.924250d0) 
c note this value is corrected from the original source
      call mass_addiso(32,71.922079d0) 
      call mass_addiso(32,72.923463d0) 
      call mass_addiso(32,73.921177d0)
      call mass_addiso(32,75.921402d0) 
      abundant(32) = 74

c As (ams 74.921d0)
      call mass_addiso(33,74.921594d0) 
      abundant(33) = 75

c Se (ams 79.916d0)
      call mass_addiso(34,75.919212d0) 
      call mass_addiso(34,76.919913d0) 
      call mass_addiso(34,77.917308d0) 
      call mass_addiso(34,79.916520d0)
      call mass_addiso(34,81.916698d0) 
      abundant(34) = 80

c Br (ams 78.914d0) ***
      call mass_addiso(35,78.918336d0)
      call mass_addiso(35,80.916289d0) 
      abundant(35) = 79

c Kr (ams 83.911d0)
      call mass_addiso(36,79.916380d0) 
      call mass_addiso(36,81.913482d0) 
      call mass_addiso(36,82.914135d0) 
      call mass_addiso(36,83.911507d0)
      call mass_addiso(36,85.910616d0) 
      abundant(36) = 84

c Rb (ams 84.911d0)
      call mass_addiso(37,84.911794d0) 
      call mass_addiso(37,86.909187d0) 
      abundant(37) = 85

c Sr (ams 87.905d0)
      call mass_addiso(38,85.909267d0) 
      call mass_addiso(38,86.908884d0) 
      call mass_addiso(38,87.905619d0)
      abundant(38) = 88

c Y (ams 88.905d0)
      call mass_addiso(39,88.905849d0)
      abundant(39) = 89

c Zr (ams 89.904d0)
      call mass_addiso(40,89.904703d0) 
      call mass_addiso(40,90.905644d0)
      call mass_addiso(40,91.905039d0) 
      call mass_addiso(40,93.906315d0) 
      call mass_addiso(40,95.908275d0) 
      abundant(40) = 90

c Nb (ams 92.906d0)
      call mass_addiso(41,92.906377d0)
      abundant(41) = 93

c Mo (ams 97.906d0)
      call mass_addiso(42,91.906809d0) 
      call mass_addiso(42,93.905085d0) 
      call mass_addiso(42,94.905841d0) 
      call mass_addiso(42,95.904679d0) 
      call mass_addiso(42,96.906021d0) 
      call mass_addiso(42,97.905407d0)
      call mass_addiso(42,99.907477d0) 
      abundant(42) = 98

c Tc (ams 98.906d0) *****
      call mass_addiso(43,97.907215d0)
      abundant(43) = 98

c Ru (ams 101.903d0)
      call mass_addiso(44,95.907599d0) 
      call mass_addiso(44,97.905287d0) 
      call mass_addiso(44,98.905939d0) 
      call mass_addiso(44,99.904219d0) 
      call mass_addiso(44,100.905582d0) 
      call mass_addiso(44,101.904349d0)	 
      call mass_addiso(44,103.905424d0) 
      abundant(44) = 102

c Rh (102.904d0)
      call mass_addiso(45,102.905500d0)
      abundant(45) = 103

c Pd (105.903d0)
      call mass_addiso(46,101.905634d0) 
      call mass_addiso(46,103.904029d0) 
      call mass_addiso(46,104.905079d0) 
      call mass_addiso(46,105.903478d0)
      call mass_addiso(46,107.903895d0) 
      call mass_addiso(46,109.905167d0) 
      abundant(46) = 106

c Ag (ams 106.905d0)
      call mass_addiso(47,106.905092d0)
      call mass_addiso(47,108.904756d0) 
      abundant(47) = 107

c Cd (ams 113.903d0)
      call mass_addiso(48,105.906461d0) 
      call mass_addiso(48,109.903005d0) 
      call mass_addiso(48,110.904182d0) 
      call mass_addiso(48,111.902757d0) 
      call mass_addiso(48,112.904400d0) 
      call mass_addiso(48,113.903357d0)
      call mass_addiso(48,115.904755d0) 
      abundant(48) = 114

c In (ams 114.904d0)
      call mass_addiso(49,112.904061d0) 
      call mass_addiso(49,114.903882d0)
      abundant(49) = 115

c Sn (ams 119.902d0)
      call mass_addiso(50,115.901747d0) 
      call mass_addiso(50,116.902956d0) 
      call mass_addiso(50,117.901609d0)
      call mass_addiso(50,118.903311d0) 
      call mass_addiso(50,119.902199d0)
      call mass_addiso(50,121.903440d0) 
      call mass_addiso(50,123.905274d0) 
      abundant(50) = 120

c Sb (ams 120.903d0)
      call mass_addiso(51,120.903821d0) 
      call mass_addiso(51,122.904216d0) 
      abundant(51) = 121

c Te (ams 127.904d0)
      call mass_addiso(52,121.903050d0) 
      call mass_addiso(52,123.902818d0) 
      call mass_addiso(52,124.904429d0)  
      call mass_addiso(52,125.903310d0) 
      call mass_addiso(52,127.904463d0) 
      call mass_addiso(52,129.906229d0)
      abundant(52) = 130

c I (ams 126.9043d0)
      call mass_addiso(53,126.904473d0)
      abundant(53) = 127

c Xe (ams 131.904d0)
      call mass_addiso(54,127.903531d0) 
      call mass_addiso(54,128.904780d0) 
      call mass_addiso(54,129.903509d0) 
      call mass_addiso(54,130.905072d0) 
      call mass_addiso(54,131.904144d0)
      call mass_addiso(54,133.905395d0) 
      call mass_addiso(54,135.907214d0) 
      abundant(54) = 132

c 55-86 added 26 Oct 1997 by C. Bayse
c ref. Blatt, Modern Physics
c Cs (ams 132.9054d0)
      call mass_addiso(55,132.905429d0)
      abundant(55)=133

c Ba (ams 137.33d0)
      call mass_addiso(56,129.906282d0)
      call mass_addiso(56,131.905042d0)
      call mass_addiso(56,133.904486d0)
      call mass_addiso(56,134.905665d0)
      call mass_addiso(56,135.904553d0)
      call mass_addiso(56,136.905812d0)
      call mass_addiso(56,137.905232d0)
      abundant(56)=138

c La (ams 138.9055d0)
      call mass_addiso(57,137.907105d0)
      call mass_addiso(57,138.906346d0)
      abundant(57)=139

c Ce (ams 140.12d0)
      call mass_addiso(58,135.907140d0)
      call mass_addiso(58,137.905985d0)
      call mass_addiso(58,139.905433d0)
      call mass_addiso(58,141.909241d0)
      abundant(58)=140

c Pr (ams 140.9077d0)
      call mass_addiso(59,140.907647d0)
      abundant(59)=141

c Nd (ams 144.24d0)
      call mass_addiso(60,141.907719d0)
      call mass_addiso(60,142.909810d0)
      call mass_addiso(60,144.912570d0)
      call mass_addiso(60,145.913113d0)
      call mass_addiso(60,147.916889d0)
      call mass_addiso(60,149.920887d0)
      abundant(60)=142

c Pm (ams 145.d0)
      call mass_addiso(61,144.912743d0)
      abundant(61)=145

c Sm (ams 150.4d0)
      call mass_addiso(62,143.911998d0)
      call mass_addiso(62,146.914895d0)
      call mass_addiso(62,147.914820d0)
      call mass_addiso(62,148.917181d0)
      call mass_addiso(62,149.917273d0)
      call mass_addiso(62,151.919729d0)
      call mass_addiso(62,153.922206d0)
      abundant(62)=152

c Eu (ams 151.96d0)
      call mass_addiso(63,150.919847d0)
      call mass_addiso(63,152.921225d0)
      abundant(63)=153

c Gd (ams 157.25d0)
      call mass_addiso(64,151.919786d0)
      call mass_addiso(64,153.920861d0)
      call mass_addiso(64,154.922618d0)
      call mass_addiso(64,155.922118d0)
      call mass_addiso(64,156.923956d0)
      call mass_addiso(64,157.924099d0)
      call mass_addiso(64,159.927049d0)
      abundant(64)=158

c Tb (ams 158.9254d0)
      call mass_addiso(65,158.925342d0)
      abundant(65)=159

c Dy (ams 162.50d0)
      call mass_addiso(66,155.925277d0)
      call mass_addiso(66,157.924403d0)
      call mass_addiso(66,159.925193d0)
      call mass_addiso(66,160.926930d0)
      call mass_addiso(66,161.926795d0)
      call mass_addiso(66,162.928728d0)
      call mass_addiso(66,163.929171d0)
      abundant(66)=164

c Ho (ams 164.9304d0)
      call mass_addiso(67,164.930319d0)
      abundant(67)=165

c Er (ams 167.26d0)
      call mass_addiso(68,161.928775d0)
      call mass_addiso(68,163.929198d0)
      call mass_addiso(68,165.930290d0)
      call mass_addiso(68,166.932046d0)
      call mass_addiso(68,167.932368d0)
      call mass_addiso(68,169.935461d0)
      abundant(68)=166

c Tm (ams 168.9342d0)
      call mass_addiso(69,168.934212d0)
      abundant(69)=169

c Yb (ams 173.04d0)
      call mass_addiso(70,167.933894d0)
      call mass_addiso(70,169.934759d0)
      call mass_addiso(70,170.936323d0)
      call mass_addiso(70,171.936378d0)
      call mass_addiso(70,172.938208d0)
      call mass_addiso(70,173.938859d0)
      call mass_addiso(70,174.942564d0)
      abundant(70)=174

c Lu (ams 174.967d0)
      call mass_addiso(71,174.940770d0)
      call mass_addiso(71,175.942679d0)
      abundant(71)=175

c Hf (ams 178.49d0)
      call mass_addiso(72,173.940044d0)
      call mass_addiso(72,175.941406d0)
      call mass_addiso(72,176.943217d0)
      call mass_addiso(72,177.943696d0)
      call mass_addiso(72,178.945812d0)
      call mass_addiso(72,179.946545d0)
      abundant(72)=180

c Ta (ams 180.9479d0)
      call mass_addiso(73,179.947462d0)
      call mass_addiso(73,180.947992d0)
      abundant(73)=181

c W (ams 183.85d0)
      call mass_addiso(74,179.946701d0)
      call mass_addiso(74,181.948202d0)
      call mass_addiso(74,182.950220d0)
      call mass_addiso(74,183.950928d0)
      call mass_addiso(74,185.954357d0)
      abundant(74)=184

c Re (ams 186.207d0)
      call mass_addiso(75,184.952951d0)
      call mass_addiso(75,186.955744d0)
      abundant(75)=187

c Os (ams 190.2d0)
      call mass_addiso(76,183.952488d0)
      call mass_addiso(76,185.953830d0)
      call mass_addiso(76,186.955741d0)
      call mass_addiso(76,187.955860d0)
      call mass_addiso(76,188.958137d0)
      call mass_addiso(76,189.958436d0)
      call mass_addiso(76,191.961467d0)
      abundant(76)=192

c Ir (ams 192.22d0)
      call mass_addiso(77,190.960584d0)
      call mass_addiso(77,192.962917d0)
      abundant(77)=193

c Pt (ams 195.09d0)
      call mass_addiso(78,189.959917d0)
      call mass_addiso(78,191.961019d0)
      call mass_addiso(78,193.962655d0)
      call mass_addiso(78,194.964766d0)
      call mass_addiso(78,195.964926d0)
      call mass_addiso(78,197.967869d0)
      abundant(78)=195

c Au (ams 196.9665d0)
      call mass_addiso(79,196.966543d0)
      abundant(79)=197

c Hg (ams 200.59d0)
      call mass_addiso(80,195.965807d0)
      call mass_addiso(80,197.966743d0)
      call mass_addiso(80,198.968254d0)
      call mass_addiso(80,199.968300d0)
      call mass_addiso(80,200.970277d0)
      call mass_addiso(80,201.970617d0)
      call mass_addiso(80,203.973467d0)
      abundant(80)=202

c Tl (ams 204.37d0)
      call mass_addiso(81,202.972320d0)
      call mass_addiso(81,204.974401d0)
      abundant(81)=205

c Pb (ams 207.2d0)
      call mass_addiso(82,203.973020d0)
      call mass_addiso(82,205.974440d0)
      call mass_addiso(82,206.975872d0)
      call mass_addiso(82,207.976627d0)
      abundant(82)=208

c Bi (ams 208.9804d0)
      call mass_addiso(83,208.980374d0)
      abundant(83)=209

c Po (ams 209.d0)
      call mass_addiso(84,208.982404d0)
      abundant(84)=209

c At (ams 210.d0)
      call mass_addiso(85,209.987126d0)
      abundant(85)=210

c Rn (ams 222.d0)
      call mass_addiso(86,222.017570d0)
      abundant(86)=222

      end
c***********************************************************************
c
c  mass_numvec:  retrieve number of atomic mass vectors currently 
c                stored, to allow loops over isotopomers in analysis
c                modules
c
c***********************************************************************

      integer function mass_numvec()
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
      mass_numvec = nmassv
      end      

c***********************************************************************
c
c mass_addiso:  add a new isotopic mass - no checking if it's
c               already there (see mass_newiso for this)
c
c***********************************************************************
      subroutine mass_addiso(nz,w)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
      integer nz
      REAL w
      if(niso .eq. maxiso)call caserr('too many isotope masses')
      niso = niso + 1
      isotpz(niso) = nz
      wisotp(niso) = w
      end

c***********************************************************************
c
c mass_rdvec: input a mass vector
c
c   the user should provide the vector in the input atom ordering, and
c   the data is stored with this ordering.
c
c   ivec    vector to be input (-1 means create new vector) The vector
c           handle is returned
c
c***********************************************************************
      function mass_rdvec(izittel)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/work)
INCLUDE(common/errcodes)
      integer ivec, iat, mass_rdvec, izittel

      ivec = izittel
      if(ivec .eq. -1)then
         nmassv = nmassv + 1
         if (nmassv .gt. maxmassv)
     &    call gamerr('attempt to define too many atomic mass vectors',
     &        ERR_NO_CODE, ERR_FIXED_DIMENSION, ERR_SYNC, ERR_NO_SYS)

         ivec = nmassv
      endif

      if(jrec .eq. jump)then
         do iat = 1,nat
            imassv(ivec,iat) = 0
         enddo
      else
         do iat = 1,nat
            call inpi(imassv(ivec,iat))
         enddo
      endif

      mass_rdvec = ivec
      end

c***********************************************************************c
c
c mass_editvec: modify an isotopic assignment
c
c  if iso .eq. -2, the mass is user defined (ie not an entry in the
c                  nucleide table) and wght must be provided. in this
c                  case ivec must equal 1
c
c***********************************************************************c
      subroutine mass_editvec(ivec,iat,iso,wght)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/infob)
      integer iat,ivec,iso
      REAL wght, dum

      external amass_get
      REAL amass_get

      if (ivec .gt. maxmassv)
     &     call caserr('attempt to edit invalid mass vector')

      if(iso .eq. -2)then
c
c input user-defined weight
c
         if(ivec .ne. 1)then
            call caserr('cant edit user-defined mass for vector .ne.1')
         endif
c
c write into infob table to account for allow for reordering later
c
         amass(iat) = wght
         imassv(ivec,iat) = -2

      else
c 
c simply subsititute the requested isotope into the table
c
         imassv(ivec,iat) = iso
c
c check we can get the value, in case we compute then second derivs 
c and then fail
c
         dum = amass_get(1,iat)

      endif

      end

c**********************************************************************
c
c amass_get : obtain the mass of an atom
c
c     (main access routine for the masses tables)
c
c  ivec     mass vector to interrogate, ivec = 1 is the default vector
c           which will be initialised to the most abundant masses
c
c  iat      atom index after reordering due to symmetry
c
c*********************************************************************

      function amass_get(ivec,iat)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/infob)
INCLUDE(common/runlab)
INCLUDE(common/errcodes)
INCLUDE(common/iofile)
      integer ivec, iat, isubst, iz, iat2
      REAL amass_get
      REAL wnucld
c
      iat2 = mapmas(iat)

      if(imassv(ivec,iat2) .eq. -2)then
c
c user-defined mass (from common/infob)
c
         amass_get = amass(iat2)
c
      else
c
c a stored isotope 
c
         iz = isubst(zaname(iat))
         if(iz .le. 0)then
            if (iz.eq.0) then
c  bq centre
             amass_get = dble(max(0,imassv(ivec,iat2)))
             return
            else
             write(iwr,*)'problem with atom tag: ',zaname(iat)
             call gamerr(
     &            "failed to decode element type from atom tag",
     &            ERR_NO_CODE, ERR_NO_CLASS, ERR_SYNC, ERR_NO_SYS)
            endif
         endif
         amass_get = wnucld(iz,imassv(ivec,iat2))
c
      endif
      end

c**********************************************************************
c  mass_print - summary print of all mass data stored
c**********************************************************************

      subroutine mass_print
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      integer iat, ivec

      external amass_get
      REAL amass_get

      external zmass_nuclab
      character*8 zmass_nuclab

      character*8 ztmp(maxmassv)

      write(iwr,100)nmassv
      write(iwr,101)(' ',ivec,ivec=1,nmassv)
      write(iwr,102)(' ',ivec=1,nmassv)

      do iat=1,nat
         do ivec=1,nmassv
            ztmp(ivec) = zmass_nuclab(ivec,iat)
         enddo
         write(iwr,103)iat,(ztmp(ivec),
     +        amass_get(ivec,iat), ivec=1,nmassv)
      enddo
 100  format(/,1x,'Number of nuclear mass vectors stored:',i3,/)
 101  format(2x,'atom',5x,a1,10(4x,'vector',i2,4x,a1))
 102  format(2x,'----',4x,a1,10(1x,'----------------',a1))
 103  format(3x,i3,5x,10(a6,1x,f10.6,1x))
      end

c**********************************************************************
c  mass_prvec - print of a single atomic mass vector
c**********************************************************************

      subroutine mass_prvec(ivec)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/runlab)

      integer iat, ivec

      external amass_get
      REAL amass_get

      external zmass_nuclab
      character*8 zmass_nuclab, ztemp

      write(iwr,*)'Nuclear masses: '
      do iat=1,nat
         ztemp = zmass_nuclab(ivec,iat)
         write(iwr,100)iat,zaname(iat),ztemp(1:6),
     &        amass_get(ivec,iat)
      enddo
 100  format(i3,3x,a8,5x,a6,f12.8)
      end

c**********************************************************************
c mass_permute : store ordering information (called once from atoms2)
c**********************************************************************

      subroutine mass_permute(map)
      implicit none
      integer map(*)
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
       integer iat
       do iat=1,nat
          mapmas(map(iat)) = iat
       enddo
       end

c**********************************************************************
c  zmass_nuclab: return a label for the nucleus for printing:
c                in one of the following forms:
c  
c   c-13  - specific nucleus
c  <ave>  - average of natural abundance
c  <user> - user defined mass
c**********************************************************************

      function zmass_nuclab(ivec,iat)
      implicit none
      character*8 zmass_nuclab
INCLUDE(common/sizes)
INCLUDE(common/masses)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/periodic)
INCLUDE(common/errcodes)
INCLUDE(common/iofile)

      integer ivec, iat

      external isubst
      integer isubst

      character*2 zzz
      character*3 zm
      integer nz, iso, first, i, iat2

      nz = isubst(zaname(iat))

      if(nz .lt. 0)then
         write(iwr,*)'problem with atom tag: ',zaname(iat)
         call gamerr("failed to decode element type from atom tag",
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_SYNC, ERR_NO_SYS)
      endif
c
c data is stored in users input ordering, translate the
c atom index
c
      iat2 = mapmas(iat)

      if(imassv(ivec,iat2) .eq. -1)then
c
c average
c
         zmass_nuclab = "<ave>"
        
      else if(imassv(ivec,iat2) .eq. -2)then
c
c user input mass
c
         zmass_nuclab = "<user>"

      else
c
c specific isotope
c
         iso = imassv(ivec,iat2)

         if (nz.eq.0.and.iso.eq.0) then
            write(iwr,*)'Problem with atom ',zaname(iat)
            write(iwr,*)'No defaults available for BQ centers'
            write(iwr,*)'Please use WEIGHTS directive to specify a mass'
            call gamerr("No default mass for BQ centers",
     &           ERR_NO_CODE, ERR_NO_CLASS, ERR_SYNC, ERR_NO_SYS)
         endif

         if(iso .eq. 0)iso = abundant(nz)
 
         zzz = zelem(nz)
         if(zzz(2:2) .eq. ' ') zzz = ' '//zelem(nz)(1:1)
 
         write(zm,'(i3)')iso
         do i=3,1,-1
            if(zm(i:i).ne.' ')first=i
         enddo

         zmass_nuclab = zzz//'-'//zm(first:3)

      endif
      return
      end
c**********************************************************************
c omass_mod - returns .true. if non-default masses were requested
c**********************************************************************
      function omass_mod()
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/masses)
      logical omass_mod
      integer iat
      omass_mod = .false.
      if(nmassv .ne. 1)then
         omass_mod = .true.
      else 
         do iat=1,nat
            if (imassv(1,iat) .ne. 0)omass_mod = .true.
         enddo
      endif
      return
      end   
c**********************************************************************
c
c Old tables for reference
c
c from mains.m (common/datprp) used in properties package
c
c      data atwt/
c     *1.0078252d0,4.0026036d0,7.016005d0,9.012186d0,11.0093051d0,
c     *12.0d0,14.0030744d0,15.9949149d0,18.998405d0,19.9924404d0,
c     *22.989773d0,23.985045d0,26.981535d0,27.976927d0,30.973763d0,
c     *31.972074d0,34.968854d0,39.962384d0,38.963714d0,39.962589d0,
c     *44.955919d0,47.947948d0,50.943978d0,52.940514d0,54.938054d0,
c     *55.93493d0,58.933189d0,57.93534d0,62.92959d0,63.929145d0/
c
c  original ams (averaged masses) - now stored as wave
c
cc    after some discussion, it has been decided to replace the
cc    averaged isotope at. weights with values for the
cc    leading isotope .. should have minor impact on calc.
cc    vib. frequencies (chap2.vtab to be modified ..)
c      data ams/1.00797d0,4.0026d0,6.939d0,9.0122d0,10.811d0,
c     *12.01115d0,14.0067d0,15.9994d0,18.9984d0,20.183d0,
c     *22.9898d0,24.312d0,26.9815d0,28.086d0,30.9738d0,
c     *32.064d0,35.453d0,39.948d0,39.102d0,40.08d0,44.956d0,
c     *47.90d0,50.942d0,51.996d0,54.9380d0,55.85d0,58.9332d0,
c     *58.71d0,63.54d0,65.37d0,69.72d0,72.59d0,74.9216d0,
c     *78.96d0,79.909d0,83.80d0,85.47d0,87.62d0,88.905d0,
c     *91.22d0,92.906d0,95.94d0,99.0d0,101.07d0,102.905d0,
c     *106.4d0,107.87d0,112.40d0,114.82d0,118.69d0,121.75d0,
c     *127.60d0,126.9044d0,131.30d0/
c
c revised ams
c
c      data ams/
c     + 1.0078252d0, 4.0026036d0,7.016005d0,9.012186d0,11.0093051d0,
c     + 12.000000d0,14.0030744d0,15.9949149d0,18.998405d0,19.9924404d0,
c     + 22.989773d0,23.985045d0,26.981535d0,27.976927d0,30.973763d0,
c     + 31.972074d0,34.968854d0,39.962384d0,38.963714d0,39.962589d0,
c     + 44.955919d0,47.947948d0,50.943978d0,51.940514d0,54.938054d0,
c     + 55.93493d0,58.933189d0, 57.93534d0,62.92959d0,63.929145d0,
c     + 68.9257d0,73.921d0,74.921d0,79.916d0,78.914d0,83.911d0,
c     + 84.911d0,87.905d0,88.905d0, 89.904d0, 92.906d0,97.906d0,
c     + 98.906d0,101.903d0,102.904d0,105.903d0,106.905d0,113.903d0,
c     + 114.904d0,119.902d0,120.903d0,127.904d0,126.9043d0,131.904d0/
c**********************************************************************
c
      subroutine ver_mass(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mass.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
