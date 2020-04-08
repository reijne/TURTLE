      module junk_dma

      save

      real(kind=8) :: ox, oy, oz, tol, tshift
      real(kind=8) d(20,20)
      real(kind=8), target :: qx(121)
      real(kind=8) binom(20,20), rt(20)
      real(kind=8) gx(112), gy(112), gz(112)

      logical lshb, linear, notscf,ladd,lcorr,lmos,nuclei

      integer nc, lmax, mindc, maxdc, kw, ida

      real(kind=8), allocatable, dimension(:,:) :: x
      integer, allocatable, dimension(:) :: limit
      real(kind=8), allocatable, dimension(:) :: radius

      character*8, allocatable, dimension(:) :: name

      real(kind=8), pointer :: 
     &           q00, q10,q11c,q11s, q20,q21c,q21s,q22c,q22s,
     &           q30,q31c,q31s,q32c,q32s,q33c,q33s,
     &           q40,q41c,q41s,q42c,q42s,q43c,q43s,q44c,q44s,
     &           q50,q51c,q51s,q52c,q52s,q53c,q53s,q54c,q54s,q55c,q55s,
     &           q60,q61c,q61s,q62c,q62s,q63c,q63s,q64c,q64s,q65c,q65s,
     &           q66c,q66s


      contains

      subroutine alloc_junk_dma(natoms)
      integer, intent(in) :: natoms
      integer istat
      allocate(x(3,natoms),stat=istat)
      allocate(radius(natoms),stat=istat)
      allocate(limit(natoms),stat=istat)
      allocate(name(natoms),stat=istat)

c provide aliases to lower multipole

      q00  => qx( 1)
      q10  => qx( 2)
      q11c => qx( 3)
      q11s => qx( 4) 
      q20  => qx( 5)
      q21c => qx( 6)
      q21s => qx( 7)
      q22c => qx( 8)
      q22s => qx( 9)
      q30  => qx(10)
      q31c => qx(11)
      q31s => qx(12)
      q32c => qx(13)
      q32s => qx(14)
      q33c => qx(15)
      q33s => qx(16)
      q40  => qx(17)
      q41c => qx(18)
      q41s => qx(19)
      q42c => qx(20)
      q42s => qx(21)
      q43c => qx(22)
      q43s => qx(23)
      q44c => qx(24)
      q44s => qx(25)
      q50  => qx(26)
      q51c => qx(27)
      q51s => qx(28)
      q52c => qx(29)
      q52s => qx(30)
      q53c => qx(31)
      q53s => qx(32)
      q54c => qx(33)
      q54s => qx(34)
      q55c => qx(35)
      q55s => qx(36)
      q60  => qx(37)
      q61c => qx(38)
      q61s => qx(39)
      q62c => qx(40)
      q62s => qx(41)
      q63c => qx(42)
      q63s => qx(43)
      q64c => qx(44)
      q64s => qx(45)
      q65c => qx(46)
      q65s => qx(47)
      q66c => qx(48)
      q66s => qx(49)

      end subroutine

      subroutine dealloc_junk_dma
      integer istat
      deallocate(x,stat=istat)
      deallocate(radius,stat=istat)
      deallocate(limit,stat=istat)
      deallocate(name,stat=istat)
      end subroutine

      end module junk_dma
