      program test_pack
      implicit none
c
c     The various packers and associated bit operations are a constant
c     source of headaches in GAMESS-UK. So I have written a little 
c     test program to be able to check whether the packers are doing
c     what they are supposed to be doing.
c
c     This program tests 
c
c       - subroutine pack
c       - subroutine unpack
c
c     The program does not produce any output if the data obtained
c     after packing and unpacking is the same as the data before (i.e.
c     the packing routines work correctly). The program will print all
c     failures though.
c
      integer maxnum      ! the maximum number of integers to pack
      parameter(maxnum = 256)
      integer iu(maxnum)  ! the original unpacked numbers
      integer iu2(maxnum) ! the unpacked numbers after packing/unpacking
      integer ip(maxnum)  ! the numbers packed
      integer i, j        ! counters
      integer ibits       ! the number of bits to pack
      logical oerror      ! was an error found already?
      integer intsize     ! the size of an integer in bits
      integer nbint       ! the number of N bit numbers that fit in 
                          ! an integer
_IF(i8)
      parameter(intsize = 64)
_ELSE
      parameter(intsize = 32)
_ENDIF
c begin old packers
c     The following 2 lines are needed for the old implementation of 
c     the packers. They are no longer needed for the new packers. So
c     once we are convinced that the new packers work properly these
c     lines can be deleted. (HvD, 2006/11/16)
c
INCLUDE(common/restrj)
      oclunp = .true.
c
c end old packers
c
      oerror = .false.
      ibits = 8
      do j = 1, 3
c
        nbint = intsize/ibits
c
        do i = 1, maxnum
          iu(i)  = (i-1)*(2**mod(i,ibits-7))
          iu2(i) = -1
          ip(i)  = -1
        enddo
c
        call pack(ip,ibits,iu,maxnum)
c
        call unpack(ip,ibits,iu2,maxnum)
c
        do i = 1, maxnum
          if (iu(i).ne.iu2(i)) then
            if (.not.oerror) then
              oerror = .true.
              write(*,'("test_pack: errors detected:")')
              write(*,*)
              write(*,'(a5,a5,a18,2a12)')"#bits","i","hex packed",
     &                                   "before","after"
            endif
            write(*,'(i5,i5,z18,2i12)')ibits,i,ip((i+nbint-1)/nbint),
     &                                 iu(i),iu2(i)
          endif
        enddo
        ibits = 2*ibits
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine caserr(string)
      implicit none
c
c     Dummy caserr to cure linking issues
c
      character*(*) string
c
      write(*,*)string
      stop 10
      end
c
c-----------------------------------------------------------------------
c
      subroutine setsto(nw,ival,iu)
      implicit none
c
c     Dummy setsto to solve linking problems with old unpack
c
c     Input:
c
      integer nw     ! the number of array elements
      integer ival   ! the value to set the array to
c
c     Output:
c
      integer iu(nw) ! the array to initialise
c
c     Local:
c
      integer i
c
      do i = 1, nw
        iu(i) = ival
      enddo
c
      end
