
c
c in this version ipsc uses generic routines wherever possible
c
c  NB - still have to check the "maximum wait count" message
c
c
c** generic walltime
c
c return wall time in seconds
c
      subroutine walltime(elapse)
c
c      Fortran90
c      using date_and_time
c      one can use tidajt as well
c      we are amazingly optimistic; a job may run for more than a year
      implicit none
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 timbeg,  tmlast,  timsec
      real*8 wtimbeg, wtmlast, wtimsec
      real*8 cpwa
      integer nodel
      common/timanb/timbeg,tmlast,timsec(14),
     *wtimbeg,wtmlast,wtimsec(14),cpwa,nodel
c
      real*8 elapse,ccpu,cpulft
      integer i,iboem,imon,ipg_nodeid,ienter,year0
      integer time_array(8),tim0_array(8),monday(24)
      integer*8 time,tim0
      save tim0,ienter,iboem,year0
      save monday
      data ienter/0/,iboem/0/
      data monday/31,28,31,30,31,30,31,31,30,31,30,31,
     1            31,28,31,30,31,30,31,31,30,31,30,31/
c
c time_array(1)
c year
c time_array(2)
c month of the year
c time_array(3)
c day of the month
c time_array(4)
c time offset with respect to UTC in minutes
c time_array(5)
c hour of the day
c time_array(6)
c minutes of the hour
c time_array(7)
c seconds of the minute
c time_array(8)
c milliseconds of the second
c  
c      call tidajt(zdate,ztime,zyear,zaccno,zanam,isecs)
c
      if (ienter.eq.0) then
         call date_and_time(values=tim0_array)
         year0=tim0_array(1)
         if ((year0/4)*4.eq.year0) then
c...       leapyear
            monday(2) = 29
         else if (((year0+1)/4)*4.eq.year0+1) then
c...      +1 leapyear
            monday(14) = 29
         end if
         elapse = 0.0d0
         ienter = 1
c...     tim0 wr, beginning of year
         tim0 = 0
         do i=1,tim0_array(2)
            tim0 = tim0+monday(i)*24*3600 
         end do
         tim0 = tim0 + (tim0_array(3)-1)*24*3600
         tim0 = tim0 + (tim0_array(5)-1)*3600
         tim0 = tim0 + (tim0_array(6)-1)*60
         tim0 = tim0 + (tim0_array(7)-1)
         tim0 = tim0*1000
         tim0 = tim0 + tim0_array(8)
       else
         call date_and_time(values=time_array)
         imon = time_array(2)
         if (time_array(1).ne.year0) imon = imon + 12

c...     time wr, beginning of year
         time = 0
         do i=1,imon
            time = time+monday(i)*24*3600 
         end do
         time = time + (time_array(3)-1)*24*3600
         time = time + (time_array(5)-1)*3600
         time = time + (time_array(6)-1)*60
         time = time + (time_array(7)-1)
         time = time*1000
         time = time + time_array(8)
c
         elapse = (time-tim0)*0.001d0
         ccpu = cpulft(1)
         if (ccpu.lt.elapse*cpwa.and.elapse.gt.500.0d0) then
            iboem = iboem + 1
c           print *, ' **** node ',ipg_nodeid(),' cpu  ',ccpu,' ****'
c           print *, ' **** node ',ipg_nodeid(),' wall ',elapse,' ****'
            if (iboem.gt.3) then
c              call caserr('cpu time less then walltime*.05')
            end if
         else
            iboem = 0
         end if
       end if
c
c..    old stuff disabled
c
      end
c We currently don't print the hostname on Windoze
c
c generic gethost
c
      subroutine gethost(host)
      implicit None

      character*44 host


      call hostcc(host,44)
      end
c     END OF GETHOST
c
c     
c**   charwall returns walltime in hours:minutes:seconds
c
      character*10 function charwall()
c     
c.. due to compiler problems we use concatenate action
c
      implicit none
      real*8 wall
      integer hours,minutes,seconds,h(4),m(2),s(2),i,n
      character*1 sym(0:10)
      data sym/'0','1','2','3','4','5','6','7','8','9',' '/
c     
      call walltime(wall)
      hours = wall/3600
      minutes = wall/60 - hours*60
      seconds = wall - hours*3600 - minutes*60 + 0.5
c
      n = 1000
      do i =1,4
         h(i) = hours/n
         hours = hours - h(i)*n
         n = n/10
         if (i.lt.4.and.h(i).eq.0) h(i) = 10
      end do
         m(1) = minutes/10
         m(2) = minutes - m(1)*10
         s(1) = seconds/10
         s(2) = seconds - s(1)*10
c
      charwall = sym(h(1))//sym(h(2))//sym(h(3))//sym(h(4))//':'//
     1           sym(m(1))//sym(m(2))//':'//sym(s(1))//sym(s(2))
c
      return
      end
c
c** generic cputime
c
c  return cputime in seconds  (user, system and total)
c
      subroutine gms_cputime(cpu)
      real*8 cpu(3)
c use C call
      call cputimecc(cpu)
      end
c
c** generic cpulft
c
      function cpulft(i)
      implicit none
      integer i
      real*8 cpulft, ak
      real*8 cpu(3)
      logical oex
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
c ------ return the clock time in seconds if i=1
c ------        the cpu time remaining if i=0
c
      call gms_cputime(cpu)
      ak=cpu(1)
      if(i.eq.1) then
         cpulft=ak
      else if(i.eq.0) then
         cpulft=timlim-ak
c
c  if file STOP exists kill GAMESS orderly
c  this inquire proves to be a fatal action on massively  machines
c  with a shared filesystem
c
         inquire(file='STOP',exist=oex)
         if (oex) cpulft = 0.5d0
      else
         write(6,*) ' called cpulft with invalid i ',i,' returning 0'
         cpulft = 0.0d0
      endif
      return
      end
c
c** generic err
c
      subroutine err(n)
      write(6,10)n
 10   format(//1x,' **** stop ',i4,' ****'//)
c
c use C exit routine
c
      call exitc(11)
      stop
      end
c
c** generic fmove
c
      subroutine fmove(a,b,n)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      do 3 i=1,n
 3    b(i)=a(i)
      return
      end
c
c** generic sptchk
c
      subroutine sptchk(text)
      implicit real*8  (a-h,o-z)
      character*(*) text
      character*16 logfil
      character*131 logbuf
      common/clogfl/logfil,logbuf(100)
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c ***
c *** write text to the file :hostchar:'logfil'
c ***
      if(.not.ologf)return
      if(logf.lt.100) then
        logf=logf+2
      else
        do 10 i=3,100
10          logbuf(i-2)=logbuf(i)
      endif
      return
      end
c
c** generic tidajt
c
      subroutine tidajt(zdate,ztime,zyear,zaccno,zanam,jtime)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer uid
      character*24 z24
      character*31 zlogg
      if(jtime.le.0) jtime=3600000
      call daycc(z24)
      call namcc(zlogg,31)
      ztime=z24(12:19)
      zdate=z24(5:11)
      zyear=z24(21:24)
      zanam=zlogg(1:8)
c
c account number
c
      call uidcc(uid)      
      write(zaccno,'(i8)') uid
      return
      end
cc_IF(cray)
cc_IF(unicos)
cc      cpud=tsecnd(cpud)
cc_ELSE
cc      cpud=second(cpud)
cc_ENDIF
cc      elapsd=timef(elapsd)*1.0d-3
cc_ENDIF
cc old FPS/Cyber calls
cc_IF1(f)      icode=sys$gettime(begin,ebegin)
cc_IF1(u)      call givtim(begin,ebegin)

cDEBUG
c     The gotcha functions are meant to help detecting explicit
c     unpacking operations (i.e. where clever tricks are played with
c     common blocks rather than calling unpack, see e.g. the checkin
c     message on version 1.19 of util3.m). Known culprits are
c
c        - sgmata
c        - gsup
c        - proc2
c        - sgmat
c
c     most likely this list is rather incomplete.
c
c     If the explicit unpacking stuff were removed the nshift tables in 
c     packers could go. Furthermore, it should then be possible to 
c     replace the 3 packers with 1 subroutine, and the same for the 
c     unpackers.
c
c     Gotcha functions take 8, 16 or 32 bit integers and change them
c     in imaginative ways. Anything goes subject to the condition that 
c     applying the function twice returns the original number i.e.
c     igotchax(igotchax(i)).eq.i must be true.
c
c     integer function igotcha8(i)
c     implicit none
c     integer i
c     igotcha8 = z'ff' - i
c     return
c     end
c     integer function igotcha16(i)
c     implicit none
c     integer i
c     igotcha16 = z'ffff' - i
c     return
c     end
c     integer function igotcha32(i)
c     implicit none
c     integer i
c     igotcha32 = z'ffffffff' - i
c     return
c     end
cDEBUG
c
c** generic pack
c
      subroutine pack(ip,nbits,iu,nw)
      implicit none
c
c     Pack takes the 8, 16, or 32 rightmost bits of the integers in
c     the unpacked array iu and packs them closely in the packed array
c     ip. The number of elements in iu does not have to be a multiple
c     of the integer size divided by nbits.
c
c     The implementation of the packers has varied considerably over
c     the years. Mostly due to the lack of standardisation of bit 
c     operations in Fortran. With Fortran 90 this has improved and a
c     sensible implementation using bit operations should now be
c     possible.
c
c     The implementation requires three bit operations which are defined
c     as [1]:
c
c     - iand(i,j): 
c
c       returns the logical AND of all the bits in i and the
c       corresponding bits in j. The arguments i and j must have the
c       same type parameter value, which is [also] the type parameter
c       value of the result.
c
c     - ior(i,j):
c
c       returns the logical inclusive OR of all the bits of i and the
c       corresponding bits of j. The arguments i and j must have the
c       same type parameter value, which is [also] the type parameter
c       value of the result.
c
c     - ishft(i,shift):
c
c       returns an integer, with the same type parameter as i, and value
c       equal to i except that the bits are shifted shift places to the
c       left (-shift places to the right if shift is negative). Zeros
c       are shifted in from the other end. The arguments shift must be
c       an integer with value satisfying the inequality
c       abs(shift) <= bit_size(i) (bit_size returns the number of bits
c       i can hold).
c
c     A final note: it seems that integer*1 is unsigned with most
c     compilers. I do not know whether this is defined in the fortran
c     standard, whether this is a defacto standard, or whether this is
c     implementation dependent. I tried looking this up but I could not
c     find any statements specific about integer*1 to decide this.
c
c     [1] "Fortran 95/2003 explained",
c         M. Metcalf, J. Reid, and M. Cohen,
c         Oxford University Press, Great Clarendon Street,
c         Oxford OX2 6DP
c         Published in 2004, ISBN 0-19-852693-8
c         pages 170-171.
c
c     Input:
c
      integer nw     ! the number of integers in iu that need packing
      integer nbits  ! the number of bits to pack from each integer
      integer iu(nw) ! the integers that need packing
c
c     Output:
c
      integer ip(*)  ! the closely packed integers
c
c     Code:
c
c ***   not the right number of words does occurr !!!
c
      if (nbits.eq.8) then
      call vipk8(iu,ip,nw)
      else if (nbits.eq.16) then
      call vipk16(iu,ip,nw)
      else if (nbits.eq.32) then
      call vipk32(iu,ip,nw)
      else
         call caserr(' *** pack called with nbits ne 8,16,32 ***')
      end if
c
      return
      end
c
c** generic vipk8
c
      subroutine vipk8(iu,ip,nw)
      implicit none
c
c     Pack the 8 rightmost bits of every integer in iu into an element
c     of ip such that iout becomes a closely packed representation
c     of iu. See subroutine pack for definitions of the required bit
c     operations.
c
c     Parameters:
c
      integer mx               ! the size of integer divided by 8
      parameter( mx = 8 )
c
c     Input:
c
      integer nw               ! the number of integers in iu to pack
      integer iu(nw)           ! the (unpacked) integers that need
                               ! packing
c
c     Output:
c
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Local:
c
      integer ibits8           ! all 8 rightmost bits set
      parameter(ibits8 = z'ff')
      integer ibits0           ! no bits set
      parameter(ibits0 = z'0')
      integer indu             ! the number of integers from iu packed
      integer indp             ! counter over integers in ip
      integer i                ! counter
      integer ibits            ! the 8 rightmost bits
      integer ishifted         ! the ibits shifted to the right position
c
cDEBUG
c     integer igotcha8
cDEBUG
      integer nshift(mx)   ! the number of bits to shift
                           ! required for explicit unpacking see gotcha
      data nshift/56,0,8,16,24,32,40,48/
c
c     Code:
c
      indu = 0
      do indp = 1, nw/mx
        ip(indp) = ibits0
        do i = 1, mx
          indu     = indu + 1
          ibits    = iand(iu(indu),ibits8)
cDEBUG
c         ibits    = igotcha8(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      enddo
c
c     if nw is not a multiple of mx then pack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        ip(indp) = ibits0
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ibits    = iand(iu(indu),ibits8)
cDEBUG
c         ibits    = igotcha8(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      endif
c       
      end
c
c** generic vipk16
c
      subroutine vipk16(iu,ip,nw)
      implicit none
c
c     Pack the 16 rightmost bits of every integer in iu into an element
c     of ip such that iout becomes a closely packed representation
c     of iu. See subroutine pack for definitions of the required bit
c     operations.
c
c     Parameters:
c
      integer mx               ! the size of integer divided by 16
      parameter( mx = 4 )
c
c     Input:
c
      integer nw               ! the number of integers in iu to pack
      integer iu(nw)           ! the (unpacked) integers that need
                               ! packing
c
c     Output:
c
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Local:
c
      integer ibits16          ! all 16 rightmost bits set
      parameter(ibits16 = z'ffff')
      integer ibits0           ! no bits set
      parameter(ibits0  = z'0')
      integer indu             ! the number of integers from iu packed
      integer indp             ! counter over integers in ip
      integer i                ! counter
      integer ibits            ! the 16 rightmost bits
      integer ishifted         ! the ibits shifted to the right position
cDEBUG
c     integer igotcha16
cDEBUG
c
      integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking
                               ! see gotcha
      data nshift/48,0,16,32/
c
      indu = 0
      do indp = 1, nw/mx
        ip(indp) = ibits0
        do i = 1, mx
          indu     = indu + 1
          ibits    = iand(iu(indu),ibits16)
cDEBUG
c         ibits    = igotcha16(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      enddo
c
c     if nw is not a multiple of mx then pack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        ip(indp) = ibits0
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ibits    = iand(iu(indu),ibits16)
cDEBUG
c         ibits    = igotcha16(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      endif
c       
      end
c
c** generic vipk32
c
      subroutine vipk32(iu,ip,nw)
      implicit none
c
c     Pack the 32 rightmost bits of every integer in iu into an element
c     of ip such that iout becomes a closely packed representation
c     of iu. See subroutine pack for definitions of the required bit
c     operations.
c
c     Parameters:
c
      integer mx               ! the size of integer divided by 32
      parameter( mx = 2 )
c
c     Input:
c
      integer nw               ! the number of integers in iu to pack
      integer iu(nw)           ! the (unpacked) integers that need
                               ! packing
c
c     Output:
c
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Local:
c
      integer ibits0           ! no bits set
      parameter(ibits0  = z'0')
      integer indu             ! the number of integers from iu packed
      integer indp             ! counter over integers in ip
      integer i                ! counter
      integer ibits            ! the 32 rightmost bits
      integer ishifted         ! the ibits shifted to the right position
cDEBUG
c     integer igotcha32
cDEBUG
c
      integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking
                               ! see gotcha

      integer ibits32          ! all 32 rightmost bits set (we hope)
      parameter(ibits32 = z'ffffffff')
      data nshift/32,0/
c
c     Code:
c
      indu = 0
      do indp = 1, nw/mx
        ip(indp) = ibits0
        do i = 1, mx
          indu     = indu + 1
          ibits    = iand(iu(indu),ibits32)
cDEBUG
c         ibits    = igotcha32(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      enddo
c
c     if nw is not a multiple of mx then pack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        ip(indp) = ibits0
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ibits    = iand(iu(indu),ibits32)
cDEBUG
c         ibits    = igotcha32(ibits)
cDEBUG
          ishifted = ishft( ibits, nshift(mod(indu,mx)+1) )
          ip(indp) = ior( ip(indp) , ishifted )
        enddo
      endif
c       
      end
c
c** generic unpack
c
      subroutine unpack(ip,nbits,iu,nw)
      implicit none
c
c     Unpack performs the inverse of the operation performed by pack.
c     Note that the number of integers to unpack does not have to be a
c     multiple of the integer size divided by nbits. See the subroutine 
c     pack for further details.
c
c     Input:
c
      integer nw     ! the number of integer to unpack
      integer nbits  ! the number of bits to unpack for every integer
      integer ip(*)  ! the closely packed integers
c
c     Output:
c
      integer iu(nw) ! the unpacked integers
c
c ***    ((nw/nbits)*nbits.ne.nw)  does happen
c
      if (nbits.eq.8) then
      call viup8(ip,iu,nw)
      else if (nbits.eq.16) then
      call viup16(ip,iu,nw)
      else if (nbits.eq.32) then
      call viup32(ip,iu,nw)
      else
         call caserr(' *** unpack called with nbits ne 8,16,32 ***')
      end if
c
      return
      end
c
c** generic viup32
c
      subroutine viup32(ip,iu,nw)
      implicit none
c
c     Unpack the 32 rightmost bits of every integer for iu from the
c     closely packed array ip. See the subroutine pack for the
c     definitions of the bit operations. The implementation of this
c     subroutine is derived of the implementation of vipk32. 
c
c     Parameters:
c
      integer mx               ! the size of integer divided by 32
      parameter( mx = 2 )
c
c     Input:
c
      integer nw               ! the number of integer to unpack into iu
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Output:
c
      integer iu(nw)           ! the unpacked integers
c
c     Local:
c
      integer ibits32          ! all 32 rightmost bits set (we hope)
      parameter(ibits32 = z'ffffffff')
      integer ibits0           ! no bits set
      parameter(ibits0 = z'0')
      integer indu             ! the number of unpacked integers
      integer indp             ! counter over elements in ip
      integer i                ! counter
      integer j                ! limit for the remaining integers
      integer ishifted         ! the bits of interest shifted to the
                               ! right position
      integer ibits            ! the bits of interest only
      integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking
                               ! see gotcha
      data nshift/32,0/
cDEBUG
c     integer igotcha32
cDEBUG
c
c     Code:
c
      indu = 0
      do indp = 1, nw/mx
        do i = 1, mx
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits32 )
cDEBUG
c         ibits    = igotcha32(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      enddo
c
c     if nw is not a multiple of mx then unpack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits32 )
cDEBUG
c         ibits    = igotcha32(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      endif
c       
      end
c
c** generic viup16
c
      subroutine viup16(ip,iu,nw)
      implicit none
c
c     Unpack the 16 rightmost bits of every integer for iu from the
c     closely packed array ip. See the subroutine pack for the
c     definitions of the bit operations. The implementation of this
c     subroutine is derived of the implementation of vipk16. 
c
c     Parameters:
c
      integer mx               ! the size of integer divided by 16
      parameter( mx = 4 )
c
c     Input:
c
      integer nw               ! the number of integer to unpack into iu
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Output:
c
      integer iu(nw)           ! the unpacked integers
c
c     Local:
c
      integer ibits16          ! all rightmost 16 bits set
      parameter(ibits16 = z'ffff')
      integer indu             ! the number of unpacked integers
      integer indp             ! counter over elements in ip
      integer i                ! counter
      integer j                ! limit for the remaining integers
      integer ishifted         ! the bits of interest shifted to the
                               ! right position
      integer ibits            ! the bits of interest only
      integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking
                               ! see gotcha
      data nshift/48,0,16,32/
cDEBUG
c     integer igotcha16
cDEBUG
c
c     Code:
c
      indu = 0
      do indp = 1, nw/mx
        do i = 1, mx
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits16 )
cDEBUG
c         ibits    = igotcha16(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      enddo
c
c     if nw is not a multiple of mx then unpack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits16 )
cDEBUG
c         ibits    = igotcha16(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      endif
c       
      end
c
c** generic viup8
c
      subroutine viup8(ip,iu,nw)
      implicit none
c
c     Unpack the 8 rightmost bits of every integer for iu from the
c     closely packed array ip. See the subroutine pack for the
c     definitions of the bit operations. The implementation of this
c     subroutine is derived of the implementation of vipk8. 
c
c     Parameters:
c
      integer mx               ! the size of integer divide by 8
      parameter( mx = 8 )
c
c     Input:
c
      integer nw               ! the number of integer to unpack into iu
      integer ip((nw+mx-1)/mx) ! the packed integers
c
c     Output:
c
      integer iu(nw)           ! the unpacked integers
c
c     Local:
c
      integer ibits8           ! all rightmost 8 bits set
      parameter(ibits8 = z'ff')
      integer indu             ! the number of unpacked integers
      integer indp             ! counter over elements in ip
      integer i                ! counter
      integer j                ! limit for the remaining integers
      integer ishifted         ! the bits of interest shifted to the
                               ! right position
      integer ibits            ! the bits of interest only
      integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking 
                               ! see gotcha
      data nshift/56,0,8,16,24,32,40,48/
cDEBUG
c     integer igotcha8
cDEBUG
c
c     Code:
c
      indu = 0
      do indp = 1, nw/mx
        do i = 1, mx
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits8 )
cDEBUG
c         ibits    = igotcha8(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      enddo
c
c     if nw is not a multiple of mx then unpack the remaining integers
c
      if (mod(nw,mx) .gt. 0) then
        indp = nw/mx+1
        do i = 1, mod(nw,mx)
          indu = indu + 1
          ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
          ibits    = iand( ishifted, ibits8 )
cDEBUG
c         ibits    = igotcha8(ibits)
cDEBUG
          iu(indu) = ibits
        enddo
      endif
c       
      end
c
c** generic strtrm
c
      subroutine strtrm( string, length)
c
c     called by:
c     calls to:
c
      character*(*)  string
      length=len(string)
      l=length
      do 5 i=1,l
      if(string(length:length).ne.' ') return
      length=length-1
 5    continue
      length=0
      return
      end
c
c** generic fget
c
      subroutine fget(gin,m,iunit)
      implicit real*8  (a-h,o-z)
      dimension gin(*)
      call find(iunit)
      call get(gin,m)
      return
      end
c
c** generic flags
c
      subroutine flags
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      write(iwr,6)
 6    format(//
     *35x,'*******************************************************',
     *'*********'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*               ===  G A M E S S - U K    ===          ',
     *'        *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'* Generalised Atomic and Molecular Electronic Structure',
     *' System *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*               ===   Opteron version  8.0  ===        ',
     *'        *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*******************************************************',
     *'*********'/)
      return
      end
c** generic initj
c
      subroutine initj(lword)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*44 host
      character*256 com
      character*10 cdate,cname,cversion
      character*5 ctime
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      character *8 zdate,ztime,zaccno,zanam
      common/jinfo/zdate,ztime,zaccno,zanam
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
c
c ----- initialise job
c
c
      call start_time_period(TP_ENTIRE)
      call inita(lword,mfree)
      call tidajt(zdate,ztime,zyear,zaccno,zanam,isecs)
c
      ti=0.0d0
      tx=0.0d0
      timlim=isecs
      call flags
      call gethost(host)
      call getarg(0,com)
      write(iwr,10) host, com
 10   format(
     *40x,'Hostname  : ',a44/
     *40x,'GAMESS-UK Executable: ',a70)
      call gms_version(cdate,ctime,cname,cversion)
      if (cversion.ne."") then
        write(iwr,12)cversion
      endif
12    format(40x,'Revision    ',a10)
      write (iwr,14) cdate,ctime,cname
14    format(40x,'Compiled on ',a10,' ',a5,' by ',a10)

      write(iwr,20)zanam,zdate(1:6),zyear(1:4),ztime,zaccno,
     &           isecs/60,lword,lword*8/(1000*1000)
 20   format(
     *40x,'job name   ',a8/
     *40x,'date    ',a6,1x,a4/
     *40x,'time       ',a8/
     *40x,'acct       ',a8//
     *40x,'job time specified ',i10,' minutes'//
     *40x,'main store requested',i10,' real*8  words',
     * ' =',i6,' mbytes'/)

c
c    data set initialisation
c
      call initb
c
      return
      end
c
c** generic iofail
c
      subroutine iofail(irw,itype)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      dimension ztext(3)
      data ztext/' input',' output',' search'/
c
      write(iwr,1)ztext(itype),yed(irw),irep
 1    format(/a8,' error on ',a4,' reply word = ',i3)
      call errors(62)
      return
      end
c
c** generic kilbuf
c
      subroutine kilbuf
      implicit real*8  (a-h,o-z)
      call clredx
      return
      end
c
c** generic gtnv
c
      subroutine gtnv(s,v)
      implicit real*8  (a-h,o-z)
      character s*(*), v*(*)
      logical ostart, ostop
      ostart = .false.
      ostop = .false.
      lens = len(s)
      lenv = len(v)
      do 10 i = 1,lens
         if(s(i:i).ne.' '.and..not.ostart)then
            istart=i
            ostart = .true.
         endif
         if(s(i:i).eq.' '.and..not.ostop)then
            istop=i-1
            ostop = .true.
         endif
10    continue
      if(.not.ostop)istop=lens
      call gtnvcc(s(istart:istop), istop-istart+1,v,lenv,ier)
      return
      end
c
c** generic quit
c
      subroutine quit(ia,n)
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *80 ia
      call whtps
      call shut
      call err(n)
      return
      end
c
c** generic whtps
c
      subroutine whtps
      implicit real*8  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/blksiz2/ns2(12),n2len,n2max
      common/blksi3/ns3(11),n3len,n3max
c...
      data ztexts,ztextp,ztextt/'sort ','psort','table'/
      do 3 k=1,maxlfn
      if(ipos(k))3,3,4
3     continue
      return
4     write(iwr,5000)
 5000  format(/1x,'file positions'/
     * 1x,'lfn','   block  length'/1x,19('=')/)
      do 1 i=k,maxlfn
      if(ipos(i))1,1,2
2     write(iwr,6000)yed(i),ipos(i),ipos(i)
1     continue
6000  format(1x,a4,2i8 )
      if(nsmax.gt.0)write(iwr,6001)ztexts,nslen,nsmax
      if(n2max.gt.0)write(iwr,6001)ztextp,n2len,n2max
      if(n3max.gt.0)write(iwr,6001)ztextt,n3len,n3max
6001  format(1x,a6,i6,i8)
      return
      end
c
c** generic print0
c
      subroutine print0
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *10 chpid, chnod
      character *132 rec132,ginput
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer pid
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical o_inclin,o_inc_err
c
c --- write a copy of the input 
c
c     
c   =====  define stream number for copy...
c
      o_inclin = .false.
      o_inc_err = .false.
      irdcpy=95

c
c  ====== define name for file
c
      call pidcc(pid)
      write(chpid,'(i10)')pid
      do 1020 i=1,10
      if(chpid(i:i).ne.' ') then
      length=i
      go to 1031
      endif
1020  continue

 1031 write(chnod,'(i10)')ipg_nodeid()

      do i=1,10
         if(chnod(i:i).ne.' ') then
            length2=i
            go to 1030
         endif
      enddo

1030  ginput='gamess_input.'//chpid(length:10)//'.'//chnod(length2:10)

c
c  spool data from ird to irdcpy
c
      open(unit=irdcpy,form='formatted',file=ginput,status='unknown')
cjvl
1060  if (.not.o_inclin) read(ird,1080,end=1090) rec132
c     include statement (allows including files in the input )
c
      call inclin(rec132,irdcpy-1,iwr,o_inclin,o_inc_err)
cjvl
      write(irdcpy,1080)rec132
      go to 1060
 1080 format(a132)
c
c ==== end of input file, rewind copy file 
c 
 1090 continue
      ird=irdcpy
      rewind ird
      return
      end
c
c** generic clredx (serial)
c
      subroutine clredx
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      return
      end
c
c** generic shut (serial)
c
      subroutine shut
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,
     *            ipos(maxlfn),nam(maxlfn),
     *  keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      call clredx
       irep = 0
      do 1 i=1,maxlfn
      if(ipos(i).gt.0 ) then
       if(.not.omem(i)) then
          if(keep(i).gt.0) then
          call closecc(nam(i))
          istat(i)=2
          else
          call closecc(nam(i))
          call strtrm(zedfil(i),length)
          call delcc(zedfil(i),length,ierrio)
          if(ierrio.ne.0)call ioerr('delete',i,zedfil(i))
          istat(i)=0
          endif
        if(irep.ne.0)then
        write(iwr,5)yed(i)
 5      format(1x,'***** error closing/deleting file',1x,a4)
           else
        if(ooprnt)write(iwr,3)yed(i)
 3      format(//' file ',a4,' closed')
        endif
       ipos(i)=0
       endif
      endif
 1    continue
c ... close fortran files
      do 10 i=1,maxfrt
      if(iposf(i).gt.0) then
      if(keepf(i).gt.0) then
         close(unit=i,iostat=irep,status='keep')
          else
         close(unit=i,iostat=irep,status='delete')
          endif
         if(irep.ne.0)then
          write(iwr,5)zftfil(i)
         else
          if(ooprnt)write(iwr,30)zftfil(i)
 30       format(//' fortran file ',a6,' closed')
         endif
      endif
 10   continue
      close(unit=95,status='delete')
      return
      end
c
c** generic find (serial)
c
      subroutine find(iunit)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      isel=iunit
      if (ooprnt) then
       write(iwr,10) isel
10     format(1x,'**** find unit', i3)
      endif
      return
      end
c
c** generic get (serial)
c
      subroutine get(c,nword)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer nchk(2)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension c(*)
c ... single-buffered  i/o
      common/bufa/abc(512)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      common/vinteg/q(1)
c
      irec=ipos(isel)
      if (irec.lt.1) call iofail(isel,1)
      if (ooprnt) then
       write(iwr,10) isel, irec
 10    format(1x,'**** get from unit', i3,
     +           ' at block',i5)
      endif
c
      if(omem(isel)) then
       iword = (irec-1)*512+iqqoff(isel)+1
       call upack2(q(iword+511),iibl,nword)
       call dcopy(nword,q(iword),1,c,1)
      else
c
c --- this is for the non multi-buffered streams
c
         call getcc(nam(isel),abc,ierrio)
         if(ierrio.ne.0)call ioerr('read',isel,' ')
c
c --- prepare to return results to the program
c
         call upack2(abc(512),iibl,nword)
         if (ooprnt) then
          write(iwr,20) iibl, nword
 20       format(/1x,'**** get: iblock,nword = ',
     +    2i6)
         endif
         call dcopy(nword,abc(1),1,c,1)
      endif
      ipos(isel)=irec+1
      return
      end
c
c** generic search (serial)
c
      subroutine search(iblk,iunit)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
c --- position a stream at the necessary point
c
      if (ooprnt) then
       write(iwr,10) iunit, iblk
 10    format(1x,'**** position unit', i3,
     +           ' at block',i5)
      endif
c
      if(iblk.le.0)call iofail(iunit,3)
      if(ipos(iunit).le.0)call rdtpnm(iunit)
      if(.not.omem(iunit)) then
        call srchcc(nam(iunit),iblk,ierrio)
        if(ierrio.ne.0)call ioerr('search',iunit,' ')
      endif
      ipos(iunit)=iblk
      return
      end
c
c** generic inita (serial)
c
      subroutine inita(lword,mfree)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      character *16 logfil
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension ytext(20),zinout(2)
      common/clogfl/logfil
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/bufa/aaaaa(8190)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/blksiz2/ns2,ns2512(10),npstat,
     *  n2len,n2max ,nport
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     *  n3len,n3max ,ntort
      common/restri/nfile(63),lda(508),isect(508),ldx(508),
     *  iacsct(508)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      dimension yfname(maxfrt)
      dimension yynam(maxlfn),zdef(maxfrt)
c dump core on caserr
      common/crdmp/odumpc
c system derived error codes
      character mess*100
      common/syscod/icode1,mess
c
      data yynam/'ed0','ed1' ,'ed2' ,'ed3' ,'ed4' ,'ed5' ,'ed6' ,
     *          'ed7' ,'ed8' ,'ed9' ,'ed10','ed11','ed12','ed13',
     *          'ed14','ed15','ed16','ed17','ed18','ed19',
     *          'mt0' ,'mt1' ,'mt2' ,'mt3' ,'mt4' ,'mt5' ,'mt6' ,
     *          'mt7' ,'mt8' ,'mt9' ,'mt10','mt11','mt12','mt13',
     *          'mt14','mt15','mt16','mt17','mt18','mt19'
     *         /
      data yfname/
     *'ft01','ft02','ft03','ft04','ft05','ft06','ft07','ft08',
     *'ft09','ft10','ft11',
     *'ft12','ft13','ft14','ft15','ft16','ft17','ft18','ft19',
     *'ft20','ft21','ft22','ft23','ft24','ft25','ft26','ft27',
     *'ft28','ft29','ft30','ft31',
     *'ft32','ft33','ft34','ft35','ft36','ft37','ft38','ft39',
     *'ft40','ft41','ft42','ft43','ft44','ft45','ft46','ft47',
     *'ft48','ft49','ft50','ft51',
     *'ft52','ft53','ft54','ft55','ft56','ft57','ft58','ft59',
     *'ft60'/
      data zdef/
     *'ftn001','ftn002','ftn003','ftn004','ftn005','ftn006','ftn007',
     *'ftn008','ftn009','ftn010','ftn011',
     *'ftn012','ftn013','ftn014','ftn015','ftn016','ftn017','ftn018',
     *'ftn019','ftn020','ftn021','ftn022','ftn023','ftn024','ftn025',
     *'ftn026','ftn027','ftn028','ftn029','ftn030','ftn031',
     *'ftn032','ftn033','ftn034','ftn035','ftn036','ftn037','ftn038',
     *'ftn039',
     *'ftn040','ftn041','ftn042','ftn043','ftn044','ftn045',
     *'ftn046','ftn047','ftn048','ftn049','ftn050','ftn051',
     *'ftn052','ftn053','ftn054','ftn055','ftn056','ftn057',
     *'ftn058','ftn059','ftn060'/
      data ytext/'chan','swit','file','seti','memo','stor','core',
     * 'lpag','time','sort','psor','driv','dcor','max','logf',
     * 'tabl','node','dire','para','noec'/
      data zblank,zkeep,zold,zlength/' ','keep','old','length'/
      data zdel/'delete'/
      data yout,yprin/'outp','prin'   /
      data  zinout/  'dataout','punch'/
c
      nslen = 0
      n2len = 0
      n3len = 0
      nsmax = -1
      n2max = -1
      n3max = 0
      nsz=10
      ns2=1
      ns3=1
      nsort=101
      nport=102
      ntort=103
      ichek=0
      isel= 0
      iselr= 0
      iselw= 0
      irep= 0
      do 4000 loop=1,maxlfn
4000  yed(loop) = yynam(loop)
      do 3100 loop = 1,maxfrt
      yft(loop)   =   yfname(loop)
      zftstat(loop)=   'unknown'
      oform(loop) = .false.
      zftn(loop) =    zdef(loop)
      zftfil(loop)=   zdef(loop)
      iposf(loop)=0
      keepf(loop)=0
3100  istatf(loop)=0
      numf = 60
c
c rwah default file settings..
c
      do 3500 loop=1,maxlfn
      iop(loop)  = 0
      keep(loop) = 0
      istat(loop) = 1
      omem(loop) = .false.
      ipos(loop) = 0
      oadd(loop) = .false.
3500  nam (loop) = numf+loop
c
c === hardwire /machin/ values
c
      m511 = 511
c
c  set default format of integral storage based on number of
c  basis functions (this controls all packing/unpacking
c  routines, using the parameters in atmblk).
c
      num2e   = 340
      num2ep  = 340
      num2ejk = 204
      o255i   = .true.
      lab1632 = 16
      lab816  =  8
      numlab  = 1360
      numlabp = 680
      numlabjk = 408
c
      nav = lenwrd()
      mach12 = nav
c
      mach(1)=1200+2/nav
      mach(2)=(7*mxshel+4)/nav
      if(nav-1)1111,1111,1112
1111  mach(3)=72
      mach(11)=6*mach(3)
      mach(12)=mach(11)+mach(3)+1
      go to 1113
1112  mach(3)=77
      mach(11)=6*mach(3)*nav
      mach(12)=1012
1113  mach(4)=-1
c     mach(4)=12*maxat+3+4/nav (disabled, cf rdrec1/server)
      mach(5)=699+(13+nav-1)/nav
      mach(6)=682/nav
      mach(7)=3*maxat+(5+nav-1)/nav
      mach(8)=maxorb+maxorb+1+6/nav
      mach(9)=mxorb3+(2*maxorb+mxorb3+nav)/nav
      mach(10)=200+600/nav
      mach(13)=(maxorb*2+65+nav)/nav
      nw=maxorb+3
      mach(14)=(nw+nw+nav-1)/nav
      mach(15)=maxorb+2 +
     *         (4*maxorb+5+nav)/nav
      mach(16)=(4*maxorb+3+nav)/nav
      mach(17)=3600 + (maxorb+7+nav)/nav
      mach(18)=(4-nav)*511
      mach(19)=0
      mach(20)=0
      mvadd = 1 + lensec(mach(8)) + lensec(mach(9))
c
c --- set default section numbers
c
      do loop=1,204
       isect(loop)=loop+200
      enddo
      do loop=205,300
       isect(loop)=0
      enddo
      do loop=301,400
       isect(loop) = loop - 300
      enddo
      do loop = 401, 508
       isect(loop)=loop
      enddo
      do loop=1,508
       iacsct(loop)=0
      enddo
c
c ... and lengths
c
      ldx(isect(501)) = 50
      lda(isect(501)) = 30 + 2674/nav
c
      lword=icoresize()
      otime=.false.
      icode1 = 0
      odumpc = .false.
      oprintm = .false.
      omemin = .false.
      call gmem_set_nanify(.false.)
      call gmem_set_debug(.false.)
      call gmem_set_quiet(.true.)
      isecs=0
      mfree=0
      call setsto(maxlfn,0,len)
      call setsto(maxfrt,0,nrecf)
      call setsto(maxfrt,0,lrecf)
      zzz=' '
      do 110 i=1,maxlfn
      zzz=yed(i)(1:4)
 110  zedfil(i)=zzz
c
c --- change the default to 'keep' if environment variable set
c
      do 115 i=1,maxlfn
       call gtnv(yed(i),filatt)
       if(filatt.ne.' ') keep(i)=1
 115  continue
c
c ---- and repeat for the fortran files
c
      do 116 i=1,maxfrt
       call gtnv(zftn(i),filatt)
       if(filatt.ne.' ') keepf(i)=1
 116  continue
c
c --- take a copy of input file
c
      call print0
      oecho = .false.
      if(opg_root()) then
       write(iwr,500)
 500   format(/1x,'**********************************'/
     +         1x,'**** GAMESS-UK Input Data Set ****'/
     +         1x,'**********************************'/)
       oecho = .true.
      endif
      opark = .false.
      oreact = .false.
      odci = .false.
      odlfi = .false.
c
c --- default output file ?
c
      call input
      call inpa4(yyy)
      if(yyy.ne.yout)go to 300
      call inpa(zinout(1))
      call inpa4(yyy)
      if(yyy.eq.yprin)ooprnt=.true.
      call input
      go to 310
 300  jrec=jrec-1
 310  call prinit(oreact,opark,odci,odlfi)
      oecho = .false.
      go to 320
 10   call input
 320  call inpa4(yyy)
      j=locatc(ytext,20,yyy)
      if(j)20,30,20
 20   go to (40,40,40,120,50,50,50,60,70,200,240,210,245,250,900,285,
     * 290,290,290,290),j
c
c  max
c
250   call inpa(ztext)
      if(ztext.ne.'off')go to 10
      omax=.false.
      go to 10
c
c  logfile '<logfile>' (up to 16 characters)
c
900   ologf=.true.
      call inpanpc(logfil)
      go to 10
c
c  drives
c
 210  call inpa(ztext)
      go to 10
c
c  sort
c
 200   call inpi(nsz)
       go to 10
c
c  psort
c
 240   call inpi(ns2)
       go to 10
c
c  table
c
 285   call inpi(ns3)
       call inpa(ztext)
       go to 10
c
c   para,node,dire
c   noecho
c
c  -- all null commands for the serial code

 290   continue
       go to 10
c
c    setio
c
 120   call inpi(nbuff)
       call inpi(ibuff)
      go to 10
c
c dumpcore
c      
 245  odumpc = .true.
      goto 10
c
c     =change/switch/file=
c
 40   if(jrec.ge.jump)go to 10
      call inpa(ztext)
c
c...  fortran sequential/vbs
c
      j=locatc(zdef,maxfrt,ztext)
      if(j.eq.0) then
         j=locatc(yft,maxfrt,ztext(1:4))
         if (j.eq.0) then
            if (ztext.eq.'punch') j = 58
         end if
         if(j.eq.0) go to 150
      endif
      iposf(j)=1
      call inpanpc(zftfil(j))
      if(zftstat(j).eq.zold)istatf(j)=2
 800  call inpa(ztext)
      if(ztext.eq.zblank)go to 10
      if(ztext.eq.zkeep) then
         keepf(j)=1
         go to 800
      endif
      if(ztext(1:4).eq.'form') then
         oform(j)=.true.
         go to 800
      endif
      if(ztext.ne.zlength)go to 800
      call inpi(nrecf(j))
      call inpi(lrecf(j))
      if(nrecf(j).le.0.or.nrecf(j).gt.999999)call caserr
     *('invalid length for fortran file')
      go to 800
c
c...  512 word files (direct access)
c
 150  j=locatc(yed,maxlfn,ztext(1:4))
      if(j.eq.0) call caserr(
     *'invalid lfn specified in file directive')
      call inpanpc(zedfil(j))
      keep(j)=1
 80   call inpa(ztext)
      if(ztext.eq.zblank)go to 10
      if(ztext.eq.zkeep) keep(j)=1
      if(ztext.eq.zdel)keep(j)=0
      if(ztext.ne.zlength)go to 80
      call inpi(len(j))
      if(len(j).le.0.or.len(j).gt.9999999)call caserr
     *('invalid length for new file')
      go to 80
c
c     =store/core/memo=
c
 50   call inpi(llll)
      if(jump.gt.2) then
 51    call inpa(ztext)
       if (ztext(1:4).eq.'mwor') then
          llll = llll*1000*1000
       else if (ztext(1:4).eq.'mbyt') then
          llll = llll*1000*1000/8
       else if (ztext(1:4).eq.'gwor') then
          llll = llll*1000*1000*1000
       else if (ztext(1:4).eq.'gbyt') then
          llll = llll*1000*1000*1000/8
       else if (ztext(1:4).eq.'tota') then
          nn = ipg_nnodes()
          if (nn.le.0) nn = 1
          llll = llll/nn
       endif
       if(ztext(1:4).eq.'dump') then
        omem(4)=.true.
       else if(ztext(1:4).eq.'scra') then
        omem(8)=.true.
       else if(ztext(1:4).eq.'inte') then
        omem(3)=.true.
       else if(ztext(1:4).eq.'prin') then
        oprintm = .true.
       else if(ztext(1:4).eq.'debu') then
        call gmem_set_debug(.true.)
       else if(ztext(1:4).eq.'norm') then
        oprintm = .false.
        call gmem_set_debug(.false.)
        call gmem_set_quiet(.false.)
       else if(ztext(1:4).eq.'quie') then
        oprintm = .false.
        call gmem_set_quiet(.true.)
       else if(ztext(1:4).eq.'nona') then
        call gmem_set_nanify(.false.)
       else if(ztext(1:4).eq.'nani') then
        call gmem_set_nanify(.true.)
       else if(ztext.eq.zblank) then
        go to 90
       else
        go to 51
       endif
       go to 51
      endif
      go to 90
c
c     =lpage=
c
 60   call inpi(llll)
      llll=llll*65536
 90   lword=llll
      omemin = .true.
      go to 10
 70   call inpf(ffmin)
      insec=nint(ffmin*60.0d0)
      otime=.true.
      go to 10
c
 30   jrec=0
      if(mfree.le.0) mfree = 40000
      if(isecs.le.0)  isecs = 3600000
      if(otime)      isecs=insec
      if(nsz.le.0.or.nsz.gt.10)nsz=10
      if(ns2.le.0.or.ns2.gt.1)ns2=1
      if(ns3.le.0.or.ns3.gt.1)ns3=1
c
       call forstr
c
      if(opun7)write(iwr,201)zinout(2)
201   format(/
     *' punchfile allocated : lfn = ',a8)
c
c
c  overwrite default 4 MW core for  semi-direct mrdci code
c  hopefully this is a short term requirement.
c  Its caused by the diag requirements with increased mxrootn
c
      if(opark.and. .not.omemin) then
       lword = 20000000
       opark = .false.
      endif
c
c  for i8 code, overwrite default 4 MW core for 
c  solvation code
c  hopefully this is a short term requirement.
c  Its caused by the absurd default integer requirements ..
c
      if(oreact.and. .not.omemin) then
       lword = 8000000
       oreact = .false.
      endif
c
      return
      end
c
c** generic forstr (serial)
c
      subroutine forstr
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      character *132 filatt
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/disc/isel(5+maxlfn*8),
     *ipos(maxfrt),keep(maxfrt),istatf(maxfrt),
     *lrecf(maxfrt),nrecf(maxfrt),oform(maxfrt)
c
      character *44 zcfsdr
      integer lcfsdr
      common/zcfsinfo/zcfsdr
      common/icfsinfo/lcfsdr
      logical ohost
      common/ipsc/ohost(maxfrt)
c
      do 1 i=1,maxfrt
      if(ipos(i).le.0.and.keep(i).le.0) go to 1
      filatt=' '
      iunit=i
      call gtnv(zftfil(iunit),filatt)
      if (filatt.eq.' ') filatt = zftfil(iunit)
      do 10 loop = len(filatt),1,-1
        if(filatt(loop:loop).ne.' ') goto 20
 10   continue
c
 20   if(oform(i)) then
      open(unit=iunit,iostat=ioerr,file=filatt,status='unknown',
     *form='formatted')
      else
      open(unit=iunit,iostat=ioerr,file=filatt,status='unknown',
     *form='unformatted')
      endif
      if(ioerr.ne.0)then
      if(ooprnt) write(iwr,3)yft(i),filatt
3     format(1x,'error opening fortran file',1x,a4,
     * ' : filename : ',a132)
      endif
      if(ooprnt) write(iwr,2)yft(i),filatt
 2    format(/' fortran stream ',a4,' allocated: file name = ',a)
 1    continue
      end
c
c** generic initb (serial)
c
      subroutine initb
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      character *44 zcfsdr
      integer lcfsdr
      common/zcfsinfo/zcfsdr
      common/icfsinfo/lcfsdr
c
      odd=.false.
c
c  print CFS directory
c
      if(lcfsdr .ne. 0)then
         write(iwr,301)zcfsdr
 301     format(//1x,'cfs directory = ',a44//)
      endif
c
c  ==== check for non-default-names edn streams ======
c
      do 100 i=1,maxlfn
         zied=yed(i)
         call gtnv(zedfil(i),filatt)
         if(zied.ne.zedfil(i).or.filatt.ne.' ')  go to 200
 100  continue
      go to 600
c
c  ==== print ed stream file assignments ====
c
 200  write(iwr,300)
 300  format(//1x,
     *'local     ',2x,'status',2x,'external file name'/
     *1x,'file name '/1x,80('*'))
      odd=.true.
      do 400 i=1,maxlfn
      zied=yed(i)
      filatt=zedfil(i)
      if(zied.eq.zedfil(i))then
      call gtnv(zedfil(i),filatt)
      if(filatt.eq.' ')  go to 400
      endif
      inquire (file=filatt,exist=ok)
      if (ok) then
         ztemp = 'old'
      else
         ztemp = 'new'
      endif
      write(iwr,500)yed(i),ztemp,filatt
 400  continue
 500  format(1x,a4,8x,a7,2x,a80)
c
c  ==== check for non-default-named fortran strams
c
 600  do 700 i=1,maxfrt
      zied=zftn(i)
      call gtnv(zied,filatt)
      if(zied.ne.zftfil(i).or.filatt.ne.' ')  go to 800
 700  continue
      return
c
c  ====  print fortran file assignments
c
 800  if(odd) then
      write(iwr,850)
 850  format(/)
      else
      write(iwr,300)
      endif
      do 900 i=1,maxfrt
      zied=zftn(i)
      filatt=zftfil(i)
      if(zied.eq.zftfil(i))then
        call gtnv(zftn(i),filatt)
        if(filatt.eq.' ') go to 900
      endif
      write(iwr,500)yft(i),zftstat(i),filatt
 900  continue
      return
      end
c
c** generic ed2mem (serial)
c
      subroutine ed2mem(nbasis,nat)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/vinteg/q(1)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
c     allocate memory for 2-electron list if available
c
      iqaddr = 0
      iqqoff(3) = 0
      if(maxset(3).gt.0) then
       mxblock = maxset(3)
      else
       nsq = nbasis*(nbasis+1)/2
       n4 = nsq * (nsq +1 )/2
       mxblock = (n4-1)/num2e + 1
c .... allow for j+k supermatrix
       if(nopk.eq.-1) mxblock = mxblock + mxblock
c .... allow for symmetry skeletonisation
       if (nat.gt.3.and.nt.gt.1) then
         if (nat.gt.5.and.nopk.eq.0) mxblock = mxblock + mxblock
         nttt = nt / 2
         mxblock = mxblock / nttt
       endif
      endif
      need = mxblock * 512
      call getmem(need, q, iqaddr, iqqoff(3))
c
      if (iqaddr.eq.0.or.iqqoff(3).eq.0) then
       if(n2file.eq.1) then
        write(iwr,300) need
        call caserr('no memory available for ed2')
        omem(3)=.false.
        go to 200
       else
        ii = 1
100     mxblock = mxblock / 2
        need = mxblock * 512
        call getmem(need, q, iqaddr, iqqoff(3))
        if (iqaddr.eq.0.or.iqqoff(3).eq.0) then
         ii = ii + 1
         if (ii.le.10.and.mxblock.gt.0) go to 100
         omem(3) = .false.
         go to 200
        endif
       endif
      endif
      write(iwr,400) mxblock,need,mxblock*4/1000.0d0
      omem(3)=.true.
      maxbl(3) = mxblock
      n2last(1)= mxblock
c
200   return
300   format(1x,'requested memory not available (',i12,' words')
400   format(1x,'sufficient memory allocated for an mfile of',i10,1x,
     +'integral blocks',/,' (',i12,' words / ',f11.3,' Mbyte )')
      end
c
c** generic ed0mem (serial)
c
      subroutine ed0mem(iunit)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      common/vinteg/q(1)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
c     allocate memory for dump or scratchfile if available
c
      iqaddr = 0
      iqqoff(iunit) = 0
      if(maxset(iunit).le.0) call caserr(
     + 'failure to allocate dump- or scratch-file memory')
      mxblock = maxset(iunit)
      need = mxblock * 512
      call getmem(need, q, iqaddr, iqqoff(iunit))
      if (iqaddr.eq.0.or.iqqoff(iunit).eq.0) then
       write(iwr,*)' requested memory not available (',need,' words)',
     1             ' so unit ',iunit,' will use disk'
       omem(iunit)=.false.
      else
        write(iwr,*)' sufficient memory allocated for ',maxset(iunit),
     + ' 512 word blocks (unit =',iunit,')'
       omem(iunit)=.true.
      endif
      return
      end
c
c** generic put (serial)
c
      subroutine put(c,nword,iunit)
      implicit real*8  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension c(*)
c ... atmol  i/o
      common/bufa/abc(512)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
      common/vinteg/q(1)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
      isel=iunit
      irec=ipos(isel)
      if(irec.lt.1) call iofail(iunit,2)
      ipos(isel)=irec+1
      if (ooprnt) then
       write(iwr,10) nword, iunit, irec
 10    format(1x,'**** put ', i4,' words to unit', i3,
     +           ' at block',i5)
      endif
      if(omem(iunit)) then
        if (irec.gt.maxbl(iunit)) then
c
c         we have exceeded the memory available for this file
c
          write(iwr,*)'*** attempting to write block ',irec,
     +                ' but have memory for only ',maxbl(iunit),
     +                ' blocks'
          call iofail(iunit,2)
        endif
        iword = (irec-1)*512+iqqoff(iunit)+1
        q(iword+511)=pack2(irec,nword)
        call dcopy(nword,c,1,q(iword),1)
      else
        call dcopy(nword,c,1,abc,1)
        abc(512)=pack2(irec,nword)
        call putcc(nam(isel),abc,ierrio)
        if(ierrio.ne.0)call ioerr('write',isel,' ')
      endif
      return
      end
c
c** generic rdtpnm (serial)
c
      subroutine rdtpnm(iunit)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      logical tell
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,
     * ipos(maxlfn),nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),
     * len(maxlfn*4)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      filatt = zedfil(iunit)
      isel=iunit
      if(yed(isel).eq.zedfil(isel)) then
       call gtnv(zedfil(isel),filatt)
       if(filatt.eq.' ') filatt = zedfil(isel)
      endif
      call strtrm(filatt,lencha)
      length = 0
c
c     open file
c
      if (iunit.ne.3.and.omem(iunit)) then
       call ed0mem(iunit)
      endif
      if(omem(iunit)) return
      inquire (file=filatt,exist=tell)
      if(tell) then
         Call length_file( filatt, length, irep )
         istat(isel)=2
         if(irep.ne.0)call iofail(isel,3)
      endif
      call opencc(filatt,lencha,nam(isel),ierrio)
      if(ierrio.ne.0)call ioerr('open',isel,filatt)
      if(ooprnt)write(iwr,1)yed(isel),length
  1    format(/' lfn ',a4,' allocated: current length =',
     * i5,' blocks')
      return
      end
c
c** generic delfil (serial)
c
      subroutine delfil(iunit)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock

      if(ipos(iunit).gt.0)then
       if(.not.omem(iunit)) then
          isel=iunit
          call clredx
           if(keep(iunit).gt.0) then
             call closecc(nam(iunit))
             istat(iunit)=2
           else
             call closecc(nam(iunit))
             call strtrm(zedfil(iunit),length)
             call delcc(zedfil(iunit),length,ierrio)
            if(ierrio.ne.0)call ioerr('delete',iunit,zedfil(iunit))
             istat(iunit)=0
           endif
          if(ooprnt)write(iwr,3)yed(iunit)
 3      format(//' file ',a4,' closed')
       endif
       ipos(iunit)=-1
      endif
      return
      end
c
c** generic shut1 (serial)
c
      subroutine shut1(i)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *            nam(maxlfn),
     *  keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      if(ipos(i).gt.0)then
        ipos(i)=0
        istat(i)=2
        if(.not.omem(i)) then
          call clredx
          call closecc(nam(i))
          if(ooprnt)write(iwr,3)yed(i)
 3        format(//' file ',a4,' closed')
        endif
      endif
      return
      end

c****************************************************************
c************ bitmanipulation routines (i8 machines)  ***********
c************ temp. fix to check out base modules     ***********
c************ note that most of these functions are defunct *****
c****************************************************************
      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines 
c      (M4 flags were opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8)
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
      integer*4 itable(15)
      parameter (zero=0,mask=15)
      equivalence(a,ii(1))
      data itable/3,2,2,1,1,1,1,8*0/
      data mm4,mm8/-4,-8/
c
      a  =  b
      leadz = 64
c
      do 100 i=1,4
         iii = ii(5-i)
         if (iii.ne.zero) then
c...  we found a non-zero 2-byte number so split it
            jj = ishft(iii,mm8)
            if (jj.ne.zero) then
               koff = 0
            else
               jj = iii
               koff = 8
            end if
c...  and split it again
            kk = and(ishft(jj,mm4),mask)
            if (kk.ne.zero) then
               leadz = (i-1)*16 + itable(kk)  + koff
            else
               leadz = (i-1)*16 + itable(jj) + koff + 4
            end if
            return
         end if
100   continue
c
      return
      end


      integer function popcnt(b)
c
c...   64-bit popcount
c...   cray-imitation for direct for 32 bit machines
c...   M4 flags user were opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      integer*4 itable(0:15),ii(2),j,jj,m1,m3,m4
      integer*4 mask
      parameter (mask=15)
      equivalence(a,ii(1))
      data itable/ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4/
      data m1,m3,m4/1,3,4/
c
      a  =  b
      kk = 0
c

      do 100 i=1,2
         do 10 j=1,8
            jj = and(ishft(ii(m3-i),-(j-m1)*m4),mask)
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk

c
      return
      end
c
c** cio/fcio ioerr
c
      subroutine ioerr(oper,iunit,file)
      character oper*(*), file*(*)
      character *80 errstr,temp
      logical ocaser
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      character mess*100
      common/syscod/icode1,mess
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      common/iotrp/ocaser
c
c  report i/o problems from cio (or fcio)
c  routines, and exit
c
      if(file.ne.' ')then
         temp = '; file = '//file
      else
         temp = '; logical file ='//yed(iunit)
      endif
      if(iunit.ne.0)then
         write(iwr,1000)oper,iunit,temp
      else
         write(iwr,1001)oper,temp
      endif
      if(icode1.ne.0)write(iwr,1002)mess(1:30)
      write(iwr,1003)
      errstr='i/o error'//temp

      call gamerr(errstr,
     &     ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_SYS)

 1000 format(//,1x,'*** i/o error during ',a8,1x,'; stream =',
     &     i3,1x,a80)
 1001 format(//,1x,'*** i/o error during ',a8,1x,a80)
 1002 format(1x,'*** system message is ',a30)
 1003 format(//)
      end
c
c** cio iolog
c
      subroutine iolog(icode,mes)
      character mes*100
      character mess*100
      common/syscod/icode1,mess
c
c  to save any useful error info from system
c  (called from c i/o routines)
c
      icode1=icode
      mess = mes
      return
      end
c
c** non-ipsc dclock
c
      function dclock()
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c use elapsed times here on non-ipsc 
      call walltime(dum)
      dclock = dum
      return
      end
      subroutine ver_machscf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/machscf.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end

      Subroutine length_file( file, blocks, error )
*
*  This routine returns the length of the named file in blocks
* Highly system dependent !
*
      Implicit none

      Character*(*) file
      Integer       blocks
      Integer       error

      Integer       length_name
      Integer       nbytes

      length_name = Len( file )
      Call statcc( file, length_name, nbytes )
      If( nbytes .EQ. -1 ) Then
         Call caserr('statcc error')
      End If
      blocks = ( nbytes / 8 - 1 ) / 512 + 1
      error  = 0

      End

c
c     Revised error handling
c     
      subroutine gamerr(desc, code, class, sync, sys)
      implicit none
      integer code, class, sync, sys, code2
      character*(*) desc

c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical opg_root
      external opg_root

      integer ipg_nodeid
      external ipg_nodeid

      integer ipg_nnodes
      external ipg_nnodes

      integer length
      character pad*80
      character desc2*80
      
      logical odumpc
      common/crdmp/odumpc

      logical ocaser
      common/iotrp/ocaser
c
c     
c     check/set flag to avoid recursive calls
c     (typically from i/o activity arising from shut call
c     below)
c
      if(ocaser)return
      ocaser = .true.
c     
c     Print out last input record
c     
      if (jrec.ge.0) then
         write (iwr,10)
 10      format(' current input line is:')
         call outrec
      endif

      call strtrm(desc,length)
      desc2 = desc
      if(length.gt.80)length=80
      pad = ' '


      if(code .ne. ERR_NO_CODE)then
          code2 = code
      else
          code2 = -1
      endif
c     
c     Output the message to stdout
c     syncronous error processing - node 0 writes the caserr banner
c     we also write this for async errors if we hit them on
c     node 0 (which covers the serial case)
c     
      if(opg_root())then
         if(sync .eq. ERR_SYNC)then
            if(code .eq. ERR_NO_CODE)then
               write(6,30)desc2(1:length)//pad(1:80-length)
 30            format(//65('*#')//50x,'error detected'/50x,
     &              '=************=',/
     &              /20x,a80//65('*#')//)
            else
               write(6,31)code,desc2(1:length)//pad(1:80-length)
 31            format(//65('*#')//50x,'error ',i5,' detected '/50x,
     &              '=*******************=',/
     &              /20x,a80//65('*#')//)

            endif
         else
            if(code .eq. ERR_NO_CODE)then
               write(6,33)desc2(1:length)//pad(1:80-length)
 33            format(//65('*#')//50x,'error detected'/50x,
     &              '**************',/
     &              /20x,a80//65('*#')//)
            else
               write(6,34)code,desc2(1:length)//pad(1:80-length)
 34            format(//65('*#')//50x,'error ',i5,' detected '/50x,
     &              '*********************',/
     &              /20x,a80//65('*#')//)
            endif
         endif
         if(class .eq. ERR_INCOMPREHENSIBLE)then
           write(6,*)'A directive or keyword was not recognised or '
           write(6,*)'was out of order. Fix input errors and try again'
         else if(class .eq. ERR_INCOMPATIBLE)then
            write(6,*)
     &     'You have provided an incompatible set of control directives'
            write(6,*)'Fix input errors and try again'
         else if(class .eq. ERR_DIMENSION)then
            write(6,*)'You will need to redimension the source code'
            write(6,*)'Edit the file common/sizes and recompile, or'
            write(6,*)'contact gamess_uk_support@dl.ac.uk'
         else if(class .eq. ERR_USER_DIMENSION)then
            write(6,*)'Modify the input file and try again'
         else if(class .eq. ERR_FIXED_DIMENSION)then
            write(6,*)'You have encountered an internal size limit'
            write(6,*)'Please contact gamess_uk_support@dl.ac.uk'
         else if(class .eq. ERR_INTERNAL)then
            write(6,*)'GAMESS_UK internal error'
            write(6,*)'Please contact gamess_uk_support@dl.ac.uk'
         else if(class .eq. ERR_SWITCHED_OUT)then
           write(6,*)'You requested an option that was not enabled when'
           write(6,*)'GAMESS-UK was built. If you have the source code'
           write(6,*)'distribution, repeat the configure command, and '
           write(6,*)'then make clean and recompile.'
         else if(class .eq. ERR_UNLUCKY)then
c
c nothing to say here really
c
         else if(class .ne. ERR_NO_CLASS)then
            write(6,*)'ERROR CALL WITH INVALID ERROR CLASS'
         endif
c     flush sequential fortran stream on unix machines
         call flushn(6)
      endif

      if( sync .eq. ERR_ASYNC .and. .not. opg_root())then
c     
c     the message has been encountered on a slave node
c     the node should output a short form of the message
c     to the log file
c     (this may or may not be visible to the user)
c     
         if(code.eq.ERR_NO_CODE)then
            write(6,100)ipg_nodeid(),desc2(1:length)//pad(1:80-length)
 100        format(' **** GAMESS-UK Error on node ',i3,': ',a80)
         else
            write(6,101)code,ipg_nodeid(),
     &           desc2(1:length)//pad(1:80-length)
 101        format(' **** GAMESS-UK Error ',i4,' on node ',i3,': ',a80)
         endif
      endif
c     
c     stderr logging (includes system error message if requested)
c     sync:  only on node 0 
c     async: independent of the node number
c     
      if(opg_root() .or. sync .eq. ERR_ASYNC)then
         call stderrc(desc2,length,code,sys,ipg_nodeid())
      endif

      call flushn(6)

      if(sync .eq. ERR_SYNC)then
c     
c     below this point, the OS may kill other nodes in a parallel
c     job, for synchronous errors wait till all output has been 
c     completed
c     
         call pg_synch(99)
      endif
c
      if(sync .eq. ERR_SYNC .or. ipg_nnodes() .eq. 1)then
c     
c     attempt to shut down I/O system
c
         call whtps
         call shut
      endif

c     
c  attempt a core dump if requested. At present this is only
c  available on the root node with the C I/O system
c
      if(odumpc)then
         if(opg_root())then
            write(iwr,*)'attempting to dump core on root node error=',
     &           desc
            call cordmp
         endif
         call exitc(-1)
      endif
c     
c     exit
c     
      call clenup2
c
      if(sync .eq. ERR_SYNC .or. ipg_nnodes() .eq. 1)then
c
c clean exit will help ensure that the output is complete
c

         call pg_end(code2)
c
      else 
c
c Async errors (may occur on a subset of nodes)
c
         call pg_err(code2)
      endif
      end
c=========================================================================
c
c
c also need a 16 bit version for leadz split algorithm
c
      function ishft2(in2,bits2)
      implicit none
      integer*2 ishft2
      integer*2 in2, bits2
      integer in,out,bits
      integer*2 buff(2)
      equivalence (buff(1),in)

      bits = bits2
c      call pr1s('inp ',in2)
     
      buff(1)=0
      buff(2)=in2

c      call pr1i('int i',in)

      call ishftcc(in,out,bits)

c      call pr1i('int o',out)

      in = out
      ishft2 = buff(2)
      end
      integer function ibset(in,bit)
      implicit none
      integer in, out, bit
      call ibsetcc(in,out,bit)
      ibset = out
      end

      integer function ibclr(in,bit)
      implicit none
      integer in, out, bit
      call ibclrcc(in,out,bit)
      ibclr = out
      end

      logical function btest(in,bit)
      implicit none
      integer in, out, bit
      call btestcc(in,out,bit)
      btest = ( out .ne. 0)
      end
c
c** generic ioupdate
c
c  route selected files (GA and in eventually in memory)
c  back to disk
c
      subroutine ioupdate

      implicit none

      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
************************************************************************
*
*  Include file for the files mapped into GA tools
*
*
*  Parameter for the file blocking. Chosen to be consistent
* with the disk based system.
*
      Integer    ga_file_block
      Parameter( ga_file_block = 511 )
*
*  Maximum number of files mapped into GA tools.
*
      Integer    max_ga_files
      Parameter( max_ga_files = 2 )
*
*  The maximum number of writes performed to a GA file 
* This is required so that the file  can be dumped
* to disk if required.
*
      Integer    max_ga_calls
      Parameter( max_ga_calls = 60000 / max_ga_files )
*
*  The largest unit number that the code will use in a call
* to the IO subystem.
*
      Integer    max_unit_number
      Parameter( max_unit_number = 99 )
*
*  The smallest unit number that the code will use in a call
* to the IO subystem.
*
      Integer    min_unit_number
      Parameter( min_unit_number = 0 )
*
*  A flag value to indicate when something is not being used.
*
      Integer    not_in_use
      Parameter( not_in_use = -999 )
*
*  Couple of Parameters to make the mode variable more comprehensible.
*
      Logical    SINGLE_NODE          , COLLECTIVE
      Parameter( SINGLE_NODE = .False., COLLECTIVE = .True. )
*
*  Flag to tell if GA tools has been initialized.
*
      Logical ga_file_initialized
*
*  The next three integers set limits on how much memory the
* GA file system will use. These are NOT parameters to allow the
* user to tune them to their requirements. However there is a default
* value for the amount of memory per processor. This is set in the
* Block Data ga_file_initialize.
*
*  The amount of memory ( in terms of 64 bit quantities )
* that we are going to allow each processor to use for the virtual
* file system. 
*
      Integer ga_file_proc_mem
*
*  The total amount of memory across the whole machine
*  that may be used.
*
      Integer ga_file_total_mem
*
*  The total space used at present by the virtual file system
* and how much left we have got.
*
      Integer ga_file_total, ga_file_unused
*
*  The total number of processes used by the GA tools
*
      Integer ga_file_processes
*
*  The id of this process
*
      Integer ga_file_proc_id
*
*  Flag to indicate if this is the root node. This is used to ensure
* coherency during writes
*
      Logical ga_file_root
*
*  The descriptor for each file. Oh for derived data types !
* But anyway they are:
* GA_FILE_HANDLE   :  The GA handle for the array
* GA_FILE_OPEN     :  Flag to tell if the file is open
* GA_FILE_SIZE     :  How much memory used internally by GA tools
*                     for this file. 
* GA_FILE_NEXT     :  Next block that will be written to in this file
* GA_FILE_MODE     :  Looks after file consistency when  OSWED3
*                     in GAMESS's nodeio common block changes.
* GA_FILE_KEEP     :  If set dump file to disk at the end of the job.
* GA_FILE_DUMP_MODE:  Controls whether only the node0 dumps its version
*                     of the file or all nodes dump versions.
* GA_FILE_ACCESSED :  If false indicates that this is the first read
*                     from write to a file. Used in determining if
*                     this isa restart file.
* GA_FILE_LENG1    :  Length of the first block. Used for restarts.
*
      Integer      ga_file_handle   ( 1:max_ga_files )
      Logical      ga_file_open     ( 1:max_ga_files )
      Integer      ga_file_size     ( 1:max_ga_files )
      Integer      ga_file_next     ( 1:max_ga_files )
      Logical      ga_file_mode     ( 1:max_ga_files )
      Logical      ga_file_keep     ( 1:max_ga_files )
      Integer      ga_file_dump_mode( 1:max_ga_files )
      Logical      ga_file_accessed ( 1:max_ga_files )
      Integer      ga_file_leng1    ( 1:max_ga_files )
*
*  The next two look after the history of the ga_file. This
* is used if the file needs to be dumped to disk.
* GA_FILE_HISTORY  :  The history of writes to the files. These are
*                     required so that the
*                     file can be dumped to disk if required.
* GA_FILE_HIST_NEXT:  The next entry in the history array
*
      Integer ga_file_history  ( 1:3, 1:max_ga_calls, 1:max_ga_files )
      Integer ga_file_hist_next( 1:max_ga_files )
*
*  An array to map from the unit number used by the code to the
* unit number used internally by these routines.
*
      Integer ga_file_unit_to_ga( min_unit_number:max_unit_number )
*
*  The common blocks for the above
*      
      Common / ga_file_status     / ga_file_initialized 
      Common / ga_file_mem_limits / ga_file_total_mem,
     +                              ga_file_proc_mem
      Common / ga_file_mem_used   / ga_file_total, ga_file_unused
      Common / ga_file_parallel   / ga_file_processes,
     +                              ga_file_proc_id, 
     +                              ga_file_root
      Common / ga_file_descriptor / ga_file_handle,
     +                              ga_file_open,
     +                              ga_file_size,
     +                              ga_file_next,
     +                              ga_file_mode,
     +                              ga_file_keep,
     +                              ga_file_dump_mode,
     +                              ga_file_accessed,
     +                              ga_file_leng1
      Common / ga_file_lifestory  / ga_file_history,
     +                              ga_file_hist_next
      Common / ga_file_units      / ga_file_unit_to_ga
*
*  Some buffers used for storing integers in double precisions and a 
* horrible equivalence. The length of this buffer is chosen to give 
* decent transfers rates on the T3D. NOTE THE LENGTH MUST BE A
* MULTIPLE 
* OF GA_FILE_BLOCK to ensure internal consistency of the files.
*
      Integer    ga_file_buf_len
      Parameter( ga_file_buf_len = 32 * ga_file_block )
*
      Integer          ga_file_integer_buf( 1:ga_file_buf_len )
      Double Precision ga_file_double_buf ( 1:ga_file_buf_len )
*
      Equivalence ( ga_file_integer_buf, ga_file_double_buf )
*
*  And stick em in common so that one space is globally accessible,
* rather than lots of different buffers building up. Use the double
* precision version as it is likely to be the longer.
*
      Common / ga_file_horror / ga_file_double_buf
*
*  This works out the biggest number of disk blocks we can
* fit in the internal buffer. Note there is a potential
* bug here if ga_file_buf_len > ga_file_block ** 2.
*
      Integer    blocks_wanted
      Parameter( blocks_wanted = ga_file_buf_len / 
     +                         ( ga_file_block + 1 ) )
*
*  The maximum number of processors is required when trying to find
* out if this processor has some part of the global array locally. 
*
      Integer    max_processes
      Parameter( max_processes = mxproc )
*
*  A couple of parameters and a safety value used internally to avoid
* excessive calls to the sychronization routine.
*
      Integer    ga_file_safe    , ga_file_risky
      Parameter( ga_file_safe = 1, ga_file_risky = 2 )
*
      Integer ga_file_safety_mode
*
      Common / ga_file_safety / ga_file_safety_mode
*
*  This is required if the GA files are accessed by the low level
* GET/PUT calls. See comments in the routines search_ga, put_ga 
* and get_ga about this.
*
      Integer ga_file_last_unit
      Common / ga_file_low_level / ga_file_last_unit
*
***********************************************************************


      end
      subroutine inclin(rec132,ifile,iw,o_inclin,inc_err)
c
c...  process an include or %include statement (lower or upper case)
c...  routine now returns next line to write 
c...  o_inclin tell's if we are currently including
c     rec132 : current card
c     ifile : number of extra unit to use
c     iw    : number of output (print) unit
c     o_inclin : true => we are including
c     o_inc_err : used for flagging errors
c
      implicit none
      character*132 rec132,filen
      integer i,if,il,iw,iread,ifile,linp,iftchr,lstchr
      logical oprint,oex,o_inclin,inc_err
      save filen,oprint,linp
c
      inc_err=.false.
      if (o_inclin) go to 10
c
      i = iftchr(rec132)
      if (i.gt.73) return
c...   check on include  consider ? a comment
      if (rec132(i:i).eq.'?') return
      if (rec132(i:i).eq.'%') i = i + 1
      if (rec132(i:i+6).ne.'include'.and.rec132(i:i+5).ne.'INCLUD')
     1   return 
c
      i = i+7
      if = iftchr(rec132(i:)) + i-1
      il = index(rec132(if:),' ') + if-2
      linp = il-if+1
      filen = rec132(if:il)
c
      inquire(file=filen(1:linp),exist=oex)
      if (.not.oex) then
         write(iw,63)filen(1:linp)
63       format(1x,/,' // trying to include ',a,' \\')
         inc_err = .true.
      call caserr('include file not found')
      end if 
      open(unit=ifile,form='formatted',file=filen(1:linp),
     +     status='unknown')
      o_inclin = .true.
c
      write(iw,62) filen(1:linp)
62    format(1x,/,' // file ',a,' included  \\\\')
c
c..   check for print or PRINT on line
c
      if = iftchr(rec132(il+1:)) + il
      if (if.gt.76) then
         oprint = .false.
      else if (rec132(if:if+4).eq.'print'.or.
     1         rec132(if:if+4).eq.'PRINT') then
            oprint = .true.
         else
         oprint = .false.
      end if
      rec132 = '? begin  file - '//filen(1:linp)
      return
c
c..    copy lines
c
10    read(ifile,1080,end=20) rec132
 1080 format(a132)
      if (oprint) write(iw,61) filen(1:linp),rec132(1:lstchr(rec132))
61    format(' -- ',a,' : ',a)
c
c...   no includes in includes are permited now (recursive)
c
c
      return
c
20    close(ifile)
      o_inclin = .false. 
      rec132 = '? end of file - '//filen(1:linp)
c
      return
      end
      integer function iftchr(a)
      character*(*) a
      integer n,len,i
c routine to return index of first non-blank character in character var
c cf. lstchr(a)
      n=len(a)
      do 10 i=1,n
      if( a(i:i).ne.' ' ) go to 20
10    continue
      i=n+1
20    iftchr=i
      return
      end

cjmht rewftn previously in dirctb.m
      subroutine rewftn(ntape)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer error_code
      character *132 filatt
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     * nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn),
     * iwhere(maxlfn),isync(maxlfn),
     * iposf(maxfrt),keepf(maxfrt),istatf(3,maxfrt),
     * oform(maxfrt)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
c     write (6,*) 'rewind fortran stream', ntape
c
      close(unit=ntape, iostat=error_code, status='keep')
      if(error_code.ne.0)then
         write(6,5)zftn(ntape)
 5      format(1x,'***** error closing file',1x,a8)
      endif
c
      filatt=' '
c      call gtnv(zftn(ntape),filatt)
c      if (filatt.eq.' ') filatt = zftn(ntape)
cjmht psh change to get file routing to work with file directive
      call gtnv(zftfil(ntape),filatt)
      if (filatt.eq.' ') filatt = zftfil(ntape)

      do 10 loop = len(filatt),1,-1
        if(filatt(loop:loop).ne.' ') goto 20
 10   continue
 20   if(oform(ntape)) then
       open(unit=ntape,iostat=error_code,
     +      form='formatted',file=filatt,status='unknown')
      else
       open(unit=ntape,iostat=error_code,
     +      form='unformatted',file=filatt,status='unknown')
      endif
      if(error_code.ne.0)then
      write(6,3)yft(ntape),filatt
3     format(1x,'error opening fortran file',1x,a4,
     * ' : filename : ',a132)
      endif
      return
      end
c
c** generic prinit (serial,parallel and cray)
c
      subroutine prinit(oreact,opark,odci,odlfi)
c
c     this routine pre-scans the input file, checking for major
c     inconsistencies in data and establishing requested
c     runtypes to pre-allocate fortran data sets
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      logical okillv
      integer isecvv, itypvv
      common /vectrn/ isecvv,itypvv,okillv
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn),keep(maxlfn),istat(maxlfn),llll(maxlfn),
     * iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     * iposf(maxfrt),keepf(maxfrt),istatf(maxfrt),
     * lrecf(maxfrt),nrecf(maxfrt),oform(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c The oblock array determines what is written out to the punchfile:
c
c Below sourced from server.m
c oblock(1)   - coords - single point or last point
c keyword: coor
      integer PUN_COORD
      parameter (PUN_COORD=1)
c oblock(2)   - coordinates at each optimisation step
c keyword: opti
      integer PUN_OPT_COORD
      parameter (PUN_OPT_COORD=2)
c oblock(3)   - starting coordinates
c keyword: init
      integer PUN_START_COORD
      parameter (PUN_START_COORD=3)
c oblock(4)   - atom connectivity
c keyword: conn
      integer PUN_ATOM_CONN
      parameter (PUN_ATOM_CONN=4)
c oblock(5)   - molecular orbital occupations
c keyword: occu
      integer PUN_MO_OCC
      parameter (PUN_MO_OCC=5)
c oblock(6)   - eigenvectors
c keyword: vect
      integer PUN_EIGV
      parameter (PUN_EIGV=6)
c oblock(7)   - run type
c keyword: type
      integer PUN_RUNT
      parameter (PUN_RUNT=7)
c oblock(8)   - gvb pair data
c keyword: gvb
      integer PUN_GVBPAIR
      parameter (PUN_GVBPAIR=8)
c oblock(9)   - force constant matrix
c keyword: forc
      integer PUN_FCMAT
      parameter (PUN_FCMAT=9)
c oblock(10)  - vibrational frequencies
c keyword: vibr
      integer PUN_VIBFREQ
      parameter (PUN_VIBFREQ=10)
c oblock(11)  - basis set
c keyword: basi
      integer PUN_BASIS
      parameter (PUN_BASIS=11)
c oblock(12)  - orbital energies
c keyword: eige
      integer PUN_ORB_ENER
      parameter (PUN_ORB_ENER=12)
c oblock(13)  - total scf energies
c keyword: scfe
      integer PUN_SCF_ENER
      parameter (PUN_SCF_ENER=13)
c oblock(14)  - job title
c keyword: titl
      integer PUN_JOB_TITLE
      parameter (PUN_JOB_TITLE=14)
c oblock(15)  - overlap matrix
c keyword: over
      integer PUN_OVLP_MAT
      parameter (PUN_OVLP_MAT=15)
c oblock(16)  - grid data from dumpfile section(s)...
c keyword: grid
      integer PUN_GRID_DATA
      parameter (PUN_GRID_DATA=16)
c oblock(17)  - mulliken analysis
c keyword: mull
      integer PUN_MULLIKEN
      parameter (PUN_MULLIKEN=17)
c oblock(18)  - lowdin analysis
c keyword: lowd
      integer PUN_LOWDIN
      parameter (PUN_LOWDIN=18)
c oblock(19)  - spin densities
c keyword: spin
      integer PUN_SPIN_DENS
      parameter (PUN_SPIN_DENS=19)
c oblock(20)  - scf type
c keyword: leve
      integer PUN_SCF_TYPE
      parameter (PUN_SCF_TYPE=20)
c oblock(21)  - normal coordinates for vibrational modes
c keyword: norm
      integer PUN_VIB_MODES
      parameter (PUN_VIB_MODES=21)
c oblock(22)  - two electron integrals
c keyword: twoe
      integer PUN_2E_INT
      parameter (PUN_2E_INT=22)
c oblock(23)  - gradients
c keyword: grad
      integer PUN_GRADIENTS
      parameter (PUN_GRADIENTS=23)
c oblock(24)  - transformation matrix
c keyword: tran
      integer PUN_TRANS_MAT
      parameter (PUN_TRANS_MAT=24)
c oblock(25)  - potential derived charges
c keyword: pdc
      integer PUN_POT_DERV_CHG
      parameter (PUN_POT_DERV_CHG=25)
c oblock(26)  - grid symmetry array
c keyword: gsym
      integer PUN_GRID_SYMM
      parameter (PUN_GRID_SYMM=26)
c oblock(27)  - cartesian hessian matrix
c keyword: secd
      integer PUN_HESSIAN
      parameter (PUN_HESSIAN=27)
c oblock(28)  - distributed multipole analysis
c keyword: dma
      integer PUN_MPOLE_ANAL
      parameter (PUN_MPOLE_ANAL=28)
c oblock(29)  - density matrix
c keyword: dens
      integer PUN_DENS_MAT
      parameter (PUN_DENS_MAT=29)
c oblock(30)  - total energy
c keyword: ener
      integer PUN_TOT_ENER
      parameter (PUN_TOT_ENER=30)
c oblock(31)  - internal coordinate hessian matrix
c keyword: hess
      integer PUN_ZMAT_HESS
      parameter (PUN_ZMAT_HESS=31)
c oblock(32)  - dipole moment
c keyword: dipo
      integer PUN_DIPOLE
      parameter (PUN_DIPOLE=32)
c oblock(33)  - infrared intensities 
c keyword: infr
      integer PUN_INFRARED
      parameter (PUN_INFRARED=33)
c oblock(34)  - raman intensities
c keyword: rama
      integer PUN_RAMAN
      parameter (PUN_RAMAN=34)
c CURRENTLY ONLY USED FOR XML
c oblock(35)  - metadata (user, compilation date etc)
c keyword: meta
      integer PUN_METADATA
      parameter (PUN_METADATA=35)
c CURRENTLY ONLY USED FOR XML
c oblock(36)  - input parameters
c keyword: inpa
      integer PUN_INPUT_PARAM
      parameter (PUN_INPUT_PARAM=36)


      integer LENPUN
      parameter(LENPUN=60)

      logical oblock, opunold
      integer nfblck, nblsec, iblsec, nblgri, iblgri
      integer itwoe, ntwoe, iblfmt
      common/blocks/nfblck,nblsec,iblsec(10),nblgri,iblgri(10),
     +              itwoe,ntwoe,iblfmt(LENPUN),oblock(LENPUN),opunold
c

      parameter (nvec=12)
      dimension yvect(nvec)
      dimension zscf(8)
      dimension index(18)
      data index/1,2,3,4,8,9,10,
     *          22,31,33,34,35,36,41,
c     additional streams for new mrdci
     +          12,11,44,42 /
      data zscf/'rhf','uhf','gvb','grhf','casscf','mcscf','multi',
     +          'direct'/
      data yvect/'ming','extg','hcor','getq','extr','alph','card',
     +           'noge','atom','atde','ator','mcar'/
      jrec=0
      omrdci=.false.
      opark=.false.
      oreact = .false.
      odci = .false.
      oci = .false.
      numbci = 14
      ocas=.false.
      otda=.false.
      orpa=.false.
      ofulci=.false.
      omech=.false.
      omopac=.false.
      oblk=.false.
      oente = .false.
      orest = .false.
      okillv = .false.
      odlfi = .false.
      isect = 0
      go to 103
 1    call input
      if(oswit) go to 4000
 103  call inpa(ztest)
      ytest=ytrunc(ztest)
      if ((ytest.eq.'util'.or.ytest.eq.'serv')
     1     .and.ztest.ne.'servec') go to 4000
      if (ztest.eq.'servec') orest = .true.
      if(ytest.eq.'punc')then
       oblk=.true.
       go to 1
      endif
      if (ytest.eq.'rest') then
       orest = .true.
       go to 1
      endif
      if (ytest.eq.'noec') then
       oecho = .false.
       if (opg_root()) then
          write(iwr,*)'**** Echo of input suppressed'
       endif
       go to 1
      endif
      i=locatc(ydd,limit1,ytest)
      if (i.gt.0) then
         if(i.eq.62) oreact = .true.
         goto 1
      endif
      i=locatc(zdd3,limit3,ztest)
      if (i.gt.0) then
         if (i.eq.15) then
c...        found weights directive, spin until 'end'
 790        call input
            call inpa(ztest)
            if (ztest(1:3).ne.'end') goto 790
         endif
         goto 1
      endif
      i=locatc(ydd(101),limit2,ytest)
      if (i.eq.0) go to 1
c     vectors (check for vectors section specification
c     in start-up job; kill this *before* dumpfile is
c     overwritten)
      if(i.eq.84.and..not.oente) then
       call inpa4(ytest)
       ii = locatc(yvect,nvec,ytest)
       if (ii.eq.0) then
        jrec = jrec - 1
        call inpi(isect)
       endif
c .... runtype
      else if(i.eq.62) then
       call inpa4(ytest)
       jrun=locatc(yrunt,mxtask,ytest)
       if(jrun.eq.3)oci = .true.
       if(jrun.eq.12)otda=.true.
       if(jrun.eq.42) then
        call inpa(ztest)
        if (ztest.eq.'rpa') orpa=.true.
       endif
       if(jrun.eq.41)omopac=.true.
c .... scftype
      else if(i.eq.85) then
       call inpa(ztest)
       ii=locatc(zscf,8,ztest)
       if(ii.eq.5)ocas=.true.
c .... direct-ci
       else if(i.eq.71)then
       odci = .true.
c .... mrdci
      else if(i.eq.78)then
       omrdci=.true.
       if(jump.eq.2) then
        numbci = 18
        opark = .true.
       endif
c .... fullci
      else if(i.eq.79)then
       ofulci=.true.
c .... model
      else if(i.eq.80) then
       omech=.true.
      else if (i.eq.61) then
       oente = .true.
c .... dl-find
      else if (i.eq.101) then
       odlfi = .true.
c        print*,"found dl-find in first run."
      endif
      go to 1
c
 4000 oswit=.false.
      if (.not.orest.and.isect.ne.0) okillv = .true.
      jrec=0
      jump=0
      oecho = .false.
      rewind ird
      call input
      call inpa4(ytest)
      if(ytest.eq.'outp') call input
      jrec=0
      if(ocas) then
      iposf(1)=1
      iposf(2)=1
      endif
      if(omrdci) then
c     otab=.true.
      do loop=1,numbci
       iposf(index(loop))=1
      enddo
       if(numbci.eq.18) then
          oform(11)=.true.
          oform(12)=.true.
       endif
      endif
      if(otda) then
      do loop=1,8
       iposf(index(loop))=1
      enddo
      endif
      if(orpa) then
c      iposf(14) = 1
c      iposf(15) = 1
      endif
      if(omopac) then
      do loop=1,4
       iposf(index(loop))=1
      enddo
      endif
      if(ofulci) then
      do loop=1,5
       iposf(index(loop))=1
      enddo
      endif
      if(omech) then
      do loop=1,4
      iposf(index(loop))=1
       oform(index(loop))=.true.
      enddo
      endif
      if(opg_root()) then
        nfblck = 58
       if(oblk)then
          iposf(nfblck)=1
          oform(nfblck)=.true.
          keepf(nfblck)=1
        endif
      endif
c
      oci = .not.(ofulci.or.omrdci).and.oci
      if(oci) odci = .true.
      return
      end
