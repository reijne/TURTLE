
_IF(win95)
_IFN(ifc)
*$pragma aux delcc "!_" parm (value,reference,reference,reference)
*$pragma aux opencc "!_" parm (value,reference,reference,reference)
*$pragma aux getdate "!_" parm (value)
*$pragma aux iolog "!_" parm(reference,value)
*$pragma aux gtnvcc "!_" parm(value,reference,value,reference,reference)
*$pragma aux statcc "!_" parm (value,reference,reference)
*$pragma aux stderrc "!_" parm (value,reference,reference,reference,reference)
      subroutine flush(n)
      integer flushunit
      external flushunit
      integer n,i
      i=flushunit(n)
      return
      end
_ELSE
      subroutine flush(n)
      integer n,i
c     call flushcc(n,i)
      call flushout
      return
      end
_ENDIF
      subroutine exit(ia)
      character*80 ia
      close(5)
      close(6)
      stop
      end
_ENDIF
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
INCLUDE(common/iofile)
INCLUDE(common/timanb)
      REAL elapse,ccpu,cpulft
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
_IF()
_IF(convex)
      REAL  secnds
      elapse=secnds(0.0d0)
_ELSEIF(t3e)
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         call timef(t)
         tst = 0.0d0
         elapse = 0.0d0
      else
         call timef(t)
         t = t * 1.0d-3
         elapse = t - tst
      endif
_ELSEIF(t3d)
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
        ienter = 1
        t = rtc() / 150.0d6
        tst = t
        elapse = 0.0d0
      else
         t = rtc() / 150.0d6
         elapse = t - tst
      endif
_ELSEIF(bluegene)
c
c     use the C-implementation as the etime function is bust
c
      call walltimecc(elapse)
_ELSEIF(rs6000)
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         t = rtc()
         tst = t
         elapse = 0.0d0
      else
         t = rtc()
         elapse = t - tst
      endif
_ELSEIF(ipsc)
      REAL dclock
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         t = dclock()
         tst = t
         elapse = 0.0d0
      else
         t = dclock()
         elapse = t - tst
      endif
_ELSE
c  C implementation
      call walltimecc(elapse)
_ENDIF
_ENDIF
      end
_IFN(win95)
c We currently don't print the hostname on Windoze
c
c generic gethost
c
      subroutine gethost(host)
      implicit None

      character*44 host

_IF(convex,sgi,sun,dec,ksr)
      integer hostnm
_ELSEIF(vms)
      integer fgethostname
_ENDIF

_IF(convex,sgi,sun,dec,hp700,hpux11,rs6000,ksr,itanium,linux)
_IF(vms)
      idum=fgethostname(%ref(host), %val(44))
_ELSEIF(hp700,hpux11,rs6000,t3d,itanium,linux)
      call hostcc(host,44)
_ELSE
      idum=hostnm(host)
_ENDIF
_ELSE
      host = 'unknown host'
_ENDIF
      end
c     END OF GETHOST
_ENDIF(win95)
c
_IF()
_IF(rs6000)
      real function rtc()
c...   imitation of rtc in case not there, jvl 1999
      real t4(2), etime
_IF(rs6000_noextname)
c in the absence of extname we include the _ here
      rtc = etime_(t4)
_ELSE
      rtc = etime(t4)
_ENDIF
      return
      end
_ENDIF
_ENDIF
c     
c**   charwall returns walltime in hours:minutes:seconds
c
      character*10 function charwall()
c     
c.. due to compiler problems we use concatenate action
c
      implicit none
      REAL wall
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
      REAL cpu(3)
_IF(bluegene)
c
c     use C call as fortran implementations are bust
c
      call cputimecc(cpu)
_ELSEIF(rs6000)
      cpu(2) = 0.0d0
      cpu(1) = dfloat(mclock() ) / 100.0d0
      cpu(3) = cpu(1)
_ELSEIF(t3d,unicos)
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         tst=tsecnd()
         cpu(1) = 0.0
      else
         cpu(1) = tsecnd()-tst
      endif
      cpu(2) = 0.0d0
      cpu(3) = cpu(1)
_ELSEIF(ksr)
      dimension t4(2)
      save ienter,tst
      dum = etime(t4)
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         tst(1) = t4(1)
         tst(2) = t4(2)
         tst(3) = dum
         cpu(1) = 0.0d0
         cpu(2) = 0.0d0
         cpu(3) = 0.0d0
      else
         cpu(1) = t4(1)-tst(1)
         cpu(2) = t4(2)-tst(2)
         cpu(3) = dum  -tst(3)
      endif
_ELSEIF(ipsc)
      REAL tst(3)
      REAL dclock
      save ienter,tst
      data ienter/0/
      if (ienter.eq.0) then
         ienter = 1
         tst(1) = dclock()
         tst(2) = 0.0d0
         tst(3) = dclock()
         cpu(1) = 0.0d0
         cpu(2) = 0.0d0
         cpu(3) = 0.0d0
      else
         cpu(1) = dclock()-tst(1)
         cpu(2) = 0.0d0
         cpu(3) = dclock()-tst(3)
      endif
_ELSE
c use C call
      call cputimecc(cpu)
_ENDIF
      end
c
c** generic cpulft
c
      function cpulft(i)
      implicit none
      integer i
      REAL cpulft, ak
      REAL cpu(3)
      logical oex
INCLUDE(common/timez)
c
c ------ return the clock time in seconds if i=1
c ------        the cpu time remaining if i=0
c
_IF(cray)
      call ustime(ak)
_ELSE
      call gms_cputime(cpu)
      ak=cpu(1)
_ENDIF
      if(i.eq.1) then
         cpulft=ak
      else if(i.eq.0) then
         cpulft=timlim-ak
_IFN(t3d,parallel)
c
c  if file STOP exists kill GAMESS orderly
c  this inquire proves to be a fatal action on massively  machines
c  with a shared filesystem
c
         inquire(file='STOP',exist=oex)
         if (oex) cpulft = 0.5d0
_ENDIF
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
_IF(parallel)
      call pg_err(n)
_ELSE
_IF(rs6000,hp700,hpux11,linux)
c
c use C exit routine
c
      call exitc(11)
_ELSE
      call exit(11)
_ENDIF
_ENDIF
      stop
      end
_IFN(cray,ibm,vax)
c
c** generic fmove
c
      subroutine fmove(a,b,n)
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
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
      implicit REAL  (a-h,o-z)
      character*(*) text
      character*16 logfil
      character*131 logbuf
      common/clogfl/logfil,logbuf(100)
INCLUDE(common/utilc)
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
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(vms)
      character *9 zdat, ztim
      integer decc$getuid
_ELSEIF(unix)
      integer uid
_ELSE
      integer getuid
_ENDIF
      character*24 z24
_IF(convex)
      character*16 zlogg
_ELSEIF(unix,ipsc)
      character*31 zlogg
_ELSEIF(vms)
      character*31 getlog,zlogg
_ELSE
      character*24 fdate
      character*31 getlog,zlogg
_ENDIF
_IF(ipsc)
      data zblnk/' '/
_ENDIF
      if(jtime.le.0) jtime=3600000
_IF(convex)
      call fdate(z24)
      call getlog(zlogg)
_ELSEIF(win95)
      zlogg = 'GAMESS  '
      call getdate(zdate,ztime,zyear)
      zaccno = 'WINDOWS '
_ELSEIF(unix)
      call daycc(z24)
      call namcc(zlogg,31)
_ELSEIF(vms)
       call date(zdat)
       zdate(1:8) = zdat(1:8)
       zyear='19'//zdat(8:9)
       call time(ztim)
       ztime(1:8) = ztim(1:8)
       zlogg = getlog()
_ELSEIF(ipsc)
      z24=zblnk
      zlogg=zblnk
_ELSE
      z24 = fdate()
      zlogg = getlog()
_ENDIF
_IFN(vms,win95)
      ztime=z24(12:19)
      zdate=z24(5:11)
      zyear=z24(21:24)
_ENDIF
      zanam=zlogg(1:8)
c
c account number
c
_IF(vms)
      write(zaccno,'(i8)') decc$getuid(x)
_ELSEIF(win95)
_ELSEIF(unix)
      call uidcc(uid)      
      write(zaccno,'(i8)') uid
_ELSEIF(ipsc)
      zaccno = zblnk
_ELSE
      write(zaccno,'(i8)') getuid(x)
_ENDIF
      return
      end
_IF(macosx)
      subroutine namcc(zlogg,i)
      character*(*) zlogg
      call getlog(zlogg)
      end
_ENDIF
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

_IFN(t3d)
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
_IF(i8)
      parameter( mx = 8 )
_ELSE
      parameter( mx = 4 )
_ENDIF
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
_IF(i8)
_IF(littleendian)
      data nshift/56,0,8,16,24,32,40,48/
_ELSE
      data nshift/0,56,48,40,32,24,16,8/
_ENDIF
_ELSE
_IF(littleendian)
      data nshift/24,0,8,16/
_ELSE
      data nshift/0,24,16,8/
_ENDIF
_ENDIF
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
_IF(i8)
      parameter( mx = 4 )
_ELSE
      parameter( mx = 2 )
_ENDIF
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
_IF(i8)
_IF(littleendian)
      data nshift/48,0,16,32/
_ELSE
      data nshift/0,48,32,16/
_ENDIF
_ELSE
_IF(littleendian)
      data nshift/16,0/
_ELSE
      data nshift/0,16/
_ENDIF
_ENDIF
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
_IF(i8)
      parameter( mx = 2 )
_ELSE
      parameter( mx = 1 )
_ENDIF
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

_IF(GFS)
CMR   the Fortran95/2003 standard for 'BOZ' constants does not allow the MSB to be set.
CMR   and gfortran does not use the postfix=typeless convention either
      integer ibits32          ! all 32 rightmost bits set
      integer*8 ibits64 
      integer*4 LR32(2)
      equivalence (ibits64,LR32(1))
_ELSE
      integer ibits32          ! all 32 rightmost bits set (we hope)
      parameter(ibits32 = z'ffffffff')
_ENDIF
_IF(i8)
_IF(littleendian)
      data nshift/32,0/
_ELSE
      data nshift/0,32/
_ENDIF
_ELSE
      data nshift/0/
_ENDIF
c
c     Code:
c
_IF(GFS)
CMR   It only hurts when I laugh
CMR   but we force the compiler to use I8 internally
      ibits64= z'00000000ffffffff'
      if (LR32(2).ne.ibits0) then
         ibits32=LR32(2)
      else if (LR32(1).ne.ibits0) then
         ibits32=LR32(1)
      else
         call caserr(' *** machscf: lost my bits from ibits64 ***')
      endif
_ENDIF
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
_IF(i8)
      parameter( mx = 2 )
_ELSE
      parameter( mx = 1 )
_ENDIF
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
_IF(GFS)
CMR   the Fortran95/2003 standard for 'BOZ' constants does not allow the MSB to be set.
CMR   and gfortran does not use the postfix=typeless convention either
      integer ibits32          ! all 32 rightmost bits set
      integer*8 ibits64 
      integer*4 LR32(2)
      equivalence (ibits64,LR32(1))
_ELSE
      integer ibits32          ! all 32 rightmost bits set (we hope)
      parameter(ibits32 = z'ffffffff')
_ENDIF
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
_IF(i8)
_IF(littleendian)
      data nshift/32,0/
_ELSE
      data nshift/0,32/
_ENDIF
_ELSE
      data nshift/0/
_ENDIF
cDEBUG
c     integer igotcha32
cDEBUG
c
c     Code:
c
_IF(GFS)
CMR   we force the compiler to use I8 internally
      ibits64= z'00000000ffffffff'
      if (LR32(2).ne.ibits0) then
         ibits32=LR32(2)
      else if (LR32(1).ne.ibits0) then
         ibits32=LR32(1)
      else
         call caserr(' *** machscf: lost my bits from ibits64 ***')
      endif
_ENDIF
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
_IF(i8)
      parameter( mx = 4 )
_ELSE
      parameter( mx = 2 )
_ENDIF
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
_IF(i8)
_IF(littleendian)
      data nshift/48,0,16,32/
_ELSE
      data nshift/0,48,32,16/
_ENDIF
_ELSE
_IF(littleendian)
      data nshift/16,0/
_ELSE
      data nshift/0,16/
_ENDIF
_ENDIF
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
_IF(i8)
      parameter( mx = 8 )
_ELSE
      parameter( mx = 4 )
_ENDIF
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
_IF(i8)
_IF(littleendian)
      data nshift/56,0,8,16,24,32,40,48/
_ELSE
      data nshift/0,56,48,40,32,24,16,8/
_ENDIF
_ELSE
_IF(littleendian)
      data nshift/24,0,8,16/
_ELSE
      data nshift/0,24,16,8/
_ENDIF
_ENDIF
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
_ENDIF(t3d)
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
_IF(cio,fcio,ipsc)
c
c** generic fget
c
      subroutine fget(gin,m,iunit)
      implicit REAL  (a-h,o-z)
      dimension gin(*)
      call find(iunit)
      call get(gin,m)
      return
      end
c
c** generic flags
c
      subroutine flags
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
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
_IF1(x)     *35x,'*               ===  CONVEX  version 8.0   ===         ',
_IF1(a)     *35x,'*               ===  ALLIANT version  8.0  ===         ',
_IF1(s)     *35x,'*               ===    SUN   version  8.0  ===         ',
_IF1(p)     *35x,'*               ===  APOLLO  version  8.0  ===         ',
_IF1(g)     *35x,'*               ===    SGI   version  8.0  ===         ',
_IF1(b)     *35x,'*               ===  HP-9000  version  8.0  ===        ',
_IF1(h)     *35x,'*              ===  IPSC/860  version 8.0   ===        ',
_IF(dec)
_IF(alpha)
_IF(osf)
     *35x,'*          ===   Compaq Tru64 version  8.0    ===      ',
_ELSE
     *35x,'*           ===    DEC AXP VMS version  8.0  ===       ',
_ENDIF
_ELSE
     *35x,'*               ===    DEC   version  8.0  ===         ',
_ENDIF
_ENDIF
_IF(macosx)
_IF1(r)     *35x,'*             ===  MACOSX-xlf version  8.0  ===        ',
_ELSE
_IF1(r)     *35x,'*               ===  RS-6000 version  8.0  ===         ',
_ENDIF
_IF1(t)     *35x,'*               ===   TITAN  version  8.0  ===         ',
_IF(t3e)
     *35x,'*                 === CRAY T3E version  8.0  ===           ',
_ELSEIF(t3d)
     *35x,'*                 === CRAY T3D version  8.0  ===           ',
_ELSE
_IF1(k)     *35x,'*               ===   KSR-2  version  8.0  ===         ',
_ENDIF
_IF(itanium)
     *35x,'*               ===   Itanium version  8.0  ===        ',
_ELSEIF(opteron)
     *35x,'*               ===   Opteron version  8.0  ===        ',
_ELSEIF(xeon)
     *35x,'*                ===   Xeon version  8.0  ===          ',
_ELSEIF(em64t)
     *35x,'*             ===   Intel EM64T version  8.0  ===      ',
_ELSE
_IF1(G)     *35x,'*               ===   Generic version  8.0  ===        ',
_ENDIF
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
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(convex,sgi,sun,dec,ksr)
      integer hostnm
_ELSEIF(vms)
      integer fgethostname
_ENDIF
_IF(convex,sgi,sun,dec,rs6000,hp700,hpux11,ksr,itanium,linux)
      character*44 host
      character*256 com
_ENDIF
      character*10 cdate,cname,cversion
      character*5 ctime
INCLUDE(common/sizes)
INCLUDE(common/timez)
INCLUDE(common/jinfo)
INCLUDE(common/iofile)
INCLUDE(common/timeperiods)
INCLUDE(common/segm)
INCLUDE(common/blksiz)
_IF(ga)
      MA_LOGICAL ga_uses_ma
_ENDIF
_IF(itanium)
       integer*4 null_4
       data null_4/0/
_ENDIF
c
c ----- initialise job
c
_IF(signals)
      call catch
_ENDIF
c
      call start_time_period(TP_ENTIRE)
      call inita(lword,mfree)
      call tidajt(zdate,ztime,zyear,zaccno,zanam,isecs)
c
      ti=0.0d0
      tx=0.0d0
      timlim=isecs
      call flags
_IFN(win95)
_IF(convex,sgi,sun,dec,hp700,hpux11,rs6000,ksr,itanium,linux)
      call gethost(host)
_IF(ksr,t3d)
      com=" "
_ELSEIF(itanium)
      call getarg(null_4,com)
c     call getarg(0,com)
_ELSE
      call getarg(0,com)
_ENDIF
      write(iwr,10) host, com
 10   format(
     *40x,'Hostname  : ',a44/
     *40x,'GAMESS-UK Executable: ',a70)
_ENDIF
_ENDIF(NOT win95)
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
     *40x,'main store requested',i10,' REAL  words',
     * ' =',i6,' mbytes'/)

_IF(ga)
      if (ga_uses_ma()) then
         write(iwr,40)
 40      format(
     &   40x,'main memory array used to hold global arrays instead of',
     &   ' shared memory'/)
      else
         write(iwr,50)
 50      format(
     &   40x,'shared memory used to hold global arrays instead of ',
     &   'main memory array'/)
      endif
_ENDIF
_IF(parallel)
c
c  print node / communications information
c
cccc
cccc  Moved down after  GA initialisation
cccc
cccc      call pg_pproc
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/iofile)
INCLUDE(common/discc)
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
      implicit REAL  (a-h,o-z)
      call clredx
      return
      end
c
c** generic gtnv
c
      subroutine gtnv(s,v)
      implicit REAL  (a-h,o-z)
      character s*(*), v*(*)
      logical ostart, ostop
      ostart = .false.
      ostop = .false.
      lens = len(s)
_IF(t3d,ipsc,unix)
      lenv = len(v)
_ENDIF
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
_IF(vms)
      call fgetenv(%ref(s(istart:istop)),%ref(v))
      v = ' '
_ELSEIF(t3d,ipsc,unix)
      call gtnvcc(s(istart:istop), istop-istart+1,v,lenv,ier)
_ELSE
      call getenv(s(istart:istop),v)
_ENDIF
      return
      end
c
c** generic quit
c
      subroutine quit(ia,n)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
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
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/iofile)
c
INCLUDE(common/blksiz)
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
_IF(ipsc)
2     write(iwr,6000)yed(i),ipos(i),iop(i)
_ELSE
2     write(iwr,6000)yed(i),ipos(i),ipos(i)
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *10 chpid, chnod
      character *132 rec132,ginput
INCLUDE(common/sizes)
      integer pid
_IF(unix)
_ELSEIF(sgi,hp700,hpux11,convex,rs6000,t3d)
      integer getpid
_ELSEIF(vms)
      integer decc$getpid
_ENDIF
INCLUDE(common/work)
INCLUDE(common/iofile)
_IF(parallel)
       dimension ibuff(132)
_ENDIF
_IF(cfs)
      character*132 zznam
INCLUDE(common/cfsinfo)
_ENDIF
_IF(taskfarm)
      integer global_nodeid
_ENDIF
_IF(mips4)
c
c PS - patch for failure to delete input file on power challenge
c
      common/inpfil/ginput
_ENDIF
      logical o_inclin,o_inc_err
c
c --- write a copy of the input 
c
c     
c   =====  define stream number for copy...
c
      o_inclin = .false.
      o_inc_err = .false.
_IF1(a)      irdcpy=59
_IF1(Grkh)      irdcpy=95
_IFN1(Garkh)      irdcpy=121

c
c  ====== define name for file
c
_IF(vms)
      pid=decc$getpid()
_ELSEIF(nec)
c
c experiment
c
      call pidcc(pid)
ccc      pid=pid+ipg_nodeid()
c
_ELSEIF(unix)
      call pidcc(pid)
_ELSEIF(sgi,hp700,hpux11,convex,rs6000,t3d)
      pid=getpid()
_ELSEIF(parallel)
_IF(cfs)
c
c all nodes will read from the same file
c
      pid=0
_ELSE
c
c fix to separate node files if on a shared file system
c
      pid=ipg_nodeid()
_ENDIF
_ELSE
      pid=0
_ENDIF
      write(chpid,'(i10)')pid
      do 1020 i=1,10
      if(chpid(i:i).ne.' ') then
      length=i
      go to 1031
      endif
1020  continue

_IF(taskfarm)
 1031 write(chnod,'(i10)')global_nodeid()
_ELSE
 1031 write(chnod,'(i10)')ipg_nodeid()
_ENDIF

      do i=1,10
         if(chnod(i:i).ne.' ') then
            length2=i
            go to 1030
         endif
      enddo

1030  ginput='gamess_input.'//chpid(length:10)//'.'//chnod(length2:10)

_IF(cfs)
c
c  This block is for parallel machines with a concurrent
c  file system (root node writes, and all nodes can read
c  so far, it has only been tested on the intel iPSC
c
c  - it relies on the pathname coming from an envir. var.
c    see routine inita
c
c ===== prepend /cfs directory, node 0 opens for writing
c
      if(lcfsdr .ne. 0)then
         zznam = zcfsdr(1:lcfsdr)//ginput
      else
         zznam = ginput
      endif
      call strtrm(zznam,length)

      if(opg_root()) then
c
c .. read from standard input
c      open(unit=5,iostat=ioerr,form='formatted',access='sequential',
c    * status='old',file='datain')
c      if(ioerr.ne.0) call caserr('error opening data input file')
c
         open(unit=irdcpy,iostat=ioerr,form='formatted',
     *        access='sequential',
     *        status='unknown',
     *        file=zznam(1:length) )
         print *,' just opened datain  for node ',minode

cjvl
 1060    if (.not.o_inclin) read(ird,1080,end=1090) rec132
c     include statement (allows including files in the input )
c
      call inclin(rec132,irdcpy-1,iwr,o_inclin,o_inc_err)
cjvl
         write(irdcpy,1080)rec132
         go to 1060
      endif
_ELSEIF(parallel)
c
c =====  broadcast input file, all write their own copy =====
c
c      write(6,*)'node',ipg_nodeid(),'opens',ginput
      open(unit=irdcpy,form='formatted',file=ginput,status='unknown')
_IF(datain)
c
c datain fix: this is needed when using some broken MPI implementations
c that do not read properly from stdin (sometimes they work for short
c files and hang for long ones)
c
       if (opg_root()) then
         open(unit=ird,iostat=ioerr,form='formatted',
     *        access='sequential',status='old',file='datain')
         if(ioerr.ne.0) call caserr('error opening data input file')
       endif
_ENDIF(datain)
 1060 continue
      if(opg_root())then
cjvl        ... include ...
         if (.not.o_inclin) read(ird,1080,iostat=nerr) rec132
         call inclin(rec132,irdcpy-1,iwr,o_inclin,o_inc_err)
cjvl
         if(nerr.eq.0)then
            call chtoi(ibuff(1),rec132)
         else
            ibuff(1) = 0
         endif
c        check if we found the include file
         if(o_inc_err) then
            ibuff(1) = -1
         endif
      endif
      liw=8/lenwrd()
      call pg_brdcst(100,ibuff,132*liw,0)

cjmht Additional synch need on hector (7/10/2010)
      call pg_synch(99)

      if(ibuff(1).eq.0)goto 1090
      if(ibuff(1).eq.-1)call caserr2('include file not found')
      if(.not.opg_root())call itoch(ibuff(1),rec132)
      write(irdcpy,1080)rec132
      go to 1060
_ELSE
c
c  spool data from ird to irdcpy
c
      open(unit=irdcpy,form='formatted',file=ginput,status='unknown')
_IF(datain)
c
c datain fix: this is needed when trying to debug gamess under
c Visual Studio .NET as there seems to be no option to redirect
c standard input.
c
       open(unit=ird,iostat=ioerr,form='formatted',access='sequential',
     *     status='old',file='datain')
       if(ioerr.ne.0) call caserr('error opening data input file')
_ENDIF(datain)
cjvl
1060  if (.not.o_inclin) read(ird,1080,end=1090) rec132
c     include statement (allows including files in the input )
c
      call inclin(rec132,irdcpy-1,iwr,o_inclin,o_inc_err)
cjvl
      write(irdcpy,1080)rec132
      go to 1060
_ENDIF
 1080 format(a132)
c
c ==== end of input file, rewind copy file 
c 
 1090 continue
      ird=irdcpy
_IF(cfs)
c
c ===== slave nodes open file written by node 0
c      
      call pg_synch(999)
      if(.not.opg_root()) then
         open(unit=95,iostat=ioerr,form='formatted',access='sequential',
     *        status='unknown',file=zznam(1:length) )
         if(ioerr.ne.0) call caserr(
     +        'error opening cfs data input file')
      endif
_ENDIF
      rewind ird
_IF(datain)
      close(unit=5)
_ENDIF
      return
      end
_IFN(parallel)
c
c** generic clredx (serial)
c
      subroutine clredx
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/disc)
      return
      end
c
c** generic shut (serial)
c
      subroutine shut
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,
     *            ipos(maxlfn),nam(maxlfn),
     *  keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/maxlen)
      call clredx
       irep = 0
      do 1 i=1,maxlfn
      if(ipos(i).gt.0 ) then
       if(.not.omem(i)) then
          if(keep(i).gt.0) then
_IF(fcio)
          close(unit=nam(i),iostat=irep,status='keep')
_ENDIF
_IF(cio)
          call closecc(nam(i))
_ENDIF
          istat(i)=2
          else
_IF(fcio)
          close(unit=nam(i),iostat=irep,status='delete')
_ENDIF
_IF(cio)
          call closecc(nam(i))
          call strtrm(zedfil(i),length)
_IFN(vms)
          call delcc(zedfil(i),length,ierrio)
_ELSE
          call delcc(%ref(zedfil(i)),length,ierrio)
_ENDIF
          if(ierrio.ne.0)call ioerr('delete',i,zedfil(i))
_ENDIF
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
_IFN1(Garkh)      close(unit=121,status='delete')
_IF1(Grkh)      close(unit=95,status='delete')
_IF1(a)      close(unit=59,status='delete')
      return
      end
c
c** generic find (serial)
c
      subroutine find(iunit)
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
INCLUDE(common/utilc)
INCLUDE(common/iofile)
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
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer nchk(2)
INCLUDE(common/sizes)
      dimension c(*)
c ... single-buffered  i/o
      common/bufa/abc(512)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 iword
_ENDIF
INCLUDE(common/iofile)
INCLUDE(common/utilc)
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
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn)
INCLUDE(common/maxlen)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      character *16 logfil
INCLUDE(common/sizes)
      dimension ytext(20),zinout(2)
      common/clogfl/logfil
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/work)
INCLUDE(common/timez)
INCLUDE(common/segm)
INCLUDE(common/gmempara)
INCLUDE(common/prints)
INCLUDE(common/disc)
INCLUDE(common/discc)
      common/bufa/aaaaa(8190)
INCLUDE(common/blksiz)
      common/blksiz2/ns2,ns2512(10),npstat,
     *  n2len,n2max ,nport
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     *  n3len,n3max ,ntort
      common/restri/nfile(63),lda(508),isect(508),ldx(508),
     *  iacsct(508)
INCLUDE(common/atmblk)
INCLUDE(common/maxlen)
      dimension yfname(maxfrt)
      dimension yynam(maxlfn),zdef(maxfrt)
c dump core on caserr
      common/crdmp/odumpc
c system derived error codes
INCLUDE(common/syscod)
_IF(chemshell)
      integer default_heap
_ENDIF
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
_IF(t3e)
      common/maskc/mask(64)
      data mask0,mask1/'0'x,'1'x /
_ENDIF
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
_IF(ksr)
      nsort=91
      nport=92
      ntort=93
_ELSE
      nsort=101
      nport=102
      ntort=103
_ENDIF
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
_IF(chemshell)
      lword = default_heap()
_ELSE
      lword=icoresize()
_ENDIF
      otime=.false.
      icode1 = 0
      odumpc = .false.
      oprintm = .false.
      omemin = .false.
_IF(debug)
      call gmem_set_nanify(.true.)
_ELSE
      call gmem_set_nanify(.false.)
_ENDIF
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
_IF(t3e)
c     set up mask (in setio routins on cray)
      mask(1)=ishft(mask1,63)
      do i=2,64
       mask(i)=ior(ishft(mask1,64-i),mask(i-1))
      enddo
_ENDIF
c
c  overwrite default 4 MW core for  semi-direct mrdci code
c  hopefully this is a short term requirement.
c  Its caused by the diag requirements with increased mxrootn
c
      if(opark.and. .not.omemin) then
       lword = 20000000
       opark = .false.
      endif
_IF(i8)
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
_ENDIF
c
      return
      end
c
c** generic forstr (serial)
c
      subroutine forstr
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character *132 filatt
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
      common/disc/isel(5+maxlfn*8),
     *ipos(maxfrt),keep(maxfrt),istatf(maxfrt),
     *lrecf(maxfrt),nrecf(maxfrt),oform(maxfrt)
c
INCLUDE(common/cfsinfo)
INCLUDE(common/ipsc)
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
_IFN(vms)
 20   if(oform(i)) then
_ELSE
 20   oform(6)=.true.
      if(oform(i)) then
_ENDIF
      open(unit=iunit,iostat=ioerr,file=filatt,status='unknown',
_IFN(vms)
     *form='formatted')
_ELSE
     *form='formatted',shared)
_ENDIF
      else
      open(unit=iunit,iostat=ioerr,file=filatt,status='unknown',
_IFN(vms)
     *form='unformatted')
_ELSE
     *form='unformatted',shared)
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/discc)
INCLUDE(common/cfsinfo)
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
INCLUDE(common/iofile)
INCLUDE(common/atmblk)
      common/vinteg/q(1)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
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
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
      common/vinteg/q(1)
INCLUDE(common/iofile)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
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
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension c(*)
c ... atmol  i/o
      common/bufa/abc(512)
_IF(i8drct)
      integer *8 ipack2, pack2
      equivalence (rpack2,ipack2)
_ENDIF
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
      common/vinteg/q(1)
INCLUDE(common/maxlen)
INCLUDE(common/utilc)
INCLUDE(common/iofile)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iword
_ENDIF
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
_IF(bits8)
        call pack2(irec,nword,q(iword+511))
_ELSEIF(i8drct)
        ipack2 = pack2(irec,nword)
        q(iword+511)= rpack2
_ELSE
        q(iword+511)=pack2(irec,nword)
_ENDIF
        call dcopy(nword,c,1,q(iword),1)
      else
        call dcopy(nword,c,1,abc,1)
_IF(bits8)
        call pack2(irec,nword,abc(512))
_ELSEIF(i8drct)
        ipack2 = pack2(irec,nword)
        abc(512)= rpack2
_ELSE
        abc(512)=pack2(irec,nword)
_ENDIF
        call putcc(nam(isel),abc,ierrio)
        if(ierrio.ne.0)call ioerr('write',isel,' ')
      endif
      return
      end
c
c** generic rdtpnm (serial)
c
      subroutine rdtpnm(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IFN(vms)
      character *132 filatt
_ELSE
      character *8 filatt
_ENDIF
      logical tell
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,
     * ipos(maxlfn),nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),
     * len(maxlfn*4)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/maxlen)
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
_IF(cio)
_IFN(vms)
      call opencc(filatt,lencha,nam(isel),ierrio)
_ELSE
      call opencc(%REF(filatt),lencha,nam(isel),ierrio)
_ENDIF
      if(ierrio.ne.0)call ioerr('open',isel,filatt)
_ENDIF
_IF(fcio)
      open(nam(isel),file=filatt,form='unformatted',
_IFN1(t)     & access='direct',status='unknown',recl=4096)
_IF1(t)     & access='direct',status='unknown',recl=1024)
_ENDIF
      if(ooprnt)write(iwr,1)yed(isel),length
  1    format(/' lfn ',a4,' allocated: current length =',
     * i5,' blocks')
      return
      end
c
c** generic delfil (serial)
c
      subroutine delfil(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/maxlen)

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
_IFN(vms)
             call delcc(zedfil(iunit),length,ierrio)
_ELSE
             call delcc(%ref(zedfil(iunit)),length,ierrio)
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *            nam(maxlfn),
     *  keep(maxlfn),istat(maxlfn),len(maxlfn),iop(maxlfn),
     *  iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/maxlen)
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
_ENDIF
_IF(parallel)
c
c  parallel I/O dependent on iwhere(iunit)
c
c      iwhere:  1 : virtual io to host / 
c               0 : concurrent filesystem (## files)
c      iwhere:  2 : concurrent filesystem  - node 0
c      iwhere:  3 : concurrent filesystem (single file)
c      iwhere:  4 : node memory (## file)
c      iwhere:  5 : node memory (node 0)
c      iwhere:  6 : GA (output to ## file)
c      iwhere:  7 : GA (output only on node 0)
c              -1 : io via host
c      isync :  1 : synchronous io  / 0 : asynchronous
c
c
c** generic clredx (parallel)
c
      subroutine clredx
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/nodinf)
INCLUDE(common/disc)
INCLUDE(common/nodeio)
INCLUDE(common/maxlen)
c
c...  clear i/o transfers on all files
c
      do 10 iunit=1,maxlfn

         if(ipos(iunit).eq.0) then
c - skip, file is not open
         else if(omem(iunit)) then
c - skip, file is in core
         else if (iwhere(iunit).eq.2 .and. .not. opg_root()) then
c - skip, file is not open on this node
         else if( iwhere(iunit).eq.2 .or.
     &           iwhere(iunit).eq.3 .or.
     &           iwhere(iunit).eq.0) then
c
c clear pending i/o
c
            call waitedn(iunit)
_IF(ga)
         Else If( iwhere( iunit ) .EQ. 6 .Or. 
     +            iwhere( iunit ) .EQ. 7 ) Then 
            Call clredx_ga
_ENDIF
         else
            call caserr('invalid iwhere in clredx')
         endif

10    continue
      return
      end
c
c** generic shut (parallel)
c
      subroutine shut
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
c
_IF(ipsc)
INCLUDE(common/ipsc)
_ENDIF
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/infoa)
c
       common /brdcst_stats/ ibrdcst(maxlfn), brdcst_l(maxlfn)
_IF(mips4)
c
c PS - patch for failure to delete input file on power challenge
c
       character*132 ginput
       common/inpfil/ginput
_ENDIF
c
c    print brdcst statistics (at present from direct-SCF)
c
      do loop = 1,maxlfn
       if(ibrdcst(loop).gt.0 ) then
        tmp = brdcst_l(loop) / 1024.0d0
        ltri = 4*(num+1)*num
        tri = ltri
        tri = brdcst_l(loop) / tri
        write(iwr,605) loop, ibrdcst(loop), tmp, tri
 605    format(' **** unit = ', i2, ' no. of 4K broadcasts = ', i8,
     +  ' amount transferred = ', f10.1, ' Kbytes (',
     +       f10.1,' triangles)')
        endif
      end do
c
c...  close all files
c
_IFN(ipsc)
      call clredx
_ENDIF
_IF(fcio,ipsc)
      irep = 0
_ENDIF
c
c close ed streams
c
      do 10 iunit=1,maxlfn

         if(ipos(iunit).ne.0 .and. isync(iunit).eq.0) then
            write(iwr,80010)iunit,maxio(iunit)
80010       format(1x,'Maximum wait count for unit',i2,' = ',i8)
         endif

         call delfil(iunit)
         ipos(iunit)=0
 10   continue

c
c ... close fortran files
c
      do 100 i=1,maxfrt
_IF(ipsc)
      if(.not.ohost(i).and.iposf(i).gt.0) then
_ELSE
      if(iposf(i).gt.0) then
_ENDIF
         if(keepf(i).gt.0) then
            close(unit=i,iostat=irep,status='keep')
         else
            close(unit=i,iostat=irep,status='delete')
         endif
         if(irep.ne.0)then
            write(iwr,5)zftfil(i)(1:80)
 5          format(1x,'***** error closing/deleting file',1x,a80)
         else
            if(ooprnt)write(iwr,110)zftn(i)(1:6)
 110        format(//' fortran file ',a6,' closed')
         endif
      endif
 100  continue
c
c delete spooled copy of input file
c
_IF(cfs)
      if(opg_root())close(unit=irdcpy,status='delete')
      close(unit=irdcpy,err=9998,status='delete')
      return
      end
_ELSEIF(mips4)
c
c PS - patch for failure to delete file
c
      close(unit=irdcpy,status='keep')
      call strtrm(ginput,length)
      call delcc(ginput,length,ierr)
      if(ierr.ne.0)call ioerr('delete',0,ginput)
      return
      end
_ELSE
      close(unit=irdcpy,err=9998,status='delete')
      return
9998  write(iwr,50)
50    format(1x,' **** problem closing copy of input file ***')
      return
      end
_ENDIF
c
c - interface as used in drv2e
c
      subroutine shut1(iunit)
      implicit none
      integer iunit
      call shutedn(iunit,0)
      end
c
c** shut shutedn 
c
c
c  - close stream on iunit
c
c  idel = 0 leave file intact
c       = 1 delete file if keep(iunit) is not set
c
      subroutine shutedn(iunit,idel)
c
      implicit none
      integer iunit, idel
c
INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
c
      integer iret
      logical odel
_IF(cio)
      integer length
_ENDIF
c
      odel = (idel .ne. 0) .and. (keep(iunit).eq.0)
_IF(fcio,ipsc)
      if(odel)then
         close(unit=nam(iunit),iostat=iret,status='delete')
         istat(iunit)=0
      else
         close(unit=nam(iunit),iostat=iret,status='keep')
         istat(iunit)=2
      endif
      if(iret.ne.0)then
         write(iwr,5)yed(iunit)
 5       format(1x,'***** error closing/deleting file',1x,a4)
      else
         if(ooprnt)write(iwr,3)yed(iunit)
 3       format(//' file ',a4,' closed')
      endif
_ELSE
      call closecc(nam(iunit))
      if(odel)then
         call strtrm(zedfil(iunit),length)
         call delcc(zedfil(iunit),length,iret)
         if(iret.ne.0)call ioerr('delete',iunit,zedfil(iunit))
      endif
_ENDIF
      if(odel)then
         istat(iunit)=0
      else
         istat(iunit)=2
      endif
      ipos(iunit)=0
      end
c
c** waitedn(iunit)
c
c  check for pending i/o on stream iunit
c
c  stream must be open
c
      subroutine waitedn(iunit)
      implicit none
INCLUDE(common/sizes)      
INCLUDE(common/nodeio)
INCLUDE(common/disc)
      integer iunit

      if (isync(iunit).eq.0) then
         if (ids48(iunit).ge.0) then
            call iowt(ids48(iunit))
            ids48(iunit)=-1
         end if
         if (idr48(iunit).ge.0) then
            call iowt(idr48(iunit))
            idr48(iunit)=-1
         end if

      end if
      return
      end
c
c** generic find (parallel)
c
      subroutine find(iunit)
_IF(ipsc)
c...
c...   asynchronous (virtual) host io
c...   asynchronous  + cfs  io
c...
_ENDIF
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/nodinf)
INCLUDE(common/timeperiods)
INCLUDE(common/nodeio)
INCLUDE(common/disc)
c
      common/bufa/abc(512*maxlfn)
c ====
      integer timflag
      dimension timflag(-1:7)
      data timflag/TP_IOFM1,TP_IOF0,TP_IOF1,TP_IOF2,
     &     TP_IOF3,TP_IOF4,TP_IOF5,TP_IOF6, TP_IOF7/
c
      isel = iunit
      irec=ipos(isel)
      nstart=(iunit-1)*512+1
cPS ipsc has fail code 3 here?
      if (irec.lt.1) call iofail(iunit,1)
c
      call start_time_period(timflag(iwhere(iunit)))

      if (opg_root().and.iwhere(iunit).eq.1) then
cjvl     handshake
         if (oputpp(iunit)) then
            itype = 999990+iunit
             call cprob(itype)
            if(infocount().le.0) then
_IF(ipsc)
            call crecv(itype,dum,0)
_ELSE
            call caserr('crecv only on ipsc')
_ENDIF
            else
               print *,' node ',ipg_nodeid(),
     +                 ' **** find error on unit ',isel
            endif
cjvl        handshake
         end if
         oputpp(iunit) = .false.
         if (ids48(iunit).ge.0) call msgw(ids48(iunit))
c
c ---   ask it from host memory (node-0)
c
         ityp48 = (maxlfn+iunit)*10000+irec
         ids48(iunit) = isend(ityp48,buf,0,mihost,mpid)
         idr48(iunit) = irecv(ityp48,abc(nstart),4096)
c
      else  if (iwhere(iunit).eq.0 .or.
     +          iwhere(iunit).eq.3 ) then
c
c...   concurrent io on each node
c
         if (isync(iunit).eq.0) then
_IFN(charmm)
c...   asynchronous
            if (ids48(iunit).ge.0) then
              call iowt(ids48(iunit))
              ids48(iunit)=-1
            end if
            idr48(iunit) = iread(nam(iunit),abc(nstart),4096)
_ELSE
            call caserr('async i/o clashes with charmm')
_ENDIF
         end if
c
      else  if (opg_root().and.iwhere(iunit).eq.2) then
c
c...   concurrent io on node 0
c
c...   asynchronous
         if (isync(iunit).eq.0) then
_IFN(charmm)
            if (ids48(iunit).ge.0) then
              call iowt(ids48(iunit))
              ids48(iunit)=-1
            end if
            idr48(iunit) = iread(nam(iunit),abc(nstart),4096)
_ELSE
            call caserr('async i/o clashes with charmm')
_ENDIF
         end if
c
      end if
c
      call end_time_period(timflag(iwhere(iunit)))

      return
      end
c
c** generic get (parallel)
c
      subroutine get(c,nword)
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension c(*)
c
INCLUDE(common/disc)
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
c
      common/vinteg/q(1)
INCLUDE(common/maxlen)
c
      common/bufa/abc(512*maxlfn)
      common /brdcst_stats/ ibrdcst(maxlfn), brdcst_l(maxlfn)
c
INCLUDE(common/timeperiods)
c
_IF(64bitpointers)
      integer*8 iword
_ENDIF
      integer timflag
      dimension timflag(-1:5)
      data timflag/TP_IOGM1,TP_IOG0,TP_IOG1,TP_IOG2,
     &     TP_IOG3,TP_IOG4,TP_IOG5/

_IF(ga)
      If( iwhere( isel ) .EQ. 6 .Or. iwhere( isel ) .EQ. 7 ) Then 
         Call get_ga( c, nword )
      Else
_ENDIF

      call start_time_period(timflag(iwhere(isel)))

      iunit = isel
      irec = ipos(iunit)
      nstart=(iunit-1)*512+1
c
      if (iwhere(iunit).eq.1) then
c...     virtual io
         if (opg_root()) then
c        root has to get from host and send
            if(ids48(iunit).ge.0)then
             call msgw(ids48(iunit))
             ids48(iunit) = -1
            end if
            call msgw(idr48(iunit))
            idr48(iunit) = -1
c...        record received  send
            itypn = (maxlfn+iunit)*10000+irec
_IF(ipsc)
            call csend(itypn,abc(nstart),4096,-1,mpid)
_ELSE
            call caserr('csend only on ipsc')
_ENDIF
         else
c...    not a node root so receive
            itypn = (maxlfn+iunit)*10000+irec
            call cprob(itypn)
            if(infocount().le.4096) then
_IF(ipsc)
               call crecv(itypn,abc(nstart),4096)
_ELSE
               call caserr('crecv only on ipsc')
_ENDIF
            else
               print *,' node ',minode,' **** get error on unit ',
     =                   isel
            endif
         end if
      else if (iwhere(iunit).eq.0 .or.
     +         iwhere(iunit).eq.3   ) then
c...   concurrent io -- same for all
         if (isync(iunit).eq.0) then
c...   asynchronous
            call iowt(idr48(iunit))
            idr48(iunit)=-1
         else
_IF1(h)            call cread(nam(iunit),abc(nstart),4096)
_IFN1(h)            call getcc(nam(isel),abc(nstart),ierrio)
_IFN1(h)            if(ierrio.ne.0)call ioerr('read',isel,' ')
         end if
      else if (iwhere(iunit).eq.2) then
         if (opg_root()) then
c        root has to get from cfs and send
          if (isync(iunit).eq.0) then
            call iowt(idr48(iunit))
            idr48(iunit)=-1
         else
_IF1(h)            call cread(nam(iunit),abc(nstart),4096)
_IFN1(h)            call getcc(nam(isel),abc(nstart),ierrio)
_IFN1(h)            if(ierrio.ne.0)call ioerr('read',isel,' ')
         end if
            if (oswed3(iunit)) then
c...        record received  send
_IF(ipsc)
              itypn = (maxlfn+iunit)*10000+irec
              call csend(itypn,abc(nstart),4096,-1,mpid)
_ELSE
              itypn = iunit*1000+irec
              call pg_brdcst(itypn,abc(nstart),4096,0)
_ENDIF
             ibrdcst(iunit) = ibrdcst(iunit) + 1
             brdcst_l(iunit) = brdcst_l(iunit) + 4096.0d0
            endif
         else
c...    not a node root so receive
            if (oswed3(iunit)) then
_IF(ipsc)
            itypn = (maxlfn+iunit)*10000+irec
            call cprob(itypn)
            if(infocount().le.4096) then
               call crecv(itypn,abc(nstart),4096)
            else
               print *,' node ',minode,' **** get error on unit ',
     +                   isel
            endif
_ELSE
            itypn = iunit*1000+irec
            call pg_brdcst(itypn,abc(nstart),4096,0)
_ENDIF
           endif
         endif
      else if (iwhere(iunit).eq.4) then
c
c     read from node memory
c
         iword = (irec-1)*512+iqqoff(iunit)+1
         call upack2(q(iword+511),iibl,nword)
         if (nword.gt.512.or.nword.lt.0.
     *        or.iibl.ne.ipos(iunit)) then
            print *,' node ',minode,' **** memory error on unit ',iunit
            print *,' node ',minode,' **** block read was',iibl
            print *,' node ',minode,' **** position was  ',ipos(iunit)
            call caserr('error in node memory access')
         end if
         call dcopy(nword,q(iword),1,c,1)
         go to 200
c
      else if (iwhere(iunit).eq.5) then
       if (opg_root()) then
c     read from node memory (node 0)
c
         iword = (irec-1)*512+iqqoff(iunit)+1
         call dcopy(512,q(iword),1,abc(nstart),1)
         call upack2(q(iword+511),iibl,nword)
         if (nword.gt.512.or.nword.lt.0.
     *   or.iibl.ne.ipos(iunit)) then
        print *,' node ',minode,' **** memory error on unit ',iunit
        print *,' node ',minode,' **** block read was',iibl
        print *,' node ',minode,' **** position was  ',ipos(iunit)
          call caserr('error in node memory access')
          end if
          if(oswed3(iunit)) then
c...         record received  send
_IF(ipsc)
           itypn = (maxlfn+iunit)*10000+irec
           call csend(itypn,abc(nstart),4096,-1,mpid)
_ELSE
           itypn = iunit*1000+irec
           call pg_brdcst(itypn,abc(nstart),4096,0)
_ENDIF
          endif
       else
c...    not a node root so receive
        if(oswed3(iunit)) then
_IF(ipsc)
          itypn = (maxlfn+iunit)*10000+irec
          call cprob(itypn)
          if(infocount().le.4096) then
             call crecv(itypn,abc(nstart),4096)
          else
          print *,' node ',minode,' **** get error on unit ',isel
            endif
_ELSE
          itypn = iunit*1000+irec
          call pg_brdcst(itypn,abc(nstart),4096,0)
_ENDIF
        endif
       endif
      else
c...   host-io synchronous only
         call synget(abc(nstart),512,nam(iunit),irec)
      end if
c...
c...    received => check
c...
      call upack2(abc(iunit*512),iibl,nword)
cjvl       call checkipos(iunit,ipos(iunit),iibl,nword)
      if (nword.gt.512.or.nword.lt.0.
     * or.iibl.ne.ipos(iunit)) then
         print *,'nword check',nword
         print *,' node ',minode,' **** error on unit ',iunit
         print *,' node ',minode,' **** block read was',iibl
         print *,' node ',minode,' **** position was  ',ipos(iunit)
         call caserr('i/o error detected')
      end if
      call dcopy(nword,abc(nstart),1,c,1)
c
 200  ipos(isel)=irec+1
c
c      write(6,*)'get done',ipg_nodeid()
c
      call end_time_period(timflag(iwhere(iunit)))

_IF(ga)
      End If
_ENDIF
c
      return
      end

_IF()
cjvl       kept for times when a parallel io-error occurs
      subroutine checkipos(iunit,ipos,iibl,nword)
c
c...   remco all reduce tells if a problem occured
c
      implicit integer (a-z)
INCLUDE(common/iofile)
c
      ii = 0
      if (ipos.ne.iibl.or.nword.gt.512.or.nword.lt.0) ii = 1
c
c...   gather individual alarm cry's
c
      call pg_igop(123456,ii,1,'+')
      call pg_synch(1234)
c
      if (ii.gt.0) then 
         write(iwr,600) ii
600      format(1X,/,' ***** IPOS errpr on ',i6,' nodes')
         if (ipos.ne.iibl.or.nword.gt.512.or.nword.lt.0) then
            print *,' error remco ',ipg_nodeid(),' : ',' ipos ',ipos,
     1              ' iibl ',iibl,' nword ',nword,' iunit ',iunit
         else
            print *,' remco ',ipg_nodeid(),' : ',' ipos ',ipos,
     1              ' iibl ',iibl,' nword ',nword,' iunit ',iunit
         end if
         call pg_synch(1234)
         if (iibl.ne.ipos) call caserr('i/o error detected')
      end if
c
      return
      end
_ENDIF
      
c
c** generic search (parallel)
c
      subroutine search(iblk,iunit)
c
c --- position a stream at the necessary point
c --- open it if necessary
c --- files start at 0 / atmol starts at block 1
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

INCLUDE(common/sizes)
INCLUDE(common/disc)
_IF(ga)
      If( iwhere( iunit ) .EQ. 6 .Or. iwhere( iunit ) .EQ. 7 ) Then 
         Call search_ga( iblk, iunit )
      Else
_ENDIF
c
      if(ipos(iunit).le.0)call rdtpnm(iunit)
c...
c...  now real search
c...
      if (iwhere(iunit).eq.0 .or.
     &    (iwhere(iunit).eq.2 .and.  opg_root()) .or. 
     &    iwhere(iunit).eq.3)then
_IF(ipsc)
        irep = lseek(nam(iunit),(iblk-1)*4096,0)
_ELSE
        call srchcc(nam(iunit),iblk,ierrio)
        if(ierrio.ne.0) call ioerr('search',iunit,' ')
_ENDIF
      endif
      ipos(iunit)=iblk
      return
_IF(ga)
      End If
_ENDIF
      end
c
c** generic inita (parallel)
c
      subroutine inita(lword,mfree)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*16 logfil
INCLUDE(common/sizes)
      dimension ytext(21),zinout(2)
      common/clogfl/logfil
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/prints)
INCLUDE(common/work)
INCLUDE(common/timez)
INCLUDE(common/segm)
INCLUDE(common/ga_file)
      common/scftim/tdiag(10)
c
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/cfsinfo)
INCLUDE(common/ipsc)
      parameter(nbuf=40)
      common/bufa/aaaaa(512*nbuf)
      common/blksiz/nsz,nsz512(10),nsstat,
     *  nslen,nsmax ,nsort,nstid(2),oasync,lensrt
      common/blksiz2/ns2,ns2512(10),npstat,
     *  n2len,n2max ,nport,nptid(2),opsync,lenprt
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     *  n3len,n3max ,ntort,nttid(2),otsync,lentrt
c
c***   ***NODE-IPSC***
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
      common /brdcst_stats/ ibrdcst(maxlfn), brdcst_l(maxlfn)
      character*11 fname
c***   ***NODE-IPSC***
INCLUDE(common/parcntl)
INCLUDE(common/fock_critic)
      common/restri/nfile(63),lda(508),isect(508),ldx(508),
     *  iacsct(508)
INCLUDE(common/atmblk)
INCLUDE(common/maxlen)
c dump core on caserr
      common/crdmp/odumpc
c system derived error codes
INCLUDE(common/syscod)
c
      logical onodid(maxlfn)
_IF(chemshell)
      integer default_heap
_ENDIF

      character*132 zztemp
      data ytext/'chan','swit','file','seti','memo','stor','core',
     * 'lpag','time','sort','psor','driv','dcor','max','logf',
     * 'tabl','node','dire','para','noec','fipc' /

      data zblank,zkeep,zold,zlength/' ','keep','old','length'/
      data zdel/'delete'/
      data yout,yprin/'outp','prin'   /
      data  zinout/  'dataout','punch'/

c
      ioupd = 0
      nslen = 0
      n2len = 0
      n3len = 0
      nsmax = -1
      n2max = -1
      n3max = 0
      lensrt = 0
      lenprt = 0
      lentrt = 0
      nsz=10
      ns2=1
      ns3=1
_IF(ksr,ipsc)
      nsort=91
      nport=92
      ntort=93
_ELSE
      nsort=101
      nport=102
      ntort=103
_ENDIF
      oasync=.true.
      opsync=.true.
      otsync=.true.
      ichek=0
      isel= 0
      iselr= 0
      iselw= 0
      irep= 0
      oprintm = .false.
_IF(debug)
      call gmem_set_nanify(.true.)
_ELSE
      call gmem_set_nanify(.false.)
_ENDIF
      call gmem_set_debug(.false.)
      call gmem_set_quiet(.true.)
c  .. clear timing accumulators
c
      do 3010 loop=1,10
3010  tdiag(loop)=0.0d0

c***   ***NODE-IPSC***
      oprint(60)=.true.
      opnode = .false.
      opkill=.true.
      do 4000 loop=1,maxlfn
      ids48(loop) = -1
      idr48(loop) = -1
_IF(parallel)
      oswed3(loop) = .true.
      ibrdcst(loop) = 0
      brdcst_l(loop) = 0.0d0
_ENDIF
      maxio(loop) = 0
      oputpp(loop) = .false.
4000  maxnod(loop)=0

c***   ***NODE-IPSC***
      do 3100 loop = 1,maxfrt
      zftstat(loop)=   'unknown'
      oform(loop) = .false.
_IF(ipsc)
      ohost(loop) = .false.
_ENDIF
      iposf(loop)=0
      keepf(loop)=0
3100  istatf(loop)=0
c
      zcfsdr = zblank
      lcfsdr = 0
_IF(cfs)
_IF(ipsc)
c
c ... determine CFS directory from nqs environment
c
c  assume /cfs/username/- very intel specific
c  this will need to be modified for the paragon
c
      call gtnv('QSUB_USER',zztemp)
      if (zztemp.ne.' ') then
         zcfsdr = '/cfs/'//zztemp
         zztemp = zcfsdr
         call strtrm(zztemp,lcfsdr)
         zcfsdr = zztemp(1:lcfsdr)//'/'
         lcfsdr = lcfsdr + 1
         if(opg_root())then
            write(iwr,*)
     &           'cfs directory set from environment var QSUB_USER - ',
     &           zcfsdr(1:lcfsdr)
         endif
      else
         write(iwr,*)
     &    '** set environment var QSUB_USER to /cfs  subdirectory **'
         call caserr('cfs directory could not be determined')
      endif
_ENDIF
_ENDIF
      numf = 60
      do 3500 loop=1,maxlfn
      iop(loop)  = 0
      keep(loop) = 0
      oform(loop) = .false.
      istat(loop) = 1
      ipos(loop) = 0
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
1113  mach(4)= -1
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
_IF(chemshell)
      lword = default_heap()
_ELSE
      lword=icoresize()
_ENDIF
      otime=.false.
      isecs=0
      mfree=0
      odumpc = .false.
      icode1 = 0
_IF(ipsc)
c
c --- open explicit input file
c
_ENDIF
      call setsto(maxlfn,0,len)
_IF(ipsc)
      call setsto(maxlfn,0,isync)
c
c ==== tpf for synchronous processing of mainfile
c      problem with PGroup = now reset
c
      isync(3) =0
c
_ELSE
c     default synchronous I/O
      call setsto(maxlfn,1,isync)
_ENDIF
c
      call setstl(maxlfn,.false.,omem)  
c
c ====
      call setsto(maxfrt,0,nrecf)
      call setsto(maxfrt,0,lrecf)
c
c ... default ed3 and ed7, ed14 to cfs disk
c ... other ed streams to concurrent io,
c     with unique file names
c
      call setsto(maxlfn,0,iwhere)
      call setstl(maxlfn,.true.,onodid)
c
c      ed0
      iwhere(1) = 2
      onodid(1) = .not. opg_root()
c     ed3
      iwhere(4) = 2
      onodid(4) = .not. opg_root()
c     ed7
      iwhere(8) = 2
      onodid(8) = .not. opg_root()
c
c     where are ed14 and ed15 used ??
      iwhere(15) = 2
      onodid(15) = .not. opg_root()
c ed15 shared file
      iwhere(16)= 3
      onodid(16) = .false.
      ogafsmode = .false.
_IF(ga)
c periodic update of GA files
c extend to in core later
c default at end of scf; more frequent dumping driven
c through iofile 
      ioupd = 1
_ENDIF


_IF(ga)
c
c ogafsmode controls the automatic routing of ed3, ed7 to GAs
c disable by default 

c
c NB previous behaviour was to enable by default for the t3d,
c t3e, sp2, and beowulf clusters
c
_IF(edgafs)
c     Put the ed files in the GAs by default when edgafs m4 keyword present
      ogafsmode = .true.
_ELSE
      ogafsmode = .false.
_ENDIF(edgafs)

_ENDIF(ga)

      ovirt = .false.
c
c --- save filenames and change the default to 'keep' if
c     environment variable set
c
      do 115 i=1,maxlfn
         call gtnv(yed(i),zztemp)
         if(zztemp.ne.' ')then
            keep(i)=1
            zedfil(i)=zztemp
            oedpr(i)=.true.
         else
            zedfil(i)=yed(i)
         endif
c...   change properties based on name (nzo for now)
       if (index(zztemp,'nzo').ne.0.or.index(zztemp,'cfs-0').ne.0.or.
     1     index(zztemp,'NZO').ne.0.or.index(zztemp,'FFS-0').ne.0) then
          iwhere(i)=2
          onodid(i) = .not. opg_root()
       end if
 115  continue
c
c ---- and repeat for the fortran files
c
      do 116 i=1,maxfrt
         call gtnv(zftn(i),zztemp)
         if(zztemp.ne.' ')then
            keepf(i)=1
            zftfil(i)=zztemp
            oftpr(i)=.true.
         else
            zftfil(i)=zftn(i)
         endif
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
300   jrec=jrec-1
 310  call prinit(oreact,opark,odci,odlfi)
      if (odci) then
c      write(6,*)' allocating GAs for ed19'
c     only set direct-ci GAs if direct-CI is invoked
c rwah ed19 common c/z file (global array) for parallel direct-ci
c
       iwhere(20) = 6
       onodid(20)= .not. opg_root()
       oadd(20) = .true.
      endif
      go to 320
 10   call input
 320  call inpa4(yyy)

      j=locatc(ytext,21,yyy)
      if(j)20,30,20
 20   go to (40,40,40,120,50,50,50,60,70,200,240,210,245,250,900,
     * 285,260,275,280,290,295),j
c
c  == fipc ==   shared memory fock building ;  for crossversion compatibility
c
295   ofipc = .false.
298   call inpa4(yyy)
      if (yyy.eq.'debu') ofipc_debug = .true.
      if (yyy(1:3).eq.'shm') ofipc_shmp  = .true.
      if (yyy.eq.'nofo') ofipc_nofock = .true.
      if (yyy.eq.'root') ofipc_root = .true.
      if (yyy.ne.' ') go to 298
      if (oroot()) then
         write(iwr,296)
         if (ofipc_debug) write(iwr,297) 'DEBUG'
         if (ofipc_shmp) write(iwr,297) 'SHMPR'
         if (ofipc_root) write(iwr,297) 'ROOT '
         if (ofipc_nofock) write(iwr,*)
     1 ' *&*& no FOCK build, energies STUPID *&*&'
296      Format(1x,/,60x,'*****************',
     1             /,60x,'*   NO FIPC     *',
     1             /,60x,'*****************')
297      Format(     60x,'*     ',a5,'     *',
     1             /,60x,'*****************')
      end if
      go to 10

c
c  == dire ==      specify default directory for all files
c
275   call inpanpc(zcfsdr)
      call strtrm(zcfsdr,lcfsdr)
      if(zcfsdr(lcfsdr:lcfsdr) .ne. '/')then
         lcfsdr = lcfsdr + 1
         zcfsdr(lcfsdr:lcfsdr) = '/'
      endif
      go to 10
c
c  == para ==      parallel control pre-directive
c
280   continue
_IF(parallel)
      call inpa(ztext)
      yyyy=ytrunc(ztext)
      if(yyyy.eq.'    ') then
         goto  10
      else if(yyyy.eq.'chun') then
c
c  set chunk setting
c
         call inpa(ztext)
         yyyy=ytrunc(ztext)
         if(yyyy.eq.'limi')then
c
c  set maximum chunk size
c
            call inpi(limchnk)
            if(limchnk.le.0)call caserr('bad chunk limit')
            write(iwr,*)'chunk size limit is ',limchnk
         else if(yyyy.eq.'task')then
c
c  set optimum tasks per node
c
            call inpi(ntchnk)
            write(iwr,*)'chunk size based on ',ntchnk,' tasks/SCF cycle'
         else
            call caserr('bad parallel chunk control directive')
         endif
      else if(yyyy.eq.'diag') then
c
c diagonalisation method
c
 100     call inpa(ztext)
         call intval(ztext,itmp,onumb)
         yyyy=ytrunc(ztext)
         if(onumb)then
            idpdiag = itmp
         else if(yyyy.eq.'alwa')then
            idpdiag = 0
         else if(yyyy.eq.'neve')then
            idpdiag = 99999999
         else if(yyyy.eq.'peig')then
            ipdiagmode = IDIAG_PEIGS
         else if(yyyy.eq.'gawr')then
            ipdiagif = IDIAG_GAWRAP
         else if(yyyy.eq.'pdsy')then
            if     (ztext.eq.'pdsyev' )then
               ipdiagmode = IDIAG_PDSYEV
            else if(ztext.eq.'pdsyevx')then
               ipdiagmode = IDIAG_PDSYEVX
            else if(ztext.eq.'pdsyevd')then
               ipdiagmode = IDIAG_PDSYEVD
            else if(ztext.eq.'pdsyevr')then
               ipdiagmode = IDIAG_PDSYEVR
            else
               call caserr('invalid ScaLAPACK diagonalizer specified')
            endif
         else if(yyyy.eq.'    ')then
            goto 110
         else 
            jrec = jrec -1
            goto 110
         endif
         goto 100
 110     continue
         if(idpdiag .lt. ipg_nnodes())then
            if(opg_root())write(iwr,*)
     &  'WARNING: Parallel diag unsafe for dimension '//
     &   '< nnodes, set to nnodes'
            idpdiag = ipg_nnodes()
         endif
      else if(yyyy.eq.'inve') then
c
c       parallel matrix inversion
c
        call inpa(ztext)
        yyyy=ytrunc(ztext)
        if (yyyy.eq.'chol') then
          ipinvmode = INV_CHOLESKY
        else if (yyyy.eq.'diag') then
          ipinvmode = INV_DIAG
        else
          call caserr('invalid inversion algorithm specified')
        endif
      else if(yyyy.eq.'diis') then
c
c parallel diis option
c
         call inpa(ztext)
         call intval(ztext,itmp,onumb)
         yyyy=ytrunc(ztext)

         if(onumb)then
            idpdiis = itmp
         else if(yyyy.eq.'alwa')then
            idpdiis = 0
         else if(yyyy.eq.'neve')then
            idpdiis = 99999999
         else if(yyyy.eq.'    ')then
            idpdiis = 0
         else 
            idpdiis = 0
            jrec = jrec -1
         endif
      else if(yyyy.eq.'mult') then
c
c parallel mult2 option
c
         call inpa(ztext)
         call intval(ztext,itmp,onumb)
         yyyy=ytrunc(ztext)

         if(onumb)then
            idpmult2 = itmp
         else if(yyyy.eq.'alwa')then
            idpmult2 = 0
         else if(yyyy.eq.'neve')then
            idpmult2 = 99999999
         else if(yyyy.eq.'    ')then
            idpmult2 = 0
         else 
            idpmult2 = 0
            jrec = jrec -1
         endif
      else if(yyyy.eq.'orth') then
c
c parallel orthog option
c
         call inpa(ztext)
         call intval(ztext,itmp,onumb)
         yyyy=ytrunc(ztext)

         if(onumb)then
            idporth = itmp
         else if(yyyy.eq.'alwa')then
            idporth = 0
         else if(yyyy.eq.'neve')then
            idporth = 99999999
         else if(yyyy.eq.'    ')then
            idporth = 0
         else 
            idporth = 0
            jrec = jrec -1
         endif
c 
      else if(yyyy.eq.'prin') then
         call inpa(ztext)
         yyyy=ytrunc(ztext)
         if(yyyy.eq.'all ')then
            iparapr = 2
         else if(yyyy.eq.'none')then
            iparapr = 0
         else if(yyyy.eq.'ga')then
            iparapr = 3
         else
            call caserr('bad parallel print keyword')
         endif
      else if(yyyy.eq.'mame')then
c
c MA memory allocation
c
         call caserr('para mememory directive obsolete - use core')
        
      else if(yyyy.eq.'iomo') then
c
c iomode
c
         call inpa(ztext)
         yyyy=ytrunc(ztext)
         if(yyyy.eq.'indi')then
c
c individual files (see code at end of
c this file to set iwhere on ed3/7
c
            ipiomode = IO_A
            if(ogafsmode) ogafsmode = .false.
         else if(yyyy.eq.'scre')then
            ipiomode = IO_NZ_S
            if(ogafsmode) ogafsmode = .false.
         else if(yyyy.eq.'nzo ')then
            ipiomode = IO_NZ
            if(ogafsmode) ogafsmode = .false.
         else if(yyyy.eq.'gafs')then
_IFN(ga)
            call caserr(
     &   'IOMODE GAFS not available in this (non-GA) build')
_ENDIF
            ipiomode = IO_NZ
            ogafsmode = .true.
            if(jrec .lt. jump)then
               call inpa4(ytemp)
               if(ytemp .eq. 'noup')then
                   ioupd = 0
               else if (ytemp.eq.'upda') then
                   ioupd = 2
               else
                   ioupd = 1
               endif
            endif
         else if(yyyy.eq.'noga')then
            ogafsmode = .false.
         else
            call caserr('bad parallel iomode keyword')
         endif
      else if(yyyy.eq.'test') then
c
c request diagnostic code for testing comms
c
         iptest = 1
      else 
         call caserr('bad parallel control parameter')
      endif
_ELSE
      write(iwr,*)"***** parallel directive ignored *****"
_ENDIF
      goto 280
c
c noecho - suppress output of input file (processed in
c          prinit)
c
 290  continue
      goto 10

c...  node
c...  special info only for other nodes (root excluded)
c...  also allows directives for possible host
c...  ****primitive***
c
260   continue
      call inpa4(yyy)
      if (yyy.eq.'    ') go to  10
      if (yyy.eq.'prin') then
         opnode = .true.
         go to 260
      end if
      if(yyy.eq.'outp') then
      opkill = .false.
      go to 260
      end if
      j=locatc(yed,maxlfn,yyy)
      if (j.le.0) go to 260
      ispec1 = j
      call inpa4(yyy)
      if (yyy.eq.'virt') ivirt1 = 1
      if (yyy.eq.'real') ivirt1 = -1
      call inpa4(yyy)
      j=locatc(yed,maxlfn,yyy)
      if (j.le.0) go to 260
      ispec2 = j
      if (yyy.eq.'virt') ivirt2 = 1
      if (yyy.eq.'real') ivirt2 = -1
      go to 260
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
      call inpan(logfil)
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
 410   call inpa4(yyy)
       if(yyy.eq.' ') go to 10
       if(yyy.eq.'sync') oasync = .false.
       if(yyy.eq.'leng') then
        call inpi(lensrt)
       endif
       go to 410
c
c  psort
c
 240   call inpi(ns2)
 510   call inpa4(yyy)
       if(yyy.eq.' ') go to 10
       if(yyy.eq.'sync') opsync = .false.
       if(yyy.eq.'leng') then
        call inpi(lenprt)
       endif
       go to 510
c
c  table
c
 285   call inpi(ns3)
 610   call inpa4(yyy)
       if(yyy.eq.' ') go to 10
       if(yyy.eq.'sync') otsync = .false.
       if(yyy.eq.'leng') then
        call inpi(lentrt)
       endif
       go to 610
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
      j=locatc(zftn,maxfrt,ztext)
      if(j.eq.0) then
         j=locatc(yft,maxfrt,ztext(1:4))
         if (j.eq.0) then
            if (ztext.eq.'punch') j = 58
         end if
         if(j.eq.0) go to 150
      endif
c
c      fortran sequential/vbs
c
      oftpr(j) = .true.
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
      if(ztext(1:4).eq.'host') then
      ohost(j) =.true.
      go to 800
      endif
      if(ztext.ne.zlength)go to 800
      call inpi(nrecf(j))
       call inpi(lrecf(j))
      if(nrecf(j).le.0.or.nrecf(j).gt.999999)call caserr
     *('invalid length for fortran file')
      go to 800
c
c      512 word files
c
 150  j=locatc(yed,maxlfn,ztext(1:4))
      if(j.eq.0) call caserr(
     *     'invalid lfn specified in file directive')
      oedpr(j) = .true.

      call inpa(ztext)
      if(ztext(1:4).eq.'memo') then
         omem(j) = .true.
c
c by default ed3,7 routed through node0, rest have
c independent memory files
c NB may reset 5 -> 4 at end of inita based on ipiomode=IO_A
c
         if(j .eq. 4 .or. j .eq. 8)then
           iwhere(j) = 5
         else
           iwhere(j) = 4
         endif
      else
         jrec=jrec-1 
         call inpanpc(zedfil(j))
      endif
      keep(j)=1
 80   call inpa(ztext)
      if(ztext.eq.zblank)then
         go to 10
      else if(ztext.eq.zkeep)then
         keep(j)=1
      else if(ztext.eq.zdel)then
         keep(j)=0
      else if(ztext.eq.'cfs') then
         iwhere(j)=0
      else if(ztext.eq.'cfs-0'.or.ztext.eq.'nzo') then
         iwhere(j)=2
         onodid(j) = .not. opg_root()
      else if(ztext.eq.'cfs-3') then
         iwhere(j)=3
      else if(ztext.eq.'virtual') then
         iwhere(j)=1
      else if(ztext(1:4).eq.'sync') then
         isync(j)=1
      else if(ztext(1:4).eq.'asyn') then
         isync(j)=0
      else if(ztext(1:4).eq.'leng' .or. (ztext(1:4).eq.'len '))then
         call inpi(len(j))
         if(len(j).le.0.or.len(j).gt.9999999)call caserr
     *        ('invalid length for new file')
_IF(ga)
      Else If( ztext( 1:4 ) .EQ. 'gafs' ) Then
         iwhere( j ) = 7
         onodid( j ) = .not. opg_root()
_ENDIF
      else 
         call caserr('bad file keyword')
      endif
      go to 80
c
c     =store/core/memo=
c
 50   call inpi(llll)
 51   call inpa4(yjunk)
      if (yjunk.eq.'mwor') then
         llll = llll*1000*1000
      else if (yjunk.eq.'mbyt') then
         llll = llll*1000*1000/8
      else if (yjunk.eq.'gwor') then
         llll = llll*1000*1000*1000
      else if (yjunk.eq.'gbyt') then
         llll = llll*1000*1000*1000/8
      else if (yjunk.eq.'tota') then
         nn = ipg_nnodes()
         if (nn.le.0) nn = 1
         llll = llll/nn
      else if (yjunk.eq.'prin') then
         oprintm = .true.
      else if(yjunk.eq.'debu') then
         call gmem_set_debug(.true.)
      else if(yjunk.eq.'norm') then
         oprintm = .false.
         call gmem_set_debug(.false.)
         call gmem_set_quiet(.false.)
      else if(yjunk.eq.'quie') then
         oprintm = .false.
         call gmem_set_quiet(.true.)
      else if(yjunk.eq.'nona') then
        call gmem_set_nanify(.false.)
      else if(yjunk.eq.'nani') then
        call gmem_set_nanify(.true.)
      else if (yjunk.eq.'    ') then
         go to 90
      else
         call caserr('bad memory directive')
      endif
      go to 51
c
c     =lpage=
c
 60   call inpi(llll)
      llll=llll*65536
 90   lword=llll
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
      if(ovirt) then
       if(iwhere(4).ne.1. or. iwhere(8).ne.1)
     *  call caserr('both ed3 and ed7 must be specified as virtual')
      endif
c
c...  set iwr to 94 to seperate output of nodes
c
      if (.not.opg_root()) then
         iwr = 94
c...      different numbers are not needed / different names are
c...      if no node output is requested send it to sink
         nn = 100+minode
         write(fname,'(a,i3)') 'out_gam_',nn
         if (opkill) oprint(60) = .false.
         if (opnode) then
          open(iwr,form='formatted',file=fname,err=9998,
     +         status='unknown')
         else
          open(iwr,form='formatted',file='/dev/null',err=9998,
     +         status='unknown')
         endif
c
      end if

_IF(dl-find)
c pass units to dl-find module
      call dlf_output(iwr,0)
_ENDIF

      igafs = 0
      do 118 i = 1,maxlfn
         if(ipiomode.eq.IO_A) then
c
c convert iwhere values corresponding to storage on node 0 only
c to replicated storage (by this will affect ed3 and ed7)
c note - file naming on node 0 is not affected (no 000 appended)

c    disk files
            if(iwhere(i).eq.2)iwhere(i)=0
c
c    in core
            if(iwhere(i).eq.5)iwhere(i)=4
c
c    GA files, only required to decide how to dump the files to disk
c    at the end of the run, if required as every node holds a distinct
c    part of these files.
            if(iwhere(i).eq.7)iwhere(i)=6

         endif
         if(ogafsmode) then
c
c  convert disk-based iwhere values to the corresponding GAFS
c  options (note true ##  are not affected)
c  allow for files defined by data as GAfs
            if(iwhere(i).eq.2.or.iwhere(i).eq.6.or.
     +          iwhere(i).eq.7)then
               igafs = igafs + 1
               if (igafs.le.max_ga_files) then
                 oedpr(i) = .true.
                 iwhere(i)=7
               endif
            endif

         endif
118   continue
c
c massage file names
c
c   prepend cfs spec to 
c
      oerr = .false.
      do i = 1, maxlfn
c
         if(lcfsdr .ne. 0 .and. zedfil(i)(1:1).ne.'/')then
c     &        .and. zedfil(i)(1:1).ne.'.')then
            zztemp = zedfil(i)
            call strtrm(zztemp,length)
            zedfil(i)=zcfsdr(1:lcfsdr)//zztemp(1:length)
         endif
c
c add a node identifier - see references to onodid above
c
         zztemp = zedfil(i)
         if(onodid(i))then
            call strtrm(zedfil(i),length)
            call intwrt(zedfil(i),length,ipg_nodeid(),3,ocrap)
         endif

         if(omem(i) .and. len(i) .eq. 0 .and. (i .ne .3) )then
            write(iwr,*)
     &    'error: you must provide a length for memory file ed', i-1
            oerr= .true.
         endif

      enddo
c
      if(oerr)call caserr('missing length spec. for memory file')

      do i = 1, maxfrt
c
         if(lcfsdr .ne. 0 .and. zftfil(i)(1:1).ne.'/')then
            zztemp = zftfil(i)
            call strtrm(zztemp,length)
            zftfil(i)=zcfsdr(1:lcfsdr)//zztemp(1:length)
         endif
c
c add a node identifier (currently all files)
c
         if(i .ne. 58)then
         zztemp = zftfil(i)
ccc         if(iwhere(i).eq.0)then
            call strtrm(zftfil(i),length)
            call intwrt(zftfil(i),length,ipg_nodeid(),3,ocrap)
ccc         endif
         endif
      enddo
c
c open fortran files
c
       call forstr
_IF1(h)       call pg_synch(999)

_IF(ga)
*
*  Set up the GA files
*
      Call ga_file_set_up( maxlfn, iwhere, len, keep, zedfil )
_ENDIF

      if(opun7)write(iwr,201)zinout(2)
 201  format(/' punchfile allocated : lfn = ',a8)

c
      return
9998  write(6,9997)
9997  format(/' *** failed to open output file *** '/)
      call caserr('failure to open output file')
      return
      end
c
c** generic forstr (parallel)
c
      subroutine forstr
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/cfsinfo)
_IF(ipsc)
INCLUDE(common/ipsc)
_ENDIF
c
      do 1 i=1,maxfrt
         if(iposf(i).le.0.and.keepf(i).le.0) go to 1
_IF(ipsc)
         if(ohost(i)) go to 1
_ENDIF
         iunit=i
         call strtrm(zftfil(iunit),length)

         if(oform(i)) then
            open(unit=iunit,iostat=ioerr,
     *           file=zftfil(iunit)(1:length),
     *           status='unknown', form='formatted')
         else
            open(unit=iunit,iostat=ioerr,
     *           file=zftfil(iunit)(1:length),
     *           status='unknown', form='unformatted')
         endif

         if(ooprnt)write(iwr,1000)zftfil(iunit)(1:80)
 1000    format('  opening ',a80)

         if(ioerr.ne.0)then
            write(iwr,3)yft(iunit),zftfil(iunit)
 3          format(1x,'error opening fortran file',1x,a4,
     *           ' : filename : ',a80)
         else
            if(ooprnt)then
               write(iwr,2)yft(iunit),zftfil(iunit)
 2             format(/' fortran stream ',a4,
     *              ' allocated: file name = ',a)
            endif
         endif
_IF(ipsc_io)
c
cPS - is this correct given all that all nodes access
c     different files ?
c
         call setiomode(iunit,0)
_ENDIF
 1    continue

_IF(ipsc)
c temp fix for QM/MM
      if(ohost(1)) call funitm(96,1,zftfil(1))
      if(ohost(3)) call funitm(97,3,zftfil(3))
_ENDIF
      return
      end
_IF(ipsc)
c
c** ipsc funitm
c
      subroutine funitm(ihost,inode,fname)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 rec132
      character *132 fname
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/nodinf)
      character *132 zznam
_IF(cfs)
INCLUDE(common/cfsinfo)
c
c --- root writes a copy of the unit ihost into unit inode on /cfs
c
      call strtrm(fname,lenhst)
      if(lcfsdr .ne. 0)then
         zznam=zcfsdr(1:lcfsdr)//fname(1:lenhst)
      else
         zznam=fname
      endif
      call strtrm(zznam,length)

      if(opg_root()) then
c 
c ... open file on /cfs/ and host
c
         open(unit=ihost,iostat=ioerr,form='formatted',
     +        access='sequential',
     *        status='unknown',file=fname(1:lenhst) )
         open(unit=inode,iostat=ioerr,form='formatted',
     +        access='sequential',
     *        status='unknown',file=zznam(1:length) )
         print *,' just opened ',fname(1:lenhst),'  for node ',minode
 1060    read(ihost,1080,end=1090) rec132
         write(inode,1080)rec132
 1080    format(a132)
         go to 1060
 1090    continue
         rewind inode
      endif
c
      call gsync
c
      if(.not.opg_root()) then
         open(unit=inode,iostat=ioerr,form='formatted',
     +        access='sequential',
     +        status='unknown',file=zznam(1:length) )
         if(ioerr.ne.0) call caserr(
     +        'error opening cfs input file')
         print *,' just opened ',fname,' for node ',minode
         rewind inode
      endif
c     
      call gsync
_ELSE
c
c this is a null routine as the cfs m4 flag was not set
c
_ENDIF
      return
      end
_ENDIF
c
c** generic initb (parallel)
c
      subroutine initb
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 zname
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/disc)
INCLUDE(common/discc)
_IF(ipsc)
INCLUDE(common/ipsc)
_ENDIF
INCLUDE(common/maxlen)

      dimension zattr(8)

INCLUDE(common/cfsinfo)
c
c  print default directory is present
c
      if(lcfsdr .ne. 0)then
         write(iwr,301)zcfsdr
 301     format(//1x,'cfs directory = ',a44//)
      endif

      if(ogafsmode)then
         write(iwr,302)
 302     format(//1x,'Node 0 files routed to GAs')
      endif

      odd=.false.
      do 100 i=1,maxlfn
         if(oedpr(i)) goto 200
 100  continue
      go to 600
 200  write(iwr,300)
 300  format(//1x,
     * 'lfn       external file name',24x,'status    attributes '/
     *     1x,72('*'))
      odd=.true.

      do 400 i=1,maxlfn
         if(oedpr(i))then
_IF(ipsc)
            ztemp = ' '
_ELSE
            inquire (file=zedfil(i),exist=ok)
            if (ok) then
               ztemp = 'old'
            else
               ztemp = 'new'
            endif
_ENDIF
c
c tabulate file attributes
c
            nattr=0
            if(isync(i).eq.0)then
               nattr=nattr+1
               zattr(nattr)='async   '
            else 
               nattr=nattr+1
               zattr(nattr)='sync    '
            endif
            if(.not. omem(i))then
               if(keep(i).eq.1)then
                  nattr=nattr+1
                  zattr(nattr)='keep    '
               else 
                  nattr=nattr+1
                  zattr(nattr)='delete  '
               endif
            endif
            if(iwhere(i).eq.-1)then
               nattr=nattr+1
               zattr(nattr)='via host'
            else if(iwhere(i).eq.0)then
               nattr=nattr+1
               zattr(nattr)='per-node'
            else if(iwhere(i).eq.1)then
               nattr=nattr+1
               zattr(nattr)='vir-host'
            else if(iwhere(i).eq.2)then
               nattr=nattr+1
               zattr(nattr)='node-0  '
            else if(iwhere(i).eq.3)then
               nattr=nattr+1
               zattr(nattr)='shared '
            else if(iwhere(i).eq.4)then
               nattr=nattr+1
               zattr(nattr)='node-mem'
            else if(iwhere(i).eq.5)then
               nattr=nattr+1
               zattr(nattr)='root-mem'
            else if(iwhere(i).eq.6)then
               nattr=nattr+1
               zattr(nattr)='node-ga'
            else if(iwhere(i).eq.7)then
               nattr=nattr+1
               zattr(nattr)='root-ga'
            endif
            
            zname = zedfil(i)
            if(omem(i))zname="<in core>"

            write(iwr,500)yed(i),zname(1:44),ztemp,(zattr(k),k=1,nattr)
 500        format(1x,a4,6x,a44,1x,a8,8(1x,a8))
 501        format(1x,a8,2x,a44,1x,a8,8(1x,a8))
         endif
 400  continue
c
 600  do 700 i=1,maxfrt
         if(oftpr(i))goto 800
 700  continue
      return
c
 800  if(odd) then
         write(iwr,850)
 850     format(/)
      else
         write(iwr,300)
      endif
      do 900 i=1,maxfrt
         if(oftpr(i))then
            
            nattr=0
            if(oform(i))then
               nattr=nattr+1
               zattr(nattr)='formattd'
            endif
_IF(ipsc)
            if(ohost(i))then
               nattr=nattr+1
               zattr(nattr)='host'
            endif
_ENDIF
            write(iwr,501)yft(i),zftfil(i)(1:44),
     &           zftstat(i),(zattr(k),k=1,nattr)
         endif
 900  continue
      return
      end
c
c** generic ed2mem (parallel)
c
      subroutine ed2mem(nbasis,nat)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
INCLUDE(common/atmblk)
INCLUDE(common/iofile)
INCLUDE(common/nodinf)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn)
      common/vinteg/q(1)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
c
c     allocate memory for 2-electron list if available
c
      iqaddr = 0
      iqqoff(3) = 0
      n2 = nbasis*(nbasis+1)/2
      n4 = n2*(n2+1)/2
      iblock = (n4-1)/num2e + 1
c
c ... assume even distribution over nodes
c
c .... allow for j+k supermatrix
      if(nopk.eq.-1) mxblock = mxblock + mxblock
      mxblock = (iblock-1)/nnodes + 1
c .... allow for symmetry skeletonisation
      if (nat.gt.3.and.nt.gt.1) then
        if (nat.gt.5.and.nopk.eq.0) mxblock = mxblock + mxblock
        nttt = nt / 2
        mxblock = mxblock / nttt
      endif

      if(len(3).gt.0.or.maxset(3).gt.0) then
         if (len(3).gt.0) mxblock = (len(3)-1)/nnodes + 1
         if (maxset(3).gt.0) mxblock = maxset(3)
         if(opg_root())write(iwr,*)'ED2 allocation per node: ',
     &        mxblock,' blocks'
      endif
c
c
c
c  added ps to check .....
c
c
      need = mxblock * 512 
      len(3) = mxblock
      call getmem(need, q, iqaddr, iqqoff(3))
c
      if (iqaddr.eq.0.or.iqqoff(3).eq.0) then
       if(n2file.eq.1) then
        if(minode.eq.0) write(iwr,300) need
        call caserr('no memory available for ed2')
        omem(3)=.false.
        iwhere(3) = 0
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
         iwhere(3) = 0
         go to 200
        endif
       endif
      endif
      if(minode.eq.0) write(iwr,400) mxblock,need,mxblock*4/1000.0d0
      omem(3)=.true.
      maxbl(3) = mxblock
      n2last(1)= mxblock
c
200   return
300   format(1x,'requested memory not available (',i12,' words)')
400   format(1x,'sufficient memory allocated for an mfile of',i10,1x,
     +'integral blocks',/,' (',i12,' words / ',f11.3,' Mbyte (/node) )')
      end
c
c** generic ed0mem (parallel)
c
      subroutine ed0mem(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn)
      common/vinteg/q(1)
INCLUDE(common/iofile)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
c
c     allocate memory for dump or scratchfile if available
c
      iqaddr = 0
      iqqoff(iunit) = 0
c      if(maxset(iunit).le.0) call caserr(
      if(len(iunit).le.0) call caserr(
     +     'failure to allocate dumpfile memory')
c     mxblock = maxset(iunit)
      mxblock = len(iunit)
      need = mxblock * 512
      call getmem(need, q, iqaddr, iqqoff(iunit))
      if (iqaddr.eq.0.or.iqqoff(iunit).eq.0) then
         if(opg_root())then
        write(iwr,*)' requested memory not available (',need,' words)',
     1              ' so unit ',iunit,' is on disk '
         endif
            omem(iunit)=.false.
            if(iwhere(iunit).eq.5)then
               iwhere(iunit) = 2
            else if(iwhere(iunit).eq.4)then
               iwhere(iunit) = 1
            else
               call caserr('invalid iwhere in ed0mem')
            endif
         else
            if(opg_root())then
         write(iwr,*)' sufficient memory allocated for ',mxblock,
     +           ' blocks on unit',iunit
               omem(iunit)=.true.
            endif
         endif
      return
      end
c
c** generic put (parallel)
c
      subroutine put(c,nword,iunit)
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
c
      common/bufa/abc(512*maxlfn)
_IF(i8drct)
      integer *8 ipack2, pack2
      equivalence (rpack2,ipack2)
_ENDIF
c
INCLUDE(common/iofile)
INCLUDE(common/discc)
INCLUDE(common/disc)
c
      common/vinteg/q(1)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 iword
_ENDIF
c
INCLUDE(common/timeperiods)
      integer timflag
      dimension timflag(-1:5)
      data timflag/TP_IOPM1,TP_IOP0,TP_IOP1,TP_IOP2,
     &     TP_IOP3,TP_IOP4,TP_IOP5/
c
      dimension c(*)

_IF(ga)
      If( iwhere( iunit ) .EQ. 6 .Or. iwhere( iunit ) .EQ. 7 ) Then 
         Call put_ga( c, nword, iunit )
      Else
_ENDIF
c
      call start_time_period(timflag(iwhere(iunit)))

      isel=iunit
      nstart=(iunit-1)*512+1
      irec=ipos(isel)
      if(irec.lt.1) call iofail(iunit,2)
      ipos(isel)=irec+1
c
      if (iwhere(iunit).eq.1) then
         if (opg_root()) then
c...   virtual io to host // root only
            if (ids48(iunit).ge.0) then
             call msgw(ids48(iunit))
             ids48(iunit)=-1
            end if
            call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
            call pack2(irec,nword,abc(isel*512))
_ELSEIF(i8drct)
            ipack2 = pack2(irec,nword)
            abc(isel*512) = rpack2
_ELSE
            abc(isel*512) = pack2(irec,nword)
_ENDIF
cjvl     handshake
            if (oputpp(iunit)) then
               itype=999990+iunit
               call cprob(itype)
               if(infocount().le.0) then
_IF(ipsc)
                  call crecv(itype,dum,0)
_ELSE
                  call caserr('crecv only on ipsc')
_ENDIF
               else
           print *,' node ',minode,' **** put error on unit ',isel
               endif
cjvl     handshake
            endif
            oputpp(isel) = .true.
c...   send to host memory
            ityp = iunit*10000+irec
            ids48(iunit) = isend(ityp,abc(nstart),4096,mihost,mpid)
         end if
      else if (iwhere(iunit).eq.0 .or.
     +         iwhere(iunit).eq.3 ) then
c...     concurrent put on all
         if (isync(iunit).eq.0) then
c...     asynchronous
            if (ids48(iunit).ge.0) then
             call iowt(ids48(iunit))
             ids48(iunit)=-1
            end if
            call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
            call pack2(irec,nword,abc(isel*512))
_ELSE
            abc(isel*512) = pack2(irec,nword)
_ENDIF
            ids48(iunit) = iwrite(nam(iunit),abc(nstart),4096)
         else
c...     synchronous
            call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
            call pack2(irec,nword,abc(isel*512))
_ELSEIF(i8drct)
            ipack2 = pack2(irec,nword)
            abc(isel*512) = rpack2
_ELSE
            abc(isel*512) = pack2(irec,nword)
_ENDIF
_IF1(h)            call cwrite(nam(iunit),abc(nstart),4096)
_IFN1(h)            call putcc(nam(isel),abc(nstart),ierrio)
_IFN1(h)            if(ierrio.ne.0)call ioerr('write',isel,' ')
         end if
c
      else if (iwhere(iunit).eq.2) then
c
         if (opg_root()) then
c...   cfs io // root only  -
          if (isync(iunit).eq.0) then
            if (ids48(iunit).ge.0) then
             call iowt(ids48(iunit))
             ids48(iunit)=-1
            end if
            call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
            call pack2(irec,nword,abc(isel*512))
_ELSEIF(i8drct)
            ipack2 = pack2(irec,nword)
            abc(isel*512) = rpack2
_ELSE
            abc(isel*512) = pack2(irec,nword)
_ENDIF
            ids48(iunit) = iwrite(nam(iunit),abc(nstart),4096)
          else
c...     synchronous
            call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
            call pack2(irec,nword,abc(isel*512))
_ELSEIF(i8drct)
            ipack2 = pack2(irec,nword)
            abc(isel*512) = rpack2
_ELSE
            abc(isel*512) = pack2(irec,nword)
_ENDIF
_IF1(h)            call cwrite(nam(iunit),abc(nstart),4096)
_IFN1(h)            call putcc(nam(isel),abc(nstart),ierrio)
_IFN1(h)            if(ierrio.ne.0)call ioerr('write',isel,' ')
          end if
         end if
c
      else if (iwhere(iunit).eq.4 .or. iwhere(iunit).eq.5) then
c
c     write to memory
c
         if(irec.gt.len(iunit))then
            if(opg_root())write(iwr,*)'attempt to write block',irec,
     &           'on stream ',yed(iunit),' length allocated is',
     &           len(iunit),' blocks'
            call caserr(
     &        'insufficient memory allocated for memory-mapped file')
         endif
         if(iwhere(iunit).eq.4.or.opg_root())then
            iword = (irec-1)*512+iqqoff(iunit)+1
_IF(bits8)
            call pack2(irec,nword,q(iword+511))
_ELSEIF(i8drct)
            ipack2 = pack2(irec,nword)
            q(iword+511)=rpack2
_ELSE
            q(iword+511)=pack2(irec,nword)
_ENDIF
            call dcopy(nword,c,1,q(iword),1)
         endif
      else
c
c back to host
c
         call dcopy(nword,c,1,abc(nstart),1)
_IF(bits8)
         call pack2(irec,nword,abc(isel*512))
_ELSEIF(i8drct)
         ipack2 = pack2(irec,nword)
         abc(isel*512) = rpack2
_ELSE
         abc(isel*512) = pack2(irec,nword)
_ENDIF
         call synput(abc(nstart),512,nam(iunit),irec)
      end if
c
      call end_time_period(timflag(iwhere(iunit)))
      return

_IF(ga)
      End If
_ENDIF
      end
c
c** generic rdtpnm (parallel) 
c
      subroutine rdtpnm(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical tell
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
INCLUDE(common/disc)
INCLUDE(common/discc)
_IF(ipsc_io)
      integer restrictvol,volist
      dimension volist(8),newlst(8)
_ENDIF
INCLUDE(common/cfsinfo)
_IF(ipsc)
INCLUDE(common/ipsc)
_ENDIF
c
c      nam contains (fortran) iunit number of the stream
c      iwhere:  1 : virtual io to host / 
c               0 : concurrent filesystem (## files)
c      iwhere:  2 : concurrent filesystem  - node 0
c      iwhere:  3 : concurrent filesystem (single file)
c      iwhere:  4 : node memory (## file)
c      iwhere:  5 : node memory (node 0)
c              -1 : io via host
*      iwhere:  6 : GA file, only node 0 dumps file at end of job
*      iwhere:  7 : GA file, all nodes dump their version at the end of the job
c      isync :  1 : synchronous io  / 0 : asynchronous
c
INCLUDE(common/timeperiods)
      integer timflag
      dimension timflag(-1:5)
      data timflag/TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,
     &     TP_IOO3,TP_IOO4,TP_IOO5/

_IF(ipsc_io)
       data volist/0,1,2,3,4,5,6,7/
_ENDIF
      isel=iunit

_IF(ga)
      If( iwhere( iunit ) .EQ. 6 .Or. iwhere( iunit ) .EQ. 7 ) Then 
         Call rdtpnm_ga( iunit )
      Else
_ENDIF

      call start_time_period(timflag(iwhere(iunit)))

c
c dont open memory mapped files
c
      if (iwhere(iunit).eq.5 .or.
     &     iwhere(iunit).eq.4) then
c
c ed2mem is done from master
c
         if(iunit .eq. 3)goto 100
c
c ed0mem allocation now
c
         call ed0mem(iunit)
         go to 1020
      endif
c
c
      if(opg_root().and.ooprnt)write(iwr,1000)zedfil(iunit)(1:80)
 1000 format('  opening ',a80)

      call strtrm(zedfil(iunit),length)      
      if (iwhere(iunit).eq.-1)then
c
c...  open direct access host-io
c
           open(unit=nam(iunit),file=zedfil(iunit)(1:length),
     *          status='unknown',
     *          form='unformatted',access='direct',recl=4096)
c
        else if (iwhere(iunit).eq.1) then
c
c  ...  send open message to host
c
_IF(ipsc)
           call csend(iunit,buf,0,mihost,mpid)
_ELSE
            call caserr('csend only on ipsc')
_ENDIF
c
        else if (iwhere(iunit).eq.2 .and. .not. opg_root()) then
c
c ... skip rest of routine - file is not needed on this node
c
           goto 100
c
        else if(iwhere(iunit).eq.0 .or.
     &          iwhere(iunit).eq.2 .or.
     &          iwhere(iunit).eq.3)then
c
c  ... file i/o 
c
_IF(ipsc)
_IF(ipsc_io)
           numvol = restrictvol(-1,8,volist)
           if(numvol.lt.0) call caserr(
     +          'error detected in restrictvol')
c
_ENDIF
           open(unit=nam(iunit),file=zedfil(iunit)(1:length),
     *          status='unknown',
     *          form='unformatted')
_IF(ipsc_io)
           if(iwhere(iunit).ne.2)call setiomode(nam(iunit),0)
c
c  define defaut size
c  PS? are these still sensible defaults?
c
           if(iwhere(iunit).eq.0)then

              if (len(iunit).eq.0) then
                 mxbl = 5000
              else
                 mxbl = len(iunit) / nnodes
              endif

           else if(iwhere(iunit).eq.2)then

              if (len(iunit).eq.0) then
                 mxbl = 1000
              else
                 mxbl = len(iunit)
              endif

           else if(iwhere(iunit).eq.3)then

              if (len(iunit).eq.0) then
                 mxbl = 5000
              else
                 mxbl = len(iunit)
              endif

           endif

           ibytes = lsize (nam(iunit),mxbl*4096,0)
           if(ibytes.lt.0) call caserr('isize call aborted')
c
           numvol = restrictvol(nam(iunit),-1,newlst)
           if(numvol.gt.0) then
              if (ooprnt) write(iwr,2000) 
     +             zedfil(iunit)(1:80),(newlst(loop),loop=1,numvol)
 2000         format(1x,'file ',a80,/,5x,'allocated to volumes ',8i6)
           else
              call caserr('error in restrictvol')
           end if
_ENDIF
_ELSE
c
c     determine status and length of file
c
           inquire (file=zedfil(iunit),exist=tell)
           if(tell) then
              Call length_file( zedfil( iunit ), lenfil, irep )
              istat(isel)=2
              if(irep.ne.0)call iofail(isel,3)
           endif
_IF(fcio)
           open(nam(isel),file=zedfil(iunit),form='unformatted',
_IFN1(t)     &          access='direct',status='unknown',recl=4096)
_IF1(t)     &          access='direct',status='unknown',recl=1024)
_ELSE
           call opencc(zedfil(iunit),length,nam(isel),ierrio)
           if(ierrio.ne.0)call ioerr('open',isel,zedfil(iunit))
_ENDIF
_ENDIF
      else
        call caserr('invalid iwhere in rdtpnm')
      end if
c
1020  continue
      if(ooprnt)write(iwr,1)yed(isel),lenfil
 1    format(/' lfn ',a4,' allocated: current length =',
     *     i5,' blocks')

 100  continue

      call end_time_period(timflag(iwhere(iunit)))

_IF(ga)
      End If
_ENDIF
      return
      end
c
c** generic delfil (parallel)
c
c   close a file, delete if keep is not set
c
      subroutine delfil(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/iofile)
INCLUDE(common/maxlen)
c
c not open
c
_IF(ga)

      If( iwhere( iunit ) .EQ. 6 .Or. iwhere( iunit ) .EQ. 7 ) Then 
         Call delfil_ga( iunit )
      Else
_ENDIF
      if (ipos(iunit).le.0)return
c
c memory file
c
      if(omem(iunit))return
c
      if (iwhere(iunit).eq.0 ) then
         call waitedn(iunit)
         call shutedn(iunit,1)
      else if (iwhere(iunit).eq.3 ) then
         call waitedn(iunit)
         if(.not. opg_root())call shutedn(iunit,0)
         call pg_synch(123)
         if(opg_root())call shutedn(iunit,1)
      else if (iwhere(iunit).eq.2) then
         if(opg_root()) then
            call waitedn(iunit)
            call shutedn(iunit,1)
         endif
      else
         call caserr('invalid iwhere in delfil')
      end if
      ipos(iunit)=0

      return
_IF(ga)
      End If
_ENDIF
      end
_ENDIF
_IF(fcio)
c
c** fcio srchcc
c
      subroutine srchcc(itape,iblock,ierr)
      common/scum/ipos(40)
      ierr=0
      ipos(itape) = iblock
      end
c
c** fcio getcc
c
      subroutine getcc(itape,b,ierr)
      REAL  b(512)
      common/scum/ipos(40)
      ierr = 0
      iblock = ipos(itape)
c     write(6,*) ' reading from ',itape,' at ',iblock
      read (itape,rec=iblock,err=999) b
      ipos(itape) = ipos(itape) + 1
      return
 999  ierr = 1
      return
      end
c
c** fcio putcc
c
      subroutine putcc(itape,b,ierr)
      REAL  b(512)
      common/scum/ipos(40)
      ierr = 0
      iblock = ipos(itape)
c     write(6,*) ' writing to ',itape,' at ',iblock
      write (itape,rec=iblock,err=999) b
      ipos(itape) = ipos(itape) + 1
      return
 999  ierr = 1
      return
      end
c
c** fcio getccn
c
      subroutine getccn(itape,b,n,ierr)
      REAL  b(512,n)
      do 10 i = 1,n
         call getcc(itape,b(1,i),ierr)
10    continue
      end
c
c** fcio closecc
c
      subroutine closecc(itape)
c     write(6,*) ' closing ',itape
      close(itape)
      end
_ENDIF
_IF(fortio)
c
c** fortio chekwr
c
      subroutine chekwr(iunit)
       implicit REAL  (a-h,o-z)
      integer     await,error_code
      integer     handle,words_io
      logical ochek
INCLUDE(common/sizes)
      common/discne/mlen(maxlfn*2),instat(maxlfn*2)
      common/disc/irep,isel(maxlfn+4),handle(maxlfn),
     * keep(3*maxlfn),io(maxlfn)
INCLUDE(common/maxlen)
c
      iostat=io(iunit)
      ochek = omem(iunit)
      if(iostat.gt.0.and..not.ochek) then
      error_code=    await(handle(iunit),words_io)
      if(error_code.ne.0)then
      irep=error_code
      call iofail(iunit,instat(iunit))
      endif
      endif
      io(iunit)=0
      irep=0
      instat(iunit)=0
      return
      end
c
c** fortio clredx
c
      subroutine clredx
      implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
      integer *8 mask0
INCLUDE(common/sizes)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     * nam(maxlfn)
      common/discne/mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      data mask0/z'0' /
c...
c... to clear up i/o on all internal streams
c...
      kunkpb=mask0
      do 1 mbuff=1,npb
1     call termbf
      return
      end
c
c** fortio fget
c
      subroutine fget(ic,nw,iunit)
      implicit REAL  (a-h,o-z)
      integer *8 ibufa,igetm,iputm,kunkpb,junjpb
      integer *8 lunlpb,mask,mask0,l,lmask
      integer *8 packi8
      integer *8 ic
      dimension ic(*)
INCLUDE(common/sizes)
      common/bufnew/ibufa(512)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     *nam(maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lunx,latmx,lblx,lbfbax,lbuffx
      common/maskc/mask(64)
      data lmask/'8000000000000000'x /
      data mask0/'0'x /
c...
c... synchronous input routine
c...
      lun=iunit
      lbl=ipos(lun)
      lpb=(lbl-1)/npbsz
c      lunlpb=insert(lpb,0,32,lun)
      lunlpb=packi8(lun,lpb)
      lat=lpb*npbsz
      latm=lbl-lat
      if(lunlpb.eq.kunkpb)goto 776
      kunkpb=lunlpb
      mbuff=locatz(junjpb,npb,lunlpb)
      if(mbuff)777,888,777
c... physical block not in buffers
888   mbuff=locat1(ipribf,npb,npbm1)
      igetm(mbuff)=mask0
      lbfbas=ibfbas(mbuff)
      call termbf
      jpbbf(mbuff)=lpb
      junibf(mbuff)=lun
      junjpb(mbuff)=lunlpb
      call revpri
      goto 333
c... physical block in buffer-check for put status
777   call revpri
776   lbfbas=ibfbas(mbuff)
      if(istabf(mbuff))666,444,444
666   mzero=0
      iat=lat
      l=iputm(mbuff)
      m=lbfbas
      call chekwr(lun)
      istabf(iobuff(lun))=0
3       if(l.eq.mask0)goto 4
       nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=ishft(l,nzero)
       mzero=leadz(not(l))
      call chekwr(lun)
      call flushv(lun,ibufa(m+1),mzero,iat)
       l=ishft(l,mzero)
      iat=iat+mzero
      goto 3
4     if(length(lun).lt.iat)length(lun)=iat
      iobuff(lun)=mbuff
      istabf(mbuff)=2
444   l=ishft(igetm(mbuff),latm-1)
      if(iand(l,lmask).eq.lmask)goto 999
c... internal block missing
333   m=length(lun)-lat-npbsz
      iat=npbsz
       istabf(iobuff(lun))=0
      call chekwr(lun)
      if(m)21,20,20
21    iat=iat+m
      if(iat)1,1,20
1     call iofail(lun,1)
20    call fillv(lun,ibufa(lbfbas+1),iat,lat)
      iobuff(lun)=mbuff
      istabf(mbuff)=1
      igetm(mbuff)=mask(iat)
999    ipos(lun)=lbl+1
      if(istabf(mbuff).eq.1)call chekwr(lun)
      m=latm*512+lbfbas
      call upaki8(ibufa(m),nw,iblok)
      if(iblok.ne.lbl)call iofail(lun,1)
      if(nw.gt.0)call vmov(ibufa(m-511),1,ic,1,nw)
      return
      end
c
c** fortio  fillv
c
      subroutine fillv(iunit,bufa,mbl,ibl)
       implicit REAL  (a-h,o-z)
      integer handle,aread
      integer error_code
INCLUDE(common/sizes)
      common/disc/isel(5),ipos(maxlfn),handle(maxlfn),
     *keep(maxlfn*3),io(maxlfn)
      dimension bufa(*)
      common/discne/
     *mlen(maxlfn*2),iostat(maxlfn),iobuff(maxlfn)
      common/vinteg/q(1)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 iword
_ENDIF
      iposw=ibl*512
      nword=mbl*512
      io(iunit)=io(iunit)+1
      if(omem(iunit)) then
       iword = iposw+iqqoff(iunit)+1
       call vmov(q(iword),1,bufa,1,nword)
      else
       iposb=iposw*8
       nbytes=nword*8
       error_code=aread(handle(iunit),iposb,bufa,nbytes)
       if(error_code.gt.0)then
       isel(1)=error_code
       call caserr('i/o error in fillv')
       endif
      endif
      iostat(iunit)=1
      return
      end
c
c** fortio find
c
      subroutine find(iunit)
       implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
      integer *8 mask,mask0,lunlpb,l,lmask
      integer *8 packi8
      integer *8 ibufa
INCLUDE(common/sizes)
      common/bufnew/ibufa(512)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     * nam(maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     * iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     * jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      common/maskc/mask(64)
      data lmask/'8000000000000000'x /
      data mask0/'0'x /
c...
c... initiates asynchronous input
c...
      lun=iunit
      lbl=ipos(lun)
      lpb=(lbl-1)/npbsz
      lunlpb=packi8(lun,lpb)
      lat=lpb*npbsz
      latm=lbl-lat
      if(lunlpb.eq.kunkpb)goto 776
      kunkpb=lunlpb
      mbuff=locatz(junjpb,npb,lunlpb)
      if(mbuff)777,888,777
c... physical block not in buffers
888   mbuff=locat1(ipribf,npb,npbm1)
      igetm(mbuff)=mask0
      lbfbas=ibfbas(mbuff)
      call termbf
      jpbbf(mbuff)=lpb
      junibf(mbuff)=lun
      junjpb(mbuff)=lunlpb
      call revpri
      goto 333
c... physical block in buffer-check for put status
777   call revpri
776   lbfbas=ibfbas(mbuff)
      if(istabf(mbuff))666,444,444
666   mzero=0
      iat=lat
      l=iputm(mbuff)
      m=lbfbas
      call chekwr(lun)
      istabf(iobuff(lun))=0
3       if(l.eq.mask0)goto 4
       nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=ishft(l,nzero)
       mzero=leadz(not(l))
      call chekwr(lun)
      call flushv(lun,ibufa(m+1),mzero,iat)
       l=ishft(l,mzero)
      iat=iat+mzero
      goto 3
4     if(length(lun).lt.iat)length(lun)=iat
      iobuff(lun)=mbuff
      istabf(mbuff)=2
444   l=ishft(igetm(mbuff),latm-1)
      if(iand(l,lmask).eq.lmask)goto 999
c... internal block missing
333   m=length(lun)-lat-npbsz
      iat=npbsz
       istabf(iobuff(lun))=0
      call chekwr(lun)
      if(m)21,20,20
21    iat=iat+m
      if(iat)1,1,20
 1     call iofail(lun,1)
20    call fillv(lun,ibufa(lbfbas+1),iat,lat)
      iobuff(lun)=mbuff
      istabf(mbuff)=1
      igetm(mbuff)=mask(iat)
999   lbuff=mbuff
       ipos(lun)=lbl+1
      return
      end
c
c** fortio flushv
c
      subroutine flushv(iunit,bufa,mbl,ibl)
      implicit REAL  (a-h,o-z)
      integer error_code,awrite
      integer handle
INCLUDE(common/sizes)
      dimension bufa(*)
      common/discne/
     *mlen(maxlfn*2),iostat(maxlfn),iobuff(maxlfn)
      common/disc/isel(5),ipos(maxlfn),handle(maxlfn),
     *keep(3*maxlfn),io(maxlfn)
      common/vinteg/q(1)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 iword
_ENDIF
      iposw=ibl*512
      nword=mbl*512
      io(iunit)=io(iunit)+1
      if(omem(iunit)) then
       iword = iposw+iqqoff(iunit)+1
       call vmov(bufa,1,q(iword),1,nword)
      else
       nbytes=nword*8
       iposb=iposw*8
       error_code=awrite(handle(iunit),iposb,bufa,nbytes)
       if(error_code.gt.0)then
       isel(1)=error_code
       call caserr('i/o error in flushv')
       endif
      endif
      iostat(iunit)=2
      return
      end
c
c** fortio get
c
      subroutine get(ic,nw)
       implicit REAL  (a-h,o-z)
      integer *8 ibufa,igetm,iputm,kunkpb,junjpb
      integer *8 ic
INCLUDE(common/sizes)
      dimension ic(*)
      common/bufnew/ibufa(512)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     * nam(maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      common/maskc/mask(64)
c...
c... completes asynchronous input
c...
      if(istabf(lbuff).eq.1)call chekwr(lun)
      m=latm*512+lbfbas
      call upaki8(ibufa(m),nw,iblok)
      if(iblok.ne.lbl)call iofail(lun,1)
      if(nw.gt.0)call vmov(ibufa(m-511),1,ic,1,nw)
      return
      end
c
c** fortio inita
c
      subroutine inita(lword,nbuff,ibuff)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      character*16 logfil
INCLUDE(common/sizes)
      dimension ytext(15),zinout(2)
      common/clogfl/logfil
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
      common/discne/mlen(maxlfn),length(maxlfn)
INCLUDE(common/work)
INCLUDE(common/timez)
INCLUDE(common/segm)
INCLUDE(common/disc)
INCLUDE(common/discc)
      common/bufnew/bufa(512,maxblo,maxbuf)
      common/blksiz/nsz,nsz512(8),nsznav,nsstat,
     *  nslen,nsmax ,nsort, nstid
      common/blksi2/ns2,ns2512(8),npznav,npstat,
     *  n2len,n2max ,nport, nptid
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     *  n3len,n3max ,ntort, nttid
      common/restri/nfile(63),lda(508),isect(508),ldx(508),
     *  iacsct(508)
INCLUDE(common/atmblk)
INCLUDE(common/maxlen)
INCLUDE(common/prints)
      dimension yfname(maxfrt)
      dimension yynam(maxlfn),zdef(maxfrt)
c
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
      data ytext/
     + 'chan','swit','file','seti','memo','stor','core',
     * 'lpag','time','sort','psor','driv','max','logf',
     * 'tabl' /
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
      oprintm = .false.
_IF(debug)
      call gmem_set_nanify(.true.)
_ELSE
      call gmem_set_nanify(.false.)
_ENDIF
      call gmem_set_debug(.false.)
      call gmem_set_quiet(.true.)
      ichek=0
      isel= 0
      iselr= 0
      iselw= 0
      irep= 0
      do 4000 loop=1,maxlfn
      omem(loop)=.false.
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
      do 3500 loop=1,maxlfn
      iop(loop)  = 0
      keep(loop) = 0
      oform(loop) = .false.
      istat(loop) = 1
      ipos(loop) = 0
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
      isecs=0
      nbuff=maxbuf
      ibuff=maxblo
c
      call setsto(maxlfn,0,length)
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
       call getenv(yed(i),filatt)
       if(filatt.ne.' ') keep(i)=1
c...   change properties based on name (nzo for now)
       if (index(filatt,'nzo').ne.0.or.index(filatt,'cfs-0').ne.0.or.
     1     index(filatt,'NZO').ne.0.or.index(filatt,'CFS-0').ne.0) then
          iwhere(i)=2
          onodid(i) = .not. opg_root()
       end if
 115  continue
c
c ---- and repeat for the fortran files
c
      do 116 i=1,maxfrt
       call getenv(zftn(i),filatt)
       if(filatt.ne.' ') keepf(i)=1
 116  continue
c
c --- take a copy of input file
c
      call print0
      oecho = .false.
      opark = .false.
      oreact = .false.
      odci = .false.
      odlfi = .false.
      if(opg_root()) then
       write(iwr,500)
 500   format(/1x,'**********************************'/
     +         1x,'**** GAMESS-UK Input Data Set ****'/
     +         1x,'**********************************'/)
       oecho = .true.
      endif
c
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
300   jrec=jrec-1
 310  call prinit(oreact,opark,odci,odlfi)
      go to 320
 10   call input
 320  call inpa4(yyy)
      j=locatc(ytext,15,yyy)
      if(j)20,30,20
 20   go to (40,40,40,120,50,50,50,60,70,200,240,210,250,900,285)
     *,j
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
      call inpan(logfil)
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
       call inpa(ztext)
       go to 10
c
c  psort
c
 240   call inpi(ns2)
       call inpa(ztext)
       go to 10
c
c  table
c
 285   call inpi(ns3)
       call inpa(ztext)
       go to 10
c
c    setio
c
 120   call inpi(nbuff)
       call inpi(ibuff)
      go to 10
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
c...  512 word files
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
      go to 10
 70   call inpf(ffmin)
      insec=nint(ffmin*60.0d0)
      otime=.true.
      go to 10
c
 30   jrec=0
      if(isecs.le.0)  isecs = 3600000
      if(otime)      isecs=insec
      if(nsz.le.0.or.nsz.gt.10)nsz=10
      if(ns2.le.0.or.ns2.gt.1)ns2=1
      if(ns3.le.0.or.ns3.gt.1)ns3=1
      if(nbuff.le.0)nbuff=maxbuf
      if(ibuff.le.0)ibuff=maxblo
c
       call forstr
c
      if(opun7)write(iwr,201)zinout(2)
201   format(/
     *' punchfile allocated : lfn = ',a8)
      return
      end
c
c** fortio initb
c
      subroutine initb
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/discc)
c
      odd=.false.
      do 100 i=1,maxlfn
      zied=yed(i)
       call getenv(zedfil(i),filatt)
      if(zied.ne.zedfil(i).or.filatt.ne.' ')  go to 200
 100  continue
      go to 600
 200  write(iwr,300)
 300  format(//1x,
     *'local     ',2x,'status',2x,'external file name'/
     *1x,'file name '/1x,80('*'))
      odd=.true.
      do 400 i=1,maxlfn
      zied=yed(i)
      filatt=zedfil(i)
      if(zied.eq.zedfil(i))then
      call getenv(zedfil(i),filatt)
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
 600  do 700 i=1,maxfrt
      zied=zftn(i)
      call getenv(zied,filatt)
      if(zied.ne.zftfil(i).or.filatt.ne.' ')  go to 800
 700  continue
      return
 800  if(odd) then
      write(iwr,850)
 850  format(/)
      else
      write(iwr,300)
      endif
      do 900 i=1,maxfrt
      zied=zftn(i)
      filatt=zftfil(i)
      if(zied.eq.zft(i))then
        call getenv(zftn(i),filatt)
        if(filatt.eq.' ') go to 900
      endif
      write(iwr,500)yft(i),zftstat(i),filatt
 900  continue
      return
      end
c
c** fortio ed2mem
c
      subroutine ed2mem(nbasis,nat)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
INCLUDE(common/atmblk)
INCLUDE(common/iofile)
      common/vinteg/q(1)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
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
       need = mxblock * 512
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
300   format(1x,'requested memory not available (',i6,' words')
400   format(1x,'sufficient memory allocated for an mfile of',i10,1x,
     +'integral blocks',/,' (',i12,' words / ',f11.3,' Mbyte )')
      end
c
c** fortio ed0mem
c
      subroutine ed0mem(iunit)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
      common/vinteg/q(1)
INCLUDE(common/iofile)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
_IF(64bitpointers)
      integer*8 iqaddr
_ENDIF
c
c     allocate memory for dump or scratchfile if available
c
      iqaddr = 0
      iqqoff(iunit) = 0
      if(maxset(iunit).le.0) call caserr(
     + 'failure to allocate dumpfile memory')
      mxblock = maxset(iunit)
      need = mxblock * 512
      call getmem(need, q, iqaddr, iqqoff(iunit))
      if (iqaddr.eq.0.or.iqqoff(iunit).eq.0) then
       write(iwr,*)' requested memory not available (',need,' words)',
     1             ' so unit ',iunit,' resides on disk'
       omem(iunit)=.false.
      else
        write(iwr,*)' sufficient memory allocated for ',maxset(iunit),
     + ' 512 word blocks (unit =',iunit,')'
       omem(iunit)=.true.
      endif
      return
      end
c
c** fortio initj
c
      subroutine initj(lword)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF(convex,sgi,sun,dec,ksr)
      integer hostnm
      character*44 host
      character*256 com
_ENDIF
_IF(hp700,hpux11)
      integer gethostname
      character*44 host
_ENDIF
_IF(rs6000)
      integer gethostname
      character*44 host
      character*256 com
_ENDIF
INCLUDE(common/timez)
INCLUDE(common/jinfo)
INCLUDE(common/iofile)
INCLUDE(common/segm)
INCLUDE(common/utilc)
INCLUDE(common/timeperiods)
INCLUDE(common/blksiz)
c
c ----- initialise job
c
_IF(signals)
c
c  ===  UNIX signals for flushing o/p (iflsh6) and
c  ===  time dumping (idumpt)
c
      call catch
_ENDIF
      call start_time_period(TP_ENTIRE)
      call inita(lword,nbuff,ibuff)
      call tidajt(zdate,ztime,zyear,zaccno,zanam,isecs)
      call getcor(lword,mfree)
      if(.not.ologf)go to 3
      close(unit=101,status='keep')
      call sptchk('process logfile directive')
3     ti=0.0d0
      tx=0.0d0
      timlim=isecs
      call flags
_IF(convex,sgi,sun,dec,ksr)
      idum=hostnm(host)
      call getarg(0,com)
      write(iwr,10) host, com
 10   format(
     *40x,'Hostname  : ',a44/
     *40x,'GAMESS-UK Executable: ',a70)
_ENDIF
_IF(hp700,hpux11)
      idum=gethostname(host)
      write(iwr,10) host
 10   format(
     *40x,'Hostname  : ',a44)
_ENDIF
_IF(rs6000)
      idum=gethostname(host)
      call getarg(0,com)
      write(iwr,10) host, com
 10   format(
     *40x,'Hostname  : ',a44/
     *40x,'GAMESS-UK Executable: ',a70)
_ENDIF
      write(iwr,20)zanam,zdate(1:6),zyear(1:4),ztime,zaccno,
     &           isecs,lword,lword*8/(1000*1000)
 20   format(
     *40x,'job name   ',a8/
     *40x,'date    ',a6,1x,a4/
     *40x,'time       ',a8/
     *40x,'acct       ',a8//
     *40x,'job time specified',i7,' seconds'//
     *40x,'main store requested',i10,' words',' =',i6,' mbytes'/)
c
c    data set initialisation
c
      call initb
      call setio(nbuff,ibuff)
c
      return
      end
c
c** fortio flags
c
      subroutine flags
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
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
     *35x,'*               ===  CONVEX  version 7.0   ===         ',
     *'        *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*******************************************************',
     *'*********'/)
      return
      end
c
c** fortio print0
c
      subroutine print0
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *10 chpid
      character *132 rec132
c      integer pid,getpid
       logical o_inclin,o_inc_err
INCLUDE(common/sizes)
INCLUDE(common/work)
INCLUDE(common/iofile)
c
c --- write a copy of the input to unit 121
c
c      pid=getpid()
c      write(chpid,'(i10)')pid
c      do 1020 i=1,10
c      if(chpid(i:i).ne.' ') then
c      length=i
c      go to 1030
c      endif
c1020  continue
c1030  continue
      o_inc_err = .false.
      open(unit=121,form='formatted',status='scratch')
1060  if (.not.o_inclin) read(ird,1080,end=1090) rec132
cjvl
c     include statement (allows including files in the input )
c
      call inclin(rec132,122,iwr,o_inclin,o_inc_err)
      write(121,1080)rec132
cjvl
1080  format(a132)
      go to 1060
1090  ird=121
      rewind ird
      return
      end
1060  read(ir,1080,end=1090) rec80
 1080 format(a80)
 1090 continue
      ir=121
      rewind ir
      return
      end
c
c** fortio forstr
c
      subroutine forstr
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character *132 filatt
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
      common/disc/isel(5+maxlfn*8),
     *ipos(maxfrt),keep(maxfrt),istatf(maxfrt),
     *lrecf(maxfrt),nrecf(maxfrt),oform(maxfrt)
c
      do 1 i=1,maxfrt
      if(ipos(i).le.0.and.keep(i).le.0) go to 1
      filatt=' '
      iunit=i
      call getenv(zftfil(iunit),filatt)
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
     * ' : filename : ',a80)
      endif
      if(ooprnt) write(iwr,2)yft(i),filatt
 2    format(/' fortran stream ',a4,' allocated: file name = ',a)
 1    continue
      end
c
c** fortio iofail
c
      subroutine iofail(irw,itype)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn)
INCLUDE(common/iofile)
INCLUDE(common/discc)
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
c** fortio kilbuf
c
      subroutine kilbuf
      implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
      integer *8 mask0
INCLUDE(common/sizes)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      data mask0/'0'x /
c...
c... to clear up all io on internal streams
c... and disconnect buffers
c...
      call clredx
      call vfill(mask0,junjpb,1,npb)
      return
      end
c
c** fortio open
c
      subroutine open(iunit,istat,handle,nrec,length,filnam)
      implicit REAL  (a-h,o-z)
      integer handle
      integer stat,statb(20)
      integer error_code,aset
      logical tell
c
      character *(*) filnam
c
c     open file
c
      inquire (file=filnam,exist=tell)
      if(.not.tell) then
      open (unit=handle,file=filnam,status='unknown')
      error_code = aset(handle,1)
      if(error_code.ne.0)then
      call caserr('error creating/opening file')
      endif
      else
      open (unit=handle,file=filnam,status='unknown')
      error_code = aset(handle,1)
      if(error_code.ne.0) then
      call caserr('error opening file')
      endif
      error_code=stat(filnam,statb)
      istat=2
      length = (statb(7)/8-1)/512+1
      if(error_code.ne.0)
     * call caserr('error in file status')
      endif
      return
      end
c
c** fortio packi8
c
      function packi8(nw,m)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *8 packi8,pack
      dimension  i4(2)
      equivalence (i4(1),pack)
      i4(1)=nw
      i4(2)=m
      packi8 = pack
      return
      end
c
c** fortio put
c
      subroutine put(ic,nw,iunit)
       implicit REAL  (a-h,o-z)
      integer *8 ibufa,igetm,iputm,kunkpb,junjpb
      integer *8 mask0,mask1,masker,l,iunipb
      integer *8 packi8
      integer *8 ic,icompl
INCLUDE(common/sizes)
      dimension ic(*)
      common/bufnew/ibufa(512)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     *nam(maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      data mask0,mask1/'0'x,'1'x /
c...
c... drives asynchronous output
c...
      iun=iunit
      ibl=ipos(iun)
      ipb=(ibl-1)/npbsz
      iunipb=packi8(iun,ipb)
      iat=ipb*npbsz
      iatm=ibl-iat
      masker=ishft(mask1,64-iatm)
      if(iunipb.eq.kunkpb)goto 776
      kunkpb=iunipb
      mbuff=locatz(junjpb,npb,iunipb)
      if(mbuff)777,888,777
c... physical block not in buffer -- find lowest priority buffer
888   mbuff=locat1(ipribf,npb,npbm1)
      igetm(mbuff)=mask0
      call revpri
      call termbf
      jpbbf(mbuff)=ipb
      junibf(mbuff)=iun
      junjpb(mbuff)=iunipb
      goto 555
777   call revpri
c... physical block in buffer - check for i/o activity
776   if(istabf(mbuff))666,555,444
444   iobuff(iun)=maxbuf+1
      call chekwr(iun)
555   iputm(mbuff)=mask0
      istabf(mbuff)=-1
666   m=ibfbas(mbuff)
      k=iatm*512+m
      ibufa(k)=packi8(nw,ibl)
      if(nw.gt.0)call vmov(ic,1,ibufa(k-511),1,nw)
      l=ior(iputm(mbuff),masker)
      iputm(mbuff)=l
      igetm(mbuff)=ior(igetm(mbuff),masker)
      if(iatm.ne.npbsz)goto 999
c... switch to out status
      mzero=0
      istabf(iobuff(iun))=0
      call chekwr(iun)
3      if(l.eq.mask0)goto 4
       nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=ishft(l,nzero)
      mzero=leadz(not(l))
      call chekwr(iun)
      call flushv(iun,ibufa(m+1),mzero,iat)
       l=ishft(l,mzero)
      iat=iat+mzero
      goto 3
4     istabf(mbuff)=2
      if(length(iun).lt.iat)length(iun)=iat
      iobuff(iun)=mbuff
999   ipos(iun)=ibl+1
      return
      end
c
c** fortio quit
c
      subroutine quit(ia,n)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *80 ia
      call whtps
      call shut
      call err(n)
      return
      end
c
c** fortio rdtpnm
c
      subroutine rdtpnm(iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character *132 filatt
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),nam(maxlfn)
     * ,keep(maxlfn),istat(maxlfn),
     * len(maxlfn*2)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
      common/discne/mlen(maxlfn),length(maxlfn)
INCLUDE(common/discc)
INCLUDE(common/maxlen)
      if(omem(iunit)) then
       if(iunit.ne.3) call ed0mem(iunit)
      else
       filatt = zedfil(iunit)
       if(yed(iunit).eq.zedfil(iunit)) then
        call getenv(zedfil(iunit),filatt)
        if(filatt.eq.' ') filatt = zedfil(iunit)
       endif
       call open(iunit,istat(iunit),nam(iunit),len(iunit),
     *  length(iunit),filatt)
       if(ooprnt)write(iwr,1)yed(iunit),length(iunit)
  1     format(/' lfn ',a4,' allocated: current length =',
     *  i5,' blocks')
      endif
      return
      end
c
c** fortio revpri
c
      subroutine revpri
      implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
INCLUDE(common/sizes)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      common/disc/irep(2),npb
       ix=ipribf(mbuff)
       do 1 loop=1,npb
       if(ix.gt.ipribf(loop))ipribf(loop)=ipribf(loop)+1
   1   continue
       ipribf(mbuff)=0
      return
      end
c
c** fortio search
c
      subroutine search(iblock,iunit)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     *name(maxlfn)
      if(iblock.le.0)call iofail(iunit,3)
      if(ipos(iunit).le.0)call rdtpnm(iunit)
      ipos(iunit)=iblock
      return
      end
c
c** fortio setio
c
      subroutine setio(lpb,lpbsz)
      implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
      integer *8 mask,mask0,mask1
INCLUDE(common/sizes)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     *nam(maxlfn),keep(6*maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      common/maskc/mask(64)
c------------internal stream data---------------------------------------
c   ipos(maxlfn)  =  block position
c   nam(maxlfn)   =  lfn  file handle
c   mlen(maxlfn)  =  max file length in 512 word blocks
c   length(maxlfn)=  current length of file in 512 word blocks
c   iostat(maxlfn)=  null/0 or input/1 or output/2 flag for physical io
c   iobuff(maxlfn)=  store buffer number or 34 if no physical io
c   itab(maxlfn)  =  list of file handles
c--------------buffer status data---------------------------------------
c   istabf(maxbuf+1)=  put/-1  nothing/0  input/1  output/2
c   ipribf(maxbuf) = buffer priority 0=maximum npbm1=minimum
c   junibf(maxbuf) = internal stream number
c   jpbbf(maxbuf)  = physical block number
c   junjpb(maxbuf) = junibf<32.or.jpbbf
c   ibfbas(maxbuf) = buffer base point = 0, npbsz*512, 2*npbsz*512,....
c   igetm(maxbuf)  = buffer status flags 1/block present 0/block absent
c   iputm(maxbuf)  = put flags       1/block destined for output
c----------------scalar variables---------------------------------------
c   irep   =reply word for physical io
c   npb    =number of store buffers
c   npbsz  =number of 512 word blocks in one buffer
c   npbm1  =npb-1
c   mbuff  =current buffer number
c   kunkpb =internal stream/physical block number of current buffer
c   lun    =internal stream in find/get
c   latm   =relative 512 word block number in find/get
c   lbl    =absolute 512 word block number in find/get
c   lbfbas =buffer base point in find/get
c   lbuff  =buffer number in find/get
c-----------------------------------------------------------------------
      data mask0,mask1/'0'x,'1'x /
      irep=0
      kunkpb=mask0
      npb=lpb
      npbsz=lpbsz
      npbm1=npb-1
      mask(1)=ishft(mask1,63)
      do 2 i=2,64
 2    mask(i)=ior(ishft(mask1,64-i),mask(i-1))
c      do 100 i=1,64
c      iii=leadz(mask(i))
      call setstz(npb,0,junjpb)
      call setsto(npb,0,istabf)
      call ibasgn(npb,npbm1,-1,ipribf)
      call ibasgn(npb,0,npbsz*512,ibfbas)
      call setsto(maxlfn,0,iostat)
      call setsto(maxlfn,0,ipos)
      call setsto(maxlfn,(maxbuf+1),iobuff)
      call icopy(maxlfn,nam,1,itab,1)
      return
      end
c
c** fortio shut
c
      subroutine shut
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer error_code
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     *nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn),
     *iwhere(maxlfn),isync(maxlfn),
     *iposf(maxfrt),keepf(maxfrt)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/maxlen)
      call clredx
      error_code = 0
      do 1 i=1,maxlfn
      if(ipos(i).gt.0)then
        if(.not.omem(i)) then
          if(keep(i).gt.0) then
          close(unit=nam(i),iostat=error_code,status='keep')
          istat(i)=2
          else
          close(unit=nam(i),iostat=error_code,status='delete')
          istat(i)=0
          endif
        endif
        ipos(i)=0
        if(error_code.ne.0)then
        write(iwr,5)yed(i)
 5      format(1x,'***** error closing/deleting file',1x,a4)
           else
        if(ooprnt)write(iwr,3)yed(i)
 3      format(//' file ',a4,' closed')
        endif
      endif
 1    continue
c ... close fortran files
      do 10 i=1,maxfrt
      if(iposf(i).gt.0) then
      if(keepf(i).gt.0) then
         close(unit=i,iostat=error_code,status='keep')
          else
         close(unit=i,iostat=error_code,status='delete')
          endif
         if(error_code.ne.0)then
          write(iwr,5)zftfil(i)
         else
          if(ooprnt)write(iwr,30)zftfil(i)
 30       format(//' fortran file ',a6,' closed')
         endif
      endif
 10   continue
      close(unit=121,status='delete')
      return
      end
c
c** fortio delfil
c
      subroutine delfil(i)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer error_code
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     *nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn)
INCLUDE(common/utilc)
INCLUDE(common/iofile)
INCLUDE(common/discc)
      common/discne/mlen(maxlfn),length(maxlfn),iostat(maxlfn)
      if(ipos(i).gt.0.and.iostat(i).eq.0)then
          call clredx
          if(keep(i).gt.0) then
          close(unit=nam(i),iostat=error_code,status='keep')
          istat(i)=2
          else
          close(unit=nam(i),iostat=error_code,status='delete')
          istat(i)=0
          endif
        if(error_code.ne.0)then
        write(iwr,5)yed(i)
 5      format(1x,'***** error closing/deleting file',1x,a4)
           else
        if(ooprnt)write(iwr,3)yed(i)
 3      format(//' file ',a4,' closed')
        endif
        ipos(i)=0
        length(i)=0
      endif
      return
      end
c
c** fortio shut1
c
      subroutine shut1(i)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer error_code
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     *nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
      common/discne/mlen(maxlfn),length(maxlfn),iostat(maxlfn)
INCLUDE(common/maxlen)
      if(ipos(i).gt.0.and.iostat(i).eq.0)then
      if(.not.omem(i)) then
        call clredx
        close(unit=nam(i),iostat=error_code,status='keep')
        if(error_code.ne.0)then
         write(iwr,5)yed(i)
 5       format(1x,'***** error closing/deleting file',1x,a4)
        else
         if(ooprnt)write(iwr,3)yed(i)
 3       format(//' file ',a4,' closed')
        endif
      endif
        ipos(i)=0
        length(i)=0
        istat(i)=2
      endif
      return
      end
c
c** fortio termbf
c
      subroutine termbf
      implicit REAL  (a-h,o-z)
      integer *8 igetm,iputm,kunkpb,junjpb
      integer *8 l,mask0
      integer *8 ibufa
INCLUDE(common/sizes)
      common/bufnew/ibufa(512)
      common/disc/irep,incr,npb,npbsz,npbm1,ipos(maxlfn),
     *nam(maxlfn)
      common/discne/
     *mlen(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn),
     *istabf(maxbuf+2),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *kunkpb,mbuff,lun,latm,lbl,lbfbas,lbuff
      data mask0/'0'x /
c...
c... to clear up all pending io for current buffer
c...
      k=istabf(mbuff)
      if(k)8,9,8
8     j=junibf(mbuff)
      call chekwr(j)
      if(k)7,6,6
c... buffer in put mode
7     mzero=0
      istabf(mbuff)=0
      l=iputm(mbuff)
      m=ibfbas(mbuff)
      k=jpbbf(mbuff)*npbsz
3     if(l.eq.mask0)goto 4
       nzero=leadz(l)
      k=k+nzero
      m=(nzero+mzero)*512+m
      l=ishft(l,nzero)
      mzero=leadz(not(l))
      call flushv(j,ibufa(m+1),mzero,k)
       l=ishft(l,mzero)
      k=k+mzero
      call chekwr(j)
      goto 3
4     if(length(j).lt.k)length(j)=k
6     istabf(iobuff(j))=0
      iobuff(j)=maxbuf+1
9     return
      end
c
c** fortio upaki8
c
      subroutine upaki8(b,left,iright)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *8 b,unpack
      dimension i4(2)
      equivalence (i4(1),unpack)
      unpack=b
      left=i4(1)
      iright=i4(2)
      return
      end
c
c** fortio whtps
c
      subroutine whtps
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/discc)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *iblksz(maxlfn)
INCLUDE(common/blksiz)
      common/blksi2/ns2(12),n2len,n2max
      common/blksi3/ns3(11),n3len,n3max
      common/discne/mlen(maxlfn),length(maxlfn)
INCLUDE(common/iofile)
c...
      data ztexts,ztextp,ztextt/'sort ','psort','table'/
      do 3 k=1,maxlfn
      if(ipos(k))3,3,4
3     continue
      return
4     write(iwr,5000)
 5000  format(/1x,'file positions'/
     * 1x,'lfn ','   block  length'/1x,20('='))
      do 1 i=k,maxlfn
      if(ipos(i))1,1,2
2     write(iwr,6000)yed(i),ipos(i),length(i)
1     continue
6000  format(1x,a4,2i8 )
      if(nsmax.gt.0)write(iwr,6001)ztexts,nslen,nsmax
      if(n2max.gt.0)write(iwr,6001)ztextp,n2len,n2max
      if(n3max.gt.0)write(iwr,6001)ztextt,n3len,n3max
6001  format(1x,a6,i7,i8)
      return
      end
_ENDIF
_ENDIF
_ENDIF
_EXTRACT(bitmanip,win95)
_IF(i8drct)
c   popcnt and leadz i8 variants for Intel ifc compiler
c   this is at present CRAP, with the ishft intrinsic
c   returning wrong values when n=64. It took a LONG time 
c   to nail problems encountered in popcnt and leadz
c   note that the littleendian flag is included to
c   permit testing on non-intel platforms
      integer function popcnt(b)
c
c...   64-bit popcount
c...   cray-imitation for direct for 32 bit machines
c...   M4 flags user were GEN_OPTIONS
c
      implicit REAL (a-h,o-z),integer(i-n)
c
      integer *8 a, b
      integer*4 itable(0:15),ii(2)
      integer*4 mask
      parameter (mask=15)
      equivalence(a,ii(1))
      data itable/ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4/
c
      a  =  b 
      kk = 0  
c
      do 100 i=1,2
         do 10 j=1,8
_IF(littleendian)
            jj = IAND32(ishft(ii(3-i),-(j-1)*4),mask)
_ELSE
            jj = IAND32(ishft(ii(i),-(j-1)*4),mask)
_ENDIF
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk
c
      return  
      end     

      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines 
c      
      implicit REAL (a-h,o-z),integer(i-n)
c
      integer*8 a,b
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
            kk = IAND32(ishft(jj,mm4),mask)
            if (kk.ne.zero) then
              leadz= (i-1)*16 + itable(kk)  + koff  
            else
              leadz= (i-1)*16 + itable(jj) + koff + 4
            end if
            return
         end if
100   continue
c
      return  
      end 
      function shiftr(x,n)
c
c...   integer *8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word !! 
c...   lshft or rshft by 32 does not empty word !! 
c
      implicit REAL  (a-h,o-z), integer(i-n)
      integer *8 shiftr
      integer *8 x,a,b
c
      integer *4 ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000) n
1000   format(' shiftr called wrongly',i20)
       call caserr('invalid argument in shiftr')
      endif   
c
      a = x   
c
      if (n.eq.0) then
         b = a
      else if (n.eq.32) then
         ib(1) = ia(2) 
         ib(2) = 0
      else if (n.lt.32) then
_IF(littleendian)
         ib(1) = IOR32(ishft(ia(1),-n),ishft(ia(2),32-n))
         ib(2) = ishft(ia(2),-n)
_ELSE
         ib(2) = IOR32(ishft(ia(2),-n),ishft(ia(1),32-n))
         ib(1) = ishft(ia(1),-n)
_ENDIF
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else    
_IF(littleendian)
         ib(1) = ishft(ia(2),32-n)
         ib(2) = 0
_ELSE
         ib(2) = ishft(ia(1),32-n)
         ib(1) = 0
_ENDIF
      end if  
c
      shiftr = b
c
      return  
      end     
      function shiftl(x,n)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !! 
c
      implicit REAL  (a-h,o-z), integer(i-n)
      integer *8 shiftl,x,a,b
c
      integer *4 ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl called wrongly')
       call caserr('invalid argument in shiftl')
      endif   
c
      a = x   
      if (n.eq.0) then
         b = a
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = ishft(ia(1),n)
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ishft(ia(2),n)
_ENDIF
      else if (n.eq.32) then
         ib(2) = ia(1) 
         ib(1) = 0
      else    
_IF(littleendian)
         ib(2) = ishft(ia(1),n-32)
         ib(1) = 0
_ELSE
         ib(1) = ishft(ia(2),n-32)
         ib(2) = 0
_ENDIF
      end if  
c
      shiftl = b
c
      return  
      end 
      function shift(x,nn)
c
c...   real*8 shift (with wraparound)
c...   n <0 shifts aborts
c...   n >64 is OK now (used in direct) 
c...   lshft or rshft by 32 does not empty word !! 
c
      implicit REAL (a-h,o-z), integer(i-n)
c
      integer*8 shift, a, b, x
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (nn.lt.0) then
       write(6,1000) nn
1000   format(' shift called wrongly',i20)
       call caserr('invalid argument in shift8')
      endif   
c
      n = IAND32(nn,63)
c
      a = x   
      if (n.eq.0.or.n.eq.64) then
         b = a
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
_ELSE  
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
_ENDIF
      else if (n.eq.32) then
         ib(1) = ia(2) 
         ib(2) = ia(1) 
      else    
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
_ELSE  
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
_ENDIF
      end if  
c
      shift = b
c
      return
      end 
_ENDIF(EXTRACT bitmanip win95)

_IFN(i8,integer64,i8drct)
************
c****************************************************************
c************ bitmanipulation routines                ***********
c****************************************************************
_IFN(bits8)
      function pad(icon)
c
c     pack2 with lefhand-side 0
c
      REAL  pp, pad
      integer ii(2),icon
      equivalence (pp,ii)
_IF(littleendian)
      data ii(2)/0/
      ii(1) = icon
_ELSE
      data ii(1)/0/
      ii(2) = icon
_ENDIF
      pad = pp
c
      return
      end
      function ipad(con)
c
c     upack2 of only righthand (cf. pad)
c
      REAL  pp,con
      integer ii(2)
      equivalence (pp,ii)
c
      pp = con
_IF(littleendian)
      ipad = ii(1)
_ELSE
      ipad = ii(2)
_ENDIF
c
      return
      end
_ENDIF(bits8)

_IFN(hitachi)
_IFN(bits8)
      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines 
c      (M4 flags were GEN_OPTIONS)
c
      implicit REAL (a-h,o-z),integer(i-n)
c
_IF(apollo)
      integer*2 ii(4),jj,kk,zero,mask
_ELSEIF(littleendian)
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
_ELSE
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
_ENDIF
      integer*4 itable(15)
_IF(linux)
      integer*2 ishft2
_ENDIF
      parameter (zero=0,mask=15)
      equivalence(a,ii(1))
      data itable/3,2,2,1,1,1,1,8*0/
_IFN1(p)      data mm4,mm8/-4,-8/
c
      a  =  b
_IF1(r)      leadzz = 64
_IFN1(r)      leadz = 64
c
      do 100 i=1,4
_IF(littleendian)
         iii = ii(5-i)
         if (iii.ne.zero) then
_ELSE
         iii = ii(i)
         if (iii.ne.zero) then
_ENDIF
c...  we found a non-zero 2-byte number so split it
_IF(apollo)
            jj = rshft(iii,8)
_ELSEIF(linux)
            jj = ishft2(iii,mm8)
_ELSEIF(littleendian)
            jj = ishft(iii,mm8)
_ELSE
            jj = ishft(iii,mm8)
_ENDIF
            if (jj.ne.zero) then
               koff = 0
            else
               jj = iii
               koff = 8
            end if
c...  and split it again
_IF(apollo)
            kk = and(rshft(jj,4),mask)
_ELSEIF(linux)
            kk = IAND32(ishft2(jj,mm4),mask)
_ELSE
            kk = IAND32(ishft(jj,mm4),mask)
_ENDIF
            if (kk.ne.zero) then
_IF1(r)               leadzz = (i-1)*16 + itable(kk)  + koff
_IFN1(r)               leadz = (i-1)*16 + itable(kk)  + koff
            else
_IF1(r)               leadzz = (i-1)*16 + itable(jj) + koff + 4
_IFN1(r)               leadz = (i-1)*16 + itable(jj) + koff + 4
            end if
            return
         end if
100   continue
c
      return
      end
_ENDIF(bits8)
_ELSE

      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines 
c      (M4 flags were GEN_OPTIONS)
c
      implicit REAL (a-h,o-z),integer(i-n)
c
_IF(apollo)
      integer*2 ii(4),jj,kk,zero,mask
_ELSEIF(littleendian)
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
_ELSE
      integer*2 ii(4),jj,kk,zero,mask,mm4,mm8
_ENDIF
      integer*4 itable(15)
_IF(linux,hitachi)
      integer*2 ishft2
_ENDIF
chitachi
cccc      parameter (zero=0,mask=15)
      equivalence(a,ii(1))
      data itable/3,2,2,1,1,1,1,8*0/
_IFN1(p)      data mm4,mm8/-4,-8/
c
      a  =  b
chitachi
      zero=0
      mask=15
_IF1(r)      leadzz = 64
_IFN1(r)      leadz = 64
c
      do 100 i=1,4
_IF(littleendian)
         iii = ii(5-i)

         if (iii.ne.zero) then
_ELSE
         if (ii(i).ne.zero) then
_ENDIF
c...  we found a non-zero 2-byte number so split it
_IF(apollo)
            jj = rshft(ii(i),8)
_ELSEIF(hitachi)
            jj = ishft2(ii(i),mm8)
_ELSEIF(linux)
            jj = ishft2(iii,mm8)
_ELSEIF(littleendian)
            jj = ishft(iii,mm8)
_ELSE
            jj = ishft(ii(i),mm8)
_ENDIF
            if (jj.ne.zero) then
               koff = 0
            else
_IF(littleendian)
               jj = iii
_ELSE
               jj = ii(i)
_ENDIF
               koff = 8
            end if
c...  and split it again
_IF(apollo)
            kk = and(rshft(jj,4),mask)
_ELSEIF(hitachi)
c
c cant do iand on integer*2s
c
        ii1=mask
        ii2=ishft2(jj,mm4)
        kk = IAND32(ii1,ii2)
c            kk = IAND32(ishft2(jj,mm4),mask)
_ELSEIF(linux)
            kk = IAND32(ishft2(jj,mm4),mask)
_ELSE
            kk = IAND32(ishft(jj,mm4),mask)
_ENDIF
            if (kk.ne.zero) then
_IF1(r)               leadzz = (i-1)*16 + itable(kk)  + koff
_IFN1(r)               leadz = (i-1)*16 + itable(kk)  + koff
            else
_IF1(r)               leadzz = (i-1)*16 + itable(jj) + koff + 4
_IFN1(r)               leadz = (i-1)*16 + itable(jj) + koff + 4
            end if
            return
         end if
100   continue
c
      return
      end

_ENDIF(hitachi)


_IF(debug-bits)
      subroutine pr2b(label,word)
      implicit none
      real*8 word
      real*8 a
      integer integ(2), i
      equivalence (a,integ(1))
      character*64 string
      character*(*) label
      logical btest

      a = word

      do i = 0,31
         string(32-i:32-i)='0'
         if(btest(integ(2),i))string(32-i:32-i)='1'
      enddo
      do i = 0,31
         string(64-i:64-i)='0'
         if(btest(integ(1),i))string(64-i:64-i)='1'
      enddo

      write(6,101)string,label
 101  format(1x,a64,1x,a20)
      end


      subroutine pr1i(label,word)
      implicit none
      integer word, i
      character*32 string
      character*(*) label
      logical btest

      do i = 0,31
         string(32-i:32-i)='0'
         if(btest(word,i))string(32-i:32-i)='1'
      enddo     

      write(6,101)string,label
 101  format(1x,a32,1x,a20)

      end
      subroutine pr1s(label,short)
      implicit none
      integer*2 short,temp(2)
      integer i
      character*16 string
      character*(*) label
      logical btest
      integer word
      temp(1)=short
      temp(2)=0

      do i = 0,15
         string(16-i:16-i)='0'
         if(btest(temp,i))string(16-i:16-i)='1'
      enddo     

      write(6,101)string,label
 101  format(1x,a16,1x,a20)
      end
_ENDIF(debug-bits)
_IFN(bits8)
      integer function popcnt(b)
c
c...   64-bit popcount
c...   cray-imitation for direct for 32 bit machines
c...   M4 flags user were GEN_OPTIONS
c
      implicit REAL (a-h,o-z),integer(i-n)
c
      integer*4 itable(0:15),ii(2)
      integer*4 mask
      parameter (mask=15)
      equivalence(a,ii(1))
      data itable/ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4/
c
      a  =  b
      kk = 0
c

      do 100 i=1,2
         do 10 j=1,8
_IF(apollo)
            jj = and(rshft(ii(i),(j-1)*4),mask)
_ELSEIF(littleendian)
            jj = IAND32(ishft(ii(3-i),-(j-1)*4),mask)
_ELSE
            jj = IAND32(ishft(ii(i),-(j-1)*4),mask)
_ENDIF
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk

c
      return
      end
      function dand(a,b)
c
c...   64-bit and function
c...   a and b are real*8 in outside world
c
      REAL  rr, dand
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IAND32(a(1),b(1))
      r(2) = IAND32(a(2),b(2))
      dand = rr
c
      return
      end
_IFN(vpp300)
      function dor(a,b)
c
c...   64-bit and function
c...   a and b are real*8 in outside world
c
      REAL  rr, dor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IOR32(a(1),b(1))
      r(2) = IOR32(a(2),b(2))
      dor = rr
c
      return
      end
      function dxor(a,b)
c
c...   64-bit xor function
c...   a and b are real*8 in outside world
c
      REAL  rr, dxor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IXOR32(a(1),b(1))
      r(2) = IXOR32(a(2),b(2))
      dxor = rr
c
      return
      end
      function SHIFT(x,nn)
c
c...   real*8 shift (with wraparound)
c...   n <0 shifts aborts
c...   n >64 is OK now (used in direct)
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (nn.lt.0) then
       write(6,1000) nn
1000   format(' shift called wrongly',i20)
       call caserr('invalid argument in shiftc')
      endif
c
      n = IAND32(nn,63)
c
      a = x
      if (n.eq.0.or.n.eq.64) then
         b = a
      else if (n.lt.32) then
_IF(apollo)
         ib(1) = or(lshft(ia(1),n),rshft(ia(2),32-n))
         ib(2) = or(lshft(ia(2),n),rshft(ia(1),32-n))
_ELSEIF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
_ENDIF
      else if (n.eq.32) then
         ib(1) = ia(2)
         ib(2) = ia(1)
      else
_IF(apollo)
         ib(1) = or(lshft(ia(2),n-32),rshft(ia(1),64-n))
         ib(2) = or(lshft(ia(1),n-32),rshft(ia(2),64-n))
_ELSEIF(littleendian)
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
_ELSE
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
_ENDIF
      end if
c
_IFN1(gh)      shift = b
_IF1(gh)      shiftc = b
c
       return
      end
      function shiftl(x,n)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl called wrongly')
       call caserr('invalid argument in shiftl')
      endif
c
      a = x
      if (n.eq.0) then
         b = a
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else if (n.lt.32) then
_IF(apollo)
         ib(1) = or(lshft(ia(1),n),rshft(ia(2),32-n))
         ib(2) = lshft(ia(2),n)
_ELSEIF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = ishft(ia(1),n)
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ishft(ia(2),n)
_ENDIF
      else if (n.eq.32) then
_IF(littleendian)
         ib(2) = ia(1)
         ib(1) = 0
_ELSE
         ib(1) = ia(2)
         ib(2) = 0
_ENDIF
      else
_IF(apollo)
         ib(1) = lshft(ia(2),n-32)
         ib(2) = 0
_ELSEIF(littleendian)
         ib(2) = ishft(ia(1),n-32)
         ib(1) = 0
_ELSE
         ib(1) = ishft(ia(2),n-32)
         ib(2) = 0
_ENDIF
      end if
c
      shiftl = b
c
      return
      end
      function shiftr(x,n)
c
c...   real*8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word !!
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000) n
1000   format(' shiftr called wrongly',i20)
       call caserr('invalid argument in shiftr')
      endif
c
      a = x
c
      if (n.eq.0) then
         b = a
      else if (n.eq.32) then
_IF(littleendian)
         ib(1) = ia(2)
         ib(2) = 0
_ELSE
         ib(2) = ia(1)
         ib(1) = 0
_ENDIF
      else if (n.lt.32) then
_IF(apollo)
         ib(2) = or(rshft(ia(2),n),lshft(ia(1),32-n))
         ib(1) = rshft(ia(1),n)
_ELSEIF(littleendian)
         ib(1) = IOR32(ishft(ia(1),-n),ishft(ia(2),32-n))
         ib(2) = ishft(ia(2),-n)
_ELSE
         ib(2) = IOR32(ishft(ia(2),-n),ishft(ia(1),32-n))
         ib(1) = ishft(ia(1),-n)
_ENDIF
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else
_IF(apollo)
         ib(2) = rshft(ia(1),n-32)
         ib(1) = 0
_ELSEIF(littleendian)
         ib(1) = ishft(ia(2),32-n)
         ib(2) = 0
_ELSE
         ib(2) = ishft(ia(1),32-n)
         ib(1) = 0
_ENDIF
      end if
c
      shiftr = b
c
      return
      end
_ENDIF(vpp300)
_ENDIF(bits8)
_IF(bits8)
_IF(linux)
c****************************************************************
c************ bitmanipulation routines                ***********
c**** subroutine rather than function driven to avoid ***********
c************* pentium FP problems                    ***********
c**** this version hardwired for littleendian and G   ***********
c**** routines/functions affected:
c**** pad       =>       pad8 
c**** dand      =>       dand8 
c**** dor       =>       dor8 
c**** dxor      =>       dxor8 
c**** shift     =>       shift8
c**** shiftl    =>       shiftl8
c**** shiftr    =>       shiftr8
c**** also included here are linux versions for:
c**** ipad
c**** popcnt
c**** leadz
c****************************************************************
      subroutine pad8(icon,ipad)
c
c     pack2 with lefhand-side 0
c
      implicit none
      integer icon,ipad(2)
      ipad(2) = 0
      ipad(1) = icon
c
      return
      end
c
      integer function ipad(ii)
c
c     upack2 of only righthand (cf. pad)
c
      implicit none
      integer ii(2)
c
      ipad = ii(1)
c
      return
      end

      subroutine dand8(a,b,dand)
c
c...   64-bit and subroutine
c...   a and b are real*8 in outside world
c
      implicit none
      integer*4 a(2),b(2),dand(2)
_IF(absoft)
      integer IAND32
_ENDIF
c
      dand(1) = IAND32(a(1),b(1))
      dand(2) = IAND32(a(2),b(2))
c
      return
      end
      subroutine dor8(a,b,dor)
c
c...   64-bit or subroutine
c...   a, b and dor are real*8 in outside world
c
      implicit none
      integer*4 a(2),b(2),dor(2)
_IF(absoft)
      integer IOR32
_ENDIF
c
      dor(1) = IOR32(a(1),b(1))
      dor(2) = IOR32(a(2),b(2))
c
      return
      end
      subroutine dxor8(a,b,dxor)
c
c...   64-bit xor subroutine
c...   a, b and dxorare real*8 in outside world
c
      implicit none
      integer*4 a(2),b(2),dxor(2)
_IF(absoft)
      integer IXOR32
_ENDIF
c
      dxor(1) = IXOR32(a(1),b(1))
      dxor(2) = IXOR32(a(2),b(2))
c
      return
      end
      subroutine shift8(ix,n,ishift)
c
c...   real*8 shift (with wraparound)
c...   n < 0 shifts abort
c...   n >64 is OK now (used in direct)
c...   lshft or rshft by 32 does not empty word !!
c
      implicit none
      integer ix,n,ishift,ia
_IF(absoft)
      integer IOR32
_ENDIF
c
      dimension ia(2),ishift(2),ix(2)
c
      if (n.lt.0) then
       write(6,1000)n
1000   format(' shift8 called wrongly',i20)
       call caserr('invalid argument in shift8')
      endif
c
c      write(6,20) ix(2), ix(1)
c20    format(' shift8: argument = ', z8,1x,z8)
      ia(1) = ix(1)
      ia(2) = ix(2)
      if (n.eq.0.or.n.eq.64) then
         ishift(1) = ia(1)
         ishift(2) = ia(2)
      else if (n.lt.32) then
         ishift(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ishift(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
      else if (n.eq.32) then
         ishift(1) = ia(2)
         ishift(2) = ia(1)
      else
         ishift(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
         ishift(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
      end if
c
      return
      end
      subroutine shiftl8(ix,n,ishiftl)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !!
c
      implicit none
      integer ix,ishiftl,n,ia
_IF(absoft)
      integer IOR32
_ENDIF
c
      dimension ishiftl(2),ix(2),ia(2)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl8 called wrongly')
       call caserr('invalid argument in shiftl8')
      endif
c
      ia(1) = ix(1)
      ia(2) = ix(2)
      if (n.eq.0) then
         ishiftl(1) = ia(1) 
         ishiftl(2) = ia(2)
      else if (n.eq.64) then
         ishiftl(1) = 0
         ishiftl(2) = 0
      else if (n.lt.32) then
         ishiftl(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ishiftl(1) = ishft(ia(1),n)
      else if (n.eq.32) then
         ishiftl(2) = ia(1)
         ishiftl(1) = 0
      else
         ishiftl(2) = ishft(ia(1),n-32)
         ishiftl(1) = 0
      end if
c
      return
      end
      subroutine shiftr8(ix,n,ishiftr)
c
c...   real*8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word 
c...   lshft or rshft by 32 does not empty word 
c
      implicit none
      integer ix,n,ishiftr,ia
_IF(absoft)
      integer IOR32
_ENDIF
c
      dimension ishiftr(2),ix(2),ia(2)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftr8 called wrongly')
       call caserr('invalid argument in shiftr8')
      endif
c
      ia(1) = ix(1)
      ia(2) = ix(2)
      if (n.eq.0) then
         ishiftr(1) = ia(1)
         ishiftr(2) = ia(2)
      else if (n.eq.32) then
         ishiftr(1) = ia(2)
         ishiftr(2) = 0
      else if (n.lt.32) then
         ishiftr(1) = IOR32(ishft(ia(1),-n),ishft(ia(2),32-n))
         ishiftr(2) = ishft(ia(2),-n)
      else if (n.eq.64) then
         ishiftr(1) = 0
         ishiftr(2) = 0
      else
         ishiftr(1) = ishft(ia(2),32-n)
         ishiftr(2) = 0
      end if
c
      return
      end
      integer function popcnt(ii)
c
c...   64-bit popcount
c...   cray-imitation for direct for 32 bit machines
c...   M4 flags user were linux,littleendian,cio,unix,doublebackslash,bits8
c
      implicit none
      integer *4 itable(0:15),ii(2),mask
      integer i,j,jj,kk
_IF(absoft)
      integer IAND32
_ENDIF
c
      parameter (mask=15)
      data itable/ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4/
c
      kk = 0
c
      do 100 i=1,2
         do 10 j=1,8
            jj = IAND32(ishft(ii(3-i),-(j-1)*4),mask)
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk

c
      return
      end
      function leadz(ii)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines
c      (M4 flags were linux,littleendian,cio,unix,doublebackslash,bits8)
c
      implicit none
_IF(absoft)
      integer IAND32
_ENDIF
c
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
      integer*4 itable(15),i,koff,leadz
      integer*2 ishft2
      parameter (zero=0,mask=15)
      data itable/3,2,2,1,1,1,1,8*0/
      data mm4,mm8/-4,-8/
c
      leadz = 64
c
      do 100 i=1,4
         iii = ii(5-i)
         if (iii.ne.zero) then
c...  we found a non-zero 2-byte number so split it
            jj = ishft2(iii,mm8)
c           write(6,*) 'i,iii,jj = ', i,iii,jj
            if (jj.ne.zero) then
               koff = 0
            else
               jj = iii
               koff = 8
            end if
c...  and split it again
            kk = IAND32(ishft2(jj,mm4),mask)
c           write(6,*) 'i,jj,kk = ', i,jj,kk
            if (kk.ne.zero) then
               leadz = (i-1)*16 + itable(kk)  + koff
            else
               leadz = (i-1)*16 + itable(jj) + koff + 4
            end if
c           write(6,*)' leadz = ', leadz
            return
         end if
100   continue
c
      return
      end
_ELSE
c****************************************************************
c************ bitmanipulation routines                ***********
c**** subroutine rather than function driven to avoid ***********
c************* pentium FP problems                    ***********
c**** routines/functions affected:
c**** pad       =>       pad8 
c**** dand      =>       dand8 
c**** dor       =>       dor8 
c**** dxor      =>       dxor8 
c**** shift     =>       shift8
c**** shiftl    =>       shiftl8
c**** shiftr    =>       shiftr8
c****************************************************************
      subroutine pad8(icon,pad)
c
c     pack2 with lefhand-side 0
c
      REAL  pp, pad
      integer ii(2),icon
      equivalence (pp,ii)
_IF(littleendian)
      data ii(2)/0/
      ii(1) = icon
_ELSE
      data ii(1)/0/
      ii(2) = icon
_ENDIF
      pad = pp
c
      return
      end

      subroutine dand8(a,b,dand)
c
c...   64-bit and subroutine
c...   a and b are real*8 in outside world
c
      REAL  rr, dand
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IAND32(a(1),b(1))
      r(2) = IAND32(a(2),b(2))
      dand = rr
c
      return
      end
      subroutine dor8(a,b,dor)
c
c...   64-bit or subroutine
c...   a and b are real*8 in outside world
c
      REAL  rr, dor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IOR32(a(1),b(1))
      r(2) = IOR32(a(2),b(2))
      dor = rr
c
      return
      end
      subroutine dxor8(a,b,dxor)
c
c...   64-bit xor subroutine
c...   a and b are real*8 in outside world
c
      REAL  rr, dxor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IXOR32(a(1),b(1))
      r(2) = IXOR32(a(2),b(2))
      dxor = rr
c
      return
      end
      subroutine shift8(x,n,shift)
c
c...   real*8 shift (with wraparound)
c...   n < 0 shift aborts
c...   n >64 is OK now (used in direct)
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.lt.0) then
       write(6,1000)n
1000   format(' shift8 called wrongly',i20)
       call caserr('invalid argument in shift8')
      endif
c

      a = x
      if (n.eq.0.or.n.eq.64) then
         b = a
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
_ENDIF
      else if (n.eq.32) then
         ib(1) = ia(2)
         ib(2) = ia(1)
      else
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
_ELSE
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
_ENDIF
      end if
c
      shift = b
      return
      end
      subroutine shiftl8(x,n,shiftl)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl8 called wrongly')
       call caserr('invalid argument in shiftl8')
      endif
c
      a = x
      if (n.eq.0) then
         b = a
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = ishft(ia(1),n)
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ishft(ia(2),n)
_ENDIF
      else if (n.eq.32) then
_IF(littleendian)
         ib(2) = ia(1)
         ib(1) = 0
_ELSE
         ib(1) = ia(2)
         ib(2) = 0
_ENDIF
      else
_IF(littleendian)
         ib(2) = ishft(ia(1),n-32)
         ib(1) = 0
_ELSE
         ib(1) = ishft(ia(2),n-32)
         ib(2) = 0
_ENDIF
      end if
c
      shiftl = b
c
      return
      end
      subroutine shiftr8(x,n,shiftr)
c
c...   real*8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word 
c...   lshft or rshft by 32 does not empty word 
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftr8 called wrongly')
       call caserr('invalid argument in shiftr8')
      endif
c
      a = x
c
      if (n.eq.0) then
         b = a
      else if (n.eq.32) then
_IF(littleendian)
         ib(1) = ia(2)
         ib(2) = 0
_ELSE
         ib(2) = ia(1)
         ib(1) = 0
_ENDIF
      else if (n.lt.32) then
_IF(littleendian)
         ib(1) = IOR32(ishft(ia(1),-n),ishft(ia(2),32-n))
         ib(2) = ishft(ia(2),-n)
_ELSE
         ib(2) = IOR32(ishft(ia(2),-n),ishft(ia(1),32-n))
         ib(1) = ishft(ia(1),-n)
_ENDIF
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else
_IF(littleendian)
         ib(1) = ishft(ia(2),32-n)
         ib(2) = 0
_ELSE
         ib(2) = ishft(ia(1),32-n)
         ib(1) = 0
_ENDIF
      end if
c
      shiftr = b
c
      return
      end
_ENDIF
_ENDIF
_ENDIF
_IF(ksr)
_IFN(t3d)
      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c
      implicit REAL (a-h,o-z),integer(i-n)
c
      integer*2 ii(4),jj,kk,zero,mask,mm4,mm8
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
         if (ii(i).ne.zero) then
c...  we found a non-zero 2-byte number so split it
            jj = ishft(ii(i),mm8)
            if (jj.ne.zero) then
               koff = 0
            else
               jj = ii(i)
               koff = 8
            end if
c...  and split it again
            kk = iand(ishft(jj,mm4),mask)
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
c
      implicit REAL (a-h,o-z),integer(i-n)
c
      integer*4 itable(0:15),ii(2)
      parameter (mask=15)
      equivalence(a,ii(1))
      data itable/ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4/
c
      a  =  b
      kk = 0
c
      do 100 i=1,2
         do 10 j=1,8
            jj = iand(ishft(ii(i),-(j-1)*4),mask)
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk
c
      return
      end
_ENDIF
_ENDIF
_IF(i8)
c****************************************************************
c************ bitmanipulation routines (i8 machines)  ***********
c************ temp. fix to check out base modules     ***********
c************ note that most of these functions are defunct *****
c****************************************************************
_IF()
      function pad(icon)
c
c     pack2 with lefhand-side 0
c
      REAL  pp, pad
      integer *4 ii(2),icon
      equivalence (pp,ii)
_IF(littleendian)
      data ii(2)/0/
      ii(1) = icon
_ELSE
      data ii(1)/0/
      ii(2) = icon
_ENDIF
      pad = pp
c
      return
      end
      function ipad(con)
c
c     upack2 of only righthand (cf. pad)
c
      REAL  pp,con
      integer *4 ii(2)
      equivalence (pp,ii)
c
      pp = con
_IF(littleendian)
      ipad = ii(1)
_ELSE
      ipad = ii(2)
_ENDIF
c
      return
      end
_ENDIF
      function leadz(b)
c
c...   count the number of leading zero's in a 64-bit word
c      imitation for 32 bit machines 
c      (M4 flags were GEN_OPTIONS)
c
      implicit REAL (a-h,o-z),integer(i-n)
c
_IF(littleendian)
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
_ELSE
      integer*2 ii(4),jj,kk,zero,mask,iii,mm4,mm8
_ENDIF
      integer*4 itable(15)
      parameter (zero=0,mask=15)
      equivalence(a,ii(1))
      data itable/3,2,2,1,1,1,1,8*0/
_IFN1(p)      data mm4,mm8/-4,-8/
c
      a  =  b
_IF1(r)      leadzz = 64
_IFN1(r)      leadz = 64
c
      do 100 i=1,4
_IF(littleendian)
         iii = ii(5-i)
         if (iii.ne.zero) then
_ELSE
         iii = ii(i)
         if (iii.ne.zero) then
_ENDIF
c...  we found a non-zero 2-byte number so split it
_IF(littleendian)
            jj = ishft(iii,mm8)
_ELSE
            jj = ishft(iii,mm8)
_ENDIF
            if (jj.ne.zero) then
               koff = 0
            else
               jj = iii
               koff = 8
            end if
c...  and split it again
            kk = IAND32(ishft(jj,mm4),mask)
            if (kk.ne.zero) then
_IF1(r)               leadzz = (i-1)*16 + itable(kk)  + koff
_IFN1(r)               leadz = (i-1)*16 + itable(kk)  + koff
            else
_IF1(r)               leadzz = (i-1)*16 + itable(jj) + koff + 4
_IFN1(r)               leadz = (i-1)*16 + itable(jj) + koff + 4
            end if
            return
         end if
100   continue
c
      return
      end

_IF(debug-bits)
      subroutine pr2b(label,word)
      implicit none
      real*8 word
      real*8 a
      integer integ(2), i
      equivalence (a,integ(1))
      character*64 string
      character*(*) label
      logical btest

      a = word

      do i = 0,31
         string(32-i:32-i)='0'
         if(btest(integ(2),i))string(32-i:32-i)='1'
      enddo
      do i = 0,31
         string(64-i:64-i)='0'
         if(btest(integ(1),i))string(64-i:64-i)='1'
      enddo

      write(6,101)string,label
 101  format(1x,a64,1x,a20)
      end


      subroutine pr1i(label,word)
      implicit none
      integer word, i
      character*32 string
      character*(*) label
      logical btest

      do i = 0,31
         string(32-i:32-i)='0'
         if(btest(word,i))string(32-i:32-i)='1'
      enddo     

      write(6,101)string,label
 101  format(1x,a32,1x,a20)

      end
      subroutine pr1s(label,short)
      implicit none
      integer*2 short,temp(2)
      integer i
      character*16 string
      character*(*) label
      logical btest
      integer word
      temp(1)=short
      temp(2)=0

      do i = 0,15
         string(16-i:16-i)='0'
         if(btest(temp,i))string(16-i:16-i)='1'
      enddo     

      write(6,101)string,label
 101  format(1x,a16,1x,a20)
      end
_ENDIF(debug-bits)

      integer function popcnt(b)
c
c...   64-bit popcount
c...   cray-imitation for direct for 32 bit machines
c...   M4 flags user were GEN_OPTIONS
c
      implicit REAL (a-h,o-z),integer(i-n)
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
_IF(littleendian)
            jj = IAND32(ishft(ii(m3-i),-(j-m1)*m4),mask)
_ELSE
            jj = IAND32(ishft(ii(i),-(j-m1)*m4),mask)
_ENDIF
            kk = kk + itable(jj)
10       continue
100   continue
c
      popcnt = kk

c
      return
      end
_IF()
      function dand(a,b)
c
c...   64-bit and function
c...   a and b are real*8 in outside world
c
      REAL  rr, dand
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IAND32(a(1),b(1))
      r(2) = IAND32(a(2),b(2))
      dand = rr
c
      return
      end
      function dor(a,b)
c
c...   64-bit and function
c...   a and b are real*8 in outside world
c
      REAL  rr, dor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IOR32(a(1),b(1))
      r(2) = IOR32(a(2),b(2))
      dor = rr
c
      return
      end
      function dxor(a,b)
c
c...   64-bit xor function
c...   a and b are real*8 in outside world
c
      REAL  rr, dxor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
_IF1(Gp)      data r/2*0/
c
      r(1) = IXOR32(a(1),b(1))
      r(2) = IXOR32(a(2),b(2))
      dxor = rr
c
      return
      end
      function SHIFT(x,nn)
c
c...   real*8 shift (with wraparound)
c...   n <0 shifts aborts
c...   n >64 is OK now (used in direct)
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      integer *4 ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (nn.lt.0) then
       write(6,1000) nn
1000   format(' shift called wrongly',i20)
       call caserr('invalid argument in shiftc')
      endif
c
      n = IAND32(nn,63)
c
      a = x
      if (n.eq.0.or.n.eq.64) then
         b = a
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
_ENDIF
      else if (n.eq.32) then
         ib(1) = ia(2)
         ib(2) = ia(1)
      else
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
_ELSE
         ib(1) = IOR32(ishft(ia(2),n-32),ishft(ia(1),n-64))
         ib(2) = IOR32(ishft(ia(1),n-32),ishft(ia(2),n-64))
_ENDIF
      end if
c
_IFN1(gh)      shift = b
_IF1(gh)      shiftc = b
c
       return
      end
      function shiftl(x,n)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      integer *4 ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl called wrongly')
       call caserr('invalid argument in shiftl')
      endif
c
      a = x
      if (n.eq.0) then
         b = a
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else if (n.lt.32) then
_IF(littleendian)
         ib(2) = IOR32(ishft(ia(2),n),ishft(ia(1),n-32))
         ib(1) = ishft(ia(1),n)
_ELSE
         ib(1) = IOR32(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ishft(ia(2),n)
_ENDIF
      else if (n.eq.32) then
_IF(littleendian)
         ib(2) = ia(1)
         ib(1) = 0
_ELSE
         ib(1) = ia(2)
         ib(2) = 0
_ENDIF
      else
_IF(littleendian)
         ib(2) = ishft(ia(1),n-32)
         ib(1) = 0
_ELSE
         ib(1) = ishft(ia(2),n-32)
         ib(2) = 0
_ENDIF
      end if
c
      shiftl = b
c
      return
      end
      function shiftr(x,n)
c
c...   real*8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word !!
c...   lshft or rshft by 32 does not empty word !!
c
      implicit REAL  (a-h,o-z), integer(i-n)
c
      integer *4 ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000) n
1000   format(' shiftr called wrongly',i20)
       call caserr('invalid argument in shiftr')
      endif
c
      a = x
c
      if (n.eq.0) then
         b = a
      else if (n.eq.32) then
_IF(littleendian)
         ib(1) = ia(2)
         ib(2) = 0
_ELSE
         ib(2) = ia(1)
         ib(1) = 0
_ENDIF
      else if (n.lt.32) then
_IF(littleendian)
         ib(1) = IOR32(ishft(ia(1),-n),ishft(ia(2),32-n))
         ib(2) = ishft(ia(2),-n)
_ELSE
         ib(2) = IOR32(ishft(ia(2),-n),ishft(ia(1),32-n))
         ib(1) = ishft(ia(1),-n)
_ENDIF
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else
_IF(littleendian)
         ib(1) = ishft(ia(2),32-n)
         ib(2) = 0
_ELSE
         ib(2) = ishft(ia(1),32-n)
         ib(1) = 0
_ENDIF
      end if
c
      shiftr = b
c
      return
      end
_ENDIF
_ENDIF
_ENDEXTRACT
_IF(cio,fcio)
c
c** cio/fcio ioerr
c
      subroutine ioerr(oper,iunit,file)
      character oper*(*), file*(*)
      character *80 errstr,temp
      logical ocaser
INCLUDE(common/sizes)
INCLUDE(common/discc)
INCLUDE(common/syscod)
INCLUDE(common/iofile)
INCLUDE(common/errcodes)
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
_ENDIF
_IF(cio)
c
c** cio iolog
c
      subroutine iolog(icode,mes)
      character mes*100
INCLUDE(common/syscod)
c
c  to save any useful error info from system
c  (called from c i/o routines)
c
      icode1=icode
      mess = mes
      return
      end
_ENDIF
_IFN(ipsc,charmm)
c
c** non-ipsc dclock
c
      function dclock()
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c use elapsed times here on non-ipsc 
      call walltime(dum)
      dclock = dum
      return
      end
_ENDIF
_IF(ipsc)
c
c  ===========    iPSC specific routines  ====================
c
c** ipsc synget
c
      subroutine synget(buf,nword,nam,irec)
c
c...  synchronous host input
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension buf(nword)
c
      read(nam,rec=irec,err=999) buf
c
c
999   call caserr(' error on synchronous input')
      return
      end
c
c** ipsc synput
c
      subroutine synput(buf,nword,nam,irec)
c
c...  synchronous host output
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension buf(nword)
c
      write(nam,rec=irec,err=999) buf
      return
c
999   call caserr(' error on synchronous output')
      return
      end
c
c** ipsc iowt
c
      subroutine iowt(itype)
c
c***   ***NODE-IPSC***
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
INCLUDE(common/nodinf)
INCLUDE(common/nodeio)
c
   20 numb = 0
   10 if(iodone(itype).ne.0) then
         maxio(isel)=max(maxio(isel),numb)
      else
         call flick()
         numb=numb+1
         if(numb.gt.10000) then
            print *,' node ',minode,' **** warning on iowait ',itype
            go to 20
         end if
         go to 10
      endif
c
      return
      end
c
c** ipsc msgw
c
      subroutine msgw(itype)
c
c***   ***NODE-IPSC***
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
INCLUDE(common/nodinf)
c
   20 numb = 0
   10 if(msgdone(itype).ne.0) then
         maxnod(isel)=max(maxnod(isel),numb)
         return
      else
         call flick()
         numb=numb+1
         if(numb.gt.10000) then
            print *,' node ',minode,' **** error on wait ',itype
            go to 20
         end if
         go to 10
      endif
c
      end
c
c** ipsc cprob
c
      subroutine cprob(itype)
c
c***   ***NODE-IPSC***
c
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c
      character *20 zinit
INCLUDE(common/disc)
INCLUDE(common/nodinf)
      common/vdisk/zinit
      common/scfopt/maxit(4),accdi(2),odiis(4),dmpc(7),iter
c
   20 numb = 0
   10 if(iprobe(itype).eq.1) then
         maxnod(isel)=max(maxnod(isel),numb)
         return
      else
         call flick()
         numb=numb+1
         if(numb.gt.1000000) then
            write(6,40)minode,itype
 40   format(' node ',i2,' ****  probe ',i8)
c           write(6,40)minode,itype,iter,zinit
c40   format(' node ',i2,' ****  probe ',i8,
c    * ' iter =',i3,
c    *' :initiated from ',a20)
            go to 20
         endif
         go to 10
      endif
c
      end
_ENDIF
_IF(alpha)
_IF(vms)
c
c** vms cputime
c
      SUBROUTINE CPUTIME(TCOMP)
*
*
*       CPU time.
*       ---------
*
      real*8  TCOMP(3), dum

      COMMON/ICPU0/ ICPU
      real*4  TARRAY(2)
*
*
*       Initialize timer (first call) or obtain elapsed time.
      IF (ICPU.EQ.0) THEN
         CALL LIB$INIT_TIMER
         TCOMP = 0.0
C          TCOMP = ETIME(TARRAY)
          ICPU = 1
      ELSE
         CALL LIB$STAT_TIMER(2,ITIME)
         dum = itime
cps - need to resolve system component
         TCOMP(1) = dum/100.0d0
         TCOMP(2) = 0.0d0
         TCOMP(3) = tcomp(1)
c
*         WRITE(0,*)' TCOMP IN SEC_ETIMEIS',TCOMP
*          TCOMP = ETIME(TARRAY)/60.0
*       Elapsed CPU time in minutes on SUN or MIPS.
      END IF
*
      ETIME=TCOMP
      RETURN
*
      END
      integer function get_for_vals
c
c        This function obtains the foreign command line from the
c     VMS operating system and parses it. It returns the number of
c     arguments (zero if there are none).
c     If there is a system error of some kind, the function prints
c     a message to the terminal and returns -1.
c        The number of arguments and indices into the command string
c     are placed in appropriate common block variables. If there is
c     some error, the number of arguments is considered to be zero.
c
c     Restrictions:
c        Number of characters on line (excluding command and any blanks
c     that follow command) = 250
c        Number of tokens = 50
c
c     John Faricelli   Device/Process Simulation
c     Created: Nov. 2, 1984
c
      include '($ssdef)'
c
      common/foreign/ number_args, token_limits, command_line

      integer number_args  ! number of command line arguments
      integer token_limits(50,2) ! limits of tokens in command_line
      character*250 command_line ! storage for command line
c
      logical have_command_line  ! if TRUE, have gotten command line
                                 ! from system.
      integer j                  ! token counter
      integer flag               ! flag bits (0)
      integer status             ! return status
      integer*2 strlen           ! length of command line string
      integer lib$get_foreign
c
      data have_command_line/.FALSE./
c
      save have_command_line
c
      if( .not. have_command_line ) then
         flag = 0
c
c        get command line from VMS
         status = lib$get_foreign(command_line, , strlen, flag)
         if( status .ne. ss$_normal ) then
            write(6,1000) status
1000        format(' ','internal error in getarg'/
     $             ' ','lib$getforeign returned ',i5)
            get_for_vals = -1
            number_args = 0
            have_command_line = .TRUE.
            return
         endif
c
         have_command_line = .TRUE.
c        test for no command line
         if( strlen .eq. 0 ) then
            get_for_vals = 0
            number_args = 0
            return
         endif
c
c        else parse input
         j = 1
         in_token = .FALSE.
         do i = 1, strlen
c
c           look for blank or comma
            if (command_line(i:i) .ne. ' ' .and.
     $          command_line(i:i) .ne. ',' ) then
               if( .not. in_token ) then
c
c              new token
                  if( j .gt. 50 ) then
c
c                    more tokens than we can handle
                     write(6,1001)
1001                 format(' ','internal error in getarg'/
     $                      ' ','more than 50 tokens on command line')
                     get_for_vals = -1
                     number_args = 0
                     return
                  endif
c
c                 record beginning of token
                  token_limits(j,1) = i
                  in_token = .TRUE.
               endif
            else if( in_token ) then
c
c              hit end of token, record end, advance token counter
               token_limits(j,2) = i-1
               j = j + 1
               in_token = .FALSE.
            endif
         end do
c
c        VMS strips any trailing blanks so we must set end of last token
         token_limits(j,2) = strlen
         get_for_vals = j
         number_args = j
c
c     else we've gotten the command line already; do nothing
      endif
c
      return
      end
c
c** vms iargc
c
      integer function iargc
c
c     This function returns the number of command line arguments
c     The function returns zero if there are no arguments
c     (unlike the UNIX f77 call, where argument zero is the
c     command name). The function returns -1 if there is some
c     system error. An informational message is written to the terminal.
c
c     Restrictions:
c        Number of characters on line (excluding command and any blanks
c     that follow command) = 250
c        Number of tokens = 50
c
c     John Faricelli   Device/Process Simulation
c     Created: Nov. 2, 1984
c
      common/foreign/ number_args, token_limits, command_line
c
      integer number_args  ! number of command line arguments
      integer token_limits(50,2) ! limits of tokens in command_line
      character*250 command_line ! storage for command line
c
      integer get_for_vals, status
c
c     get arguments (get_for_vals only parses the command line
c     once). get_for_val does error messages.
      status = get_for_vals()
      iargc = number_args
      return
      end
c
c** vms getarg
c
      subroutine getarg(k, arg)
c
c     This subroutine gets the kth string on the command line.
c     The function tests for invalid token numbers
c     (<= 0 or > number_args ) and an informational message is
c     written to the terminal.
c
c     Restrictions:
c        Number of characters on line (excluding command and any blanks
c     that follow command) = 250
c        Number of tokens = 50
c
c     John Faricelli   Device/Process Simulation
c     Created: Nov. 2, 1984
c
      common/foreign/ number_args, token_limits, command_line
c
      integer number_args  ! number of command line arguments
      integer token_limits(50,2) ! limits of tokens in command_line
      character*250 command_line ! storage for command line
c
      character*(*) arg   ! argument to get
      integer k           ! which one
c
      integer status
      integer ilen
c
c     Get tokens (just in case). Don't worry. We only parse the
c     command line once.
      status = get_for_vals()
c
c     valid argument number ?
      if (k.eq.0) return
      if( k .lt. 0 .or. k .gt. number_args ) then
         write(6,1000) k, number_args
1000     format(' ','error in use of subroutine getarg'/
     $      ' ','program requested argument ',i4/
     $      ' ','arguments on command line were in range 1 to ',i4)
         return
      endif
c
c     sockit toem. But first, check to see if arg can store entire string.
c
      ilen = LEN(arg)
      if( ilen .lt. (token_limits(k,2) - token_limits(k,1) + 1) ) then
         write(6,1100) k, (token_limits(k,2) - token_limits(k,1) + 1),
     1                 ilen
1100     format(' ','Error in call to getarg:'/
     1          ' ','Token ',i5,' on command line was of length ',i5/
     2          ' ','Input parameter string was of length ',i5/
     3          ' ','Token will be truncated')
      endif
      arg = command_line(token_limits(k,1):token_limits(k,2))
      return
      end
c
c** vms getlog
c
          character*31 function getlog ()
          character*8 name
          call fgetenv('USER',name)
          getlog = name
          return
          end
c
c** vms flush
c
          SUBROUTINE FLUSH (LU)
c
          include '($SYSSRVNAM)'
          include '($SSDEF)'
c
c
          INTEGER*4 FOR$RAB
          INTEGER STATUS
c
          STATUS = SYS$FLUSH(%VAL(FOR$RAB(LU)))
          RETURN
          END
_ENDIF
_ENDIF
_IF(hp700)
c 
c** hp700 getarg  (emulator for hp-ux)
c** note this is NOT required in HP-UX 11.00 onwards
      subroutine getarg(i,s)
      integer i
      character s*(*)
      l=len(s)
      call igetarg(i,s,l)
      end
c
c** hp700 flush dummy 
c
      subroutine flush
      return
      end
c
c** hp700 exit replacement (for extname compilation)
c
      subroutine exit
      stop
      end
_ENDIF      
_IF(convex)
c
c** convex usage
c
      subroutine usage
      implicit REAL  (a-h,o-z)
      integer*4  iruse(20)
      iw = 0
      call cgetru ( %val(iw) , iruse )
      usert =dfloat(iruse(1)) + dfloat(iruse(2))/1.d6
      systt =dfloat(iruse(3)) + dfloat(iruse(4))/1.d6
      write(6,101) usert,systt,(iruse(j),j=5,18)
101   format(////,1x,129('*'),/
     *1x'* resource usage'48x'|'63x'*'/
     *1x'*'127('-')'*'/,
     *1x'* user-time                         : 'f10.2' seconds'8x'|'
     *,' system-time                       : 'f10.2' seconds'8x'*'/
     *1x'* maximum resident set size utilized:'5xi6' kb             |'
     *,' integral shared memory size       :'5xi6' mb*clock-ticks *'/
     *1x'* integral unshared data size       :'5xi6' mb*clock-ticks |'
     *,' integral unshared stack size      :'5xi6' mb*clock-ticks *'/
     *1x'* page reclaims                     :'5xi6'                |'
     *,' page faults                       :'5xi6'                *'/
     *1x'* swaps                             :'5xi6'                |'
     *,' block input operations            :'5xi6'                *'/
     *1x'* block output operations           :'5xi6'                |'
     *,' number of messages sent           :'5xi6'                *'/
     *1x'* number of messages received       :'5xi6'                |'
     *,' number of signals received        :'5xi6'                *'/
     *1x'* voluntary context switches        :'5xi6'                |'
     *,' involuntary context switches      :'5xi6'                *'/
     *1x,129('*'))
      return
      end
_ENDIF
_IF1()
_IF1() =============== start of CRAY section =================
_IF1()
_IF(cray,unicos)
c ******************************************************
c ******************************************************
c             =   machscf  =
c ******************************************************
c ******************************************************
c
c** cray setio
c
      subroutine setio(lpb,lpbsz)
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(common/sizes)
INCLUDE(common/discc)
      common/maskc/mask(64)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
     +  ,keep(maxlfn),keepf(3,maxfrt)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),
     *junibf(maxbuf),jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c------------ stream data------------------------------------------
c   ipos(maxlfn)  =  block position
c   yed(maxlfn)  =  dsname
c   i40(maxlfn)   =  dsp base point
c   length(maxlfn)=  current length of file in atmol blocks
c   iostat(maxlfn)=  input/1  or  output/2 flag for physical io
c   iobuff(maxlfn)=  store buffer number or maxbuf+1 if no physical io
c   itab(maxlfn*40) =  list of 40 word dsps
c--------------buffer status data---------------------------------------
c   istabf(maxbuf+1)=  put/-1  nothing/0  input/1  output/2
c   ipribf(maxbuf) = buffer priority 0=maximum npbm1=minimum
c   junibf(maxbuf) = atmol stream number
c   jpbbf(maxbuf)  = physical block number
c   junjpb(maxbuf) = junibf<32.or.jpbbf
c   ibfbas(maxbuf) = buffer base point = 0, npbsz*512, 2*npbsz*512,....
c   igetm(maxbuf)  =buffer status flags 1/block present 0/block absent
c   iputm(maxbuf)  =put flags       1/block destined for output
c----------------scalar variables---------------------------------------
c   irep   =reply word for physical io
c   incr   =base point for next active dsp
c   npb    =number of store buffers
c   npbsz  =number of atmol blocks in one buffer
c   npbm1  =npb-1
c   mbuff  =current buffer number
c   kunkpb =atmol stream/physical block number of current buffer
c   lun    =atmol stream in find/get
c   l40    =dsp address in find/get
c   latm   =relative atmol block number in find/get
c   lbl    =absolute atmol block number in find/get
c   lbfbas =buffer base point in find/get
c   lbuff  =buffer buffer in find/get
c-----------------------------------------------------------------------
      irep=0
      kunkpb=0
      n=0
      incr=1
      npb=lpb
      npbsz=lpbsz
      npbm1=npb-1
      mask(1)=shiftl(1,63)
      do 4 i=2,64
 4    mask(i)=or(shiftl(1,64-i),mask(i-1))
      nn=npbm1
      m=npbsz*512
      do 1 loop=1,npb
      ipribf(loop)=nn
      ibfbas(loop)=n
      nn=nn-1
1     n=n+m
      do 2 loop=1,npb
      junjpb(loop)=0
2     istabf(loop)=0
      numf = 50
      do 3 loop=1,maxlfn
      iobuff(loop)=maxbuf+1
      nam(loop)=loop+numf
3     ipos(loop)=0
      return
      end
c
c** cray flushw
c
      subroutine flushw(tab,buff,nblock,iblock,iun)
      dimension tab(*),buff(*)
INCLUDE(common/sizes)
INCLUDE(common/discc)
      common/disc/isppp(3),irep,ispppp,ipos(maxlfn),
     * nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     * iobuff(maxlfn),itab(maxlfn*40),
     * istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     * jpbbf(maxbuf),junjpb(maxbuf),ibfbas(maxbuf),
     * igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
      len=length(iun)
      ibl=len+1
      nfill=iblock-ibl
5     nfill=nfill-nblock
      if(nfill)4,3,3
3     call flushv(tab,buff,nblock,ibl,iun)
      ibl=ibl+nblock
      call chekwr(tab)
      if(irep)11,5,11
11    call iofail(yed(iun),2)
4     nfill=nfill+nblock
      if(nfill)77,77,6
6     call flushv(tab,buff,nfill,ibl,iun)
      call chekwr(tab)
      if(irep)11,77,11
77    call flushv(tab,buff,nblock,iblock,iun)
      ibl=iblock+nblock
      if(ibl.gt.len)length(iun)=ibl-1
      return
      end
c
c** cray initj
c
      subroutine initj(lword)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(unicos)
      character*256 com
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/timez)
INCLUDE(common/jinfo)
INCLUDE(common/iofile)
INCLUDE(common/timeperiods)
INCLUDE(common/segm)
      common/blksiz/nsz,nsz51,nszz(51)
c
c...   set up signal-catching (unix signals)
c
      call catch
c
c ----- initialise job
c
      call start_time_period(TP_ENTIRE)
      call inita(lword,nbuff,ibuff)
      call tidajt(zdate,ztime,zaccno,zanam,isecs)
cpsh
cc      call getcor(lword)
      ti=0.0d0
      tx=0.0d0
      call flags
_IF(unicos)
      call getarg(0,com)
      write(iwr,'('' GAMESS-UK Executable: '',a70)') com
_ENDIF
      write(iwr,7)zanam,zdate,ztime,zaccno,isecs,lword,
     #            lword*8/(1000*1000)
 7    format(
     *40x,'job name   ',a8/
     *40x,'date       ',a8/
     *40x,'time       ',a8/
     *40x,'acct       ',a8//
     *40x,'job time specified ',i5,' seconds'//
     *40x,'main store requested',i10,' words',' =',i6,' mbytes'/)
c
c    data set initialisation
c
      call setio(nbuff,ibuff)
      call initb
c
      return
      end
c
c** cray flags
c
      subroutine flags
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      write(iwr,6)
 6    format(//
     *35x,'*******************************************************',
     *'*********'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*                ===  G A M E S S - U K   ===          ',
     *'        *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'* Generalised Atomic and Molecular Electronic Structure',
     *' System *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*               ===  UNICOS  version 7.0   ===         ',
     *'        *'/
     *35x,'*                                                      ',
     *'        *'/
     *35x,'*******************************************************',
     *'*********'/)
      return
      end
c
c** cray inita
c
      subroutine inita(lword,nbuff,ibuff)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer getenv
      character *132 filatt
INCLUDE(common/sizes)
_IF(parallel)
      character*11 fname
INCLUDE(common/nodinf)
      logical oswed3
      common/nodeio/oswed3(maxlfn)
      common/scftim/tdiag(10)
_ENDIF
INCLUDE(common/iofile)
INCLUDE(common/prints)
INCLUDE(common/utilc)
      common/discne/mlen(maxlfn),length(maxlfn)
INCLUDE(common/work)
INCLUDE(common/machin)
INCLUDE(common/timez)
INCLUDE(common/segm)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     +  nam(maxlfn),keep(maxlfn),
     +  keepf(maxfrt),iposf(maxfrt),oform(maxfrt)
INCLUDE(common/discc)
INCLUDE(common/discc1)
      common/bufnew/ibufa(512,maxblo,maxbuf)
      common/blksiz/nsz,nsz512(8),nsznav,nsstat,
     *  nslen,nsmax ,nsort(40)
      common/blksi2/ns2,ns2512(8),npznav,npstat,
     *  n2len,n2max ,nport(40)
      common/blksi3/ns3,ns3512(8),npdnav,npdtat,
     *  n3len,n3max ,ntort(40)
      common/restri/nfile(63),lda(508),isect(508),ldx(508),
     *  iacsct(508)
INCLUDE(common/atmblk)
INCLUDE(common/maxlen)
      dimension  ytext(15),zinout(2),yynam(maxlfn)
c
      data yynam/'ed0' ,'ed1' ,'ed2' ,'ed3' ,'ed4' ,'ed5' ,'ed6' ,
     *          'ed7' ,'ed8' ,'ed9' ,'ed10','ed11','ed12','ed13',
     *          'ed14','ed15','ed16','ed17','ed18','ed19',
     *          'mt0' ,'mt1' ,'mt2' ,'mt3' ,'mt4' ,'mt5' ,'mt6' ,
     *          'mt7' ,'mt8' ,'mt9' ,'mt10','mt11','mt12','mt13',
     *          'mt14','mt15','mt16','mt17','mt18','mt19'
     *         /
      data ytext/'chan','swit','file','seti','memo','stor','core',
     * 'lpag','time','sort','psor','driv','****','logf',
     * 'tabl'/
      data zblank/' '/
      data yout,yprin/'outp','prin'   /
      data  zinout/  'dataout','punch'/
c
      ichek=0
      isel=0
      iselr=0
      iselw=0
      irep=0
      nslen=0
      n2len=0
      n3len=0
      nsmax=0
      n2max=0
      n3max=0
      nsz=10
      ns2=1
      ns3=1
      nav = lenwrd()
      call setsto(40,1,nsort)
      call setsto(40,1,nport)
      call setsto(40,1,ntort)
      nsort(1)=91
      nport(1)=92
      ntort(1)=93
      ooprnt = .false.
      do 3200 loop=1,maxlfn
      yed(loop)=yynam(loop)
_IF(parallel)
      oswed3(loop) = .true.
_ENDIF
      omem(loop) = .false.
3200  keep(loop) = 0
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
      do 3201 loop=1,maxfrt
      oform(loop) = .false.
      iposf(loop)=0
       do 3202 i=1,maxfrt
          if (i.lt.10) then
           write(zzzz1(i),'(3hft0,i1)') i
           write(zzzz2(i),'(5hftn00,i1)') i
           write(zzzz3(i),'(5hfort.,i1)') i
          else
           write(zzzz1(i),'(2hft,i2)') i
           write(zzzz2(i),'(4hftn0,i2)') i
           write(zzzz3(i),'(5hfort.,i2)') i
          end if
3202   continue
      zftn(loop) = zzzz3(loop)
      zftfil(loop)= zftn(loop)
3201  keepf(loop) = 0
      mach12 = nav
c
c === hardwire /machin/ values
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
      oprintm = .false.
_IF(debug)
      call gmem_set_nanify(.true.)
_ELSE
      call gmem_set_nanify(.false.)
_ENDIF
      call gmem_set_debug(.false.)
      call gmem_set_quiet(.true.)
      otime=.false.
      lword=icoresize()
      nbuff=maxbuf
      ibuff=maxblo
      isecs=0
      ift=0
_IF(parallel)
      do 3010 loop=1,10
3010  tdiag(loop)=0.0d0
_ENDIF
      do 110 i=1,maxlfn
      zzz=yed(i)(1:4)
 110  zedfil(i)=zzz
c
c..     global setenv's (do not delete any file of a class)
c
      ii = getenv('ed*',filatt)
      if (ii.ne.0) call setsto(maxlfn,1,keep)
      ii = getenv('fort*',filatt)
      if (ii.ne.0) call setsto(maxfrt,1,keepf)
c
c --- change the default to 'keep' if environment variable set
c
      do 115 i=1,maxlfn
       filatt = ' '
       ii  =  getenv(yed(i),filatt)
       if(ii.ne.0) keep(i)=1
 115  continue
c
c      and repeat for fortran files
c
       do 101 i=1,50
       filatt = ' '
       ii = getenv(zzzz1(i),filatt) + 
     +      getenv(zzzz2(i),filatt) + 
     +      getenv(zzzz3(i),filatt)
       if (ii.ne.0) keepf(i) = 1
101    continue
c
c --- take a copy of input file
c
      call print0
      oecho = .false.
      opark = .false.
      oreact = .true.
      odci = .false.
      odlfi = .false.
      if(opg_root()) then
       write(iwr,500)
 500   format(/1x,'**********************************'/
     +         1x,'**** GAMESS-UK Input Data Set ****'/
     +         1x,'**********************************'/)
       oecho = .true.
      endif
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
300   jrec=jrec-1
310   call prinit(oreact,opark,odci,odlfi)
      go to 320
 10   call input
 320  call inpa4(yyy)
      j=locatc(ytext,15,yyy)
      if(j)20,30,20
 20   go to (40,40,40,120,50,50,50,60,70,200,240,10,10,925,
     *       285), j
c
c logfile
c
 925  ologf=.true.
      go to 10
c
c width
c
c920  call inpi(jwidth)
c     call inpwid(jwidth)
c     go to 10
c
c  sort
c
 200   call inpi(nsz)
       call inpa(ztext)
       go to 10
c
c  psort
c
 240   call inpi(ns2)
       call inpa(ztext)
       go to 10
c
c  table
c
 285   call inpi(ns3)
       call inpa(ztext)
       go to 10
c
c    setio
c
 120   call inpi(nbuff)
       call inpi(ibuff)
      go to 10
c
c     =change/switch/file=
c
 40   if(jrec.ge.jump)go to 10
      call inpa(ztext)
c
c...  fortran sequential/vbs
c
      j=locatc(zzzz1,maxfrt,ztext)
      if(j.eq.0) then
         j=locatc(zzzz2,maxfrt,ztext)
         if(j.eq.0) then
            j=locatc(zzzz3,maxfrt,ztext)
            if(j.eq.0) go to 150
         endif
      endif
      iposf(j)=1
      keepf(j)=1
      call inpanpc(zftfil(j))
      call inpa(ztext)
      if(ztext(1:4).eq.'form') oform(j)=.true.
      go to 10
c
c...  512 word files
c
150   j=locatc(yed,maxlfn,ztext(1:4))
      if(j.eq.0) call caserr(
     *'invalid lfn specified in file directive')
      keep(j) = 1
      call inpa(zedfil(j))
      go to 10
c
c     =store/core/memo=
c
 50   call inpi(llll)
      if(jump.gt.2) then
 51    call inpa(ztext)
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
      go to 10
 70   call inpf(ffmin)
      insec=nint(ffmin*60.0d0)
      otime=.true.
      go to 10
c
 30   jrec=0
      if(nsz.le.0.or.nsz.gt.10)nsz=10
      if(ns2.le.0.or.ns2.gt.1)ns2=1
      if(ns3.le.0.or.ns3.gt.1)ns3=1
      if(otime) isecs=insec
      if(isecs.le.0)isecs=36000
      if(nbuff.le.0)nbuff=maxbuf
      if(ibuff.le.0)ibuff=maxblo
c
       call forstr
c
_IF(parallel)
c
c... set iwr to 94 to seperate output of nodes
c
      if (.not.opg_root()) then
         iwr=94
c...      different numbers are not needed/different names are
         nn = 100 + ipg_nodeid()
         write(fname,'(a,i3)')'out_gam_',nn
         open(iwr,form='formatted',file=fname,err=9998,
     +         status='unknown')
      end if
_ENDIF
      if(opun7)write(iwr,201)zinout(2)
201   format(/
     *' punchfile allocated : lfn = ',a8)
      return
_IF(parallel)
9998  write(6,9997)
9997  format(/' *** failed to open output file *** '/)
      call caserr('failure to open output file')
      return
_ENDIF
      end
c
c** cray initb
c
      subroutine initb
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      character*34 filpr
      integer getenv
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/discc1)
c
      odd=.false.
      do 100 i=1,maxlfn
      zied=yed(i)
      filatt = ' '
      ii = getenv(zedfil(i),filatt)
      if(zied.ne.zedfil(i).or.ii.ne.0 ) go to 200
 100  continue
      go to 600
 200  write(iwr,300)
 300  format(//1x,
     *'local     ',2x,'status',2x,'external file name'/
     *1x,'file name '/1x,80('*'))
      odd=.true.
      do 400 i=1,maxlfn
      zied=yed(i)
      filatt=zedfil(i)
      if(zied.eq.zedfil(i))then
       ii = getenv(zedfil(i),filatt)
       if(ii.eq.0)  go to 400
      endif
      inquire (file=filatt,exist=ok)
      if (ok) then
         zstat = 'old'
      else
         zstat = 'new'
      endif
c...    pretty print
      filpr = ' '
      ii = index(filatt,' ')
      if (ii.le.0) ii = 80
      if (ii.le.34) then
         filpr(34-ii+1:) = filatt
      else
         iii = index(filatt(ii-34+3:),'/') + ii-34+2 
         iii = max(iii,1)
         filpr(34-(ii-iii+3):) = '...'//filatt(iii:ii)
      end if
      write(iwr,500)yed(i),zstat,filpr
 400  continue
 500  format(1x,a4,8x,a7,2x,a80)
 600  do 700 i=1,maxfrt
      zied=zftn(i)
      filatt=' '
      ii = getenv(zzzz1(i),filatt)
      jj = 0
      kk = 0
      if (ii.eq.0) then
         jj = getenv(zzzz2(i),filatt)
         if (jj.eq.0) then
            kk = getenv(zzzz3(i),filatt)
            if (kk.eq.0) filatt = zftn(i)
         endif
      endif
      ijk = ii + jj + kk
      if(zied.ne.zftn(i).or.ijk.gt.0 )  go to 800
 700  continue
      return
 800  if(odd) then
       write(iwr,850)
 850   format(/)
      else
       write(iwr,300)
      endif
      do 900 i=1,maxfrt
      zied=zzzz3(i)
      filatt=zftn(i)
      if(zied.eq.zftn(i))then
        ii = getenv(zzzz1(i),filatt)
        jj = 0
        kk = 0
         if (ii.eq.0) then
            jj = getenv(zzzz2(i),filatt)
            if (jj.eq.0) then
               kk = getenv(zzzz3(i),filatt)
               if (kk.eq.0) filatt = zftn(i)
            endif
         endif
        ijk = ii + jj + kk
        if(ijk.eq.0) go to 900
      endif
      inquire (file=filatt,exist=ok)
      if (ok) then
         zstat = 'old'
      else
         zstat = 'new'
      endif
      write(iwr,501)zied,zstat,filatt
 501  format(1x,a8,9x,a25,2x,a4)
 900  continue
      return
      end
c
c** cray ed2mem
c
      subroutine ed2mem(nbasis,nat)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
INCLUDE(common/atmblk)
INCLUDE(common/iofile)
      common/vinteg/q(1)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
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
      write(iwr,400) mxblock,mxblock*512,mxblock*4/1000.0d0
      omem(3)=.true.
      maxbl(3) = mxblock
      n2last(1)= mxblock
c
200   return
300   format(1x,'requested memory not available (',i6,' words')
400   format(1x,'sufficient memory allocated for an mfile of',i10,1x,
     +'integral blocks',/,' (',i12,' words / ',f11.3,' Mbyte )')
      end
c
c** cray ed0mem
c
      subroutine ed0mem(iunit)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/segm)
INCLUDE(common/maxlen)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
     +  ,keep(maxlfn),keepf(3,maxfrt)
      common/vinteg/q(1)
INCLUDE(common/iofile)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/files)
c
c     allocate memory for dump or scratchfile if available
c
      iqaddr = 0
      iqqoff(iunit) = 0
      if(maxset(iunit).le.0) call caserr(
     + 'failure to allocate dumpfile memory')
      mxblock = maxset(iunit)
      need = mxblock * 512
      call getmem(need, q, iqaddr, iqqoff(iunit))
      if (iqaddr.eq.0.or.iqqoff(iunit).eq.0) then
       write(iwr,*)' requested memory not available (',need,' words)',
     1             ' so unit ',iunit,' resides on disk'
       omem(iunit)=.false.
      else
        write(iwr,*)' sufficient memory allocated for ',maxset(iunit),
     + ' 512 word blocks (unit =',iunit,')'
       omem(iunit)=.true.
      endif
      return
      end
c
c** cray termbf
c
      subroutine termbf
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/discc)
      common/bufnew/ibufa(512)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),ibfbas(maxbuf),igetm(maxbuf),
     *iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... to clear up all pending io for current buffer
c...
      k=istabf(mbuff)
      if(k)8,9,8
8     j=junibf(mbuff)
      j40=i40(j)
      call chekwr(itab(j40))
      if(irep)1,2,1
2     if(k)7,6,6
c... buffer in put mode
7     iostat(j)=2
      mzero=0
      istabf(mbuff)=0
      l=iputm(mbuff)
      m=ibfbas(mbuff)
      k=jpbbf(mbuff)*npbsz
3     if(l)99,4,99
99     nzero=leadz(l)
      k=k+nzero
      m=(nzero+mzero)*512+m
      l=shiftl(l,nzero)
      mzero=leadz(compl(l))
      call flushw(itab(j40),ibufa(m+1),mzero,k+1,j)
       l=shiftl(l,mzero)
      k=k+mzero
      call chekwr(itab(j40))
      if(irep)1,3,1
1     call iofail(yed(j),iostat(j))
4     if(length(j).lt.k)length(j)=k
6     istabf(iobuff(j))=0
      iobuff(j)=maxbuf+1
9     return
      end
c
c** cray clredx
c
      subroutine clredx
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     &iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... to clear up i/o on all atmol streams
c...
      kunkpb=0
      do 1 mbuff=1,npb
      call termbf
1     continue
      return
      end
c
c** cray put
c
      subroutine put(ic,nw,iunit)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      dimension ic(*)
      common/bufnew/ibufa(512)
INCLUDE(common/discc)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     * iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... drives asynchronous output
c...
      iun=iunit
      ibl=ipos(iun)
      ipb=(ibl-1)/npbsz
      iunipb=shiftl(iun,32).or.ipb
      j40=i40(iun)
      iat=ipb*npbsz
      iatm=ibl-iat
      masker=shiftl(1,64-iatm)
      if(iunipb.eq.kunkpb)goto 776
      kunkpb=iunipb
      mbuff=locat(junjpb,iunipb)
      if(mbuff)777,888,777
c... physical block not in buffer -- find lowest priority buffer
888   mbuff=locat(ipribf,npbm1)
      igetm(mbuff)=0
      call revpri
      call termbf
      jpbbf(mbuff)=ipb
      junibf(mbuff)=iun
      junjpb(mbuff)=iunipb
      goto 555
777   call revpri
c... physical block in buffer - check for i/o activity
776   if(istabf(mbuff))666,555,444
444   iobuff(iun)=maxbuf+1
      call chekwr(itab(j40))
      if(irep)1,555,1
1     call iofail(yed(iun),iostat(iun))
555   iputm(mbuff)=0
      istabf(mbuff)=-1
666   m=ibfbas(mbuff)
      k=iatm*512+m
      ibufa(k)=shiftl(nw,32).or.ibl
      if(nw.gt.0)call fmove(ic,ibufa(k-511),nw)
      l=iputm(mbuff).or.masker
      iputm(mbuff)=l
      igetm(mbuff)=igetm(mbuff).or.masker
      if(iatm.ne.npbsz)goto 999
c... switch to out status
      mzero=0
      istabf(iobuff(iun))=0
      call chekwr(itab(j40))
      if(irep)1,77,1
77    iostat(iun)=2
3      if(l)99,4,99
99     nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=shiftl(l,nzero)
      mzero=leadz(compl(l))
      call chekwr(itab(j40))
      if(irep)1,5,1
5     call flushw(itab(j40),ibufa(m+1),mzero,iat+1,iun)
       l=shiftl(l,mzero)
      iat=iat+mzero
      goto 3
4     istabf(mbuff)=2
      if(length(iun).lt.iat)length(iun)=iat
      iobuff(iun)=mbuff
999   ipos(iun)=ibl+1
      return
      end
c
c** cray find
c
      subroutine find(iunit)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      common/bufnew/ibufa(512)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
INCLUDE(common/discc)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... initiates asynchronous input
c...
      lun=iunit
      lbl=ipos(lun)
      lpb=(lbl-1)/npbsz
      l40=i40(lun)
      lunlpb=shiftl(lun,32).or.lpb
      lat=lpb*npbsz
      latm=lbl-lat
      if(lunlpb.eq.kunkpb)goto 776
      kunkpb=lunlpb
      mbuff=locat(junjpb,lunlpb)
      if(mbuff)777,888,777
c... physical block not in buffers
888   mbuff=locat(ipribf,npbm1)
      igetm(mbuff)=0
      lbfbas=ibfbas(mbuff)
      call termbf
      jpbbf(mbuff)=lpb
      junibf(mbuff)=lun
      junjpb(mbuff)=lunlpb
      call revpri
      goto 333
c... physical block in buffer-check for put status
777   call revpri
776   lbfbas=ibfbas(mbuff)
      if(istabf(mbuff))666,444,444
666   mzero=0
      iat=lat
      l=iputm(mbuff)
      m=lbfbas
      call chekwr(itab(l40))
      if(irep)1,77,1
77    istabf(iobuff(lun))=0
      iostat(lun)=2
3       if(l)99,4,99
99     nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=shiftl(l,nzero)
       mzero=leadz(compl(l))
      call chekwr(itab(l40))
      if(irep)1,5,1
1     call iofail(yed(lun),iostat(lun))
5     call flushw(itab(l40),ibufa(m+1),mzero,iat+1,lun)
       l=shiftl(l,mzero)
      iat=iat+mzero
      goto 3
4     if(length(lun).lt.iat)length(lun)=iat
      iobuff(lun)=mbuff
      istabf(mbuff)=2
444   l=shiftl(igetm(mbuff),latm-1)
      if(l)999,333,333
c... atmol block missing
333   m=length(lun)-lat-npbsz
      iat=npbsz
       istabf(iobuff(lun))=0
      call chekwr(itab(l40))
      if(irep)1,19,1
19     iostat(lun)=1
      if(m)21,20,20
21    iat=iat+m
      if(iat)1,1,20
20    call fillv(lun,itab(l40),ibufa(lbfbas+1),iat,lat+1)
      iobuff(lun)=mbuff
      istabf(mbuff)=1
      igetm(mbuff)=mask(iat)
999   lbuff=mbuff
       ipos(lun)=lbl+1
      return
      end
c
c** cray get
c
      subroutine get(ic,nw)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      dimension ic(*)
      common/bufnew/ibufa(512)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
INCLUDE(common/discc)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... completes asynchronous input
c...
      if(istabf(lbuff).ne.1)goto 2
      call chekrd(itab(l40),latm)
      if(irep)1,2,1
1     call iofail(yed(lun),1)
2     m=latm*512+lbfbas
      mword=ibufa(m)
      iblok=mword.and.37777777777b
      if(iblok.ne.lbl)goto 1
      nw=shiftr(mword,32)
      if(nw.gt.0)call fmove(ibufa(m-511),ic,nw)
      return
      end
c
c** cray fget
c
      subroutine fget(ic,nw,iunit)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      common/bufnew/ibufa(512)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
INCLUDE(common/discc)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
      dimension ic(*)
c...
c... synchronous input routine
c...
      lun=iunit
      lbl=ipos(lun)
      lpb=(lbl-1)/npbsz
      l40=i40(lun)
      lunlpb=shiftl(lun,32).or.lpb
      lat=lpb*npbsz
      latm=lbl-lat
      if(lunlpb.eq.kunkpb)goto 776
      kunkpb=lunlpb
      mbuff=locat(junjpb,lunlpb)
      if(mbuff)777,888,777
c... physical block not in buffers
888   mbuff=locat(ipribf,npbm1)
      igetm(mbuff)=0
      lbfbas=ibfbas(mbuff)
      call termbf
      jpbbf(mbuff)=lpb
      junibf(mbuff)=lun
      junjpb(mbuff)=lunlpb
      call revpri
      goto 333
c... physical block in buffer-check for put status
777   call revpri
776   lbfbas=ibfbas(mbuff)
      if(istabf(mbuff))666,444,444
666   mzero=0
      iat=lat
      l=iputm(mbuff)
      m=lbfbas
      call chekwr(itab(l40))
      if(irep)1,77,1
77    istabf(iobuff(lun))=0
      iostat(lun)=2
3       if(l)99,4,99
99     nzero=leadz(l)
      iat=iat+nzero
      m=(nzero+mzero)*512+m
      l=shiftl(l,nzero)
       mzero=leadz(compl(l))
      call chekwr(itab(l40))
      if(irep)1,5,1
1     call iofail(yed(lun),iostat(lun))
5     call flushw(itab(l40),ibufa(m+1),mzero,iat+1,lun)
       l=shiftl(l,mzero)
      iat=iat+mzero
      goto 3
4     if(length(lun).lt.iat)length(lun)=iat
      iobuff(lun)=mbuff
      istabf(mbuff)=2
444   l=shiftl(igetm(mbuff),latm-1)
      if(l)999,333,333
c... atmol block missing
333   m=length(lun)-lat-npbsz
      iat=npbsz
       istabf(iobuff(lun))=0
      call chekwr(itab(l40))
      if(irep)1,19,1
19     iostat(lun)=1
      if(m)21,20,20
21    iat=iat+m
      if(iat)1,1,20
20    call fillv(lun,itab(l40),ibufa(lbfbas+1),iat,lat+1)
      iobuff(lun)=mbuff
      istabf(mbuff)=1
      igetm(mbuff)=mask(iat)
999    ipos(lun)=lbl+1
      if(istabf(mbuff).ne.1)goto 220
      call chekrd(itab(l40),latm)
      if(irep)210,220,210
210   call iofail(yed(lun),1)
220   m=latm*512+lbfbas
      mword=ibufa(m)
      iblok=mword.and.37777777777b
      if(iblok.ne.lbl)goto 210
      nw=shiftr(mword,32)
      if(nw.gt.0)call fmove(ibufa(m-511),ic,nw)
      return
      end
c
c** cray shut
c
      subroutine shut
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer unlink
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
     + ,keep(maxlfn),
     +  keepf(maxfrt),iposf(maxfrt),oform(maxfrt)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
      common/blksiz/nsz(13),nsort(40)
      common/blksi2/ns2(13),nport(40)
      logical oex
       call clredx
      do 1 i=1,maxlfn
         call delfil(i)
 1    continue
c
c     delete fortran-files (if any)
c
      do 2 i=1,maxfrt
          if (keepf(i).ne.0) go to 2
          inquire(file=zftn(i),exist=oex)
          if (oex) then
            ii = unlink(zftn(i))
            if(ooprnt) then
              write(iwr,3) zftn(i)
 3            format(' file ',a8,' unlinked')
              if (ii.ne.0) write(iwr,4)zftn(i)
 4            format(' ******* error condition on unlink **** ',a)
            end if
         end if
 2    continue
c
c...   delete sort and psort files
c
      inquire(file='sort',exist=oex)
      if (oex) then
       write(zunlnk,'(5hfort.,i2)')nsort(1)
       iii = unlink(zunlnk)
       ii = unlink('sort')
      endif
      inquire(file='psort',exist=oex)
      if (oex) then
       write(zunlnk,'(5hfort.,i2)')nport(1)
       iii = unlink(zunlnk)
       ii = unlink('psort')
      endif
c
      return
      end
c
c** cray delfil
c
      subroutine delfil (iunit)
c...   more perm fix delete iunit (after closing) for cray (jvl)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      integer  unlink
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
     +           ,keep(maxlfn),keepf(maxfrt)
INCLUDE(common/maxlen)
      call shut1(iunit)
      if (omem(iunit).or.keep(iunit).ne.0) return
       filatt = zedfil(iunit)
       if(yed(iunit).eq.zedfil(iunit)) then
        ii = getenv(zedfil(iunit),filatt)
        if(ii.eq.0) filatt = zedfil(iunit)
      endif
      inquire(file=filatt,exist=oex)
      if (oex) then
      ii = unlink(filatt)
       if(ooprnt) then
           write(iwr,3)yed(iunit)
 3         format(' file ',a4,' unlinked')
           if (ii.ne.0) write(iwr,4)
 4         format(' ******* error condition on unlink *****')
       end if
      end if
      return
      end
c
c** cray shut1
c
      subroutine shut1(iunit)
c...   shut 1 file (added for cray - jvl)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/maxlen)
      i=iunit
      if(ipos(i).gt.0) then
       call clredx
       if(.not.omem(i)) then
         ii=i40(i)
         call vclose(itab(ii))
         write(zunlnk,'(5hfort.,i2)')nam(iunit)
         iii = unlink(zunlnk)
         if (ooprnt.and.iii.ne.0) write(iwr,44)zunlnk
 44      format(' ******* error condition on unlink ***** ',a)
       endif
       ipos(i) = 0
       iostat(i) = 0
       length(i) = 0
       if(ooprnt)write(iwr,3)yed(i)
 3     format(//' file ',a4,' closed')
      endif
      return
      end
c
c** cray iofail
c
      subroutine iofail(itext,itype)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character *(*) itext
      dimension ztext(3)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),
     * nam(maxlfn)
INCLUDE(common/iofile)
      data ztext/' input',' output',' search'/
      write(iwr,1)ztext(itype),itext,irep
 1    format(/a8,' error on ',a4,' reply word = ',o4)
      call errors(62)
      return
      end
c
c** cray search
c
      subroutine search(iblock,iunit)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 filatt
      integer getenv
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/utilc)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
INCLUDE(common/discc)
INCLUDE(common/maxlen)
c
      if(iblock.le.0) then
       write(iwr,4444) iblock,zedfil(iunit)
       call iofail(yed(iunit),3)
      endif
      if(iunit.le.0.or.iunit.gt.maxlfn) then
       write(iwr,5555) iunit
       call iofail(yed(iunit),3)
      endif
      if(ipos(iunit).eq.0) then
        if(omem(iunit)) then
         if(iunit.ne.3) call ed0mem(iunit)
        else
         filatt = zedfil(iunit)
         if(yed(iunit).eq.zedfil(iunit)) then
          ii = getenv(zedfil(iunit),filatt)
          if(ii.eq.0) filatt = zedfil(iunit)
         endif
         incr = (iunit-1)*40
         itab(incr) = nam(iunit)
         call vopen(itab(incr),filatt,length(iunit))
         if(irep.ne.0)  call iofail(yed(iunit),3)
         i40(iunit)=incr
         if(ooprnt)write(iwr,1)yed(iunit),length(iunit)
        endif
      endif
      ipos(iunit)=iblock
      return
  1    format(/' lfn ',a4,' allocated: current length =',
     * i5,' blocks')
4444   format(/10x,
     +        '**** illegal search for block ',i6,' on file ',a8)
5555   format(/10x,'**** illegal unit ',i6,' in search ')
      end
c
c** cray getcor
c
      subroutine getcor(nmaxly)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c  ==== NULL routine under unicos ==
      return
      end
c
c** cray tidajt
c
      subroutine tidajt(zdate,ztime,zaccno,zanam,jtime)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/timez)
      dimension ifac(50)
       data ifac/50*0/
      call date(idate)
      call clock(itime)
      timlim=jtime+0.1d0
      cpusec = timlim
      call timef(elapl)
      if(jtime.ne.0)go to 10
      timlim=cpusec
      jtime=nint(timlim)
 10   continue
      write(ztime,20)itime
      write(zdate,20)idate
      zanam  = '        '
      zaccno = '        '
 20   format(a8)
      return
      end
c
c** cray whtps
c
      subroutine whtps
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
INCLUDE(common/blksiz)
      common/blksi2/ns2(11),n2len,n2max
      common/blksi3/ns3(11),n3len,n3max
INCLUDE(common/discc)
INCLUDE(common/iofile)
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
2     write(iwr,6000)yed(i),ipos(i),length(i)
1     continue
6000  format(1x,a4,2i8 )
      if(nsmax.gt.0)write(iwr,6001)ztexts,nslen,nsmax
      if(n2max.gt.0)write(iwr,6001)ztextp,n2len,n2max
      if(n3max.gt.0)write(iwr,6001)ztextt,n3len,n3max
6001  format(1x,a6,i6,i8)
      return
      end
c
c** cray quit
c
      subroutine quit(ia,n)
      implicit REAL  (a-h,p-w),integer  (i-n),logical   (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *80 ia
      character *72 ia72
      ia72=ia(1:72)
      call shut
      call whtps
      call remark2(ia72)
      call abort
      return
      end
c
c** cray kilbuf
c
      subroutine kilbuf
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character  *4 (y)
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40),
     *istabf(maxbuf+1),ipribf(maxbuf),junibf(maxbuf),
     *jpbbf(maxbuf),junjpb(maxbuf),
     *ibfbas(maxbuf),igetm(maxbuf),iputm(maxbuf),
     *mbuff,kunkpb,lun,l40,latm,lbl,lbfbas,lbuff
c...
c... to clear up all io on atmol streams
c... and disconnect buffers
c...
      call clredx
      call setsto(npb,0,junjpb)
      return
      end
c
c** cray sptchk
c
      subroutine sptchk(text)
      implicit REAL  (a-h,o-z), integer (i-n)
      character*(*) text
INCLUDE(common/utilc)
      if(ologf) then
      l72=min(lstchr(text),72)
      call remark2 (text(1:l72))
      endif
      return
      end
c
c** cray strtrm
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
c** cray vopen
c
      subroutine vopen(ftable,filnam,nblock)
      implicit integer (a-z)
      character *(*) filnam
      character *8 zunlnk
      common /disc/incr,npb,npbsz,irep
INCLUDE(common/iofile)
INCLUDE(common/utilc)
      integer ftable(*)
c...
c....   stat : info in /usr/include/sys/stat.h
c....   assumed to be max 53 (ok for 6.0 + ...?)
c...
      integer buf(53),stat
c...
      logical ok
c
      inquire (file=filnam,exist=ok)
      if(ooprnt)write(iwr,*) ' opening ',filnam,ok
      if (.not. ok) then
        nblock = 0
      else
        icrap = stat(filnam,buf)
        nblock = buf(8) / 4096 - 1
        if(ooprnt)write(iwr,*) ' nblock ',nblock
      endif
      iblock = npbsz
      write(zunlnk,'(5hfort.,i2)')ftable(1)
      iii = unlink(zunlnk)
      call linkcc(ftable(1),filnam,80)
      call wopen(ftable(1),iblock,0,ierr)
c...   because wopen cannot handle file-names > 8 chars
      if (ierr.ne.0) then
        write(iwr,*) ' ierr in vopen from wopen ',ierr
        call caserr('wopen failed')
      endif
      if(ooprnt)write(iwr,1) filnam,iblock,nblock
1     format(1x,' vopen : ',a,2x,2i8)
      return
      end
c
c** cray vclose
c
      subroutine vclose(ftable)
      implicit integer (a-h,o-z)
INCLUDE(common/iofile)
      dimension ftable(*)
      call wclose(ftable(1),ierr)
      if (ierr.ne.0) then
        write(iwr,*) ' ierr in vclose from wclose ',ierr
        call caserr('wclose failed')
      endif
      return
      end
c
c** cray fillv
c
      subroutine fillv(iunit,ftable,buffer,mblock,iblock)
      implicit integer (a-h,o-z)
INCLUDE(common/sizes)
      common /disc/incr,npb,npbsz,irep
      common/vinteg/q(1)
INCLUDE(common/maxlen)
      dimension ftable(1),buffer(*)
      nword = mblock*512
      if(omem(iunit)) then
       iword = (iblock-1)*512 + 1
       iword = iword + iqqoff(iunit)
       call fmove(q(iword),buffer,nword)
       irep = 0
      else
       iword = iblock*512 + 1
       call seek(ftable(1),iword,nword,ierrs)
       if (ierrs.ne.0) then
         write(6,*) ' error from seek in fillv ',ierrs
         write(6,*) ' iblock, iword, nword ',iblock,iword,nword
         call caserr('seek failed')
       endif
       call getwa(ftable(1),buffer,iword,nword,ierrg)
       if (ierrg.ne.0) then
         write(6,*) ' error from getwa in fillv ',ierrg
         write(6,*) ' iblock, iword, nword ',iblock,iword,nword
         call caserr('getwa failed')
       endif
       irep = ierrs + ierrg
      endif
      return
      end
c
c** cray flushv
c
      subroutine flushv(ftable,buffer,mblock,iblock,iunit)
      implicit integer (a-h,o-z)
INCLUDE(common/sizes)
      common/vinteg/q(1)
INCLUDE(common/maxlen)
      common /disc/incr,npb,npbsz,irep
      dimension ftable(1),buffer(*)
      nword = mblock*512
      if(omem(iunit)) then
       iword = (iblock-1)*512 + 1
       iword = iword + iqqoff(iunit)
       call fmove(buffer,q(iword),nword)
       irep = 0
      else
       iword = iblock*512 + 1
       call putwa(ftable(1),buffer,iword,nword,ierrp)
       if (ierrp.ne.0) then
         write(6,*) ' error from putwa in flushv ',ierrp
         write(6,*) ' iblock, iword, nword ',iblock,iword,nword
         call caserr('putwa failed')
       endif
       irep = ierrp
      endif
      return
      end
      subroutine chekrd(f,i)
      return
      end
      subroutine chekwr(f)
      return
      end
c
c** cray catch
c
      subroutine catch
c
c     set up signal-catching  (cray-ymp unicos)
c
c  ===  UNIX signals for
c  ===  flushing o/p (sigusr1(49))  (not effective due to SSD)
c  ===  time dumping (sigusr2(50))
c
      external iflsh6,idumpt,caught
      call fsigctl('REGISTER','SIGUSR1',iflsh6)
      call fsigctl('REGISTER','SIGUSR2',idumpt)
c
c...  catch all kind of system stuff
c
      do 1 i=1,17
1     call fsigctl('REGISTER',i,caught)
      call fsigctl('REGISTER',22,caught)
      do 2 i=24,26
2     call fsigctl('REGISTER',i,caught)
c
      return
      end
c
c** cray caught
c
      subroutine caught(flop)
c
c...   to get a bit more info from system after an error
c...   (the error-number is passed by value / use loc to get it)
c
      character*60 text
      write(text,1) loc(flop)
1     format(' ****** UNIX error signal ',i5,' caught ******')
      call sync()
      call qquit(text)
      end
      subroutine qquit(ia)
c...    quit for recovered errors (e.g. io)
      character *(*) ia
      print *,' ********************************************'
      print *,ia
      print *,' ********************************************'
c...   first print file positions and traceback
      call whtps
      call tracebk
c...   try to flush dumpfile
      call revind
c...   now error-exit
      call remark2(ia)
      call exit(999)
      return
      end
c
c** cray cleanci
c
       subroutine cleanci(phase)
c
c...   cleanup files for big CI's
c...   reacts to setenv cleanci
c...   phase denotes the phase in the program
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
      common/disc/incr,npb,npbsz,irep,npbm1,ipos(maxlfn),nam(maxlfn)
      common/discne/i40(maxlfn),length(maxlfn),iostat(maxlfn),
     *iobuff(maxlfn),itab(maxlfn*40)
c...   sort-file
      common/blksiz/nsz(11),nslen,nsmax,nsort(40)
c...   psort-file
      common/blksi2/ns2(11),n2len,n2max,nport(40)
INCLUDE(common/iofile)
INCLUDE(common/discc)
c
       character*(*) phase
       integer getenv
c
      parameter (nphas=2)
      dimension zphase(nphas),numph(nphas),nedph(nphas,3)
      dimension ztexsp(2)
c...
      data ztexsp/'sort ','psort'/
c....  close and delete :
c....  after tran : ed2 ed4 sort
      data zphase(1)/'tran'/
      data numph(1)/3/
      data (nedph(1,i),i=1,3)/3,5,100/
c....  after sorting in direct : ed6 sort psort
      data zphase(2)/'dirsort'/
      data numph(2)/3/
      data (nedph(2,i),i=1,3)/7,100,101/
c
       if (getenv('cleanci',ztext).eq.0) return
c
      iphase = 0
      do 10 i=1,nphas
10    if (phase.eq.zphase(i)) iphase = i
      if (iphase.eq.0) call caserr(' program error1 cleanci')
c
      write(iwr,600)
      do 20 k=1,numph(iphase)
         i = nedph(iphase,k)
         if (i.le.maxlfn) then
            if (ipos(i).le.0) call caserr(' program error2 cleanci')
            write(iwr,601) yed(i),zedfil(i),length(i)
            ii=i40(i)
            call vclose(itab(ii))
            call unlink(zedfil(i))
            ipos(i) = 0
         else if (i.eq.100) then
            if(nsmax.le.0) call caserr(' program error3 cleanci')
            write(iwr,601) '    ',ztexsp(1),nsmax
            call vclose(nsort)
            call unlink(ztexsp(1))
            nsmax = 0
         else if (i.eq.101) then
            if(n2max.le.0) call caserr(' program error4 cleanci')
            write(iwr,601) '    ',ztexsp(2),n2max
            call vclose(nport)
            call unlink(ztexsp(2))
            n2max = 0
         else
            call caserr(' program error5 cleanci')
         end if
20    continue
      write(iwr,602)
c
600   format(/' **************************',
     *       /1x,' direct access files removed '/,
     *        3x,'lfn','   name   length'/1x,24('='))
601   format(2x,a4,2x,a8,i8 )
602   format(' **************************'/)
c
      return
      end
c
c** cray forstr
c
      subroutine forstr
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      character *132 filatt
      integer getenv
INCLUDE(common/iofile)
INCLUDE(common/utilc)
INCLUDE(common/discc)
INCLUDE(common/discc1)
      common/disc/isel(5+maxlfn*3),
     * keep(maxfrt),ipos(maxfrt),oform(maxfrt)
c
      do 1 i=1,maxfrt
      if(ipos(i).le.0.and.keep(i).le.0) go to 1
      iunit=i
      filatt=' '
      ii = getenv(zzzz1(iunit),filatt)
      if (ii.eq.0) then
         jj = getenv(zzzz2(iunit),filatt)
         if (jj.eq.0) then
            kk = getenv(zzzz3(iunit),filatt)
            if (kk.eq.0) filatt = zftn(iunit)
         endif
      endif
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
      if(ooprnt) write(iwr,3)zftn(i),filatt
3     format(1x,'error opening fortran file',1x,a8,
     * ' : filename : ',a)
      endif
      if(ooprnt) write(iwr,2)zftn(i),filatt
 2    format(/' fortran stream ',a8,' allocated: file name = ',a)
 1    continue
      return
      end
c
c** cray print0
c
      subroutine print0
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *132 rec132
INCLUDE(common/sizes)
INCLUDE(common/work)
INCLUDE(common/iofile)
_IF(parallel)
      dimension ibuff(132)
_ENDIF
      logical o_inclin,o_inc_err
      o_inclin = .false.
      o_inc_err = .false.
c
c --- write a copy of the input to unit 99
c
      open(unit=99,form='formatted',status='scratch')
1060  continue
_IF(parallel)
      if (opg_root()) then
        read(ird,1080,end=4790) rec132
        call chtoi(ibuff(1),rec132)
        goto 4791
4790    ibuff(1)=0
4791    continue
      endif
      call pg_brdcst(100,ibuff,132*8,0)
      if (ibuff(1).eq.0) goto 1090
      if (.not.opg_root()) call itoch(ibuff(1),rec132)
_ELSE
      if (.not.o_inclin) read(ird,1080,end=1090) rec132
_ENDIF
cjvl
c     include statement (allows including files in the input )
c
      call inclin(rec132,98,iwr,o_inclin,o_inc_err)
cjvl
      write(99,1080)rec132
1080  format(a132)
      go to 1060
1090  ird=99
      rewind ird
      return
      end
_ENDIF
_IF1()
_IF1() =============== end of CRAY section =================
_IF1()
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

_IF(titan,rs6000,hp700,hpux11,t3d,unix)
      Integer       length_name
_ENDIF
_IF(rs6000,hp700,hpux11,t3d,unix)
      Integer       nbytes
_ENDIF

_IF(titan)
      Integer stat, statb( 1:20 )
      Call strtrm( file, length_name )
      error  = stat( file( 1:length_name ), statb )
      blocks = ( statb( 8 ) / 8 - 1 ) / 512 + 1
_ELSEIF(rs6000,hp700,hpux11,t3d,unix)
      length_name = Len( file )
      Call statcc( file, length_name, nbytes )
      If( nbytes .EQ. -1 ) Then
         Call caserr('statcc error')
      End If
      blocks = ( nbytes / 8 - 1 ) / 512 + 1
      error  = 0
_ELSEIF(vms)
      Integer decc$stat,statb( 1:20 )
      error  = decc$stat( file, statb )
      blocks = ( statb( 8 ) / 8 - 1 ) / 512 + 1
_ELSEIF(mips4)
      Integer stat, statb( 1:23 ),statt(2)
      integer*8 off_t
      equivalence (off_t,statt)
      
      error  = stat( file, statb )
      statt(1) = statb(7)
      statt(2) = statb(8)
      blocks = ( off_t / 8 - 1 ) / 512 + 1
_ELSE
      Integer stat, statb( 1:20 )
      error  = stat( file, statb )
      blocks = ( statb( 8 ) / 8 - 1 ) / 512 + 1
_ENDIF

      End

c
c     Revised error handling
c     
      subroutine gamerr(desc, code, class, sync, sys)
      implicit none
      integer code, class, sync, sys, code2
      character*(*) desc

INCLUDE(common/errcodes)
INCLUDE(common/workc)
INCLUDE(common/work)
INCLUDE(common/iofile)
_IF(charmm,chemshell)
INCLUDE(common/chemshell)
_ELSEIF(taskfarm)
INCLUDE(common/taskfarm)
      character*132 infile, outfile
_ENDIF
_IF(xml)
c Need to include to see if we are outputting xml
INCLUDE(common/blocks)
INCLUDE(common/xmlout)
_ENDIF(xml)
      logical opg_root
      external opg_root

      integer ipg_nodeid
      external ipg_nodeid

      integer ipg_nnodes
      external ipg_nnodes

      integer length
      character pad*80
      character desc2*80
      
_IF(cio)
      logical odumpc
      common/crdmp/odumpc
_ENDIF(cio)

      logical ocaser
      common/iotrp/ocaser
c
_IF(xml)
      if (oxml) call xmlAddError( desc, code )
_ENDIF
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

_IF(cio)
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
_ENDIF(cio)
c     
c     exit
c     
      call clenup2
c
      if(sync .eq. ERR_SYNC .or. ipg_nnodes() .eq. 1)then
c
c clean exit will help ensure that the output is complete
c
_IF(altix)
c
c including to avoid hangs on altix   
c we shouldnt actually need this here but will leave it for now.
c
      if (ipg_nnodes().gt.1)then
        call abort('Aborting by FORCE - have a nice day ...')
      endif
_ENDIF

_IF(charmm)
         CALL WRNDIE(-5,'<GUKINI>',desc2)
_ELSE
         call pg_end(code2)
_ENDIF
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
_IF(linux)
_IF(f2c)
c
c linux (f2c) replacements for bit intrinsics 
c
      integer function ishft(in,bits)
      implicit none
      integer in,out,bits
      call ishftcc(in,out,bits)
      ishft = out
      end
_ENDIF
_ENDIF
_IF(linux,hitachi)
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
_ENDIF
_IF(linux)
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
_ENDIF
c
c** generic ioupdate
c
c  route selected files (GA and in eventually in memory)
c  back to disk
c
      subroutine ioupdate

      implicit none

INCLUDE(common/sizes)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/maxlen)
INCLUDE(common/ga_file)

_IF(ga)
c
c start with the dumpfile only
c
      integer unit, ga_to_use
      logical opg_root

      unit = idaf

      if(ioupd .ne. 0) then

      If( iwhere( unit ) .EQ. 6 .Or. iwhere( unit ) .EQ. 7 ) Then 

         ga_to_use = ga_file_unit_to_ga( unit )

          If( ga_file_keep( ga_to_use ) ) Then
              if (opg_root())then

                 write(iwr,100)zedfil(unit)(1:50)
 100  format(1x,'writing GA file to :',a50)

              endif
              Call dump_ga_file( unit, zedfil(unit) )
          End If

       endif
      endif

_ENDIF
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
_IF(parallel)
         return
_ELSE
      call caserr('include file not found')
_ENDIF
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

_IF(signals)
c
c** generic idumpt
c
      function idumpt()
      implicit REAL  (a-h,p-w),integer(i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/timez)
      integer*4 output_file
      parameter (output_file=6)
      timlim=-1.0d8
      write(6,20)
 20   format(/10x,' ****************************************'/
     1        10x,' *        a kill signal has been        *'/
     2        10x,' *               received               *'/
     2        10x,' *  attempt to initiate dumping process *'/
     3        10x,' ****************************************'/)
_IF1(c)      call sync()
      call flush(output_file)
      idumpt=0
      return
      end
c
c** generic iflsh6
c
      function iflsh6()
      implicit REAL  (a-h,o-z)
      character*8 zdate,ztime,zyear,zaccno,zanam
      integer*4 output_file
      parameter (output_file=6)
      call tidajt(zdate,ztime,zyear,zaccno,zanam,isecs)
      write(6,100) ztime,zdate(1:6),zyear,cpulft(1)
 100  format(/10x,' ***********************************************'/
     1        10x,' *    a kill -USR1 signal has been received    *'/
     2        10x,' * at ',a8,  ' on ',a6,2x,a4,' cpu',f12.2,   ' *'/
     3        10x,' *           flush output to unit 6            *'/
     4        10x,' ***********************************************'/)
_IF1(c)      call sync()
      call flush(output_file)
      iflsh6=0
      return
      end
c
c** generic icpu
c
      function icpu()
      integer*4 output_file
      parameter (output_file=6)
_IF1(c)      call sync()
      call flush(output_file)
      call caserr('CPU-timelimit EXCEEDED')
c
      icpu=0
      return
      end

_IFN(cray)
      subroutine catch

_IF(vms)
      external iflsh6,idumpt
      integer signal
      integer *4 isig4,min1
      character*(*) act
      data min1/-1/
      isig4 = 30
      icode=decc$signal(isig4,iflsh6,min1)
      isig4 = 31
      icode=decc$signal(isig4,idumpt,min1)
_ELSE
c
c  UNIX signal handling, see c.m
c  Actual numbers used are system dependent so use
c  kill -USR1 for iflsh6
c  kill -USR2 for idumpt
      call catchcc
_ENDIF
      return
      end
_ENDIF
_ENDIF
cjmht rewftn previously in dirctb.m
      subroutine rewftn(ntape)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(linux)
      integer error_code
      character *132 filatt
INCLUDE(common/sizes)
      common/disc/irep,iout,iin,iun,iblkk,ipos(maxlfn),
     * nam(maxlfn),keep(maxlfn),istat(maxlfn*2),io(maxlfn),
     * iwhere(maxlfn),isync(maxlfn),
     * iposf(maxfrt),keepf(maxfrt),istatf(3,maxfrt),
     * oform(maxfrt)
INCLUDE(common/discc)
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
_ELSE
      rewind ntape
_ENDIF
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
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/direc)
INCLUDE(common/machin)
INCLUDE(common/work)
INCLUDE(common/vectrn)
INCLUDE(common/prints)
_IF(ipsc)
INCLUDE(common/ipsc)
_ENDIF
c
_IF(parallel)
INCLUDE(common/disc)
_ELSEIF(cray,unicos)
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     + nam(maxlfn),keep(maxlfn),
     + keepf(maxfrt),iposf(maxfrt),oform(maxfrt)
_ELSE
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     * nam(maxlfn),keep(maxlfn),istat(maxlfn),llll(maxlfn),
     * iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     * iposf(maxfrt),keepf(maxfrt),istatf(maxfrt),
     * lrecf(maxfrt),nrecf(maxfrt),oform(maxfrt)
_ENDIF
INCLUDE(common/iofile)
INCLUDE(common/blocks)
_IF(xml)
INCLUDE(common/xmlin)
INCLUDE(common/xmlout)
_ENDIF
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
_IF(xml)
      if(ytest.eq.'xmlo')then
       oxml=.true.
       go to 1
      endif
      if(ytest.eq.'xmli')then
ccc       oxml_in=.true.
       go to 1
      endif
_ENDIF(xml)
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
_IF(ipsc)
      ohost(index(1)) = .true.
      ohost(index(3)) = .true.
_ENDIF
      endif
      if(opg_root()) then
_IF(cray)
        nfblck = 94
_ELSE
        nfblck = 58
_ENDIF
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
