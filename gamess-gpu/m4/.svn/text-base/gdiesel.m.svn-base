      subroutine dieselinp
c
      implicit NONE
c
      character*256 line
c
INCLUDE(common/diesel)
INCLUDE(common/iofile)
      integer len_trim
c
      open(unit=idiesop,file='diesel.options',form='formatted',
_IFN(mips4)
     &     status = 'replace')
_ELSE
     &     status = 'unknown')
_ENDIF
10    read (ird,'(a)', end = 1000) line
      if (line(1:9).ne.'dieselend') then
         write(idiesop,'(a)') line(1:len_trim(line))
         go to 10
      end if
      close(idiesop)
      return
1000  write (iwr,'(a)') 
     &      'unexpected end of file while reading diesel input'      
      call caserr('close diesel input section with dieselend')
      end

      subroutine diesel
c
      implicit NONE
c
      character*256 line, program
c
INCLUDE(common/diesel)
INCLUDE(common/iofile)
      integer len_trim
c
      open(unit=idiesop,file='diesel.options',form='formatted', 
     &     status='old')
10    read (idiesop,'(a)', end=1000) line 
      if (line(1:3).ne.'end') then
         program=line
         open(unit=idiesin,file='diesel.input',form='formatted',
_IFN(mips4)
     &        status='replace')
_ELSE
     &        status='unknown')
_ENDIF
20       read(idiesop,'(a)', end=1010) line
         if (line(1:3).ne.'end') then
            write(idiesin,'(a)') line(1:len_trim(line))
            go to 20
         end if      
         close(idiesin)
         call startdiesel(program)
      end if
      goto 10
1000  write (iwr,'(a)') 'finished running all specified diesel code'
      return
1010  call caserr('unexpected end of diesel input')
      end

      subroutine startdiesel(program)
c
      implicit NONE
c
      character*256 program, pcall, path
      integer l, iret, ip
      logical ofound
c
INCLUDE(common/iofile)
      integer len_trim
c
_IF(itanium)
       integer*4 null_4
       data null_4/0/
_ENDIF
_IF(diesel)
_IFN(macosx)
       external system
       integer system
_ENDIF
_ENDIF
      path=' '
_IFN(win95)
_IF(convex,sgi,sun,dec,hp700,hpux11,rs6000,ksr,itanium,linux)
_IF(itanium)
      call getarg(null_4,path)
c     call getarg(0,path)
_ELSE
      call getarg(0,path)
_ENDIF
_ENDIF
_ENDIF
      print *,'GAMESS is', path
      ip = len(path)+1
      ofound = .false.
      do while (.not.ofound.and.(ip.gt.1))
         ip = ip - 1
         if (path(ip:ip).eq."/") ofound=.true.
      end do

      l = len_trim(program)
      if (ofound) then 
         write (pcall,'(4a)') path(1:ip),'diesel/',program(1:l),
     &                        ' < diesel.input' 
      else 
         write (pcall,'(2a)') program(1:l),
     &                        ' < diesel.input'
      end if
      write (iwr,'(2a)') 'will call ',pcall
_IF(diesel)
_IF(macosx) 
      call system(pcall,iret)
_ELSE
      iret = system(pcall)
_ENDIF
      print *,'Return value ',iret
_ELSE
      call caserr('DIESEL not available')
_ENDIF
      return 
      end
_IF(mips4)
      integer function len_trim(string)
c...  fortran 90 intrinsic
      character*(*) string
      ll = len(string) + 1
10    ll = ll - 1
      if (string(ll:ll).eq.' '.and.ll.gt.1) go to 10
      if (string(ll:ll).eq.' ') ll = 0
      len_trim = ll
      return
      end
_ENDIF
