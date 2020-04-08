      subroutine dieselinp
c
      implicit NONE
c
      character*256 line
c
c
      integer idiesin,idiesop
      parameter (idiesin=37,idiesop=78)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer len_trim
c
      open(unit=idiesop,file='diesel.options',form='formatted',
     &     status = 'replace')
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
c
      integer idiesin,idiesop
      parameter (idiesin=37,idiesop=78)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer len_trim
c
      open(unit=idiesop,file='diesel.options',form='formatted', 
     &     status='old')
10    read (idiesop,'(a)', end=1000) line 
      if (line(1:3).ne.'end') then
         program=line
         open(unit=idiesin,file='diesel.input',form='formatted',
     &        status='replace')
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer len_trim
c
      path=' '
      call getarg(0,path)
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
      call caserr('DIESEL not available')
      return 
      end
