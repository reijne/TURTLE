      program gamessfarm
      implicit none
      include 'mpif.h'
INCLUDE(common/iofile)
INCLUDE(common/errcodes)
INCLUDE(common/taskmain)

      integer root, tag
      parameter (root=0, tag=10)

      character*132 infile,outfile
      integer comm, ntasks, i, j, k, l, group_done
      integer ierr, ioerr,status(MPI_STATUS_SIZE)
      integer init, gamret, evengrps, length, smry
      integer rtn_code, maxjobs, njob
      parameter (maxjobs=5000)
      integer secs
      logical looping, rootloop, groupsize_one
      logical finished,reading
      intrinsic len, mod
      integer nsuccess, nfail, ncrash, jobdone, nbadfiles
      character*132 filename, addsuffix
      character*132 filelog(maxjobs),badfiles(maxjobs)
      character*132 tmplist(maxjobs),tasklist(maxjobs)
      character*132 gamessin, gamessout, msg, grp2file(maxgroupsize)
      character*80 errmsg
      character *8 zyear,zdate,ztime,zaccno,zanam, job_start(maxjobs)
      common /filenames/ gamessin, gamessout
      REAL totaltime
      REAL  rgroup_done, rngroups
      integer idleg
      character*4 ytext

INCLUDE(common/taskfarm)
c     To use the GAMESS-UK the input reader
INCLUDE(common/workc)
INCLUDE(common/work)

      data filename /' '/
      data ntasks /0/
      data idleg /100/
      data length /0/
      data totaltime /0.0d0/
      data smry /10/
      data nsuccess /0/
      data ncrash /0/
      data nfail /0/

c     REM: root is the overall master, the processors are then grouped into
c     'satellites' - a group of processors under the control of their own
c     local master.
c     There is a separate 'master' group that links the root and 
c     the masters


c     Start MPI
      call MPI_INIT(ierr)

c     Determine number of processors and groups
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

c     The rank of the particular process in question
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

c     Check if we're running on > maxproc processors
      if (nprocs .gt. maxproc) then
         msg = 'Too many processors! Please change maxproc!'
         call taskerr(msg)
      elseif (nprocs .le. 1) then
         msg = 'Taskfarming on 1 processor?!? I think not...'
         call taskerr(msg)
      endif

c     Open the task.list file with the taskfarm directives and tasks
      if (rank .eq. root) then
         ioerr=0
         open(unit=ird, file='task.list', status='old', iostat=ioerr)
         if (ioerr .ne. 0) then
            msg = 'Taskfarm needs a file called task.list with a list \n
     +of the tasks to be processed !'
            call taskerr(msg)
         endif
         
c     read in the groupsize
         call input
         call inpa4(ytext)
         if (ytext .eq. 'grou')then
            call inpi(groupsize)
            if ((groupsize .gt. (nprocs-1)).or.(groupsize .lt.1)) then
               write(unit=msg,fmt=113)groupsize,nprocs
               call taskerr(msg)
            endif
         else
            msg = "First line in task.list specifys the group size\n
     +as in: \"groupsize x\" where x is the integer size."
            call taskerr(msg)
         endif

c     Check if there is  an idlegroup field, else assume its the first file
         call input
         call inpa4(ytext)
         if (ytext .eq. 'idle')then
            call inpi(idleg)
            if ( ( idleg .gt. 100 ).or.( idleg .le. 0 ) ) then
               write(unit=msg,fmt=115)idleg
               call taskerr(msg)
            endif
         elseif( ytext .eq. '*****' ) then
            msg = "2nd line in task.list is **** indicating EOF!"
            call taskerr(msg)
         else
            filename = char1c
c     Warn about lengthy filenames
            call strtrm(filename,length)
            if (length .gt. 124) then
               write(iwr,118)
            endif
            tmplist(1)=filename(1:length)
            ntasks=1
         endif

c     Now just loop over all  lines in the file until we are done
         call input
10171    continue
c     check oswit to see if we've reached the end of the file
         if( oswit ) goto 20
         if(jrec .ge. jump) call input
         call inpa4(ytext)
         if(ytext .eq. '****')then
            goto 20
         else
            filename = char1c
c     Warn about lengthy filenames
            call strtrm(filename,length)
            if (length .gt. 124) then
               write(iwr,118)
            endif
            tmplist(ntasks+1)=filename(1:length)
            ntasks=ntasks+1
            goto 10171
         endif
         
 20      continue
         close(unit=ird)

         if (ntasks .eq. 0) then
            msg = ("**Hey! Am I supposed to guess what to do here?\n
     &         Please ensure that task.list has a list of jobs in it.")
            call taskerr(msg)
         endif

c     The filenames have been read into the array tmplist, we now go through
c     this array and check we can open the files. The good files are added to
c     the array tasklist and the dodgyones to badfiles, which is printed to the
c     taskfarm summary file.

         j=0
         nbadfiles=0
         do i=1, ntasks
            filename = addsuffix( tmplist(i), '.in' )
            open(unit=smry, file=filename, status='old', iostat=ioerr )
            if ( ioerr .ne. 0 ) then
               write(*,*)"*** Could not open file ",filename
               nbadfiles = nbadfiles + 1
               badfiles(nbadfiles)= filename
            else
               j = j +1
               tasklist(j) = tmplist(i)
            endif
            close(unit=smry)
         enddo

         ntasks = j

      endif

c     Now sort every1 out with how big/many groups we have
      call MPI_BCAST(groupsize, 1, MPI_INTEGER,
     +     0, MPI_COMM_WORLD,ierr)
      ngroups=(nprocs-1) / groupsize

c     split the processors into groups under control of a root server
      call split_processors()
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     Commands for the Root node
      if (rank .eq. root) then
         call walltime(totaltime)
         open(unit=smry,file='taskfarm.summary',status='unknown'
     +        ,iostat=ioerr)
         if (ioerr .ne. 0) then
            write(*,*)"Could not open summary for writing!"
         endif

c     Write summary header and proc & groupsize info
         write(smry,201)
c     State which files we could not open
         if ( nbadfiles > 0 ) then
            write(smry,215)
            do i=1, nbadfiles
               write(smry,216) badfiles(i)
            enddo
         endif

c        Now write out the start time, nprocs & groupsize
         write(iwr,114) nprocs,groupsize
         write(smry,114) nprocs,groupsize
         call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
         write(smry,203)ztime,zdate(1:6),zyear(1:4)
         write(smry,204)

c        Loop to send out the first job to each of the masters
         njob=0
         group_done=0
         finished = .false.
         do i=1, ngroups
            if (finished) then
               filename(1:5)="*****"
            else
               call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
               filename=tasklist(njob+1)
               write(iwr,119)i,filename(1:20),ztime,zdate(1:6),
     &                 zyear(1:4)
               job_start(njob+1)=ztime
            endif

            call MPI_SEND(filename, 132, MPI_CHARACTER, i,
     +           tag, root_comm, ierr)

            if (finished) then
               group_done=group_done+1
            else
               njob=njob+1
c     grp2file: array for taskname indexed master rank in rgroup
               grp2file(njob) = filename
               if (njob .eq. ntasks) then
                  write(iwr,120)ztime,zdate(1:6),zyear(1:4)
                  finished = .true.
                  filename(1:5)="*****"
               endif
            endif
         enddo
c        end of loop sending out the first inputs


c     Now the loop to receive the result, determine source and send out new filename
c     to the masters - this will run until either all the jobs have finished and all 
c     the masters have recieved the kill signal or the percentage of idle groups has
c     increased above the value specified for idleg

         rootloop=.true.
c        jobdone to count completed jobs for summary
         jobdone = 0
         do while (rootloop)
            call MPI_RECV ( rtn_code, 1 , MPI_INTEGER, MPI_ANY_SOURCE,
     &           tag, root_comm, status, ierr)
            k = status(MPI_SOURCE)
            jobdone = jobdone + 1
            
c           state when and how things went
            do l=1,njob
               if (tasklist(l) == grp2file(k)) then
                  call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
                  if (rtn_code == 0) then
                     nsuccess=nsuccess+1
                     write(smry,205)tasklist(l),job_start(l),ztime
                     write(iwr,205)tasklist(l),job_start(l),ztime
                  elseif (rtn_code == 1) then
                     nfail=nfail+1
                     write(smry,206)tasklist(l),job_start(l),ztime
                     write(iwr,206)tasklist(l),job_start(l),ztime
                  elseif (rtn_code  == 2) then
                     ncrash=ncrash+1
                     write(smry,207)tasklist(l),job_start(l),ztime
                     write(iwr,207)tasklist(l),job_start(l),ztime
                  else
                     write(iwr,*)"Missed a result"
                  endif

               endif
            enddo

            if (.not. finished) then
               filename=tasklist(njob+1)
               call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
               write(iwr,119)k,filename(1:20),ztime,
     &              zdate(1:6),zyear(1:4)
               job_start(njob+1)=ztime

            endif

            call MPI_SEND(filename, 132, MPI_CHARACTER,k,
     +           tag, root_comm,ierr)

            if (finished) then
c     make sure all masters get the kill signal
               group_done=group_done+1
               rgroup_done=real(group_done)
               rngroups=real(ngroups)
               if ( ANINT( ( rgroup_done/rngroups)*100 ).ge.idleg) then
c       Finished so write out the summary
                  call walltime(totaltime)
                  write(smry,208)totaltime
                  write(smry,209)jobdone,nsuccess,nfail,ncrash
                  call flushn(smry)
                  close(unit=smry)

c     Can only finish cleanly if all jobs are allowed to finish - otherwise messy abort
                  if ( idleg .eq. 100 ) then
                     rootloop =.false.
                  else
                     call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
                  endif
               endif
            else
               njob=njob+1
               grp2file(k)=filename
               if (njob .eq. ntasks) then
                  write(iwr,120)ztime,zdate(1:6),zyear(1:4)
                  finished = .true.
                  filename="*****"
               endif
            endif

         enddo
      endif
c     end of root if

c     Loop for the masters and slaves
      if (rank .ne. root) then
         looping = .true.
         do while (looping)
            if (masters(rank)) then
               call MPI_RECV( filename, 132, MPI_CHARACTER, root,
     +              tag, root_comm, status, ierr)
            endif

            if (.not. groupsize_one()) then
               call MPI_BCAST(filename, 132, MPI_CHARACTER,
     +              0, group_comm, ierr)
            endif

c           run the gamess subroutine
            if (filename(1:5) .eq. '*****') then
               looping = .false.
            else
               gamessin=addsuffix(filename,'.in')
               gamessout=addsuffix(filename,'.out')
c              Set the names of ed's 2,3, & 7
               call setenv_tskfm(filename)
               init=1
               call rungamessf(init,gamret)
c               write(*,*)"ran job ",filename
c               gamret=0

               if (.not. groupsize_one()) then
                  call MPI_BARRIER(group_comm, ierr)
                  call MPI_REDUCE(gamret, rtn_code, 1, MPI_INTEGER,
     +                 MPI_MAX, root, group_comm, ierr)
               else
                  rtn_code = gamret
               endif

c     Delete the edfiles
               if (masters(rank)) then
                  call del_edfiles(filename,odscf1)
               endif

               if (masters(rank)) then
                  call MPI_SEND( rtn_code, 1, MPI_INTEGER, root,
     +                 tag, root_comm, ierr)
              endif
            endif
         enddo
      endif


      call MPI_FINALIZE(ierr)

c     FORMAT Statements
 110  format(a10,i4)
 111  format(a7,i4)
 112  format(132a)
 113  format("Groupsize of ",i5," incompatible with ",i6,
     +        " processors.")
 114  format(/5x,"Taskfarm running on ",i5,
     +     " procs with ",i4," procs per group.")
 115  format("Invalid idleg of ",i4,"  must be 0 < idleg <= 100")
 118  format(/"**Warning! - filename is > 124 characters."/
     &     "Filenames will be trimmed and files could get overwritten.")
 119  format(/"Group ",i4," processing task: ",a20," at ",a8,
     &              " on ",a6,1x,a4)
 120  format(/"Taskfarming stopped reading inputs at: ",a8,
     +                 " on ",a6,1x,a4)
C  201  format(/38x,23('*')/40x,'Taskfarming Summary',
C      &     //38x,23('*')
C      &     //5x,'Taskfarming started at : ',a8,' on ',a6,1x,a4,
C      &     //5x,'Job Name',32x,5x,'Started',6x,'Finished',5x,'Result',
C      &     /5x,8('-'),37x,7('-'),6x,8('-'),5x,6('-'))
 201  format(/38x,23('*')/40x,'Taskfarming Summary',
     &     //38x,23('*'),// )

 203  format( /5x,'Taskfarming started at : ',a8,' on ',a6,1x,a4,/ )
 204  format(/5x,'Job Name',32x,5x,'Started',6x,'Finished',5x,'Result',
     &     /5x,8('-'),37x,7('-'),6x,8('-'),5x,6('-'))
 215  format( /5x,'*** WARNING! The following file(s) in task.list', 
     &     ' could not be opened!',/)
 216  format(10x,a40)



 205  format(5x,a40,5x,a8,5x,a8,5x,'Succeeded')
 206  format(5x,a40,5x,a8,5x,a8,5x,'Failed')
 207  format(5x,a40,5x,a8,5x,a8,5x,'**Crashed**')


 208  format(/5x,'Total runtime: ',f7.2,' seconds')

 209  format(/5x,i4,' jobs processed'/5x,19('-')/
     +     5x,'Succeeded : ',i4/
     +     5x,'Failed    : ',i4/
     +     5x,'Crashed   : ',i4)


      end

      subroutine split_processors()
      implicit none
      include 'mpif.h'
INCLUDE(common/iofile)
INCLUDE(common/errcodes)
INCLUDE(common/taskmain)

      integer ierr, i, j, k,workerprocs, evengrps, colour

      integer tmp, proc_start, proc_end, system_context
      integer procs_in_context( 1:maxgroupsize )


c     Work out if we've even or uneven groups
      workerprocs=(nprocs-1)
      extraprocs = mod(workerprocs,groupsize)

      if (extraprocs .ge. 1) then
         ngroups = (workerprocs/groupsize) +1
         evengrps = ngroups-1
      else
         ngroups = (workerprocs/groupsize)
         evengrps = ngroups
      endif

c     Create the array defining the satellites
c     rank2group: array linking the rank of a processor to its group 
c     rgroup: array to define the root group communicator this is
c             1 if true or MPI_UNDEFINED if not.
c     masters: true if we are not a master, false if not

_IF(newscf)
      call blacs_get( tmp, 0, system_context )
_ENDIF
      k=1
      rank2group(0)=MPI_UNDEFINED
      do i=1, evengrps
         proc_start = k
         do j=1, groupsize
            rank2group(k)=i
            procs_in_context( j ) = k
            k=k+1
         enddo
         proc_end = k-1
_IF(newscf)
         tmp = system_context
         Call blacs_gridmap( tmp, procs_in_context, groupsize, 
     +        groupsize, 1 )
         If( rank >= proc_start .And. rank <= proc_end ) Then
            grp_context = tmp
         End If
_ENDIF
      enddo
      extragroup=0
      if (extraprocs .ne. 0) then
         extragroup=evengrps+1
c     Tak on extra group using last k value as start
         do i=1, extraprocs
            proc_start = k
            rank2group(k)=evengrps+1
            procs_in_context( i ) = k
            k=k+1
         enddo
         proc_end = k - 1
_IF(newscf)
         tmp = system_context
         Call blacs_gridmap( tmp, procs_in_context,
     &        extraprocs, extraprocs, 1 )
         If( rank >= proc_start .And. rank <= proc_end ) Then
            grp_context = tmp
         End If
_ENDIF
      endif
         
c     Create the array defining the root/masters group
c     and the array defining the masters.
      rgroup(0)=1
      k=1
      do i=1, evengrps
         do j=1, groupsize
            if (j .eq. 1) then
               rgroup(k) = 1
               masters(k)=.true.
            else
               rgroup(k)=MPI_UNDEFINED
               masters(k)=.false.
            endif
            k=k+1
         enddo
      enddo
      if (extraprocs .ne. 0) then
c     Tack on extra group using last k value as start
         do i=1, extraprocs 
            if (i .eq. 1) then
               rgroup(k)=1
               masters(k)=.true.
            else
               rgroup(k)=MPI_UNDEFINED
               masters(k)=.false.
            endif
            k=k+1
         enddo
      endif

c     Create the communicators for the satellites
      colour=rank2group(rank)
      call MPI_Comm_Split( MPI_COMM_WORLD, colour,
     +     rank, group_comm, ierr )

c     Create the communicator for the root/masters group
      colour=rgroup(rank)
      call MPI_Comm_Split( MPI_COMM_WORLD, colour,
     +        rank, root_comm, ierr )

      end

      logical function groupsize_one()
      implicit none
INCLUDE(common/taskmain)

      if ((groupsize .eq. 1) .or. ((extraprocs .eq. 1) .and.
     &        (rank2group(rank) .eq. extragroup))) then
            groupsize_one=.true.
      else
         groupsize_one=.false.
      endif
      return
      end
            
      integer function group_id()
      implicit none
      include 'mpif.h'
INCLUDE(common/taskmain)
      
      group_id = rank2group(rank)
      if(group_id .eq.MPI_UNDEFINED)group_id = 0
      return
      end

      subroutine get_communicator(comm)
      implicit none
INCLUDE(common/taskmain)

      integer comm

      comm=group_comm
      return
      end

      subroutine get_context(context)
      implicit none
INCLUDE(common/taskmain)

      integer context
     
      context=grp_context
      return
      end

      subroutine getfilenames(infile,outfile,lenfile)
      implicit none
      character*132 infile,outfile
      integer comm
      common /filenames/ gamessin,gamessout
      character*132 gamessin,gamessout
      integer lenfile

      infile=gamessin
      outfile=gamessout
      return
      end

      function addsuffix(filename, suffix)
      implicit none
      character*132 filename, addsuffix
      character*(*) suffix
      integer i, length, len_suff
      intrinsic len
c     strip off trailing blanks from filename and add the suffix

c     strtrm (machscf) returns the length excluding blanks
      call strtrm(filename,length)

c     need to make sure there is enough room for the suffix and to
c     add "./" for del_edfiles
      if (length .gt. 126) then
         length = 126
      endif
      len_suff=len(suffix)
      addsuffix=filename(1:length)//suffix(1:len_suff)
      return
      end

      subroutine rungamessf(init,gamret)
c     This calls the c rungamess routine that uses the long_jmp
      implicit none
      integer init, gamret
      call rungamessc(init, gamret)
c      write(*,*) "fortran rungamessf returned with ",gamret
      return
      end

      subroutine taskerr(msg)
c     Bail out if somethings gone horribly wrong.
      implicit none
      include 'mpif.h'
INCLUDE(common/taskmain)

      integer l, ierr, err_rtn
      intrinsic len
      character*(*) msg

      err_rtn=1

      if (rank .eq. 0) then
         l = len(msg)
         call stderrc(msg,l,0,0,rank)
      endif
      call MPI_ABORT(MPI_COMM_WORLD,err_rtn,ierr)
      stop
      return
      end

      subroutine setenv_tskfm(filename)
      implicit none
      character*(*)filename
      character*3 ed2_envname,ed3_envname,ed7_envname
      character*6 pun_envname
      character*132 ed2_value,ed3_value,ed7_value,pun_value,addsuffix
      integer led2_envname,led3_envname,led7_envname,lpun_envname
      integer led2_value,led3_value,led7_value,lpun_value
      intrinsic len

      ed2_envname="ed2"
      ed3_envname="ed3"
      ed7_envname="ed7"
      pun_envname="ftn058"
      led2_envname=3
      led3_envname=3
      led7_envname=3
      lpun_envname=6

      ed2_value = addsuffix(filename,'.ed2')
      ed3_value = addsuffix(filename,'.ed3')
      ed7_value = addsuffix(filename,'.ed7')
      pun_value = addsuffix(filename,'.pun')
      led2_value = len(ed2_value)
      led3_value = len(ed3_value)
      led7_value = len(ed7_value)
      lpun_value = len(pun_value)

      call setenvc(ed2_envname,led2_envname,ed2_value,led2_value)
      call setenvc(ed3_envname,led3_envname,ed3_value,led3_value)
      call setenvc(ed7_envname,led7_envname,ed7_value,led7_value)
      call setenvc(pun_envname,lpun_envname,pun_value,lpun_value)

      return
      end

      subroutine del_edfiles(filename,odscf)
      implicit none
INCLUDE(common/iofile)
ccccINCLUDE(common/cslosc)
INCLUDE(common/sizes)
INCLUDE(common/maxlen)
      character*132 filename
      character*132 ed2name,ed3name,ed7name,addsuffix
      integer led2name, led3name, led7name, ierr, length
      integer i, ipg_nnodes, ipg_nodeid
      character num*3
c      intrinsic len
      logical odscf
c     Delete the ed files

      if(.not. odscf .and. .not. omem(3) )then
         do i = 1,ipg_nnodes()
            ed2name = addsuffix(filename,'.ed2')
            write(num,'(i3.3)')i-1
            ed2name = addsuffix(ed2name,num)
            call strtrm(ed2name,led2name)
cc            write(6,*) "Deleting ",ed2name
            call delcc(ed2name,led2name,ierr)
         enddo
      endif

      ed3name = addsuffix(filename,'.ed3')
      call strtrm(ed3name,led3name)
cc      write(6,*) "Deleting ",ed3name
      call delcc(ed3name,led3name,ierr)

      ed7name = addsuffix(filename,'.ed7')
      call strtrm(ed7name,led7name)
cc      write(6,*) "Deleting ",ed7name
      call delcc(ed7name,led7name,ierr)

      return
      end

      subroutine freaderr(ioerr)
c     set the filename to ***** if we hit errors on reading
      implicit none
      integer ioerr
INCLUDE(common/iofile)
      
      if (ioerr .lt.0) then
         write(iwr,*)"End of task.list but no file called *****!"
         write(iwr,*)"Setting filename to *****"
      elseif (ioerr .gt. 0) then
         write(iwr,*)"Error trying to read task.list!"
         write(iwr,*)"Setting filename to *****"
      endif

      return
      end

      subroutine print_jobresult(rtn_code,k,grp2file)
c     print a summary saying how and when the job went
      implicit none
INCLUDE(common/iofile)
      integer maxgroupsize,k,secs,rtn_code
      parameter (maxgroupsize=1024)
      character*132 grp2file(maxgroupsize)
      character *8 zyear,zdate,ztime,zaccno,zanam

      call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
      if (rtn_code .eq. 0) then
         write(iwr,117) k,grp2file(k)(1:20),ztime,zdate(1:6), 
     &        zyear(1:4)
      elseif (rtn_code .eq. 1) then
         write(iwr,119) k,grp2file(k)(1:20),ztime,zdate(1:6),
     &        zyear(1:4)
      else
         write(iwr,113) k,grp2file(k)(1:20),ztime,zdate(1:6),
     &        zyear(1:4)
      endif
            
 117  format(/"Group ",i4," returned successfully from task: ",a20,
     &        "at ",a8," on ",a6,1x,a4)
 119  format(/"** Group ",i4," was unsuccessful for task: ",a20,
     &        "at ",a8," on ",a6,1x,a4)
 113  format(/"** There was an error on group ",i4," processing task: ",
     &     a20," at: ",a8," on ",a6,1x,a4)

      return
      end


      subroutine write_summary(njob,filelog,rtncode_log,totaltime)
      implicit none
cINCLUDE(common/iofile)
      integer maxjobs,secs,njob,ioerr,iwrt,i,access
      integer nsuccess,nfail,ncrash
      parameter(maxjobs=5000,iwrt=9)
      character *132 filelog(maxjobs)
      integer rtncode_log(maxjobs)
      character *8 zyear,zdate,ztime,zaccno,zanam
      REAL totaltime
      logical endoffile

      save access
      data access/0/

      if (access .eq. 0) then
c     First time through just print the header
         open(unit=iwrt,file='taskfarm.summary',status='new'
     +        ,iostat=ioerr)
         if (ioerr .ne. 0) then
            write(*,*)"Could not open summary for writing!"
            return
         endif
         
         call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)

         write(iwrt,001)
         write(iwrt,002)ztime,zdate(1:6),zyear(1:4)
         close(unit=iwrt)
         access=1
         return
      else
         open(unit=iwrt,file='taskfarm.summary',status='old'
     +        ,iostat=ioerr)
         if (ioerr .ne. 0) then
            write(*,*)"Could not re-open summary for writing!"
            return
         endif

c     Read to the end of the file then back one
         endoffile = .false.
         do while (.not. endoffile)
            read(unit=9,fmt=*,iostat=ioerr)
            if (ioerr .lt. 0) then
               endoffile = .true.
               backspace(iwrt)
            elseif (ioerr .gt. 0) then
               write(*,*)"** Taskfarm.summary EOF error!"
               return
            endif
         enddo

         call tidajt(zdate,ztime,zyear,zaccno,zanam,secs)
         write(iwrt,003)ztime,zdate(1:6),zyear(1:4)
         write(iwrt,004)totaltime

         nsuccess = 0
         nfail = 0
         ncrash = 0
         
         write(iwrt,005)
         
         do i=1,njob
            if (rtncode_log(i) == 0) then
               nsuccess=nsuccess+1
               write(iwrt,006)filelog(i)
            elseif (rtncode_log(i) == 1) then
               nfail=nfail+1
               write(iwrt,007)filelog(i)
            elseif (rtncode_log(i) == 2) then
               ncrash=ncrash+1
               write(iwrt,008)filelog(i)
            else
               write(iwrt,*)"Missed a result"
            endif
         enddo
         
         write(iwrt,009)njob,nsuccess,nfail,ncrash
         close(iwrt)

         return
      endif

c 001  format(/38x,23('*')/40x,'Taskfarming Summary',
c     &     //38x,'job started at')
c 002  format(/38x,a8,' on ',a6,1x,a4/38x,23('*')//)

 001  format(/38x,23('*')/40x,'Taskfarming Summary',
     &     //38x,23('*'))
 002  format(//5x,'Taskfarming started at : ',a8,' on ',a6,1x,a4)

 003  format(5x,'Taskfarming finished at: ',a8,' on ',a6,1x,a4)
 004  format(/5x,'Total runtime: ',f7.2,' seconds')

 005  format(/5x,'Job Name',32x,5x,'Result'/5x,8('-'),37x,11('-'))
 006  format(5x,a40,5x,'Succeeded')
 007  format(5x,a40,5x,'Failed')
 008  format(5x,a40,5x,'**Crashed**')

 009  format(/5x,i4,' jobs processed'/5x,19('-')/
     +     5x,'Succeeded : ',i4/
     +     5x,'Failed    : ',i4/
     +     5x,'Crashed   : ',i4)
      end
