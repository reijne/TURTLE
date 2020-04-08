Module mpi_global_counter

  ! MPI2 Global counter module originally from A.G.Sunderland,
  ! slightly modified and extended by IJB. 
  ! 08/2008 Modified again by AGS: Deleted IJB's original modifications for IBM.
  ! 08/2008 Modified by HvD: Getting the MPI_ADDRESS_KIND arguments in the right
  !         places.
  ! 10/2008 Modified by HvD: Replaced nested locks (which did not work properly)
  !         with non-overlapping communications. This seems to work but
  !         certainly will not scale well on large processor counts.

  Use newscf_modules ! This module is currently used to access the 
                     ! MPI_COMM_GAMESS variable for the communicator.
                     ! A cleaner approach would be to add a subroutine to set
                     ! the communicator (e.g. set_communicator).
                     ! If there is a need for multiple counters (for example
                     ! across different subsets of processors) then the
                     ! communicator, task_ctr, win, myrank and com_size could
                     ! all be wrapped up in a user defined data type and passed
                     ! into the appropriate routines. For use within GAMESS-UK
                     ! these modifications are currently not necessary but 
                     ! a generally useful module could be derived from these
                     ! routines in that way. (HvD, 08/2008)

  Use mpi            ! Should be the prefered way to include fortran interfaces
                     ! not supported widely at present though. So if this fails
                     ! comment out the use statement and uncomment the mpif.h
                     ! include statement.

  Implicit None

! Include 'mpif.h'

  Public :: set_tasks, get_task, end_tasks, reset_tasks

  Private

  Integer, Save, Allocatable :: task_ctr(:)  ! task counter
  Integer, Save, Allocatable :: buffer(:)    ! task counter local buffer
  Integer, Save              :: win          ! window handle
  Integer, Save              :: myrank       ! processor number
  Integer, Save              :: com_size     ! # processors
  Integer, Save              :: my_ctr       ! the number of tasks I have done
  Integer, Save              :: size_int     ! stride size

Contains

  Subroutine set_tasks

    ! set up shared memory area defining the current task - BLOCKING

    Integer                   :: err, info, size_int
    Integer(MPI_ADDRESS_KIND) :: size_addr, lb_addr

    ! find rank and number of processors:

    Call MPI_COMM_RANK (MPI_COMM_GAMESS, myrank, err)
    Call MPI_COMM_SIZE (MPI_COMM_GAMESS, com_size, err)

    Call MPI_TYPE_GET_EXTENT (MPI_INTEGER, lb_addr, size_addr, err)

    ! create local buffer of task counters

    Allocate ( buffer(0:com_size-1) )

    ! create shared memory window on processor 0:

    size_int = size_addr
    If (myrank /= 0) size_addr = 0

    If (myrank == 0 ) then
      Allocate ( task_ctr(0:com_size-1) )
    Else
      Allocate ( task_ctr(0:0) ) ! allocate 1 integer to create a valid address
    Endif

    Call MPI_WIN_CREATE (task_ctr, com_size*size_addr, size_int, &
         MPI_INFO_NULL, MPI_COMM_GAMESS, win, err)

    Call reset_tasks

  End Subroutine set_tasks

  Subroutine get_task (new_task)

    ! Get and update task counter

    ! This implementation was suggested by David Henty. It addresses problems
    ! with the ordering of communications during a single epoch in that they
    ! are not guaranteed to complete in the order they were issued. Attempts
    ! to use multiple epochs to guarantee the order led to problems with locks.
    ! This implementation avoids the issue by have a task counter for every
    ! process. During get_task all elements of rank lower than myrank and 
    ! higher than myrank are read, the element of rank myrank is written. Due
    ! to the non-overlapping nature of accessing the elements the order problem
    ! goes away. The price of this is having to communicate and add com_size
    ! integers instead of just one. (HvD, 01/10/2008)

    Integer, Intent(out)           :: new_task
    Integer                        :: err, assert, m1, p1, root, i
    Integer(KIND=MPI_ADDRESS_KIND) :: size_addr

    new_task = -777
    assert = 0
    m1 = -1
    p1 = 1
    root = 0

    Call MPI_WIN_LOCK (MPI_LOCK_EXCLUSIVE, root, assert, win, err)

    size_addr = 0 ! size_addr trick to pass correct integer type to MPI
    Call MPI_GET (buffer(0:myrank-1), myrank, MPI_INTEGER, &
         root, size_addr, myrank, MPI_INTEGER, win, err)

    size_addr = myrank+1
    Call MPI_GET (buffer(myrank+1:com_size-1), com_size-myrank-1, MPI_INTEGER, &
         root, size_addr, com_size-myrank-1, MPI_INTEGER, win, err)

    size_addr = myrank
    buffer(myrank) = my_ctr
    my_ctr = my_ctr + p1
    Call MPI_PUT(my_ctr, 1, MPI_INTEGER, &
         root, size_addr, 1, MPI_INTEGER, win, err)

    Call MPI_WIN_UNLOCK (root, win, err)

    new_task = 0
    Do i = 0, com_size-1
      new_task = new_task + buffer(i)
    Enddo

  End Subroutine get_task

  Subroutine reset_tasks

    ! Reset counter - BLOCKING

    Integer :: err

    Call MPI_BARRIER( MPI_COMM_GAMESS, err )
    my_ctr = 0
    If( myrank == 0 ) Then
       task_ctr = 0
    End If
    Call MPI_BARRIER( MPI_COMM_GAMESS, err )

  End Subroutine reset_tasks

  Subroutine end_tasks

    ! Kill MPI windows - BLOCKING

    Integer   :: err

    Call MPI_BARRIER (MPI_COMM_GAMESS, err)

    Call MPI_WIN_FREE (win , err)

    Deallocate( buffer)
    Deallocate( task_ctr )

  End Subroutine end_tasks

End Module mpi_global_counter
