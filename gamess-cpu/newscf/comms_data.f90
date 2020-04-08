Module comms_data

  Use allocation

  Implicit None

  Private

  Integer, Public :: comms_data_all_comm
  Integer, Public :: comms_data_uhf_k_comm

  Integer, Public :: comms_data_all_context
  Integer, Public :: comms_data_uhf_k_context

  Integer, Public :: comms_data_all_row_comm
  Integer, Public :: comms_data_all_col_comm

  Integer, Public :: comms_data_uhf_k_row_comm
  Integer, Public :: comms_data_uhf_k_col_comm

  Integer, Public :: comms_data_all_nproc
  Integer, Public :: comms_data_all_me

  Integer, Public :: comms_data_uhf_k_nproc
  Integer, Public :: comms_data_uhf_k_me

  Integer, Public :: comms_data_all_nrow
  Integer, Public :: comms_data_all_ncol
  Integer, Public :: comms_data_all_my_row
  Integer, Public :: comms_data_all_my_col

  Integer, Public :: comms_data_uhf_k_nrow
  Integer, Public :: comms_data_uhf_k_ncol
  Integer, Public :: comms_data_uhf_k_my_row
  Integer, Public :: comms_data_uhf_k_my_col

  Integer, Public :: comms_data_uhf_k_repl_comm
  Integer, Public :: comms_data_uhf_k_repl_nproc
  Integer, Public :: comms_data_uhf_k_repl_me

  ! Jens
  Logical, Public :: taskfarm_alreadywarned

  Logical, Dimension( : ), Allocatable, Public :: comms_data_i_own

  Logical, Public :: comms_data_uhf_k_is_all

  Integer, Parameter, Public :: comms_data_block_default = 96
  Integer, Parameter, Public :: comms_data_not_valid     = -999

  Public :: comms_data_setup
  Public :: comms_data_free

Contains

  Subroutine comms_data_setup( base_comm, base_context, nspin, nk )

    Use newscf_modules, Only : is_farm, iwr

    Include 'mpif.h'

    Integer, Intent( In ) :: base_comm
    Integer, Intent( In ) :: base_context
    Integer, Intent( In ) :: nspin
    Integer, Intent( In ) :: nk

    Integer, Dimension( :, : ), Allocatable :: proc_map

    Integer :: dummy
    Integer :: tmp
    Integer :: nproc_per_job
    Integer :: error
    Integer :: this_proc, base_proc
    Integer :: repl_flag
    Integer :: i, j, job

    ! Set base communicator for matrix parallel ops
    Call mpi_comm_dup( base_comm, comms_data_all_comm, error )

    ! Set up the 'all' stuff
    Call mpi_comm_size( comms_data_all_comm, comms_data_all_nproc, error )
    Call mpi_comm_rank( comms_data_all_comm, comms_data_all_me   , error )

    ! Factor the grid as square as possible
    Call factor( comms_data_all_nproc, comms_data_all_nrow, comms_data_all_ncol )

    If( .Not. is_farm ) Then
        ! Set up the process map for BLACS
       Allocate( proc_map( 0:comms_data_all_nrow - 1, 0:comms_data_all_ncol - 1 ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in COMMS_DATA_SETUP', 'abort' )
       End If
       Call alloc_add( ALLOC_INTEGER, Size( proc_map ) )
       this_proc = 0
       Do i = 0, comms_data_all_ncol - 1
          Do j = 0, comms_data_all_nrow - 1
             proc_map( j, i ) = this_proc
             this_proc = this_proc + 1
          End Do
       End Do
       
       ! Set up the blacs proc grid and related stuff
       tmp = base_context
       Call blacs_gridmap( tmp, proc_map, comms_data_all_nrow, comms_data_all_nrow, &
            comms_data_all_ncol )
       comms_data_all_context = tmp
       Call alloc_free( ALLOC_INTEGER, Size( proc_map ) )
       Deallocate( proc_map )
    Else
       comms_data_all_context = base_context
    End If

    Call blacs_gridinfo( comms_data_all_context, comms_data_all_nrow  , &
                                                 comms_data_all_ncol  , &
                                                 comms_data_all_my_row, &
                                                 comms_data_all_my_col )

    ! Now set up the row and column communicators for the 'all' communicator
    Call mpi_comm_split( comms_data_all_comm, comms_data_all_my_col, &
                         comms_data_all_my_row, comms_data_all_col_comm, error )
    Call mpi_comm_split( comms_data_all_comm, comms_data_all_my_row, &
                         comms_data_all_my_col, comms_data_all_row_comm, error )

    ! Make life easy - if number of procs is not a multiple
    ! of the number of spin * k DO NOT use split comms
    If( nspin * nk == 1 ) Then
       comms_data_uhf_k_is_all = .True.
    Else If( comms_data_all_nproc == 1 ) Then
       comms_data_uhf_k_is_all = .True.
    Else If( is_farm ) Then
       If( comms_data_all_me == 0 ) Then
          If( taskfarm_alreadywarned .eqv. .True. ) Then
              continue
          Else
             Write( iwr, * ) 'For taskfarms split communicators are NOT yet implemented'
          End If
       End If
       taskfarm_alreadywarned = .True.
       comms_data_uhf_k_is_all = .True.
    Else If( Mod( comms_data_all_nproc, nspin * nk ) /= 0 ) Then
       If( comms_data_all_me == 0 ) Then
          Write( iwr, * ) 'Number of procs is not a multiple of the number of jobs'
          Write( iwr, * ) 'Split communicators NOT being used'
       End If
       comms_data_uhf_k_is_all = .True.
    Else
       comms_data_uhf_k_is_all = .False.
    End If

    Allocate( comms_data_i_own( 1:nspin * nk ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in COMMS_DATA_SETUP', 'abort' )
    End If
    Call alloc_add( ALLOC_LOGICAL, Size( comms_data_i_own ) )

    If( comms_data_uhf_k_is_all ) Then

       ! Easy case - not using split communicators
       comms_data_uhf_k_comm      = comms_data_all_comm
       comms_data_uhf_k_context   = comms_data_all_context
       comms_data_uhf_k_nproc     = comms_data_all_nproc
       comms_data_uhf_k_me        = comms_data_all_me
       comms_data_uhf_k_nrow      = comms_data_all_nrow
       comms_data_uhf_k_ncol      = comms_data_all_ncol
       comms_data_uhf_k_my_row    = comms_data_all_my_row
       comms_data_uhf_k_my_col    = comms_data_all_my_col
       comms_data_uhf_k_col_comm  = comms_data_all_col_comm
       comms_data_uhf_k_row_comm  = comms_data_all_row_comm

       comms_data_i_own = .True.

    Else

       ! Tricky case - using split communicators

       ! Split the MPI base communicator
       ! Put first job on the first NPROC_PER_JOB procs etc.
       nproc_per_job = comms_data_all_nproc / ( nspin * nk )
       Call mpi_comm_split( comms_data_all_comm, comms_data_all_me / nproc_per_job, &
                            comms_data_all_me, comms_data_uhf_k_comm, error )

       Call mpi_comm_size( comms_data_uhf_k_comm, comms_data_uhf_k_nproc, error )
       Call mpi_comm_rank( comms_data_uhf_k_comm, comms_data_uhf_k_me   , error )

       
       ! Again factor
       Call factor( comms_data_uhf_k_nproc, comms_data_uhf_k_nrow, comms_data_uhf_k_ncol )

       If( .Not. is_farm ) Then
          ! BLACS is really, really, really, really icky for splitting contexts !
          Allocate( proc_map( 0:comms_data_uhf_k_nrow - 1, 0:comms_data_uhf_k_ncol - 1 ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in COMMS_DATA_SETUP', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( proc_map ) )
          this_proc = 0
          comms_data_i_own = .False.
          Do job = 1, nspin * nk
             base_proc = this_proc
             Do i = 0, comms_data_uhf_k_ncol - 1
                Do j = 0, comms_data_uhf_k_nrow - 1
                   proc_map( j, i ) = this_proc
                   this_proc = this_proc + 1
                End Do
             End Do
             tmp = base_context
             Call blacs_gridmap( tmp, proc_map, comms_data_uhf_k_nrow, comms_data_uhf_k_nrow, &
                  comms_data_uhf_k_ncol )
             If( comms_data_all_me >= base_proc .And. comms_data_all_me < this_proc ) Then
                comms_data_uhf_k_context = tmp
                comms_data_i_own( job ) = .True.
             Else
                comms_data_i_own( job ) = .False.
             End If
             base_proc = this_proc
          End Do
          Call alloc_free( ALLOC_INTEGER, Size( proc_map ) )
          Deallocate( proc_map )
       End If

       Call blacs_gridinfo( comms_data_uhf_k_context, comms_data_uhf_k_nrow  , &
                                                      comms_data_uhf_k_ncol  , &
                                                      comms_data_uhf_k_my_row, &
                                                      comms_data_uhf_k_my_col )

       ! Now set up the row and column communicators for the 'uhf_k' communicator
       Call mpi_comm_split( comms_data_uhf_k_comm, comms_data_uhf_k_my_col, &
                            comms_data_uhf_k_my_row, comms_data_uhf_k_col_comm, error )
       Call mpi_comm_split( comms_data_uhf_k_comm, comms_data_uhf_k_my_row, &
                            comms_data_uhf_k_my_col, comms_data_uhf_k_row_comm, error )

       ! Finally set up a communicator that can ne used to replicate UHF_K
       ! arrays of matrices
       Call mpi_comm_split( comms_data_all_comm, comms_data_uhf_k_me, &
                            comms_data_all_me, comms_data_uhf_k_repl_comm, error )
       Call mpi_comm_size( comms_data_uhf_k_repl_comm, comms_data_uhf_k_repl_nproc, error )
       Call mpi_comm_rank( comms_data_uhf_k_repl_comm, comms_data_uhf_k_repl_me   , error )

    End If

  Contains
    
    Subroutine factor( n, nr, nc )

      ! Factor n into two factors which as as close to each other
      ! as possible

      Integer, Intent( In    ) :: n
      Integer, Intent(   Out ) :: nr
      Integer, Intent(   Out ) :: nc

      Integer :: top
      Integer :: i

      top = Nint( Sqrt( Real( n ) ) )

      ! Do it this way so that in the worst case we can always ensure
      ! nc = n and nr = 1 as factors
      Do i = top, 1, -1
         nr = i
         nc = n / nr
         If( nc * nr == n ) Then
            Exit
         End If
      End Do

    End Subroutine factor

  End Subroutine comms_data_setup

  Subroutine comms_data_free

    Use newscf_modules, Only : is_farm

    Integer :: error

    If( .Not. comms_data_uhf_k_is_all ) Then
       Call mpi_comm_free ( comms_data_uhf_k_repl_comm, error )
       Call mpi_comm_free ( comms_data_uhf_k_row_comm , error )
       Call mpi_comm_free ( comms_data_uhf_k_col_comm , error )
       If( .Not. is_farm ) Then
          Call blacs_gridexit( comms_data_uhf_k_context          )
       End If
       Call mpi_comm_free ( comms_data_uhf_k_comm     , error )
    End If

    Call alloc_free( ALLOC_LOGICAL, Size( comms_data_i_own ) )
    Deallocate( comms_data_i_own )

    Call mpi_comm_free( comms_data_all_row_comm, error )
    Call mpi_comm_free( comms_data_all_col_comm, error )

    If( .Not. is_farm ) Then
       Call blacs_gridexit( comms_data_all_context     )
    End If
    Call mpi_comm_free ( comms_data_all_comm, error )

  End Subroutine comms_data_free

End Module comms_data
