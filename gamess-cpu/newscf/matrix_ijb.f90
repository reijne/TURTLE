Module distributed_matrices

  Use newscf_numbers
  Use allocation

  Implicit None

  Public :: matrix
  Public :: matrix_create
  Public :: matrix_destroy
  Public :: matrix_add
  Public :: matrix_add_diagonal
  Public :: matrix_assign
  Public :: matrix_set
  Public :: matrix_set_element
  Public :: matrix_set_random
  Public :: matrix_set_from_triangle
  Public :: matrix_set_from_square
  Public :: matrix_get
  Public :: matrix_get_element
  Public :: matrix_get_to_triangle
  Public :: matrix_get_to_square
  Public :: matrix_set_column
  Public :: matrix_set_row
  Public :: matrix_get_column
  Public :: matrix_get_row
  Public :: matrix_length
  Public :: matrix_dot_product
  Public :: matrix_daxpy
  Public :: matrix_copy
  Public :: matrix_copy_cols
  Public :: matrix_scale
  Public :: matrix_dimension
  Public :: matrix_inquire
  Public :: matrix_multiply
  Public :: matrix_dgemm
  Public :: matrix_mult_by_diag
  Public :: matrix_mult2
  Public :: matrix_mult2t
  Public :: matrix_combine
  Public :: matrix_abs
  Public :: matrix_absmax
  Public :: matrix_diagonalise
  Public :: matrix_invert
  Public :: matrix_absmaxloc_col
  Public :: matrix_check_abs_ov
  Public :: matrix_write
  Public :: matrix_print
  Public :: matrix_print_titled
  Public :: matrix_read_rectangle
  Public :: matrix_read_triangle
  Public :: matrix_zero_patch
  Public :: matrix_set_patch
  Public :: matrix_orthogonalize
  Public :: matrix_comms_init
  Public :: matrix_comms_finalize
  Public :: matrix_write_rectangle
  Public :: matrix_write_triangle
  Public :: matrix_tdown_setup
  Public :: matrix_tdown_free
  Public :: matrix_tdown
  Public :: matrix_stretch
  Public :: matrix_direct_product
  Public :: matrix_get_col_biggest

  Integer, Parameter, Public :: MATRIX_REAL    = 1
  Integer, Parameter, Public :: MATRIX_COMPLEX = 2

  Integer, Parameter, Public :: MATRIX_DISTRIB_REPL  = 1
  Integer, Parameter, Public :: MATRIX_DISTRIB_ALL   = 2
  Integer, Parameter, Public :: MATRIX_DISTRIB_UHF_K = 3
  Integer, Parameter, Public :: MATRIX_DISTRIB_DEFAULT = MATRIX_DISTRIB_REPL

  Private

  Type matrix_distribution
     Private
     Integer                   :: type
     Integer                   :: communicator
     Integer                   :: row_comm
     Integer                   :: col_comm
     Integer                   :: n_proc
     Integer                   :: me
     Integer                   :: context
     Integer                   :: n_proc_row, n_proc_col
     Integer                   :: my_proc_row, my_proc_col
     Integer, Dimension( 1:9 ) :: descriptor
  End Type matrix_distribution

  Type matrix_indexing
     Private
     Integer                          :: n
     Integer                          :: blocking_factor
     Integer                          :: local_n
     Integer                          :: local_ld
     Integer, Dimension( : ), Pointer :: global_to_local
     Integer, Dimension( : ), Pointer :: local_to_global
     Integer, Dimension( : ), Pointer :: who_owns
  End Type matrix_indexing

  Type matrix
     Private
     Type( matrix_distribution ), Pointer      :: distribution
     Type( matrix_indexing     ), Pointer      :: row
     Type( matrix_indexing     ), Pointer      :: col
     Logical                                   :: i_own
     Logical                                   :: inited = .False.
     Logical                                   :: wierd
     Integer                                   :: genus
     Integer                                   :: n, m
     Integer                                   :: local_n , local_m
     Integer                                   :: local_ld, local_sd
     Character( Len = 32 )                     :: name
     Real   ( wp ), Dimension( :, : ), Pointer :: real_data    => Null()
     Complex( wp ), Dimension( :, : ), Pointer :: complex_data => Null()
  End Type matrix

  Integer, Parameter :: MATRIX_UNKNOWN_GENUS   = -1
  Integer, Parameter :: MATRIX_NOT_OWN         = -1

  Type( matrix_distribution ), Target :: distrib_repl
  Type( matrix_distribution ), Target :: distrib_all
  Type( matrix_distribution ), Target :: distrib_uhf_k

  Type( matrix_indexing ), Target :: index_repl_row
  Type( matrix_indexing ), Target :: index_all_row
  Type( matrix_indexing ), Target :: index_uhf_k_row

  Type( matrix_indexing ), Target :: index_repl_col
  Type( matrix_indexing ), Target :: index_all_col
  Type( matrix_indexing ), Target :: index_uhf_k_col

  Integer, Parameter :: MATRIX_SHIFT_COL = 1
  Integer, Parameter :: MATRIX_SHIFT_ROW = 2

  Type tdown_ops
     Integer, Dimension( : ), Pointer :: j_vec_add
     Integer, Dimension( : ), Pointer :: j_tran    
     Integer, Dimension( : ), Pointer :: j_vec_base
  End Type tdown_ops
  
  Type( tdown_ops ), Dimension( : ), Allocatable :: remote_ops
  
  Type( tdown_ops ) :: all_ops

  ! General workspace buffer, used in such things as replications
  ! and I/O
  Integer, Parameter :: n_repl_buff = ( 1024 * 1024 ) / 8
  Real( wp ), Dimension( 1:n_repl_buff ) :: buff
  Integer, Parameter :: io_block = 511
  Integer, Parameter :: n_io_buff = ( n_repl_buff / io_block ) * io_block
  ! Unfortunately need a seperate buffer for the writes
  Real( wp ), Dimension( 1:n_repl_buff ) :: buff_write

  Interface matrix_create
     Module Procedure matrix_create_one
     Module Procedure matrix_create_multi
     Module Procedure matrix_create_template_one
     Module Procedure matrix_create_template_multi
  End Interface

  Interface matrix_destroy
     Module Procedure matrix_destroy_one
     Module Procedure matrix_destroy_multi
  End Interface

  Interface matrix_copy
     Module Procedure matrix_copy_one
     Module Procedure matrix_copy_multi
  End Interface

  Interface matrix_copy_cols
     Module Procedure matrix_copy_cols_one
     Module Procedure matrix_copy_cols_multi
  End Interface

  Interface matrix_assign
     Module Procedure matrix_assign_one
     Module Procedure matrix_assign_multi
  End Interface

  Interface matrix_set
     Module Procedure matrix_assign_one
     Module Procedure matrix_assign_multi
     Module Procedure matrix_set_element
     Module Procedure matrix_set_1d
     Module Procedure matrix_set_random
  End Interface

  Interface matrix_get
     Module Procedure matrix_get_element
     Module Procedure matrix_get_1d
  End Interface

  Interface matrix_write
     Module Procedure matrix_print
     Module Procedure matrix_print_titled
     Module Procedure matrix_print_to_file
  End Interface

  Interface matrix_add_diagonal
     Module Procedure matrix_add_diagonal_one
     Module Procedure matrix_add_diagonal_multi
  End Interface

  Interface matrix_daxpy
     Module Procedure matrix_daxpy_one
     Module Procedure matrix_daxpy_multi
  End Interface

  Interface matrix_dot_product
     Module Procedure matrix_dot_product_one
     Module Procedure matrix_dot_product_multi
  End Interface

  Interface matrix_length
     Module Procedure matrix_length_one
     Module Procedure matrix_length_multi
  End Interface

  Interface matrix_abs
     Module Procedure matrix_abs_one
     Module Procedure matrix_abs_multi
  End Interface

  Interface matrix_absmax
     Module Procedure matrix_absmax_one
     Module Procedure matrix_absmax_multi
  End Interface

  Interface matrix_diagonalise
     Module Procedure matrix_diagonalise_one
     Module Procedure matrix_diagonalise_multi
  End Interface

  Interface matrix_invert
     Module Procedure matrix_invert_one
     Module Procedure matrix_invert_multi
  End Interface

  Interface matrix_multiply
     Module Procedure matrix_multiply_one
     Module Procedure matrix_multiply_multi
  End Interface

  Interface matrix_dgemm
     Module Procedure matrix_dgemm_one
     Module Procedure matrix_dgemm_multi
  End Interface

  Interface matrix_mult_by_diag
     Module Procedure matrix_mult_by_diag_one
     Module Procedure matrix_mult_by_diag_multi
  End Interface

  Interface matrix_mult2
     Module Procedure matrix_mult2_one
     Module Procedure matrix_mult2_multi
  End Interface

  Interface matrix_mult2t
     Module Procedure matrix_mult2t_one
     Module Procedure matrix_mult2t_multi
  End Interface

  Interface matrix_zero_patch
     Module Procedure matrix_zero_patch_one
     Module Procedure matrix_zero_patch_multi
  End Interface

  Interface matrix_set_patch
     Module Procedure matrix_set_patch_one
     Module Procedure matrix_set_patch_multi
  End Interface

  Interface matrix_absmaxloc_col
     Module Procedure matrix_absmaxloc_col_one
     Module Procedure matrix_absmaxloc_col_multi
  End Interface

  Interface matrix_check_abs_ov
     Module Procedure matrix_check_abs_ov_one
     Module Procedure matrix_check_abs_ov_multi
  End Interface

  Interface matrix_orthogonalize
     Module Procedure matrix_orthogonalize_one
     Module Procedure matrix_orthogonalize_multi
  End Interface

  Interface matrix_stretch
     Module Procedure matrix_stretch_one
     Module Procedure matrix_stretch_multi
  End Interface

  Interface matrix_direct_product
     Module Procedure matrix_direct_product_one
     Module Procedure matrix_direct_product_multi
  End Interface

  Interface matrix_get_col_biggest
     Module Procedure matrix_get_col_biggest_one
     Module Procedure matrix_get_col_biggest_multi
  End Interface

  Interface prmat
     Module Procedure prmat_real
     Module Procedure prmat_complex
  End Interface

  Interface Operator( == )
     Module Procedure matrix_compare_distrib
  End Interface

Contains

  Subroutine matrix_create_one( n, m, a, name, genus, distrib, i_own )

    Use comms_data

    Include 'mpif.h'

    Integer              , Intent( In )           :: n
    Integer              , Intent( In )           :: m
    Type( matrix )                                :: a
    Integer              , Intent( In )           :: genus
    Character( Len = * ) , Intent( In )           :: name
    Integer              , Intent( In ), Optional :: distrib
    Logical              , Intent( In ), Optional :: i_own

    Integer, External :: numroc

    Logical, External :: scalapack_avail

    Integer :: this_distrib
    Integer :: local_n , local_m
    Integer :: local_ld, local_sd
    Integer :: block
    Integer :: tmp
    Integer :: error

    Logical :: this_i_own


    If( Present( distrib ) ) Then
       this_distrib = distrib
    Else
       this_distrib = MATRIX_DISTRIB_DEFAULT
    End If

    If( Present( i_own ) ) Then
       this_i_own = i_own
    Else
       this_i_own = .True.
    End If

    If( comms_data_all_nproc == 1 .or. .not. scalapack_avail()) Then
       this_distrib =  MATRIX_DISTRIB_REPL
    End If

    Select Case( this_distrib )

    Case( MATRIX_DISTRIB_REPL )
       a%distribution => distrib_repl
       a%row          => index_repl_row
       a%col          => index_repl_col
       this_i_own     =  .True.
       
    Case( MATRIX_DISTRIB_ALL )
       a%distribution => distrib_all
       a%row          => index_all_row
       a%col          => index_all_col
       this_i_own     =  .True.
       
    Case( MATRIX_DISTRIB_UHF_K )
       a%distribution => distrib_uhf_k
       a%row          => index_uhf_k_row
       a%col          => index_uhf_k_col
       
    Case Default
       Call dlc_error( 'Unknown matrix distribution in MATRIX_CREATE', 'abort' )
       
    End Select

    a%i_own    = this_i_own
    ! Check if matrix distrib corresponds to one of the defaults
    If( a%row%n == n .And. a%col%n == m ) Then

       a%wierd    = .False.

    Else

       ! It doesn't have to set by hand for this case
       a%wierd    = .True.

       Allocate( a%distribution, Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, 1 )
       Allocate( a%row         , Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, 1 )
       Allocate( a%col         , Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, 1 )

       ! First set up the distribution data. Most of this can be inherited from
       ! the default.
       ! The one exception is the descriptor for which the leading local
       ! dim will have to change. Might as well get all the local dims here.

       Select Case( this_distrib )

       Case( MATRIX_DISTRIB_REPL )
          a%distribution = distrib_repl
          block = index_repl_row%blocking_factor

       Case( MATRIX_DISTRIB_ALL )
          a%distribution = distrib_all
          block = index_all_row%blocking_factor

       Case( MATRIX_DISTRIB_UHF_K )
          a%distribution = distrib_uhf_k
          block = index_uhf_k_row%blocking_factor

       Case Default
          Call dlc_error( 'Unknown matrix distribution in MATRIX_CREATE', 'abort' )

       End Select

       If( this_distrib /= MATRIX_DISTRIB_REPL ) Then

          local_n = numroc( n, block, a%distribution%my_proc_row, 0, &
               a%distribution%n_proc_row )
          local_m = numroc( m, block, a%distribution%my_proc_col, 0, &
               a%distribution%n_proc_col )

          tmp = local_n
          Call mpi_allreduce( tmp, local_ld, 1, mpi_integer, mpi_max, &
               a%distribution%communicator, error )
          tmp = local_m
          Call mpi_allreduce( tmp, local_sd, 1, mpi_integer, mpi_max, &
               a%distribution%communicator, error )
          
          ! Know sizes, can now set up descriptor
          Call descinit( a%distribution%descriptor, n, m, &
               block, block, 0, 0, &
               a%distribution%context, local_ld, error )

       Else
          
          local_n  = n
          local_m  = m
          local_ld = n
          local_sd = m

          Call descinit( a%distribution%descriptor, n, m, &
               block, block, 0, 0, &
               a%distribution%context, local_ld, error )

       End If

       If( a%i_own ) Then

          a%row%n               = n
          a%row%blocking_factor = block
          a%row%local_n         = local_n
          a%row%local_ld        = local_ld
          Allocate( a%row%global_to_local( 1:n ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%row%global_to_local ) )
          Call set_g_to_l( block, a%distribution%n_proc_row, &
               a%distribution%my_proc_row, a%row%global_to_local )
          Allocate( a%row%local_to_global( 1:local_n ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%row%local_to_global ) )
          Call set_l_to_g( n, block, a%distribution%n_proc_row, &
               a%distribution%my_proc_row, a%row%local_to_global )
          Allocate( a%row%who_owns( 1:n ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%row%who_owns ) )
          Call set_who_owns( n, block, a%distribution%n_proc_row, a%row%who_owns )

          a%col%n               = m
          a%col%blocking_factor = block
          a%col%local_n         = local_m
          a%col%local_ld        = local_sd
          Allocate( a%col%global_to_local( 1:m ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%col%global_to_local ) )
          Call set_g_to_l( block, a%distribution%n_proc_col, &
               a%distribution%my_proc_col, a%col%global_to_local )
          Allocate( a%col%local_to_global( 1:local_m ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%col%local_to_global ) )
          Call set_l_to_g( m, block, a%distribution%n_proc_col, &
               a%distribution%my_proc_col, a%col%local_to_global )
          Allocate( a%col%who_owns( 1:m ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( a%col%who_owns ) )
          Call set_who_owns( m, block, a%distribution%n_proc_col, a%col%who_owns )

       End If

    End If

    a%n        = n
    a%m        = m
    If( .Not. a%wierd .Or. a%i_own ) Then
       a%local_n  = a%row%local_n
       a%local_m  = a%col%local_n
       a%local_ld = a%row%local_ld
       a%local_sd = a%col%local_ld
    End If
    a%name     = name
    a%genus    = genus

    If( a%i_own ) Then

       Select Case( a%genus )
          
          Case( MATRIX_REAL )
             Allocate( a%real_data( 1:a%local_ld, 1:a%local_sd ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( a%real_data ) )

          Case( MATRIX_COMPLEX )
             Allocate( a%complex_data( 1:a%local_ld, 1:a%local_sd ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Matrix allocation error in MATRIX_CREATE', 'abort' )
             End If
             Call alloc_add( ALLOC_COMPLEX, Size( a%complex_data ) )

          Case Default
             Call dlc_error( 'Unknown matrix type in MATRIX_CREATE', 'abort' )

       End Select

    End If

    a%inited = .True.

  End Subroutine matrix_create_one

  Subroutine matrix_create_multi( n, m, a, name, genus, distrib )

    Use comms_data

    Integer              ,                 Intent( In )           :: n
    Integer              ,                 Intent( In )           :: m
    Type( matrix )       , Dimension( : )                         :: a
    Integer              , Dimension( : ), Intent( In )           :: genus
    Character( Len = * ) , Dimension( : ), Intent( In )           :: name
    Integer                              , Intent( In ), Optional :: distrib

    Integer :: this_distrib
    Integer :: i

    If( Present( distrib ) ) Then
       this_distrib = distrib
    Else
       this_distrib = MATRIX_DISTRIB_DEFAULT
    End If

    Do i = 1, Size( a )
       If( this_distrib == MATRIX_DISTRIB_UHF_K ) Then
          Call matrix_create( n, m, a( i ), name( i ), genus( i ), this_distrib, &
               comms_data_i_own( i ) )
       Else
          Call matrix_create( n, m, a( i ), name( i ), genus( i ), this_distrib, &
               .True. )
       End If
    End Do

  End Subroutine matrix_create_multi

  Subroutine matrix_create_template_one( a, b, name )

    Type( matrix )                               :: a
    Type( matrix )                               :: b
    Character( Len = * ), Intent( In ), Optional :: name

    Character( Len = 32 ) :: this_name

    If( .Not. a%inited ) Then
       Call dlc_error( 'Matrix not initialized in MATRIX_CREATE_TEMPLATE', 'abort' )
    End If

    If( Present( name ) ) Then
       this_name = name
    Else
       this_name = a%name
    End If

    Call matrix_create( a%n, a%m, b, this_name, a%genus, a%distribution%type )

  End Subroutine matrix_create_template_one

  Subroutine matrix_create_template_multi( a, b, name )

    Type( matrix )      , Dimension( : )                         :: a
    Type( matrix )      , Dimension( : )                         :: b
    Character( Len = * ), Dimension( : ), Intent( In ), Optional :: name

    Character( Len = 32 ), Dimension( : ), Allocatable :: this_name

    Integer :: error
    Integer :: i

    Allocate( this_name( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_CREATE_TEMPLATE', 'abort' )
    End If

    Do i = 1, Size( a )
       If( .Not. a( i )%inited ) Then
          Call dlc_error( 'Matrix not initialized in MATRIX_CREATE_TEMPLATE', 'abort' )
       End If
    End Do

    If( Present( name ) ) Then
       this_name = name
    Else
       this_name = a%name
    End If

    Call matrix_create( a( 1 )%n, a( 1 )%m, b, this_name, a%genus, a( 1 )%distribution%type )

    Deallocate( this_name )

  End Subroutine matrix_create_template_multi

  Subroutine matrix_destroy_one( a )

    Type( matrix ) :: a

    Character( Len = 14 ), Parameter :: routine = 'MATRIX_DESTROY' 

    If( a%distribution%type /= MATRIX_DISTRIB_UHF_K ) Then
       Call matrix_check_valid( routine, a )
    End If

    If( a%i_own ) Then

       Select Case( a%genus )

          Case( MATRIX_REAL )
             Call alloc_free( ALLOC_REAL, Size( a%real_data ) )
             Deallocate( a%real_data )
             Nullify( a%real_data )

          Case( MATRIX_COMPLEX )
             Call alloc_free( ALLOC_COMPLEX, Size( a%complex_data ) )
             Deallocate( a%complex_data )
             Nullify( a%complex_data )

       End Select

    End If

    a%genus    = MATRIX_UNKNOWN_GENUS
    a%n        = 0
    a%m        = 0
    a%local_n  = 0
    a%local_m  = 0
    a%local_ld = 0
    a%local_sd = 0
    a%name     = ''

    If( .Not. a%wierd ) Then

       Nullify( a%distribution )
       Nullify( a%row          )
       Nullify( a%col          )

    Else

       If( a%i_own ) Then

          Call alloc_free( ALLOC_INTEGER, Size( a%col%who_owns ) )
          Deallocate( a%col%who_owns )
          Call alloc_free( ALLOC_INTEGER, Size( a%col%local_to_global ) )
          Deallocate( a%col%local_to_global )
          Call alloc_free( ALLOC_INTEGER, Size( a%col%global_to_local ) )
          Deallocate( a%col%global_to_local )
          
          Call alloc_free( ALLOC_INTEGER, Size( a%row%who_owns ) )
          Deallocate( a%row%who_owns )
          Call alloc_free( ALLOC_INTEGER, Size( a%row%local_to_global ) )
          Deallocate( a%row%local_to_global )
          Call alloc_free( ALLOC_INTEGER, Size( a%row%global_to_local ) )
          Deallocate( a%row%global_to_local )

       End If

       Call alloc_free( ALLOC_DERIVED, 1 )
       Deallocate( a%col )
       Call alloc_free( ALLOC_DERIVED, 1 )
       Deallocate( a%row )
       Call alloc_free( ALLOC_DERIVED, 1 )
       Deallocate( a%distribution )

    End If

    a%i_own    = .False.
    a%inited = .False.

  End Subroutine matrix_destroy_one

  Subroutine matrix_destroy_multi( a )

    Type( matrix ), Dimension( : ) :: a

    Integer :: i

    Do i = Size( a ), 1, -1
       Call matrix_destroy( a( i ) )
    End Do

  End Subroutine matrix_destroy_multi

  Subroutine matrix_add( c, a, b )

    Type( matrix ) :: c
    Type( matrix ) :: a
    Type( matrix ) :: b

    Character( Len = 10 ), Parameter :: routine = 'MATRIX_ADD' 

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )
    Call matrix_check_valid( routine, c )

    Call matrix_check_conform( routine, a, b )
    Call matrix_check_conform( routine, a, c )

    Select Case( c%genus )

       Case( MATRIX_REAL )
          c%real_data( 1:c%local_n, 1:c%local_m ) = a%real_data( 1:a%local_n, 1:a%local_m ) + &
                                                    b%real_data( 1:b%local_n, 1:b%local_m )

       Case( MATRIX_COMPLEX )
          c%complex_data( 1:c%local_n, 1:c%local_m ) = a%complex_data( 1:a%local_n, 1:a%local_m ) + &
                                                       b%complex_data( 1:b%local_n, 1:b%local_m )

    End Select

  End Subroutine matrix_add

  Subroutine matrix_add_diagonal_one( a, b )

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: b

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_ADD_DIAGONAL' 

    Integer :: row, col
    Integer :: i

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = 1, a%n
             row = a%row%global_to_local( i )
             col = a%col%global_to_local( i )
             If( row /= MATRIX_NOT_OWN .And. col /= MATRIX_NOT_OWN ) Then
                a%real_data( row, col ) = a%real_data( row, col ) + b( i )
             End If
          End Do

       Case( MATRIX_COMPLEX )
          Do i = 1, a%n
             row = a%row%global_to_local( i )
             col = a%col%global_to_local( i )
             If( row /= MATRIX_NOT_OWN .And. col /= MATRIX_NOT_OWN ) Then
                a%complex_data( row, col ) = a%complex_data( row, col ) + b( i )
             End If
          End Do

    End Select

  End Subroutine matrix_add_diagonal_one

  Subroutine matrix_add_diagonal_multi( a, b )

    Type( matrix ), Dimension( :    )               :: a
    Real( wp )    , Dimension( :, : ), Intent( In ) :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_add_diagonal( a( i ), b( :, i ) )
       End If
    End Do

  End Subroutine matrix_add_diagonal_multi

  Subroutine matrix_assign_one( a, value )

    Type( matrix )           :: a
    Real( wp ), Intent( In ) :: value

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_ASSIGN'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          a%real_data = value

       Case( MATRIX_COMPLEX )
          a%complex_data = value

    End Select

  End Subroutine matrix_assign_one

  Subroutine matrix_assign_multi( a, value )

    Type( matrix ), Dimension( : )               :: a
    Real( wp )                    , Intent( In ) :: value

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_ASSIGN'

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_assign( a( i ), value )
       End If
    End Do

  End Subroutine matrix_assign_multi

  Subroutine matrix_set_element( a, value, j, i )

    Type( matrix )               :: a
    Real( wp )    , Intent( In ) :: value
    Integer       , Intent( In ) :: j
    Integer       , Intent( In ) :: i

    Integer :: local_i, local_j

    Character( Len = 18 ), Parameter :: routine = 'MATRIX_SET_ELEMENT'

    Call matrix_check_valid( routine, a )

    local_i = a%col%global_to_local( i )
    local_j = a%row%global_to_local( j )

    If( local_i /= MATRIX_NOT_OWN .And. local_j /= MATRIX_NOT_OWN ) Then

       Select Case( a%genus )
          
          Case( MATRIX_REAL )
             a%real_data( local_j, local_i ) = value
             
          Case( MATRIX_COMPLEX )
             a%complex_data( local_j, local_i ) = value

       End Select

    End If

  End Subroutine matrix_set_element

  Subroutine matrix_set_1d( a, data )

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: data

    Integer :: icount
    Integer :: i, j

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_SET_1D'

    Call matrix_check_valid( routine, a )

    icount = 0

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do  i = 1, a%local_m
             Do  j = 1, a%local_n
                icount = icount + 1
                a%real_data( j, i ) = data( icount )
             End Do
          End Do

       Case( MATRIX_COMPLEX )
          Do  i = 1, a%local_m
             Do  j = 1, a%local_n
                icount = icount + 1
                a%complex_data( j, i ) = data( icount )
             End Do
          End Do

    End Select

  End Subroutine matrix_set_1d

  Subroutine matrix_set_from_triangle( a, data )

    ! Even in a UHF_K case this can be called by processors that do NOT
    ! own any data for a - simplifies the interface.

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: data

    Integer :: icount
    Integer :: row_i, col_i
    Integer :: row_j, col_j
    Integer :: i, j

    Character( Len = 24 ), Parameter :: routine = 'MATRIX_SET_FROM_TRIANGLE'

    If( a%i_own ) Then

!!$    Call matrix_check_valid( routine, a )

       Select Case( a%genus )

       Case( MATRIX_REAL )

          icount = 0
          Do  i = 1, a%m
             col_i = a%col%global_to_local( i )
             If( col_i /= MATRIX_NOT_OWN ) Then
                Do  j = 1, i
                   row_j = a%row%global_to_local( j )
                   icount = icount + 1
                   If( row_j /= MATRIX_NOT_OWN ) Then
                      a%real_data( row_j, col_i ) = data( icount )
                   End If
                End Do
             Else
                icount = icount + i
             End If
          End Do

          icount = 0
          Do  i = 1, a%m
             row_i = a%row%global_to_local( i )
             If( row_i /= MATRIX_NOT_OWN ) Then
                Do  j = 1, i
                   col_j = a%col%global_to_local( j )
                   icount = icount + 1
                   If( col_j /= MATRIX_NOT_OWN ) Then
                      a%real_data( row_i, col_j ) = data( icount )
                   End If
                End Do
             Else
                icount = icount + i
             End If
          End Do

          ! Complex screwed anyway
!!$       Case( MATRIX_COMPLEX )
!!$          Do  i = 1, a%local_m
!!$             Do  j = 1, i
!!$                icount = icount + 1
!!$                a%complex_data( j, i ) = data( icount )
!!$                a%complex_data( i, j ) = data( icount )
!!$             End Do
!!$          End Do

       End Select

    End If

  End Subroutine matrix_set_from_triangle

  Subroutine matrix_set_from_square( a, data )

    ! Even in a UHF_K case this can be called by processors that do NOT
    ! own any data for a - simplifies the interface.

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: data

    Integer :: icount
    Integer :: row_i, col_i
    Integer :: row_j, col_j
    Integer :: i, j

    Character( Len = 24 ), Parameter :: routine = 'MATRIX_SET_FROM_TRIANGLE'

    If( a%i_own ) Then

!!$    Call matrix_check_valid( routine, a )

       Select Case( a%genus )

       Case( MATRIX_REAL )

          icount = 0
          Do  i = 1, a%m
             col_i = a%col%global_to_local( i )
             If( col_i /= MATRIX_NOT_OWN ) Then
                Do  j = 1, a%n
                   row_j = a%row%global_to_local( j )
                   icount = icount + 1
                   If( row_j /= MATRIX_NOT_OWN ) Then
                      a%real_data( row_j, col_i ) = data( icount )
                   End If
                End Do
             Else
                icount = icount + a%n
             End If
          End Do

          ! Complex screwed anyway
!!$       Case( MATRIX_COMPLEX )
!!$          Do  i = 1, a%local_m
!!$             Do  j = 1, i
!!$                icount = icount + 1
!!$                a%complex_data( j, i ) = data( icount )
!!$                a%complex_data( i, j ) = data( icount )
!!$             End Do
!!$          End Do

       End Select

    End If

  End Subroutine matrix_set_from_square

  Subroutine matrix_set_random( a )

    Type( matrix ) :: a

    Character( Len = 17 ), Parameter :: routine = 'MATRIX_SET_RANDOM'

    Real( wp ), Dimension( :, : ), Allocatable :: real_part
    Real( wp ), Dimension( :, : ), Allocatable :: imag_part

    Integer :: error

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Call Random_number( a%real_data( 1:a%local_n, 1:a%local_m ) )
          
       Case( MATRIX_COMPLEX )
          Allocate( real_part( 1:a%local_n, 1:a%local_m ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
          End If
          Call alloc_add( ALLOC_REAL, Size( real_part ) )
          Allocate( imag_part( 1:a%local_n, 1:a%local_m ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
          End If
          Call alloc_add( ALLOC_REAL, Size( imag_part ) )
          Call Random_number( real_part )
          Call Random_number( imag_part )
          a%complex_data( 1:a%local_n, 1:a%local_m ) = Cmplx( real_part, imag_part, Kind = wp )
          Call alloc_free( ALLOC_REAL, Size( imag_part ) )
          Deallocate( imag_part )
          Call alloc_free( ALLOC_REAL, Size( real_part ) )
          Deallocate( real_part )

    End Select

  End Subroutine matrix_set_random

  Subroutine matrix_get_element( a, value, j, i )

    Type( matrix )                  :: a
    Real( wp )    , Intent(   Out ) :: value
    Integer       , Intent( In    ) :: j
    Integer       , Intent( In    ) :: i

    Character( Len = 18 ), Parameter :: routine = 'MATRIX_GET_ELEMENT'

    Call matrix_check_valid( routine, a )

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

       Select Case( a%genus )

       Case( MATRIX_REAL )
          value = a%real_data( j, i )

       Case( MATRIX_COMPLEX )
          value = a%complex_data( j, i )

       End Select

    Else

       Call pdelget( 'A', ' ', value, a%real_data, j, i, a%distribution%descriptor )

    End If

  End Subroutine matrix_get_element

  Subroutine matrix_get_1d( a, data )

    Type( matrix )                            :: a
    Real( wp ), Dimension( : ), Intent( Out ) :: data

    Integer :: icount
    Integer :: i, j

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_GET_1D'

    Call matrix_check_valid( routine, a )

    icount = 0

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do  i = 1, a%local_m
             Do  j = 1, a%local_n
                icount = icount + 1
                data( icount ) = a%real_data( j, i )
             End Do
          End Do

       Case( MATRIX_COMPLEX )
          Do  i = 1, a%local_m
             Do  j = 1, a%local_n
                icount = icount + 1
                data( icount ) = a%complex_data( j, i )
             End Do
          End Do

    End Select

  End Subroutine matrix_get_1d

  Subroutine matrix_get_to_triangle( a, data )

    ! To simplify the interface any processor can call this, even
    ! if it does not own the data for a. AS A RESULT THIS ROUTINE IS
    ! GLOBALLY BLOCKING

    Use comms_data

    Type( matrix )                            :: a
    Real( wp ), Dimension( : ), Intent( Out ) :: data

    Integer :: icount
    Integer :: col_i
    Integer :: row_j
    Integer :: i, j

    Character( Len = 22 ), Parameter :: routine = 'MATRIX_GET_TO_TRIANGLE'

!!$    Call matrix_check_valid( routine, a )

    If( a%i_own ) Then

       icount = 0

       Select Case( a%genus )

       Case( MATRIX_REAL )

          Do i = 1, a%m
             col_i = a%col%global_to_local( i )
             If( col_i /= MATRIX_NOT_OWN ) Then
                Do j = 1, i
                   row_j = a%row%global_to_local( j )
                   icount = icount + 1
                   If( row_j /= MATRIX_NOT_OWN ) Then
                      data( icount ) = a%real_data( row_j, col_i )
                   Else
                      data( icount ) = 0.0_wp
                   End If
                End Do
             Else
                data( icount + 1:icount + i ) = 0.0_wp
                icount = icount + i
             End If
          End Do

          ! Complex buggered anyway
!!$       Case( MATRIX_COMPLEX )
!!$          Do  i = 1, a%local_m
!!$             Do  j = 1, a%local_n
!!$                icount = icount + 1
!!$                data( icount ) = a%complex_data( j, i )
!!$             End Do
!!$          End Do

       End Select

    Else

       data = 0.0_wp

    End If

    If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
       Call replicate_array( ( a%n * ( a%n + 1 ) ) / 2, &
            comms_data_all_comm, data )
    End If

  End Subroutine matrix_get_to_triangle

  Subroutine matrix_get_to_square( a, data )

    ! To simplify the interface any processor can call this, even
    ! if it does not own the data for a. AS A RESULT THIS ROUTINE IS
    ! GLOBALLY BLOCKING

    Use comms_data

    Type( matrix )                            :: a
    Real( wp ), Dimension( : ), Intent( Out ) :: data

    Integer :: icount
    Integer :: col_i
    Integer :: row_j
    Integer :: i, j

    Character( Len = 22 ), Parameter :: routine = 'MATRIX_GET_TO_TRIANGLE'

!!$    Call matrix_check_valid( routine, a )

    If( a%i_own ) Then

       icount = 0

       Select Case( a%genus )

       Case( MATRIX_REAL )

          Do i = 1, a%m
             col_i = a%col%global_to_local( i )
             If( col_i /= MATRIX_NOT_OWN ) Then
                Do j = 1, a%n
                   row_j = a%row%global_to_local( j )
                   icount = icount + 1
                   If( row_j /= MATRIX_NOT_OWN ) Then
                      data( icount ) = a%real_data( row_j, col_i )
                   Else
                      data( icount ) = 0.0_wp
                   End If
                End Do
             Else
                data( icount + 1:icount + a%n ) = 0.0_wp
                icount = icount + a%n
             End If
          End Do

!!$          If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
!!$             Call replicate_array( a%n * a%m, a%distribution%communicator, data )
!!$          End If

          ! Complex buggered anyway
!!$       Case( MATRIX_COMPLEX )
!!$          Do  i = 1, a%local_m
!!$             Do  j = 1, a%local_n
!!$                icount = icount + 1
!!$                data( icount ) = a%complex_data( j, i )
!!$             End Do
!!$          End Do

       End Select

    Else

       data = 0.0_wp

    End If

    If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
       Call replicate_array( a%n * a%m, comms_data_all_comm, data )
    End If

  End Subroutine matrix_get_to_square

  Subroutine matrix_set_column( a, data, i  )

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: data
    Integer                   , Intent( In ) :: i

    Character( Len = 17 ), Parameter :: routine = 'MATRIX_SET_COLUMN'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          a%real_data( 1:a%local_n, i ) = data

       Case( MATRIX_COMPLEX )
          a%complex_data( 1:a%local_n, i ) = data

    End Select

  End Subroutine matrix_set_column

  Subroutine matrix_set_row( a, data, i  )

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: data
    Integer                   , Intent( In ) :: i

    Character( Len = 14 ), Parameter :: routine = 'MATRIX_SET_ROW'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          a%real_data( i, 1:a%local_m ) = data

       Case( MATRIX_COMPLEX )
          a%complex_data( i, 1:a%local_m ) = data

    End Select

  End Subroutine matrix_set_row

  Subroutine matrix_get_column( a, data, i  )

    Type( matrix )                              :: a
    Real( wp ), Dimension( : ), Intent(   Out ) :: data
    Integer                   , Intent( In )    :: i

    Character( Len = 17 ), Parameter :: routine = 'MATRIX_GET_COLUMN'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
           data = a%real_data( 1:a%local_n, i )

       Case( MATRIX_COMPLEX )
          data = a%complex_data( 1:a%local_n, i )

    End Select

  End Subroutine matrix_get_column

  Subroutine matrix_get_row( a, data, i  )

    Type( matrix )                              :: a
    Real( wp ), Dimension( : ), Intent(   Out ) :: data
    Integer                   , Intent( In )    :: i

    Character( Len = 14 ), Parameter :: routine = 'MATRIX_GET_ROW'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          data = a%real_data( i, : )

       Case( MATRIX_COMPLEX )
          data = a%complex_data( i, : )

    End Select

  End Subroutine matrix_get_row

  Subroutine matrix_length_one( a, result )
    
    Include 'mpif.h'

    Type( matrix )            :: a
    Real( wp ), Intent( Out ) :: result

    Integer    :: error
    Real( wp ) :: work

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_LENGTH'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          work = Sum( a%real_data( 1:a%local_n, 1:a%local_m ) ** 2 )

       Case( MATRIX_COMPLEX )
          work = Sum( a%complex_data( 1:a%local_n, 1:a%local_m ) * &
               Conjg( a%complex_data( 1:a%local_n, 1:a%local_m ) ) )

    End Select

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then
       result = work
    Else
       Call mpi_allreduce( work, result, 1, mpi_double_precision, &
            mpi_sum, a%distribution%communicator, error )
    End If

    result = Sqrt( result )

  End Subroutine matrix_length_one

  Subroutine matrix_length_multi( a, result )

    Use comms_data
    
    Include 'mpif.h'

    Type( matrix ), Dimension( : )                :: a
    Real( wp )    , Dimension( : ), Intent( Out ) :: result

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_LENGTH_MULTI'

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer :: error
    Integer :: i

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work ) )

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then

          Select Case( a( i )%genus )

             Case( MATRIX_REAL )
                work( i ) = Sum( a( i )%real_data( 1:a( i )%local_n, 1:a( i )%local_m ) ** 2 )
                
             Case( MATRIX_COMPLEX )
                work( i ) = Sum( a( i )%complex_data( 1:a( i )%local_n, 1:a( i )%local_m ) * &
                     Conjg( a( i )%complex_data( 1:a( i )%local_n, 1:a( i )%local_m ) ) )

          End Select

       Else
          
          work( i ) = 0.0_wp

       End If
    End Do

    If( a( 1 )%distribution%type == MATRIX_DISTRIB_REPL ) Then
       result = work
    Else
       Call mpi_allreduce( work, result, Size( work ), mpi_double_precision, &
            mpi_sum, comms_data_all_comm, error )
    End If

    Call alloc_free( ALLOC_REAL, Size( work ) )
    Deallocate( work )

    result = Sqrt( result )

  End Subroutine matrix_length_multi

  Subroutine matrix_dot_product_one( a, b, result, dont_sum )
    
    ! Complex case screwed !

    Use comms_data

    Include 'mpif.h'

    Type( matrix )                     :: a
    Type( matrix )                     :: b
    Real( wp )       , Intent(   Out ) :: result
    Logical, Optional, Intent( In    ) :: dont_sum

    Character( Len = 18 ), Parameter :: routine = 'MATRIX_DOT_PRODUCT'

    Real( wp ) :: work

    Integer :: error

    Logical :: do_sum

    Call matrix_check_valid( routine, a )

    Call matrix_check_conform( routine, a, b )

    If( .Not. Present( dont_sum ) ) Then
       do_sum = .True.
    Else
       do_sum = .Not. dont_sum
    End If

    Select Case( a%genus )

       Case( MATRIX_REAL )
          result = Sum( a%real_data( 1:a%local_n, 1:a%local_m ) * &
                        b%real_data( 1:b%local_n, 1:b%local_m ) )

       Case( MATRIX_COMPLEX )
          result = Sum( a%complex_data( 1:a%local_n, 1:a%local_m ) * & 
                 Conjg( b%complex_data( 1:b%local_n, 1:b%local_m ) ) )

    End Select

    If( do_sum ) Then
       If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
          work = result
          Call mpi_allreduce( work, result, 1, mpi_double_precision, &
               mpi_sum, comms_data_all_comm, error )
       End If
    End If


  End Subroutine matrix_dot_product_one

  Subroutine matrix_dot_product_multi( a, b, result )
    
    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( : )                :: a
    Type( matrix ), Dimension( : )                :: b
    Real( wp )    , Dimension( : ), Intent( Out ) :: result

    Character( Len = 24 ), Parameter :: routine = 'MATRIX_DOT_PRODUCT_MULTI'

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer :: error
    Integer :: i

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work ) )

    work = 0.0_wp

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_dot_product( a( i ), b( i ), work( i ), .True. )
       End If
    End Do

    If( a( 1 )%distribution%type == MATRIX_DISTRIB_REPL ) Then
       result = work
    Else
       Call mpi_allreduce( work, result, Size( work ), mpi_double_precision, &
            mpi_sum, comms_data_all_comm, error )
    End If

    Call alloc_free( ALLOC_REAL, Size( work ) )
    Deallocate( work )

  End Subroutine matrix_dot_product_multi

  Subroutine matrix_daxpy_one( alpha, a, b )

    Real( wp ), Intent( In ) :: alpha
    Type( matrix )           :: a
    Type( matrix )           :: b

    Character( Len = 12 ), Parameter :: routine = 'MATRIX_DAXPY'

    Call matrix_check_valid( routine, a )

    Call matrix_check_conform( routine, a, b )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          b%real_data( 1:b%local_n, 1:b%local_m ) = b%real_data( 1:b%local_n, 1:b%local_m ) + &
               alpha * a%real_data( 1:a%local_n, 1:a%local_m )

       Case( MATRIX_COMPLEX )
          b%complex_data( 1:b%local_n, 1:b%local_m ) = b%complex_data( 1:b%local_n, 1:b%local_m ) + &
               alpha * a%complex_data( 1:a%local_n, 1:a%local_m )

    End Select

  End Subroutine matrix_daxpy_one

  Subroutine matrix_daxpy_multi( alpha, a, b )

    Real( wp )                    , Intent( In ) :: alpha
    Type( matrix ), Dimension( : )               :: a
    Type( matrix ), Dimension( : )               :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_daxpy( alpha, a( i ), b( i ) )
       End If
    End Do

  End Subroutine matrix_daxpy_multi

  Subroutine matrix_copy_one( a, b )

    Use comms_data

    Include 'mpif.h'
    
    Type( matrix ) :: a
    Type( matrix ) :: b

    Complex( wp ), Dimension( :, : ), Allocatable :: cwork

    Real( wp ), Dimension( :, : ), Allocatable :: rwork

    Integer i, j
    Integer :: error

    Logical :: a_repl, b_repl

    Character( Len = 12 ), Parameter :: routine = 'MATRIX_COPY'

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )

    Call matrix_check_conform( routine, a, b, .False. )

    a_repl = a%distribution%type == MATRIX_DISTRIB_REPL
    b_repl = b%distribution%type == MATRIX_DISTRIB_REPL

    If( (       a_repl .And.       b_repl ) .Or. &
        ( .Not. a_repl .And. .Not. b_repl ) .Or. &
        (       a_repl .And. .Not. b_repl ) ) Then

       ! Both distributed, both replicated, or A replicated and B distributed
       ! No comms required

       ! At present can distrib to/from all and uhf_k
       If( .Not. a_repl .And. .Not. b_repl ) Then
          If( a%distribution%type /= b%distribution%type ) Then
             Call dlc_error( 'Disparate matrix distributions in ' // routine, 'abort' )
          End If
       End If

       If( a_repl .Neqv. b_repl ) Then

          Select Case( a%genus )

             Case( MATRIX_REAL )
                Do i = 1, b%local_m
                   Do j = 1, b%local_n
                      b%real_data( j, i ) = a%real_data( b%row%local_to_global( j ), &
                           b%col%local_to_global( i ) )
                   End Do
                End Do

             Case( MATRIX_COMPLEX )
                Do i = 1, b%local_m
                   Do j = 1, b%local_n
                      b%complex_data( j, i ) = a%complex_data( b%row%local_to_global( j ), &
                           b%col%local_to_global( i ) )
                   End Do
                End Do

          End Select

       Else

          Select Case( a%genus )

             Case( MATRIX_REAL )
                Do i = 1, b%local_m
                   Do j = 1, b%local_n
                      b%real_data( j, i ) = a%real_data( j, i )
                   End Do
                End Do

             Case( MATRIX_COMPLEX )
                Do i = 1, b%local_m
                   Do j = 1, b%local_n
                      b%complex_data( j, i ) = a%complex_data( j, i )
                   End Do
                End Do

          End Select

       End If

    Else

       ! A distributed, B replicated - Needs comms to replicate B
       
       Select Case( a%genus )

          Case( MATRIX_REAL )
             Allocate( rwork( 1:b%local_n, b%local_m ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Matrix allocation error in MATRIX_COPY', 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( rwork ) )
             rwork = 0.0_wp
             Do i = 1, a%local_m
                Do j = 1, a%local_n
                   rwork( a%row%local_to_global( j ), a%col%local_to_global( i ) ) = &
                        a%real_data( j, i )
                End Do
             End Do
             Call mpi_allreduce( rwork, b%real_data, Size( rwork ), mpi_double_precision, &
                  mpi_sum, a%distribution%communicator, error )
             Call alloc_free( ALLOC_REAL, Size( rwork ) )
             Deallocate( rwork )

          Case( MATRIX_COMPLEX )
             Allocate( cwork( 1:b%local_n, b%local_m ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Matrix allocation error in MATRIX_COPY', 'abort' )
             End If
             Call alloc_add( ALLOC_COMPLEX, Size( cwork ) )
             cwork = 0.0_wp
             Do i = 1, a%local_m
                Do j = 1, a%local_n
                   cwork( a%row%local_to_global( j ), a%col%local_to_global( i ) ) = &
                        a%complex_data( j, i )
                End Do
             End Do
             Call mpi_allreduce( cwork, b%complex_data, 2 * Size( cwork ), mpi_double_precision, &
                  mpi_sum, a%distribution%communicator, error )
             Call alloc_free( ALLOC_COMPLEX, Size( cwork ) )
             Deallocate( cwork )

       End Select

    End If

  End Subroutine matrix_copy_one

  Subroutine matrix_copy_multi( a, b )

    Use comms_data
    
    Type( matrix ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own .And. b( i )%i_own ) Then
          Call matrix_copy( a( i ), b( i ) )
       End If
    End Do

    ! If result matrix replicated need to, duh, replicate it
    If( a( 1 )%distribution%type == MATRIX_DISTRIB_UHF_K ) Then
       If( b( 1 )%distribution%type == MATRIX_DISTRIB_REPL ) Then
          Call matrix_replicate_uhf_k( b )
       End If
    End If

  End Subroutine matrix_copy_multi

  Subroutine matrix_copy_cols_one( col_start_a, col_end_a, a, &
                                   col_start_b, col_end_b, b )

    Include 'mpif.h'
    
    Integer, Intent( In ) :: col_start_a
    Integer, Intent( In ) :: col_end_a
    Type( matrix )        :: a
    Integer, Intent( In ) :: col_start_b
    Integer, Intent( In ) :: col_end_b
    Type( matrix )        :: b

    Character( Len = 16 ), Parameter :: routine = 'MATRIX_COPY_COLS'

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )

    If( a%distribution%type /= b%distribution%type ) Then
       Call dlc_error( 'Inconsistent matrix distributions in ' // routine, 'abort' )
    End If

    If( col_end_a - col_start_a /= col_end_b - col_start_b ) Then
       Call dlc_error( 'Inconsistent dimensions in ' // routine, 'abort' )
    End If


    Call pdgemr2d( a%n, col_end_a - col_start_a + 1,             &
         a%real_data, 1, col_start_a, a%distribution%descriptor, &
         b%real_data, 1, 1,           b%distribution%descriptor, &
         a%distribution%context )

  End Subroutine matrix_copy_cols_one

  Subroutine matrix_copy_cols_multi( col_start_a, col_end_a, a, &
                                     col_start_b, col_end_b, b )

    Include 'mpif.h'
    
    Integer                       , Intent( In ) :: col_start_a
    Integer                       , Intent( In ) :: col_end_a
    Type( matrix ), Dimension( : )               :: a
    Integer                       , Intent( In ) :: col_start_b
    Integer                       , Intent( In ) :: col_end_b
    Type( matrix ), Dimension( : )               :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_copy_cols( col_start_a, col_end_a, a( i ), &
                                 col_start_b, col_end_b, b( i ) )
       End If
    End do

  End Subroutine matrix_copy_cols_multi

  Subroutine matrix_scale( a, alpha )

    Real( wp ), Intent( In ) :: alpha
    Type( matrix )           :: a

    Character( Len = 12 ), Parameter :: routine = 'MATRIX_SCALE'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          a%real_data( 1:a%local_n, 1:a%local_m ) = alpha * a%real_data( 1:a%local_n, 1:a%local_m )

       Case( MATRIX_COMPLEX )
          a%complex_data( 1:a%local_n, 1:a%local_m ) = alpha * a%complex_data( 1:a%local_n, 1:a%local_m )

    End Select

  End Subroutine matrix_scale

  Integer Function matrix_dimension( a, dimension )

    ! RETURNS GLOBAL SIZES !!!!!!!!!!
    ! Further .....
    ! NOT RECOMMENDED !!!!!!!!!!
    ! Use MATRIX_INQUIRE instead

    Type( matrix )        :: a
    Integer, Intent( In ) :: dimension

    Character( Len = 16 ), Parameter :: routine = 'MATRIX_DIMENSION'

    Call matrix_check_valid( routine, a )

    Select Case( dimension )

       Case( 1 )
          matrix_dimension = a%n

       Case( 2 )
          matrix_dimension = a%m
          
       Case Default
          Call dlc_error( 'Invalid dimension in ' // routine, 'abort' )

    End Select

  End Function matrix_dimension

  Subroutine matrix_inquire( a, i_own, genus, distrib, global_n, global_m, &
                                                       local_n , local_m , &
                                                       local_ld, local_sd, &
                                                       name )
    
    Type( matrix )                                     :: a
    Integer             , Optional, Intent( Out )      :: distrib
    Logical             , Optional, Intent( Out )      :: i_own
    Integer             , Optional, Intent( Out )      :: genus
    Integer             , Optional, Intent( Out )      :: global_n
    Integer             , Optional, Intent( Out )      :: global_m
    Integer             , Optional, Intent( Out )      :: local_n
    Integer             , Optional, Intent( Out )      :: local_m
    Integer             , Optional, Intent( Out )      :: local_ld
    Integer             , Optional, Intent( Out )      :: local_sd
    Character( Len = * ), Optional, Intent( Out )      :: name
    
    Character( Len = 14 ), Parameter :: routine = 'MATRIX_INQUIRE'

    If( Present( distrib ) ) Then
       distrib = a%distribution%type
    End If
    
    If( Present( i_own ) ) Then
       i_own = a%i_own
    End If
    
    If( Present( genus ) ) Then
       genus = a%genus
    End If
    
    If( Present( global_n ) ) Then
       global_n = a%n
    End If
    
    If( Present( global_m ) ) Then
       global_m = a%m
    End If
    
    If( Present( local_n ) ) Then
       local_n = a%local_n
    End If
    
    If( Present( local_m ) ) Then
       local_m = a%local_m
    End If
    
    If( Present( local_ld ) ) Then
       local_ld = a%local_ld
    End If
    
    If( Present( local_sd ) ) Then
       local_sd = a%local_sd
    End If
    
    If( Present( name ) ) Then
       name = a%name
    End If
    
  End Subroutine matrix_inquire

  Subroutine matrix_multiply_one( alpha, a, b, beta, c )
  
    Real( wp ), Intent( In ) :: alpha
    Type( matrix )           :: a
    Type( matrix )           :: b
    Real( wp ), Intent( In ) :: beta
    Type( matrix )           :: c

    Call matrix_dgemm( 'n', 'n', alpha, a, b, beta, c )

  End Subroutine matrix_multiply_one

  Subroutine matrix_multiply_multi( alpha, a, b, beta, c )

    Real( wp )                    , Intent( In ) :: alpha
    Type( matrix ), Dimension( : )               :: a
    Type( matrix ), Dimension( : )               :: b
    Real( wp )                    , Intent( In ) :: beta
    Type( matrix ), Dimension( : )               :: c

    Integer :: i

    Call matrix_dgemm( 'n', 'n', alpha, a, b, beta, c )

  End Subroutine matrix_multiply_multi

  Subroutine matrix_dgemm_one( t1, t2, alpha, a, b, beta, c )
  
    Character , Intent( In ) :: t1
    Character , Intent( In ) :: t2
    Real( wp ), Intent( In ) :: alpha
    Type( matrix )           :: a
    Type( matrix )           :: b
    Real( wp ), Intent( In ) :: beta
    Type( matrix )           :: c

    Character( Len = 15 ), Parameter :: routine = 'MATRIX_DGEMM'

    Character :: t1_here
    Character :: t2_here

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )
    Call matrix_check_valid( routine, c )

    If( a%genus /= b%genus .Or. a%genus /= c%genus ) Then
       Call dlc_error( 'Mismatched matrix types in ' // routine, 'abort' )
    End If

    If( .Not. a%i_own .Or. .Not. b%i_own .Or. .Not. c%i_own ) Then
       Call dlc_error( 'Disparate matrix distributions in ' // routine, 'abort' )
    End If

    t1_here = matrix_to_lower( t1 )
    If( t1_here == 'c' .And. a%genus == MATRIX_COMPLEX ) Then
       t1_here = 't'
    End If
    t2_here = matrix_to_lower( t2 )
    If( t2_here == 'c' .And. a%genus == MATRIX_COMPLEX ) Then
       t2_here = 't'
    End If

    Select Case( t1_here )

       Case( 'n' )
          Select Case( t2_here )
             Case( 'n' )
                Call dgemm_nn
             Case( 't' )
                Call dgemm_nt
             Case Default
                Call dlc_error( 'Invalid transpose option in ' // routine, 'abort' )
          End Select

       Case( 't')
          Select Case( t2_here )
             Case( 'n' )
                Call dgemm_tn
             Case( 't' )
                Call dgemm_tt
             Case Default
                Call dlc_error( 'Invalid transpose option in ' // routine, 'abort' )
          End Select

       Case Default
          Call dlc_error( 'Invalid transpose option in ' // routine, 'abort' )

    End Select

  Contains

    Subroutine dgemm_nn

      If( a%m /= b%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 2 ) and Dim( b, 1 )', &
              'abort' ) 
      End If

      If( a%n /= c%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 1 ) and Dim( c, 1 )', &
              'abort' ) 
      End If

      If( b%m /= c%m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( b, 2 ) and Dim( c, 2 )', &
              'abort' ) 
      End If
      
      If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

         Select Case( a%genus )
            
            Case( MATRIX_REAL )
               Call dgemm( 'n', 'n', c%n, c%m, b%n, alpha, a%real_data, a%local_ld, &
                                                           b%real_data, b%local_ld, &
                                                    beta , c%real_data, c%local_ld )

            Case( MATRIX_COMPLEX )
               Call zgemm( 'n', 'n', c%n, c%m, b%n, alpha, a%complex_data, a%local_ld, &
                                                           b%complex_data, b%local_ld, &
                                                    beta , c%complex_data, c%local_ld )

         End Select

      Else

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call pdgemm( 'n', 'n', c%n, c%m, b%n, alpha, a%real_data, 1, 1, a%distribution%descriptor, &
                                                            b%real_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%real_data, 1, 1, c%distribution%descriptor )
            Case( MATRIX_COMPLEX )
               Call pzgemm( 'n', 'n', c%n, c%m, b%n, alpha, a%complex_data, 1, 1, a%distribution%descriptor, &
                                                            b%complex_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%complex_data, 1, 1, c%distribution%descriptor )

         End Select

      End If

    End Subroutine dgemm_nn

    Subroutine dgemm_nt

      If( a%m /= b%m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 2 ) and Dim( b, 1 )', &
              'abort' ) 
      End If

      If( a%n /= c%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 1 ) and Dim( c, 1 )', &
              'abort' ) 
      End If

      If( b%n /= c%m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( b, 2 ) and Dim( c, 2 )', &
              'abort' ) 
      End If

      If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call dgemm( 'n', 't', c%n, c%m, b%m, alpha, a%real_data, a%local_ld, &
                                                           b%real_data, b%local_ld, &
                                                    beta , c%real_data, c%local_ld )

            Case( MATRIX_COMPLEX )
               Call zgemm( 'n', 'c', c%n, c%m, b%m, alpha, a%complex_data, a%local_ld, &
                                                           b%complex_data, b%local_ld, &
                                                    beta , c%complex_data, c%local_ld )

         End Select

      Else

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call pdgemm( 'n', 't', c%n, c%m, b%m, alpha, a%real_data, 1, 1, a%distribution%descriptor, &
                                                            b%real_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%real_data, 1, 1, c%distribution%descriptor )
            Case( MATRIX_COMPLEX )
               Call pzgemm( 'n', 't', c%n, c%m, b%m, alpha, a%complex_data, 1, 1, a%distribution%descriptor, &
                                                            b%complex_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%complex_data, 1, 1, c%distribution%descriptor )

         End Select

      End If

    End Subroutine dgemm_nt

    Subroutine dgemm_tn

      If( a%n /= b%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 2 ) and Dim( b, 1 )', &
              'abort' ) 
      End If

      If( a%m /= c%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 1 ) and Dim( c, 1 )', &
              'abort' ) 
      End If

      If( b%m /= c%m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( b, 2 ) and Dim( c, 2 )', &
              'abort' ) 
      End If

      If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call dgemm( 't', 'n', c%n, c%m, b%n, alpha, a%real_data, a%local_ld, &
                                                           b%real_data, b%local_ld, &
                                                    beta , c%real_data, c%local_ld )

            Case( MATRIX_COMPLEX )
               Call zgemm( 'c', 'n', c%n, c%m, b%n, alpha, a%complex_data, a%local_ld, &
                                                           b%complex_data, b%local_ld, &
                                                    beta , c%complex_data, c%local_ld )

         End Select

      Else

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call pdgemm( 't', 'n', c%n, c%m, b%n, alpha, a%real_data, 1, 1, a%distribution%descriptor, &
                                                            b%real_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%real_data, 1, 1, c%distribution%descriptor )
            Case( MATRIX_COMPLEX )
               Call pzgemm( 'c', 'n', c%n, c%m, b%n, alpha, a%complex_data, 1, 1, a%distribution%descriptor, &
                                                            b%complex_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%complex_data, 1, 1, c%distribution%descriptor )

         End Select

      End If

    End Subroutine dgemm_tn

    Subroutine dgemm_tt

      If( a%n /= b%m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 2 ) and Dim( b, 1 )', &
              'abort' ) 
      End If

      If( a%m /= c%n ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( a, 1 ) and Dim( c, 1 )', &
              'abort' ) 
      End If

      If( b%local_n /= c%local_m ) Then
         Call dlc_error( 'Mismatched dimensions in ' // routine // ' : Dim( b, 2 ) and Dim( c, 2 )', &
              'abort' ) 
      End If

      If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call dgemm( 't', 't', c%n, c%m, b%m, alpha, a%real_data, a%local_ld, &
                                                           b%real_data, b%local_ld, &
                                                    beta , c%real_data, c%local_ld )

            Case( MATRIX_COMPLEX )
               Call zgemm( 'c', 'c', c%n, c%m, b%m, alpha, a%complex_data, a%local_ld, &
                                                           b%complex_data, b%local_ld, &
                                                    beta , c%complex_data, c%local_ld )
         End Select

      Else

         Select Case( a%genus )

            Case( MATRIX_REAL )
               Call pdgemm( 't', 't', c%n, c%m, b%m, alpha, a%real_data, 1, 1, a%distribution%descriptor, &
                                                            b%real_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%real_data, 1, 1, c%distribution%descriptor )
            Case( MATRIX_COMPLEX )
               Call pzgemm( 'c', 'c', c%n, c%m, b%m, alpha, a%complex_data, 1, 1, a%distribution%descriptor, &
                                                            b%complex_data, 1, 1, b%distribution%descriptor, & 
                                                      beta, c%complex_data, 1, 1, c%distribution%descriptor )

         End Select

      End If


    End Subroutine dgemm_tt

  End Subroutine matrix_dgemm_one

  Subroutine matrix_dgemm_multi( t1, t2, alpha, a, b, beta, c )
  
    Character ,                    Intent( In ) :: t1
    Character ,                    Intent( In ) :: t2
    Real( wp ),                    Intent( In ) :: alpha
    Type( matrix ), Dimension( : )              :: a
    Type( matrix ), Dimension( : )              :: b
    Real( wp ),                    Intent( In ) :: beta
    Type( matrix ), Dimension( : )              :: c

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_dgemm( t1, t2, alpha, a( i ), b( i ), beta, c( i ) )
       End If
    End Do

  End Subroutine matrix_dgemm_multi

  Subroutine matrix_mult_by_diag_one( a, diag, b )

    Type( matrix )                           :: a
    Real( wp ), Dimension( : ), Intent( In ) :: diag
    Type( matrix )                           :: b

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_MULT_BY_DIAG'

    Integer :: col
    Integer :: i

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )

    Call matrix_check_conform( routine, a, b )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = 1, a%local_m
             col = a%col%local_to_global( i )
             b%real_data( 1:b%local_n, i ) = diag( col ) * a%real_data( 1:a%local_n, i )
          End Do

       Case( MATRIX_COMPLEX )
          Do i = 1, a%local_m
             col = a%col%local_to_global( i )
             b%complex_data( 1:b%local_n, i ) = diag( col ) * a%complex_data( 1:a%local_n, i )
          End Do

    End Select

  End Subroutine matrix_mult_by_diag_one

  Subroutine matrix_mult_by_diag_multi( a, diag, b )

    Type( matrix ), Dimension( :    )               :: a
    Real( wp )    , Dimension( :, : ), Intent( In ) :: diag
    Type( matrix ), Dimension( :    )               :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_mult_by_diag( a( i ), diag( :, i ), b( i ) )
       End If
    End Do

  End Subroutine matrix_mult_by_diag_multi

  Subroutine matrix_mult2_one( a, b, c )

    Type( matrix ) :: a
    Type( matrix ) :: b
    Type( matrix ) :: c

    Type( matrix ) :: work

    Call matrix_create( b%m, b%n, work, 'mult2 work', a%genus, &
         a%distribution%type, a%i_own )

    Call matrix_dgemm( 't', 'n', 1.0_wp, b   , a, 0.0_wp, work )
    Call matrix_dgemm( 'n', 'n', 1.0_wp, work, b, 0.0_wp, c    )

    Call matrix_destroy( work )

  End Subroutine matrix_mult2_one

  Subroutine matrix_mult2_multi( a, b, c )

    Type( matrix ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: b
    Type( matrix ), Dimension( : ) :: c

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_mult2( a( i ), b( i ), c( i ) )
       End If
    End Do

  End Subroutine matrix_mult2_multi

  Subroutine matrix_mult2t_one( a, b, c )

    Type( matrix ) :: a
    Type( matrix ) :: b
    Type( matrix ) :: c

    Type( matrix ) :: work

    Call matrix_create( b%n, b%m, work, 'mult2t work', a%genus, &
         a%distribution%type, a%i_own )

    Call matrix_dgemm( 'n', 'n', 1.0_wp, b   , a, 0.0_wp, work )
    Call matrix_dgemm( 'n', 't', 1.0_wp, work, b, 0.0_wp, c    )

    Call matrix_destroy( work )

  End Subroutine matrix_mult2t_one

  Subroutine matrix_mult2t_multi( a, b, c )

    Type( matrix ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: b
    Type( matrix ), Dimension( : ) :: c

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_mult2t( a( i ), b( i ), c( i ) )
       End If
    End Do
    
  End Subroutine matrix_mult2t_multi

  Subroutine matrix_combine( alpha, a, beta, b, c )

    Real( wp ), Intent( In ) :: alpha
    Type( matrix )           :: a
    Real( wp ), Intent( In ) :: beta
    Type( matrix )           :: b
    Type( matrix )           :: c

    Character( Len = 14 ), Parameter :: routine = 'MATRIX_COMBINE' 

    Call matrix_check_valid( routine, a )
    Call matrix_check_valid( routine, b )
    Call matrix_check_valid( routine, c )

    Call matrix_check_conform( routine, a, b )
    Call matrix_check_conform( routine, a, c )

    Select Case( c%genus )

       Case( MATRIX_REAL )
          c%real_data( 1:c%local_n, 1:c%local_m ) = alpha * a%real_data( 1:a%local_n, 1:a%local_m ) + &
                                                    beta  * b%real_data( 1:b%local_n, 1:b%local_m )

       Case( MATRIX_COMPLEX )
          c%complex_data( 1:c%local_n, 1:c%local_m ) = alpha * a%complex_data( 1:a%local_n, 1:a%local_m ) + &
                                                       beta  * b%complex_data( 1:b%local_n, 1:b%local_m )

    End Select

  End Subroutine matrix_combine

  Subroutine matrix_abs_one( a )

    Type( matrix ) :: a

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_ABS'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          a%real_data( 1:a%local_n, 1:a%local_m ) = Abs( a%real_data( 1:a%local_n, 1:a%local_m ) ) 

       Case( MATRIX_COMPLEX )
          a%complex_data( 1:a%local_n, 1:a%local_m ) = Abs( a%complex_data( 1:a%local_n, 1:a%local_m ) )

    End Select

  End Subroutine matrix_abs_one

  Subroutine matrix_abs_multi( a )

    Type( matrix ), Dimension( : ) :: a

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_abs( a( i ) )
       End If
    End Do

  End Subroutine matrix_abs_multi

  Subroutine matrix_absmax_one( a, absmax )

    Type( matrix )            :: a
    Real( wp ), Intent( Out ) :: absmax

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_ABSMAX'

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          absmax = Maxval( Abs( a%real_data( 1:a%local_n, 1:a%local_m ) ) ) 

       Case( MATRIX_COMPLEX )
          absmax = Maxval( Abs( a%complex_data( 1:a%local_n, 1:a%local_m ) ) ) 

    End Select

  End Subroutine matrix_absmax_one

  Subroutine matrix_absmax_multi( a, result )

    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( : )                :: a
    Real( wp ),     Dimension( : ), Intent( Out ) :: result

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_ABSMAX_MULTI'

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer :: error
    Integer :: i

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work ) )

    work = -1.0_wp

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_absmax( a( i ), work( i ) )
       End If
    End Do

    If( a( 1 )%distribution%type == MATRIX_DISTRIB_REPL ) Then
       result = work
    Else
       Call mpi_allreduce( work, result, Size( work ), mpi_double_precision, &
            mpi_max, comms_data_all_comm, error )
    End If

    Call alloc_free( ALLOC_REAL, Size( work ) )
    Deallocate( work )

  End Subroutine matrix_absmax_multi

  Subroutine matrix_diagonalise_one( a, evecs, evals )

    Use comms_data

    Include 'mpif.h'

    Type( matrix )                            :: a
    Type( matrix )                            :: evecs
    Real( wp ), Dimension( : ), Intent( Out ) :: evals

    Integer, External :: ilaenv
    Logical, External :: divide_and_conquer_requested

    Character( Len = 18 ), Parameter :: routine = 'MATRIX_DIAGONALISE' 

    Complex( wp ), Dimension( : ), Allocatable  :: imag_work

    Real( wp ), Dimension( : ), Allocatable  :: real_work

    Real( wp ) :: r_work_size

    Integer :: nb
    Integer :: work_size, i_work_size, itmp
    Integer :: error
    Integer :: i, j

!    !Hack for workspace - only for divide and conquer - see below
    Integer, Dimension( : ), Allocatable :: iwork

    Call matrix_check_conform( routine, a, evecs )

    Call matrix_check_valid( routine, a     )
    Call matrix_check_valid( routine, evecs )

    If( a%n /= a%m ) Then
       Call dlc_error( 'In ' // routine // ' trying to diag a rectangular matrix.', 'abort' )
    End If

    If( Size( evals ) < a%n ) Then
       Call dlc_error( 'In ' // routine // ' insufficient space for evals.', 'abort' )
    End If

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

       Call matrix_copy( a, evecs )
       
       Select Case( a%genus )
          
          Case( MATRIX_REAL )
             nb = ilaenv( 1, 'DSYTRD', 'U', a%local_n, -1, -1, -1 )
             work_size = ( nb + 2 ) * a%local_n
             Allocate( real_work( 1:work_size ), stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( real_work ) )
             evecs%real_data = a%real_data
             Call dsyev( 'V', 'U', a%local_n, evecs%real_data, a%local_ld, evals, &
                  real_work, work_size, error )
             Call alloc_free( ALLOC_REAL, Size( real_work ) )
             Deallocate( real_work )

          Case( MATRIX_COMPLEX )
             nb = ilaenv( 1, 'ZHETRD', 'U', a%local_n, -1, -1, -1 )
             work_size = ( nb + 2 ) * a%local_n
             Allocate( imag_work( 1:work_size ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_COMPLEX, Size( imag_work ) )
             Allocate( real_work( 1:3 * a%local_n ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( real_work ) )
             evecs%complex_data = a%complex_data
             Call zheev( 'V', 'U', a%local_n, evecs%complex_data, a%local_ld, evals, &
                  imag_work, work_size, error )
             Call alloc_free( ALLOC_REAL, Size( real_work ) )
             Deallocate( real_work )
             Call alloc_free( ALLOC_COMPLEX, Size( imag_work ) )
             Deallocate( imag_work )

       End Select

    Else

! Interface for divide and conquer routines - but hit a bug ! So go
! back to basic driver
!!$       Call pdsyevd( 'V', 'U', a%n, &
!!$                a%real_data, 1, 1,     a%distribution%descriptor, evals, &
!!$            evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
!!$            work, Size( work ), iwork, Size( iwork ), error )

       Select Case( a%genus )
          
          Case( MATRIX_REAL )

             If (divide_and_conquer_requested()) Then

                ! Workspace inquiry
                i_work_size = 1
                Call pdsyevd( 'V', 'U', a%n, &
                     a%real_data, 1, 1,     a%distribution%descriptor, evals, &
                     evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
                     r_work_size, -1, i_work_size, 1, error )


                work_size = Nint( r_work_size )

!!$                Write( *, * ) work_size, i_work_size
!!$                itmp = work_size
!!$                Call mpi_allreduce( itmp, work_size, 1, mpi_integer, mpi_max, &
!!$                     a%distribution%communicator, error )
                !!!!!!!   THE FACTOR OF FOUR IS A BODGE TO HOPEFULLY GET AROUND A BUG IN PDSYEVD !!!!!!!!!!!!!!
                work_size = 4 * work_size
!!$                i_work_size = 1000000

                Allocate( real_work( 1:work_size ), Stat = error )
                If( error /= 0 ) Then
                   Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
                End If
                Call alloc_add( ALLOC_REAL, Size( real_work ) )

                Allocate( iwork( 1:i_work_size ), Stat = error )
                If( error /= 0 ) Then
                   Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
                End If
                Call alloc_add( ALLOC_INTEGER, Size( iwork ) )

                Call pdsyevd( 'V', 'U', a%n, &
                     a%real_data, 1, 1,     a%distribution%descriptor, evals, &
                     evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
                     real_work, Size( real_work ), iwork, Size( iwork ), error )
!!$                If( comms_data_uhf_k_me == 0 ) Then
!!$                   Write( *, * ) 'Parallel Diag status for psyevd: ', error, ' Work size: ', Size( real_work )
!!$                End If
                If( error /= 0 ) Then
                   Call dlc_error( 'Internal error: Parallel diag failed in ' // routine, 'abort' )
                End If

                Call alloc_free( ALLOC_INTEGER, Size( iwork ) )
                Deallocate( iwork )

                Call alloc_free( ALLOC_REAL, Size( real_work ) )
                Deallocate( real_work )

             Else
             ! Workspace inquiry
             Call pdsyev( 'V', 'U', a%n, &
                  a%real_data, 1, 1,     a%distribution%descriptor, evals, &
                  evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
                  r_work_size, -1, error )
!!$             Call pdsyevd( 'V', 'U', a%n, &
!!$                  a%real_data, 1, 1,     a%distribution%descriptor, evals, &
!!$                  evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
!!$                  r_work_size, -1, i_work_size, 1, error )

             work_size = Nint( r_work_size )
             Allocate( real_work( 1:work_size ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( real_work ) )
!!$             Allocate( iwork( 1:i_work_size ), Stat = error )
!!$             If( error /= 0 ) Then
!!$                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
!!$             End If
!!$             Call alloc_add( ALLOC_INTEGER, Size( iwork ) )

             Call pdsyev( 'V', 'U', a%n, &
                  a%real_data, 1, 1,     a%distribution%descriptor, evals, &
                  evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
                  real_work, Size( real_work ), error )
!!$             Call pdsyevd( 'V', 'U', a%n, &
!!$                  a%real_data, 1, 1,     a%distribution%descriptor, evals, &
!!$                  evecs%real_data, 1, 1, evecs%distribution%descriptor,        &
!!$                  real_work, Size( real_work ), iwork, Size( iwork ), error )
!!$             If( comms_data_uhf_k_me == 0 ) Then
!!$                Write( *, * ) 'Parallel Diag status for psyevd: ', error, ' Work size: ', Size( real_work )
!!$             End If
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: Parallel diag failed in ' // routine, 'abort' )
             End If
!!$             Call alloc_free( ALLOC_INTEGER, Size( iwork ) )
!!$             Deallocate( iwork )
             Call alloc_free( ALLOC_REAL, Size( real_work ) )
             Deallocate( real_work )
             Endif

          Case( MATRIX_COMPLEX )

             Call dlc_error( 'Internal error: Hermitian diag not implemented in ' // routine, 'abort' )

       End Select

    End If

  End Subroutine matrix_diagonalise_one

  Subroutine matrix_diagonalise_multi( a, evecs, evals )

    Use comms_data

    Type( matrix ), Dimension(   :  )                :: a
    Type( matrix ), Dimension(   :  )                :: evecs
    Real( wp )    , Dimension( :, : ), Intent( Out ) :: evals

    Integer :: i, j

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_diagonalise( a( i ), evecs( i ), evals( :, i ) )
       End If
    End Do

  End Subroutine matrix_diagonalise_multi

  Subroutine matrix_invert_one( a )

    Type( matrix ) :: a

    Character( Len = 13 ), Parameter :: routine = 'MATRIX_INVERT'

    Complex( wp ), Dimension( : ), Allocatable :: cwork

    Real( wp ), Dimension( : ), Allocatable :: rwork

    Integer, Dimension( : ), Allocatable :: ipiv
    Integer, Dimension( : ), Allocatable :: iwork

    Real( wp ) :: rlwork

    Integer :: lwork, liwork
    Integer :: error

    Call matrix_check_valid( routine, a )

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

       Select Case( a%genus )
          
          Case( MATRIX_REAL )
             Allocate( ipiv( 1:Min( a%n, a%m ) ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_INTEGER, Size( ipiv ) )
             Call dgetrf( a%n, a%m, a%real_data, a%local_ld, ipiv, error )
             If( error > 0 ) Then
                Call dlc_error( 'Internal error: Matrix singular in ' // routine, 'abort' )
             End If
             Allocate( rwork( 1:a%n ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( rwork ) )
             Call dgetri( a%n, a%real_data, a%local_ld, ipiv, rwork, Size( rwork ), error )
             If( error > 0 ) Then
                Call dlc_error( 'Internal error: Matrix singular in ' // routine, 'abort' )
             End If
             Call alloc_free( ALLOC_REAL, Size( rwork ) )
             Deallocate( rwork )
             Call alloc_free( ALLOC_INTEGER, Size( ipiv ) )
             Deallocate( ipiv )
             
          Case( MATRIX_COMPLEX )
             Allocate( ipiv( 1:Min( a%n, a%m ) ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_INTEGER, Size( ipiv ) )
             Call zgetrf( a%n, a%m, a%complex_data, a%local_ld, ipiv, error )
             If( error > 0 ) Then
                Call dlc_error( 'Internal error: Matrix singular in ' // routine, 'abort' )
             End If
             Allocate( cwork( 1:a%n ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_COMPLEX, Size( cwork ) )
             Call zgetri( a%n, a%complex_data, a%local_ld, ipiv, cwork, Size( cwork ), error )
             If( error > 0 ) Then
                Call dlc_error( 'Internal error: Matrix singular in ' // routine, 'abort' )
             End If
             Call alloc_free( ALLOC_COMPLEX, Size( cwork ) )
             Deallocate( cwork )
             Call alloc_free( ALLOC_INTEGER, Size( ipiv ) )
             Deallocate( ipiv )

       End Select

    Else

       Select Case( a%genus )

          Case( MATRIX_REAL )
             Allocate( ipiv( 1:a%n + a%row%blocking_factor ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_INTEGER, Size( ipiv ) )

             Call pdgetrf( a%n, a%m, a%real_data, 1, 1, a%distribution%descriptor, ipiv, error )
             If( error > 0 ) Then
                Call dlc_error( 'Internal error: Matrix singular in ' // routine, 'abort' )
             End If

             Call pdgetri( a%n, a%real_data, 1, 1, a%distribution%descriptor, ipiv, rlwork, -1, &
                  liwork, -1, error )
             lwork = Nint( rlwork )
             Allocate( rwork( 1:lwork ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( rwork ) )
             Allocate( iwork( 1:liwork ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
             End If
             Call alloc_add( ALLOC_INTEGER, Size( iwork ) )

             Call pdgetri( a%n, a%real_data, 1, 1, a%distribution%descriptor, ipiv, rwork, lwork, &
                  iwork, liwork, error )

             Call alloc_free( ALLOC_INTEGER, Size( iwork ) )
             Deallocate( iwork )
             Call alloc_free( ALLOC_REAL, Size( rwork ) )
             Deallocate( rwork )
             Call alloc_free( ALLOC_INTEGER, Size( ipiv ) )
             Deallocate( ipiv )

          Case( MATRIX_COMPLEX )

             Call dlc_error( 'Internal error: complex case not implemented in ' // routine, 'abort' )

       End Select

    End If

  End Subroutine matrix_invert_one

  Subroutine matrix_invert_multi( a )

    Type( matrix ), Dimension( : ) :: a

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_invert( a( i ) )
       End If
    End Do

  End Subroutine matrix_invert_multi

  Subroutine matrix_absmaxloc_col_one( a, big_loc )

    Include 'mpif.h'

    Type( matrix )                         :: a
    Integer, Dimension( : ), Intent( Out ) :: big_loc

    ! It's done this way to remind me of potential
    ! pitfalls when parallelizing - what if two procs
    ! have identical values ?
    Type biggest
       Sequence
       Real( wp ) :: value
       Integer    :: index
       Integer    :: who_owns
    End Type biggest

    Type( biggest ), Dimension( : ), Allocatable :: biggest_value

    Real( wp ), Dimension( : ), Allocatable :: biggest_value_all

    Integer, Dimension( : ), Allocatable :: big_loc_local

    Integer :: lg
    Integer :: error
    Integer :: i, j

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_MAXLOC_PRINT' 

    Call matrix_check_valid( routine, a )

    Allocate( biggest_value( 1:a%m ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( biggest_value ) )
    Allocate( biggest_value_all( 1:a%m ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( biggest_value_all ) )
    Allocate( big_loc_local( 1:a%m ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( big_loc_local ) )

    biggest_value%value = 0.0_wp
    biggest_value%index = -1

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = 1, a%local_m
             lg = a%col%local_to_global( i )
             biggest_value( lg )%value    = Abs( a%real_data( 1, i ) )
             biggest_value( lg )%index    = a%row%local_to_global( 1 )
             biggest_value( lg )%who_owns = a%distribution%me
             Do j = 2, a%local_n
                If( Abs( a%real_data( j, i ) ) > biggest_value( lg )%value  ) Then
                   biggest_value( lg )%value    = Abs( a%real_data( j, i ) )
                   biggest_value( lg )%index    = a%row%local_to_global( j )
                   biggest_value( lg )%who_owns = a%distribution%me
                End If
             End Do
          End Do

       Case( MATRIX_COMPLEX )
          Do i = 1, a%local_m
             lg = a%col%local_to_global( i )
             biggest_value( lg )%value    = Abs( a%complex_data( 1, i ) )
             biggest_value( lg )%index    = a%row%local_to_global( 1 )
             biggest_value( lg )%who_owns = a%distribution%me
             Do j = 2, a%local_n
                If( Abs( a%complex_data( j, i ) ) > Abs( biggest_value( lg )%value ) ) Then
                   biggest_value( lg )%value    = Abs( a%complex_data( j, i ) )
                   biggest_value( lg )%index    = a%row%local_to_global( j )
                   biggest_value( lg )%who_owns = a%distribution%me
                End If
             End Do
          End Do

    End Select

    ! This is potentially bugy if two different procs have the same value
    ! in the same col - however keep it simple for the moment
    ! Could also improve effciency by doing reduces over columns
    ! of proc grid, but again KISS.

    If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then

       Call mpi_allreduce( biggest_value%value, biggest_value_all, &
            Size( biggest_value_all ), mpi_double_precision, mpi_max, &
            a%distribution%communicator, error )

       Do i = 1, a%m
          If( biggest_value_all( i ) == biggest_value( i )%value ) Then
             big_loc_local( i ) = biggest_value( i )%index
          Else
             big_loc_local( i ) = -1
          End If
       End Do       

       Call mpi_allreduce( big_loc_local, big_loc, &
            Size( big_loc ), mpi_integer, mpi_max, &
            a%distribution%communicator, error )

    Else
       
       big_loc = biggest_value%index

    End If

    Call alloc_free( ALLOC_INTEGER, Size( big_loc_local ) )
    Deallocate( big_loc_local )
    Call alloc_free( ALLOC_REAL, Size( biggest_value_all ) )
    Deallocate( biggest_value_all )
    Call alloc_free( ALLOC_DERIVED, Size( biggest_value ) )
    Deallocate( biggest_value )

  End Subroutine matrix_absmaxloc_col_one

  Subroutine matrix_absmaxloc_col_multi( a, big_loc )

    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( :    )                :: a
    Integer       , Dimension( :, : ), Intent( Out ) :: big_loc

    Integer :: error
    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_absmaxloc_col( a( i ), big_loc( :, i ) )
       End If
    End Do
    
    ! If UHF_K need to redistribute to all procs for all spins/k
    If( a( 1 )%distribution%type == MATRIX_DISTRIB_UHF_K ) Then
       Do i = 0, comms_data_uhf_k_repl_nproc - 1
          Call mpi_bcast( big_loc( :, i + 1 ), Size( big_loc, Dim = 1 ), &
               mpi_integer, i, comms_data_uhf_k_repl_comm, error )
       End Do
    End If

  End Subroutine matrix_absmaxloc_col_multi

  Subroutine matrix_check_abs_ov_one( a, occs, ov )

    Type( matrix )                           :: a
    Logical, Dimension( : ), Intent( In    ) :: occs
    Real( wp )             , Intent(   Out ) :: ov

    Integer :: i, j

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_CHECK_ABS_OV' 

    Call matrix_check_valid( routine, a )

    ov = -1.0_wp

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = 1, a%local_m
             If( occs( a%col%local_to_global( i ) ) ) Then
                Do j = 1, a%local_n
                   If( .Not. occs( a%row%local_to_global( j ) ) ) Then
                      ov = Max( ov, Abs( a%real_data( j, i ) ) )
                   End If
                End Do
             End If
          End Do

       Case( MATRIX_COMPLEX )
          Do i = 1, a%local_m
             If( occs( a%col%local_to_global( i ) ) ) Then
                Do j = 1, a%local_n
                   If( .Not. occs( a%row%local_to_global( j ) ) ) Then
                      ov = Max( ov, Abs( a%complex_data( j, i ) ) )
                   End If
                End Do
             End If
          End Do

    End Select

  End Subroutine matrix_check_abs_ov_one

  Subroutine matrix_check_abs_ov_multi( a, occs, result )

    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( :    )                  :: a
    Logical       , Dimension( :, : ), Intent( In    ) :: occs
    Real( wp )    , Dimension( :    ), Intent(   Out ) :: result

    Character( Len = 25 ), Parameter :: routine = 'MATRIX_CHECK_ABS_OV_MULTI'

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer :: error
    Integer :: i

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error in ' // routine // ' : failed to alloc memory ', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work ) )

    work = -1.0_wp

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_check_abs_ov( a( i ), occs( :, i ), work( i ) )
       End If
    End Do

    ! If UHF_K need to replicate for all spins/k
    If( a( 1 )%distribution%type == MATRIX_DISTRIB_REPL ) Then
       result = work
    Else
       Call mpi_allreduce( work, result, Size( work ), mpi_double_precision, &
            mpi_max, comms_data_all_comm, error )
    End If

    Call alloc_free( ALLOC_REAL, Size( work ) )
    Deallocate( work )

  End Subroutine matrix_check_abs_ov_multi

  Subroutine matrix_orthogonalize_one( a, metric, have_metric )
    
    !TMP HACK!!
    Use comms_data

    Include 'mpif.h'

    ! No complex case

    Type( matrix )        :: a
    Type( matrix )        :: metric
    Logical, Intent( In ) :: have_metric

    Logical, External :: opg_root

    Type( matrix ) :: metric_prod

    Real( wp ), Dimension( :, : ), Allocatable :: work1_par
    Real( wp ), Dimension( :, : ), Allocatable :: work2_par

    Real( wp ), Dimension( : ), Allocatable :: work1
    Real( wp ), Dimension( : ), Allocatable :: work2

    Real( wp ), Dimension( 1:2 ) :: tmp
    Real( wp ), Dimension( 1:2 ) :: tmp2

    Real( wp ) :: si, scale
    Real( wp ) :: sum_local
    Real( wp ) :: sum_global
    Real( wp ) :: work_i

    Integer, Parameter :: max_passes = 3

    Integer, Dimension( : ), Allocatable :: who_holds_i

    Integer :: n, m
    Integer :: npass, npass_tot
    Integer :: gl
    Integer :: error
    Integer :: i, j

    Character( Len = 20 ), Parameter :: routine = 'MATRIX_ORTHOGONALIZE'

    Call matrix_check_valid( routine, a      )
    Call matrix_check_valid( routine, metric )

!!$    Call matrix_check_conform( routine, a, metric )

    n = a%n
    m = a%m

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then
       Allocate( work1( 1:a%local_ld ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( work1 ) )
       Allocate( work2( 1:a%local_ld ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( work2 ) )
    Else
       Allocate( work1_par( 1:a%local_ld, 1:a%local_sd ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( work1_par ) )
       Allocate( work2_par( 1:a%local_ld, 1:a%local_sd ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in ' // routine, 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( work2_par ) )
       Call matrix_create( a, metric_prod )
       If( have_metric ) Then
          Call matrix_multiply( 1.0_wp, metric, a, 0.0_wp, metric_prod )
       End If
       work1_par = 0.0_wp
    End If

    npass_tot = 0
    Do i = 1, m

       gl = a%col%global_to_local( i )

       npass = 0
       Do While( npass < max_passes )
          npass = npass + 1
          
          If( have_metric ) Then
             If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then
                Call dgemv( 'n', n, n, 1.0_wp, metric%real_data   , n, &
                                               a%real_data( :, i ), 1, &
                                       0.0_wp, work1              , 1 )
             Else
                If( npass /= 1 ) Then
                   Call pdgemv( 'n', n, n, &
                        1.0_wp, metric%real_data, 1, 1, metric%distribution%descriptor, &
                        a%real_data, 1, i, a%distribution%descriptor, 1, &
                        0.0_wp, work1_par, 1, i, a%distribution%descriptor, 1 )
                Else
                   If( gl /= MATRIX_NOT_OWN ) Then
                      work1_par( 1:metric_prod%local_n, gl ) = metric_prod%real_data( 1:metric_prod%local_n, gl )
                   End If
                End If
             End If
          Else
             ! BUG HERE in parallel - not used but need to fix.
             work1 = a%real_data( 1:a%local_ld, gl )
          End If

          If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then
             Call dgemv( 't', n, i, 1.0_wp, a%real_data, n, &
                                            work1      , 1, &
                                    0.0_wp, work2      , 1 )

             If( i /= 1 ) Then
                Call dgemv( 'n', n, i - 1, -1.0_wp, a%real_data,         n, &
                                                    work2      ,         1, &
                                            1.0_wp, a%real_data( :, i ), 1 )
             End If
             si = work2( i ) - Sum( work2( 1:i - 1 ) * work2( 1:i - 1 ) )
             work_i = work2( i )
          Else
             Call pdgemv( 't', n, i, &
                  1.0_wp, a%real_data, 1, 1, a%distribution%descriptor, &
                  work1_par, 1, i, a%distribution%descriptor, 1, &
                  0.0_wp, work2_par, 1, i, a%distribution%descriptor, 1 )

             If( i /= 1 ) Then
                Call pdgemv( 'n', n, i - 1, &
                     -1.0_wp, a%real_data, 1, 1, a%distribution%descriptor, &
                     work2_par, 1, i, a%distribution%descriptor, 1, &
                     1.0_wp, a%real_data, 1, i, a%distribution%descriptor, 1 )
             End If
             If( gl /= MATRIX_NOT_OWN ) Then
                sum_local = 0.0_wp
                Do j = 1, a%local_n
                   If( a%row%local_to_global( j ) >= i ) Then
                      Exit
                   End If
                   sum_local = sum_local + work2_par( j, gl  ) ** 2
                End Do
                If( a%row%who_owns( i ) == a%distribution%my_proc_row ) Then
                   work_i = work2_par( a%row%global_to_local( i ), gl )
                Else
                   work_i = 0.0_wp
                End If
                tmp( 1 ) = work_i
                tmp( 2 ) = sum_local
                Call mpi_allreduce( tmp, tmp2, 2, mpi_double_precision, &
                     mpi_sum, a%distribution%col_comm, error )
                work_i     = tmp2( 1 )
                sum_global = tmp2( 2 )
                si         = work_i - sum_global
             Else
                si     = 0.0_wp
                work_i = 0.0_wp
             End If
          End If
             

          If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
             tmp( 1 ) = si
             tmp( 2 ) = work_i
             Call mpi_allreduce( tmp, tmp2, 2, mpi_double_precision, &
                  mpi_sum, a%distribution%row_comm, error )
             si     = tmp2( 1 )
             work_i = tmp2( 2 )
          End If

          If( i > 1 ) Then
             scale = si / work_i
             If( scale < 0.9_wp ) Then
                If( npass < max_passes ) Then
                   Cycle
                Else
                   Call dlc_error( routine // ': Failed to orthog vector', 'abort' )
                End If
             End If
          End If
             
          If( si < 0.0_wp ) Then
             Call dlc_error( routine // ': negative', 'abort' )
          Else If( si == 0.0_wp ) Then
             Call dlc_error( routine // ': Hard zero', 'abort' )
          Else
             scale = 1.0_wp / Sqrt( si )
          End If
             
          If( gl /= MATRIX_NOT_OWN ) Then
             a%real_data( 1:a%local_n, gl ) = a%real_data( 1:a%local_n, gl ) * scale
          End If

          npass_tot = npass_tot + npass

          Exit

       End Do

    End Do

!!$    If( a%distribution%me == 0 ) Then
!!$       Write( *, * ) 'Average number of passes in orthog: ', Real( npass_tot ) / m
!!$    End If

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then
       Call alloc_free( ALLOC_REAL, Size( work2 ) )
       Deallocate( work2 )
       Call alloc_free( ALLOC_REAL, Size( work1 ) )
       Deallocate( work1 )
    Else
       Call matrix_destroy( metric_prod )
       Call alloc_free( ALLOC_REAL, Size( work2_par ) )
       Deallocate( work2_par )
       Call alloc_free( ALLOC_REAL, Size( work1_par ) )
       Deallocate( work1_par )
    End If

  End Subroutine matrix_orthogonalize_one

  Subroutine matrix_orthogonalize_multi( a, metric, have_metric )
    
    Include 'mpif.h'

    ! 1) No complex case
    ! 2) No harmonc available

    Type( matrix ), Dimension( : )               :: a
    Type( matrix ), Dimension( : )               :: metric
    Logical                       , Intent( In ) :: have_metric

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_orthogonalize( a( i ), metric( i ), have_metric )
       End If
    End Do

  End Subroutine matrix_orthogonalize_multi

  Subroutine matrix_zero_patch_one( a, rows, cols )

    Type( matrix )                        :: a
    Logical, Dimension( : ), Intent( In ) :: rows
    Logical, Dimension( : ), Intent( In ) :: cols

    Character( Len = 21 ), Parameter :: routine = 'MATRIX_ZERO_PATCH_ONE'

    Integer :: i, j

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = 1, a%local_m
             If( cols( a%col%local_to_global( i ) ) ) Then
                Do j = 1, a%local_n
                   If( rows( a%row%local_to_global( j ) ) ) Then
                      a%real_data( j, i ) = 0.0_wp
                   End If
                End Do
             End If
          End Do

       Case( MATRIX_COMPLEX )
          Do i = 1, a%local_m
             If( cols( a%col%local_to_global( i ) ) ) Then
                Do j = 1, a%local_n
                   If( rows( a%row%local_to_global( j ) ) ) Then
                      a%complex_data( j, i ) = 0.0_wp
                   End If
                End Do
             End If
          End Do

    End Select

  End Subroutine matrix_zero_patch_one

  Subroutine matrix_zero_patch_multi( a, rows, cols )

    Type( matrix ), Dimension( :    )               :: a
    Logical       , Dimension( :, : ), Intent( In ) :: rows
    Logical       , Dimension( :, : ), Intent( In ) :: cols

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_zero_patch( a( i ), rows( :, i ), cols( :, i ) )
       End If
    End Do

  End Subroutine matrix_zero_patch_multi

  Subroutine matrix_set_patch_one( a, values, first_row, first_col )

    Type( matrix )                              :: a
    Real( wp ), Dimension( :, : ), Intent( In ) :: values
    Integer,                       Intent( In ) :: first_row
    Integer,                       Intent( In ) :: first_col

    Character( Len = 20 ), Parameter :: routine = 'MATRIX_SET_PATCH_ONE'

    Integer :: local_i, local_j
    Integer :: i, j

    Call matrix_check_valid( routine, a )

    Select Case( a%genus )

       Case( MATRIX_REAL )
          Do i = first_col, first_col + Size( values, Dim = 2 ) - 1
             local_i = a%col%global_to_local( i )
             If( local_i /= MATRIX_NOT_OWN ) Then
                Do j = first_row, first_row + Size( values, Dim = 1 ) - 1
                   local_j = a%row%global_to_local( j )
                   If( local_j /= MATRIX_NOT_OWN ) Then
                      a%real_data( local_j, local_i ) = &
                           values( j - first_row + 1, i - first_col + 1 ) 
                   End If
                End Do
             End If
          End Do

       Case( MATRIX_COMPLEX )
          Do i = first_col, first_col + Size( values, Dim = 2 ) - 1
             local_i = a%col%global_to_local( i )
             If( local_i /= MATRIX_NOT_OWN ) Then
                Do j = first_row, first_row + Size( values, Dim = 1 ) - 1
                   local_j = a%row%global_to_local( j )
                   If( local_j /= MATRIX_NOT_OWN ) Then
                      a%complex_data( local_j, local_i ) = &
                           values( j - first_row + 1, i - first_col + 1 ) 
                   End If
                End Do
             End If
          End Do

    End Select

  End Subroutine matrix_set_patch_one

  Subroutine matrix_set_patch_multi( a, values, first_row, first_col )

    Type( matrix ), Dimension( :       )               :: a
    Real( wp )    , Dimension( :, :, : ), Intent( In ) :: values
    Integer,        Dimension( :       ), Intent( In ) :: first_row
    Integer,        Dimension( :       ), Intent( In ) :: first_col

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_set_patch( a( i ), values( :, :, i ), first_row( i ), first_col( i ) )
       End If
    End Do

  End Subroutine matrix_set_patch_multi

  Subroutine matrix_print( a )

    Type( matrix ) :: a

    Character( Len = 12 ), Parameter :: routine = 'MATRIX_PRINT' 

    Logical, External :: opg_root

    Real( wp ) :: aij

    Integer :: start, finish
    Integer :: i, j

    Call matrix_check_valid( routine, a )

    If( a%distribution%type == MATRIX_DISTRIB_REPL ) Then

       If( opg_root() ) Then

          Write( *, * )
          Write( *, '( a, a ) ' ) 'Contents of matrix: ', a%name
          Write( *, '( a )' )     '---------------------------------------'
    
          Select Case( a%genus )

          Case( MATRIX_REAL )
             call prmat( a%real_data, a%local_n, a%local_m, .False. )

          Case( MATRIX_COMPLEX )
             call prmat( a%complex_data, a%local_n, a%local_m, .False. )

          End Select

       End If

    Else

       finish = 0

       Do While( finish < a%m )
          start = finish + 1
          finish = Min( a%n, finish + 12 )
          If( a%distribution%me == 0 ) Then
             Write( *, * )
             Write( *, '( 6x, 12( 3x, i3, 3x ) )' ) ( i, i = start, finish )
          End If
          Do j = 1, a%n
             If( a%distribution%me == 0 ) Then
                Write( *, '( i5, 1x )', Advance = 'No' ) j
             End If
             Do i = start, finish
                Call pdelget( 'A', ' ', aij, a%real_data, j, i, a%distribution%descriptor )
                If( a%distribution%me == 0 ) Then
                   Write( *, '( f9.5 )', Advance = 'No' ) aij
                End If
             End Do
             If( a%distribution%me == 0 ) Then
                Write( *, * )
             End If
          End Do
       End Do
    End If

  End Subroutine matrix_print

  Subroutine matrix_print_titled( a, title )

    Type( matrix )                     :: a
    Character( Len = * ), Intent( In ) :: title

    Character( Len = 19 ), Parameter :: routine = 'MATRIX_PRINT_TITLED'

    Integer :: i

    Call matrix_check_valid( routine, a )

    Write( 6, * )
    Write( 6, '( a )' ) title 
    Write( 6, '( 120a1 )' ) ( '-', i, i = 1, Len_trim( title ) )
    
    Select Case( a%genus )

       Case( MATRIX_REAL )
          call prmat( a%real_data, a%local_n, a%local_m, .False. )

       Case( MATRIX_COMPLEX )
          call prmat( a%complex_data, a%local_n, a%local_m, .False. )

    End Select

  End Subroutine matrix_print_titled

  Subroutine matrix_print_to_file( a, unit )

    Type( matrix )        :: a
    Integer, Intent( In ) :: unit

    Character( Len = 12 ), Parameter :: routine = 'MATRIX_PRINT' 

    If( a%i_own ) Then

    Call matrix_check_valid( routine, a )

    write( unit, * )
    write( unit, '( a, a ) ' ) 'Contents of matrix: ', a%name
    write( unit, '( a )' )     '---------------------------------------'
    
    Select Case( a%genus )

       Case( MATRIX_REAL )
          call prmat( a%real_data, a%local_n, a%local_m, .False., unit )

       Case( MATRIX_COMPLEX )
          call prmat( a%complex_data, a%local_n, a%local_m, .False., unit )

    End Select

    End If

  End Subroutine matrix_print_to_file

  Subroutine prmat_real( a, m, n, high_precision, unit )

    Real( wp ), Dimension( :, : ), Intent( In )           :: a
    Integer                      , Intent( In )           :: m
    Integer                      , Intent( In )           :: n
    Logical                      , Intent( In )           :: high_precision
    Integer                      , Intent( In ), Optional :: unit

    Integer :: this_unit
    Integer :: start, finish, step
    Integer :: i, j

    If( Present( unit ) ) Then
       this_unit = unit
    Else
       this_unit = 6
    End If

    If( high_precision ) Then
       step = 7
    Else
       step = 12
    End If

    finish = 0

    Do While( finish < n )
       start  = finish + 1
       finish = Min( n, finish + step )
       Write( this_unit, * )
       If( high_precision ) Then
          Write( this_unit, '( 6x,  7( 6x, i3, 6x ) )' ) ( i, i = start, finish )
       Else
          Write( this_unit, '( 6x, 12( 3x, i3, 3x ) )' ) ( i, i = start, finish )
       End If
       Do j = 1, m
          If( high_precision ) Then
             Write( this_unit, '( i5, 1x, 7f15.10 )' ) j, ( a( j, i ), i = start, finish )
          Else
             Write( this_unit, '( i5, 1x, 12f9.5  )' ) j, ( a( j, i ), i = start, finish )
          End If
       End Do
    End Do

  End Subroutine prmat_real
          
  Subroutine prmat_complex( a, m, n, high_precision, unit )

    Complex( wp ), Dimension( :, : ), Intent( In )           :: a
    Integer                         , Intent( In )           :: m
    Integer                         , Intent( In )           :: n
    Logical                         , Intent( In )           :: high_precision
    Integer                         , Intent( In ), Optional :: unit

    Integer :: this_unit
    Integer :: start, finish, step
    Integer :: i, j

    If( Present( unit ) ) Then
       this_unit = unit
    Else
       this_unit = 6
    End If

    If( high_precision ) Then
       step = 3
    Else
       step = 6
    End If

    finish = 0

    Do While( finish < n )
       start  = finish + 1
       finish = Min( n, finish + step )
       Write( this_unit, * )
       If( high_precision ) Then
          Write( this_unit, '( 6x,  7( 6x, i3, 6x ) )' ) ( i, i = start, finish )
       Else
          Write( this_unit, '( 6x, 12( 3x, i3, 3x ) )' ) ( i, i = start, finish )
       End If
       Do j = 1, m
          If( high_precision ) Then
             Write( this_unit, '( i5, 1x, 3( "(", f15.10, ",", f15.10 ) )' ) j, ( a( j, i ), i = start, finish )
          Else
             Write( this_unit, '( i5, 1x, 6( "(", f9.5  , ",", f9.5   ) )' ) j, ( a( j, i ), i = start, finish )
          End If
       End Do
    End Do

  End Subroutine prmat_complex
          
  Subroutine matrix_read_rectangle( a, where, idaf )

    Type( matrix )        :: a
    Integer, Intent( In ) :: where
    Integer, Intent( In ) :: idaf

    Integer :: left
    Integer :: length, blocks, location
    Integer :: col_i, row_j
    Integer :: n, m
    Integer :: i, j, k
    
    Call matrix_inquire( a, global_n = n )
    Call matrix_inquire( a, global_m = m )

    left     = n * m
    location = where

    If( a%i_own ) Then
       a%real_data = 0.0_wp
    End If

    i = 1
    j = 1

    Do While( left > 0 )

       length = Min( left, n_io_buff )
       blocks = length / io_block

       Call rdedx_prec( buff, length, location, idaf )

       If( a%i_own ) Then

          Do k = 1, length

             row_j = a%row%global_to_local( j )
             col_i = a%col%global_to_local( i )

             If( row_j /= MATRIX_NOT_OWN .And. col_i /= MATRIX_NOT_OWN ) Then
                a%real_data( row_j, col_i ) = buff( k )
             End If

             j = j + 1
             If( j > n ) Then
                j = 1
                i = i + 1
             End If

          End Do

       End If

       left     = left - length
       location = location + blocks

    End Do

  End Subroutine matrix_read_rectangle

  Subroutine matrix_read_triangle( a, where, idaf )

    Type( matrix )        :: a
    Integer, Intent( In ) :: where
    Integer, Intent( In ) :: idaf

    Integer :: left
    Integer :: length, blocks, location
    Integer :: col_i, row_j
    Integer :: row_i, col_j
    Integer :: n
    Integer :: i, j, k
    
    Call matrix_inquire( a, global_n = n )

    left     = ( n * ( n + 1 ) ) / 2
    location = where

    If( a%i_own ) Then
       a%real_data = 0.0_wp
    End If

    i = 1
    j = 1

    Do While( left > 0 )

       length = Min( left, n_io_buff )
       blocks = length / io_block

       Call rdedx_prec( buff, length, location, idaf )

       If( a%i_own ) Then

          Do k = 1, length

             row_j = a%row%global_to_local( j )
             col_i = a%col%global_to_local( i )

             If( row_j /= MATRIX_NOT_OWN .And. col_i /= MATRIX_NOT_OWN ) Then
                a%real_data( row_j, col_i ) = buff( k )
             End If

             col_j = a%col%global_to_local( j )
             row_i = a%row%global_to_local( i )

             If( col_j /= MATRIX_NOT_OWN .And. row_i /= MATRIX_NOT_OWN ) Then
                a%real_data( row_i, col_j ) = buff( k )
             End If

             j = j + 1
             If( j > i ) Then
                j = 1
                i = i + 1
             End If

          End Do

       End If

       left     = left - length
       location = location + blocks

    End Do

  End Subroutine matrix_read_triangle

  Subroutine matrix_write_rectangle( a, where, idaf )

    Use comms_data

    ! Needs finishing off.

    Type( matrix )        :: a
    Integer, Intent( In ) :: where
    Integer, Intent( In ) :: idaf

    Integer, External :: lensec

    Integer :: left
    Integer :: length, blocks, location
    Integer :: col_i, row_j
    Integer :: n, m
    Integer :: i, j, k
    
    Call matrix_inquire( a, global_n = n )
    Call matrix_inquire( a, global_m = m )

    left     = n * m
    location = where

    i = 1
    j = 1

    Do While( left > 0 )

       length = Min( left, n_io_buff )
       blocks = lensec( length )

       buff_write = 0.0_wp

       If( a%i_own ) Then

          Do k = 1, length

             row_j = a%row%global_to_local( j )
             col_i = a%col%global_to_local( i )

             If( row_j /= MATRIX_NOT_OWN .And. col_i /= MATRIX_NOT_OWN ) Then
                buff_write( k ) = a%real_data( row_j, col_i )
             End If

             j = j + 1
             If( j > n ) Then
                j = 1
                i = i + 1
             End If

          End Do

       End If

       If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
          Call replicate_array( length, comms_data_all_comm, buff_write )
       End If

       Call wrt3( buff_write, length, location, idaf )

       left     = left - length
       location = location + blocks

    End Do
       
  End Subroutine matrix_write_rectangle

  Subroutine matrix_write_triangle( a, where, idaf )

    Use comms_data

    ! Needs finishing off.

    Type( matrix )        :: a
    Integer, Intent( In ) :: where
    Integer, Intent( In ) :: idaf

    Integer :: left
    Integer :: length, blocks, location
    Integer :: col_i, row_j
    Integer :: n
    Integer :: i, j, k
    Logical, External :: opg_root
    
    Call matrix_inquire( a, global_n = n )

    left     = ( n * ( n + 1 ) ) / 2
    location = where

    i = 1
    j = 1

    Do While( left > 0 )

       length = Min( left, n_io_buff )
       blocks = length / io_block

       buff_write = 0.0_wp

       If( a%i_own ) Then

          Do k = 1, length

             row_j = a%row%global_to_local( j )
             col_i = a%col%global_to_local( i )

             If( row_j /= MATRIX_NOT_OWN .And. col_i /= MATRIX_NOT_OWN ) Then
                buff_write( k ) = a%real_data( row_j, col_i )
             End If

             j = j + 1
             If( j > i ) Then
                j = 1
                i = i + 1
             End If

          End Do

       End If

       If( a%distribution%type /= MATRIX_DISTRIB_REPL ) Then
          Call replicate_array( length, comms_data_all_comm, buff_write )
       End If

       Call wrt3( buff_write, length, location, idaf )

       left     = left - length
       location = location + blocks

    End Do

  End Subroutine matrix_write_triangle

  Subroutine matrix_check_valid( routine, a )

    Character( Len = * ), Intent( In ) :: routine
    Type( matrix )                     :: a

    If( .Not. a%inited ) Then
       Call dlc_error( 'Uninitialized matrix in ' // routine, 'abort' )
    End If

    If( .Not. a%i_own ) Then
       Call dlc_error( 'Unowned matrix in ' // routine, 'abort' )
    End If

  End Subroutine matrix_check_valid
    
  Subroutine matrix_check_conform( routine, a, b, check_distrib )

    Character( Len = * ), Intent( In )           :: routine
    Type( matrix )                               :: a
    Type( matrix )                               :: b
    Logical             , Intent( In ), Optional :: check_distrib

    Logical :: this_check_distrib
    
    If( Present( check_distrib ) ) Then
       this_check_distrib = check_distrib
    Else
       this_check_distrib = .True.
    End If

    If( a%n /= b%n ) Then
       Call dlc_error( 'Incompatible matrix dimensions in ' // routine, 'abort' )
    End If

    If( a%m /= b%m ) Then
       Call dlc_error( 'Incompatible matrix dimensions in ' // routine, 'abort' )
    End If

    If( a%genus /= b%genus ) Then
       Call dlc_error( 'Incompatible matrix types in ' // routine, 'abort' )
    End If
    
    If( this_check_distrib ) Then
       If( .Not. ( a%distribution == b%distribution ) ) Then
          Call dlc_error( 'Incompatible matrix distributions in ' // routine, 'abort' )
       End If
    End If

  End Subroutine matrix_check_conform

  Logical Function matrix_compare_distrib( a, b )

    Type( matrix_distribution ), Intent( In ) :: a
    Type( matrix_distribution ), Intent( In ) :: b

    Logical :: result

    result = a%communicator == b%communicator

    result = result .And. a%context      == b%context
    result = result .And. All( a%descriptor == b%descriptor )

    matrix_compare_distrib = result

  End Function matrix_compare_distrib

  Character Function matrix_to_lower( char )

    Character, Intent( In ) :: char

    Character( Len = 26 ), Parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character( Len = 26 ), Parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'

    Integer :: where

    If( char >= 'A' .And. char <= 'Z' ) Then
       where = Index( upper, char )
       matrix_to_lower = lower( where:where )
    Else
       matrix_to_lower = char
    End If

  End Function matrix_to_lower

  Subroutine matrix_replicate_uhf_k( a )

    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( : ) :: a

    Integer :: error
    Integer :: i

    Do i = 0, comms_data_uhf_k_repl_nproc - 1

       Select Case( a( i + 1 )%genus )

          Case( MATRIX_REAL )
             Call mpi_bcast( a( i + 1 )%real_data, Size( a( i + 1 )%real_data ), &
                  mpi_double_precision, i, comms_data_uhf_k_repl_comm, error )

          Case( MATRIX_COMPLEX )
             Call mpi_bcast( a( i + 1 )%complex_data, 2 * Size( a( i + 1 )%complex_data ), &
                  mpi_double_precision, i, comms_data_uhf_k_repl_comm, error )

       End Select

    End Do

  End Subroutine matrix_replicate_uhf_k

  Subroutine matrix_comms_init( n_basis, base_comm, base_context, nspin, nk )

    Use comms_data

    Include 'mpif.h'

    Integer, Intent( In ) :: n_basis
    Integer, Intent( In ) :: base_comm
    Integer, Intent( In ) :: base_context
    Integer, Intent( In ) :: nspin
    Integer, Intent( In ) :: nk

    Integer, External :: numroc

    Integer, Dimension( 1:9 ) :: desc_repl, desc_all, desc_uhf_k

    Integer :: block_repl
    Integer :: block_all
    Integer :: block_uhf_k
    Integer :: local_n_all
    Integer :: local_m_all
    Integer :: local_ld_all
    Integer :: local_sd_all
    Integer :: local_n_uhf_k
    Integer :: local_m_uhf_k
    Integer :: local_ld_uhf_k
    Integer :: local_sd_uhf_k
    Integer :: file
    Integer :: error
    Integer :: i

    ! Set up the various communicators/contexts etc.
    Call comms_data_setup( base_comm, base_context, nspin, nk )

    ! Set up the default blocking factors for the 3 distributions
    ! Need error condition here for too many procs for too small a case
    block_repl  = n_basis
    block_all   = Min( comms_data_block_default, n_basis / Max( &
         comms_data_all_nrow, comms_data_all_ncol ) )
    If( block_all == 0 ) Then
       Call dlc_error( 'Blocking factor dropped to zero in MATRIX_COMMS_INIT. ' // &
            'You are running too small a problem on too many processors.', 'abort' )
    End If
    block_uhf_k = Min( comms_data_block_default, n_basis / Max( &
         comms_data_uhf_k_nrow, comms_data_uhf_k_ncol ) )
    If( block_uhf_k == 0 ) Then
       Call dlc_error( 'Blocking factor dropped to zero in MATRIX_COMMS_INIT. ' // &
            'You are running too small a problem on too many processors.', 'abort' )
    End If

    ! Set up the descriptors
    ! First need to know how much local storage we will need for the array.
    local_n_all = numroc( n_basis, block_all, comms_data_all_my_row, 0, &
         comms_data_all_nrow )
    local_m_all = numroc( n_basis, block_all, comms_data_all_my_col, 0, &
         comms_data_all_ncol )

    local_n_uhf_k = numroc( n_basis, block_uhf_k, comms_data_uhf_k_my_row, 0, &
         comms_data_uhf_k_nrow )
    local_m_uhf_k = numroc( n_basis, block_uhf_k, comms_data_uhf_k_my_col, 0, &
         comms_data_uhf_k_ncol )

    ! All implementations I know require that each proc has the same leading dim
    ! on each proc, be extra safe and make same size on each proc.
    Call mpi_allreduce( local_n_all, local_ld_all, 1, mpi_integer, mpi_max, &
         comms_data_all_comm, error )
    Call mpi_allreduce( local_m_all, local_sd_all, 1, mpi_integer, mpi_max, &
         comms_data_all_comm, error )

    Call mpi_allreduce( local_n_uhf_k, local_ld_uhf_k, 1, mpi_integer, mpi_max, &
         comms_data_uhf_k_comm, error )
    Call mpi_allreduce( local_m_uhf_k, local_sd_uhf_k, 1, mpi_integer, mpi_max, &
         comms_data_uhf_k_comm, error )

    ! Set up the descriptors
    Call descinit( desc_repl, n_basis, n_basis, block_repl, block_repl, 0, 0, &
         comms_data_all_context, n_basis, error )
    Call descinit( desc_all, n_basis, n_basis, block_all, block_all, 0, 0, &
         comms_data_all_context, local_ld_all, error )
    Call descinit( desc_uhf_k, n_basis, n_basis, block_uhf_k, block_uhf_k, 0, 0, &
         comms_data_uhf_k_context, local_ld_uhf_k, error )

    ! Set up the distribution info for the 3 cases
    distrib_repl%type         = MATRIX_DISTRIB_REPL
    distrib_repl%communicator = comms_data_all_comm
    distrib_repl%row_comm     = comms_data_all_row_comm
    distrib_repl%col_comm     = comms_data_all_col_comm
    distrib_repl%n_proc       = 1
    distrib_repl%me           = 0
    distrib_repl%context      = comms_data_not_valid
    distrib_repl%n_proc_row   = 1
    distrib_repl%n_proc_col   = 1
    distrib_repl%my_proc_row  = 0
    distrib_repl%my_proc_col  = 0
    distrib_repl%descriptor   = desc_repl

    distrib_all%type         = MATRIX_DISTRIB_ALL
    distrib_all%communicator = comms_data_all_comm
    distrib_all%row_comm     = comms_data_all_row_comm
    distrib_all%col_comm     = comms_data_all_col_comm
    distrib_all%n_proc       = comms_data_all_nproc
    distrib_all%me           = comms_data_all_me
    distrib_all%context      = comms_data_all_context
    distrib_all%n_proc_row   = comms_data_all_nrow
    distrib_all%n_proc_col   = comms_data_all_ncol
    distrib_all%my_proc_row  = comms_data_all_my_row
    distrib_all%my_proc_col  = comms_data_all_my_col
    distrib_all%descriptor   = desc_all

    If( .Not. comms_data_uhf_k_is_all ) Then
       distrib_uhf_k%type         = MATRIX_DISTRIB_UHF_K
    Else
       distrib_uhf_k%type         = MATRIX_DISTRIB_ALL
    End If
    distrib_uhf_k%communicator = comms_data_uhf_k_comm
    distrib_uhf_k%row_comm     = comms_data_uhf_k_row_comm
    distrib_uhf_k%col_comm     = comms_data_uhf_k_col_comm
    distrib_uhf_k%n_proc       = comms_data_uhf_k_nproc
    distrib_uhf_k%me           = comms_data_uhf_k_me
    distrib_uhf_k%context      = comms_data_uhf_k_context
    distrib_uhf_k%n_proc_row   = comms_data_uhf_k_nrow
    distrib_uhf_k%n_proc_col   = comms_data_uhf_k_ncol
    distrib_uhf_k%my_proc_row  = comms_data_uhf_k_my_row
    distrib_uhf_k%my_proc_col  = comms_data_uhf_k_my_col
    distrib_uhf_k%descriptor   = desc_uhf_k

    ! Now set up the indexing arrays and descriptor
    ! First for the replicated case

    index_repl_row%n               = n_basis
    index_repl_row%blocking_factor = block_repl
    index_repl_row%local_n         = n_basis
    index_repl_row%local_ld        = n_basis
    Allocate( index_repl_row%global_to_local( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_repl_row%global_to_local ) )
    index_repl_row%global_to_local = (/ ( i, i = 1, n_basis ) /)

    Allocate( index_repl_row%local_to_global( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_repl_row%global_to_local ) )
    index_repl_row%local_to_global = (/ ( i, i = 1, n_basis ) /)

    Allocate( index_repl_row%who_owns( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_repl_row%who_owns ) )
    index_repl_row%who_owns = comms_data_all_me

    index_repl_col = index_repl_row

    ! Now for the distributed arrays
    ! Set up the basic indexing stuff
    index_all_row%n               = n_basis
    index_all_row%blocking_factor = block_all
    index_all_row%local_n         = local_n_all
    index_all_row%local_ld        = local_ld_all
    
    index_all_col%n               = n_basis
    index_all_col%blocking_factor = block_all
    index_all_col%local_n         = local_m_all
    index_all_col%local_ld        = local_sd_all
    
    index_uhf_k_row%n               = n_basis
    index_uhf_k_row%blocking_factor = block_uhf_k
    index_uhf_k_row%local_n         = local_n_uhf_k
    index_uhf_k_row%local_ld        = local_ld_uhf_k
     
    index_uhf_k_col%n               = n_basis
    index_uhf_k_col%blocking_factor = block_uhf_k
    index_uhf_k_col%local_n         = local_m_uhf_k
    index_uhf_k_col%local_ld        = local_sd_uhf_k

    ! Now need to set up the 'translation' arrays
    Allocate( index_all_row%global_to_local( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_row%global_to_local ) )

    Call set_g_to_l( index_all_row%blocking_factor, &
                     distrib_all%n_proc_row,        &
                     distrib_all%my_proc_row,       &
                     index_all_row%global_to_local )

    Allocate( index_all_col%global_to_local( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_col%global_to_local ) )
    Call set_g_to_l( index_all_col%blocking_factor, &
                     distrib_all%n_proc_col,        &
                     distrib_all%my_proc_col,       &
                     index_all_col%global_to_local )

    Allocate( index_all_row%local_to_global( 1:local_n_all ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_row%local_to_global ) )
    Call set_l_to_g( n_basis,                       &
                     index_all_row%blocking_factor, &
                     distrib_all%n_proc_row,        &
                     distrib_all%my_proc_row,       &
                     index_all_row%local_to_global )

    Allocate( index_all_col%local_to_global( 1:local_m_all ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_col%local_to_global ) )
    Call set_l_to_g( n_basis,                       &
                     index_all_col%blocking_factor, &
                     distrib_all%n_proc_col,        &
                     distrib_all%my_proc_col,       &
                     index_all_col%local_to_global )

    Allocate( index_uhf_k_row%global_to_local( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_row%global_to_local ) )
    Call set_g_to_l( index_uhf_k_row%blocking_factor, &
                     distrib_uhf_k%n_proc_row,        &
                     distrib_uhf_k%my_proc_row,       &
                     index_uhf_k_row%global_to_local )

    Allocate( index_uhf_k_col%global_to_local( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_col%global_to_local ) )
    Call set_g_to_l( index_uhf_k_col%blocking_factor, &
                     distrib_uhf_k%n_proc_col,        &
                     distrib_uhf_k%my_proc_col,       &
                     index_uhf_k_col%global_to_local )

    Allocate( index_uhf_k_row%local_to_global( 1:local_n_uhf_k ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_row%local_to_global ) )
    Call set_l_to_g( n_basis,                       &
                     index_uhf_k_row%blocking_factor, &
                     distrib_uhf_k%n_proc_row,        &
                     distrib_uhf_k%my_proc_row,       &
                     index_uhf_k_row%local_to_global )


    Allocate( index_uhf_k_col%local_to_global( 1:local_m_uhf_k ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_col%local_to_global ) )
    Call set_l_to_g( n_basis,                       &
                     index_uhf_k_col%blocking_factor, &
                     distrib_uhf_k%n_proc_col,        &
                     distrib_uhf_k%my_proc_col,       &
                     index_uhf_k_col%local_to_global )

    ! Finally set up the distributed matrix ownership arrays
    Allocate( index_all_row%who_owns( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_row%who_owns ) )
    Call set_who_owns( n_basis, index_all_row%blocking_factor, distrib_all%n_proc_row, &
         index_all_row%who_owns )

    Allocate( index_all_col%who_owns( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_all_col%who_owns ) )
    Call set_who_owns( n_basis, index_all_col%blocking_factor, distrib_all%n_proc_col, &
         index_all_col%who_owns )

    Allocate( index_uhf_k_row%who_owns( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_row%who_owns ) )
    Call set_who_owns( n_basis, index_uhf_k_row%blocking_factor, distrib_uhf_k%n_proc_row, &
         index_uhf_k_row%who_owns )

    Allocate( index_uhf_k_col%who_owns( 1:n_basis ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( index_uhf_k_col%who_owns ) )
    Call set_who_owns( n_basis, index_uhf_k_col%blocking_factor, distrib_uhf_k%n_proc_col, &
         index_uhf_k_col%who_owns )

  End Subroutine matrix_comms_init

  Subroutine set_g_to_l( block, nproc, me, g_to_l )
    
    Integer,                 Intent( In    ) :: block
    Integer,                 Intent( In    ) :: nproc
    Integer,                 Intent( In    ) :: me
    Integer, Dimension( : ), Intent(   Out ) :: g_to_l
    
    Integer :: n
    Integer :: start, base, skip
    Integer :: j
    Integer :: i
    
    n = Size( g_to_l )
    
    g_to_l = MATRIX_NOT_OWN
    
    start = me * block
    skip  = nproc * block
    
    base = start
    j    = 0
    
    Do While( base <= n )
       Do i = base + 1, Min( base + block, n )
          j = j + 1
          g_to_l( i ) = j
       End Do
       base = base + skip
    End Do
    
  End Subroutine set_g_to_l
  
  Subroutine set_l_to_g( n, block, nproc, me, l_to_g )
    
    Integer,                 Intent( In    ) :: n
    Integer,                 Intent( In    ) :: block
    Integer,                 Intent( In    ) :: nproc
    Integer,                 Intent( In    ) :: me
    Integer, Dimension( : ), Intent(   Out ) :: l_to_g
    
    Integer :: start, base, skip
    Integer :: j
    Integer :: i
    
    start = me * block
    skip  = nproc * block
    
    base = start
    j    = 0
    
    Do While( base <= n )
       Do i = base + 1, Min( base + block, n )
          j = j + 1
          l_to_g( j ) = i
       End Do
       base = base + skip
    End Do
    
  End Subroutine set_l_to_g
  
  Subroutine set_who_owns( n, block, nproc, who_owns )

    Integer,                 Intent( In    ) :: n
    Integer,                 Intent( In    ) :: block
    Integer,                 Intent( In    ) :: nproc
    Integer, Dimension( : ), Intent(   Out ) :: who_owns

    Integer, External :: indxg2p
    
    Integer :: i

    Do i = 1, n
       who_owns( i ) = indxg2p( i, block, 0, 0, nproc ) 
    End Do
    
  End Subroutine set_who_owns

  Subroutine matrix_comms_finalize

    Use comms_data

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_col%who_owns ) )
    Deallocate( index_uhf_k_col%who_owns )

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_row%who_owns ) )
    Deallocate( index_uhf_k_row%who_owns )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_col%who_owns ) )
    Deallocate( index_all_col%who_owns )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_row%who_owns ) )
    Deallocate( index_all_row%who_owns )

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_col%local_to_global ) )
    Deallocate( index_uhf_k_col%local_to_global )

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_row%local_to_global ) )
    Deallocate( index_uhf_k_row%local_to_global )

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_col%global_to_local ) )
    Deallocate( index_uhf_k_col%global_to_local )

    Call alloc_free( ALLOC_INTEGER, Size( index_uhf_k_row%global_to_local ) )
    Deallocate( index_uhf_k_row%global_to_local )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_col%local_to_global ) )
    Deallocate( index_all_col%local_to_global )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_row%local_to_global ) )
    Deallocate( index_all_row%local_to_global )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_col%global_to_local ) )
    Deallocate( index_all_col%global_to_local )

    Call alloc_free( ALLOC_INTEGER, Size( index_all_row%global_to_local ) )
    Deallocate( index_all_row%global_to_local )

    Call alloc_free( ALLOC_INTEGER, Size( index_repl_row%who_owns ) )
    Deallocate( index_repl_row%who_owns )

    Call alloc_free( ALLOC_INTEGER, Size( index_repl_row%local_to_global ) )
    Deallocate( index_repl_row%local_to_global )

    Call alloc_free( ALLOC_INTEGER, Size( index_repl_row%global_to_local ) )
    Deallocate( index_repl_row%global_to_local )

    Call comms_data_free

  End Subroutine matrix_comms_finalize

  Subroutine replicate_array( n, comm, a )

    Include 'mpif.h'

    Integer                   , Intent( In    ) :: n
    Integer                   , Intent( In    ) :: comm
    Real( wp ), Dimension( : ), Intent( InOut ) :: a

    Integer :: start, finish, length
    Integer :: error

    start = 1

    Do While( start <= n )
       finish = Min( n, start + n_repl_buff - 1 )
       length = finish - start + 1
       buff( 1:length ) = a( start:finish )
       Call mpi_allreduce( buff, a( start ), length, mpi_double_precision, &
            mpi_sum, comm, error )
       start = finish + 1
    End Do

  End Subroutine replicate_array

  Subroutine matrix_tdown_setup( a, ctran, ntran, itran, ilifc, otran )

    Type( matrix )                               :: a
    Real( wp )    , Dimension( : ), Intent( In ) :: ctran
    Integer       , Dimension( : ), Intent( In ) :: ntran
    Integer       , Dimension( : ), Intent( In ) :: itran
    Integer       , Dimension( : ), Intent( In ) :: ilifc
    Logical                       , Intent( In ) :: otran

    Integer, External :: indxg2l
    Integer, External :: indxg2p

    Integer :: proc_from
    Integer :: proc_shift
    Integer :: n_all_ops_local, n_this_ops_local
    Integer :: global_j, local_j
    Integer :: error
    Integer :: i, j, k, l, m, n

    If( .Not. otran ) Then

       ! First of all work out which of the possible ops in tdown
       ! will act on local data for local processing. First count
       ! them
       n_all_ops_local = 0
       Do j = 1, a%n
          n = ntran( j )
          Do k = 1, n
             l = ilifc( j ) + k
             If( a%row%global_to_local( itran( l ) ) /= MATRIX_NOT_OWN ) Then
                n_all_ops_local = n_all_ops_local + 1
             End If
          End Do
       End Do

       Allocate( all_ops%j_vec_add ( 1:n_all_ops_local ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
       End If
       Call alloc_add( ALLOC_INTEGER, Size(all_ops%j_vec_add  ) )
       Allocate( all_ops%j_tran    ( 1:n_all_ops_local ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
       End If
       Call alloc_add( ALLOC_INTEGER, Size( all_ops%j_tran ) )
       Allocate( all_ops%j_vec_base( 1:n_all_ops_local ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
       End If
       Call alloc_add( ALLOC_INTEGER, Size( all_ops%j_vec_base ) )

       ! Now store list of all possible ops that work on local data
       n_all_ops_local = 0
       Do j = 1, a%n
          n = ntran( j )
          Do k = 1, n
             l = ilifc( j ) + k
             If( a%row%global_to_local( itran( l ) ) /= MATRIX_NOT_OWN ) Then
                n_all_ops_local = n_all_ops_local + 1
                all_ops%j_vec_add ( n_all_ops_local ) = j
                all_ops%j_tran    ( n_all_ops_local ) = l
                all_ops%j_vec_base( n_all_ops_local ) = a%row%global_to_local( itran( l ) )
             End If
          End Do
       End Do

       ! Now loop over procs we will recv from and work out which
       ! of the above ops correspond to which proc

       Allocate( remote_ops( 0:a%distribution%n_proc_row - 1 ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( remote_ops ) )

       Do proc_shift = 0, a%distribution%n_proc_row - 1

          proc_from = Modulo( a%distribution%my_proc_row + proc_shift, &
               a%distribution%n_proc_row )

          ! Look for ops that include the appropriate indexing for the
          ! remote data. First count how many.
          n_this_ops_local = 0
          Do i = 1, n_all_ops_local
             global_j = all_ops%j_vec_add( i )
             If( indxg2p( global_j, a%row%blocking_factor, 0, 0, &
                  a%distribution%n_proc_row ) == proc_from ) Then
                n_this_ops_local = n_this_ops_local + 1
             End If
          End Do

          Allocate( remote_ops( proc_shift )%j_vec_add ( 1:n_this_ops_local ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_vec_add ) )
          Allocate( remote_ops( proc_shift )%j_tran    ( 1:n_this_ops_local ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_tran ) )
          Allocate( remote_ops( proc_shift )%j_vec_base( 1:n_this_ops_local ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Matrix allocation error in MATRIX_COMMS_INIT', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_vec_base ) )

          ! Now generate the ops list
          n_this_ops_local = 0
          Do i = 1, n_all_ops_local
             global_j = all_ops%j_vec_add( i )
             If( indxg2p( global_j, a%row%blocking_factor, 0, 0, &
                  a%distribution%n_proc_row ) == proc_from ) Then
                n_this_ops_local = n_this_ops_local + 1
                local_j = indxg2l( global_j, a%row%blocking_factor, 0, 0, &
                     a%distribution%n_proc_row )
                remote_ops( proc_shift )%j_vec_add ( n_this_ops_local ) = local_j
                remote_ops( proc_shift )%j_tran    ( n_this_ops_local ) = all_ops%j_tran    ( i )
                remote_ops( proc_shift )%j_vec_base( n_this_ops_local ) = all_ops%j_vec_base( i )
             End If
          End Do

       End Do

    End If

    Call alloc_free( ALLOC_INTEGER, Size( all_ops%j_vec_base ) ) 
    Deallocate( all_ops%j_vec_base )
    Call alloc_free( ALLOC_INTEGER, Size( all_ops%j_tran ) ) 
    Deallocate( all_ops%j_tran )
    Call alloc_free( ALLOC_INTEGER, Size( all_ops%j_vec_add ) ) 
    Deallocate( all_ops%j_vec_add )

  End Subroutine matrix_tdown_setup

  Subroutine matrix_tdown_free

    Integer :: proc_shift

    Do proc_shift = Ubound( remote_ops, Dim = 1 ), 0, -1
       Call alloc_free( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_vec_base ) )
       Deallocate( remote_ops( proc_shift )%j_vec_base )
       Call alloc_free( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_tran ) )
       Deallocate( remote_ops( proc_shift )%j_tran     )
       Call alloc_free( ALLOC_INTEGER, Size( remote_ops( proc_shift )%j_vec_add ) )
       Deallocate( remote_ops( proc_shift )%j_vec_add  )
    End Do

    Call alloc_free( ALLOC_DERIVED, Size( remote_ops ) )
    Deallocate( remote_ops )

  End Subroutine matrix_tdown_free

  Subroutine matrix_tdown( a, b, ctran, n_harmonic, otran )

    Type( matrix ), Dimension( : )               :: a
    Type( matrix ), Dimension( : )               :: b
    Real( wp )    , Dimension( : ), Intent( In ) :: ctran
    Integer                       , Intent( In ) :: n_harmonic
    Logical                       , Intent( In ) :: otran

    Type( matrix ) :: c

    Integer :: n_harmonic_local
    Integer :: ja, jb, jt
    Integer :: spin
    Integer :: proc_shift
    Integer :: i, op

    Do spin = 1, Size( a )
       If( a( spin )%i_own ) Then

          n_harmonic_local = 0
          Do i = n_harmonic, 1, -1
             If( a( spin )%col%global_to_local( i ) /= MATRIX_NOT_OWN ) Then
                n_harmonic_local = a( spin )%col%global_to_local( i )
                Exit
             End If
          End Do

          If( otran ) Then

             ! No symm adapt
             b( spin )%real_data( :, 1:n_harmonic_local ) = &
                  a( spin )%real_data( :, 1:n_harmonic_local )
             b( spin )%real_data( :, n_harmonic_local + 1: ) = 0.0_wp

          Else

             b( spin )%real_data = 0.0_wp

             Call matrix_create( a( spin ), c )

             Do proc_shift = 0, a( spin )%distribution%n_proc_row - 1

                ! Recv into c data from PROC_SHIFT procs further down the
                ! proc column ( with a mod or two )
                Call matrix_proc_cshift( a( spin ), c, proc_shift, MATRIX_SHIFT_COL )

                ! And apply the transformation
                Do i = 1, n_harmonic_local
                   Do op = 1, Size( remote_ops( proc_shift )%j_vec_add )
                      ja = remote_ops( proc_shift )%j_vec_add ( op )
                      jt = remote_ops( proc_shift )%j_tran    ( op )
                      jb = remote_ops( proc_shift )%j_vec_base( op )
                      b( spin )%real_data( jb, i ) = b( spin )%real_data( jb, i ) + &
                           ctran( jt ) * c%real_data( ja, i )
                   End Do
                End Do

             End Do

             Call matrix_destroy( c )

          End If

       End If
    End Do


  End Subroutine matrix_tdown

  Subroutine matrix_proc_cshift( a, b, shift, row_col )

    Include 'mpif.h'

    Type( matrix )               :: a
    Type( matrix )               :: b
    Integer       , Intent( In ) :: shift
    Integer       , Intent( In ) :: row_col

    Integer, Dimension( 1:mpi_status_size ) :: status

    Integer :: proc_to, proc_from
    Integer :: comm
    Integer :: tag
    Integer :: me
    Integer :: error

    Select Case( row_col )

    Case( MATRIX_SHIFT_COL )

       proc_to   = Modulo( a%distribution%my_proc_row - shift, &
            a%distribution%n_proc_row )
       proc_from = Modulo( a%distribution%my_proc_row + shift, &
            a%distribution%n_proc_row )
       comm = a%distribution%col_comm
       me   = a%distribution%my_proc_row

    Case( MATRIX_SHIFT_ROW )

       proc_to   = Modulo( a%distribution%my_proc_col - shift, &
            a%distribution%n_proc_col )
       proc_from = Modulo( a%distribution%my_proc_col + shift, &
            a%distribution%n_proc_col )
       comm = a%distribution%row_comm
       me   =  a%distribution%my_proc_col

    Case Default
       Call dlc_error( 'Illegal shift operation', 'abort' )

    End Select

    If( proc_to /= me ) Then
       Call mpi_irecv( b%real_data, Size( b%real_data ), mpi_double_precision, &
            proc_from, 10, comm, tag, error )
       Call mpi_send ( a%real_data, Size( a%real_data ), mpi_double_precision, &
            proc_to  , 10, comm, error )
       Call mpi_wait(  tag, status, error )
    Else
       b%real_data = a%real_data
    End If

  End Subroutine matrix_proc_cshift

  Subroutine matrix_stretch_one( a, b )

    Type( matrix ) :: a
    Type( matrix ) :: b

    b%real_data = 0.0_wp
    b%real_data( 1:a%local_n, 1:a%local_m ) = a%real_data( 1:a%local_n, 1:a%local_m )

  End Subroutine matrix_stretch_one

  Subroutine matrix_stretch_multi( a, b )

    Type( matrix ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: b

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_stretch( a( i ), b( i ) )
       End If
    End Do

  End Subroutine matrix_stretch_multi

  Subroutine matrix_direct_product_one( a, b, c )

    Type( matrix ) :: a
    Type( matrix ) :: b
    Type( matrix ) :: c

    c%real_data = a%real_data * b%real_data

  End Subroutine matrix_direct_product_one

  Subroutine matrix_direct_product_multi( a, b, c )

    Type( matrix ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: b
    Type( matrix ), Dimension( : ) :: c

    Integer :: i

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_direct_product( a( i ), b( i ), c( i ) )
       End If
    End Do

  End Subroutine matrix_direct_product_multi

  Subroutine matrix_get_col_biggest_one( a, biggest, ind_biggest )

    Include 'mpif.h'

    Type( matrix )                               :: a
    Real( wp ), Dimension( :, : ), Intent( Out ) :: biggest
    Integer   , Dimension( :, : ), Intent( Out ) :: ind_biggest

    Real( wp ), Dimension( :    ), Allocatable :: work_real_1d
    Real( wp ), Dimension( :, : ), Allocatable :: work_real_2d

    Real( wp ), Dimension( : ), Allocatable :: biggest_local
    Real( wp ), Dimension( : ), Allocatable :: col_local

    Integer, Dimension( :    ), Allocatable :: work_int_1d
    Integer, Dimension( :, : ), Allocatable :: work_int_2d

    Integer, Dimension( : ), Allocatable :: ind_biggest_local

    Real( wp ) :: big
    Real( wp ) :: norm, tmp

    Integer :: n_tot_pops, n_tot_pops_loc
    Integer :: loc_tot, loc_start
    Integer :: glob_m
    Integer :: loc, loc_j
    Integer :: error
    Integer :: i, j, k

    biggest     = 0.0_wp
    ind_biggest = 0

    n_tot_pops     = Size( biggest, Dim = 1 )
    n_tot_pops_loc = Min( Size( biggest, Dim = 1 ), a%local_n )

    loc_tot = n_tot_pops_loc * a%distribution%n_proc_row

    Allocate( col_local( 1:a%local_n ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( col_local ) )
    Allocate( biggest_local( 1:loc_tot ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( biggest_local ) )
    Allocate( ind_biggest_local( 1:loc_tot ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( ind_biggest_local ) )
    Allocate( work_real_1d( 1:loc_tot ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work_real_1d ) )
    Allocate( work_int_1d( 1:loc_tot ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( work_int_1d ) )

    Do i = 1, a%local_m

       glob_m = a%col%local_to_global( i )

       col_local         = a%real_data( 1:a%local_n, i )
       biggest_local     = 0.0_wp
       ind_biggest_local = 0

       loc_start = n_tot_pops_loc * a%distribution%my_proc_row

       Do k = 1, n_tot_pops_loc

          big   = 0.0_wp
          loc   = 1
          loc_j = 1

          Do j = 1, a%local_n
             If( Abs( col_local( j ) ) > Abs( big ) ) Then
                big = col_local( j )
                loc   = j
                loc_j = j
             End If
          End Do

          loc_start = loc_start + 1
          biggest_local    ( loc_start ) = big
          ind_biggest_local( loc_start ) = a%row%local_to_global( loc )
          col_local        ( loc_j     ) = 0.0_wp

       End Do

       Call mpi_allreduce( biggest_local, work_real_1d, Size( biggest_local ), &
            mpi_double_precision, mpi_sum, &
            a%distribution%col_comm, error )
       biggest_local = work_real_1d

       Call mpi_allreduce( ind_biggest_local, work_int_1d, Size( ind_biggest_local ), &
            mpi_integer, mpi_sum, &
            a%distribution%col_comm, error )
       ind_biggest_local = work_int_1d

       Do k = 1, n_tot_pops

          big = 0.0_wp
          loc = 1
          loc_j = 1
       
          Do j = 1, loc_tot
             If( Abs( biggest_local( j ) ) > Abs( big ) ) Then
                big   = biggest_local( j )
                loc   = ind_biggest_local( j )
                loc_j = j
             End If
          End Do

          biggest      ( k, glob_m ) = big
          ind_biggest  ( k, glob_m ) = loc
          biggest_local( loc_j     ) = 0.0_wp

       End Do

    End Do

    Call alloc_free( ALLOC_INTEGER, Size( work_int_1d ) )
    Deallocate( work_int_1d )
    Call alloc_free( ALLOC_REAL, Size( work_real_1d ) )
    Deallocate( work_real_1d )
    Call alloc_free( ALLOC_INTEGER, Size( ind_biggest_local ) )
    Deallocate( ind_biggest_local )
    Call alloc_free( ALLOC_REAL, Size( biggest_local ) )
    Deallocate( biggest_local )
    Call alloc_free( ALLOC_REAL, Size( col_local ) )
    Deallocate( col_local )

    Allocate( work_real_2d( 1:n_tot_pops, 1:a%m ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( work_real_2d ) )
    Call mpi_allreduce( biggest, work_real_2d, Size( biggest ), mpi_double_precision, mpi_sum, &
         a%distribution%row_comm, error )
    biggest = work_real_2d
    Call alloc_free( ALLOC_REAL, Size( work_real_2d ) )
    Deallocate( work_real_2d )

    Allocate( work_int_2d( 1:n_tot_pops, 1:a%m ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( work_int_2d ) )
    Call mpi_allreduce( ind_biggest, work_int_2d, Size( ind_biggest ), mpi_integer, mpi_sum, &
         a%distribution%row_comm, error )
    ind_biggest = work_int_2d
    Call alloc_free( ALLOC_INTEGER, Size( work_int_2d ) )
    Deallocate( work_int_2d )

  End Subroutine matrix_get_col_biggest_one

  Subroutine matrix_get_col_biggest_multi( a, biggest, ind_biggest )

    Use comms_data

    Include 'mpif.h'

    Type( matrix ), Dimension( : )                      :: a
    Real( wp )    , Dimension( :, :, : ), Intent( Out ) :: biggest
    Integer       , Dimension( :, :, : ), Intent( Out ) :: ind_biggest

    Real( wp ), Dimension( :, :, : ), Allocatable :: work_real_3d

    Integer, Dimension( :, :, : ), Allocatable :: work_int_3d

    Integer :: n1, n2, n3
    Integer :: error
    Integer :: i

    n1 = Size( biggest, Dim = 1 )
    n2 = Size( biggest, Dim = 2 )
    n3 = Size( biggest, Dim = 3 )

    biggest     = 0.0_wp
    ind_biggest = 0

    Do i = 1, Size( a )
       If( a( i )%i_own ) Then
          Call matrix_get_col_biggest( a( i ), biggest( :, :, i ), ind_biggest( :, :, i ) )
       End If
    End Do

    If( a( 1 )%distribution%type == MATRIX_DISTRIB_UHF_K ) Then

       Allocate( work_real_3d( 1:n1, 1:n2, 1:n3 ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( work_real_3d ) )
       Call mpi_allreduce( biggest, work_real_3d, Size( biggest ), mpi_double_precision, mpi_sum, &
            comms_data_uhf_k_repl_comm , error )
       biggest = work_real_3d
       Call alloc_free( ALLOC_REAL, Size( work_real_3d ) )
       Deallocate( work_real_3d )
       
       Allocate( work_int_3d( 1:n1, 1:n2, 1:n3 ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Matrix allocation error in MATRIX_GET_COL_BIGGEST_ONE', 'abort' )
       End If
       Call alloc_add( ALLOC_INTEGER, Size( work_int_3d ) )
       Call mpi_allreduce( ind_biggest, work_int_3d, Size( ind_biggest ), mpi_integer, mpi_sum, &
            comms_data_uhf_k_repl_comm, error )
       ind_biggest = work_int_3d
       Call alloc_free( ALLOC_INTEGER, Size( work_int_3d ) )
       Deallocate( work_int_3d )

    End If
       
  End Subroutine matrix_get_col_biggest_multi

End Module distributed_matrices


   
