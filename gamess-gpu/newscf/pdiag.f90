Subroutine pdiag_f90( n, a, q, e, error )

  Use newscf_numbers
  Use newscf_modules, Only : mpi_comm_workers, blacs_context
  Use distributed_matrices

  Implicit None

  Integer                                         , Intent( In    ) :: n
  Real( wp ), Dimension( 1:( n * ( n + 1 ) ) / 2 ), Intent( In    ) :: a
  Real( wp ), Dimension( 1:n * n  )               , Intent(   Out ) :: q
  Real( wp ), Dimension( 1:n      )               , Intent(   Out ) :: e
  Integer                                         , Intent(   Out ) :: error

  Type( matrix ) :: a_mat
  Type( matrix ) :: q_mat

  error = 0

  Call matrix_comms_init( n, mpi_comm_workers, blacs_context, 1, 1 )

  Call matrix_create( n, n, a_mat, 'diag a', MATRIX_REAL, MATRIX_DISTRIB_ALL )
  Call matrix_create( a_mat, q_mat )

  Call matrix_set_from_triangle( a_mat, a )

  Call matrix_diagonalise( a_mat, q_mat, e )

  Call matrix_get_to_square( q_mat, q )

  Call matrix_destroy( q_mat )
  Call matrix_destroy( a_mat )

  Call matrix_comms_finalize

End Subroutine pdiag_f90

Subroutine mult2_f90( q, a, h, n1, n2, n )

  Use newscf_numbers
  Use newscf_modules, Only : mpi_comm_workers, blacs_context
  Use distributed_matrices
  Use allocation

  Implicit None

  Integer                                           , Intent( In    ) :: n
  Integer                                           , Intent( In    ) :: n1
  Integer                                           , Intent( In    ) :: n2
  Real( wp ), Dimension( 1:n * n1  )                , Intent( In    ) :: q
  Real( wp ), Dimension( 1:( n1 * ( n1 + 1 ) ) / 2 ), Intent(   Out ) :: a
  Real( wp ), Dimension( 1:( n  * ( n  + 1 ) ) / 2 ), Intent( In    ) :: h

  Type( matrix ) :: a_mat
  Type( matrix ) :: q_mat
  Type( matrix ) :: h_mat

  Integer :: me, error

  Call matrix_comms_init( n, mpi_comm_workers, blacs_context, 1, 1 )

  Call matrix_create( n1, n1, a_mat, 'diag a', MATRIX_REAL, MATRIX_DISTRIB_ALL )
  Call matrix_create( n , n1, q_mat, 'diag q', MATRIX_REAL, MATRIX_DISTRIB_ALL )
!!$  Call matrix_create( a_mat, q_mat )
  Call matrix_create( n , n , h_mat, 'diag h', MATRIX_REAL, MATRIX_DISTRIB_ALL )
!!$  Call matrix_create( a_mat, h_mat )

  Call matrix_set_from_triangle( h_mat, h )
  Call matrix_set_from_square  ( q_mat, q )

  Call matrix_mult2( h_mat, q_mat, a_mat )

  Call matrix_get_to_triangle( a_mat, a )

  Call matrix_destroy( h_mat )
  Call matrix_destroy( q_mat )
  Call matrix_destroy( a_mat )

  Call matrix_comms_finalize

End Subroutine mult2_f90

Subroutine tfsqc_f90( n_h, n_c, a, b )

  Use newscf_numbers
  Use newscf_modules, Only : mpi_comm_workers, blacs_context
  Use distributed_matrices

  Implicit None

  Integer                             , Intent( In    ) :: n_h
  Integer                             , Intent( In    ) :: n_c
  Real( wp ), Dimension( 1:n_c * n_h ), Intent( In    ) :: a
  Real( wp ), Dimension( 1:n_c * n_h ), Intent( InOut ) :: b

  Type( matrix ) :: a_mat
  Type( matrix ) :: b_mat
  Type( matrix ) :: c_mat

  Call matrix_comms_init( n_c, mpi_comm_workers, blacs_context, 1, 1 )

  Call matrix_create( n_c, n_h, a_mat, 'tfsqc a', MATRIX_REAL, MATRIX_DISTRIB_ALL )
  Call matrix_create( n_h, n_h, b_mat, 'tfsqc b', MATRIX_REAL, MATRIX_DISTRIB_ALL )
  Call matrix_create( n_c, n_h, c_mat, 'tfsqc c', MATRIX_REAL, MATRIX_DISTRIB_ALL )
!!$  Call matrix_create( a_mat, b_mat )
!!$  Call matrix_create( a_mat, c_mat )

  Call matrix_set_from_square( a_mat, a )
  Call matrix_set_from_square( b_mat, b )

  Call matrix_multiply( 1.0_wp, a_mat, b_mat, 0.0_wp, c_mat )

  b = 0.0_wp
  Call matrix_get_to_square( c_mat, b )

  Call matrix_destroy( c_mat )
  Call matrix_destroy( b_mat )
  Call matrix_destroy( a_mat )

  Call matrix_comms_finalize

End Subroutine tfsqc_f90

Subroutine dgemm_f90( ta, tb, n, alpha, a, b, beta, c )

  Use newscf_numbers
  Use newscf_modules, Only : mpi_comm_workers, blacs_context
  Use distributed_matrices

  Implicit None

  Character( Len = 1 )            , Intent( In    ) :: ta
  Character( Len = 1 )            , Intent( In    ) :: tb
  Integer                         , Intent( In    ) :: n
  Real( wp )                      , Intent( In    ) :: alpha
  Real( wp ), Dimension( 1:n * n ), Intent( In    ) :: a
  Real( wp ), Dimension( 1:n * n ), Intent( In    ) :: b
  Real( wp )                      , Intent( In    ) :: beta
  Real( wp ), Dimension( 1:n * n ), Intent(   Out ) :: c

  Type( matrix ) :: a_mat
  Type( matrix ) :: b_mat
  Type( matrix ) :: c_mat

  Call matrix_comms_init( n, mpi_comm_workers, blacs_context, 1, 1 )

  Call matrix_create( n, n, a_mat, 'diag a', MATRIX_REAL, MATRIX_DISTRIB_ALL )
  Call matrix_create( a_mat, b_mat )
  Call matrix_create( a_mat, c_mat )

  Call matrix_set_from_square( a_mat, a )
  Call matrix_set_from_square( b_mat, b )

  Call matrix_dgemm( ta, tb, alpha, a_mat, b_mat, 0.0_wp, c_mat )

  Call matrix_get_to_square( c_mat, c )

  Call matrix_destroy( c_mat )
  Call matrix_destroy( b_mat )
  Call matrix_destroy( a_mat )

  Call matrix_comms_finalize

End Subroutine dgemm_f90
