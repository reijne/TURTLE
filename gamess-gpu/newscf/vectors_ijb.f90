Module distributed_vectors

  Use newscf_numbers
  Use distributed_matrices
  Use allocation

  Implicit None

  Public :: vector
  Public :: vector_create
  Public :: vector_destroy
  Public :: vector_copy
  Public :: vector_make_density
  Public :: vector_assign_by_energy
  Public :: vector_assign_by_overlap
  Public :: vector_assign_by_fermi
  Public :: vector_get_gap
  Public :: vector_read
  Public :: vector_write
  Public :: vector_similarity
  Public :: vector_diagonalise
  Public :: vector_back_transform
  Public :: vector_level_shift
  Public :: vector_check_ov
  Public :: vector_save_orbitals
  Public :: vector_diis_initialize
  Public :: vector_diis_reset
  Public :: vector_diis_solve
  Public :: vector_diis_free_scratch
  Public :: vector_orthogonalize
  Public :: vector_tdown_setup
  Public :: vector_tdown
  Public :: vector_initialize
  Public :: vector_get_biggest_orb_pops

  Type vector
!Early INTEL compiler, i.e. my laptop, falls over on this - remember to uncomment !!!!!!!!!!!!!
     Private
     Integer                             :: spin
     Integer                             :: k
     Integer                             :: n_basis
     Integer                             :: n_vecs
     Type( matrix )                      :: evecs
     Real( wp ), Dimension( : ), Pointer :: evals
     Real( wp ), Dimension( : ), Pointer :: occupations
  End Type vector

  Private

  Interface vector_create
     Module Procedure vector_create_one
     Module Procedure vector_create_multi
  End Interface

  Interface vector_destroy
     Module Procedure vector_destroy_one
     Module Procedure vector_destroy_multi
  End Interface

  Interface vector_copy
     Module Procedure vector_copy_one
     Module Procedure vector_copy_multi
  End Interface

  Interface vector_diagonalise
     Module Procedure vector_diagonalise_one
     Module Procedure vector_diagonalise_multi
  End Interface

  Interface vector_initialize
     Module Procedure vector_initialize_one
     Module Procedure vector_initialize_multi
  End Interface

  ! DIIS data

  Type( matrix ), Dimension( :, :, : ), Allocatable, Save :: fock_matrices

  Integer, Parameter :: maxdiis = 8
  Integer, Parameter :: DIIS_AO = 1
  Integer, Parameter :: DIIS_OV = 2

  Real( wp ), Dimension( 1:maxdiis + 1, 1:maxdiis + 1 ), Save :: diis_matrix
  Real( wp ), Dimension( 1:maxdiis + 1, 1:maxdiis + 1 ), Save :: diis_matrix_copy
  Real( wp ), Dimension( 1:maxdiis + 1                )       :: diis_soln
  Real( wp ), Dimension( 1:maxdiis + 1                )       :: diis_scale

  Integer, Dimension( 1:maxdiis + 1 ) :: diis_pivot

  Integer :: nstore
  Logical :: ondiis
  Logical :: declared = .False.

Contains

  Subroutine vector_create_one( a, k, spin, nbas, nvecs, name, genus, distribution )

    Integer             , Intent( In )           :: k
    Integer             , Intent( In )           :: spin
    Type( vector )                               :: a
    Integer             , Intent( In )           :: nbas
    Integer             , Intent( In )           :: nvecs
    Character( Len = * ), Intent( In )           :: name
    Integer             , Intent( In )           :: genus
    Integer             , Intent( In ), Optional :: distribution

    Integer :: this_distribution
    Integer :: error

    If( Present( distribution ) ) Then
       this_distribution = distribution
    Else
       this_distribution = MATRIX_DISTRIB_DEFAULT
    End If

    Call matrix_create( nbas, nvecs, a%evecs, name, genus, this_distribution )

    a%k       = k
    a%spin    = spin
    a%n_basis = nbas
    a%n_vecs  = nvecs

    Allocate( a%evals( 1:nvecs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Failed to allocate EVALS in vector_create', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( a%evals ) )

    Allocate( a%occupations( 1:nvecs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Failed to allocate OCCUPATIONS in vector_create', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( a%occupations ) )

  End Subroutine vector_create_one

  Subroutine vector_create_multi( a, k, spin, nbas, nvecs, name, genus, distribution )

    Type( vector )      , Dimension( : )                         :: a
    Integer             , Dimension( : ), Intent( In )           :: k
    Integer             , Dimension( : ), Intent( In )           :: spin
    Integer             ,                 Intent( In )           :: nbas
    Integer             ,                 Intent( In )           :: nvecs
    Character( Len = * ), Dimension( : ), Intent( In )           :: name
    Integer             , Dimension( : ), Intent( In )           :: genus
    Integer                             , Intent( In ), Optional :: distribution

    Integer :: this_distribution
    Integer :: error
    Integer :: i

    If( Present( distribution ) ) Then
       this_distribution = distribution
    Else
       this_distribution = MATRIX_DISTRIB_DEFAULT
    End If

    Call matrix_create( nbas, nvecs, a%evecs, name, genus, this_distribution ) 

    Do i = 1, Size( a )

       a( i )%k       = k( i )
       a( i )%spin    = spin( i )
       a( i )%n_basis = nbas
       a( i )%n_vecs  = nvecs
       
       Allocate( a( i )%evals( 1:nvecs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Failed to allocate EVALS in vector_create', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( a( i )%evals ) )
       
       Allocate( a( i )%occupations( 1:nvecs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Failed to allocate OCCUPATIONS in vector_create', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( a( i )%occupations ) )

    End do

  End Subroutine vector_create_multi

  Subroutine vector_destroy_one( a )

    Type( vector ) :: a

    Call alloc_free( ALLOC_REAL, Size( a%occupations ) )
    Deallocate( a%occupations )
    Nullify   ( a%occupations )
    Call alloc_free( ALLOC_REAL, Size( a%evals ) )
    Deallocate( a%evals )
    Nullify   ( a%evals )

    Call matrix_destroy( a%evecs )

  End Subroutine vector_destroy_one

  Subroutine vector_destroy_multi( a )

    Type( vector ), Dimension( : ) :: a

    Integer :: i

    Do i = Size( a ), 1, -1
       Call vector_destroy( a( i ) )
    End Do

  End Subroutine vector_destroy_multi

  Subroutine vector_copy_one( a, b )

    Type( vector ) :: a
    Type( vector ) :: b

    Call matrix_copy( a%evecs, b%evecs )
    
    b%evals = a%evals
    b%occupations = a%occupations

    b%spin    = a%spin
    b%k       = a%k
    b%n_basis = a%n_basis
       
  End Subroutine vector_copy_one

  Subroutine vector_copy_multi( a, b )

    Type( vector ), Dimension( : ) :: a
    Type( vector ), Dimension( : ) :: b

    Integer :: a_distrib, b_distrib
    Integer :: i

    Logical :: a_own, b_own

    Call matrix_copy( a%evecs, b%evecs )

    Do i = 1, Size( a )

       b( i )%evals = a( i )%evals
       b( i )%occupations = a( i )%occupations
       
       b( i )%spin    = a( i )%spin
       b( i )%k       = a( i )%k
       b( i )%n_basis = a( i )%n_basis
       
    End Do
    
  End Subroutine vector_copy_multi

  ! PRINITNG ROUTINES to be done

  Subroutine vector_make_density( density, a )

    Type( matrix ), Dimension( : ) :: density
    Type( vector ), Dimension( : ) :: a

    Type( matrix ), Dimension( : ), Allocatable :: work

    Real( wp ), Dimension( :, : ), Allocatable :: occs

    Integer :: error
    Integer :: i

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_MAKE_DENSITY', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( work ) )
    Allocate( occs( 1:Size( a( 1 )%occupations ), 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_MAKE_DENSITY', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( occs ) )

    Do i = 1, Size( a )
       occs( :, i ) = a( i )%occupations
    End Do

    Call matrix_create( a%evecs, work )

    Call matrix_mult_by_diag( a%evecs, occs, work )
    Call matrix_dgemm( 'n', 't', 1.0_wp, work, a%evecs, 0.0_wp, density )

    Call matrix_destroy( work )
    
    Call alloc_free( ALLOC_REAL, Size( occs ) )
    Deallocate( occs )
    Call alloc_free( ALLOC_DERIVED, Size( work ) )
    Deallocate( work )

  End Subroutine vector_make_density

  Subroutine vector_assign_by_energy( a, n_occ )

    ! HAVEN'T CONSIDERED DEGENERACY IN RHF CASE
    ! Also weights not v. well done.

    Type( vector ), Dimension( : )               :: a
    Integer       , Dimension( : ), Intent( In ) :: n_occ

    Real( wp ) :: weight

    Integer :: i

    Logical :: uhf

    uhf = Size( a ) /= 1
    
    If( uhf ) Then
       weight = 1.0_wp
    Else
       weight = 2.0_wp
    End If
    
    do i = 1, Size( a )
       a( i )%occupations( :n_occ( i )     ) = weight
       a( i )%occupations( n_occ( i ) + 1: ) = 0.0_wp
    End do
    
  End Subroutine vector_assign_by_energy

  Subroutine vector_assign_by_fermi( a, t, n_occ, e_mermin )

    Type( vector ), Dimension( : )                  :: a
    Real( wp )                    , Intent( In    ) :: t
    Integer       , Dimension( : ), Intent( In    ) :: n_occ
    Real( wp )                    , Intent(   Out ) :: e_mermin

    Real( wp ) :: etmp

    Integer :: nfocc, noccmx
    Integer :: i

    e_mermin = 0.0_wp

    Do i = 1, Size( a )
       Call fermi_smear( a( i )%evals, a( i )%occupations, etmp, nfocc, noccmx, &
            t, n_occ( i ), Size( a( i )%evals ), 6, .False. )
       e_mermin = e_mermin  + etmp
    End Do

    If( Size( a ) == 1 ) Then
       a( 1 )%occupations = a( 1 )%occupations * 2.0_wp
       e_mermin = e_mermin * 2.0_wp
    End If

  End Subroutine vector_assign_by_fermi

  Subroutine vector_get_gap( a, n_occ, egap )

    Type( vector ), Dimension( : )                  :: a
    Integer       , Dimension( : ), Intent( In    ) :: n_occ
    Real( wp )    ,                 Intent(   Out ) :: egap

    Integer :: i

    egap = -1.0_wp
    Do i = 1, Size( a )
       egap = Max( egap, a( i )%evals( n_occ( i + 1 ) ) - a( i )%evals( n_occ( i ) ) )
    End Do

  End Subroutine vector_get_gap

  Subroutine vector_assign_by_overlap( a, old_a, n_occ, check )

    Type( vector ), Dimension( : )                  :: a
    Type( vector ), Dimension( : )                  :: old_a
    Integer       , Dimension( : ), Intent( In    ) :: n_occ
    Logical                       , Intent(   Out ) :: check

    Type( matrix ), Dimension( : ), Allocatable :: work
    Type( matrix ), Dimension( : ), Allocatable :: overlap

    Logical, External :: opg_root

    Integer, Dimension( :, : ), Allocatable :: similar
    
    Integer :: genus
    Integer :: electrons
    Integer :: error
    Integer :: i, j

    Allocate( work( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_ASSIGN_BY_OVERLAP', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( work ) )
    Allocate( overlap( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_ASSIGN_BY_OVERLAP', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( overlap ) )
    Allocate( similar( 1:a( 1 )%n_basis, 1:Size( a ) ) )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_ASSIGN_BY_OVERLAP', 'abort' )
    End If
    Call alloc_add( ALLOC_INTEGER, Size( similar ) )

    Call matrix_create( a%evecs, work )
    Call matrix_create( a%evecs, overlap )

    Call matrix_copy( old_a%evecs, work )
    Call matrix_invert( work )

    Call matrix_dgemm( 'n', 'n', 1.0_wp, work, a%evecs, 0.0_wp, overlap )
    Call matrix_absmaxloc_col( overlap, similar )
    check = .True.
    Do i = 1, Size( a )
       electrons = 0
       Do j = 1, a( 1 )%n_basis
          a( i )%occupations( j ) = old_a( i )%occupations( similar( j, i ) )
          ! What about fractional occupations ?
          If( a( i )%occupations( j ) > 1.0E-10_wp ) Then
             electrons = electrons + 1
          End If
       End Do
       check = check .And. electrons == n_occ( i )
    End Do

    Call matrix_destroy( overlap )
    Call matrix_destroy( work )

    Call alloc_free( ALLOC_INTEGER, Size( similar ) )
    Deallocate( similar )
    Call alloc_free( ALLOC_DERIVED, Size( overlap ) )
    Deallocate( overlap )
    Call alloc_free( ALLOC_DERIVED, Size( work ) )
    Deallocate( work )

  End Subroutine vector_assign_by_overlap

  Subroutine vector_read( a, where_vecs, where_vals, idaf )

    Type( vector ), Dimension( : )                  :: a
    Integer       , Dimension( : ), Intent( In    ) :: where_vecs
    Integer       , Dimension( : ), Intent( In    ) :: where_vals
    Integer                       , Intent( In    ) :: idaf

    Integer :: n, m
    Integer :: i

    Do i = 1, Size( a )
       Call matrix_inquire( a( i )%evecs, global_m = m )
       Call rdedx_prec( a( i )%evals, m, where_vals( i ),  idaf )
       Call matrix_read_rectangle( a( i )%evecs, where_vecs( i ), idaf )
    End Do

  End Subroutine vector_read

  Subroutine vector_write( a, where_vecs, where_vals, idaf )

    Type( vector ), Dimension( : )                  :: a
    Integer       , Dimension( : ), Intent( In    ) :: where_vecs
    Integer       , Dimension( : ), Intent( In    ) :: where_vals
    Integer                       , Intent( In    ) :: idaf

    Integer :: n, m
    Integer :: i

    Do i = 1, Size( a )
       Call matrix_inquire( a( i )%evecs, global_m = m )
       Call wrt3( a( i )%evals, m, where_vals( i ),  idaf )
       Call matrix_write_rectangle( a( i )%evecs, where_vecs( i ), idaf )
    End Do

  End Subroutine vector_write

  Subroutine vector_similarity( fock_ao, a, fock_mo )

    Type( matrix ), Dimension( : ) :: fock_ao
    Type( vector ), Dimension( : ) :: a
    Type( matrix ), Dimension( : ) :: fock_mo

    Call matrix_mult2( fock_ao, a%evecs, fock_mo )

  End Subroutine vector_similarity

  Subroutine vector_diagonalise_one( mat, a )

    Use comms_data

    Type( matrix ) :: mat
    Type( vector ) :: a

    Integer :: num
    Integer :: n_jobs
    Integer :: error
    Integer :: file
    Integer :: i, j

    Call matrix_diagonalise( mat, a%evecs, a%evals )

  End Subroutine vector_diagonalise_one

  Subroutine vector_diagonalise_multi( mat, a )

    Use comms_data

    !Needs getting rid of comms ! - move to matrix
    Include 'mpif.h'

    Type( matrix ), Dimension( : ) :: mat
    Type( vector ), Dimension( : ) :: a

    !Hack - only need one when all distrib
    Real( wp ), Dimension( :, : ), Allocatable :: evals_work
    Real( wp ), Dimension( :, : ), Allocatable :: evals_work_2

    Integer :: num
    Integer :: n_jobs
    Integer :: error
    Integer :: file
    Integer :: i, j

    Logical :: i_own, i_own_e

    n_jobs = Size( mat )

    Call matrix_inquire( mat( 1 ), global_n = num )

    Allocate( evals_work( 1:num, 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_DIAGONALISE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( evals_work ) )
    Allocate( evals_work_2( 1:num, 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: Unable to allocate memory in ' // &
            'VECTOR_DIAGONALISE', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( evals_work ) )
    evals_work = 0.0_wp

    Call matrix_diagonalise( mat, a%evecs, evals_work )

    ! Temp hack - Get rid of when all distrib
    If( .Not. comms_data_uhf_k_is_all ) Then
       Call mpi_allreduce( evals_work, evals_work_2, Size( evals_work ), mpi_double_precision, &
            mpi_sum, comms_data_uhf_k_repl_comm, error )
       Do i = 1, n_jobs
          a( i )%evals = evals_work_2( :, i )
       End Do
    Else
       Do i = 1, n_jobs
          a( i )%evals = evals_work( :, i )
       End Do
    End If

    Call alloc_free( ALLOC_REAL, Size( evals_work_2 ) )
    Deallocate( evals_work_2 )
    Call alloc_free( ALLOC_REAL, Size( evals_work ) )
    Deallocate( evals_work   )

  End Subroutine vector_diagonalise_multi

  Subroutine vector_back_transform( evecs_mo, a, evecs_ao )

    Type( vector ), Dimension( : ) :: evecs_mo
    Type( vector ), Dimension( : ) :: a
    Type( vector ), Dimension( : ) :: evecs_ao

    Type( matrix ), Dimension( : ), Allocatable :: work

    Integer :: n
    Integer :: distribution
    Integer :: n_jobs
    Integer :: error
    Integer :: i

    Call matrix_multiply( 1.0_wp, a%evecs, evecs_mo%evecs, 0.0_wp, evecs_ao%evecs )

    Do i = 1, Size( evecs_mo )
       evecs_ao( i )%evals = evecs_mo( i )%evals
    End Do

  End Subroutine vector_back_transform

  Subroutine vector_level_shift( fock, a, shift_vals )

    Type( matrix ), Dimension( : )               :: fock
    Type( vector ), Dimension( : )               :: a
    Real( wp )    , Dimension( : ), Intent( In ) :: shift_vals

    Real( wp ), Dimension( :, : ), Allocatable :: shift_vec

    Integer :: n
    Integer :: error
    Integer :: i, j

    Call matrix_inquire( fock( 1 ), global_n = n )

    Allocate( shift_vec( 1:n, 1:Size( fock ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in&
            & VECTOR_LEVEL_SHIFT', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( shift_vec ) )

    Do i = 1, Size( fock )

       Do j = 1, n
          If( a( i )%occupations( j ) <= 0.001_wp ) Then
             shift_vec( j, i ) = shift_vals( i )
          Else
             shift_vec( j, i ) = 0.0_wp
          End If
       End Do

    End Do

    Call matrix_add_diagonal( fock, shift_vec )

    Call alloc_free( ALLOC_REAL, Size( shift_vec ) )
    Deallocate( shift_vec )

  End Subroutine vector_level_shift

  Real( wp ) Function vector_check_ov( fock, a )
    
    Type( matrix ), Dimension( : ) :: fock
    Type( vector ), Dimension( : ) :: a

    Real( wp ) :: max_ov

    Real( wp ), Dimension( : ), Allocatable :: ov

    Logical, Dimension( :, : ), Allocatable :: these_occs

    Integer :: error
    Integer :: i, j

    Allocate( ov( 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_CHECK_OV', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( ov ) )
    Allocate( these_occs( 1:Size( a( 1 )%occupations ), 1:Size( a ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_CHECK_OV', 'abort' )
    End If
    Call alloc_add( ALLOC_LOGICAL, Size( these_occs ) )
    
    max_ov = -1.0_wp

    Do i = 1, Size( a )
       Do j = 1, Size( these_occs, Dim = 1 )
          these_occs( j, i ) = a( i )%occupations( j ) > 0.001_wp
       End Do
    End Do
    Call matrix_check_abs_ov( fock, these_occs, ov )
    max_ov = Maxval( ov )

    Call alloc_free( ALLOC_LOGICAL, Size( these_occs ) )
    Deallocate( these_occs )
    Call alloc_free( ALLOC_REAL, Size( ov ) )
    Deallocate( ov )

    vector_check_ov = max_ov

  End Function vector_check_ov

  Subroutine vector_save_orbitals( vectors, density )

    Use newscf_modules

    Type( vector ), Dimension( : ) :: vectors
    Type( matrix ), Dimension( : ) :: density

    Type( matrix ), Dimension( : ), Allocatable :: vectors_tmp

    Integer :: error
    Integer :: l1, l2
    Integer :: l0
    Integer :: n_spin
    Integer :: ndaf
    Integer :: distrib
    Integer :: i

    Logical :: uhf

    Call matrix_inquire( vectors( 1 )%evecs, global_n = l1 )
    Call matrix_inquire( vectors( 1 )%evecs, global_m = l0 )
    Call matrix_inquire( vectors( 1 )%evecs, distrib  = distrib )
    n_spin = Size( vectors )
    uhf = n_spin == 2


    If( l0 == l1 ) Then

       ! Non-Harmonic case

       ndaf = mouta
       Call scfsav_new( vectors( 1 )%evecs, density( 1 ), vectors( 1 )%evals, vectors( 1 )%occupations, &
            ndaf, l1, l2, ibl3pa, ibl3ea )
       
       If( uhf ) Then
          
          ndaf = moutb
          Call scfsav_new( vectors( 2 )%evecs, density( 2 ), vectors( 2 )%evals, vectors( 2 )%occupations, &
               ndaf, l1, l2, ibl3pb, ibl3eb )
          
       End If

    Else

       ! Harmonic case
       Allocate( vectors_tmp( 1:Size( vectors ) ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_SAVE_ORBITALS', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( vectors_tmp ) )
       Call matrix_create( l1, l1, vectors_tmp, (/ ( 'vectors_tmp', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), distrib )

       Call matrix_stretch( vectors%evecs, vectors_tmp )

       ndaf = mouta
       Call scfsav_new( vectors_tmp( 1 ), density( 1 ), vectors( 1 )%evals, vectors( 1 )%occupations, &
            ndaf, l1, l2, ibl3pa, ibl3ea )
       
       If( uhf ) Then
          
          ndaf = moutb
          Call scfsav_new( vectors_tmp( 2 ), density( 2 ), vectors( 2 )%evals, vectors( 2 )%occupations, &
               ndaf, l1, l2, ibl3pb, ibl3eb )
          
       End If

       Call matrix_destroy( vectors_tmp )
       Call alloc_free( ALLOC_DERIVED, Size( vectors_tmp ) )
       Deallocate( vectors_tmp )

    End If

  End Subroutine vector_save_orbitals

  Subroutine vector_diis_initialize( n_cart, n_harm, uhf )

    Integer, Intent( In ) :: n_cart
    Integer, Intent( In ) :: n_harm
    Logical, Intent( In ) :: uhf

    Integer :: error
    Integer :: nspin
    Integer :: i, j

    nstore = 0
    ondiis = .False.

    If( .Not. declared ) Then

       If( .Not. uhf ) Then
          nspin = 1
       Else
          nspin = 2
       End If

       Allocate( fock_matrices( 1:maxdiis, 1:2, 1:nspin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_DIIS_INITIALIZE', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( fock_matrices ) )

       Do i = 1, maxdiis
          Call matrix_create( n_cart, n_cart, fock_matrices( i, DIIS_AO, : ), &
               (/ ( 'DIIS AO', j = 1, nspin ) /),     &
               (/ ( MATRIX_REAL, j = 1, nspin ) /), MATRIX_DISTRIB_UHF_K )
          Call matrix_create( n_cart, n_cart, fock_matrices( i, DIIS_OV, : ), &
               (/ ( 'DIIS MO', j = 1, nspin ) /),     &
               (/ ( MATRIX_REAL, j = 1, nspin ) /), MATRIX_DISTRIB_UHF_K )
       End Do

       declared = .True.

    End If

  End Subroutine vector_diis_initialize

  Subroutine vector_diis_free_scratch

    Integer :: i

    If( declared ) Then

       Do i = maxdiis, 1, -1
          Call matrix_destroy( fock_matrices( i, DIIS_OV, : ) )
          Call matrix_destroy( fock_matrices( i, DIIS_AO, : ) )
       End Do
       
       Call alloc_free( ALLOC_DERIVED, Size( fock_matrices ) )
       Deallocate( fock_matrices )

       declared = .False.

    End If

  End Subroutine vector_diis_free_scratch

  Subroutine vector_diis_reset

    nstore = 0
    ondiis = .False.

  End Subroutine vector_diis_reset

  Subroutine vector_diis_solve( fock, vectors, overlap, diis1, diis_on, diis_dimension, &
       tester, diis_error, accdi1 )

    Use newscf_modules

    Type( matrix ), Dimension( : ) :: fock
    Type( vector ), Dimension( : ) :: vectors
    Type( matrix ), Dimension( : ) :: overlap
    Logical      , Intent( In    ) :: diis1
    Logical      , Intent(   Out ) :: diis_on
    Integer      , Intent(   Out ) :: diis_dimension
    Real( wp )   , Intent(   Out ) :: tester
    Real( wp )   , Intent(   Out ) :: diis_error
    Real( wp )   , Intent( In    ) :: accdi1

    Logical, External :: opg_root
    Integer, External :: ipg_nodeid

    Type( matrix ), Dimension( : ), Allocatable :: scratch
    Type( matrix ), Dimension( : ), Allocatable :: scratch2

    Real( wp ), Dimension( : ), Allocatable :: dot_prod_array

    Real( wp ) :: dot_prod

    Integer :: ispace, nmin, nmax
    Integer :: ipos
    Integer :: n, n_vecs, nspin
    Integer :: ndim
    Integer :: error
    Integer :: i, j

    Logical, Dimension( :, : ), Allocatable :: occupied

    diis_dimension = 0

    ! Should we start storing stuff ?
    nmin = 3
    nmax = maxdiis
    ispace = Sqrt( Real( nx, Kind = wp ) ) + 0.01_wp
    If( nmax >= ispace ) Then
       nmax = ispace - 1
    End If
    If( nmin >= nmax ) then
       nmin = nmax - 1
    End If
    If( nmin < 1 ) Then
       ! Too few degrees of freedom - skip diis
       ondiis = .False.
    Else

       n      = vectors( 1 )%n_basis
       n_vecs = vectors( 1 )%n_vecs
       nspin  = Size( fock )

       ! Store info and do DIIS if required

       ipos = Mod( nstore, nmax ) + 1

       ! Store AO rep of fock matrix
       Call matrix_copy ( fock, fock_matrices( ipos, DIIS_AO, : ) )

       ! Transform to MO 
       Allocate( scratch( 1:nspin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_DIIS_SOLVE', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( scratch ) )
       Call matrix_create( n_vecs, n_vecs, scratch, (/ ( 'diis scratch', i = 1, nspin ) /), &
            (/ ( MATRIX_REAL, i = 1, nspin ) /), MATRIX_DISTRIB_UHF_K )
       Call matrix_mult2( fock, vectors%evecs, scratch )

       ! Zero OO and VV blocks
       Allocate( occupied( 1:n_vecs, 1:nspin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_DIIS_SOLVE', 'abort' )
       End If
       Call alloc_add( ALLOC_LOGICAL, Size( occupied ) )
       Do i =1, nspin
          occupied( :, i ) = vectors( i )%occupations > 0.001_wp
       End Do
       Call matrix_zero_patch( scratch,       occupied,       occupied )
       Call matrix_zero_patch( scratch, .Not. occupied, .Not. occupied )
       Call alloc_free( ALLOC_LOGICAL, Size( occupied ) )
       Deallocate( occupied )

       ! Calc tester
       tester = vector_check_ov( scratch, vectors )

       ! If tester too big, or converged, switch off DIIS
       If( tester > accdi1 .Or. tester == 0.0_wp ) Then
          nstore = 0
          diis_dimension = 0
          ondiis = .False.
          If( opg_root() ) Then
             Write( iwr, * ) 'Diis suppressed (tester too high)', tester
          End If

          Call matrix_destroy( scratch )
          Call alloc_free( ALLOC_DERIVED, Size( scratch ) )
          Deallocate( scratch )

       Else

          ! Transform back to the AO basis
          Allocate( scratch2( 1:nspin ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_DIIS_SOLVE', 'abort' )
          End If
          Call alloc_add( ALLOC_DERIVED, Size( scratch2 ) )
          Call matrix_create( overlap, scratch2 )
          Call matrix_mult2t( scratch, vectors%evecs, scratch2 )
          Call matrix_mult2 ( scratch2, overlap, fock_matrices( ipos, DIIS_OV, : ) )
          Call matrix_destroy( scratch2 )
          Call alloc_free( ALLOC_DERIVED, Size( scratch2 ) )
          Deallocate( scratch2 )
          Call matrix_destroy( scratch )
          Call alloc_free( ALLOC_DERIVED, Size( scratch ) )
          Deallocate( scratch )

          ! Need to think more about this bit - what about
          ! wrap around when stored > nmax ?
          ! O.K. - I get v. confused about the indexing here so this is how it 
          ! should work:
          ! IPOS  - which row/column of the DIIS matrix the dot prods will be 
          !         written to. This wraps around
          ! NDIM  - the actual dimension of the problem, including the 
          !         constraint
          ! NMAX  - the biggest prob we can handle NOT including the constraint
          ! Note also the the constraint row/column is always written directly 
          ! after the last row `in use' - so for the first few calls it is 
          ! IPOS + 1 and then just before we start wrapping around it is in 
          ! row/column NDIM for ever more.

          ndim = ipos + 1
          If( nstore >= nmax ) Then
             ndim = nmax + 1
          End If
          nstore = nstore + 1

          Allocate( dot_prod_array( 1:Size( fock_matrices, Dim = 3 ) ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_DIIS_SOLVE', 'abort' )
          End If
          Call alloc_add( ALLOC_REAL, Size( dot_prod_array ) )

          Do i = 1, ndim - 1
             Call matrix_dot_product( fock_matrices( i, DIIS_OV, : ), &
                  fock_matrices( ipos, DIIS_OV, : ), dot_prod_array )
             dot_prod = Sum( dot_prod_array )
             diis_matrix( i, ipos ) = dot_prod
             diis_matrix( ipos, i ) = dot_prod
             diis_matrix( i, ndim ) = -1.0_wp
             diis_matrix( ndim, i ) = -1.0_wp
             diis_soln( i )         = 0.0_wp
          End Do

          Call alloc_free( ALLOC_REAL, Size( dot_prod_array ) )
          Deallocate( dot_prod_array )

          diis_matrix( ndim, ndim ) = 0.0_wp
          diis_soln( ndim ) = -1.0_wp

          If ( nstore <= nmin .Or. .Not. diis1 ) Then
             ondiis = .False.
          Else
             ondiis = .True.
             Do i = 1, ndim - 1
                diis_scale( i ) = Sqrt( diis_matrix( i, i ) )
             End Do
             diis_scale( ndim ) = 1.0_wp

             Do i = 1, ndim
                Do j = 1, ndim
                   diis_matrix_copy( j, i ) = diis_matrix( j, i ) / &
                        ( diis_scale( i ) * diis_scale( j ) )
                End Do
             End Do

             diis_scale( ndim ) = Abs( diis_matrix_copy( ndim, 1 ) )
             diis_matrix_copy( :, ndim ) = diis_matrix_copy( :, ndim ) / diis_scale( ndim )
             diis_matrix_copy( ndim, : ) = diis_matrix_copy( ndim, : ) / diis_scale( ndim )
             diis_soln( ndim ) = diis_soln( ndim ) / diis_scale( ndim )

             ! And actually solve the eqns - IN SERIAL - too small to bother about parallelism
             Call dgesv( ndim, 1, diis_matrix_copy, maxdiis + 1, &
                  diis_pivot, diis_soln, maxdiis + 1, error )

             If( error /= 0 ) Then

                ! Eqn solve failed - reset DIIS but store current position

                Call matrix_copy( fock_matrices( ipos, DIIS_AO, : ), fock_matrices( 1, DIIS_AO, : ) )
                Call matrix_copy( fock_matrices( ipos, DIIS_OV, : ), fock_matrices( 1, DIIS_OV, : ) )

                nstore = 1
                diis_matrix( 1, 1 ) = diis_matrix( ipos, ipos )

                ondiis = .False.

             Else

                diis_soln( 1:ndim ) = diis_soln( 1:ndim ) / diis_scale( 1:ndim )
                diis_error = diis_soln( ndim )

                Call matrix_set( fock, 0.0_wp )
                Do i = 1, ndim - 1
                   Call matrix_daxpy( diis_soln( i ), fock_matrices( i, DIIS_AO, : ), fock )
                End Do

                diis_dimension = ndim

             End If

          End If

       End If

    End If

    diis_on = ondiis

  End Subroutine vector_diis_solve

  Subroutine vector_orthogonalize( vectors, overlap, have_metric )

    Type( vector ), Dimension( : ) :: vectors
    Type( matrix ), Dimension( : ) :: overlap
    Logical, Intent( In )          :: have_metric

    Integer :: n, n_spin
    Integer :: i
    
    ! Do we require to transform the overlap matrix to the MO basis ???????
    ! Looks like we do not have to, but I'll leave this comment
    ! to remind me in case something funny occurs ....

    Call matrix_orthogonalize( vectors%evecs, overlap, have_metric )

  End Subroutine vector_orthogonalize

  Subroutine vector_tdown_setup( a, ctran, ntran, itran, ilifc, otran )

    Type( vector ), Dimension( : )               :: a
    Real( wp )    , Dimension( : ), Intent( In ) :: ctran
    Integer       , Dimension( : ), Intent( In ) :: ntran
    Integer       , Dimension( : ), Intent( In ) :: itran
    Integer       , Dimension( : ), Intent( In ) :: ilifc
    Logical                       , Intent( In ) :: otran

    Integer :: use

    Logical :: i_own

    Call matrix_inquire( a( 1 )%evecs, i_own = i_own )
    If( i_own ) Then
       use = 1
    Else
       use = 2
    End If

    Call matrix_tdown_setup( a( use )%evecs, ctran, ntran, itran, ilifc, otran )

  End Subroutine vector_tdown_setup

  Subroutine vector_tdown( a, b, ctran, n_harmonic, otran )

    Type( vector ), Dimension( : )               :: a
    Type( vector ), Dimension( : )               :: b
    Real( wp )    , Dimension( : ), Intent( In ) :: ctran
    Integer                       , Intent( In ) :: n_harmonic
    Logical                       , Intent( In ) :: otran

    Call vector_copy( a, b )
    Call matrix_tdown( a%evecs, b%evecs, ctran, n_harmonic, otran )

  End Subroutine vector_tdown

  Subroutine vector_initialize_one( n_cart, n_harm, tol_warn, tol_project, overlap, n_use, evecs )

    Integer       , Intent( In    ) :: n_cart
    Integer       , Intent( In    ) :: n_harm
    Real( wp )    , Intent( In    ) :: tol_warn
    Real( wp )    , Intent( In    ) :: tol_project
    Type( matrix )                  :: overlap
    Integer       , Intent(   Out ) :: n_use
    Type( vector )                  :: evecs

    Logical, External :: opg_root

    Type( vector ) :: work1
    Type( vector ) :: work2

    Integer :: start
    Integer :: error
    Integer :: i

    Call vector_create( work1, 0, 1, n_cart, n_cart, 'Work1', &
         MATRIX_REAL, MATRIX_DISTRIB_ALL )
       
    Call vector_diagonalise( overlap, work1 )
       
    n_use = Count( work1%evals >= tol_project )
    
    Call vector_create( evecs, 0, 1, n_cart, n_use, 'evecs', &
         MATRIX_REAL, MATRIX_DISTRIB_ALL )
    
    Call vector_create( work2, 0, 1, n_cart, n_use, 'Work2', &
         MATRIX_REAL, MATRIX_DISTRIB_ALL )
    
    start = n_cart - n_use + 1
    
    evecs%evals( 1:n_use ) = work1%evals( start:n_cart )
    work2%evals = 1.0_wp / Sqrt( work1%evals( start:n_cart ) )
    
    Call matrix_copy_cols( start, n_cart, work1%evecs, &
         1    , n_use , work2%evecs )
    
    Call matrix_mult_by_diag( work2%evecs, work2%evals, evecs%evecs )
    
    Call vector_destroy( work2 )
    Call vector_destroy( work1 )

  End Subroutine vector_initialize_one

  Subroutine vector_initialize_multi( n_cart, n_harm, tol_warn, tol_project, overlap, n_use, evecs )

    Integer                       , Intent( In    ) :: n_cart
    Integer                       , Intent( In    ) :: n_harm
    Real( wp )                    , Intent( In    ) :: tol_warn
    Real( wp )                    , Intent( In    ) :: tol_project
    Type( matrix ), Dimension( : )                  :: overlap
    Integer                       , Intent(   Out ) :: n_use
    Type( vector ), Dimension( : )                  :: evecs

    Type( vector ), Dimension( : ), Allocatable :: work1
    Type( vector ), dimension( : ), Allocatable :: work2

    Real( wp ), Dimension( :, : ), Allocatable :: diag

    Integer :: start
    Integer :: error    
    Integer :: n_spin
    Integer :: i

!!$    Do i = 1, Size( evecs )
!!$       Call vector_initialize( n_cart, n_harm, tol_warn, tol_project, overlap( i ), n_use, evecs( i ) )
!!$    End Do

    n_spin = Size( evecs )

    Allocate( work1( 1:n_spin ) )
    Call vector_create( work1, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_cart, n_cart, (/ ( 'Work1 ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Call vector_diagonalise( overlap, work1 )

    n_use = Count( work1( 1 )%evals >= tol_project )

    Call vector_create( evecs, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_cart, n_use, (/ ( 'evecs ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( work2( 1:n_spin ) )
    Call vector_create( work2, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_cart, n_use, (/ ( 'work2 ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    
    start = n_cart - n_use + 1
    

    Do i = 1, n_spin
       evecs( i )%evals( 1:n_use ) = work1( i )%evals( start:n_cart )
       work2( i )%evals = 1.0_wp / Sqrt( work1( i )%evals( start:n_cart ) )
    End Do
    
    Call matrix_copy_cols( start, n_cart, work1%evecs, &
         1    , n_use , work2%evecs )
    

    Allocate( diag( 1:n_use, 1:2 ) )
    Do i = 1, n_spin
       diag( :, i ) = work2( i )%evals
    End Do
    Call matrix_mult_by_diag( work2%evecs, diag, evecs%evecs )
    Deallocate( diag )
    
    Call vector_destroy( work2 )
    Deallocate( work2 )
    Call vector_destroy( work1 )
    Deallocate( work1 )

  End Subroutine vector_initialize_multi

  Subroutine vector_get_biggest_orb_pops( evecs, overlap, n_elec, biggest, ind_biggest, evals )

    Type( vector ), Dimension( :       )                  :: evecs
    Type( matrix ), Dimension( :       )                  :: overlap
    Integer       , Dimension( :       ), Intent( In    ) :: n_elec
    Real( wp )    , Dimension( :, :, : ), Intent(   Out ) :: biggest
    Integer       , Dimension( :, :, : ), Intent(   Out ) :: ind_biggest
    Real( wp )    , Dimension( :, :    ), Intent(   Out ) :: evals

    Type( matrix ), Dimension( : ), Allocatable :: s_times_q
    Type( matrix ), Dimension( : ), Allocatable :: orbital_pops

    Integer :: n, m
    Integer :: n_spin
    Integer :: distrib
    Integer :: error
    Integer :: i

    n_spin = Size( evecs )

    Call matrix_inquire( evecs( 1 )%evecs, global_n = n       )
    Call matrix_inquire( evecs( 1 )%evecs, global_m = m       )
    Call matrix_inquire( evecs( 1 )%evecs, distrib  = distrib )

    Allocate( s_times_q( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_GET_BIGGEST_ORB_POPS', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( s_times_q ) )
    Call matrix_create( n, m, s_times_q, (/ ( 's_times_q', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), distrib )

    Call matrix_multiply( 1.0_wp, overlap, evecs%evecs, 0.0_wp, s_times_q )

    Allocate( orbital_pops( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in VECTOR_GET_BIGGEST_ORB_POPS', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( orbital_pops ) )
    Call matrix_create( n, m, orbital_pops, (/ ( 'orbital_pops', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), distrib )

    Call matrix_direct_product( evecs%evecs, s_times_q, orbital_pops )

    Call matrix_get_col_biggest( orbital_pops, biggest, ind_biggest )

    Call matrix_destroy( orbital_pops )
    Call alloc_free( ALLOC_DERIVED, Size( orbital_pops ) )
    Deallocate( orbital_pops )

    Call matrix_destroy( s_times_q )
    Call alloc_free( ALLOC_DERIVED, Size( s_times_q ) )
    Deallocate( s_times_q )

    Do i = 1, n_spin
       evals( :, i ) = evecs( i )%evals
    End Do

  End Subroutine vector_get_biggest_orb_pops

End Module distributed_vectors
