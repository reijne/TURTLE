Module newscf_routines

  Public :: newscf
  Public :: initial_evecs

  Private

Contains

  Subroutine newscf( core, uhf, direct, enuclear, ehf, etot, sz, s2, ek, vir, &
       rdmat_work, prefac_work, focka_work, fockb_work, densa_work, densb_work, &
       accdi1 )
 

    Use newscf_numbers
    Use newscf_modules
    Use allocation
    Use distributed_matrices
    Use distributed_vectors

    Implicit None

    Real( wp ), Dimension( * ), Intent( InOut ) :: core
    Logical                   , Intent( In    ) :: uhf
    Logical                   , Intent( In    ) :: direct
    Real( wp )                , Intent(   Out ) :: enuclear
    Real( wp )                , Intent( InOut ) :: ehf
    Real( wp )                , Intent(   Out ) :: etot
    Real( wp )                , Intent(   Out ) :: sz
    Real( wp )                , Intent(   Out ) :: s2
    Real( wp )                , Intent(   Out ) :: ek
    Real( wp )                , Intent(   Out ) :: vir
    Real( wp ), Dimension( : ), Intent(   Out ) :: rdmat_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: prefac_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: focka_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: fockb_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: densa_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: densb_work
    Real( wp )                , Intent( In    ) :: accdi1

    Integer, Parameter :: orthog_freq = 5
!!$    Integer, Parameter :: front_print_evals = 32
    Integer, Parameter :: front_print_comps = 8

    Logical   , External :: opg_root
    Integer, External :: lensec
    Real( wp ), External :: nuclear_energy

    Type( vector ), Dimension( : ), Allocatable :: old_vectors_distrib
    Type( vector ), Dimension( : ), Allocatable :: best_vectors_distrib
    Type( vector ), Dimension( : ), Allocatable :: tmp_vectors_distrib
    Type( vector ), Dimension( : ), Allocatable :: vectors_distrib
    Type( vector ), Dimension( : ), Allocatable :: vectors_no_symm_distrib
    Type( vector ), Dimension( : ), Allocatable :: vectors_fock

    Type( matrix ), Dimension( : ), Allocatable :: Fock

    Type( matrix ), Dimension( : ), Allocatable :: hcore_distrib

    Type( matrix ), Dimension( : ), AllocAtable :: fock_distrib
    Type( matrix ), Dimension( : ), Allocatable :: fock_m1_distrib
    Type( matrix ), Dimension( : ), Allocatable :: fock_m2_distrib
    Type( matrix ), Dimension( : ), Allocatable :: tfock_distrib
    Type( matrix ), Dimension( : ), Allocatable :: tvmat_distrib
    Type( matrix ), Dimension( : ), Allocatable :: density_distrib
    Type( matrix ), Dimension( : ), Allocatable :: overlap_distrib
    Type( matrix ), Dimension( : ), Allocatable :: overlap_noadapt
    Type( matrix ), Dimension( : ), Allocatable :: scratch_2_distrib
    Type( matrix ), Dimension( : ), Allocatable :: scratch_1_distrib

    Real( wp ), Dimension( :, :, : ), Allocatable :: big_pops

    Real( wp ), Dimension( :, : ), Allocatable :: evals

    Real( wp ), Dimension( : ), Allocatable :: shift_cyc
    Real( wp ), Dimension( : ), Allocatable :: maxfok_array
    Real( wp ), Dimension( : ), Allocatable :: fok_lng_array
    Real( wp ), Dimension( : ), Allocatable :: fac1
    Real( wp ), Dimension( : ), Allocatable :: fac2
    Real( wp ), Dimension( : ), Allocatable :: len1
    Real( wp ), Dimension( : ), Allocatable :: len2
    Real( wp ), Dimension( : ), Allocatable :: iso

    Real( wp ) :: coef1
    Real( wp ) :: elow
    Real( wp ) :: ehf0
    Real( wp ) :: maxfok, rmsfok
    Real( wp ) :: maxden, rmsden
    Real( wp ) :: energy_change
    Real( wp ) :: tester, diis_tester
    Real( wp ) :: alpha_quad, beta_quad
    Real( wp ) :: ehf1, ehf2, edft, emermin, esmear, egap
    Real( wp ) :: diis_error 
    Real( wp ) :: extrap_fac, extrap_len1, extrap_len2_o, extrap_len2_n

    Integer, Dimension( :, :, : ), Allocatable :: ind_big_pops

    Integer, Dimension( : ), Allocatable :: occs

    Integer :: n_cart, n_harm
    Integer :: n_jobs, n_all_jobs
    Integer :: idum
    Integer :: l1, l2, l3
    Integer :: niter
    Integer :: maxcyc
    Integer :: iphase, niter_phase, new_phase
    Integer :: check_conv
    Integer :: nquad
    Integer :: diis_dimension
    Integer :: num_extrap_cyc
    Integer :: occs_start, occs_end, print_comps
    Integer :: error
    Integer :: i, j, k
    Integer :: linfo !jens

    Logical, Dimension( 1:maxphase ) :: o1
    Logical, Dimension( 1:maxphase ) :: o2
    Logical, Dimension( 1:maxphase ) :: o3
    Logical, Dimension( 1:maxphase ) :: o4
    Logical, Dimension( 1:maxphase ) :: o5

    Logical :: dft
    Logical :: use_extrap
    Logical :: converged
    Logical :: diis1, diis_on
    Logical :: extrap1, extrap_on
    Logical :: lock1
    Logical :: check

    Logical :: this_iter_best

    Character( Len = 5 ), Dimension( 1:2 ), Parameter :: spin_labels = &
         (/ 'Alpha', 'Beta ' /)

    Character( Len = 6 ) :: flags

    Integer :: iblkv, iblk, lenv, lenq, length1, length2
    Integer :: jjfile, lds, isect, ldsect, iacsct
    common/restri/jjfile(63),lds(508),isect(508),ldsect(508),iacsct(508)

!!$    If( opg_root() ) Then
!!$       Call _dump_allocated( %val( 0 ) )
!!$    End If

    Call start_time_period( TP_NEWSCF )
    Call start_time_period( TP_F90_START )

    Call alloc_set_verbosity( .False. )

    n_cart = num
    n_harm = newbas0

    Allocate( iso( 1:nw196( 5 ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call rdedx( iso, nw196( 5 ), ibl196( 5 ), idaf )
    Call put_iso( nw196( 5 ), iso )
    Deallocate( iso )


    length1 = lensec( n_cart + 1 )
    lenv = length1 + lensec( n_cart * ( n_cart + 1 ) / 2 )
    lenq = lenv + lensec( n_cart * n_cart )
    Call secput( isecqm, 167, lenq, iblkv )
    ibl3qs = iblkv + lenv
    ! vectors section
    length1 = lensec( mach( 8 ) )
    lenv = lensec( n_cart * n_cart )
    length2 = lensec( mach( 9 ) )
    Call secput( mouta, 3, length1 + lenv + length2 + 1, iblk )
    ibl3qa = iblk + length2 + length1 + 1
    If( uhf ) Then
       Call secput( moutb, 3, length1 + lenv + length2 + 1, iblk )
       ibl3qb = iblk + length2 + length1 + 1
    End If

    If( uhf ) Then
       Allocate( occs( 1:2 ) )
       occs = (/ na, nb /)
    Else
       Allocate( occs( 1:1 ) )
       occs = (/ na /)
    End If

!!$    Call mp_disableintr( error )

!    Call alloc_report

    Call get_n_jobs( uhf, n_jobs, n_all_jobs )
    Call matrix_comms_init( n_cart, mpi_comm_workers, blacs_context, n_jobs, 1 )

    Allocate( tfock_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( tfock_distrib ) )
    Call matrix_create( n_harm, n_harm, tfock_distrib, (/ ( 'tfock distrib' // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( vectors_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_distrib ) )
    Call vector_create( vectors_distrib, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_cart, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )

    Call vector_tdown_setup( vectors_distrib, ctran, ntran, itran, ilifc, otran )

    ! Allocate and create Fock, transformed fock, transformed vector and
    ! density matrices 
    Allocate( fock_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( fock_distrib ) )
    Allocate( tvmat_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( tvmat_distrib ) )
    Call matrix_create( n_cart, n_cart, fock_distrib, (/ ( 'tfock distrib' // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )
    Call matrix_create( tfock_distrib, tvmat_distrib )

    ! Scratch space for temporary fock matrices
    Allocate( scratch_1_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( scratch_1_distrib ) )
    Allocate( scratch_2_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( scratch_2_distrib ) )
    Call matrix_create( fock_distrib, scratch_1_distrib )
    Call matrix_create( fock_distrib, scratch_2_distrib )

    dft = CD_active()
    If( dft ) Then
       idum = CD_update_geom( c )
    End If

    l1 = num
    l2 = l1 * ( l1 + 1 ) / 2
    l3 = l1 * l1

    sz  = 0.0_wp
    s2  = 0.0_wp
    ek  = 0.0_wp
    vir = 0.0_wp

    num_extrap_cyc = 0

    Call default_conv( uhf )
    Call print_conv  ( uhf )

    ! Allocate DIIS
    Call vector_diis_initialize( n_cart, n_harm, uhf )

    enuclear = nuclear_energy()

    ! Allocate Vectors

    Allocate( old_vectors_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( old_vectors_distrib ) )
    Allocate( best_vectors_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( best_vectors_distrib ) )
    Allocate( tmp_vectors_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( tmp_vectors_distrib ) )
    Allocate( vectors_no_symm_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_no_symm_distrib ) )
    Allocate( vectors_fock( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_fock ) )
    Call vector_create( old_vectors_distrib, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_cart, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )
    Call vector_create( best_vectors_distrib, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_cart, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )
    Call vector_create( tmp_vectors_distrib, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_cart, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )
    Call vector_create( vectors_no_symm_distrib, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_cart, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )
    Call vector_create( vectors_fock, (/ ( 0, i = 1, n_jobs ) /), (/ ( i, i = 1, n_jobs ) /), &
         n_harm, n_harm, (/ ( 'Vectors '   // spin_labels( i ), i = 1, n_jobs ) /), &
         (/ ( MATRIX_REAL, i = 1, n_jobs ) /), MATRIX_DISTRIB_UHF_K )

    ! AO overlap matrix
    Allocate( overlap_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( overlap_distrib ) )
    Call matrix_create( fock_distrib, overlap_distrib )

    Allocate( overlap_noadapt( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( overlap_noadapt ) )
    Call matrix_create( fock_distrib, overlap_noadapt )

    Allocate( hcore_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( hcore_distrib ) )
    Call matrix_create( fock_distrib, hcore_distrib )

    Call read_initial( uhf, vectors_distrib, overlap_distrib, hcore_distrib, overlap_noadapt )

!!$    ! Load guess vectors
!!$    Call vector_read( vectors_distrib, (/ ibl3qa, ibl3qb /), (/ ibl3ea, ibl3eb /), idaf ) 
!!$
!!$    !HACKY !!!!!!!!!!!!
!!$    Call matrix_read_triangle( overlap_distrib( 1 ), ibl7st, num8 )
!!$    If( uhf ) Then
!!$       Call matrix_read_triangle( overlap_distrib( 2 ), ibl7st, num8 )
!!$    End If
!!$
!!$    ! Core Hamiltonian
!!$    ! HACKY !!!!!!!!!!!!
!!$    Call matrix_read_triangle( hcore_distrib( 1 ), ibl7f, num8 )
!!$    If( uhf ) Then
!!$       Call matrix_read_triangle( hcore_distrib( 2 ), ibl7f, num8 )
!!$    End If

    ! Need to orthog 
    Call start_time_period( TP_F90_ORTHOG )
    Call vector_orthogonalize( vectors_distrib, overlap_distrib, .True. )
    Call end_time_period( TP_F90_ORTHOG )

    ! Set up intial density matrix
    If( smear( 1 ) == SMEAR_OFF ) Then
       Call vector_assign_by_energy( vectors_distrib, occs )
    Else
       flags( 6:6 ) = 'F'
       If( smear( 1 ) == SMEAR_ENERGY ) Then
          esmear = new_esmear_start( 1 )
       Else
          Call vector_get_gap( vectors_distrib, occs, egap )
          esmear = egap * egap_scale( 1 )
       End If
       Call vector_assign_by_fermi( vectors_distrib, esmear, occs, emermin )
    End If
    Allocate( density_distrib( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( density_distrib ) )
    Call matrix_create( fock_distrib, density_distrib )
    Call vector_tdown( vectors_distrib, vectors_no_symm_distrib, ctran, newbas0, otran )
    Call vector_make_density( density_distrib, vectors_no_symm_distrib )


    ! Start of main SCF loop.

    niter       = 0
    niter_phase = 0
    maxcyc      = maxcycp
    ehf0        = 0.0_wp
    converged   = .False.
    iphase      = 1
    elow        = 1.0_wp
    coef1       = 0.0_wp

    ! Set up extrapolation stuff if required
    use_extrap = Any( extrap( 1:maxphase ) )
    If( use_extrap ) Then
       Allocate( fock_m1_distrib( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( fock_m1_distrib ) )
       Allocate( fock_m2_distrib( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( fock_m2_distrib ) )
       Call matrix_create( fock_distrib, fock_m1_distrib )
       Call matrix_create( fock_distrib, fock_m2_distrib )
       Allocate( fac1( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( fac1 ) )
       Allocate( fac2( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( fac2 ) )
       Allocate( len1( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( len1 ) )
       Allocate( len2( 1:n_jobs ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       Call alloc_add( ALLOC_REAL, Size( len2 ) )
    End If

    Allocate( shift_cyc    ( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( shift_cyc ) )
    Allocate( maxfok_array ( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( maxfok_array ) )
    Allocate( fok_lng_array( 1:n_jobs ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( fok_lng_array ) )

    ! No fitting at present as memory stuff a bit tricky.

    ! Header for SCF output
    If( .Not. oscfprint( PR_FULL ) .And. opg_root() ) Then
       Write( iwr, "( 1x, 'niter', 5x, 'Energy', 4x, 'Energy Change', 4x, 'Tester' )" )
    End If

    ! Need to zero the fock matrix so the delta calcs work on the first cycle.
    Call matrix_assign( fock_distrib, 0.0_wp )

    Call end_time_period( TP_F90_START )

    Call start_time_period( TP_F90_SCF )

    If( direct ) Then
       dlntol = tolitr(3)
       dlnmxd = 0.0_wp
       Call secput( isect( 471 ), 171, 2 * lensec( ikyp( nshell ) ), ibl171 )
       Call rdmake( prefac_work )
       Call wrt3 ( prefac_work, ikyp( nshell ), ibl171, idaf )
    End If


    Do While( .Not. converged .And. niter < maxcyc .And. iphase /= -1 )

       ! Update cycle counters
       niter = niter + 1
       niter_phase = niter_phase + 1

       ! flags for output
       flags = ' '

       ! Save previous energy
       ehf0 = ehf

       ! Set convergance control parameters
       shift_cyc( 1:n_jobs ) = shift( iphase, 1:n_jobs )
       diis1                 = diis( iphase )
       extrap1               = extrap( iphase )
       lock1                 = lock_vec( iphase )

       ! Save current fock matrices
       Call matrix_copy( fock_distrib, scratch_1_distrib )

       ! Save current Vectors
       Call vector_copy( vectors_distrib, tmp_vectors_distrib )

       ! Need to call fock_build here - probs: need all densities
       ! for all possible spin/k-point cases - THINK
       ! HACK - just call modified old routine at present.
       Call start_time_period( TP_F90_BUILD )
       If( direct )Then
!!$          If(delta)Then
!!$             ! -- note delta density not implemented
!!$          End If
          Call start_time_period( TP_F90_RDMAT )
          Call make_rdmat( rdmat_work, uhf, density_distrib, densa_work, densb_work )
          Call end_time_period( TP_F90_RDMAT )
          dlntol = tolitr(3)
          dlnmxd = Maxval( rdmat_work )
          If( opg_root() ) Then
             Write( iwr, "( 5x, 'dlntol =  ', f15.10 / 5x, 25( '*' ) )" ) dlntol
          End If
       End If

       Call fock_build( density_distrib, direct, &
            fock_distrib, hcore_distrib, &
            uhf, dft, ehf, ehf1, ehf2, edft, &
            nquad, alpha_quad, beta_quad, &
            rdmat_work, prefac_work, focka_work, fockb_work, densa_work, densb_work, &
            core )
       Call end_time_period( TP_F90_BUILD )

       If( smear( iphase ) /= SMEAR_OFF ) Then
          ehf = ehf + emermin
       End If

       ! Check difference fock matrices
       Call start_time_period( TP_F90_DELTA_EVAL )
       Call matrix_copy  ( fock_distrib, scratch_2_distrib )
       Call matrix_daxpy ( -1.0_wp, scratch_1_distrib, scratch_2_distrib )
       Call matrix_absmax( scratch_2_distrib, maxfok_array )
       maxfok = Maxval( maxfok_array )
       Call matrix_length( scratch_2_distrib, fok_lng_array )
       rmsfok = Sqrt( Sum( fok_lng_array ) ) / ( n_all_jobs * l3 )
       Call end_time_period( TP_F90_DELTA_EVAL )

       ! Option to print fock matrices not implemented

       energy_change = ehf - ehf0

       ! Save the current vectors
       Call vector_copy( vectors_distrib, old_vectors_distrib )

       ! Transform to harmonic/symm adapted
       Call vector_tdown( vectors_distrib, vectors_no_symm_distrib, ctran, newbas0, otran )

       ! DIIS
       Call start_time_period( TP_F90_DIIS )
       Call vector_diis_solve( fock_distrib, vectors_no_symm_distrib, overlap_distrib, &
            diis1, diis_on, diis_dimension, &
            tester, diis_error, accdi1 )
       If( diis_on ) Then
          flags( 1:1 ) = 'd'
       End If
       Call end_time_period( TP_F90_DIIS )

       ! Transform Fock to MO basis
       Call start_time_period( TP_F90_SIMIL )
       Call vector_similarity( fock_distrib, vectors_no_symm_distrib, tfock_distrib )
       Call end_time_period( TP_F90_SIMIL )

       ! Fock extrapolation
       extrap_on = .False.
       If( extrap1 ) Then
          ! First check unextrapolated tester
          tester = vector_check_ov( tfock_distrib, vectors_distrib )
          ! Now if enough cycles try extrapolation
          If( num_extrap_cyc == 2 ) Then
             ! Calculate 'vectors' between points
             Call matrix_copy ( fock_m1_distrib, scratch_2_distrib )
             Call matrix_daxpy( -1.0_wp, fock_m2_distrib, scratch_2_distrib )
             Call matrix_copy ( fock_distrib, scratch_1_distrib )
             Call matrix_daxpy( -1.0_wp, fock_m1_distrib, scratch_1_distrib )
             ! Now calculate angle between the 'vectors'
             Call matrix_length( scratch_2_distrib, len2 )
             Call matrix_length( scratch_1_distrib, len1 )
             Call matrix_dot_product( scratch_2_distrib, scratch_1_distrib, fac1 )
             fac1 = fac1 / ( len1 * len2 )
             fac2 = len1 / len2
             ! If going in a straight enough line extrapolate
             extrap_fac    = fac1( 1 )
             extrap_len1   = len2( 1 )
             extrap_len2_o = len1( 1 )
             If( fac1( 1 ) > extrap_tol( iphase ) .And. &
                 fac2( 1 ) > 0.7_wp               .And. &
                 fac2( 1 ) < 1.2_wp + coef1 ) Then
                coef1 = extrap_coef( iphase )
                Call matrix_daxpy( coef1, scratch_1_distrib, fock_distrib )
                Call matrix_copy ( fock_distrib, scratch_1_distrib )
                Call matrix_daxpy( -1.0_wp, fock_m1_distrib, scratch_1_distrib )
                Call matrix_length( scratch_1_distrib, len1 )
                Call matrix_dot_product( scratch_2_distrib, scratch_1_distrib, fac1 )
                fac1 = fac1 / ( len1 * len2 )
                extrap_len2_n = len1( 1 )
                flags( 2:2 ) = 'E'
                extrap_on    = .True.
                Call start_time_period( TP_F90_SIMIL )
                Call vector_similarity( fock_distrib, vectors_no_symm_distrib, tfock_distrib )
                Call end_time_period( TP_F90_SIMIL )
             End If
          End If
       End If
       If( .Not. extrap_on ) Then
          coef1 = 0.0_wp
       End If

       num_extrap_cyc = Min( num_extrap_cyc + 1, 2 )

       ! Calculate tester if applicable
       Call start_time_period( TP_F90_TESTER_EVAL )
       If( extrap_on ) Then
       Else If( diis_on ) Then
          diis_tester = vector_check_ov( tfock_distrib, vectors_distrib )
       Else
          tester = vector_check_ov( tfock_distrib, vectors_distrib )
       End If
       If( smear( iphase ) /= SMEAR_OFF ) Then
          If( smear( iphase ) == SMEAR_ENERGY ) Then
             esmear = Max( new_esmear_final( iphase ), Min( tester, esmear ) )
          Else
             esmear = egap * egap_scale( iphase )
          End If
       End If
       Call end_time_period( TP_F90_TESTER_EVAL )

       ! Level Shift
       Call start_time_period( TP_F90_LEV_SHIFT )
       Call vector_level_shift( tfock_distrib, vectors_distrib, shift_cyc )
       If( Any( shift_cyc > 0.0_wp ) ) Then
          flags( 3:3 ) = 'S'
       End If
       Call end_time_period( TP_F90_LEV_SHIFT )

       ! Diagonalise and back tranform the vectors ( at damn last !! )
       Call start_time_period( TP_F90_DIAG )
       Call start_time_period( TP_DIAG )
       Call vector_diagonalise( tfock_distrib, vectors_fock )
       Call end_time_period( TP_DIAG )

       Call end_time_period( TP_F90_DIAG )
       Call start_time_period( TP_F90_BACK )
       Call vector_back_transform( vectors_fock, old_vectors_distrib, vectors_distrib )
       Call end_time_period( TP_F90_BACK )

       ! Printing of vectors not implemented

       ! Assign occupations
       Call start_time_period( TP_F90_ASSIGN )
       If( smear( iphase ) /= SMEAR_OFF ) Then
          Call vector_assign_by_fermi( vectors_distrib, esmear, occs, emermin )
          flags( 6:6 ) = 'F'
       Else If( lock1 ) Then
          Call vector_assign_by_overlap( vectors_distrib, old_vectors_distrib, occs, &
               check )
          If( .Not. check ) Then
             Call vector_assign_by_energy( vectors_distrib, occs )
          End If
          flags( 4:4 ) = 'L'
       Else
          Call vector_assign_by_energy( vectors_distrib, occs )
       End If
       Call end_time_period( TP_F90_ASSIGN )

       If( Mod( niter, orthog_freq ) == 0 ) Then
          Call start_time_period( TP_F90_ORTHOG )
          Call vector_orthogonalize( vectors_distrib, overlap_distrib, .True. )
          Call end_time_period( TP_F90_ORTHOG )
       End If

       ! Create new density matrix
       ! New density now in SCRATCH_1
       Call start_time_period( TP_F90_TDOWN )
       Call vector_tdown( vectors_distrib, vectors_no_symm_distrib, ctran, newbas0, otran )
       Call end_time_period( TP_F90_TDOWN )

       Call start_time_period( TP_F90_MAKE_DENS )
       Call vector_make_density( scratch_1_distrib, vectors_no_symm_distrib )
       Call end_time_period( TP_F90_MAKE_DENS )

       Call start_time_period( TP_F90_DELTA_EVAL )
       ! Check difference density
       Call matrix_copy( scratch_1_distrib, scratch_2_distrib )
       Call matrix_daxpy ( -1.0_wp, density_distrib, scratch_2_distrib )
       Call matrix_absmax( scratch_2_distrib, maxfok_array )
       maxden = Maxval( maxfok_array )
       Call matrix_length( scratch_2_distrib, fok_lng_array )
       rmsden = Sqrt( Sum( fok_lng_array ) ) / ( n_all_jobs * l3 )
       Call matrix_copy  ( scratch_1_distrib, density_distrib )
       Call end_time_period( TP_F90_DELTA_EVAL )

       ! Printing of density not implemented

       ! Check for convergance, or switch of converg params
       new_phase = check_conv( uhf, iphase, tester, energy_change, &
            niter_phase, niter, o1, o2, o3, o4, o5 )
       ! The energy
       etot = ehf + enuclear

       ! If using extrapolation store previous fock matrices
       If( use_extrap ) Then
          If( niter >=  2 ) Then
             Call matrix_copy( fock_m1_distrib, fock_m2_distrib )
          End If
          If( niter >=  1 ) Then
             Call matrix_copy( fock_distrib, fock_m1_distrib )
          End If
       End If

       ! Frontier orbital analysis
       If( oscfprint( PR_FRONTIER ) ) Then
          print_comps = Min( n_cart, front_print_comps )
          Allocate( big_pops( 1:print_comps, 1:n_harm, 1:Size( occs ) ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
          End If
          Call alloc_add( ALLOC_REAL, Size( big_pops ) )
          Allocate( ind_big_pops( 1:print_comps, 1:n_harm, 1:Size( occs ) ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
          End If
          Call alloc_add( ALLOC_INTEGER, Size( ind_big_pops ) )
          Allocate( evals( 1:n_harm, 1:size( occs ) ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
          End If
          Call alloc_add( ALLOC_REAL, Size( evals ) )
          Call vector_get_biggest_orb_pops( vectors_no_symm_distrib, overlap_noadapt, &
               occs, big_pops, ind_big_pops, evals )
          If( opg_root() ) Then
             Write( iwr, * )
             Write( iwr, * ) 'Frontier orbital population analysis:'
             Write( iwr, * )
             If( .Not. uhf ) Then
                Write( iwr, '( 1x, "Number of electrons: ", i5 )' ) occs( 1 )
             Else
                Write( iwr, '( 1x, 2( "Number of ", a5, " electrons: ", i5, 4x ) )' ) &
                     ( spin_labels( j ), occs( j ), j = 1, 2 )
             End If
             Write( iwr, * )
             occs_start = Max( Minval( occs - front_print_evals / 2 + 1 ), 1 )
             occs_end   = Min( Maxval( occs + front_print_evals / 2 ), n_harm )
             If( uhf ) Then
                Write( iwr, '( 1x, 5x, 2( 16x, a5, 17x ) )' ) &
                     ( spin_labels( j ), j = 1, 2 )
             End If
             Write( iwr, '( 1x, 1x, "M.O.", 2( 1x, a9, 1x, a6, 1x, a10, 1x, a5, 4x ) )' ) &
                  ( '    Eval ', '  Pop ', 'Bfunc type', 'Bfunc numb', &
                           j = 1, Size( occs ) )
             Write( iwr, * )
             Do k = occs_start, occs_end
                Do i = 1, print_comps
                   If( i == 1 ) Then
                      Write( iwr, '( 1x, i5, 2( 1x, f9.3, 1x, f6.3, 1x, a10, 1x, i5, 4x ) )' ) &
                           k, ( evals( k, j ), big_pops( i, k, j ), &
                           Adjustr( zbflab( ind_big_pops( i, k, j ) ) ), ind_big_pops( i, k, j ), &
                           j = 1, Size( occs ) )
                   Else
                      Write( iwr, '( 1x, 5x, 2( 1x, 9x  , 1x, f6.3, 1x, a10, 1x, i5, 4x ) )' ) &
                           (                   big_pops( i, k, j ), &
                           Adjustr( zbflab( ind_big_pops( i, k, j ) ) ), ind_big_pops( i, k, j ), &
                           j = 1, Size( occs ) )
                   End If
                End Do
                Write( iwr, * )
             End Do
          End If
          Call alloc_free( ALLOC_REAL, Size( evals ) )
          Deallocate( evals )
          Call alloc_free( ALLOC_INTEGER, Size( ind_big_pops ) )
          Deallocate( ind_big_pops )
          Call alloc_free( ALLOC_REAL, Size( big_pops ) )
          Deallocate( big_pops )
       End If

       ! Print intermediate results
!!$       Call print_intermediate
       If( opg_root() ) Then
          If( oscfprint( PR_FULL ) ) Then
             Write( iwr, '( /1x,                                 &
                  &        "Cycle ",i4,1x,                       &
                  &        "E       = ",f16.8,1x,                &
                  &        "dE      = ",g14.6,3x,                &
                  &        "Tester  = ",f16.8,3x,i4,1x,a,/,12x,  &
                  &        "EHF(1e) = ",f16.8,1x,                &
                  &        "EHF(2e) = ",f16.8,1x,                &
                  &        "EDFT    = ",f16.8 )' )               &
                  niter, etot, energy_change, tester, iphase, flags, ehf1, ehf2, edft
             If( dft ) Then
                If( uhf ) Then
                   Write( iwr, '( 12x, &
                        &        "Nquad    = ", i14, 1x,  &
                        &        "Nelec alpha = ", f14.8, &
                        &        " beta = ", f14.8 )' )   &
                        nquad, alpha_quad, beta_quad
                Else
                   Write( iwr, '( 12x, &
                        &        "Nquad    = ", i14, 1x,    & 
                        &        "Nelec    = ", f14.8 ) ' ) &
                        nquad, alpha_quad
                End If
             End If
             Write( iwr, '( 12x, "Shift   = ", f10.4, 2x, f10.4 )' ) shift_cyc 
             If( diis_on ) Then
                Write( iwr, '( 12x, "Diis On Dim = ", i3,                 &
                     &              " Error = ",f15.9,                    &
                     &              " Extrapolated tester = ",f14.8 )' )  &
                     diis_dimension, diis_error, diis_tester
             End If
             if(extrap1 .and. niter > 2 )then
                if(extrap_on)then
                   write(iwr, "(12x,'Extrapolation  test=',f14.8,' scal ',f12.4, &
                        &        ' len1',f12.4, &
                        &        ' len2(old) = ',f12.4,' len2(ext) ',f12.4)" ) &
                        extrap_fac,coef1, &
                        extrap_len1,                 &
                        extrap_len2_o, extrap_len2_n 
                else
                   write(iwr,"(12x,'Extrapolation  test=',f14.8,' suppressed',10x, &
                        &        ' len1 ',f12.4, &
                        &        ' len2 ',f12.4)" ) &
                        extrap_fac,extrap_len1,extrap_len2_o
                endif
             endif
             If( smear( iphase ) /= SMEAR_OFF ) Then
                If( smear( iphase ) == SMEAR_ENERGY ) Then
                   ! Magic number below converts Hartrees to Kelvin
                   Write( iwr, '( 12x, "Smearing Energy = ", f8.6, " ( ", f10.2, "K )", 2x, &
                        &   "Mermin Energy = ", f16.8 )' ) esmear, &
                        esmear * 315716.93_wp, emermin
                Else
                   Write( iwr, '( 12x, "Smearing Energy = ", f8.6, " ( ", f10.2, "K )", 2x, &
                        &   "Band Gap = ", f9.5, 2x,                                        &
                        &   "Mermin Energy = ", f16.8 )' ) esmear, &
                        esmear * 315716.93_wp, egap, emermin
                End If
             End If
             Write( iwr, '( 12x,                               &
                  &         "Changes - Density Max = ", f14.8, &
                  &         " RMS = ", f14.8,                  &
                  &         " Fock Max = ", f14.8,             &
                  &         " RMS = ", f14.8 )' )              &
                  maxden, rmsden, maxfok, rmsfok
             Do i = 1, nnext( iphase )
                If( nextphase( iphase, i ) == 0 ) Then
                   Write( iwr, '(  12x, "                                   &
                        &          Criteria for convergence: tester: ", l1, &
                        &               ", dE: ", l1,                       &
                        &               ", |dE|: ", l1,                     &
                        &               ", Ncyc: ", l1 )' )                 &
                        o1( i ), o2( i ), o3( i ), o4( i )
                Else
                   Write( iwr, '( 12x,                    &
                        &      "Criteria for phase",  i2, &
                        &      ": tester: ", l1,          &
                        &      ", dE: ", l1,              &
                        &      ", |dE|: ", l1,            &
                        &      ", Ncyc: ", l1 )' )        &
                        nextphase( iphase, i ), o1( i ), o2( i ), o3( i ), o4( i )
                End If
             End Do
          Else
             Write( iwr, '( 1x, i4, 2x, f14.8, 2x, g12.6, 2x, f12.6 )' ) &
                  niter, ehf + enuclear, energy_change, tester
          End If
       End If
       
       ! Update converg options
!!$       Call converg_update
       this_iter_best = iphase /= -1 .And. ehf < elow

       If( this_iter_best ) Then
          Call vector_save_orbitals( tmp_vectors_distrib, density_distrib )
          Call vector_copy  ( tmp_vectors_distrib, best_vectors_distrib )
          elow = ehf
       End If
       
       If( iphase /= new_phase ) Then
          
          If( opg_root() ) Then
             Write( iwr, * )
             If( new_phase == -1 ) Then
                Write( iwr, "( '**** Calculation aborted ' )" )
             Else If( new_phase /= 0 ) Then
                Write( iwr, "( '**** Switch to convergence phase ', i0, /, / )" ) new_phase
             End If
          End If
         
          niter_phase = 0
          iphase      = new_phase
          converged   = iphase == 0
         
         ! Restore best vecs if required
          If( iphase > 0 ) Then
             If( restore_vec( iphase ) .And. .Not. this_iter_best ) Then
                Call vector_copy ( best_vectors_distrib, vectors_distrib )
                Call start_time_period( TP_F90_TDOWN )
                Call vector_tdown( vectors_distrib, vectors_no_symm_distrib, ctran, newbas0, otran )
                Call end_time_period( TP_F90_TDOWN )
                Call start_time_period( TP_F90_MAKE_DENS )
                Call vector_make_density( density_distrib, vectors_no_symm_distrib )
                Call end_time_period( TP_F90_MAKE_DENS )
             End If
          End If
          
          ! Set up smearing based on gap if required
          If( smear( new_phase ) == SMEAR_GAP ) Then
             Call vector_get_gap( vectors_distrib, occs, egap )
          End If
          
          ! Reset DIIS if required
          If( iphase > 0 ) Then
             If( new_diis( iphase ) ) Then
                If( opg_root() ) Then
                   Write( iwr, "( 'Resetting DIIS' )" )
                End If
                Call vector_diis_reset
             End If
          End If
          
          ! Reset extrapolation
          if(iphase .gt. 0 .and. new_extrap(iphase))then
             write(iwr,*)'Resetting Extrap'
             num_extrap_cyc = 0
          endif
          
       End If
       
    End Do

!    Call alloc_report

    Call end_time_period( TP_F90_SCF )

    Call start_time_period( TP_F90_END )

    ! No dealloc of fitting for DFT 'cos no fitting

    ! So tell the world what happened !
    If( converged ) Then
       If( opg_root() ) Then
          Write( iwr, * ) '*************'
          Write( iwr, * ) 'SCF converged'
          Write( iwr, * ) '*************'
          call strtrm(pchange_info,linfo)
          If ( linfo .gt. 1 ) Write( iwr, 8999 ) pchange_info ! Write out user-defined information
          Write( iwr, 9000 ) niter, ehf, edft, enuclear, etot
       End If
8999   Format( /'Convergence Info: ',a72 )
9000   Format(/10x,14('-')/10x,'final energies after',i4,' cycles'/, &
            &          10x,14('-')/                                  &
            &          10x,'hartree-fock energy         ',f18.10/    &
            &          10x,'exchange-correlation energy ',f18.10/    &  
            &          10x,'nuclear repulsion energy    ',f18.10/    &
            &          10x,'total energy                ',f18.10)

       Call vector_save_orbitals( vectors_distrib, density_distrib )
    Else
       If( opg_root() ) Then
          write( iwr, * ) '********************'
          write( iwr, * ) 'SCF did not converge'
          write( iwr, * ) '********************'
       End If
       
!!$       If( .Not. omaxcyc ) Then
!!$          If( hardfail ) Then
!!$             Call gamerr( 'SCF convergence failure', ERR_NO_CODE, ERR_UNLUCKY, ERR_SYNC, ERR_NO_SYS )
!!$          End If
!!$       Else
          Write( iwr, * ) 'Too many iterations'
!!$       End If
       irest = 3
       nindmx = 0
    End If

    ! Destroy in opposite order to create - good practice

    Call alloc_free( ALLOC_REAL, Size( fok_lng_array ) )
    Deallocate( fok_lng_array )
    Call alloc_free( ALLOC_REAL, Size( maxfok_array ) )
    Deallocate( maxfok_array  )
    Call alloc_free( ALLOC_REAL, Size( shift_cyc ) )
    Deallocate( shift_cyc     )

    If( use_extrap ) Then
       Call alloc_free( ALLOC_REAL, Size( len2 ) )
       Deallocate( len2 )
       Call alloc_free( ALLOC_REAL, Size( len1 ) )
       Deallocate( len1 )
       Call alloc_free( ALLOC_REAL, Size( fac2 ) )
       Deallocate( fac2 )
       Call alloc_free( ALLOC_REAL, Size( fac1 ) )
       Deallocate( fac1 )
       Call alloc_free( ALLOC_DERIVED, Size( fock_m2_distrib ) )
       Call matrix_destroy( fock_m2_distrib )
       Call alloc_free( ALLOC_DERIVED, Size( fock_m1_distrib ) )
       Call matrix_destroy( fock_m1_distrib )
    End If

    Call matrix_destroy( density_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( density_distrib ) )
    Deallocate( density_distrib )

    Call matrix_destroy( scratch_2_distrib )
    Call matrix_destroy( scratch_1_distrib )

    Call alloc_free( ALLOC_DERIVED, Size( scratch_2_distrib ) )
    Deallocate( scratch_2_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( scratch_1_distrib ) )
    Deallocate( scratch_1_distrib )

    Call matrix_destroy( tvmat_distrib )
    Call matrix_destroy( fock_distrib  )

    Call alloc_free( ALLOC_DERIVED, Size( tvmat_distrib ) )
    Deallocate( tvmat_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( fock_distrib ) )
    Deallocate( fock_distrib )

    Call matrix_destroy( hcore_distrib   )
    Call alloc_free( ALLOC_DERIVED, Size( hcore_distrib ) )
    Deallocate( hcore_distrib )

    Call matrix_destroy( overlap_noadapt )
    Call alloc_free( ALLOC_DERIVED, Size( overlap_noadapt ) )
    Deallocate( overlap_noadapt )
    Call matrix_destroy( overlap_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( overlap_distrib ) )
    Deallocate( overlap_distrib )

    Call vector_destroy( vectors_fock  )
    Call vector_destroy( vectors_no_symm_distrib  )
    Call vector_destroy( tmp_vectors_distrib  )
    Call vector_destroy( best_vectors_distrib )
    Call vector_destroy( old_vectors_distrib  )

    Call alloc_free( ALLOC_DERIVED, Size( vectors_fock ) )
    Deallocate( vectors_fock )
    Call alloc_free( ALLOC_DERIVED, Size( vectors_no_symm_distrib ) )
    Deallocate( vectors_no_symm_distrib  )
    Call alloc_free( ALLOC_DERIVED, Size( tmp_vectors_distrib ) )
    Deallocate( tmp_vectors_distrib  )
    Call alloc_free( ALLOC_DERIVED, Size( best_vectors_distrib ) )
    Deallocate( best_vectors_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( old_vectors_distrib ) )
    Deallocate( old_vectors_distrib  )

    Call vector_diis_free_scratch

    Call matrix_tdown_free

    Call vector_destroy( vectors_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( vectors_distrib ) )
    Deallocate( vectors_distrib )

    Call matrix_destroy( tfock_distrib )
    Call alloc_free( ALLOC_DERIVED, Size( tfock_distrib ) )
    Deallocate( tfock_distrib )

    Call matrix_comms_finalize

!    Call alloc_report

!!$    Call mp_enableintr( error )

    Deallocate( occs )

    Call end_time_period( TP_F90_END )
    Call end_time_period( TP_NEWSCF )

  Contains
    
    Subroutine print_intermediate

      Real( wp ) :: tmp

      If( opg_root() ) Then
         If( oscfprint( PR_FULL ) ) Then
            Write( iwr, '( /1x,                                 &
                 &        "Cycle ",i4,1x,                       &
                 &        "E       = ",f17.8,1x,                &
                 &        "dE      = ",g14.6,3x,                &
                 &        "Tester  = ",f16.8,3x,i4,1x,a5,/,12x, &
                 &        "EHF(1e) = ",f17.8,1x,                &
                 &        "EHF(2e) = ",f16.8,1x,                &
                 &        "EDFT    = ",f16.8 )' )               &
                 niter, etot, energy_change, tester, iphase, flags, ehf1, ehf2, edft
            If( dft ) Then
               If( uhf ) Then
                  Write( iwr, '( 12x, &
                       &        "Nquad    = ", i14, 1x,  &
                       &        "Nelec alpha = ", f14.8, &
                       &        " beta = ", f14.8 )' )   &
                       nquad, alpha_quad, beta_quad
               Else
                  Write( iwr, '( 12x, &
                       &        "Nquad    = ", i14, 1x,    & 
                       &        "Nelec    = ", f14.8 ) ' ) &
                       nquad, alpha_quad
               End If
            End If
            Write( iwr, '( 12x, "Shift   = ", f10.4, 2x, f10.4 )' ) shift_cyc 
            If( diis_on ) Then
               Write( iwr, '( 12x, "Diis On Dim = ", i3,                 &
                    &              " Error = ",f15.9,                    &
                    &              " Extrapolated tester = ",f14.8 )' )  &
                    diis_dimension, diis_error, diis_tester
            End If
            if(extrap1 .and. niter > 2 )then
               if(extrap_on)then
                  write(iwr, "(12x,'Extrapolation  test=',f14.8,' scal ',f12.4, &
                       &        ' len1',f12.4, &
                       &        ' len2(old) = ',f12.4,' len2(ext) ',f12.4)" ) &
                       extrap_fac,coef1, &
                       extrap_len1,                 &
                       extrap_len2_o, extrap_len2_n 
               else
                  write(iwr,"(12x,'Extrapolation  test=',f14.8,' suppressed',10x, &
                       &        ' len1 ',f12.4, &
                       &        ' len2 ',f12.4)" ) &
                       extrap_fac,extrap_len1,extrap_len2_o
               endif
            endif
            Write( iwr, '( 12x,                               &
                 &         "Changes - Density Max = ", f14.8, &
                 &         " RMS = ", f14.8,                  &
                 &         " Fock Max = ", f14.8,             &
                 &         " RMS = ", f14.8 )' )              &
                 maxden, rmsden, maxfok, rmsfok
            Do i = 1, nnext( iphase )
               If( nextphase( iphase, i ) == 0 ) Then
                  Write( iwr, '(  12x,                                     &
                       &         "Criteria for convergence: tester: ", l1, &
                       &               ", dE: ", l1,                       &
                       &               ", |dE|: ", l1,                     &
                       &               ", Ncyc: ", l1,                     &
                       &               ", Totcyc: ", l1 )' )               &
                       o1( i ), o2( i ), o3( i ), o4( i ), o5( i )
               Else
                  Write( iwr, '( 12x,                    &
                       &      "Criteria for phase",  i2, &
                       &      ": tester: ", l1,          &
                       &      ", dE: ", l1,              &
                       &      ", |dE|: ", l1,            &
                       &      ", Ncyc: ", l1 ,           &
                       &      ", Totcyc: ", l1 )' )      &
                       nextphase( iphase, i ), o1( i ), o2( i ), o3( i ), o4( i ), o5( i )
               End If
            End Do
         Else
            Write( iwr, '( 1x, i4, 2x, f14.8, 2x, g12.6, 2x, f12.6 )' ) &
                 niter, ehf + enuclear, energy_change, tester
         End If
      End If

    End Subroutine print_intermediate

    Subroutine converg_update

      Logical :: this_iter_best
      Integer :: linfo

      this_iter_best = iphase /= -1 .And. ehf < elow

      If( this_iter_best ) Then
         Call vector_save_orbitals( tmp_vectors_distrib, density_distrib )
         Call vector_copy  ( tmp_vectors_distrib, best_vectors_distrib )
         elow = ehf
      End If

      If( iphase /= new_phase ) Then

         If( opg_root() ) Then
            Write( iwr, * )
            If( new_phase == -1 ) Then
               Write( iwr, "( '**** Calculation aborted ' )" )
            Else If( new_phase /= 0 ) Then
!               Write( iwr, "( '**** Switch to convergence phase ', i0, /, / )" ) new_phase
               Write( iwr, "( '**** Switch to convergence phase ', i0 )" ) new_phase
               call strtrm(pchange_info,linfo)
               If ( linfo .gt. 0 ) then
                  Write( iwr, "( '**** ', a72)") pchange_info 
               End If
               Write( iwr, "( / )" )
            End If
         End If
         
         niter_phase = 0
         iphase      = new_phase
         converged   = iphase == 0
         
         ! Restore best vecs if required
         If( iphase > 0 ) Then
            If( restore_vec( iphase ) .And. .Not. this_iter_best ) Then
               Call vector_copy ( best_vectors_distrib, vectors_distrib )
               Call start_time_period( TP_F90_TDOWN )
               Call vector_tdown( vectors_distrib, vectors_no_symm_distrib, ctran, newbas0, otran )
               Call end_time_period( TP_F90_TDOWN )
               Call start_time_period( TP_F90_MAKE_DENS )
               Call vector_make_density( density_distrib, vectors_no_symm_distrib )
               Call end_time_period( TP_F90_MAKE_DENS )
            End If
         End If
         
         ! Reset DIIS if required
         If( iphase > 0 ) Then
            If( new_diis( iphase ) ) Then
               If( opg_root() ) Then
                  Write( iwr, "( 'Resetting DIIS' )" )
               End If
               Call vector_diis_reset
            End If
         End If
         
         ! Reset extrapolation
         if(iphase .gt. 0 .and. new_extrap(iphase))then
            write(iwr,*)'Resetting Extrap'
            num_extrap_cyc = 0
         endif

      End If
      
     End Subroutine converg_update

     Subroutine read_initial( uhf, vectors, overlap, hcore, overlap_noadapt )

       Logical       , Intent( In )   :: uhf
       Type( vector ), Dimension( : ) :: vectors
       Type( matrix ), Dimension( : ) :: overlap
       Type( matrix ), Dimension( : ) :: hcore
       Type( matrix ), Dimension( : ) :: overlap_noadapt

       Integer :: block

       ! Load guess vectors
       Call vector_read( vectors, (/ ibl3qa, ibl3qb /), (/ ibl3ea, ibl3eb /), idaf ) 
       
       !HACKY !!!!!!!!!!!!
       Call matrix_read_triangle( overlap( 1 ), ibl7st, num8 )
       If( uhf ) Then
          Call matrix_read_triangle( overlap( 2 ), ibl7st, num8 )
       End If
       
       ! Core Hamiltonian
       ! HACKY !!!!!!!!!!!!
       Call matrix_read_triangle( hcore( 1 ), ibl7f, num8 )
       If( uhf ) Then
          Call matrix_read_triangle( hcore( 2 ), ibl7f, num8 )
       End If

       Call secget( ionsec, 2, block )
       block = block + 1
       Call matrix_read_triangle( overlap_noadapt( 1 ), block, idaf )
       If( uhf ) Then
          Call matrix_read_triangle( overlap_noadapt( 2 ), block, idaf )
       End If
 
     End Subroutine read_initial

  End Subroutine newscf
    
  Subroutine get_n_jobs( uhf, n_jobs, n_all_jobs )

    Logical, Intent( In    ) :: uhf
    Integer, Intent(   Out ) :: n_jobs

    If( uhf ) Then
       n_jobs = 2
       n_all_jobs = 2
    Else
       n_jobs = 1
       n_all_jobs = 1
    End If

  End Subroutine get_n_jobs

  Subroutine fock_build( density, direct,    &
       fock, hcore, uhf, do_dft, &
       ehf, ehf1,ehf2, edft, &
       nquad, alpha_quad, beta_quad, &
       rdmat_work, prefac_work, focka_work, fockb_work, densa_work, densb_work, &
       core )
    !************************************************************************
    !
    !  construct fock matrix and components of the hf/ks energy
    !
    !   alphadensity  matrix tag           in
    !   betadensity   matrix tag           in
    !   direct        logical              in    direct-scf
    !   alphafock     matrix tag           out
    !   betafock      matrix tag           out
    !   ehf1          double precision     out
    !   ehf2          double precision     out
    !   edft          double precision     out
    !   nquad         integer              out dft results for diagostics 
    !   alpha_quad    real                 "
    !   beta_quad     real                 "
    !
    !************************************************************************
    !

    Use newscf_numbers
    Use allocation
    Use distributed_matrices
    Use distributed_vectors
    Use newscf_modules

    Implicit None

    Type( matrix ), Dimension( : )              :: density
    Logical,                    Intent( In    ) :: direct
    Type( matrix ), Dimension( : )              :: fock
    Type( matrix ), Dimension( : )              :: hcore
    Real( wp ),                 Intent(   Out ) :: ehf
    Real( wp ),                 Intent(   Out ) :: ehf1
    Real( wp ),                 Intent(   Out ) :: ehf2
    Real( wp ),                 Intent(   Out ) :: edft
    Integer   ,                 Intent(   Out ) :: nquad
    Real( wp ),                 Intent(   Out ) :: alpha_quad
    Real( wp ),                 Intent(   Out ) :: beta_quad
    Logical                   , Intent( In    ) :: uhf
    Logical                   , Intent( In    ) :: do_dft
    Real( wp ), Dimension( : ), Intent( In    ) :: rdmat_work
    Real( wp ), Dimension( : ), Intent( In    ) :: prefac_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: focka_work
    Real( wp ), Dimension( : ), Intent(   Out ) :: fockb_work
    Real( wp ), Dimension( : ), Intent( InOut ) :: densa_work
    Real( wp ), Dimension( : ), Intent( InOut ) :: densb_work
    Real( wp ), Dimension( * ), Intent(   Out ) :: core

    !
    ! control parameters direct scf
    !
    ! needed for correct functioning of integrals
    ! (o = output only)
    !
    ! n dlnmxd
    ! y dlntol    (dlntol is usually set to one of the tolitr values)
    ! n delfac   - delta factor for tolerance (delfac) when dlnmxd is small than 
    ! n tolitr(3) 
    !   itrtol   - iteraction count to control dlntol
    ! 
    !
    ! y oimag    - analyse integral magnitudes
    ! o intmag   - holds integral magnitudes
    ! o intcut   - integral test counts (numbers of integrals deleted for
    !              various reasons)
    ! y odscf    - run direct
    !
    !
    ! y irdmat
    ! y iprefa
    ! y iexch    - address of exchange integrals  ? workspace
    ! y ifock    - address for fock build 
    !
    !   nhamd   ???   only for gvb
    !   nsheld
    !
    ! -- not used in this version
    !
    !  iof171  - section 171 offsets
    !  m171t .. just 171
    !  odnew

    !
    !      real*8 dlnmxd, dlntol, tolitr, deltol, delfac
    !      integer ifock, idmat, irdmat, iprefa, intcut, intmag
    !      integer ibl171, lentri, itrtol, iexch, iof171, m171t
    !      integer nhamd, nsheld
    !      logical odscf, oimag, odnew, opref

    ! external functions

    Logical, External :: opg_root

    ! local variables

    Real( wp ), Dimension( : ), Allocatable :: scratch
    Real( wp ), Dimension( : ), Allocatable :: e_work

    Real( wp ) :: tracep, dum
    Real( wp ) :: etmp, dft_accu
    Real( wp ) :: facex

    Integer :: l1, l2, l3, l22
    Integer :: idum, i
    Integer :: idens,  iscr
    Integer :: idensb, ifockb
    Integer :: error

    Logical :: o2e, odft, ocoul, oexch
    Logical :: delta
    Logical :: out, outon

    l1 = num
    l2 = l1 *( l1 + 1 ) / 2
    l3 = l1 * l1

    l22 = 2 * l2

    out = nprint == 5
    outon = nprint /= -5

    edft = 0.0_wp

    !
    ! modify fock builder options according to
    ! required dft scheme (hf exchange, fitted coulomb etc
    !
    If( do_dft ) Then
       idum = cd_set_2e()
    End If

    Call matrix_get_to_triangle( density( 1 ), densa_work )
    If(uhf)Then
       Call matrix_get_to_triangle( density( 2 ), densb_work )
    End If

    focka_work = 0.0_wp

    !
    If( cd_2e() .And. do_dft ) Then
       odft = .True.
       If( uhf ) Then
          idum = cd_uks()
       End If
       !
       !     for simplicity, only implement full coulomb version
       !     here for the moment
       !
       If( cd_hf_exchange() ) Then
          facex = cd_hf_exchange_weight()
          oexch = .True.
       Else
          oexch = .False.
       End If
    Else
       odft = .False.
       facex = 1.0_wp
       oexch = .True.
    End If
    !
    ! coulomb
    !
    ocoul = cd_hf_coulomb() .Or. .Not. odft
    o2e   = .Not.( odft .And. .Not. oexch .And. .Not. ocoul )
    !
    
    Call start_time_period( TP_F90_INTS )

!!$    Call mp_enableintr( error )
    If( uhf ) Then
       fockb_work = 0.0_wp
       !
       If( o2e ) Then
          If( direct ) Then
             irest = 0
             Call dhstaru( core, densa_work( 1 ), densb_work( 1 ), &
                                 focka_work( 1 ), fockb_work( 1 ), &
                                 prefac_work, rdmat_work, irest )
          Else
             Call hstaru( densa_work, densb_work, &
                          focka_work, fockb_work, &
                          nopk )
          End If
       End If
    Else
       !
       ! rhf case
       If( o2e )Then
          If( direct )Then
             irest = 0
             Call dhstar( core, focka_work( 1 ), densa_work( 1 ), &
                                prefac_work, rdmat_work, irest )
          Else
             Allocate( scratch( 1:l2 ), Stat = error )
             If( error /= 0 ) Then
                Call dlc_error( 'Internal error: failed to alloc memory in FOCK_BUILD', 'abort' )
             End If
             Call alloc_add( ALLOC_REAL, Size( scratch ) )
             scratch = 0.0_wp
             ! problem if on cray, convex or titan - m4 removed
             Call hstar( densa_work, focka_work, scratch, nopk )
             Call alloc_free( ALLOC_REAL, Size( scratch ) )
             Deallocate( scratch )
          End If
       End If
    End If

!!$    Call mp_disableintr( error )
    Call end_time_period( TP_F90_INTS )

    !
    ! restore fock builder options
    !
    If( do_dft ) Then
       idum = cd_reset_2e()
    End If

    ! - skip zora for now

    !     ----- symmetrize skeleton fock matrix -----
    !
    !     -h- at x(i10)
    !     scratch area at x(i20)
    !
    ! IJB to save memory use densa_work as scratch and then
    ! reget ouyt of distributed array. Bit hacky really.
    Call symh( focka_work, densa_work, iky, 0, 0 )
    If( uhf )Call symh( fockb_work, densa_work, iky, 0, 0 )
    Call matrix_get_to_triangle( density( 1 ), densa_work )

    If( out ) Then
       Write( iwr, 9088 )
       If( uhf ) Then
          Write( iwr, 9108 )
       End If
       Call prtril( focka_work, l1 )
       If( uhf ) Then
          Write( iwr, 9128 )
          Call prtril( fockb_work, l1 )
       End If
    End If

9088 Format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23( &
         '-'))
9108 Format(//1x,129('-')/// &
         50x,'----- alpha set -----')
9128 Format(//1x,129('-')/// &
         50x,'----- beta set -----')

    Call matrix_set_from_triangle( fock( 1 ), focka_work )
    If( uhf ) &
         Call matrix_set_from_triangle( fock( 2 ), fockb_work )

    !
    !     ----- read in core hamiltonian matrix
    !           and calculate hf energy -----
    !
    !     -h0- at x(i20)
    !     - h- at x(i10)
    !     - e- at x(i30)
    !

    !      call vadd(q(i10),1,q(jblkf),1,q(i10),1,nx)

    !HACKY !!!!!!
    Call matrix_daxpy( 1.0_wp, hcore, fock )
    !
    !  ehf1 = trp(p*h_1)
    !
    Allocate( e_work( 1:Size( hcore ) ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in FOCK_BUILD', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( e_work ) )
    Call matrix_dot_product( hcore, density, e_work )
    ehf1 = Sum( e_work )
    Call alloc_free( ALLOC_DERIVED, Size( e_work ) )
    Deallocate( e_work )

    If( cd_active() .And. do_dft )Then
       !
       ! ccpdft: evaluate kohn sham energy expression
       !
       If( cd_hf_exchange() .Or. cd_hf_coulomb() )Then
          !
          ! coulomb or exchange operator is in q(ifock) augmented by h_1, 
          ! compute energy using hf expression (maybe without exchange)
          !
          Allocate( e_work( 1:Size( hcore ) ), Stat = error )
          If( error /= 0 ) Then
             Call dlc_error( 'Internal error: failed to alloc memory in FOCK_BUILD', 'abort' )
          End If
          Call alloc_add( ALLOC_DERIVED, Size( e_work ) )
          Call matrix_dot_product( density, fock, e_work )
          ehf2 = Sum( e_work )
          Call alloc_free( ALLOC_DERIVED, Size( e_work ) )
          Deallocate( e_work )
          etmp = ( ehf1 + ehf2 ) * 0.5_wp
          !
       Else
          !
          ! coulomb operator has not yet been evaluated
          ! (j energy will be part of edft)
          !
          etmp = ehf1
       End If

       !
       ! update kohn-sham matrix and compute fitted/integrated 
       ! energy terms
       !
       !c         if (diff.ne.0.0_wp) dft_accu = min(dft_accu,abs(diff))

       Call matrix_get_to_triangle( fock( 1 ), focka_work )
       If(uhf) &
            Call matrix_get_to_triangle( fock( 2 ), fockb_work )

       dft_accu = 0.0_wp

       idum = cd_energy_ao( c, focka_work( 1 ), fockb_work( 1 ), &
            densa_work( 1 ), densb_work( 1 ),                    &
            edft,core,core,outon,dft_accu,iwr )


       Call matrix_set_from_triangle( fock( 1 ), focka_work )
       If(uhf) &
            Call matrix_set_from_triangle( fock( 2 ), fockb_work )

       ehf = etmp+edft
       Call cd_get_dftresults( nquad, alpha_quad, beta_quad )

    Else
       !
       !  hartree-fock energy expression
       !  e = 1/2 [ trp(p*h_1) + trp(p*(h_1 + h_2)) ]
       !                ehf1            ehf2
       !
       Allocate( e_work( 1:Size( hcore ) ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in FOCK_BUILD', 'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( e_work ) )
       Call matrix_dot_product( density, fock, e_work )
       ehf2 = Sum( e_work )
       Call alloc_free( ALLOC_DERIVED, Size( e_work ) )
       Deallocate( e_work )
       ehf = ( ehf1 + ehf2 ) * 0.5_wp
    End If

    ! - nb skip code for crystal field
  End Subroutine fock_build

  Subroutine make_rdmat( rdmat, uhf, density, densa_work, densb_work )

    Use newscf_numbers
    Use allocation
    Use distributed_matrices

    Real( wp ), Dimension( : ), Intent(   Out ) :: rdmat
    Logical                   , Intent( In    ) :: uhf
    Type( matrix ), Dimension( : )              :: density
    Real( wp ), Dimension( : ), Intent( InOut ) :: densa_work
    Real( wp ), Dimension( : ), Intent( InOut ) :: densb_work

    Integer :: n, l2
    Integer :: error

    Call matrix_inquire( density( 1 ), global_n = n )
    l2 = ( n * ( n + 1 ) ) / 2

    If( uhf ) Then
       Call matrix_get_to_triangle( density( 1 ), densa_work )
       Call matrix_get_to_triangle( density( 2 ), densb_work )
       Call mkrdmt( 'uhf', rdmat, densa_work( 1 ), l2, 0 )
    Else
       Call matrix_get_to_triangle( density( 1 ), densa_work )
       Call mkrdmt( 'rhf', rdmat, densa_work( 1 ), l2, 0 )
    End If

  End Subroutine make_rdmat


  Subroutine initial_evecs( uhf, uhfatoms, direct, rdmat, prefac, focka, fockb, densa, densb, core )

    Use newscf_numbers
    Use allocation
    Use distributed_matrices
    Use distributed_vectors
    Use newscf_modules

    Implicit None

    Logical                   , Intent( In    ) :: uhf
    Logical                   , Intent( In    ) :: uhfatoms
    Logical,                    Intent( In    ) :: direct
    Real( wp ), Dimension( : ), Intent(   Out ) :: rdmat
    Real( wp ), Dimension( : ), Intent(   Out ) :: prefac
    Real( wp ), Dimension( : ), Intent(   Out ) :: focka
    Real( wp ), Dimension( : ), Intent(   Out ) :: fockb
    Real( wp ), Dimension( : ), Intent(   Out ) :: densa
    Real( wp ), Dimension( : ), Intent(   Out ) :: densb
    Real( wp ), Dimension( * ), Intent( InOut ) :: core

    Integer, External :: lensec
    Logical, External :: opg_root

    Type( matrix ), Dimension( : ), Allocatable :: density
    Type( matrix ), Dimension( : ), Allocatable :: fock_ao
    Type( matrix ), Dimension( : ), Allocatable :: fock_mo
    Type( matrix ), Dimension( : ), Allocatable :: hcore
    Type( matrix ), Dimension( : ), Allocatable :: overlap
    Type( vector ), Dimension( : ), Allocatable :: vectors_no_symm
    Type( vector ), Dimension( : ), Allocatable :: vectors_mo
    Type( vector ), Dimension( : ), Allocatable :: vectors_copy
    Type( vector ), Dimension( : ), Allocatable :: vectors

    Real( wp ), Parameter :: tol1 = 1.0e-5_wp

    Real( wp ) :: alpha_quad, beta_quad
    Real( wp ) :: ehf1, ehf2, edft, ehf
    Real( wp ) :: tol_depen

    Integer, Dimension( : ), Allocatable :: occs

    Integer :: n_spin
    Integer :: n_cart, n_harm
    Integer :: n_use
    Integer :: nquad
    Integer :: error
    Integer :: len1, len2, lenv, lenq
    Integer :: iblkv, iblk
    Integer :: i

    Character( Len = 8 ) :: scf_type

    Integer :: jjfile, lds, isect, ldsect, iacsct
    common/restri/jjfile(63),lds(508),isect(508),ldsect(508),iacsct(508)

    Call start_time_period( TP_DENSCF )

    n_cart = num
    n_harm = newbas0

    ! If the original guess was open shell need to get
    ! how many spin states
    If( uhfatoms ) Then
       n_spin = 2
       Allocate( occs( 1:n_spin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       occs = (/ na, nb /)
    Else
       n_spin = 1
       Allocate( occs( 1:n_spin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
       End If
       occs = (/ na /)
    End If
    Call alloc_add( ALLOC_INTEGER, Size( occs ) )

    Call matrix_comms_init( n_cart, mpi_comm_workers, blacs_context, n_spin, 1 )

    ! First of all set up I/O
    ! qmat section
    len1 = lensec( n_cart + 1 )
    lenv = len1 + lensec( n_cart * ( n_cart + 1 ) / 2 )
    lenq = lenv + lensec( n_cart * n_cart )
    Call secput( isecqm, 167, lenq, iblkv )
    ibl3qs = iblkv + lenv
    ! vectors section
    len1 = lensec( mach( 8 ) )
    lenv = lensec( n_cart * n_cart )
    len2 = lensec( mach( 9 ) )
    Call secput( mouta, 3, len1 + lenv + len2 + 1, iblk )
    ibl3qa = iblk + len2 + len1 + 1
    If( uhf ) Then
       Call secput( moutb, 3, len1 + lenv + len2 + 1, iblk )
       ibl3qb = iblk + len2 + len1 + 1
    End If

    If( direct ) Then
       scf_type = zscftp
       If( uhfatoms ) Then
          zscftp   = 'uhf'
       Else
          zscftp   = 'rhf'
       End If
    End If

    Allocate( fock_ao( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( fock_ao ) )
    Call matrix_create( n_cart, n_cart, fock_ao, (/ ( 'fock_ao', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( density( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( density ) )
    Call matrix_create( n_cart, n_cart, density, (/ ( 'density', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( hcore( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( hcore ) )
    Call matrix_create( n_cart, n_cart, hcore, (/ ( 'hcore', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Call matrix_read_triangle( hcore( 1 ), ibl7f, num8 )
    Call matrix_read_triangle( density( 1 ), ibl3pa, idaf )
    If( n_spin == 2 ) Then
       Call matrix_read_triangle( hcore( 2 ), ibl7f, num8 )
       Call matrix_read_triangle( density( 2 ), ibl3pb, idaf )
    End If

    ! Generate Fock matrix
    Call start_time_period( TP_DENSCF_BUILD )
    If( direct )Then
       Call start_time_period( TP_F90_RDMAT )
       Call make_rdmat( rdmat, uhfatoms, density, densa, densb )
       Call end_time_period( TP_F90_RDMAT )
       dlntol = tolitr( 1 ) - Min( Max( dlnmxd, delfac ), 0.0_wp )
       If( opg_root() ) Then
          Write( iwr, "( 5x, 'dlntol =  ', f15.10 / 5x, 25( '*' ) )" ) dlntol
       End If
    End If
    Call fock_build( density, direct, &
         fock_ao, hcore, &
         uhfatoms, .False., ehf, ehf1, ehf2, edft, &
         nquad, alpha_quad, beta_quad, &
         rdmat, prefac, focka, fockb, densa, densb, &
         core )
    Call end_time_period( TP_DENSCF_BUILD )

    ! OK done with hcore
    Call matrix_destroy( hcore )
    Call alloc_free( ALLOC_DERIVED, Size( hcore ) )
    Deallocate( hcore )

    If( direct ) Then
       ! Need to save prefactors and reduced density matrix to disk
       Call secput( isect( 471 ), 171, 2 * lensec( ikyp( nshell ) ), ibl171 )
       Call wrt3 ( prefac, ikyp( nshell ), ibl171, idaf )
       Call wrt3s( rdmat , ikyp( nshell ), idaf )
       Call clredx
    End If

    ! Now need to generate orthonormalizing transformation

    ! First read in the overlap matrix
    Allocate( overlap( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( overlap ) )
    Call matrix_create( n_cart, n_cart, overlap, (/ ( 'overlap', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )
    Call matrix_read_triangle( overlap( 1 ), ibl7st, num8 )
    Call matrix_write_triangle( overlap( 1 ), iblkv + lensec( n_cart + 1 ), idaf )
    If( n_spin == 2 ) Then
       Call matrix_read_triangle( overlap( 2 ), ibl7st, num8 )
    End If

    If( i_depen_check < 0 ) Then
       tol_depen = 1.0e-7_wp
    Else If( i_depen_check == 0 ) Then
       tol_depen = 1.0e-20_wp
    Else
       tol_depen = 10.0 ** ( - i_depen_check )
    End If

    ! Now diag and scale with evecs. Note VECTORS is initialized in this routine
    Call start_time_period( TP_DENSCF_DIAG )
    Allocate( vectors( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors ) )
    Call vector_initialize( n_cart, n_harm, tol1, tol_depen, overlap, &
         n_use, vectors )
    Call end_time_period( TP_DENSCF_DIAG )

    If( n_use /= n_harm ) Then
       newbash = newbas0
       newbas0 = n_use
       n_harm  = n_use
       odepen  = .True.
    End If

    ! O.K., have orthonorming transformation + a fock matrix. Can move toward
    ! doing diag
    Allocate( vectors_no_symm( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_no_symm ) )
    Call vector_create( vectors_no_symm, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_cart, n_harm, (/ ( 'Vectors no symm ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( vectors_mo( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_mo ) )
    Call vector_create( vectors_mo, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_harm, n_harm, (/ ( 'Vectors MO ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( vectors_copy( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( vectors_copy ) )
    Call vector_create( vectors_copy, (/ ( 0, i = 1, n_spin ) /), (/ ( 1, i = 1, n_spin ) /), &
         n_cart, n_harm, (/ ( 'Vectors copy ', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Allocate( fock_mo( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to alloc memory in NEWSCF', 'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( fock_mo ) )
    Call matrix_create( n_harm, n_harm, fock_mo, (/ ( 'fock mo', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_UHF_K )

    Call vector_tdown_setup( vectors_no_symm, ctran, ntran, itran, ilifc, otran )

    Call vector_copy( vectors, vectors_copy )

    ! And generate the unsymmetrized version
    Call start_time_period( TP_DENSCF_TDOWN )
    Call vector_tdown( vectors, vectors_no_symm, ctran, n_harm, otran )
    Call end_time_period( TP_DENSCF_TDOWN )

    ! Transform fock AO -> MO
    Call start_time_period( TP_DENSCF_SIMIL )
    Call vector_similarity( fock_ao, vectors_no_symm, fock_mo )
    Call end_time_period( TP_DENSCF_SIMIL )

    ! Diag
    Call start_time_period( TP_DENSCF_DIAG )
    Call vector_diagonalise( fock_mo, vectors_mo )
    Call end_time_period( TP_DENSCF_DIAG )

    ! Transform vectors MO -> AO
    Call start_time_period( TP_DENSCF_BACK )
    Call vector_back_transform( vectors_mo, vectors_copy, vectors )
    Call end_time_period( TP_DENSCF_BACK )

    ! Occupancies
    Call vector_assign_by_energy( vectors, occs )

    ! And generate a density
    Call start_time_period( TP_DENSCF_TDOWN )
    Call vector_tdown( vectors, vectors_no_symm, ctran, n_harm, otran )
    Call end_time_period( TP_DENSCF_TDOWN )
    Call start_time_period( TP_DENSCF_MAKE_DENS )
    Call vector_make_density( density, vectors_no_symm )
    Call end_time_period( TP_DENSCF_MAKE_DENS )

    ! And save
    Call vector_write( (/ vectors( 1 ) /), (/ ibl3qa /), (/ ibl3ea /), idaf )
    If( uhf .Or. n_spin == 2 ) Then
       Call vector_write( (/ vectors( n_spin ) /), (/ ibl3qb /), (/ ibl3eb /), idaf )
    End If

    Call matrix_tdown_free

    Call matrix_destroy( fock_mo )
    Call alloc_free( ALLOC_DERIVED, Size( fock_mo ) )
    Deallocate( fock_mo )

    Call vector_destroy( vectors_copy )
    Call alloc_free( ALLOC_DERIVED, Size( vectors_copy ) )
    Deallocate( vectors_copy )

    Call vector_destroy( vectors_mo )
    Call alloc_free( ALLOC_DERIVED, Size( vectors_mo ) )
    Deallocate( vectors_mo )

    Call vector_destroy( vectors_no_symm )
    Call alloc_free( ALLOC_DERIVED, Size( vectors_no_symm ) )
    Deallocate( vectors_no_symm )

    Call vector_destroy( vectors )
    Call alloc_free( ALLOC_DERIVED, Size( vectors ) )
    Deallocate( vectors )

    Call matrix_destroy( overlap )
    Call alloc_free( ALLOC_DERIVED, Size( overlap ) )
    Deallocate( overlap )

    Call matrix_destroy( density )
    Call alloc_free( ALLOC_DERIVED, Size( density ) )
    Deallocate( density )

    Call matrix_destroy( fock_ao )
    Call alloc_free( ALLOC_DERIVED, Size( fock_ao ) )
    Deallocate( fock_ao )

    If( direct ) Then
       zscftp = scf_type
    End If

    Call alloc_free( ALLOC_INTEGER, Size( occs ) )
    Deallocate( occs )

    Call matrix_comms_finalize

    Call end_time_period( TP_DENSCF )

  End Subroutine initial_evecs

End Module newscf_routines

Function nuclear_energy()

  Use newscf_numbers
  Use newscf_modules, Only : nat, czan, c

  Implicit None

  Real( wp ) :: nuclear_energy

  Real( wp ), External :: enucf

  nuclear_energy = enucf( nat, czan, c )

End Function nuclear_energy

