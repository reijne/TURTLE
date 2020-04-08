Module analysis

  Use newscf_numbers
  Use distributed_matrices

  Implicit None

  Public :: analysis_driver

  Private

  Type( matrix ), Save                        :: overlap
  Type( matrix ), Dimension( : ), Allocatable :: evecs
  Type( matrix ), Dimension( : ), Allocatable :: density
  Type( matrix ), Dimension( : ), Allocatable :: dipole

  Integer :: n_cart
  Integer :: n_spin

Contains

  Subroutine analysis_driver( scf_type )

    Use newscf_modules

    Character( Len = 8 ), Intent( In ) :: scf_type

    n_cart = num

    If( scf_type /= 'rhf' ) Then
       n_spin = 2
    Else
       n_spin = 1
    End If

    Call analysis_init

    Call analysis_dipole

    Call analysis_population

    Call analysis_finalize

  End Subroutine analysis_driver

  Subroutine analysis_init

    Use newscf_modules

    ! Initialise the matrix stuff
    Call matrix_comms_init( n_cart, mpi_comm_workers, blacs_context, 1, 1 )

    ! Allocate the matrices
    Call analysis_allocate

    ! Read in the matrices
    Call analysis_read

    ! Unadapt the evecs
    Call analysis_tdown

    ! Transform the d funcs ( if any )
    Call analysis_transd

  End Subroutine analysis_init

  Subroutine analysis_allocate

    Use allocation

    Integer :: i
    Integer :: error

    Call matrix_create( n_cart, n_cart, overlap, 'overlap', MATRIX_REAL, MATRIX_DISTRIB_ALL )

    Allocate( evecs( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_ALLOCATE', &
            'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( evecs ) )
    Call matrix_create( n_cart, n_cart, evecs, (/ ( 'evecs', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_ALL )

    Allocate( density( 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_ALLOCATE', &
            'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( density ) )
    Call matrix_create( n_cart, n_cart, density, (/ ( 'density', i = 1, n_spin ) /), &
         (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_ALL )

    Allocate( dipole( 1:3 ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_ALLOCATE', &
            'abort' )
    End If
    Call alloc_add( ALLOC_DERIVED, Size( dipole ) )
    Call matrix_create( n_cart, n_cart, dipole, (/ ( 'dipole', i = 1, 3 ) /), &
         (/ ( MATRIX_REAL, i = 1, 3 ) /), MATRIX_DISTRIB_ALL )

  End Subroutine analysis_allocate

  Subroutine analysis_finalize

    Use allocation

    Call matrix_destroy( dipole )
    Call alloc_free( ALLOC_DERIVED, Size( dipole ) )
    Deallocate( dipole )

    Call matrix_destroy( density )
    Call alloc_free( ALLOC_DERIVED, Size( density ) )
    Deallocate( density )

    Call matrix_destroy( evecs )
    Call alloc_free( ALLOC_DERIVED, Size( evecs ) )
    Deallocate( evecs )

    Call matrix_destroy( overlap )

    Call matrix_comms_finalize

  End Subroutine analysis_finalize

  Subroutine analysis_read

    Use newscf_modules

    Integer, External :: lensec

    Integer :: block
    Integer :: len
    Integer :: i

    Integer :: jjfile, lds, isect, ldsect, iacsct
    Common/restri/jjfile(63),lds(508),isect(508),ldsect(508),iacsct(508)

    len = n_cart * ( n_cart + 1 ) / 2

    ! Read in overlap
    Call secget( ionsec, 2, block )
    block = block + 1
    Call matrix_read_triangle( overlap, block, idaf )

    ! Read in dipole matrices
    block = block + 3 * lensec( len )
    Do i = 1, 3
       Call matrix_read_triangle( dipole( i ), block, idaf )
       block = block + lensec( len )
    End Do

    ! Read in density matrices
    Call secget( isect( 497 ), 19, block )
    Call matrix_read_triangle( density( 1 ), block, idaf )
    If( n_spin == 2 ) Then
       block = block + lensec( n_cart ) + lensec( len )
       Call matrix_read_triangle( density( 2 ), block, idaf )
    End If

    ! Read in evecs
    Call secget( mouta, 3, block )
    block = block + 1 + lensec( mach( 8 ) ) + lensec( mach( 9 ) )
    Call matrix_read_rectangle( evecs( 1 ), block, idaf )
    If( n_spin == 2 ) Then
       Call secget( moutb, 3, block )
       block = block + 1 + lensec( mach( 8 ) ) + lensec( mach( 9 ) )
       Call matrix_read_rectangle( evecs( 2 ), block, idaf )
    End If

  End Subroutine analysis_read

  Subroutine analysis_tdown

    Use newscf_modules
    Use allocation

    Integer :: error
    Integer :: i

    Type( matrix ), Dimension( : ), Allocatable :: evecs_no_symm

    If( .Not. otran ) Then

       Allocate( evecs_no_symm( 1:n_spin ), Stat = error )
       If( error /= 0 ) Then
          Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_ALLOCATE', &
               'abort' )
       End If
       Call alloc_add( ALLOC_DERIVED, Size( evecs_no_symm ) )
       Call matrix_create( n_cart, n_cart, evecs_no_symm, (/ ( 'evecs', i = 1, n_spin ) /), &
            (/ ( MATRIX_REAL, i = 1, n_spin ) /), MATRIX_DISTRIB_ALL )

       Call matrix_tdown_setup( evecs( 1 ), ctran, ntran, itran, ilifc, otran )
       Call matrix_tdown( evecs, evecs_no_symm, ctran, n_cart, otran )
       Call matrix_tdown_free

       Call matrix_destroy( evecs_no_symm )
       Call alloc_free( ALLOC_DERIVED, Size( evecs_no_symm ) )
       Deallocate( evecs_no_symm )

    End If

  End Subroutine analysis_tdown

  Subroutine analysis_transd

    Use newscf_modules

    Real( wp ), Dimension( 1:3, 1:3 ), Parameter :: dtr = Reshape( (/             &
         0.4472135954999579_wp,  0.4472135954999579_wp, 0.4472135954999579_wp,    &
         0.8660254037844386_wp, -0.8660254037844386_wp, 0.0000000000000000_wp,    &
        -0.5000000000000000_wp, -0.5000000000000000_wp, 1.0000000000000000_wp /), &
         (/ 3, 3 /) )
    Real( wp ), Dimension( 1:3, 1:3 ), Parameter :: dtr_inv = Reshape( (/         &
         0.7453559924999299_wp,  0.7453559924999299_wp, 0.7453559924999299_wp,    &
         0.5773502691896258_wp, -0.5773502691896258_wp, 0.0000000000000000_wp,    &
        -0.3333333333333333_wp, -0.3333333333333333_wp, 0.6666666666666667_wp /), &
         (/ 3, 3 /) )

    Type( matrix ) :: trans
    Type( matrix ) :: work1

    Integer :: block
    Integer :: i

    Logical :: d_in_basis

    d_in_basis = Any( ktype( 1:nshell ) == 3 )

    If( d_in_basis ) Then

       Call matrix_create( n_cart, n_cart, trans    , 'trans    ', MATRIX_REAL, MATRIX_DISTRIB_ALL )

       ! Initialize the `forward' transformation
       Call trans_init( nshell, ktype, kloc, dtr, trans )

       Call matrix_create( n_cart, n_cart, work1    , 'work1    ', MATRIX_REAL, MATRIX_DISTRIB_ALL )

       ! Transform the overlap matrix
       Call matrix_mult2( overlap, trans, work1 )
       Call matrix_copy( work1, overlap )

       ! Transform the dipole matrix elements
       Do i = 1, 3
          Call matrix_mult2( dipole( i ), trans, work1 )
          Call matrix_copy( work1, dipole( i ) )
       End Do

       ! And now for the `back' transformations
       Call trans_init( nshell, ktype, kloc, dtr_inv, trans )

       ! Transform the density matrices
       Do i = 1, n_spin
          Call matrix_mult2( density( i ), trans, work1 )
          Call matrix_copy( work1, density( i ) )
       End Do
          

       ! Transform the vectors
       Do i = 1, n_spin
          Call matrix_dgemm( 't', 'n', 1.0_wp, trans, evecs( i ), 0.0_wp, work1 )
          Call matrix_copy( work1, evecs( i ) )
       End Do

       Call matrix_destroy( work1     )

       Call matrix_destroy( trans     )

    End If

  Contains

    Subroutine trans_init( n_shell, shell_type, shell_loc, dtr, trans )
      
      Integer                          , Intent( In ) :: n_shell
      Integer   , Dimension( : )       , Intent( In ) :: shell_type
      Integer   , Dimension( : )       , Intent( In ) :: shell_loc
      Real( wp ), Dimension( 1:3, 1:3 ), Intent( In ) :: dtr
      Type( matrix )                                  :: trans

      Character( Len = 8 ), Dimension( 1:3 ), Parameter :: labels = (/ &
           '     rr ', '   xx-yy', '   zz-rr' /)

      Integer :: n
      Integer :: i

      Call matrix_inquire( trans, global_n = n )

      Call matrix_set( trans, 0.0_wp )

      Do i = 1, n
         Call matrix_set( trans, 1.0_wp, i, i )
      End Do

      Do i = 1, n_shell
         If( shell_type( i ) == 3 ) Then
            Call matrix_set_patch( trans, dtr, shell_loc( i ), shell_loc( i ) )
            zbflab( shell_loc( i ): shell_loc( i ) + 2 ) = labels
         End If
      End Do

    End Subroutine trans_init

  End Subroutine analysis_transd

  Subroutine analysis_population

    Use allocation
    Use newscf_modules

    Logical, External :: opg_root

    Real( wp ), Dimension( :, : ), Allocatable :: mulliken_orb_pops
    Real( wp ), Dimension( :, : ), Allocatable :: lowdin_orb_pops
    Real( wp ), Dimension( :, : ), Allocatable :: mulliken_atom_pops
    Real( wp ), Dimension( :, : ), Allocatable :: lowdin_atom_pops

    Integer :: error

    Allocate( mulliken_orb_pops( 1:n_cart, 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
            'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( mulliken_orb_pops ) )

    Allocate( lowdin_orb_pops( 1:n_cart, 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
            'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( lowdin_orb_pops ) )

    Allocate( mulliken_atom_pops( 1:nat, 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
            'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( mulliken_atom_pops ) )

    Allocate( lowdin_atom_pops( 1:nat, 1:n_spin ), Stat = error )
    If( error /= 0 ) Then
       Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
            'abort' )
    End If
    Call alloc_add( ALLOC_REAL, Size( lowdin_atom_pops ) )

    Call analysis_pops_mulliken

    Call analysis_pops_lowdin

    Call analysis_pops_orbs_to_atoms

    Call analysis_pops_print

    Call alloc_free( ALLOC_REAL, Size( lowdin_atom_pops  ) )
    Deallocate( lowdin_atom_pops )

    Call alloc_free( ALLOC_REAL, Size( mulliken_atom_pops  ) )
    Deallocate( mulliken_atom_pops )

    Call alloc_free( ALLOC_REAL, Size( lowdin_orb_pops  ) )
    Deallocate( lowdin_orb_pops )

    Call alloc_free( ALLOC_REAL, Size( mulliken_orb_pops  ) )
    Deallocate( mulliken_orb_pops )

  Contains

    Subroutine analysis_pops_mulliken

      Type( matrix ) :: work

      Integer :: i, j
      Integer :: error

      Call matrix_create( n_cart, n_cart, work   , 'work   ', MATRIX_REAL, MATRIX_DISTRIB_ALL )
      
      Do i = 1, n_spin
         Call matrix_multiply( 1.0_wp, density( i ), overlap, 0.0_wp, work )
         Do j = 1, n_cart
            Call matrix_get( work, mulliken_orb_pops( j, i ), j, j )
         End Do
      End Do
      
      Call matrix_destroy( work )
      
    End Subroutine analysis_pops_mulliken

    Subroutine analysis_pops_lowdin

      Type( matrix ) :: s_sqrt
      Type( matrix ) :: s_evecs
      Type( matrix ) :: work

      Real( wp ), Dimension( : ), Allocatable :: s_evals

      Integer :: i, j
      Integer :: error

      ! Generate S**(1/2)
      Call matrix_create( n_cart, n_cart, work   , 'work   ', MATRIX_REAL, MATRIX_DISTRIB_ALL )
      Call matrix_create( n_cart, n_cart, s_sqrt , 's_sqrt ', MATRIX_REAL, MATRIX_DISTRIB_ALL )
      Call matrix_create( n_cart, n_cart, s_evecs, 's_evecs', MATRIX_REAL, MATRIX_DISTRIB_ALL )

      Allocate( s_evals( 1:n_cart ), Stat = error )
      If( error /= 0 ) Then
         Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
              'abort' )
      End If
      Call alloc_add( ALLOC_REAL, Size( s_evals ) )

      ! Need to save the overlap matrix as the diag destroys it.
      Call matrix_copy( overlap, work )
      Call matrix_diagonalise( overlap, s_evecs, s_evals )
      Call matrix_copy( work, overlap )

      s_evals = Sqrt( s_evals )

      Call matrix_mult_by_diag( s_evecs, s_evals, work )

      Call matrix_dgemm( 'n', 't', 1.0_wp, work, s_evecs, 0.0_wp, s_sqrt )

      ! Now generate lowdin pops - use s_evecs as workspace
      Do i = 1, n_spin
         Call matrix_multiply( 1.0_wp, density( i ), s_sqrt, 0.0_wp, s_evecs )
         Call matrix_multiply( 1.0_wp, s_sqrt, s_evecs, 0.0_wp, work )
         Do j = 1, n_cart
            Call matrix_get( work, lowdin_orb_pops( j, i ), j, j )
         End Do
      End Do

      Call alloc_free( ALLOC_REAL, Size( s_evals ) )
      Deallocate( s_evals )

      Call matrix_destroy( s_evecs )
      Call matrix_destroy( s_sqrt  )
      Call matrix_destroy( work )
      
    End Subroutine analysis_pops_lowdin

    Subroutine analysis_pops_orbs_to_atoms

      Integer, Dimension( : ), Allocatable :: orb_to_atom

      Integer :: this_atom
      Integer :: error
      Integer :: i, j

      Allocate( orb_to_atom( 1:n_cart ), Stat = error )
      If( error /= 0 ) Then
         Call dlc_error( 'Internal error: failed to allocate memory in ANALYSIS_POPULATION', &
              'abort' )
      End If
      Call alloc_add( ALLOC_INTEGER, Size( orb_to_atom ) )

      Do i = 1, nshell
         this_atom = katom( i )
         Do j = kloc( i ), kloc( i ) + kmax( i ) - kmin( i )
            orb_to_atom( j ) = this_atom
         End Do
      End Do

      this_atom = orb_to_atom( 1 )

      mulliken_atom_pops( this_atom, : ) = mulliken_orb_pops( 1, : )
      lowdin_atom_pops  ( this_atom, : ) = lowdin_orb_pops  ( 1, : )

      Do i = 2, n_cart

         If( orb_to_atom( i ) == this_atom ) Then

            mulliken_atom_pops( this_atom, : ) = mulliken_atom_pops( this_atom, : ) + &
                 mulliken_orb_pops( i, : )
            lowdin_atom_pops( this_atom, : ) = lowdin_atom_pops( this_atom, : ) + &
                 lowdin_orb_pops( i, : )

         Else

            this_atom = orb_to_atom( i )
            mulliken_atom_pops( this_atom, : ) = mulliken_orb_pops( i, : )
            lowdin_atom_pops  ( this_atom, : ) = lowdin_orb_pops  ( i, : )

         End If

      End Do

      Call alloc_free( ALLOC_INTEGER, Size( orb_to_atom ) )
      Deallocate( orb_to_atom )

    End Subroutine analysis_pops_orbs_to_atoms

    Subroutine analysis_pops_print

      Character( Len = 8 ), dimension( 1:2 ), Parameter :: spin_types = &
           (/ '   alpha', '    beta' /)

      Integer :: i, j

      Character( Len = 8 ) :: print_type

      If( n_spin == 2 ) Then

         ! Alpha and beta terms
         Do i = 1, n_spin

            print_type = spin_types( i )

            Write( iwr, "( // 1x, 80( '=' ) // 10x, 39( '-' ) / 10x,                &
                 &        'mulliken and lowdin population ', 'analyses', 10x, a8,   &
                 &        ' electrons' / 10x, 39 ( '-' ) ) " ) print_type

            Write( iwr, "( / 10x,'----- total gross population in aos ------' / )" )
            Do j = 1, n_cart
               Write( iwr, "( 10x, i5, 2x, a10, 2f12.5 )" ) j, zbflab( j ), &
                    mulliken_orb_pops( j, i ), lowdin_orb_pops( j, i )
            End Do
            
            Write( iwr, "( / 10x,'----- condensed to atoms -----' / )" )
            
            Write( iwr, "( / 10x,'----- total gross population on atoms ----' / )" )
            Do j = 1, nat
               Write( iwr, "( 10x, i5, 2x, a8, 2x, f6.1, 2f12.5 )" ) j, zaname( j ), czan( j ), &
                    mulliken_atom_pops( j, i ), lowdin_atom_pops( j, i )
            End Do
            
            ! Punch files
            Call blkao( mulliken_orb_pops ( :, i ), lowdin_orb_pops ( :, i ), print_type )
            Call blkat( mulliken_atom_pops( :, i ), lowdin_atom_pops( :, i ), print_type )
            
         End Do

         ! Spin terms

         Write( iwr, "( // 1x, 80( '=' ) // 10x, 39( '-' ) / 10x,                &
              &        'mulliken and lowdin population ', 'analyses', 10x, a8,   &
              &        ' electrons' / 10x, 39 ( '-' ) ) " ) '    spin'

         Write( iwr, "( / 10x,'----- aos spin population ------' / )" )
         Do j = 1, n_cart
            Write( iwr, "( 10x, i5, 2x, a10, 2f12.5 )" ) j, zbflab( j ), &
                 mulliken_orb_pops( j, 1 ) - mulliken_orb_pops( j, 2 ), &
                 lowdin_orb_pops  ( j, 1 ) - lowdin_orb_pops  ( j, 2 )
         End Do

         Write( iwr, "( / 10x,'----- condensed to atoms -----' / )" )

         Write( iwr, "( / 10x,'----- atomic spin population -----' / )" )
         Do j = 1, nat
            Write( iwr, "( 10x, i5, 2x, a8, 2x, f6.1, 2f12.5 )" ) j, zaname( j ), czan( j ), &
                 mulliken_atom_pops( j, 1 ) - mulliken_atom_pops( j, 2 ), &
                 lowdin_atom_pops  ( j, 1 ) - lowdin_atom_pops  ( j, 2 )
         End Do

         ! Punch files
         Call blkao( mulliken_orb_pops( :, 1 ) - mulliken_orb_pops( :, 2 ), &
                     lowdin_orb_pops  ( :, 1 ) - lowdin_orb_pops  ( :, 2 ), &
                     '    spin' )
         Call blkat( mulliken_atom_pops( :, 1 ) - mulliken_atom_pops( :, 2 ), &
                     lowdin_atom_pops  ( :, 1 ) - lowdin_atom_pops  ( :, 2 ), &
                     '    spin' )

      End If

      ! Charge terms - only thing required for RHF
      Write( iwr, "( // 1x, 80( '=' ) // 10x, 39( '-' ) / 10x,          &
           &        'mulliken and lowdin population ', 'analyses', 10x, a8,   &
           &        ' electrons' / 10x, 39 ( '-' ) ) " ) '     all'

      Write( iwr, "( / 10x,'----- total gross population in aos ------' / )" )
      Do j = 1, n_cart
         Write( iwr, "( 10x, i5, 2x, a10, 2f12.5 )" ) j, zbflab( j ), &
              Sum( mulliken_orb_pops( j, : ) ), Sum( lowdin_orb_pops( j, : ) )
      End Do

      Write( iwr, "( / 10x,'----- condensed to atoms -----' / )" )
      
      Write( iwr, "( / 10x,'----- total gross population on atoms ----' / )" )
      Do j = 1, nat
         Write( iwr, "( 10x, i5, 2x, a8, 2x, f6.1, 2f12.5 )" ) j, zaname( j ), czan( j ), &
              Sum( mulliken_atom_pops( j, : ) ), Sum( lowdin_atom_pops( j, : ) )
      End Do
      
      ! Punch files
      Call blkao( Sum( mulliken_orb_pops , Dim = 2 ), Sum( lowdin_orb_pops , dim = 2 ), '     all' )
      Call blkat( Sum( mulliken_atom_pops, Dim = 2 ), Sum( lowdin_atom_pops, Dim = 2 ), '     all' )

    End Subroutine analysis_pops_print

  End Subroutine analysis_population

  Subroutine analysis_dipole

    Use newscf_modules

    Real( wp ), Parameter :: dipole_tol = 1.0e-16_wp
    Real( wp ), Parameter :: au_to_debye = 2.54158059_wp

    Real( wp ), Dimension( 1:3 ) :: dipole_ele
    Real( wp ), Dimension( 1:3 ) :: dipole_nuc
    Real( wp ), Dimension( 1:3 ) :: dipole_tot

    Real( wp ) :: dipole_au
    Real( wp ) :: dipole_debye
    Real( wp ) :: work

    Integer :: i, j

    Character, Dimension( 1:3 ), Parameter :: components = (/ 'x', 'y', 'z' /)

    ! Electronic term
    dipole_ele = 0.0_wp
    Do i = 1, n_spin
       Do j = 1, 3
          Call matrix_dot_product( density( i ), dipole( j ), work )
          dipole_ele( j ) = dipole_ele( j ) + work
       End Do
    End Do

    ! Nuclear term
    dipole_nuc = 0.0_wp
    Do i = 1, 3
       dipole_nuc( i ) = Sum( czan( 1:nat ) * c( i, 1:nat ), &
            Mask = zaname( 1:nat )( 1:2 ) /= 'bq' .Or. .Not. ocryst )
    End Do

    dipole_tot = dipole_ele + dipole_nuc

    dipole_au = Dot_product( dipole_tot, dipole_tot )
    If( dipole_au > dipole_tol ) Then
       dipole_au = Sqrt( dipole_au )
    End If

    dipole_debye = dipole_au * au_to_debye

    If( ocryst ) Then
       Write( iwr, '( 17x, "dipole moments excluding bq centres" )' )
    Else
       Write( iwr, '( 17x, "dipole moments" )' )
    End If

    Write( iwr, '( //, 11x, "nuclear", 6x, "electronic", 11x, "total", / )' )
    Do i = 1, 3
       Write( iwr, '( 1x, a1, 3f16.7 )' ) components( i ), &
            dipole_nuc( i ), dipole_ele( i ), dipole_tot( i )
    End Do

    Write( iwr, '( /, " total dipole moment = ", f16.7, " (a.u.)", /, &
         & 23x, f16.7, " (debye)" )' ) dipole_au, dipole_debye

    Call blkdip( dipole_tot )

  End Subroutine analysis_dipole

End Module analysis

Subroutine analysis_wrap( scf_type )

  Use analysis

  Character( Len = 8 ), Intent( In ) :: scf_type

  Call analysis_driver( scf_type )

End Subroutine analysis_wrap
