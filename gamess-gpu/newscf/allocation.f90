Module allocation

  Implicit None

  Integer, Parameter, Public :: ALLOC_COMPLEX = 1
  Integer, Parameter, Public :: ALLOC_REAL    = 2
  Integer, Parameter, Public :: ALLOC_INTEGER = 3
  Integer, Parameter, Public :: ALLOC_LOGICAL = 4
  Integer, Parameter, Public :: ALLOC_DERIVED = 5

  Logical, Parameter, Public :: ALLOC_VERBOSE    = .True.
  Logical, Parameter, Public :: ALLOC_NO_VERBOSE = .False.

  Public :: alloc_set_verbosity
  Public :: alloc_add
  Public :: alloc_free
  Public :: alloc_report
  Public :: alloc_reset

  Private

  Integer, Dimension( 1:5 ) :: current = 0
  Integer, Dimension( 1:5 ) :: top     = 0
  Integer, Dimension( 1:5 ) :: active  = 0
  Integer, Dimension( 1:5 ) :: max_active  = 0

  Logical :: verbose = .False.

  Character( Len = 7 ), Dimension( 1:5 ), Parameter :: names = &
       (/ 'complex', 'real   ', 'integer', 'logical', 'derived' /)

Contains

  Subroutine alloc_set_verbosity( value )

    Logical, Intent( In ) :: value

    verbose = value

  End Subroutine alloc_set_verbosity

  Subroutine alloc_add( what, eles )

    Integer, Intent( In ) :: what
    Integer, Intent( In ) :: eles

    Logical, External :: opg_root

    current   ( what ) = current( what ) + eles
    top       ( what ) = Max( top( what ), current( what ) )
    active    ( what ) = active( what ) + 1
    max_active( what ) = Max( max_active( what ), active( what ) )

    If( verbose ) Then
       If( opg_root() ) Then
          Write( *, '( "Allocating ", i0, 1x, a, "s" )' ) eles, &
               names( what )( 1:Len_trim( names( what ) ) )
          Call alloc_report
       End If
    End If

  End Subroutine alloc_add

  Subroutine alloc_free( what, eles )

    Integer, Intent( In ) :: what
    Integer, Intent( In ) :: eles

    Logical, External :: opg_root

    current( what ) = current( what ) - eles
    active ( what ) = active( what ) - 1

    If( verbose ) Then
       If( opg_root() ) Then
          Write( *, '( "Freeing ", i0, 1x, a, "s" )' ) eles, &
               names( what )( 1:Len_trim( names( what ) ) )
          Call alloc_report
       End If
    End If

  End Subroutine alloc_free

  Subroutine alloc_report

    Logical, External :: opg_root

    Integer :: i

    If( opg_root() ) Then
       Write( *, * )
       Write( *, '( "Current F90 allocation status" )' )
       Write( *, '( 3x, "Type", 2x, 2x, "Active", 3x, 2x, &
            & "Current", 2x, "Max Active", 1x, 3x, "Top" )' )
       Do i = 1, Size( current )
          Write( *, '( a7, ":", 1x, i10, 1x, i10, 1x, i10, 1x, i10 )' ) &
               names( i ), active( i ), current( i ), max_active( i ), top( i )
       End Do
       Write( *, * )
    End If

  End Subroutine alloc_report

  Subroutine alloc_reset

    current = 0
    top     = 0
    active  = 0

  End Subroutine alloc_reset

End Module allocation
