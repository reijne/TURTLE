! Rather hacky file to avoid I/O in the direct ints for the newscf_f90
! driver. The iso stuff is particularly bad, as is the allocation
! methodology
Module f90_io_data

  Use newscf_numbers

  Implicit None

  Real( wp ), Dimension( : ), Allocatable :: schwartz_data
  Real( wp ), Dimension( : ), Allocatable :: iso_data

End Module f90_io_data

Subroutine put_schwartz( n, schwartz )


  Use newscf_numbers
  Use f90_io_data

  Implicit None

  Integer                     , Intent( In ) :: n
  Real( wp ), Dimension( 1:n ), Intent( In ) :: schwartz

  Integer :: i

  If( Allocated( schwartz_data ) ) Then
     Deallocate( schwartz_data )
  End If
  Allocate( schwartz_data( 1:n ) )

  Do i = 1, Size( schwartz_data )
     schwartz_data( i ) = schwartz( i )
  End Do

End Subroutine put_schwartz

Subroutine put_iso( n, iso )

  Use newscf_numbers
  Use f90_io_data

  Implicit None

  Integer                     , Intent( In ) :: n
  Real( wp ), Dimension( 1:n ), Intent( In ) :: iso

  Integer :: i

  If( Allocated( iso_data ) ) Then
     Deallocate( iso_data )
  End If
  Allocate( iso_data( 1:n ) )

  Do i = 1, Size( iso_data )
     iso_data( i ) = iso( i )
  End Do

End Subroutine put_iso

Subroutine get_schwartz( schwartz )

  Use newscf_numbers
  Use f90_io_data

  Implicit None

  Real( wp ), Dimension( * ), Intent( Out ) :: schwartz

  Integer :: i

  Do i = 1, Size( schwartz_data )
     schwartz( i ) = schwartz_data( i )
  End Do

End Subroutine get_schwartz

Subroutine get_iso( iso )

  Use newscf_numbers
  Use f90_io_data

  Implicit None

  Real( wp ), Dimension( * ), Intent( Out ) :: iso

  Integer :: i

  Do i = 1, Size( iso_data )
     iso( i ) = iso_data( i )
  End Do

End Subroutine get_iso


  

