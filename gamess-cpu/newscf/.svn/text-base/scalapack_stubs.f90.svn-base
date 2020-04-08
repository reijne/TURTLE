! Now the BLACS stuff

Logical Function scalapack_avail()

  ! Indicate if Scalapack is available in this build

  scalapack_avail = .False.
  
End Function scalapack_avail

Subroutine blacs_get( context, what, val )

  Implicit None

  Integer, Intent( In    ) :: context
  Integer, Intent( In    ) :: what
  Integer, Intent(   Out ) :: val

  val = 1

End Subroutine blacs_get

Subroutine blacs_gridexit( context )

  Integer, Intent( In ) :: context

End Subroutine blacs_gridexit

Subroutine blacs_gridinfo( context, nprow, npcol, myrow, mycol )

  Implicit None

  Integer, Intent( In    ) :: context
  Integer, Intent(   Out ) :: nprow
  Integer, Intent(   Out ) :: npcol
  Integer, Intent(   Out ) :: myrow
  Integer, Intent(   Out ) :: mycol

  nprow = 1
  npcol = 1
  myrow = 0
  mycol = 0

End Subroutine blacs_gridinfo

Subroutine blacs_gridmap( context, map, n, nprow, npcol )

  Implicit None

  Integer,                            Intent( In ) :: context
  Integer,                            Intent( In ) :: n
  Integer,                            Intent( In ) :: npcol
  Integer, Dimension( 1:n, 1:npcol ), Intent( In ) :: map
  Integer,                            Intent( In ) :: nprow

End Subroutine blacs_gridmap

! Finally the SCALAPACK stuff
! Firstly the computational routines

Subroutine pdgemm( ta, tb, n, m, k, alpha, a, npa, nra, desca, &
                                           b, npb, nrb, descb, &
                                     beta, c, npc, nrc, descc )

  Use newscf_numbers

  Character( Len = * )        , Intent( In    ) :: ta
  Character( Len = * )        , Intent( In    ) :: tb
  Integer                     , Intent( In    ) :: n
  Integer                     , Intent( In    ) :: m
  Integer                     , Intent( In    ) :: k
  Real( wp )                  , Intent( In    ) :: alpha
  Real( wp ), Dimension( 1:* ), Intent( In    ) :: a
  Integer                     , Intent( In    ) :: npa
  Integer                     , Intent( In    ) :: nra
  Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
  Real( wp ), Dimension( 1:* ), Intent( In    ) :: b
  Integer                     , Intent( In    ) :: npb
  Integer                     , Intent( In    ) :: nrb
  Integer   , Dimension( 1:9 ), Intent( In    ) :: descb
  Real( wp )                  , Intent( In    ) :: beta
  Real( wp ), Dimension( 1:* ), Intent( InOut ) :: c
  Integer                     , Intent( In    ) :: npc
  Integer                     , Intent( In    ) :: nrc
  Integer   , Dimension( 1:9 ), Intent( In    ) :: descc

  Call dgemm( ta, tb, n, m, k, alpha, a, desca( 9 ), &
                                      b, descb( 9 ), &
                                beta, c, descc( 9 ) )

End Subroutine pdgemm

Subroutine pzgemm( ta, tb, n, m, k, alpha, a, npa, nra, desca, &
                                           b, npb, nrb, descb, &
                                     beta, c, npc, nrc, descc )

  Use newscf_numbers

  Implicit None

  Character( Len = * )           , Intent( In    ) :: ta
  Character( Len = * )           , Intent( In    ) :: tb
  Integer                        , Intent( In    ) :: n
  Integer                        , Intent( In    ) :: m
  Integer                        , Intent( In    ) :: k
  Complex( wp )                  , Intent( In    ) :: alpha
  Complex( wp ), Dimension( 1:* ), Intent( In    ) :: a
  Integer                        , Intent( In    ) :: npa
  Integer                        , Intent( In    ) :: nra
  Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
  Complex( wp ), Dimension( 1:* ), Intent( In    ) :: b
  Integer                        , Intent( In    ) :: npb
  Integer                        , Intent( In    ) :: nrb
  Integer      , Dimension( 1:9 ), Intent( In    ) :: descb
  Complex( wp )                  , Intent( In    ) :: beta
  Complex( wp ), Dimension( 1:* ), Intent( InOut ) :: c
  Integer                        , Intent( In    ) :: npc
  Integer                        , Intent( In    ) :: nrc
  Integer      , Dimension( 1:9 ), Intent( In    ) :: descc

  Call zgemm( ta, tb, n, m, k, alpha, a, desca( 9 ), &
                                      b, descb( 9 ), &
                                beta, c, descc( 9 ) )

End Subroutine pzgemm

Subroutine pdgemv( ta, n, m, alpha, a, npa, nra, desca,       &
                                    x, npx, nrx, descx, incx, &
                              beta, y, npy, nry, descy, incy )

  Use newscf_numbers

  Implicit None

  Character( Len = * )        , Intent( In    ) :: ta
  Integer                     , Intent( In    ) :: n
  Integer                     , Intent( In    ) :: m
  Real( wp )                  , Intent( In    ) :: alpha
  Real( wp ), Dimension( 1:* ), Intent( In    ) :: a
  Integer                     , Intent( In    ) :: npa
  Integer                     , Intent( In    ) :: nra
  Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
  Real( wp ), Dimension( 1:* ), Intent( In    ) :: x
  Integer                     , Intent( In    ) :: npx
  Integer                     , Intent( In    ) :: nrx
  Integer   , Dimension( 1:9 ), Intent( In    ) :: descx
  Integer                     , Intent( In    ) :: incx
  Real( wp )                  , Intent( In    ) :: beta
  Real( wp ), Dimension( 1:* ), Intent( InOut ) :: y
  Integer                     , Intent( In    ) :: npy
  Integer                     , Intent( In    ) :: nry
  Integer   , Dimension( 1:9 ), Intent( In    ) :: descy
  Integer                     , Intent( In    ) :: incy

  Call dgemv( ta, n, m, alpha, a, desca( 9 ), &
                               x, incx,       &
                         beta, y, incy )

End Subroutine pdgemv

! Not used at present as locking broken
Subroutine pdgetrf( m, n, a, ia, ja, desca, ipiv, info )

  Use newscf_numbers

  Implicit None

  Integer                     , Intent( In    ) :: m
  Integer                     , Intent( In    ) :: n
  Real( wp ), Dimension( 1:* ), Intent( InOut ) :: a
  Integer                     , Intent( In    ) :: ia
  Integer                     , Intent( In    ) :: ja
  Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
  Integer   , Dimension( 1:* ), Intent( In    ) :: ipiv
  Integer                     , Intent(   Out ) :: info

  info = 0

End Subroutine pdgetrf

! Not used at present as locking broken
Subroutine pdgetri( n, a, ia, ja, desca, ipiv, work, lwork, iwork, liwork, info )

  Use newscf_numbers

  Implicit None

  Integer                          , Intent( In    ) :: n
  Real( wp ), Dimension( 1:*      ), Intent( InOut ) :: a
  Integer                          , Intent( In    ) :: ia
  Integer                          , Intent( In    ) :: ja
  Integer   , Dimension( 1:9      ), Intent( In    ) :: desca
  Integer   , Dimension( 1:*      ), Intent( In    ) :: ipiv
  Integer                          , Intent( In    ) :: lwork
  Real( wp ), Dimension( 1:lwork  ), Intent( InOut ) :: work
  Integer                          , Intent( In    ) :: liwork
  Integer   , Dimension( 1:liwork ), Intent( In    ) :: iwork
  Integer                          , Intent(   Out ) :: info

  info = 0

End Subroutine pdgetri

Subroutine pdsyev( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, &
     work, lwork, info )

  Use newscf_numbers

  Implicit None

  Character( Len = * )             , Intent( In    ) :: jobz
  Character( Len = * )             , Intent( In    ) :: uplo
  Integer                          , Intent( In    ) :: n
  Real( wp ), Dimension( 1:*      ), Intent( InOut ) :: a
  Integer                          , Intent( In    ) :: ia
  Integer                          , Intent( In    ) :: ja
  Integer   , Dimension( 1:9      ), Intent( In    ) :: desca
  Real( wp ), Dimension( 1:*      ), Intent(   Out ) :: w
  Real( wp ), Dimension( 1:*      ), Intent(   Out ) :: z
  Integer                          , Intent( In    ) :: iz
  Integer                          , Intent( In    ) :: jz
  Integer   , Dimension( 1:9      ), Intent( In    ) :: descz
  Integer                          , Intent( In    ) :: lwork
  Real( wp ), Dimension( 1:lwork  ), Intent( InOut ) :: work
  Integer                          , Intent(   Out ) :: info

  Call dsyev( jobz, uplo, n, a, desca( 9 ), w, work, lwork, info )

End Subroutine pdsyev

Subroutine pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, &
     work, lwork, iwork, liwork, info )

  Use newscf_numbers

  Implicit None

  Character( Len = * )             , Intent( In    ) :: jobz
  Character( Len = * )             , Intent( In    ) :: uplo
  Integer                          , Intent( In    ) :: n
  Real( wp ), Dimension( 1:*      ), Intent( InOut ) :: a
  Integer                          , Intent( In    ) :: ia
  Integer                          , Intent( In    ) :: ja
  Integer   , Dimension( 1:9      ), Intent( In    ) :: desca
  Real( wp ), Dimension( 1:*      ), Intent(   Out ) :: w
  Real( wp ), Dimension( 1:*      ), Intent(   Out ) :: z
  Integer                          , Intent( In    ) :: iz
  Integer                          , Intent( In    ) :: jz
  Integer   , Dimension( 1:9      ), Intent( In    ) :: descz
  Integer                          , Intent( In    ) :: lwork
  Real( wp ), Dimension( 1:lwork  ), Intent( InOut ) :: work
  Integer                          , Intent( In    ) :: liwork
  Integer   , Dimension( 1:liwork ), Intent( In    ) :: iwork
  Integer                          , Intent(   Out ) :: info

  Call pdsyev( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, &
     work, lwork, info )

End Subroutine pdsyevd

Subroutine pdelget( scope, top, alpha, a, ia, ja, desca )

  Use newscf_numbers

  Implicit None

  Character( Len = * )             , Intent( In    ) :: scope
  Character( Len = * )             , Intent( In    ) :: top
  Real( wp )                       , Intent(   Out ) :: alpha
  Real( wp ), Dimension( 1:*      ), Intent( In    ) :: a
  Integer                          , Intent( In    ) :: ia
  Integer                          , Intent( In    ) :: ja
  Integer   , Dimension( 1:9      ), Intent( In    ) :: desca

  alpha = a( ( ja - 1 ) * desca( 9 ) + ia )

End Subroutine pdelget

Subroutine pdgemr2d( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt )

  Use newscf_numbers

  Implicit None

  Integer                          , Intent( In    ) :: m
  Integer                          , Intent( In    ) :: n
  Real( wp ), Dimension( 1:*      ), Intent( In    ) :: a
  Integer                          , Intent( In    ) :: ia
  Integer                          , Intent( In    ) :: ja
  Integer   , Dimension( 1:9      ), Intent( In    ) :: desca
  Real( wp ), Dimension( 1:*      ), Intent(   Out ) :: b
  Integer                          , Intent( In    ) :: ib
  Integer                          , Intent( In    ) :: jb
  Integer   , Dimension( 1:9      ), Intent( In    ) :: descb
  Integer                          , Intent( In    ) :: ictxt

  Integer :: lda, ldb
  Integer :: posa, posb
  Integer :: starta, startb
  Integer :: i, j

  lda = desca( 9 )
  ldb = descb( 9 )

  posa = ( ja - 1 ) * lda + ia - 1
  posb = ( jb - 1 ) * ldb + ib - 1

  Do i = 1, n
     starta = posa
     startb = posb
     Do j = 1, m
        posa = posa + 1
        posb = posb + 1
        b( posb ) = a( posa )
     End Do
     posa = starta + lda
     posb = startb + ldb
  End Do

End Subroutine pdgemr2d

Integer Function indxg2p( indxglob, nb, iproc, isrcproc, nprocs )

  Use newscf_numbers

  Implicit None

  Integer, Intent( In ) :: indxglob
  Integer, Intent( In ) :: nb
  Integer, Intent( In ) :: iproc
  Integer, Intent( In ) :: isrcproc
  Integer, Intent( In ) :: nprocs

  indxg2p = 0

End Function indxg2p

Integer Function indxg2l( indxglob, nb, iproc, isrcproc, nprocs )

  Use newscf_numbers

  Implicit None

  Integer, Intent( In ) :: indxglob
  Integer, Intent( In ) :: nb
  Integer, Intent( In ) :: iproc
  Integer, Intent( In ) :: isrcproc
  Integer, Intent( In ) :: nprocs

  indxg2l = indxglob

End Function indxg2l

Subroutine descinit( desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info )

  Use newscf_numbers

  Implicit None

  Integer, Dimension( 1:9 ), Intent(   Out ) :: desc
  Integer                  , Intent( In    ) :: m
  Integer                  , Intent( In    ) :: n
  Integer                  , Intent( In    ) :: mb
  Integer                  , Intent( In    ) :: nb
  Integer                  , Intent( In    ) :: irsrc
  Integer                  , Intent( In    ) :: icsrc
  Integer                  , Intent( In    ) :: ictxt
  Integer                  , Intent( In    ) :: lld
  Integer                  , Intent(   Out ) :: info

  info = 0

  desc( 1 ) = 1
  desc( 2 ) = m
  desc( 3 ) = n
  desc( 4 ) = mb
  desc( 5 ) = nb
  desc( 6 ) = 1
  desc( 7 ) = 1
  desc( 8 ) = ictxt
  desc( 9 ) = lld

End Subroutine descinit

Integer function numroc( n, nb, iproc, isrcproc, nprocs )

  Use newscf_numbers

  Implicit None

  Integer, Intent( In ) :: n
  Integer, Intent( In ) :: nb
  Integer, Intent( In ) :: iproc
  Integer, Intent( In ) :: isrcproc
  Integer, Intent( In ) :: nprocs

  numroc = n

End function numroc
