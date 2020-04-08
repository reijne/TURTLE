      subroutine `dcopy( n, dx, incx, dy, incy )'
      integer           incx, incy, n
      REAL  dx( 1 ), dy( 1 )
      call scopy( n, dx, incx, dy, incy )
      end

      subroutine `daxpy( n, a, x, incx, y, incy )'
      integer           incx, incy, n
      REAL  a
      REAL  x( 1 ), y( 1 )
      call saxpy(  n, a, x, incx, y, incy )
      end

      subroutine `dscal( n, a, x, incx )'
      integer           incx, n
      REAL  a
      REAL  x( 1 )
      call sscal( n, a, x, incx )
      end

      subroutine `dgemm( transa,transb,m,n,k,alpha,a,lda,b,ldb,
     &     beta,c,ldc )'
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      REAL   alpha, beta
      REAL   a( lda,*), b(ldb,*), c(ldc,*)
      call sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc )
      end

      REAL function  `ddot( n, x, incx, y, incy )'
      integer           incx, incy, n
      REAL  x( 1 ), y( 1 ), sdot
      `ddot' = sdot( n, x, incx, y, incy )
      end

      subroutine `drot( n, dx, incx, dy, incy, c, s )'
      integer          incx, incy, n
      REAL c, s
      REAL dx( 1 ), dy( 1 )
      call srot ( n, dx, incx, dy, incy, c, s )
      end

      REAL  function `dnrm2( n, x, incx )'
      integer           incx, n
      REAL  x( 1 ), snrm2
      `dnrm2' = snrm2( n, x, incx )
      end

      subroutine `dswap( n, dx, incx, dy, incy )'
      integer           incx, incy, n
      REAL  dx( 1 ), dy( 1 )
      call sswap (n, dx, incx, dy, incy )
      end

      integer function `idamax( n, x, incx )'
      integer           incx, n
      integer isamax
      REAL  x( 1 )
      `idamax' = isamax( n, x, incx )
      end
c
      REAL function `dlamch(cmach)'
      character *1  cmach
      REAL slamch
      `dlamch' = slamch(cmach)
      end
