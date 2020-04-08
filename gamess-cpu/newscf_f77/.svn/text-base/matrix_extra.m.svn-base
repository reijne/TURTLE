_IF(ma)
_MACRO(QM,qq(ivoff+`$1'))
_ENDIF
c
c Assignment to/from a matrix object using triangular
c arrays
c
      integer function matrix_set_from_triangle(matrix, data)
      implicit none

INCLUDE(matrix_internal.fh)

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount, i1, i2

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      icount = 0
      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            i1 = min(i,j)
            i2 = max(i,j)
            QM(ixa+(i-1) + lda*(j-1)) = data(iky(i2) + i1)
         enddo
      enddo
      matrix_set_from_triangle=0
      end

      integer function matrix_get_to_triangle(matrix, data)
      implicit none

INCLUDE(matrix_internal.fh)

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            if(i.le.j)then
               data(iky(j) + i) = QM(ixa+(i-1) + lda*(j-1))
            endif
         enddo
      enddo

      matrix_get_to_triangle=0

      end



c
c  dgemm -      C := alpha*A*B + beta*C
c 
c      integer function matrix_multiply(alpha,a,b,beta,c)

      integer function matrix_dgemm(x1,x2,
     &     alpha,
     &     a, b, beta, c)

      implicit none

      integer a, b, c
      double precision alpha, beta
      character*1 x1, x2

INCLUDE(matrix_internal.fh)

      integer nlink, ixa, ixb, ixc, ii

      ixa = hp(matrix_handle(a))
      ixb = hp(matrix_handle(b))
      ixc = hp(matrix_handle(c))

ccc
ccc   Note - quite a lot of work to do here if matrices are 
ccc   not all square!
ccc
ccc      if(x1 .eq. 'n' .and. x2 .eq. 'n')then

      nlink = dimensions(2,a)

      if(nlink .ne. dimensions(1,b))then
         write(6,'(a,i8,i8)')'Mismatch dimensions - nlink',
     &        dimensions(2,a),dimensions(1,b)
         call dlc_error('Matrix error','abort')
      endif

      if(dimensions(1,a) .ne. dimensions(1,c))then
         write(6,'(a,i8,i8)')'Mismatch dimensions - nrow(a) nrow(c)',
     &        dimensions(1,a),dimensions(1,c)
         call dlc_error('Matrix error','abort')
      endif

      if(dimensions(2,b) .ne. dimensions(2,c))then
         write(6,'(a,i8,i8)')'Mismatch dimensions - ncol(b) ncol(c)',
     &        dimensions(2,b),dimensions(2,c)
         call dlc_error('Matrix error','abort')
      endif

ccc      endif

      call dgemm(x1,x2
     &     ,dimensions(1,c),dimensions(2,c), nlink
     &     ,alpha, QM(ixa),dimensions(1,a)
     &     ,QM(ixb),dimensions(1,b)
     &     ,beta, QM(ixc),dimensions(1,c) )

      matrix_dgemm = 0
      end


c
c Similarity transform
c
      integer function matrix_mult2(Operator,Vectors,
     &     Result, Scratch)

      implicit none

c matrix tags
      integer Operator,Vectors,Result,Scratch
      integer idum1, idum2
      integer matrix_dgemm

      character*1 xn,xt
      data xn,xt/'n','t'/

      idum1 = matrix_dgemm(xt,xn,1.0d0,
     &     Vectors, Operator, 0.0d0, Scratch)

      idum2 = matrix_dgemm(xn,xn,1.0d0,
     &     Scratch,Vectors, 0.0d0, Result)

      matrix_mult2=0
      
      end

      integer function matrix_mult2t(Operator,Vectors,
     &     Result, Scratch)

      implicit none

c matrix tags
      integer Operator,Vectors,Result,Scratch
      integer idum1, idum2
      integer matrix_dgemm
      character*1 xn,xt
      data xn,xt/'n','t'/

      idum1= matrix_dgemm(xn,xn,1.0d0,
     &     Vectors, Operator, 0.0d0, Scratch)

      idum2 = matrix_dgemm(xn,xt,1.0d0,
     &     Scratch,Vectors, 0.0d0, Result)

      matrix_mult2t=idum1 + idum2

      end
