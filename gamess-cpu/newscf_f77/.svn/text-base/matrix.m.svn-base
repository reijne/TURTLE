_IF(ma)
_MACRO(QM,qq(ivoff+`$1'))
_ENDIF
c
c matrix object implementation
c pure f77 version, calling BLAS where available
c

c
c limitations
c    number of matrices allocated at any one time (max_matrices)
c    total memory available
c
c
c  Notes - fastest varying index of fortran workspace corresponds
c          to the row index (ie traversing a column)
c          a(n,m) has n rows and m cols
c          dimensions(1,ix) is n rows
c          dimensions(2,ix) is m cols
c
c     block data matrix_data
c     integer dum
c     include 'matrix_internal.fh'
c     parameter (dum=max_matrices*2)
c     data dimensions/dum * 0/
c     end
c
c returns handle of matrix, negative value denotes 
c failure
c

      integer function matrix_create(rows,columns,tag,type)
      implicit none
      integer rows, columns,handle_alloc
      external handle_alloc
      character tag*(*)
      integer type

_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF

INCLUDE(matrix_internal.fh)

      integer index
      integer elements
      integer i
c
c Locate the next free index
c
      index = 1

      do index = 1,max_matrices
         if (dimensions(1,index) .eq. 0) goto 10
      enddo
      
      write(6,'(a,i5)')
     &   'Too many matrices allocated',number_of_matrices
      write(6,'(a)')'Matrix names are:'
      do i=1,max_matrices
         write(6,'(a)')matrix_tag(i)
      enddo
      call dlc_error('Matrix error ','abort')

 10   continue

      number_of_matrices = number_of_matrices + 1

      dimensions(1,index)=rows
      dimensions(2,index)=columns
      elements = rows*columns
      matrix_tag(index) = tag
      matrix_type(index)= type
c---c     
c---c     Allocate memory for matrix
c---c     
c---c     igmem_alloc etc will allocate relative to the
c---c     common /vcore/qq
c---c
c---c     since we want to be able to include both matrix and array
c---c     files together, we cant use the same common block.
c---c     We therefore need an integer offset between /vcore/ and 
c---c     /matrix_core/
c---c
c---      if(offset.eq.0)then
c---         call gmem_c_pointer_diff(qq(1),mq(1),offset)
c---      endif
c---
         if (type .eq. MATRIX_REAL) then
c---
c---c         call mallocc(imq(1),imq(2),elements,
c---c     &        memory_offset(index),
c---c     &        memory_handle(index))
c---
c---         iaddr = igmem_alloc(elements)
c---c we'll use handle to release the memory
c---         memory_handle(index) = iaddr
c---c and offset to access it
c---         memory_offset(index) = iaddr + offset
c---         
c---c
c---c to check the offset is correct
c---c         call chkadr(qq(memory_handle(index)))
c---c         call chkadr(mq(memory_offset(index)))

         matrix_handle(index) = handle_alloc(SIZE_DOUBLE*elements)

c         write(6,*)'check address of matrix memory'
c         call chkadr(Q(hp(matrix_handle(index))))
c         call chkadr(QM(hp(matrix_handle(index))))

      else
         write(6,'(a,i8,a)')'Unrecognised matrix type for',tag,type
         call dlc_error('Matrix error ','abort')
      endif

c      write(6,*)'matrix_create',matrix_handle(index),tag

      matrix_create = index
      end

c-----------------------------------------------------------------------
c  matrix_delete
c  returns 0 success
c         -1 failure
c-----------------------------------------------------------------------

      integer function matrix_destroy(matrix)
      implicit none
INCLUDE(matrix_internal.fh)
      integer idum
      integer matrix
cc      integer handle_free
cc      external handle_free
      call handle_free(matrix_handle(matrix))
      number_of_matrices = number_of_matrices - 1
      dimensions(1,matrix)=0
      matrix_destroy=0
      end

c
c  simple matrix addition
c
      integer function matrix_add(ic,ia,ib)
      implicit none

INCLUDE(matrix_internal.fh)
      integer ia
      integer ib
      integer ic

      integer ixa, ixb, ixc, lda, ldb, ldc, i, j

      ixa = hp(matrix_handle(ia))
      ixb = hp(matrix_handle(ib))
      ixc = hp(matrix_handle(ic))

      lda = dimensions(1,ia)
      ldb = dimensions(1,ib)
      ldc = dimensions(1,ic)

      do j = 1, dimensions(2,ic)
         do i = 1, dimensions(1,ic)
            QM(ixc+(i-1) + ldc*(j-1)) = 
     &           QM(ixa+(i-1) + lda*(j-1)) +
     &           QM(ixb+(i-1) + ldb*(j-1)) 
         enddo
      enddo
      matrix_add = 0
      end

      integer function matrix_assign_unit(matrix)
      implicit none

INCLUDE(matrix_internal.fh)
      integer matrix
      integer matrix_assign_zero
      integer ixa, lda, i, j, idum

      idum = matrix_assign_zero(matrix)

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)

      do j = 1, dimensions(2,matrix)
         QM(ixa+(j-1) + lda*(j-1)) =  1.0
      enddo
      matrix_assign_unit = 0
      end
c

      integer function matrix_assign_zero(matrix)
      implicit none

INCLUDE(matrix_internal.fh)
      integer matrix

      integer ixa, lda, i, j

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)

      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            QM(ixa+(i-1) + lda*(j-1)) =  0.0
         enddo
      enddo
      matrix_assign_zero = 0
      end

c
c Assignment
c
      integer function matrix_set(matrix, data)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      icount = 0
      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            icount = icount + 1
            QM(ixa+(i-1) + lda*(j-1)) =  data(icount)
         enddo
      enddo
      matrix_set=0
      end

      integer function matrix_set_element(matrix, data, row, col)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data
      integer ixa, lda, row, col

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      QM(ixa+(row-1) + lda*(col-1)) =  data
      matrix_set_element=0
      end

      integer function matrix_get(matrix, data)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      icount = 0
      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            icount = icount + 1
            data(icount) = QM(ixa+(i-1) + lda*(j-1))
         enddo
      enddo

      matrix_get=0

      end

      integer function matrix_get_element(matrix, data, row, col)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data
      integer ixa, lda, row, col

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      data = QM(ixa+(row-1) + lda*(col-1))
      matrix_get_element=0
      end


c
c     Row by row and column by column for efficient transfer
c
      integer function matrix_set_column(matrix, data,col)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount,col

      ixa = hp(matrix_handle(matrix)) 
      lda = dimensions(1,matrix)
      ixa = ixa+lda*(col-1)-1
      do i = 1, dimensions(1,matrix)
         QM(ixa+i) = data(i)
      enddo

      matrix_set_column=0

      end

      integer function matrix_set_row(matrix, data,row)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount,row

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      ixa = row-1+ixa
      icount = 0
      do i = 0, dimensions(2,matrix)-1
         icount=icount+1
         QM(ixa+i*lda) = data(icount)
      enddo

      matrix_set_row=0

      end

      integer function matrix_get_row(matrix, data,row)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i, j, icount,row

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      ixa = row-1+ixa
      icount = 0
      do i = 0, dimensions(2,matrix)-1
         icount=icount+1
         data(icount) = QM(ixa+i*lda) 
      enddo

      matrix_get_row=0

      end

      integer function matrix_get_column(matrix, data,col)
      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      double precision data(*)
      integer ixa, lda, i ,col

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      ixa = ixa+lda*(col-1)-1
      do i = 1, dimensions(1,matrix)
         data(i)= QM(ixa+i) 
      enddo

      matrix_get_column=0

      end

c
c printing 
c
      integer function matrix_print(matrix)
      implicit none

INCLUDE(matrix_internal.fh)
      integer matrix

      integer ixa, lda, i, j
      logical ohi

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)

      write(6,'(a)') ' '
      write(6,'(a,a)')'Contents of matrix: ',matrix_tag(matrix)
      write(6,'(a)')   '---------------------------------------'

      ohi = .false.
      call prmat(QM(ixa),
     &     dimensions(1,matrix),dimensions(2,matrix),
     &     dimensions(1,matrix), ohi)

c      write(6,'(a,i8)')'Handle of matrix: ',matrix_handle(matrix)
      matrix_print = 0
      end

      integer function matrix_print_titled(matrix,title)
      implicit none

INCLUDE(matrix_internal.fh)
      integer matrix
      character title*(*)
      character t*120

      integer ixa, lda, i, j, l
      logical ohi

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)

      t = title
      l = 0
      do i = 120,1,-1
         if (t(i:i) .ne. ' ' .and. l .eq. 0)then
            l=i
         endif
      enddo

      write(6,'(a)') ' '
      write(6,'(a)')title
      write(6,'(120a1)')('-',i=1,l)

      ohi = .false.
      call prmat(QM(ixa),
     &     dimensions(1,matrix),dimensions(2,matrix),
     &     dimensions(1,matrix), ohi)

c      write(6,'(a,i8)')'Handle of matrix: ',matrix_handle(matrix)
      matrix_print_titled = 0
      end


      subroutine prmat(v,m,n,ndim,ohi)
c
c     ----- print out a matrix -----
c
c m = number of columns
c n = number of rows
c ndim = leading dimension
c ohi print high precision
c
      implicit none
c
      integer n,m,ndim
      double precision v(ndim,*)
      logical ohi
      integer imin, imax, max, i, j, iwr
      iwr = 6
      ohi = .false.

      max = 12
      if (ohi) max=7
      imax = 0

 100  imin = imax+1
      imax = imax+max
      if (imax .gt. n) imax = n
      write (iwr,9008)
      if(.not. ohi) write (iwr,8028) (i,i = imin,imax)
      if( ohi) write (iwr,9028) (i,i = imin,imax)
      write (iwr,9008)
      do 120 j = 1,m
      if( ohi) write (iwr,9048) j,(v(j,i),i = imin,imax)
      if(.not.ohi) write (iwr,8048) j,(v(j,i),i = imin,imax)
120   continue
      if (imax .lt. n) go to 100
      return
 9008 format(1x)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.5)
      end
c
c actually computes square root of sum of squares 
c of all elements of a 2d matrix
c
      double precision function matrix_length(matrix)

      implicit none

INCLUDE(matrix_internal.fh)

      integer matrix
      integer ixa, lda, i, j, icount
      double precision sum, fac, mmq

      ixa = hp(matrix_handle(matrix))
      lda = dimensions(1,matrix)
      icount = 0
      sum = 0.0d0
      do j = 1, dimensions(2,matrix)
         do i = 1, dimensions(1,matrix)
            icount = icount + 1
            mmq=QM(ixa+(i-1) + lda*(j-1))
            sum=sum+mmq*mmq
         enddo
      enddo
c      fac = icount
      matrix_length= sqrt(sum)

      end
c
c  implements ddot:  x.y
c
      double precision function matrix_dot_product(x,y)

      implicit none

      integer x, y

INCLUDE(matrix_internal.fh)

      integer ixx, ixy, length
_IF(t3e,cray)
      REAL `sdot'
_ELSE
      REAL ddot
_ENDIF

      ixx = hp(matrix_handle(x))
      ixy = hp(matrix_handle(y))
      length = dimensions(1,x) * dimensions(2,x)
      matrix_dot_product = ddot(length
     &     ,QM(ixx),1,QM(ixy),1)
      end

c
c implements daxpy : y <- y+a.x
c warning - assumes matrices are the same shape
c
      integer function matrix_daxpy(a,x,y)
      implicit none

      double precision a
      integer x,y

INCLUDE(matrix_internal.fh)

      integer ixx, ixy, length

      ixx = hp(matrix_handle(x))
      ixy = hp(matrix_handle(y))

      length = dimensions(1,x) * dimensions(2,x)
      call daxpy(length,a
     &     ,QM(ixx),1,QM(ixy),1)

      matrix_daxpy = 0

      end

c
c implements dcopy : y <- x
c warning - assumes matrices are the same shape
c
      integer function matrix_copy(x,y)
      implicit none

      integer x,y

INCLUDE(matrix_internal.fh)

      integer length, ixx, ixy
      
      ixx = hp(matrix_handle(x))
      ixy = hp(matrix_handle(y))
      length = dimensions(1,x) * dimensions(2,x)
      call dcopy(length
     &     ,QM(ixx),1,QM(ixy),1)

      matrix_copy = 0

      end

c
c implements dscal:  x <- a.x
c
      integer function matrix_scale(x,a)
      implicit none

      integer x
      double precision a

INCLUDE(matrix_internal.fh)

      integer length, ixx
      
      ixx = hp(matrix_handle(x))
      length = dimensions(1,x) * dimensions(2,x)
      call dscal(length,a,QM(ixx),1)

      matrix_scale = 0
      end



      integer function matrix_dimension(x,index)

      implicit none

      integer x
      integer index 

INCLUDE(matrix_internal.fh)

      matrix_dimension = dimensions(index,x)

      end
c
c  dgemm -      C := alpha*A*B + beta*C
c 
      integer function matrix_multiply(alpha,a,b,beta,c)
      implicit none

      integer a, b, c
      double precision alpha, beta

INCLUDE(matrix_internal.fh)

      integer nlink, ixa, ixb, ixc, ii
      character*1 xn
      data xn/'n'/

      ixa = hp(matrix_handle(a))
      ixb = hp(matrix_handle(b))
      ixc = hp(matrix_handle(c))

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

      call dgemm(xn,xn
     &     ,dimensions(1,c),dimensions(2,c), nlink
     &     ,alpha, QM(ixa),dimensions(1,a)
     &     ,QM(ixb),dimensions(1,b)
     &     ,beta, QM(ixc),dimensions(1,c) )

      matrix_multiply = 0
      end
c
c  implement c <- alpha.a + beta.b    for alpha,beta scalar
c
      integer function matrix_combine(alpha,a,beta,b,c)
      implicit none

      integer a, b, c
      double precision alpha, beta

INCLUDE(matrix_internal.fh)

      integer ixa, ixb, ixc, length, i

      ixa = hp(matrix_handle(a))
      ixb = hp(matrix_handle(b))
      ixc = hp(matrix_handle(c))

      length = dimensions(1,a)*dimensions(2,a)
      do i = 1,length
         QM(ixc + (i-1) ) = alpha * 
     &        QM(ixa + (i-1) ) + 
     &        beta * QM(ixb + (i-1) )
      enddo

      matrix_combine = 0

      end

      double precision function matrix_absmax(x)

      implicit none

      integer x

INCLUDE(matrix_internal.fh)

      integer ixx, ldx, i, j
      double precision t

      ixx = hp(matrix_handle(x))
      ldx = dimensions(1,x)

      t = abs(QM(ixx))

      do j = 1, dimensions(2,x)
         do i = 1, dimensions(1,x)
            t = max(t,abs(QM(ixx+(i-1) + ldx*(j-1))))
         enddo
      enddo

      matrix_absmax  = t

      end

      integer function matrix_transpose(x)

      implicit none

      integer itag


INCLUDE(matrix_internal.fh)

      integer idum,i,j,x,ip
      integer t1, t2, ixt1, ixt2, ixx, ldax, istart, n

      double precision rdum

      integer matrix_create
      external matrix_create

      integer matrix_copy
      external matrix_copy

      integer matrix_print
      external matrix_print

      integer matrix_destroy
      external matrix_destroy

      integer stack_alloc
      external stack_alloc

      integer stack_free 
      external stack_free 

      integer matrix_get_row, matrix_get_column
      external matrix_get_row, matrix_get_column

      integer matrix_set_column, matrix_set_row
      external matrix_set_column, matrix_set_row

      if( (dimensions(1,x) .eq. 1) .or.
     &    (dimensions(2,x) .eq. 1) )then
c
c swap around the dimensions
c
         idum = dimensions(1,x)
         dimensions(1,x) = dimensions(2,x)
         dimensions(2,x) = idum
      else if(dimensions(1,x) .eq. dimensions(2,x)) then
c
c square case, swap elements using a scratch row
c 
         n      = dimensions(1,x)   
         ip     = stack_alloc(SIZE_DOUBLE*n)
         ixx    = hp(matrix_handle(x))
         istart = 1
         do i = 1, n
            istart = istart + 1
c load a column
            do j = istart, n
               QM(ip +(j-1)) =  
     &              QM(ixx + (i-1)*n + (j-1) )
            enddo
c replace with a row
            do j = istart, n
               QM(ixx + (i-1)*n +(j-1)) = 
     &              QM(ixx +(j-1)*n +(i-1))
            enddo
c re store saved column row-wise
            do j = istart, n
               QM(ixx + (j-1)*n + (i-1) ) = 
     &              QM(ip +(j-1)) 
            enddo
         enddo
         ip=  stack_free(ip)
      else
c
c replace dimensions and reshape matrix. There is probably a better
c way !
c
         t1   = matrix_create(dimensions(2,x),1,'scratch v',MATRIX_REAL)
         t2   = matrix_create(dimensions(2,x),dimensions(1,x),
     &        'scratch m',MATRIX_REAL)
         ixt1 = hp(matrix_handle(t1))
         do i = 1, dimensions(1,x)
            idum=matrix_get_row   (x, QM(ixt1),i) 
            idum=matrix_set_column(t2,QM(ixt1),i) 
         enddo
c
c Reshape input matrix and copy transposed data
c
         idum = dimensions(1,x)
         dimensions(1,x) = dimensions(2,x)
         dimensions(2,x) = idum

         idum = matrix_copy(t2,x)

         idum = matrix_destroy(t2)
         idum = matrix_destroy(t1)
      endif
      matrix_transpose=0
      end
c
c  Random fill of a matrix
c
      integer function matrix_set_random(itag)
      implicit none
      integer itag
INCLUDE(matrix_internal.fh)

      integer r1, r2, i, j, ixa, lda
      double precision rand, dum

      call srand(1234)

      ixa = hp(matrix_handle(itag))
      lda = dimensions(1,itag)

      do j = 1, dimensions(2,itag)
         do i = 1, dimensions(1,itag)
            QM(ixa + (i-1) + lda*(j-1)) = rand()
         enddo
      enddo
      matrix_set_random = 0
      end

c
c     Delete the record of a matrix but don't deallocate
c     return handle to matrix memory
c
      integer function matrix_to_array(itag)
      implicit none
INCLUDE(matrix_internal.fh)
      integer idum,itag
      number_of_matrices = number_of_matrices - 1
      dimensions(1,itag)=0
      matrix_to_array=matrix_handle(itag)
      end

c
c     Create a new matrix, but give its data the handel 
c     supplied value to it directly picts up contents
c     from an array
c
c
      integer function
     &     matrix_from_array(rows,columns,tag,type,point)
      implicit none
      integer rows, columns
      character tag*(*)
      integer type,point

INCLUDE(matrix_internal.fh)

      integer index
      integer elements
      integer iaddr
c
c Locate the next free index
c
      index = 1

      do index = 1,max_matrices
         if (dimensions(1,index) .eq. 0) goto 10
      enddo

      call dlc_error('Too many matrices allocate','abort')

 10   continue

      number_of_matrices = number_of_matrices + 1

      dimensions(1,index)=rows
      dimensions(2,index)=columns
      elements = rows*columns
      matrix_tag(index) = tag
      matrix_type(index)= type
      if (type .eq. MATRIX_REAL) then
         matrix_handle(index) = point
      else
         write(6,'(a,i8,i8)')'Unrecognised matrix type for',tag,type
         call dlc_error('Matrix','abort')
      endif

      matrix_from_array = index
      end

c
c matrix_f77b
c

c###################################################################
c
c                  matrix diagonalisation module
c

c
c     Driver function interfaces rsp with Quasi environemt     
c

      integer function matrix_diagonalise
     1(itag,ivect,ivalues,nvect,increasing)
      implicit none
INCLUDE(matrix_internal.fh)
      logical increasing
      integer itag,ivect,ivalues,nvect,ierr,iivect,iroot
      integer pi,pr,n,ipi,ipr,j,k,l,m,ipa,iitag,pa,i,eight
      integer matrix_create,
     &   matrix_destroy,stack_alloc,zero
      external matrix_create, matrix_destroy,stack_alloc
      integer stack_free
      external stack_free
c
c     zero is set to two as the eign vectors in dlc_managera are upside down
c     this is a nasty hack and must bee cleaned up!
c
      parameter(eight=8,zero=0)

      if(.not.increasing)stop 'matrix_diagonalise'

c     Make memory availible
      n=dimensions(1,itag)
      if(dimensions(2,itag).ne.n)then
         write(6,'(a)')'Matrix not square in matrix_diagonalise'
         matrix_diagonalise=-1
         return 
      endif

      ipi=stack_alloc(SIZE_DOUBLE*n)
      ipr=stack_alloc(SIZE_DOUBLE*n*eight)
      ipa=stack_alloc(SIZE_DOUBLE*n*n)
      iitag= hp(matrix_handle(itag))
      iivect=hp(matrix_handle(ivect))
      iroot= hp(matrix_handle(ivalues))

c     Pack the square input matrix into a lower
c     triangle matrix
      k=-1
      do i=1,n
         l=iitag+i-1
         do j=1,i
           k=k+1
           QM(ipa+k)=QM((j-1)*n+l)
         enddo
      enddo
      
c     do calculation
      i=0
      call evvrsp(i,n,nvect,n*n,n,QM(ipa),QM(ipr)
     &,QM(ipi),QM(iroot),QM(iivect),zero,ierr)

c     clear working room
      ipa= stack_free(ipa)
      ipr= stack_free(ipr)
      ipi= stack_free(ipi)

      if(ierr.ne.0)
     &   write(6,'(A,I5)')
     &'Matrix_diagonalise failed, returning ',ierr
      matrix_diagonalise=ierr
      return
      end

c
c                       end module
c
c###################################################################



c###################################################################
      double precision function matrix_symmetrise(itag)
c
c     11 nov 1996 - ajt
c     Converted to 'quasi' environemt' Oct 1997
c     symmeterise a square matrix and return
c     (rms asymetery)/(rms magnitude)
c
c     on input
c     --------
c     itag=tag of matrix

c     on output
c     ---------
c     return value is the error
c
      implicit none
INCLUDE(matrix_internal.fh)
      double precision mag,asy
      integer  n,r,i,j,x,y,itag,ix,nel
      r=dimensions(itag,1)
      n=dimensions(itag,2)
      if(r.ne.n)then
         write(6,'(a)')'Matrix not square in matrix_symmeterise'
         matrix_symmetrise=-1.0d0
         return
      endif
 
      mag=0.0d0
      asy=0.0d0
      ix = hp(matrix_handle(itag))
      do 100 i=1,n
         do 100 j=1,i
            x=ix + (j-1)*r + (i-1)
            y=ix + (i-1)*r + (j-1)
            asy = asy + (QM(x)-QM(y))**2
            QM(x)=(QM(x)+QM(y))*5d-1
            QM(y)=QM(x)
            mag=mag+QM(x)**2
 100  continue
c
c from paul - are you sure we need the n*(n+1)/2 factor
c I would have expected it cancels from top and bottom?
c
      matrix_symmetrise= sqrt(asy / mag)
      return
      end
c###################################################################



c###################################################################
c
c      module to invert a general square matrix
c
c      contains:
c           matrix_invert - driver
c              dgedi      - inverter
c              dgefa      - finds factors

      integer function matrix_invert(itag,d)
c
c     AJT 1997 - drives matrix inversion
c     runs in 'Quasi' environemt
c     input
c       itag = tag for memory offset of tag
c       det  = dble varible to which the determinant
c              will be returned
c
c     rotine dimanically allocates memory
c     for the various arrays used by the inverter
c     and compresses two double precision determinant to 
c     a simgle double precision determinent
c

      implicit none
INCLUDE(matrix_internal.fh)
      double precision  d, dd(2)
      integer   n,ni,r,job,ns,itag,info
      integer   p1,p2,of,i,stack_alloc
      external stack_alloc
      integer stack_free
      external stack_free

      integer lenwrd
      external lenwrd

c      call timePrint('Start-invert')
      job=11

      n=dimensions(1,itag)
      if(dimensions(2,itag).ne.n)then
         call dlc_error('Matrix not square in invert','abort')
         matrix_invert=-1
         return
      endif
      ns=n*n
      r=n
      
c     get memory
      p1=stack_alloc(SIZE_DOUBLE*n)
      ni = 1+(n -1)/lenwrd()
      p2=stack_alloc(SIZE_DOUBLE*ni)

      if((p1.eq.-1).or.(p2.eq.-1)) stop
      of=hp(matrix_handle(itag))
      call minvs(QM(of),r,n,QM(p2),
     &     dd,QM(p1),job)

c     return memory     
      p2 = stack_free(p2)
      p1=  stack_free(p1)

c     recompute determinant
c
c from paul - i copied this from below
c       determinant = det(1) * 10.0**det(2)
      d=dd(1) * 10.0d0**dd(2)

      matrix_invert=0
c      call timePrint('End-invert')
      return
      end

      subroutine minvs(m,r,n,ipvt,dd,work,job)
      double precision m(*),dd(*),work(*)
      integer r,n,ipvt(*),job,info
c     get factors
      call dgefa(m,r,n,ipvt,info)
c     test for singularity
      if (info.ne.0) then
        dd(1)=0d0
        dd(2)=0d0
        write(6,'(a)')
     1'warning => attempt to invert (probably) sigular matrix'
        write(6,'(a)')
     1'matrix is left unchanged and determinent set to zero'
        return
      end if
      call dgedi(m,r,n,ipvt,dd,work,job)
      return
      end
c
c                end module
c
c###################################################################


c  lapack.f

      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      implicit double precision(a-h,o-z)
      dimension a(lda,1),det(2),work(1),ipvt(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer*4
c                the leading dimension of the array  a .
c
c        n       integer*4
c                the order of the matrix  a .
c
c        ipvt    integer*4(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer*4
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. abs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran abs,mod
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
      det(1) = 1.0d+00
      det(2) = 0.0d+00
      ten = 10.0d+00
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = a(i,i)*det(1)
c        ...exit
         if (det(1) .eq. 0.0d+00) go to 60
   10    if (abs(det(1)) .ge. 1.0d+00) go to 20
         if (det(2).lt.-1000.0d0) then
             call dlc_error('Problems getting determinant','warn')
             det(1)=0.0d0
             det(2)=0.0d0
             return 
         endif
         det(1) = ten*det(1)
         det(2) = det(2) - 1.0d+00
         go to 10
   20    continue
   30    if (abs(det(1)) .lt. ten) go to 40
         det(1) = det(1)/ten
         det(2) = det(2) + 1.0d+00
         go to 30
   40    continue
   50 continue
   60 continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
      do 100 k = 1, n
         a(k,k) = 1.0d+00/a(k,k)
         t = -a(k,k)
         call dscal(k-1,t,a(1,k),1)
         kp1 = k + 1
         if (n .lt. kp1) go to 90
         do 80 j = kp1, n
            t = a(k,j)
            a(k,j) = 0.0d+00
            call daxpy(k,t,a(1,k),1,a(1,j),1)
   80    continue
   90    continue
  100 continue
c
c        form inverse(u)*inverse(l)
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 140
      do 130 kb = 1, nm1
         k = n - kb
         kp1 = k + 1
         do 110 i = kp1, n
            work(i) = a(i,k)
            a(i,k) = 0.0d+00
  110    continue
         do 120 j = kp1, n
            t = work(j)
            call daxpy(n,t,a(1,j),1,a(1,k),1)
  120    continue
         l = ipvt(k)
         if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130 continue
  140 continue
  150 continue
      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
      implicit double precision(a-h,o-z)
      dimension a(lda,1),ipvt(1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer*4
c                the leading dimension of the array  a .
c
c        n       integer*4
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer*4(n)
c                an integer*4 vector of pivot indices.
c
c        info    integer*4
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d+00) go to 40
c
c           interchange if necessary
c
         if (l .eq. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
   10    continue
c
c           compute multipliers
c
         t = -1.0d+00/a(k,k)
         call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
         do 30 j = kp1, n
            t = a(l,j)
            if (l .eq. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
   20       continue
            call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30    continue
         go to 50
   40    continue
         info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d+00) info = n
      return
      end

      subroutine matrix_list
      implicit none
INCLUDE(matrix_internal.fh)
      logical opg_root
      integer index
      if(opg_root())then
         write(6,*)'Current list of matrices'
         do index = 1,max_matrices
            if (dimensions(1,index) .ne. 0) then
               write(6,*)matrix_tag(index)
            endif
         enddo
      endif
      end

      subroutine dlc_error(s1,s2)
      character s1*(*), s2*(*)
INCLUDE(../m4/common/errcodes)
      call gamerr(s1,
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_SYNC, ERR_NO_SYS)
      end
