_IF(ma)
_MACRO(QM,qq(ivoff+`$1'))
_ENDIF
c************************************************************************
      integer function allocate_vectors(nbas,uhf)
c************************************************************************
      implicit none
      integer nbas
      logical uhf

INCLUDE(vectors_internal.fh)
INCLUDE(matrix_internal.fh)

      include 'matrix_interface.fh'

      external vec_data

_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF

c Local
      integer idum, i, ivec

      if(n_vec.eq.maxvec)then
        call dlc_error('too many sets of vectors','abort')
c        allocate_vectors = -1
c        return
      endif

      do i = 1,maxvec
        if(Alph_vec(i) .eq. -1)then
           ivec = i
           n_vec = n_vec + 1
           goto 100
        endif
      enddo

      call dlc_error('internal error - failure to allocate vectors',
     &     'abort')

 100  continue

      Alph_vec(ivec) = matrix_create(nbas,nbas,'Alpha',MATRIX_REAL)
      Occa_vec(ivec) = matrix_create(nbas,1,'Alpha Occ',MATRIX_REAL)
      Eiga_vec(ivec) = matrix_create(nbas,1,'Alpha Eig',MATRIX_REAL)

      if(uhf)then
      Beta_vec(ivec) = matrix_create(nbas,nbas,'Beta',MATRIX_REAL)
      Occb_vec(ivec) = matrix_create(nbas,1,'Beta Occ',MATRIX_REAL)
      Eigb_vec(ivec) = matrix_create(nbas,1,'Beta Eig',MATRIX_REAL)
      endif

      itop_vec(ivec) = nbas
      ncol_vec(ivec) = nbas
      nrow_vec(ivec) = nbas
c
c Zero vectors of occupations and eigenvalues
c
      idum = matrix_assign_zero(Occa_vec(ivec))
      idum = matrix_assign_zero(Eiga_vec(ivec))

      if(uhf)then
         idum = matrix_assign_zero(Occb_vec(ivec))
         idum = matrix_assign_zero(Eigb_vec(ivec))
      endif

      uhf_vec(ivec) = uhf

      allocate_vectors=ivec

      end

c************************************************************************
      integer function destroy_vectors(Vectors)
c************************************************************************
      implicit none

      integer Vectors

INCLUDE(vectors_internal.fh)
INCLUDE(matrix_internal.fh)

      include 'matrix_interface.fh'
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF

c Local
      integer idum
      logical uhf

      n_vec = n_vec  - 1

      uhf = uhf_vec(Vectors)

      idum = matrix_destroy(Alph_vec(Vectors))
      idum = matrix_destroy(Occa_vec(Vectors))
      idum = matrix_destroy(Eiga_vec(Vectors))

      if(uhf)then
         idum = matrix_destroy(Beta_vec(Vectors))
         idum = matrix_destroy(Occb_vec(Vectors))
         idum = matrix_destroy(Eigb_vec(Vectors))
      endif

      itop_vec(Vectors) = 0
      ncol_vec(Vectors) = 0
      nrow_vec(Vectors) = 0

      Alph_vec(Vectors) = -1

      destroy_vectors = 0

      end


c************************************************************************
      integer function copy_vectors(From,To)
c************************************************************************
c
c  Make a copy of a set of vectors
c  (currently Eigenvectors, Occupations and Eigenvalues)
c  Only for RHF at present
c
c************************************************************************
      implicit none
      integer From, To

INCLUDE(vectors_internal.fh)
INCLUDE(matrix_internal.fh)

      include 'matrix_interface.fh'

      integer idum
c
      idum = matrix_copy(Alph_vec(From),Alph_vec(To))
      idum = matrix_copy(Occa_vec(From),Occa_vec(To))
      idum = matrix_copy(Eiga_vec(From),Eiga_vec(To))
      if(uhf_vec(From))then
         idum = matrix_copy(Beta_vec(From),Beta_vec(To))
         idum = matrix_copy(Occb_vec(From),Occb_vec(To))
         idum = matrix_copy(Eigb_vec(From),Eigb_vec(To))
      endif

      copy_vectors = 0

      end
	

c***********************************************************************
      integer function print_vectors(Vect)
c***********************************************************************
      implicit none

      integer Vect

INCLUDE(vectors_internal.fh)
INCLUDE(../m4/common/iofile)
INCLUDE(matrix_internal.fh)

      include 'matrix_interface.fh'

      integer Alpha, nbas, ncol, ivec, ieig
      logical uhf
      
      Alpha = Alph_vec(Vect)
      uhf = uhf_vec(Vect)

c Harmonic alert
      nbas = matrix_dimension(Alpha,1)
      ncol = matrix_dimension(Alpha,2)

      ivec = hp(matrix_handle(Alpha))
      ieig = hp(matrix_handle(Eiga_vec(Vect)))
      if(uhf)write(iwr,*)'Eigenvectors - Alpha set'
      if(.not.uhf)write(iwr,*)'Eigenvectors'
      call prev(QM(ivec),QM(ieig),ncol,nbas,nbas)
      if(uhf)then
         write(iwr,*)'Eigenvectors - Beta set'
         ivec = hp(matrix_handle(Beta_vec(Vect)))
         ieig = hp(matrix_handle(Eigb_vec(Vect)))
         call prev(QM(ivec),QM(ieig),ncol,nbas,nbas)
      endif
      print_vectors = 0
      end

c***********************************************************************
      logical function vec_unrestricted(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of alpha coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      vec_unrestricted = uhf_vec(Vect)
      end

c***********************************************************************
      integer function vec_alpha_coefficients(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of alpha coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      vec_alpha_coefficients = Alph_vec(Vect)
      end

c***********************************************************************
      integer function vec_alpha_occupations(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of alpha coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      vec_alpha_occupations = Occa_vec(Vect)
      end

c***********************************************************************
      integer function vec_alpha_eigenvalues(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of alpha coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      vec_alpha_eigenvalues = Eiga_vec(Vect)
      end

c***********************************************************************
      integer function vec_beta_coefficients(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of beta coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      if(.not. uhf_vec(Vect))call caserr('No beta vectors')
      vec_beta_coefficients = Beta_vec(Vect)
      end

c***********************************************************************
      integer function vec_beta_occupations(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of beta coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      if(.not. uhf_vec(Vect))call caserr('No beta vectors')
      vec_beta_occupations = Occb_vec(Vect)
      end

c***********************************************************************
      integer function vec_beta_eigenvalues(Vect)
c***********************************************************************
c  Given a vectors handle, returns the matrix handle used to store the
c  set of beta coefficients
c***********************************************************************
      implicit none
      integer Vect
INCLUDE(vectors_internal.fh)      
      if(.not. uhf_vec(Vect))call caserr('No beta vectors')
      vec_beta_eigenvalues = Eigb_vec(Vect)
      end


c***********************************************************************
      integer function orthogonalise_vectors(Vectors,
     &	Overlap,Scratch,Otran)
c***********************************************************************
c
c Orthogonalise a set of vectors
c
c   - this version uses the orfog routine
c
c   Vectors - vectors object 
c   Overlap - matrix object
c   Scratch - matrix object  (workspace)
c   Otran   - matrix object  (scratch)
c
c***********************************************************************

      implicit none
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
      integer Vectors, Overlap, Scratch, Otran

INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'
      include 'vectors_interface.fh'

c External
      integer igmem_alloc
      external igmem_alloc

c Local
      integer Vmat, l0, l1, l2
      integer ixv, io, iscr, idum
      logical uhf

      Vmat = vec_alpha_coefficients(Vectors)
      uhf = vec_unrestricted(Vectors)
c
c  harmonic alert
c
      l0 = matrix_dimension(Vmat,2) 
      l1 = matrix_dimension(Vmat,1)
      l2 = (l1+1)*l1/2
      io   = igmem_alloc(l2)
      iscr = igmem_alloc(l2)
c
c Convert overlap to eigenvector basis
c
c          call mult2(q(i30),q(i10),q(i50),l0,l0,l1)
c          call orfog(q(i30),q(i30),q(i10),q(i20),iky,
c     +    ilifq,l0,l1,1)
c
c      write(6,*)'Orthog'
c      write(6,*)'overlap before transform'
c      call matrix_print(Overlap)
c      call matrix_print(Vmat)

      idum=matrix_mult2(Overlap,Vmat,Otran,Scratch)

c      write(6,*)'overlap after transform'
c      call matrix_print(Otran)

c      subroutine orfog(q,qp,b,c,iky,ilifq,newbas,nrow,iop)
c      dimension q(*),qp(*),b(*),c(*),iky(*),ilifq(*)
c      common/blkin/p(1)
c... to orthogonalize the cols of q - result to qp
c... overlap matrix supplied in b - destroyed on exit
c... scratch triangle in c
c... iop=1 normal mode iop=2 qp set as if q was given as identity
c... qp can overwrite q

      ixv  = hp(matrix_handle(Vmat))

      idum=matrix_get_to_triangle(Otran,Q(io))
 
ccc      call prtri(Q(io),l1)
ccc      idum = matrix_print(Vmat)

c      call prsq(qm(ixv),l1,l1,l1)
   
      call orfog(QM(ixv),QM(ixv),
     &    Q(io),Q(iscr),iky,
     &    ilifq,l0,l1,1)

      if(uhf)then

      Vmat = vec_beta_coefficients(Vectors)

      idum=matrix_mult2(Overlap,Vmat,Otran,Scratch)

      ixv  = hp(matrix_handle(Vmat))

      idum=matrix_get_to_triangle(Otran,Q(io))

      call orfog(QM(ixv),QM(ixv),
     &    Q(io),Q(iscr),iky,
     &    ilifq,l0,l1,1)

      endif

      call gmem_free(iscr)
      call gmem_free(io)

      orthogonalise_vectors = 0

      end

c************************************************************************
      subroutine make_density(AlphaDensity,BetaDensity,Vectors)
c************************************************************************
c
c  Construct Density Matrix
c
c   Density    Matrix tag           Out   Result density
c   Vectors    Vectors tag          In    Input vectors
c
c************************************************************************

      implicit none

_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
INCLUDE(vectors_internal.fh)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)

      include 'vectors_interface.fh'
      include 'matrix_interface.fh'
INCLUDE(matrix_internal.fh)

      integer Vectors, AlphaDensity, BetaDensity

      integer i, j, k, n, m
      integer Occ, Vmat
      integer idum
      integer id, iv, io
      logical uhf

      uhf = vec_unrestricted(Vectors)
      n = matrix_dimension(AlphaDensity,1)
      m = matrix_dimension(AlphaDensity,1)

      Occ  = vec_alpha_occupations(Vectors)
      Vmat  = vec_alpha_coefficients(Vectors)

      id = hp(matrix_handle(AlphaDensity))
      iv = hp(matrix_handle(Vmat))
      io = hp(matrix_handle(Occ))

c      write(6,*)'input to density build'
c      call matrix_print(Vmat)
c      call matrix_print(Occ)

      call dmtx3(QM(id),QM(iv),QM(io),iky,m,n,n)

      if(uhf)then

         Occ  = vec_beta_occupations(Vectors)
         Vmat  = vec_beta_coefficients(Vectors)

         id = hp(matrix_handle(BetaDensity))
         iv = hp(matrix_handle(Vmat))
         io = hp(matrix_handle(Occ))

         call dmtx3(QM(id),QM(iv),QM(io),iky,m,n,n)

      endif

      end

      subroutine dmtx3(d,v,p,ia,m,n,ndim)
_IF1(a)cvd$r noconcur
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ia(*),d(*),v(ndim,*),p(*)
      data dzero /0.0d0/

      call vclr(d,1,n*n)
      if(m.eq.0) return
      ii=0
      do 120 i = 1,n
        do 130 k = 1,m
        dum = p(k)*v(i,k)
        if(dum.ne.dzero) then
           call daxpy(n,dum,v(1,k),1,d(ii+1),1)
        endif
  130   continue
  120 ii=ii+n

      return
      end

      subroutine assign_occupations_by_energy(Vectors,nalpha,nbeta)
      implicit none

      integer Vectors
      integer nalpha, nbeta

INCLUDE(vectors_internal.fh)
INCLUDE(matrix_internal.fh)
      include 'vectors_interface.fh'
      include 'matrix_interface.fh'

      integer Occa, Occb, iocca, ioccb
      integer Eiga, Eigb, ieiga, ieigb, l1, i
      integer na, nocmx, nocmxa, nocmxb
      logical uhf

      uhf = vec_unrestricted(Vectors)

      Occa  = vec_alpha_occupations(Vectors)
      iocca  = hp(matrix_handle(Occa))

      Eiga = vec_alpha_eigenvalues(Vectors)
      ieiga = hp(matrix_handle(Eiga))

      l1 = matrix_dimension(Occa,1)

      if(uhf)then

         Occb  = vec_beta_occupations(Vectors)
         ioccb  = hp(matrix_handle(Occb))

         Eigb = vec_beta_eigenvalues(Vectors)
         ieigb = hp(matrix_handle(Eigb))

         call llvmo(QM(ieiga),QM(iocca),
     &              nalpha,nocmxa,l1)
         call llvmo(QM(ieigb),QM(ioccb),
     &              nbeta,nocmxb,l1)

         itop_vec(Vectors) = max(nocmxa,nocmxb)

      else

         call llvmo(QM(ieiga),QM(iocca),
     &              nalpha,nocmx,l1)
         call dscal(l1,2.0d0,QM(iocca),1)

         itop_vec(Vectors) = nocmx
      endif

      end

      subroutine assign_occupations_by_overlap(Vectors,OldVectors,
     &     nalpha, nbeta, check)

      implicit none

      integer Vectors, OldVectors
      integer nalpha, nbeta
      logical check
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
INCLUDE(vectors_internal.fh)

      include 'vectors_interface.fh'
INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'

      integer Vnew, Vold, Occnew, Occold, Eig
      integer l0, l1, l3
      integer ioldvi, iscr, iscr1, iscr2
      integer iocc, ioccold, iv, ivold, ifpvi, ieig
      integer nocmx, len

      logical checkb
      logical uhf

c External
      integer lenwrd
      external lenwrd

      integer igmem_alloc
      external igmem_alloc

      uhf = vec_unrestricted(Vectors)

      Vnew = vec_alpha_coefficients(Vectors)
      iv = hp(matrix_handle(Vnew))

      Vold = vec_alpha_coefficients(OldVectors)
      ivold = hp(matrix_handle(Vold))

      Occnew = vec_alpha_occupations(Vectors)
      iocc = hp(matrix_handle(Occnew))

      Occold = vec_alpha_occupations(OldVectors)
      ioccold = hp(matrix_handle(Occold))

      Eig = vec_alpha_eigenvalues(Vectors)
      ieig = hp(matrix_handle(Eig))
c
c harmonic alert
c
      l0 = matrix_dimension(Vnew,1)
      l1 = matrix_dimension(Vnew,2)
      len = (l1-1)/lenwrd() + 1
      l3 = l1*l1

      iscr = igmem_alloc(l3)
      ifpvi = igmem_alloc(l3)
_IF(unicos)
      iscr1 = igmem_alloc(len+len)
_ELSE
      iscr1 = igmem_alloc(len)
_ENDIF
      iscr2 = igmem_alloc(len)

      call assign0(l0,l1,QM(iv),QM(ivold),
     &     QM(iocc),QM(ioccold),
     &     nocmx,nalpha,check,
     &     QM(ieig),
     &     Q(ifpvi),Q(iscr),
     &     Q(iscr1),Q(iscr2))

      if(uhf)then

         Vnew = vec_beta_coefficients(Vectors)
         iv = hp(matrix_handle(Vnew))

         Vold = vec_beta_coefficients(OldVectors)
         ivold = hp(matrix_handle(Vold))

         Occnew = vec_beta_occupations(Vectors)
         iocc = hp(matrix_handle(Occnew))

         Occold = vec_beta_occupations(OldVectors)
         ioccold = hp(matrix_handle(Occold))

         Eig = vec_beta_eigenvalues(Vectors)
         ieig = hp(matrix_handle(Eig))

         call assign0(l0,l1,QM(iv),QM(ivold),
     &        QM(iocc),QM(ioccold),
     &        nocmx,nbeta,checkb,
     &        QM(ieig),
     &        Q(ifpvi),Q(iscr),
     &        Q(iscr1),Q(iscr2))

         check = check .and. checkb

      endif

      call gmem_free(iscr2)
      call gmem_free(iscr1)
      call gmem_free(ifpvi)
      call gmem_free(iscr)

      itop_vec(Vectors) = nocmx

      end

      subroutine assign0(l0,l1,v,oldv,occ,oldocc,
     &     nocmx,na,assign,
     &     e,oldvi,scr,scr1,scr2)

      implicit none

      integer l0, l1, nocmx, na
      logical assign, oprint
      REAL v(*), oldv(*), occ(*), oldocc(*)
      REAL oldvi(*), scr(*),scr1(*),scr2(*)
      REAL e(*)

INCLUDE(../m4/common/iofile)

      logical odbg, nochange

      integer i, j, lprnt, iforb
      REAL qmax, deter

      odbg = .false.
      oprint = .true.
c
c try and set occs based on orbitals of previous cycle
c
      call dcopy(l1*l1,oldv,1,oldvi,1)
      deter=-1.0d0
      call minvrt(oldvi,l1,deter,scr1,scr2)

c      if(odbg)write(6,*)'deter pv',deter

c     use orbitals of previous cycle
      call tfsqc(v,oldvi,scr,l0,l1,l1)

      lprnt = l0
      if(odbg )then
         write (iwr,91482) 
91482    format(//1x,100('-')//
     +        50x,12('-')/
     +        50x,'eigenvectors in old MO basis'/
     +        50x,12('-'))
         call prev(v,e,lprnt,l1,l1)
      endif

      if(oprint)then
         write(iwr,*)
         write(iwr,*)'Mapping to previous vectors by overlap:'
         write(iwr,*)'New MO   Old MO    Coef       Occ'
      endif

      nochange = .true.
      do i = 1, l0
         qmax = 0.0d0
         do j = 1, l1
            if( abs(v((i-1)*l1 + j)) .gt. abs(qmax))then
               qmax = v((i-1)*l1 + j)
               iforb = j
            endif
         enddo
         occ(i) = oldocc(iforb)
         if(oprint .and. 
     &        ((i .ne. iforb) .or. (abs(qmax) .lt. 0.95d0)))then
            write(iwr,100)i,iforb, qmax, occ(i)
            nochange = .false.
         endif
 100     format(2i8,f10.4,2x,f6.4)
      enddo

      if(oprint .and. nochange)
     &     write(iwr,*)' - no significant changes - '

      j = 0
      do i = 1, l0
         if(occ(i) .gt. 1.0d-10)then
            nocmx = i
            j = j + 1
         endif
      enddo

      if(j .ne. na)then
         write(iwr,*)'problem setting occupancies'
         write(iwr,*)'counts',j,na
         assign = .false.
      else
         assign = .true.
      endif
c
c Restore vectors to AO basis
c
      call tfsqc(v,oldv,scr,l0,l1,l1)

      end

      subroutine summarise_frontier_orbitals(Vectors)

      implicit none

      integer Vectors
      integer nalpha, nbeta

INCLUDE(../m4/common/iofile)
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
INCLUDE(../m4/common/sizes)
INCLUDE(vectors_internal.fh)

      include 'vectors_interface.fh'
INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'

      integer igmem_alloc
      external igmem_alloc

      integer Vmat, iv, Occ, iocc, Eig, ieig, l1
      integer lo, hi, off, ncont, iscr

      integer j
      integer ix(maxorb)

      logical uhf

      uhf = vec_unrestricted(Vectors)

      Occ  = vec_alpha_occupations(Vectors)
      iocc  = hp(matrix_handle(Occ))

      Eig = vec_alpha_eigenvalues(Vectors)
      ieig  = hp(matrix_handle(Eig))

      Vmat = vec_alpha_coefficients(Vectors)
      iv  = hp(matrix_handle(Vmat))

      l1 = matrix_dimension(Vmat,1)
      iscr = igmem_alloc(l1)

      write(iwr,*)
      if(uhf)then
         write(iwr,*)'Intermediate Frontier orbitals (Alpha set):'
      else
         write(iwr,*)'Intermediate Frontier orbitals:'
      endif
      write(iwr,*)

      lo=max(1,itop_vec(Vectors)-5)
      hi=min(l1,itop_vec(Vectors)+5)
      ncont= (min(5,l1))

      write(iwr,*)'Orbital     Energy     Occ.   Contributions'
      do off = lo,hi

         do j = 1,l1
            Q(iscr + j -1) = QM(iv + (off-1)*l1 + (j - 1))
            ix(j)=j
         enddo

         call sortit(l1,Q(iscr),ix)
         
         write(iwr,100)off,QM(ieig -1 + off),
     &        QM(iocc -1 + off),(ix(j),
     &        Q(iscr + j - 1),j=1,ncont)
      enddo

      if(uhf)then
         write(iwr,*)'Intermediate Frontier orbitals (Beta set):'

      Occ  = vec_beta_occupations(Vectors)
      iocc  = hp(matrix_handle(Occ))

      Eig = vec_beta_eigenvalues(Vectors)
      ieig  = hp(matrix_handle(Eig))

      Vmat = vec_beta_coefficients(Vectors)
      iv  = hp(matrix_handle(Vmat))

      write(iwr,*)'Orbital     Energy     Occ.   Contributions'
      do off = lo,hi

         do j = 1,l1
            Q(iscr + j -1) = QM(iv + (off-1)*l1 + (j - 1))
            ix(j)=j
         enddo

         call sortit(l1,Q(iscr),ix)
         
         write(iwr,100)off,QM(ieig -1 + off),
     &        QM(iocc -1 + off),(ix(j),
     &        Q(iscr + j - 1),j=1,ncont)
      enddo

      endif

      call gmem_free(iscr)

 100  format(1x,i7,f14.8,f8.4,5(i4,':',f10.4,3x))

      end

      SUBROUTINE SORTIT(N,D,IX)

      REAL D(N)
      INTEGER IX(N)

C...SORTS THE ELEMENTS OF D(N) INTO DESCENDING absolute ORDER AND MOVES THE
C   CORRESPONDING ELEMENTS OF IX

      REAL S
      INTEGER ITEMP

      IF (N.EQ.1) RETURN
      NM1 = N - 1
      DO 3 I = 1,NM1
         K=I
         S = D(I)
         IP1 = I + 1
         DO 1 J = IP1,N
            IF (abs(D(J)) .LE. abs(S)) GO TO 1
            K = J
            S = D(J)
1           CONTINUE
         IF (K .LE. I) GO TO 3
         D(K) = D(I)
         D(I) = S
         ITEMP = IX(I)
         IX(I) = IX(K)
         IX(K) = ITEMP
3        CONTINUE
      RETURN
      END

      block data vec_data
INCLUDE(vectors_internal.fh)    
      data n_vec/0/
      data alph_vec/maxvec*-1/
      end
