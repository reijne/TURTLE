_IF(ma)
_MACRO(QM,qq(ivoff+`$1'))
_ENDIF

c
c  Developmental SCF program, invoke by adding the
c  additional input block
c
c    newscf
c    <control directives here>
c    end
c
c conventions
c    initial capital denotes a matrix or vector pointer
c    (integer type in F77 version)
c    a.g. Alpha
c
c Warnings:
c
c No account is made of symmetry adaption at present
c   use adapt off
c
c Don't use harmonic (l0/l1 confusion still reigns)
c
c No direct mode
c
_IF(ccpdft)
      subroutine mem_newscf(l1,l2,l3,uhf,direct,imemnew)
      implicit none
c
c     This subroutine computes the maximum number of words that are
c     needed during the main loop
c
c...  Parameters
c
      integer MT_DBL
      parameter (MT_DBL=1013)
c
c...  In variables
c
      integer l1, l2, l3
      logical uhf, direct
c
c...  Out variables
c
      integer imemnew
c
c...  Local variables
c
      logical diis1, lock1
      integer n_rdmat, n_scr, n_dens, n_fock, n_dhstaru
      integer n_scr1, n_scr2, n_scr3, n_ipi, n_ipr, n_ipa
      integer n_fpvi, n_io, i
      integer imemcur, imemohd
c
c...  Functions
c
      integer  lenwrd
      external lenwrd
_IF(ma)
      integer  MA_sizeof_overhead
      external MA_sizeof_overhead
_ENDIF
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/ccpdft.hf77)
c
      imemcur = 0
      imemnew = 0
_IF(ma)
      imemohd = MA_sizeof_overhead(MT_DBL)
_ELSEIF(dynamic_memory)
      imemohd = 0
_ELSE
      imemohd = 4
_ENDIF
c
c     call fock_build
c
         if (direct) then
            n_rdmat = 2*ikyp(nshell)+imemohd
         else
            n_rdmat = 0
         endif
         imemcur = imemcur + n_rdmat
         imemnew = max(imemnew,imemcur)
c
c        call make_rdmat
c
            if (uhf) then
               n_scr = 2*l2+imemohd
            else
               n_scr = l2+imemohd
            endif
            imemcur = imemcur + n_scr
            imemnew = max(imemnew,imemcur)
            imemcur = imemcur - n_scr
c
c        leave make_rdmat
c
         if (uhf) then
            n_dens = 2*l2+imemohd
            n_fock = 2*l2+imemohd
         else
            n_dens = l2+imemohd
            n_fock = l2+imemohd
         endif
         n_scr = l2+imemohd
         if (CD_2e().and.(.not.CD_HF_exchange())
     +              .and.(.not.CD_HF_coulomb())) then
            n_dhstaru = 0
         else
            call mem_dhstaru(n_dhstaru)
         endif
         imemcur = imemcur + n_dens
         imemcur = imemcur + n_fock
         imemcur = imemcur + n_scr
         imemcur = imemcur + n_dhstaru
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_dhstaru
         imemcur = imemcur - n_scr
         imemcur = imemcur - n_fock
         imemcur = imemcur - n_dens
         imemcur = imemcur - n_rdmat
c
c     leave fock_build
c
      diis1 = .false.
      lock1 = .false.
      do i = 1, nphase
         diis1 = diis1.or.diis(i)
         lock1 = lock1.or.lock_vec(i)
      enddo
c
c     call solve_diis
c
      if (diis1) then
         n_scr1 = l2 + imemohd
         n_scr2 = l2 + imemohd
         n_scr3 = l2 + imemohd
         imemcur = imemcur + n_scr1
         imemcur = imemcur + n_scr2
         imemcur = imemcur + n_scr3
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_scr3
         imemcur = imemcur - n_scr2
         imemcur = imemcur - n_scr1
      endif
c
c     leave solve_diis
c
c     call matrix_diagonalise
c
         n_ipi = l1 + imemohd
         n_ipr = 8*l1 + imemohd
         n_ipa = l1*l1 + imemohd
         imemcur = imemcur + n_ipi
         imemcur = imemcur + n_ipr
         imemcur = imemcur + n_ipa
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_ipa
         imemcur = imemcur - n_ipr
         imemcur = imemcur - n_ipi
c
c     leave matrix_diagonalise
c
c     call assign_occupations_by_overlap
c
      if (lock1) then
         n_scr   = l3+imemohd
         n_fpvi  = l3+imemohd
         n_scr1  = (l1+lenwrd()-1)/lenwrd() + imemohd
         n_scr2  = (l1+lenwrd()-1)/lenwrd() + imemohd
         imemcur = imemcur + n_scr
         imemcur = imemcur + n_fpvi
         imemcur = imemcur + n_scr1
         imemcur = imemcur + n_scr2
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_scr2
         imemcur = imemcur - n_scr1
         imemcur = imemcur - n_fpvi
         imemcur = imemcur - n_scr
      endif
c
c     leave assign_occupations_by_overlap
c
c     call orthogonalise_vectors
c
         n_io  = l2 + imemohd
         n_scr = l2 + imemohd
         imemcur = imemcur + n_io
         imemcur = imemcur + n_scr
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_scr
         imemcur = imemcur - n_io
c
c     leave orthogonalise_vectors
c
c     call summarise_frontier_orbitals
c
         n_scr = l1 + imemohd
         imemcur = imemcur + n_scr
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_scr
c
c     leave summarise_frontier_orbitals
c
c     call save_orbitals
c
         n_dens = l2 + imemohd
         imemcur = imemcur + n_dens
         imemnew = max(imemnew,imemcur)
         imemcur = imemcur - n_dens
c
c     leave save_orbitals
c
      if (imemcur.ne.0) then
         call caserr('mem_newscf: messed up!!!')
      endif
      end
_ENDIF
c
c
c
      subroutine newscf(core,uhf,direct,
     &     enuc,ehf,etot,sz,s2,ek,vir)
      implicit none

      REAL core(*)
      logical uhf, direct
      REAL enuc,ehf,etot,sz,s2,ek,vir

INCLUDE(../m4/common/sizes)
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/errcodes)

INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/dump3)
INCLUDE(../m4/common/scra7)
INCLUDE(../m4/common/ccpdft.hf77)

c External
      logical opg_root
      integer igmem_alloc
      include 'vectors_interface.fh'
INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'

      REAL check_ov

c Local
      logical dft
      integer l0, l1, l2, l3
      REAL ehf1, ehf2, edft, ehf0, elow
      REAL energy_change
      integer iscr, ieig, iocc, ix

c - vector handles
      integer Vectors, OldVectors, BestVectors, TmpVectors

c - matrix handles
      integer Alpha, Beta
      integer St, Scratch, Scratch2, Scratch3
      integer FockM2, FockM0, FockM1
      integer FockM2Beta, FockM0Beta, FockM1Beta
      integer AlphaTFock, BetaTFock
      integer AlphaTrVec, BetaTrVec
      integer AlphaTVmat, BetaTVMat
      integer Svec, Sval
      integer Occa, Eiga, Occb, Eigb
      integer AlphaDensity, BetaDensity
      integer Overlap, HCore, AlphaFock, BetaFock

      integer niter, niter_phase, maxcyc
      integer idum
      integer iphase, newphase
      logical converged, check, o1(maxphase), o2(maxphase), 
     &     o3(maxphase), o4(maxphase)
      REAL shifta, shiftb, tester
      integer nquad
      REAL alpha_quad, beta_quad
      REAL len1, len2, coef1, coef2
      REAL ener0, ener1, ener2
      integer num_extrap_cyc

      REAL extrap_len1
      REAL extrap_len2_o
      REAL extrap_len2_n
      REAL extrap_fac
      REAL fac2
      logical use_extrap, use_diis

      logical lock1 
      logical diis1, diis_on
      logical extrap1, extrap_on
      REAL diis_tester, diis_error
      REAL maxden, rmsden, maxfok,rmsfok, fac
      REAL tcyc_cpu, ttot_cpu, tcyc_wall, ttot_wall
      REAL t0_cpu(3), t0_wall, t1_cpu(3), t1_wall
      REAL t2_cpu(3), t2_wall
      integer diis_dimension
      character *1 xn
      integer i
      logical restore

      integer  imemdhstaru, imemspare, imemreq, imemfree, ierror
      integer  igmem_max_memory
      external igmem_max_memory

      integer check_conv

      character flags*5

      logical ofirst
      data ofirst/.true./

      data xn/'n'/

      save ofirst

      dft = CD_active()

_IF(ccpdft)
      if(dft)then
         idum = CD_update_geom(c)
      endif
_ENDIF

      l1 = num
      l2 = l1*(l1+1)/2
      l3 = l1*l1

      sz = 0.0d0
      s2 = 0.0d0
      ek = 0.0d0
      vir = 0.0d0

      num_extrap_cyc = 0

c
c     check extrapolation and DIIS request
c   
      use_extrap=.false.
      use_diis  =.false.
      do iphase = 0, maxphase
         if(extrap(iphase)) use_extrap = .true.
         if(diis(iphase))   use_diis = .true.
      enddo


      call gms_cputime(t0_cpu)
      call walltime(t0_wall)

      call gms_cputime(t1_cpu)
      call walltime(t1_wall)
c
c Set default convergence approach, and print out
c 
      call default_conv(uhf)
      call print_conv(uhf)
c
c Initialise heap storage
c
      if(ofirst)call memory_create
      ofirst = .false.

c
c  Allocate Fock/Error matrices for DIIS
c
      call diis_initialise(l1,uhf)

      enuc = nuclear_energy()
c
c Load guess vectors
c      
      Vectors = allocate_vectors(num,uhf)
      Alpha = vec_alpha_coefficients(Vectors)
      if(uhf)Beta = vec_beta_coefficients(Vectors)
c
c previous set (used for lock)
c best set (can be replaced with I/O later)
c
      OldVectors = allocate_vectors(num,uhf)
      BestVectors = allocate_vectors(num,uhf)
      TmpVectors = allocate_vectors(num,uhf)

      iscr = igmem_alloc(l3)
      call rdedx(Q(iscr),l3,ibl3qa,idaf)
      idum=matrix_set(Alpha,Q(iscr))

      if(uhf)then
         call rdedx(Q(iscr),l3,ibl3qb,idaf)
         idum=matrix_set(Beta,Q(iscr))
      endif

      if(oscfprint(PR_GUESS))then
         idum=matrix_print_titled(Alpha,'Guess eigenvectors')
         if(uhf)idum=matrix_print_titled(Beta,'Guess eigenvectors')
      endif
c
c Transformed overlap
c
      call rdedx(Q(iscr),l2,ibl7st,num8)
      St = matrix_create(num,num,'St',MATRIX_REAL)
      idum=matrix_set_from_triangle(St,Q(iscr))
c
c Overlap eigenvectors
c (currently unused, in rhfclm etc they are used 
c  as the orthogonalising transformation)
c
      Sval = matrix_create(num,1,'S eigenvals',MATRIX_REAL)
      Svec = matrix_create(num,num,'S vectors',MATRIX_REAL)
      idum =  matrix_diagonalise(St, SVec, Sval, l1, .true.)
c
c AO Overlap matrix (used by diis solver)
c
      call rdedx(Q(iscr),l2,ibl7s,num8)
      Overlap = matrix_create(num,num,'Overlap',MATRIX_REAL)
      idum=matrix_set_from_triangle(Overlap,Q(iscr))
c
c Core Hamiltonian
c
      call rdedx(Q(iscr),l2,ibl7f,num8)
      Hcore = matrix_create(num,num,'HCore',MATRIX_REAL)
      idum=matrix_set_from_triangle(HCore,Q(iscr))

      if(oscfprint(PR_FOCK))then
         idum=matrix_print(Hcore)
      endif

      call gmem_free(iscr)
c
c Scratch matrices
c
      Scratch = matrix_create(num,num,'Scratch',MATRIX_REAL)
      Scratch2 = matrix_create(num,num,'Scratch2',MATRIX_REAL)
      Scratch3 = matrix_create(num,num,'Scratch3',MATRIX_REAL)
c
c allocate Fock, transformed fock, transformed vector
c  and Density Matrices
c
      AlphaFock = matrix_create(l1,l1,'Fock',MATRIX_REAL)
      AlphaTFock = matrix_create(l1,l1,'Transformed Fock',MATRIX_REAL)
      AlphaTVmat = matrix_create(l1,l1,'Tr Fock Vecs',MATRIX_REAL)
      AlphaDensity = matrix_create(num,num,'Alpha Density',MATRIX_REAL)

      if(uhf)then
         BetaFock = matrix_create(l1,l1,'Fock',MATRIX_REAL)
         BetaTFock = matrix_create(l1,l1,'Transformed Fock',MATRIX_REAL)
         BetaTVmat = matrix_create(l1,l1,'Tr Fock Vecs',MATRIX_REAL)
         BetaDensity = matrix_create(num,num,'Beta Density',MATRIX_REAL)
      endif

c
c these are zero-ed so that the the initial differences
c can be computed without floating pt errors
c
      idum=matrix_assign_zero(AlphaFock)
      if(uhf)idum=matrix_assign_zero(BetaFock)

c
c These are used for extrapolation
c
      if(use_extrap)then
         FockM0 = matrix_create(l1,l1,'Fock 0',MATRIX_REAL)
         FockM1 = matrix_create(l1,l1,'Fock -1',MATRIX_REAL)
         FockM2 = matrix_create(l1,l1,'Fock -2',MATRIX_REAL)
         if(uhf)then
            FockM0Beta = matrix_create(l1,l1,'Fock 0',MATRIX_REAL)
            FockM1Beta = matrix_create(l1,l1,'Fock -1',MATRIX_REAL)
            FockM2Beta = matrix_create(l1,l1,'Fock -2',MATRIX_REAL)
         endif
      endif
c
c Orthog guess vectors
c
      call orthogonalise_vectors(Vectors,St,Scratch,Scratch2)

cc      call wrt3(q(i30),l3,ibl3qa,idaf)
cc      call tdown(q(i30),ilifq,q(i30),ilifq,l1)

c
c Load guess vector eigenvalues
c Maybe a problem for cases where input occupancies are non-aufbau
c   (perhaps read pops as well)
c
      Eiga = vec_alpha_eigenvalues(Vectors)
      ieig  = hp(matrix_handle(Eiga))
      call rdedx(QM(ieig),l1,ibl3ea,idaf)

      if(uhf)then
         Eigb = vec_beta_eigenvalues(Vectors)
         ieig  = hp(matrix_handle(Eigb))
         call rdedx(QM(ieig),l1,ibl3eb,idaf)
      endif
      
      call assign_occupations_by_energy(Vectors,na,nb)
c
c  Diagnostic print of frontier orbitals
c
      if(oscfprint(PR_FRONTIER))then
         write(iwr,*)'Frontier orbitals of Guess vectors'
         call summarise_frontier_orbitals(Vectors)
      endif

      call make_density(AlphaDensity,BetaDensity,Vectors)

      if(oscfprint(PR_DENSITY))then
         idum=matrix_print_titled(AlphaDensity,'Alpha guess density')
         if(uhf)
     &        idum=matrix_print_titled(BetaDensity,'Beta guess density')
      endif

c      call wrt3(q(i50),l2,ibl3pa,idaf)
c store copy of density
c      call dcopy(l2,q(i50),1,q(jblkpa),1)

c
c  ----==========  Main SCF Loop =========-----
c

      niter = 0
      niter_phase = 0
      maxcyc = maxcycp
      ehf0 = 0.0d0
      converged = .false.
      iphase = 1
      elow = 1.0d20
      coef1 = 0.0d0

      if(.not. oscfprint(PR_FULL) .and. opg_root())then
         write(iwr,100)
 100     format(1x,'niter',5x,'Energy',4x,'Energy Change',4x,'Tester'
     &        )
      endif
_IF(ccpdft)
      if (CD_active()) then
         if (uhf) then
            idum = CD_uks()
         else
            idum = CD_rks()
         endif
         call retrieve_spare(imemspare)
         idum = CD_set_2e()
         call mem_newscf(l1,l2,l3,uhf,direct,imemdhstaru)
         idum = CD_reset_2e()
         imemfree = igmem_max_memory() - imemspare - imemdhstaru
         imemreq = CD_memreq_energy_ao(q0,q0,iwr)
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            write(iwr,600)ierror
            call caserr('Not enough memory for incore coulomb fit')
         endif
 600     format('*** Need ',i10,' more words to store any ',
     +          '3-center integrals in core')
      endif
_ENDIF
      do while (.not. converged .and. niter .le. maxcyc 
     &     .and. iphase .ne. -1 )



         niter = niter + 1
         niter_phase = niter_phase + 1
c
c Flags for output
c
         flags=" "

c     
c     save previous total SCF energy
c
         ehf0 = ehf
c
c set convergence control parameters depending on which
c convergence phase is active
c
         shifta = shift(iphase,1)
         if(uhf)then
           shiftb = shift(iphase,2)
         else
           shiftb = 0.0d0
         endif
         diis1 = diis(iphase)
         extrap1 = extrap(iphase)
         lock1 = lock_vec(iphase)

         idum = matrix_copy(AlphaFock,Scratch)
         if(uhf)idum = matrix_copy(BetaFock,Scratch3)

         call copy_vectors(Vectors,TmpVectors)

         call fock_build(AlphaDensity,BetaDensity,
     &        direct,
     &        AlphaFock,BetaFock,HCore,uhf,
     &        ehf, ehf1, ehf2, edft, 
     &        nquad, alpha_quad,beta_quad,
     &        core)

c
c  Check difference fock matrices
c
         idum=matrix_copy(AlphaFock,Scratch2)
         idum=matrix_daxpy(-1.0d0,Scratch,Scratch2)
         maxfok = matrix_absmax(Scratch2)
         fac = l3
         rmsfok = sqrt(matrix_dot_product(Scratch2,Scratch2)/fac)

         if(uhf)then
            idum=matrix_copy(BetaFock,Scratch2)
            idum=matrix_daxpy(-1.0d0,Scratch3,Scratch2)
            maxfok = matrix_absmax(Scratch2)
            fac = l3
            rmsfok = 0.5d0*(rmsfok + 
     &      sqrt(matrix_dot_product(Scratch2,Scratch2)/fac))
         endif

         if(oscfprint(PR_FOCK))then
            idum=matrix_print(AlphaFock)
            if(uhf)idum=matrix_print(BetaFock)
         endif

         energy_change = ehf-ehf0

         if(use_diis)then

            call solve_diis(AlphaFock,BetaFock,Vectors,Overlap,
     &           Scratch, Scratch2,
     &           diis1,diis_on,diis_dimension,tester,diis_error)

            if(diis_on.and.diis1)flags(1:1) = 'd'

         else
            diis_on = .false.
         endif

         extrap_on = .false.

         if(extrap1)then

c
c Compute unextrapolated tester
c
            AlphaTrVec = Alpha
            idum=matrix_mult2(AlphaFock,AlphaTrVec,
     &           AlphaTFock,Scratch)
            if(uhf)then
               BetaTrVec = Beta
               idum=matrix_mult2(BetaFock,BetaTrVec,
     &              BetaTFock,Scratch)
            endif

            tester = check_ov(AlphaTfock,BetaTfock,Vectors)
         endif
c
c Extrapolation phase
c
c set up two old fock matrices
c
         if (use_extrap) then
            if(num_extrap_cyc .gt. 1)then
               idum=matrix_copy(FockM1,FockM2)
               if(uhf)idum=matrix_copy(FockM1Beta,FockM2Beta)
               ener2 = ener1
            endif

            if(num_extrap_cyc .gt. 0)then
               idum=matrix_copy(FockM0,FockM1)
               if(uhf)idum=matrix_copy(FockM0Beta,FockM1Beta)
               ener1 = ener0
            endif

            idum=matrix_copy(AlphaFock,FockM0)
            if(uhf)idum=matrix_copy(BetaFock,FockM0Beta)
            ener0 = ehf
         endif

         if (extrap1) then
            if(num_extrap_cyc.eq.2)then
c
c Find 2 displacements
c (note the test is only done for alpha density)
c 
               idum=matrix_copy(FockM1,Scratch)
               idum=matrix_daxpy(-1.0d0,FockM2,Scratch)

               idum=matrix_copy(FockM0,Scratch2)
               idum=matrix_daxpy(-1.0d0,FockM1,Scratch2)
c
c Now check the angle between the two steps
c
               len1 = sqrt(matrix_dot_product(Scratch,Scratch))
               len2 = sqrt(matrix_dot_product(Scratch2,Scratch2))
               fac  = matrix_dot_product(Scratch,Scratch2) /
     &              (len1*len2)

               extrap_fac = fac
               extrap_len1 = len1
               extrap_len2_o = len2

               fac2 = len2 / len1

c              write(6,*)'Extrapolation fac,fac2,len1,len2,coef = ',
c    &              fac,fac2,len1,len2,coef1
               
               if(fac .gt. extrap_tol(iphase)  .and.
     &          (fac2 .gt. 0.7) .and. (fac2 .lt. (1.2 + coef1)) ) then
c
                  coef1 = extrap_coef(iphase)

                  idum=matrix_daxpy(coef1,Scratch2,AlphaFock)
                  if(uhf)then
                     idum=matrix_copy(FockM0Beta,Scratch2)
                     idum=matrix_daxpy(-1.0d0,FockM1Beta,Scratch2)
                     idum=matrix_daxpy(coef1,Scratch2,BetaFock)
                  endif
c
c  Compute fac again as a check
c
                  idum=matrix_copy(AlphaFock,Scratch2)
                  idum=matrix_daxpy(-1.0d0,FockM1,Scratch2)
                  len2 = sqrt(matrix_dot_product(Scratch2,Scratch2))
                  fac  = matrix_dot_product(Scratch,Scratch2) /
     &                 (len1*len2)
c
c printing diagnostics
c
                  flags(2:2) = 'E'
                  extrap_len2_n = len2
                  extrap_on = .true.

               endif
            endif
         endif
         if (use_extrap) then
            num_extrap_cyc = min(2, num_extrap_cyc + 1)
         endif

c
c Used to estimate valid step length on next cycle
c
         if(.not. extrap_on)coef1 = 0.0d0


c     if(oshift) then
c     call tdown(q(i30),ilifq,q(jblkqa),ilifq,l0)
c     else
c     call tdown(q(i30),ilifq,q(jblkqs),ilifq,l0)
c     endif

c     
c     Tranform Fock to orthonormal basis
c     (use of previous vectors is hardwired)
c
         AlphaTrVec = Alpha
         idum=matrix_mult2(AlphaFock,AlphaTrVec,
     &        AlphaTFock,Scratch)
         if(uhf)then
            BetaTrVec = Beta
            idum=matrix_mult2(BetaFock,BetaTrVec,
     &           BetaTFock,Scratch)
         endif

         if(extrap_on)then
c - skip tester (see above)
         elseif(diis_on)then
c
c Current Fock is extrapolated - can't really trust 
c the magnitude of the ov block, print out in case it is 
c a useful diagnostic for DIIS operation
c
            diis_tester = check_ov(AlphaTfock,BetaTFock,Vectors)
         else 
            tester = check_ov(AlphaTfock,BetaTfock,Vectors)
         endif
         
c      ocvged = (damp.lt.dmplim) .and. (diff.lt.acurcy) .and. (iter.gt.1)

         call level_shift(AlphaTFock,BetaTFock,Vectors,shifta,shiftb)
         if(shifta .gt. 0 .or. shiftb .gt. 0)then
            flags(3:3) = 'S'
         endif

         if(oscfprint(PR_FOCK))then
            idum=matrix_print_titled(AlphaTFock,
     &           'Transformed fock after Level shift')
            if(uhf)then
            idum=matrix_print_titled(BetaTFock,
     &              'Beta Transformed fock after Level shift')
            endif
         endif

         call copy_vectors(Vectors,OldVectors)
c
c  Diagonalise and Back-transform the vectors 
c
         idum=matrix_diagonalise(AlphaTFock,AlphaTVmat,Eiga,l1,.true.)
         idum=matrix_dgemm(xn,xn,1.0d0,AlphaTrVec,
     &        AlphaTVmat,0.0d0,Scratch)
         idum=matrix_copy(Scratch,Alpha)

         if(uhf)then
            idum=matrix_diagonalise(BetaTFock,BetaTVmat,Eigb,l1,.true.)
            idum=matrix_dgemm(xn,xn,1.0d0,BetaTrVec,
     &           BetaTVmat,0.0d0,Scratch)
            idum=matrix_copy(Scratch,Beta)
         endif

         if(oscfprint(PR_VECTORS))then
            call print_vectors(Vectors)
         endif
c
         if(lock1)then
c            call overlap_check(Vectors, OldVectors)
            call assign_occupations_by_overlap(Vectors,OldVectors,
     &           na,nb,check)
            if(.not.check)then
               call assign_occupations_by_energy(Vectors,na,nb) 
           endif

           flags(4:4) = 'L'
         else
            call assign_occupations_by_energy(Vectors,na,nb)
         endif
c
c Re-orthog vectors
c
      call orthogonalise_vectors(Vectors,St,Scratch,Scratch2)
c
c  Diagnostic print of frontier orbitals
c
         if(oscfprint(PR_FRONTIER))then
            call summarise_frontier_orbitals(Vectors)
         endif
c 
c  Construct new density matrix
c         
         call make_density(Scratch,Scratch3,Vectors)
c
c  Check difference density
c
         idum=matrix_copy(Scratch,Scratch2)
         idum=matrix_daxpy(-1.0d0,AlphaDensity,Scratch2)
         maxden = matrix_absmax(Scratch2)
         fac = l3
         rmsden = sqrt(matrix_dot_product(Scratch2,Scratch2)/fac)
         idum=matrix_copy(Scratch,AlphaDensity)

         if(uhf)then
            idum=matrix_copy(Scratch3,Scratch2)
            idum=matrix_daxpy(-1.0d0,BetaDensity,Scratch2)
            maxden = matrix_absmax(Scratch2)
            fac = l3
            rmsden = rmsden + 
     &           sqrt(matrix_dot_product(Scratch2,Scratch2)/fac)
            idum=matrix_copy(Scratch3,BetaDensity)
         endif

         if(oscfprint(PR_DENSITY))then
            idum=matrix_print(AlphaDensity)
            if(uhf)idum=matrix_print(BetaDensity)
         endif
c
c  Check for SCF convergence, or switch of convergence 
c  parameters
c
         newphase = check_conv(uhf,iphase,tester,energy_change,
     &        niter_phase,o1,o2,o3,o4)

c
c Print intermediate results
c
         etot = ehf+enuc

         call gms_cputime(t2_cpu)
         call walltime(t2_wall)

         tcyc_cpu  = t2_cpu(3) - t1_cpu(3)
         tcyc_wall = t2_wall   - t1_wall
         ttot_cpu  = t2_cpu(3) - t0_cpu(3)
         ttot_wall = t2_wall   - t0_wall

         call gms_cputime(t1_cpu)
         call walltime(t1_wall)

         if (opg_root())then
         if(oscfprint(PR_FULL))then

            write(iwr,101)niter, etot, energy_change, tester,
     &           iphase, flags,ehf1,ehf2,edft
            if(dft)then
               if(uhf)then
                  write(iwr,10111)nquad,alpha_quad,beta_quad
               else
                  write(iwr,1011)nquad,alpha_quad
               endif
            endif
            if(uhf)then
               if(shifta .ne. 0.0d0 .or. shiftb .ne. 0.0d0)
     &              write(iwr,1012)shifta,shiftb
            else
               if(shifta .ne. 0.0d0)write(iwr,1012)shifta
            endif
            if(diis_on)write(iwr,1013)diis_dimension,diis_error,
     &           diis_tester

            if(extrap1 .and. (num_extrap_cyc.eq.2) )then
               if(extrap_on)then
                  write(iwr,10132)extrap_fac,coef1,
     &                 extrap_len1,
     &                 extrap_len2_o, extrap_len2_n
               else
                  write(iwr,10133)extrap_fac,extrap_len1,extrap_len2_o
               endif
            endif

            write(iwr,10131)maxden, rmsden, maxfok,rmsfok
            write(iwr,1014)iphase,niter_phase
            do i = 1,nnext(iphase)
               if(nextphase(iphase,i).eq.0)then
                  write(iwr,1015)o1(i),o2(i),o3(i),o4(i)
               else
                  write(iwr,1016)nextphase(iphase,i),o1(i),o2(i),
     &                 o3(i),o4(i)
               endif
            enddo
            if(oscfprint(PR_TIMINGS))then
               write(iwr,1017)tcyc_cpu, ttot_cpu,tcyc_wall, ttot_wall
            endif
         else
            if(oscfprint(PR_TIMINGS))then
               write(iwr,102)niter, ehf+enuc, energy_change, tester,
     &              ttot_cpu
            else
               write(iwr,102)niter, ehf+enuc, energy_change, tester
            endif
         endif
         endif

 101     format(/1x,'Cycle ',i4,1x,
     &        'E        = ',f14.8,1x,
     &        'dE       = ',g12.6,3x,
     &        'Tester   = ',f14.8,3x,i4,1x,a5,/,12x,
     &        'EHF(1e)  = ',f14.8,1x,
     &        'EHF(2e)  = ',f14.8,1x,
     &        'EDFT     = ',f14.8)
 1011    format(12x,
     &        'Nquad    = ',i14,1x,
     &        'Nelec    = ',f14.8)
10111    format(12x,
     &        'Nquad    = ',i14,1x,
     &        'Nelec alpha = ',f14.8,' beta = ',f14.8)
 1012    format(12x,'Shift   = ',f10.4,2x,f10.4)
 1013    format(12x,'Diis On Dim = ',i3,' Error = ',f15.9,
     &        ' Extrapolated tester = ',f14.8)

10132    format(12x,'Extrapolation  test=',f14.8,' scal ',f12.4,
     &        ' len1',f12.4,
     &        ' len2(old) = ',f12.4,' len2(ext) ',f12.4)
10133    format(12x,'Extrapolation  test=',f14.8,' suppressed',10x,
     &        ' len1 ',f12.4,
     &        ' len2 ',f12.4)

10131    format(12x,'Changes - Density Max = ',f14.8,' RMS = ',f14.8,
     &        ' Fock Max = ',f14.8,' RMS = ',f14.8)
 1014    format(12x,'Convergence phase = ',i3,' step ',i3)
 1015    format(12x,'Criteria for convergence: tester: ',l1,', dE: ',
     &        l1,', |dE|: ',l1,', Ncyc: ',l1)
 1016    format(12x,'Criteria for phase', i2,': tester: ',l1,', dE: ',
     &        l1,', |dE|: ',l1,', Ncyc: ',l1)

 1017    format(12x,'CPU Times: cycle ',f10.2,' total ',f10.2,
     &        ' Wall times: cycle ',f10.2,' total ',f10.2)
 102     format(1x,i4,2x,f14.8,2x,g12.6,2x,f12.6,f10.2)
 
c
c Need to save orbitals that went into this fock build
c Not the ones that resulted from Extrap/Diag
c Save them before the restore, perhaps???
c
         if(iphase .ne. -1 .and. ehf .lt. elow) then
            if(oscfprint(PR_FULL) .and. opg_root())
     &           write(iwr,*)'Saving orbitals'
            call save_orbitals(TmpVectors,AlphaDensity,BetaDensity)
            call copy_vectors(TmpVectors,BestVectors)
            elow = ehf
         endif

         if(iphase .ne. newphase)then
            write(iwr,*)
            if(newphase.eq.-1)then
               write(iwr,*)'**** Calculation aborted '
            else if(newphase.ne.0) then
               write(iwr,*)'**** Switch to convergence phase ',
     &              newphase
               write(iwr,*)
            endif

            niter_phase = 0
            iphase = newphase
            converged = (iphase .eq. 0)
c
c restore best set of vectors if requested
c
            if(iphase .gt. 0 .and. restore_vec(iphase))then
               if (ehf .gt. elow) then
                  write(iwr,*)'Restoring orbitals'
                  call copy_vectors(BestVectors,Vectors)
                  call make_density(AlphaDensity,BetaDensity,Vectors)
               else
                  write(iwr,*)'Restoring orbitals implicitly'
               endif
            endif
c
c reset DIIS matrices if requested
c

            if(iphase .gt. 0 .and. new_diis(iphase))then
               write(iwr,*)'Resetting DIIS'
               call reset_diis	
            endif
c
c reset extrapolation matrices if requested
c
            if(iphase .gt. 0 .and. new_extrap(iphase))then
               write(iwr,*)'Resetting Extrap'
               num_extrap_cyc = 0
            endif

         endif

      enddo
_IF(ccpdft)
      if (CD_active()) then
         ierror = CD_jfit_clean2()
         if (ierror.ne.0) then
            call caserr('Memory failure in duhfop')
         endif
      endif
_ENDIF
      
      if(converged)then
         write(iwr,*)'*************'
         write(iwr,*)'SCF converged'
         write(iwr,*)'*************'
         write(iwr,9000)niter,ehf,edft,enuc,etot
 9000    format(/10x,14('-')/10x,'final energies after',i4,' cycles'/,
     &          10x,14('-')/
     &          10x,'hartree-fock energy         ',f18.10/
     &          10x,'exchange-correlation energy ',f18.10/
     &          10x,'nuclear repulsion energy    ',f18.10/
     &          10x,'total energy                ',f18.10)
         call save_orbitals(Vectors,AlphaDensity,BetaDensity)
      else
         write(iwr,*)'********************'
         write(iwr,*)'SCF did not converge'
         write(iwr,*)'********************'
         write(iwr,9000)niter,ehf,edft,enuc,0.0d0

         if(hardfail)then
            call gamerr(
     &           'SCF convergence failure',
     &           ERR_NO_CODE, ERR_UNLUCKY, ERR_SYNC, ERR_NO_SYS)
         endif
      endif

      idum = destroy_vectors(Vectors)
      idum = destroy_vectors(BestVectors)
      idum = destroy_vectors(OldVectors)
      idum = destroy_vectors(TmpVectors)

      idum = matrix_destroy(AlphaDensity)
      idum = matrix_destroy(AlphaFock)
      idum = matrix_destroy(AlphaTFock)
      idum = matrix_destroy(AlphaTVMat)

      if(uhf)then
         idum = matrix_destroy(BetaDensity)
         idum = matrix_destroy(BetaFock)
         idum = matrix_destroy(BetaTFock)
         idum = matrix_destroy(BetaTVMat)
      endif

      idum = matrix_destroy(Scratch)
      idum = matrix_destroy(Scratch2)
      idum = matrix_destroy(Scratch3)

      idum = matrix_destroy(Hcore)
      idum = matrix_destroy(Overlap)
      idum = matrix_destroy(Sval)
      idum = matrix_destroy(Svec)
      idum = matrix_destroy(St)

      if(use_extrap)then
         idum = matrix_destroy(FockM0)
         idum = matrix_destroy(FockM1)
         idum = matrix_destroy(FockM2)
         if(uhf)then
            idum = matrix_destroy(FockM0Beta)
            idum = matrix_destroy(FockM1Beta)
            idum = matrix_destroy(FockM2Beta)
         endif
      endif

      call deallocate_diis_scratch(uhf)

c      call matrix_list

      end
c

c************************************************************************
      REAL function nuclear_energy()
c************************************************************************
c
c  Compute nuclear energy of the system
c
c************************************************************************
      implicit none
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
      REAL enucf
      external enucf
      nuclear_energy = enucf(nat,czan,c)
      end

c************************************************************************
      subroutine fock_build(AlphaDensity,BetaDensity,
     &     direct,   
     &     AlphaFock,BetaFock,HCore,uhf,
     &     ehf, ehf1,ehf2,edft,
     &     nquad, alpha_quad, beta_quad,
     &     core)
c************************************************************************
c
c  Construct Fock matrix and components of the HF/KS energy
c
c   AlphaDensity  Matrix tag           In
c   BetaDensity   Matrix tag           In
c   direct        logical              In    Direct-SCF
c   AlphaFock     Matrix tag           Out
c   BetaFock      Matrix tag           Out
c   ehf1          double precision     Out
c   ehf2          double precision     Out
c   edft          double precision     Out
c   nquad         integer              Out dft results for diagostics 
c   alpha_quad    REAL                 "
c   beta_quad     REAL                 "
c
c************************************************************************
c
      implicit none

      integer AlphaDensity,BetaDensity
      logical direct
      integer AlphaFock,BetaFock,HCore
      REAL ehf, ehf1, ehf2, edft
      integer nquad
      REAL alpha_quad, beta_quad
      logical uhf
      REAL core(*)

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/vcore)
INCLUDE(../m4/common/ccpdft.hf77)
INCLUDE(../m4/common/iofile)

c
c Control parameters direct SCF
c
c needed for correct functioning of integrals
c (o = output only)
c
c n dlnmxd
c y dlntol    (dlntol is usually set to one of the tolitr values)
c n delfac   - delta factor for tolerance (delfac) when dlnmxd is small than 
c n tolitr(3) 
c   itrtol   - iteraction count to control dlntol
c 
c
c y oimag    - analyse integral magnitudes
c o intmag   - holds integral magnitudes
c o intcut   - integral test counts (numbers of integrals deleted for
c              various reasons)
c y odscf    - run direct
c
c
c y irdmat
c y iprefa
c y iexch    - address of exchange integrals  ? workspace
c y ifock    - address for fock build 
c
c   nhamd   ???   only for GVB
c   nsheld
c
c -- not used in this version
c
c  iof171  - section 171 offsets
c  m171t .. just 171
c  odnew

c
c      real*8 dlnmxd, dlntol, tolitr, deltol, delfac
c      integer ifock, idmat, irdmat, iprefa, intcut, intmag
c      integer ibl171, lentri, itrtol, iexch, iof171, m171t
c      integer nhamd, nsheld
c      logical odscf, oimag, odnew, opref

INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/timeperiods)

c External functions

      include 'matrix_interface.fh'

      logical opg_root
      integer igmem_alloc
      REAL tracep, dum

c Local variables

      integer l1, l2, l3, l22
      logical out, outon
      integer idum, i
      integer ifock,  idmat, irdmat, iprefa
      integer idens,  iscr
      integer idensb, ifockb
      REAL etmp, dft_accu

      logical o2e, odft, ocoul, oexch
      REAL facex

      logical delta
      integer nshtri
      integer kkk

_IF(debug_S)
      integer isma, ismb
_ENDIF

      l1 = num
      l2 = l1*(l1+1)/2
      l3 = l1*l1

      l22 = 2*l2

      out = nprint .eq. 5
      outon = nprint .ne. -5

      call pg_synch(4444)

      edft = 0.0d0


      delta = .false.

      if(direct)then

         if(delta)then
c -- note delta density not implemented
         endif

         dlnmxd=0.0d0
         dlntol=tolitr(3)

c         if(opg_root())then
c	    write(iwr,*)'Direct fock build dlntol: ',dlntol
c         endif

c
c  build reduced integral prefactor matrix
c  - can probably skip this on subsequent cycles
c
         call start_time_period(TP_RDMAT)
         nshtri=ikyp(nshell)

         irdmat = igmem_alloc(2*nshtri)
         call make_rdmat(Q(irdmat),uhf,AlphaDensity,BetaDensity)

c         write(6,*)'Reduced density'
c         call prtri(Q(irdmat),nshell)

c         iprefa = igmem_alloc(nshtri)
         iprefa = irdmat + nshtri
         call rdmake(Q(iprefa))

c         write(6,*)'Reduced prefactor'
c         call prtri(Q(iprefa),nshell)

         dlnmxd=-9999999.0d0
         do  kkk=1,nshtri
            if(Q(irdmat+kkk-1).gt.dlnmxd)
     &           dlnmxd=Q(irdmat+kkk-1)
         enddo

         call end_time_period(TP_RDMAT)

      endif

_IF(ccpdft)

c
c Modify fock builder options according to
c required DFT scheme (HF exchange, fitted coulomb etc
c
      idum = CD_set_2e()
_ENDIF

c
c allocate triangles for hstar etc
c hstaru assumes matrices are contiguous
c
      if(uhf)then
         idens = igmem_alloc(l22)
         idensb = idens + l2
         ifock = igmem_alloc(l22)
         ifockb = ifock + l2
      else
         idens = igmem_alloc(l2)
         ifock = igmem_alloc(l2)
      endif

c - read in fock builders
      idmat = idens

      idum = matrix_get_to_triangle(AlphaDensity,Q(idens))
      if(uhf)then
         idum = matrix_get_to_triangle(BetaDensity,Q(idensb))
      endif

      iscr = igmem_alloc(l2)

      do i = 1,l2
         Q(ifock + i -1) = 0.0d0
         Q(iscr + i -1) = 0.0d0
      enddo
c
      if (CD_2e())then
         odft = .true.
         if(uhf)idum = CD_uks()
c
c     for simplicity, only implement full coulomb version
c     here for the moment
c
         if(CD_HF_exchange())then
            facex = CD_HF_exchange_weight()
            oexch = .true.
         else
            oexch = .false.
         endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
c
      ocoul = CD_HF_coulomb() .or. .not. odft
      o2e = .not.(odft.and..not.oexch.and..not.ocoul)
c
      if(uhf)then
         do i = 1,l2
            Q(ifockb + i -1) = 0.0d0
         enddo
c
         if(o2e)then
            if(direct)then
               irest=0
               call dhstaru(core,Q(idens),Q(idensb),
     &              Q(ifock),Q(ifockb),
     &              Q(iprefa),Q(irdmat),irest)
            else
               call hstaru(Q(idens),Q(idensb),
     &              Q(ifock),Q(ifockb),nopk,
     &              facex,ocoul,oexch)
c            else
c               call vclr(Q(ifock),1,l2)
c               call vclr(Q(ifocka),1,l2)
            endif
         endif
      else
c
c rhf case
         if(o2e)then
            if(direct)then
               irest = 0
               call dhstar(core,Q(ifock),Q(idens),
     &                     Q(iprefa),Q(irdmat),
     &                     irest)
            else
_IF(cray,convex,titan)
               call hstar(Q(idens),Q(ifock),
     &                    Q(iscr),nopk,nss)
_ELSE
               call hstar(Q(idens),Q(ifock),
     &                    Q(iscr),nopk)
_ENDIF
            endif
         endif
      endif

_IF(ccpdft)
c
c restore fock builder options
c
      idum = CD_reset_2e()
_ENDIF

C - skip Zora for now

c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at x(i10)
c     scratch area at x(i20)
c
      call symh(Q(ifock),Q(iscr),iky,0,0)
      if(uhf)call symh(Q(ifockb),Q(iscr),iky,0,0)
      call gmem_free(iscr)

      if (out) then
         write (iwr,9088)
         if(uhf)write (iwr,9108)
         call prtril(Q(ifock),l1)
         if(uhf)then
            write (iwr,9128)
            call prtril(Q(ifockb),l1)
         endif
      endif

 9088 format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23(
     +     '-'))
 9108 format(//1x,129('-')///
     *50x,'----- alpha set -----')
 9128 format(//1x,129('-')///
     *50x,'----- beta set -----')

      idum = matrix_set_from_triangle(AlphaFock,Q(ifock))
      if(uhf)
     &  idum = matrix_set_from_triangle(BetaFock,Q(ifockb))

c
c     ----- read in core hamiltonian matrix
c           and calculate hf energy -----
c
c     -h0- at x(i20)
c     - h- at x(i10)
c     - e- at x(i30)
c

c      call vadd(q(i10),1,q(jblkf),1,q(i10),1,nx)

       idum = matrix_add(AlphaFock,AlphaFock,Hcore)
       if(uhf)idum = matrix_add(BetaFock,BetaFock,Hcore)
c
c  ehf1 = TrP(P*H_1)
c
       idum = matrix_get_to_triangle(AlphaDensity,Q(idens))
       if(uhf)idum = 
     &      matrix_get_to_triangle(BetaDensity,Q(idensb))

       idum = matrix_get_to_triangle(Hcore,Q(ifock))
       ehf1 = tracep(Q(idens),Q(ifock),l1)
       if(uhf)ehf1 = ehf1 + tracep(Q(idensb),Q(ifock),l1)

_IF(ccpdft)
      if(CD_active())then
c
c ccpdft: evaluate Kohn Sham energy expression
c
         if(CD_HF_exchange() .or. CD_HF_coulomb())then
c
c Coulomb or exchange operator is in Q(ifock) augmented by H_1, 
c compute energy using HF expression (maybe without exchange)
c
            idum = matrix_get_to_triangle(AlphaFock,Q(ifock))
            ehf2 = tracep(Q(idens),Q(ifock),l1)
            if(uhf)then
               idum = matrix_get_to_triangle(BetaFock,Q(ifockb))
               ehf2 = ehf2 + 
     &              tracep(Q(idensb),Q(ifockb),l1)
            endif
            etmp = (ehf1+ehf2)*0.5d0
c
         else
c
c Coulomb operator has not yet been evaluated
c (J energy will be part of edft)
c
            etmp = ehf1
         endif

c
c Update Kohn-Sham matrix and compute fitted/integrated 
c energy terms
c
cc         if (diff.ne.0.0d0) dft_accu = min(dft_accu,abs(diff))

         idum = matrix_get_to_triangle(AlphaFock,Q(ifock))
         if(uhf)
     &      idum = matrix_get_to_triangle(BetaFock,Q(ifockb))

         dft_accu = 0.0d0

_IF(debug_S)
         isma = igmem_alloc(l2)
         ismb = igmem_alloc(l2)
         call vclr(Q(isma),1,l2)
_ENDIF
         idum = CD_energy_ao(c,Q(ifock),Q(ifockb),
     &        Q(idens),Q(idensb),
     &        edft,core,core,outon,dft_accu,iwr
_IF(debug_S)
     +        ,Q(isma)
_ENDIF
     &	)

_IF(debug_S)
         call compare_S(Q(isma),Q(ismb),l2,num)
         call gmem_free(ismb)
         call gmem_free(isma)
_ENDIF


         idum = matrix_set_from_triangle(AlphaFock,Q(ifock))
         if(uhf)
     &      idum = matrix_set_from_triangle(BetaFock,Q(ifockb))

         ehf = etmp+edft
         call CD_get_dftresults(nquad,alpha_quad,beta_quad)

      else
c
c  Hartree-Fock energy expression
c  E = 1/2 [ TrP(P*H_1) + TrP(P*(H_1 + H_2)) ]
c                ehf1            ehf2
c
         idum = matrix_get_to_triangle(AlphaFock,Q(ifock))
         ehf2 = tracep(Q(idens),Q(ifock),l1)
         if(uhf)then
            idum = 
     &        matrix_get_to_triangle(BetaFock,Q(ifockb))
            ehf2 = ehf2 + 
     &        tracep(Q(idensb),Q(ifockb),l1)
         endif
         ehf = (ehf1+ehf2)*0.5d0
      endif
_ELSE
      idum = matrix_get_to_triangle(AlphaFock,Q(ifock))
      ehf2 = tracep(Q(idens),Q(ifock),l1)
      if(uhf)then
         idum = 
     &     matrix_get_to_triangle(BetaFock,Q(ifockb))
         ehf2 = ehf2 + 
     &        tracep(Q(idensb),Q(ifockb),l1)
      endif
      ehf = (ehf1+ehf2)*0.5d0
_ENDIF

      call gmem_free(ifock)
      call gmem_free(idens)

      if(direct)then
         call gmem_free(irdmat)
      endif

c - NB skip code for crystal field
      end

c***********************************************************************

      subroutine reset_diis
c***********************************************************************
c
c  Reset DIIS parameters
c
c***********************************************************************
      implicit none

      integer ndim
      logical beta

INCLUDE(diis.fh)
INCLUDE(../m4/common/iofile)

      REAL derror, scale, st, ct, r
      integer iposit, nstore, mp
      logical ondiis
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis
c
c  Reset control parameters
c
      nstore = 0
      ondiis = .false.
      return
      end


      subroutine diis_initialise(ndim,beta)
c************************************************************************
c
c  Allocate workspace for DIIS solver
c
c************************************************************************
      implicit none

      integer ndim
      logical beta

INCLUDE(diis.fh)
INCLUDE(../m4/common/iofile)


      REAL derror, scale, st, ct, r
      integer iposit, nstore, mp
      logical ondiis
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis

INCLUDE(matrix_internal.fh)

      character*2 l1(8), l2(8)
      data l1 / 'f1','f2','f3','f4','f5','f6','f7','f8'/
      data l2 / 'e1','e2','e3','e4','e5','e6','e7','e8'/

      integer matrix_create

      integer i, ltri, ibuff1, ibuff2

c
c  Reset control parameters
c
      nstore = 0
      ondiis = .false.

      if (declared) return
c
c maxdiis pairs of buffers for fock and error vectors
c
      do i = 1, maxdiis

       ihandle(i,1) = matrix_create(ndim,ndim,l1(i),
     &   MATRIX_REAL)
       if(ihandle(i,1).le.0)
     &      call pg_error('failed to create fock ga ',i)

       ihandle(i,2) = matrix_create(ndim,ndim,l2(i),
     &   MATRIX_REAL)
       if(ihandle(i,2).le.0)
     &      call pg_error('failed to create error ga ',i)

         if (beta) then
c
c maxdiis pairs of buffers for beta fock and error vectors
c
       ihandle(i,3) = matrix_create(ndim,ndim,l1(i),
     &   MATRIX_REAL)
       if(ihandle(i,3).le.0)
     &      call pg_error('failed to create beta fock ga ',i)

       ihandle(i,4) = matrix_create(ndim,ndim,l2(i),
     &   MATRIX_REAL)
       if(ihandle(i,4).le.0)
     &      call pg_error('failed to create beta error ga ',i)

         end if

      enddo
c
c  remaining matrices are now handled outside the DIIS
c  solver
c
      declared = .true.

      end

      subroutine deallocate_diis_scratch (uhf)
INCLUDE(diis.fh)

      integer  o
      logical uhf
      if (supressed) return
      o = 0
      do i=1,maxdiis
         o = o + matrix_destroy(ihandle(i,1))
         o = o + matrix_destroy(ihandle(i,2))
         if (uhf) then
            o = o + matrix_destroy(ihandle(i,3))
            o = o + matrix_destroy(ihandle(i,4))
         endif
      enddo
      if(o.gt.0)call caserr('error destrying diis matrices')
      declared = .false.
      end


      block data block_diis_storage_2
INCLUDE(diis.fh)
      data declared /.false./
      data supressed /.false./
      end

c************************************************************************
c
c --- sets up and solves the DIIS equations
c       ( direct inversion of iterated subspace )
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
c   This version is based on matrix and vector objects
c
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written to ga in square form
c
c arguments:
c
c   AlpahaFock  latest fock matrix (overwritten with extrapolated)
c   BetaFock    latest fock matrix (overwritten with extrapolated)
c   Vectors     Current vectors
c
c   h0          scratch (triangle?)
c   h1          scratch (triangle)
c   h2          scratch (square)
c
c   ondiis1     logical  Out  flag specifying whether DIIS was used
c   nstore1     integer  Out  number of matrices used in extrapolation
c   diff0       REAL     Out  tester 
c   derror1     REAL     Out  diis error
c
c************************************************************************

      subroutine solve_diis(AlphaFock,BetaFock,Vectors,Overlap,
     &     Scratch, Scratch2,
     &     diis1,ondiis1,nvec,diff0,derror1)
      implicit none

      integer AlphaFock, BetaFock
      integer Vectors,Overlap,Scratch, Scratch2
      REAL diff0, derror1
      integer nvec
      logical diis1,ondiis1
INCLUDE(../m4/common/vcore)
INCLUDE(../m4/common/newscf)

      REAL derror, scale, st, ct, r
      integer iposit, nstore, mp
      logical ondiis
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis

c External
      include 'vectors_interface.fh'
      include 'matrix_interface.fh'

      integer igmem_alloc

c Local
      integer i, iq, iscr1, iscr2, iscr3
      integer AlphaOcc, BetaOcc, AlphaVmat, BetaVmat
      integer l1, l2
      logical uhf

      uhf = vec_unrestricted(Vectors)
      AlphaVmat = vec_alpha_coefficients(Vectors)
      AlphaOcc  = vec_alpha_occupations(Vectors)

      if(uhf)then
         BetaVmat = vec_beta_coefficients(Vectors)
         BetaOcc  = vec_beta_occupations(Vectors)
      endif
c
c harmonic alert
c
      l1 = matrix_dimension(AlphaVmat,1)
      l2 = (l1+1)*l1/2
      iscr1 = igmem_alloc(l2)
      iscr2 = igmem_alloc(l2)
      iscr3 = igmem_alloc(l2)

      call solve_diis0(AlphaFock,AlphaVmat,AlphaOcc,
     &     BetaFock, BetaVmat, BetaOcc,
     &     Overlap,
     &     Scratch, Scratch2,
     &     Q(iscr1),Q(iscr2),Q(iscr3),
     &     nvec,diff0,oscfprint(PR_DIIS),uhf,diis1)
   
      call gmem_free(iscr3)
      call gmem_free(iscr2)
      call gmem_free(iscr1)

      ondiis1 = ondiis
      derror1 = derror

      end

      subroutine solve_diis0(AlphaFock,AlphaVmat,AlphaOcc,
     &     BetaFock, BetaVmat, BetaOcc,
     &     Overlap,
     &     Scratch, Scratch2,
     &     h0,h1,h2,
     &     itot,diff0,out,uhf,diis_requested)
c
      implicit none
      
INCLUDE(../m4/common/sizes)
      REAL h0(*),h1(*),h2(*)
      integer AlphaFock,AlphaVmat,AlphaOcc
      integer BetaFock, BetaVmat, BetaOcc
      integer Overlap,Scratch,Scratch2
      REAL diff0
      logical out
      integer itot
      logical uhf
      logical diis_requested

INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/scra7)
INCLUDE(../m4/common/scfopt)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/harmon)
INCLUDE(diis.fh)
INCLUDE(../m4/common/prints)

c storage in /diisd/:
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c        is to be used in the scf program
c accdi1  diis comes in only when diff is less than this
c
      REAL derror, scale, st, ct, r
      integer iposit, nstore, mp
      logical ondiis
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis

INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'

c External

      logical opg_root

c Local
      integer num0, l1, l3, nmin, nmax, ndim, ipos, idum
      integer i, j, k
      integer ispace, iq, ixo, ifail, ipos1
      integer mp1, k1, k2, itmp
      REAL al2, ci, cij, sc
      logical occupied(maxat)
      REAL diff0b

      if (supressed) call caserr('diis - no matrices')

      derror=0.0d0
      num0 = newbas0
      l3 = num*num

      nmin=3
      nmax=maxdiis
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1)then
c
c too few degrees of freedom to run diis
c
         ondiis = .false.
         return
      endif
c
c ----sort out the vectors, store in ga
c
cc      call rdedx(h2,l3,iblkqa,idaf)
cc      call tdown(h2,ilifq,h2,ilifq,num0)
cc      call load_ga_from_square(ih_vec,h2,num)

c
c store input fock matrix to ga location (ipos,1)
c
      ipos=mod(nstore,nmax)+1

      idum = matrix_copy(AlphaFock,ihandle(ipos,1))
c
c ----- calculate the error vector, 
c       transform to mo basis
c
      idum=matrix_mult2(AlphaFock,AlphaVmat,ihandle(ipos,2),Scratch)

c      idum=matrix_print(AlphaFock)
c      idum=matrix_print_titled(ihandle(ipos,2),"error vector")

c
c       zero oo and vv blocks, and find absmax of ov
c

c      call pg_dscal_patch(ihandle(ipos,2),1,na,1,na,0.0d0)
c      call pg_dscal_patch(ihandle(ipos,2),na+1,num0,na+1,num0,0.0d0)
c      diff0 = gax_absmax_patch(ihandle(ipos,2),1,na,na+1,num0)

       l1 = matrix_dimension(ihandle(ipos,2),1)
       iq = hp(matrix_handle(ihandle(ipos,2)))

       ixo = hp(matrix_handle(AlphaOcc))

      do i = 1,num0
         occupied(i) = (QM(ixo+i-1).gt.0.001d0)
      enddo

      do i = 1,num0
         if(occupied(i))then
            do j = 1,num0
               if(occupied(j) )then
                  QM(iq + (i-1)*l1 + (j-1)) = 0.0d0
               endif
            enddo
         endif
      enddo

      do i = 1,num0
         if(.not. occupied(i))then
            do j = 1,num0
               if(.not. occupied(j) )then
                  QM(iq + (i-1)*l1 + (j-1)) = 0.0d0
               endif
            enddo
         endif
      enddo

      diff0 = 0.0d0
      do i = 1,num0
         if(occupied(i))then
            do j = 1,num0
               if(.not. occupied(j) )then
                  diff0 = max(diff0,
     &               abs(QM(iq + (i-1)*l1 + (j-1))))
               endif
            enddo
         endif
      enddo


      if(uhf)then

cc      call rdedx(h2,l3,iblkqb,idaf)
cc      call tdown(h2,ilifq,h2,ilifq,num0)
cc      call load_ga_from_square(ih_vec,h2,num)

      idum = matrix_copy(BetaFock,ihandle(ipos,3))
c
c ----- calculate the error vector, 
c       transform to mo basis
c
      idum=matrix_mult2(BetaFock,BetaVmat,ihandle(ipos,4),Scratch)
c
c       zero oo and vv blocks, and find absmax of ov
c

c      call pg_dscal_patch(ihandle(ipos,2),1,na,1,na,0.0d0)
c      call pg_dscal_patch(ihandle(ipos,2),na+1,num0,na+1,num0,0.0d0)
c      diff0 = gax_absmax_patch(ihandle(ipos,2),1,na,na+1,num0)

       l1 = matrix_dimension(ihandle(ipos,4),1)
       iq = hp(matrix_handle(ihandle(ipos,4)))

       ixo = hp(matrix_handle(BetaOcc))

      do i = 1,num0
         occupied(i) = (QM(ixo+i-1).gt.0.001d0)
      enddo

      do i = 1,num0
         if(occupied(i))then
            do j = 1,num0
               if(occupied(j) )then
                  QM(iq + (i-1)*l1 + (j-1)) = 0.0d0
               endif
            enddo
         endif
      enddo

      do i = 1,num0
         if(.not. occupied(i))then
            do j = 1,num0
               if(.not. occupied(j) )then
                  QM(iq + (i-1)*l1 + (j-1)) = 0.0d0
               endif
            enddo
         endif
      enddo

      diff0b = 0.0d0
      do i = 1,num0
         if(occupied(i))then
            do j = 1,num0
               if(.not. occupied(j) )then
                  diff0b = max(diff0b,
     &                 abs(QM(iq + (i-1)*l1 + (j-1))))
               endif
            enddo
         endif
      enddo

      diff0 = max(diff0,diff0b)

      endif

      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) then
c
c components of error vector > accdi1, turn off diis
c
         nstore=0
         itot=0
         mp=0
         ondiis=.false.
         if(opg_root())write(iwr,*)
     &        'DIIS suppressed (tester too high)', diff0
         goto 9999
      endif
c
c ------ transform back to the ao basis
c
cc      call mult2t_ga(ihandle(ipos,2),ih_vec, ih_scr, num)
cc      call mult2_ga(ihandle(ipos,2),ih_ov, ih_scr, num)

      idum = matrix_mult2t(ihandle(ipos,2),AlphaVmat, Scratch2,Scratch)
      idum = matrix_mult2(Scratch2,Overlap,ihandle(ipos,2),Scratch)

      if(uhf)then
         idum = 
     &     matrix_mult2t(ihandle(ipos,4),BetaVmat, Scratch2,Scratch)
         idum = 
     &     matrix_mult2(Scratch2,Overlap,ihandle(ipos,4),Scratch)
      endif

      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1

      ipos1=ipos-1
      mp=iky(ipos)

      if(ipos1.ge.1) then
         do 50 i=1,ipos1
            mp=mp+1
            st(mp)=matrix_dot_product(ihandle(i,2),ihandle(ipos,2))
            if(uhf)then
               st(mp)=st(mp) + 
     &              matrix_dot_product(ihandle(i,4),ihandle(ipos,4))
            endif
 50      continue
      endif
c
c store diagonal term
      mp=mp+1
      st(mp)=matrix_dot_product(ihandle(ipos,2),ihandle(ipos,2))
      if(uhf)then
         st(mp)=st(mp) + 
     &        matrix_dot_product(ihandle(ipos,4),ihandle(ipos,4))
      endif
      itot=nstore
      if(nstore.gt.nmax)itot=nmax

      ipos1=ipos+1
      if(ipos1.le.itot) then
         do 110 i=ipos1,itot
            mp=mp+i-1
            st(mp)=matrix_dot_product(ihandle(i,2),ihandle(ipos,2))
            if(uhf)then
               st(mp)=st(mp) + 
     &              matrix_dot_product(ihandle(i,4),ihandle(ipos,4))
            endif
 110     continue
      endif

      if (.not.diis_requested) then
c
c        stored all the required data for potential use in
c        a following phase but do not actually do anything
c
         if (opg_root()) write(iwr,*)
     &      "DIIS not active, store data only"
         goto 9999
      endif
c
      if(nstore.le.nmin)then
c
c not enough triangles yet
c
         ondiis = .false.
         if (opg_root()) write(iwr,*)
     &      "DIIS supressed (subspace too small)"
         goto 9999

      endif
c
c --- now solve the diis equations
c
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)
c
c --- scale the matrix
c
      do 125 i=1,ndim-1
 125  scale(i)=dsqrt(st(ikyp(i)))
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1

      if(out)write(iwr,126)(st(k),k=1,mp)
 126  format(/1x,'diis matrix'/(1x,10e11.4))

      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h2,h2(ndim+1),ifail)
      if(ifail.ne.0) then
c
c diis failure - exit, but retain the current point
c
         if(out)write(iwr,129)ifail
 129     format(//1x,'diis failure in f04atf .... ifail= ',i2)

c swap matrix handles to bring current h,e to ipos=1
c this has not been checked
c
         itmp = ihandle(ipos,1)
         ihandle(ipos,1) = ihandle(1,1)
         ihandle(1,1) = itmp

         itmp = ihandle(ipos,2)
         ihandle(ipos,2) = ihandle(1,2)
         ihandle(1,2) = itmp

         itmp = ihandle(ipos,3)
         ihandle(ipos,3) = ihandle(1,3)
         ihandle(1,3) = itmp

         itmp = ihandle(ipos,4)
         ihandle(ipos,4) = ihandle(1,4)
         ihandle(1,4) = itmp

         nstore=1
         itot=1
         mp=iky(ipos)+ipos
         st(1)=st(mp)

         ondiis = .false.
         goto 9999

      endif

      do 114 i=1,ndim
 114     ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(out)write(iwr,128)(ct(k),k=1,ndim)
 128  format(/1x,'solution'/(1x,10e11.4))
c
c construct extrapolated fock matrix
c
      idum=matrix_assign_zero(AlphaFock)
      do k=1,itot
         idum=matrix_daxpy(ct(k),ihandle(k,1),AlphaFock)
      enddo
      if(uhf)then
         idum=matrix_assign_zero(BetaFock)
         do k=1,itot
            idum=matrix_daxpy(ct(k),ihandle(k,3),BetaFock)
         enddo
      endif
c
      if(.not.ondiis) kcount=0
      ondiis=.true.
c
c exit diis 
c
9999  continue

      return
      end

c
      REAL function check_ov(AlphaTFock,BetaTFock,Vectors)
      implicit none
      integer AlphaTFock,BetaTFock, Vectors

INCLUDE(../m4/common/sizes)

      include 'vectors_interface.fh'
      include 'matrix_interface.fh'
INCLUDE(matrix_internal.fh)

      logical occupied(maxorb)
      integer i, j, n
      integer ixo, ixt
      integer Occ, Vmat
      logical uhf
      REAL check_ov_beta

      Occ = vec_alpha_occupations(Vectors)
      Vmat = vec_alpha_coefficients(Vectors)
      uhf = vec_unrestricted(Vectors)

c Harmonic alert
      n = matrix_dimension(Vmat,1)

      ixo = hp(matrix_handle(Occ))
      ixt = hp(matrix_handle(AlphaTfock))

      do i = 1,n
         occupied(i) = (QM(ixo+i-1).gt.0.001d0)
      enddo

      check_ov = 0.0d0

c      write(iwr,*)'check_ov occs',(occupied(i),i=1,n)

      do i = 1,n
         if(occupied(i))then
            do j = 1,n
               if( .not. occupied(j) )then
                  check_ov = 
     &         max(abs(QM(ixt + (i-1)*n + (j-1))),check_ov)
               endif
            enddo
         endif
      enddo

      if(uhf)then

         Occ = vec_beta_occupations(Vectors)
         Vmat = vec_beta_coefficients(Vectors)

         ixo = hp(matrix_handle(Occ))
         ixt = hp(matrix_handle(BetaTfock))

         do i = 1,n
            occupied(i) = (QM(ixo+i-1).gt.0.001d0)
         enddo

         check_ov_beta = 0.0d0

         do i = 1,n
            if(occupied(i))then
               do j = 1,n
                  if( .not. occupied(j) )then
                     check_ov_beta = 
     &            max(abs(QM(ixt + (i-1)*n + (j-1))),check_ov)
                  endif
               enddo
            endif
         enddo
c
c Note this a change to the GAMESS-UK definition
c
         check_ov = max(check_ov,check_ov_beta)

      endif

      end


      subroutine level_shift(AlphaTFock,BetaTFock,
     &   Vectors,shifta,shiftb)
      implicit none

      integer AlphaTFock, BetaTFock, Vectors
      REAL shifta, shiftb

INCLUDE(../m4/common/sizes)

      include 'vectors_interface.fh'
      include 'matrix_interface.fh'
INCLUDE(matrix_internal.fh)

      logical occupied(maxorb)
      integer i, j, n
      integer ixo, ixf
      integer Occ, Vmat
      logical uhf

      Occ = vec_alpha_occupations(Vectors)
      Vmat = vec_alpha_coefficients(Vectors)
      uhf = vec_unrestricted(Vectors)

c Harmonic alert
      n = matrix_dimension(Vmat,1)

      ixo = hp(matrix_handle(Occ))
      ixf = hp(matrix_handle(AlphaTFock))

      do i = 1,n
         occupied(i) = (QM(ixo+i-1).gt.0.001d0)
      enddo

      do i = 1,n
         if(.not. occupied(i))then
            QM(ixf + (i-1)*n + (i-1)) =
     &           QM(ixf + (i-1)*n + (i-1)) + shifta
         endif
      enddo

      if(uhf)then
         Occ = vec_beta_occupations(Vectors)
         Vmat = vec_beta_coefficients(Vectors)

         ixo = hp(matrix_handle(Occ))
         ixf = hp(matrix_handle(BetaTFock))

         do i = 1,n
            occupied(i) = (QM(ixo+i-1).gt.0.001d0)
         enddo

         do i = 1,n
            if(.not. occupied(i))then
               QM(ixf + (i-1)*n + (i-1)) =
     &              QM(ixf + (i-1)*n + (i-1)) + shiftb
            endif
         enddo

      endif
      end


c***********************************************************************
      subroutine save_orbitals(Vectors, AlphaDensity, BetaDensity)
c***********************************************************************
c
c interface to scfsav: save mo*s + density + orbital energies 
c                      + occupations
c
c  Vectors   vectors handle   In
c  Density   matrix handle    In
c
c***********************************************************************

      implicit none

      integer Vectors, AlphaDensity, BetaDensity

INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/atmol3)
_IFN(ma)
INCLUDE(../m4/common/vcore)
_ENDIF
INCLUDE(../m4/common/dump3)

INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'
      include 'vectors_interface.fh'

      integer Occ, iocc, Eig, ieig, Vec, ivec
      integer l1, l2
      integer idens
      integer na, ndaf, idum
      
      logical uhf

      integer igmem_alloc
      external igmem_alloc

      uhf = vec_unrestricted(Vectors)

      Occ  = vec_alpha_occupations(Vectors)
      iocc  = hp(matrix_handle(Occ))
      Eig = vec_alpha_eigenvalues(Vectors)
      ieig  = hp(matrix_handle(Eig))
      Vec = vec_alpha_coefficients(Vectors)
      ivec  = hp(matrix_handle(Vec))

      l1 = matrix_dimension(Occ,1)
      l2 = (l1+1)*l1/2
c
c Triangulate density
c
      idens = igmem_alloc(l2)
      idum = matrix_get_to_triangle(AlphaDensity,Q(idens))

c      subroutine scfsav(q,p,e,pop,ndaf,l1,l2,iblkp,iblke)

      ndaf = mouta
      call scfsav(QM(ivec),Q(idens),
     &     QM(ieig),QM(iocc),
     &     ndaf,l1,l2,ibl3pa,ibl3ea)

      if(uhf)then

c Save Beta set

         Occ  = vec_beta_occupations(Vectors)
         iocc  = hp(matrix_handle(Occ))
         Eig = vec_beta_eigenvalues(Vectors)
         ieig  = hp(matrix_handle(Eig))
         Vec = vec_beta_coefficients(Vectors)
         ivec  = hp(matrix_handle(Vec))

         idum = matrix_get_to_triangle(BetaDensity,Q(idens))

         ndaf = moutb

         call scfsav(QM(ivec),Q(idens),
     &        QM(ieig),QM(iocc),
     &        ndaf,l1,l2,ibl3pb,ibl3eb)

      endif
      call gmem_free(idens)
      end


      subroutine default_conv(uhf)

      implicit none

      logical uhf

INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/iofile)

      integer iphase

      if(nphase .ne. 0)then
         write(iwr,*)'Using user-defined convergence scheme'
         return
      endif

      nphase = 2

      iphase=1
      lock_vec(iphase) = lock_vec(0)
      diis(iphase) = diis(0)
      shift(iphase,1) = shift(0,1)
      shift(iphase,2) = shift(0,2)
      nnext(iphase) = 1
      tester_chk(iphase,1) = CONV_BELOW
      tester_val(iphase,1) = 0.1d0
      nextphase(iphase,1) = 2
c
c  As phase 1, but turn off level shifter
c
      iphase=2
      lock_vec(iphase) = lock_vec(0)
      diis(iphase) = diis(0)
      shift(iphase,1) = 0.0d0
      shift(iphase,2) = 0.0d0
      nnext(iphase) = 1
      tester_chk(iphase,1) = CONV_BELOW
      tester_val(iphase,1) = tester_val(0,1)
      abs_dele_chk(iphase,1) = CONV_BELOW
      abs_dele_val(iphase,1) = abs_dele_val(0,1)
      nextphase(iphase,1) = 0

      end


      subroutine print_conv(uhf)

      implicit none

      logical uhf

INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/iofile)

      logical opg_root

      integer i,j

      write(iwr,*)
      write(iwr,*)'Convergence procedure'
      write(iwr,*)'*********************'

      do i =1,nphase
         write(iwr,*)
         write(iwr,*)'Phase ',i 
         if(uhf)then
            write(iwr,*)'  Level ',
     &           shift(i,1),shift(i,2)
         else
            write(iwr,*)'  Level ',shift(i,1)
         endif
         if(diis(i))write(iwr,*)'  DIIS'
         if(lock_vec(i))write(iwr,*)'  Lock'
         if(restore_vec(i))write(iwr,*)'  Restore'
         if(new_diis(i))write(iwr,*)'  NewDIIS'
         if(new_extrap(i))write(iwr,*)'  NewExtrapolate'
         if(extrap(i))write(iwr,*)'  Extrapolate',
     &        extrap_tol(i),extrap_coef(i)
         do j = 1,nnext(i)

            if(nextphase(i,j).eq.-1)then
               write(iwr,*)'#  Abort calculation'
            else if(nextphase(i,j).eq.0)then
               write(iwr,*)'#  Converge calculation'
            else
               write(iwr,*)'#  Switch to phase ',nextphase(i,j)
            endif

            write(iwr,*)'   next ',nextphase(i,j)

            if(tester_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'     Tester below ',tester_val(i,j)
            else if (tester_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'     Tester above ',tester_val(i,j)
            endif

            if(dele_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#    Energy change'
               write(iwr,*)' dE below ',
     &              dele_val(i,j)
            else if (dele_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Energy change '
               write(iwr,*)' dE above ',
     &              dele_val(i,j)
            endif

            if(abs_dele_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#     Absolute energy change'
               write(iwr,*)'    dEabs below',
     &              abs_dele_val(i,j)
            else if (abs_dele_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Absolute energy change'
               write(iwr,*)'    dEabs above',
     &              abs_dele_val(i,j)
            endif

            if(ncyc_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#     Cycle count'
               write(iwr,*)'   ncyc below ',ncyc_val(i,j)
            else if (ncyc_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Cycle count'
               write(iwr,*)'    ncyc above ',ncyc_val(i,j)
            endif

         enddo
      enddo

      end

c
c Check convergence
c
      integer function check_conv(uhf,iphase,tester,dele,niter,
     &     o1,o2,o3,o4)

      implicit none

      logical uhf
      integer iphase
      REAL tester
      REAL dele
      integer niter


INCLUDE(../m4/common/newscf)      
INCLUDE(../m4/common/iofile)

      logical o1(maxphase), o2(maxphase), o3(maxphase), o4(maxphase)
      integer i
c
c  default is no change yet
c
      check_conv = iphase

      do i = 1, nnext(iphase)
c
c perform convergence/abort tests
c if more than one test is passed, the last one takes
c effect
c
         o1(i) = (tester_chk(iphase,i) .eq. CONV_INACTIVE)
         if(tester_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        tester .gt. tester_val(iphase,i))o1(i)=.true.
         if(tester_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        tester .lt. tester_val(iphase,i))o1(i)=.true.

         o2(i) = (dele_chk(iphase,i) .eq. CONV_INACTIVE)
         if(dele_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        dele .gt. dele_val(iphase,i))o2(i)=.true.
         if(dele_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        dele .lt. dele_val(iphase,i))o2(i)=.true.

         o3(i) = (abs_dele_chk(iphase,i) .eq. CONV_INACTIVE)
         if(abs_dele_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        abs(dele) .gt. abs_dele_val(iphase,i))o3(i)=.true.
         if(abs_dele_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        abs(dele) .lt. abs_dele_val(iphase,i))o3(i)=.true.

         o4(i) = (ncyc_chk(iphase,i) .eq. CONV_INACTIVE)
         if(ncyc_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        niter .gt. ncyc_val(iphase,i))o4(i)=.true.
         if(ncyc_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        niter .lt. ncyc_val(iphase,i))o4(i)=.true.

         if(o1(i) .and. o2(i) .and. o3(i) .and. o4(i))then
            check_conv = nextphase(iphase,i)
         endif

      enddo

      end

_IFN(ga)
      subroutine pg_error(text,code)
      character text*(*)
      integer code
      call dlc_error(text,'abort')
      end
_ENDIF
      subroutine overlap_check(NewVec, OldVec)
      integer NewVec, OldVec

INCLUDE(../m4/common/sizes)

INCLUDE(matrix_internal.fh)
      include 'matrix_interface.fh'
      include 'vectors_interface.fh'
      
      REAL fnew(maxorb)
      REAL fold(maxorb)
      REAL sum, sumnew, sumold
      REAL res(maxorb,maxorb)

      integer i, j, k, New, inew, Old, iold

      New = vec_alpha_coefficients(NewVec)
      inew  = hp(matrix_handle(New))

      Old = vec_alpha_coefficients(OldVec)
      iold  = hp(matrix_handle(Old))

      l1 = matrix_dimension(New,1)

      do i = 1,l1
         sumnew = 0.0d0
         sumold = 0.0d0
         do j = 1,l1
            sumold = sumold + QM( iold + (i-1)*l1 + (j-1)) * 
     &           QM( iold + (i-1)*l1 + (j-1))
            sumnew = sumnew + QM( inew + (i-1)*l1 + (j-1)) * 
     &           QM( inew + (i-1)*l1 + (j-1))
         enddo
         fnew(i) = sqrt(sumnew)
         fold(i) = sqrt(sumold)
      enddo

      do i = 1,l1
         do j = 1,l1
            sum = 0.0d0
            do k = 1,l1
               sum = sum + QM( iold + (i-1)*l1 + (k-1)) * 
     &              QM( inew + (j-1)*l1 + (k-1))
            enddo
            res(i,j) = sum / (fold(i) * fnew(j))
         enddo
      enddo

      call prsq(res,l1,l1,maxorb)

      end

      subroutine make_rdmat(rdmat,uhf,AlphaDensity,BetaDensity)
      implicit none

      include 'matrix_interface.fh'      
INCLUDE(../m4/common/vcore)

c Arguments
      REAL rdmat(*)
      logical uhf 
      integer AlphaDensity,BetaDensity

c External function
      integer igmem_alloc
      external igmem_alloc

c Local
      integer idum
      integer iscr, l2, n
      character*8  zscftp

      n = matrix_dimension(AlphaDensity,1)

      l2 = ((n+1)*n)/2
      if(uhf)then
         zscftp = 'uhf'
         iscr = igmem_alloc(2*l2)
         idum=matrix_get_to_triangle(AlphaDensity,
     &        Q(iscr))
         idum=matrix_get_to_triangle(BetaDensity,
     &        Q(iscr+l2))
      else
         zscftp = 'rhf'
         iscr = igmem_alloc(l2)
         idum=matrix_get_to_triangle(AlphaDensity,
     &        Q(iscr))
      endif

      call mkrdmt(zscftp,rdmat,Q(iscr),l2,0)

      call gmem_free(iscr)

      end

