Subroutine fock_build(alphadensity,betadensity, &
     direct,    &
     alphafock,betafock,hcore,uhf, &
     ehf, ehf1,ehf2,edft, &
     nquad, alpha_quad, beta_quad, &
     core)
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
  Use distributed_matrices
  Use distributed_vectors
  Use newscf_modules

  Implicit None

  !$$$      integer alphadensity,betadensity
  Type( matrix ) :: alphadensity,betadensity
  Logical direct
  !$$$      integer alphafock,betafock,hcore
  Type( matrix ) :: alphafock,betafock,hcore
  Real( wp ) ehf, ehf1, ehf2, edft
  Integer nquad
  Real( wp ) alpha_quad, beta_quad
  Logical uhf
  Real( wp ) core(*)

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

  !$$$      include 'matrix_interface.fh'

  Logical opg_root
  Integer igmem_alloc
  Real( wp ) tracep, dum

  ! local variables

  Integer l1, l2, l3, l22
  Logical out, outon
  Integer idum, i
  Integer idens,  iscr
  Integer idensb, ifockb
  Real( wp ) etmp, dft_accu

  Logical o2e, odft, ocoul, oexch
  Real( wp ) facex

  Logical delta
  Integer nshtri
  Integer kkk

  l1 = num
  l2 = l1*(l1+1)/2
  l3 = l1*l1

  l22 = 2*l2

  out = nprint == 5
  outon = nprint /= -5

  Call pg_synch(4444)

  edft = 0.0_wp

  delta = .False.

  If(direct)Then

     If(delta)Then
        ! -- note delta density not implemented
     End If

     dlnmxd=0.0_wp
     dlntol=tolitr(3)

     !         if(opg_root())then
     !	    write(iwr,*)'direct fock build dlntol: ',dlntol
     !         End If

     !
     !  build reduced integral prefactor matrix
     !  - can probably skip this on subsequent cycles
     !
     Call start_time_period(tp_rdmat)
     nshtri=ikyp(nshell)

     irdmat = igmem_alloc(2*nshtri)
     Call make_rdmat(q(irdmat),uhf,alphadensity,betadensity)

     !         write(6,*)'reduced density'
     !         call prtri(q(irdmat),nshell)

     !         iprefa = igmem_alloc(nshtri)
     iprefa = irdmat + nshtri
     Call rdmake(q(iprefa))

     !         write(6,*)'reduced prefactor'
     !         call prtri(q(iprefa),nshell)

     dlnmxd=-9999999.0_wp
     Do  kkk=1,nshtri
        If(q(irdmat+kkk-1).Gt.dlnmxd) &
             dlnmxd=q(irdmat+kkk-1)
     End Do

     Call end_time_period(tp_rdmat)

  End If

  !
  ! modify fock builder options according to
  ! required dft scheme (hf exchange, fitted coulomb etc
  !
  idum = cd_set_2e()

  !
  ! allocate triangles for hstar etc
  ! hstaru assumes matrices are contiguous
  !
  If(uhf)Then
     idens = igmem_alloc(l22)
     idensb = idens + l2
     ifock = igmem_alloc(l22)
     ifockb = ifock + l2
  Else
     idens = igmem_alloc(l2)
     ifock = igmem_alloc(l2)
  End If

  ! - read in fock builders
  idmat = idens

  Call matrix_get_to_triangle(alphadensity,q(idens))
  If(uhf)Then
     Call matrix_get_to_triangle(betadensity,q(idensb))
  End If

  iscr = igmem_alloc(l2)

  Do i = 1,l2
     q(ifock + i -1) = 0.0_wp
     q(iscr + i -1) = 0.0_wp
  End Do
  !
  If (cd_2e())Then
     odft = .True.
     If(uhf)idum = cd_uks()
     !
     !     for simplicity, only implement full coulomb version
     !     here for the moment
     !
     If(cd_hf_exchange())Then
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
  o2e = .Not.(odft.And..Not.oexch.And..Not.ocoul)
  !
  If(uhf)Then
     Do i = 1,l2
        q(ifockb + i -1) = 0.0_wp
     End Do
     !
     If(o2e)Then
        If(direct)Then
           irest=0
           Call dhstaru(core,q(idens),q(idensb), &
                q(ifock),q(ifockb),q(iprefa),q(irdmat),irest)
        Else
           Call hstaru(q(idens),q(idensb), &
                q(ifock),q(ifockb),nopk, &
                facex,ocoul,oexch)
           !            else
           !               call vclr(q(ifock),1,l2)
           !               call vclr(q(ifocka),1,l2)
        End If
     End If
  Else
     !
     ! rhf case
     If(o2e)Then
        If(direct)Then
           irest = 0
           Call dhstar(core,q(ifock),q(idens),q(iprefa),q(irdmat),irest)
        Else
           ! problem if on cray, convex or titan - m4 removed
           Call hstar(q(idens),q(ifock), &
                q(iscr),nopk)
        End If
     End If
  End If

  !
  ! restore fock builder options
  !
  idum = cd_reset_2e()

  ! - skip zora for now

  !     ----- symmetrize skeleton fock matrix -----
  !
  !     -h- at x(i10)
  !     scratch area at x(i20)
  !
  Call symh(q(ifock),q(iscr),iky,0,0)
  If(uhf)Call symh(q(ifockb),q(iscr),iky,0,0)
  Call gmem_free(iscr)

  If (out) Then
     Write (iwr,9088)
     If(uhf)Write (iwr,9108)
     Call prtril(q(ifock),l1)
     If(uhf)Then
        Write (iwr,9128)
        Call prtril(q(ifockb),l1)
     End If
  End If

9088 Format(/20x,23('-')/20x,'symmetrized fock matrix'/20x,23( &
       '-'))
9108 Format(//1x,129('-')/// &
       50x,'----- alpha set -----')
9128 Format(//1x,129('-')/// &
       50x,'----- beta set -----')

  Call matrix_set_from_triangle(alphafock,q(ifock))
  If(uhf) &
       Call matrix_set_from_triangle(betafock,q(ifockb))

  !
  !     ----- read in core hamiltonian matrix
  !           and calculate hf energy -----
  !
  !     -h0- at x(i20)
  !     - h- at x(i10)
  !     - e- at x(i30)
  !

  !      call vadd(q(i10),1,q(jblkf),1,q(i10),1,nx)

  Call matrix_add(alphafock,alphafock,hcore)
  If(uhf) Call matrix_add(betafock,betafock,hcore)
  !
  !  ehf1 = trp(p*h_1)
  !
  Call matrix_get_to_triangle(alphadensity,q(idens))
  If(uhf) Call matrix_get_to_triangle(betadensity,q(idensb))

  Call matrix_get_to_triangle(hcore,q(ifock))
  ehf1 = tracep(q(idens),q(ifock),l1)
  If(uhf)ehf1 = ehf1 + tracep(q(idensb),q(ifock),l1)

  If(cd_active())Then
     !
     ! ccpdft: evaluate kohn sham energy expression
     !
     If(cd_hf_exchange() .Or. cd_hf_coulomb())Then
        !
        ! coulomb or exchange operator is in q(ifock) augmented by h_1, 
        ! compute energy using hf expression (maybe without exchange)
        !
        Call matrix_get_to_triangle(alphafock,q(ifock))
        ehf2 = tracep(q(idens),q(ifock),l1)
        If(uhf)Then
           Call matrix_get_to_triangle(betafock,q(ifockb))
           ehf2 = ehf2 +  &
                tracep(q(idensb),q(ifockb),l1)
        End If
        etmp = (ehf1+ehf2)*0.5_wp
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

     Call matrix_get_to_triangle(alphafock,q(ifock))
     If(uhf) &
          Call matrix_get_to_triangle(betafock,q(ifockb))

     dft_accu = 0.0_wp

     idum = cd_energy(c,q(ifock),q(ifockb), &
          q(idens),q(idensb), &
          edft,core,core,outon,dft_accu,iwr &
          )


     Call matrix_set_from_triangle(alphafock,q(ifock))
     If(uhf) &
          Call matrix_set_from_triangle(betafock,q(ifockb))

     ehf = etmp+edft
     Call cd_get_dftresults(nquad,alpha_quad,beta_quad)

  Else
     !
     !  hartree-fock energy expression
     !  e = 1/2 [ trp(p*h_1) + trp(p*(h_1 + h_2)) ]
     !                ehf1            ehf2
     !
     Call matrix_get_to_triangle(alphafock,q(ifock))
     ehf2 = tracep(q(idens),q(ifock),l1)
     If(uhf)Then
        Call matrix_get_to_triangle(betafock,q(ifockb))
        ehf2 = ehf2 +  &
             tracep(q(idensb),q(ifockb),l1)
     End If
     ehf = (ehf1+ehf2)*0.5_wp
  End If

  Call gmem_free(ifock)
  Call gmem_free(idens)

  If(direct)Then
     Call gmem_free(irdmat)
  End If

  ! - nb skip code for crystal field
End Subroutine fock_build

