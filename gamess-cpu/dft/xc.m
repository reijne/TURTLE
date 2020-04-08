c
c  $Author: hvd $ $Revision: 6206 $ $Date: 2010-10-31 18:05:53 +0100 (Sun, 31 Oct 2010) $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/dft/xc.m,v $
c
C  *********************************************************************
c  * exquad:
C  * Generate exchange/correlation contribution using numerical 
c  * quadrature  - unified energy/ksmatrix/gradient version
C  *********************************************************************
c
c---- memory counting routines -----------------------------------------
c     
      subroutine memreq_exquad(memory_fp,memory_int,
     &     igrid,
     &     ao_tag,nao,nvec,naocc,nbocc,npert,
     &     geometric_pert_sw,
     &     prpt,prwt,
     &     e_sw,ks_sw,grad_sw,
     &     dksm_exp_sw,rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &     mxp0,ao_in_sw,mo_in_sw,accuracy,iout)

      implicit none
c
c     Parameters  
c
INCLUDE(common/dft_parameters)
c
c     In variables  
c
INCLUDE(common/dft_api)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_xc)
INCLUDE(common/dft_dij)
      integer igrid             ! which grid to use
      integer ao_tag
      integer nao               ! the number of AOs
      integer nvec              ! the number of MOs
      integer naocc             ! the number of occupied alpha orbitals
      integer nbocc             ! the number of occupied beta  orbitals
      integer npert             ! the number of perturbations
      logical geometric_pert_sw ! true if geometric perturbations are
                                ! considered, otherwise the source of
                                ! the perturbation is an external field.
      REAL accuracy
      logical e_sw              ! calculate the energy
      logical ks_sw             ! calculate Kohn-Sham matrix elements
      logical grad_sw           ! calculate energy gradient contribution
      logical dksm_exp_sw       ! calculate the derivative of the Kohn-
                                ! Sham matrix with respect to the 
                                ! explicit dependencies on the nuclear
                                ! coordinates
      logical lhs_sw            ! calculate CHF "left-hand-side" 
                                ! contributions
      logical rhs_sw            ! calculate CHF "right-hand-side" 
                                ! contributions
      logical dksm_sw           ! calculate the "wavefunction" term to
                                ! the derivative of the Kohn-Sham matrix
      logical hess_sw           ! calculate energy 2nd derivative 
                                ! contribution
      logical ao_in_sw          ! true if adens and bdens present
      logical mo_in_sw          ! true if avec and bvec present
      logical extwr_sw
      integer mxp0
      integer iout
c
c     Out variables  
c
      REAL kma   ,kmb                       
      REAL grad                       
      REAL ba_ao                         ! the alpha b vector (Au=b)
      REAL bb_ao                         ! the beta  b vector (Au=b)
      REAL ba_mo                         ! the alpha b vector (Au=b)
      REAL bb_mo                         ! the beta  b vector (Au=b)
      REAL ga_ao                         ! the alpha Au matrix (Au=b)
      REAL gb_ao                         ! the beta  Au matrix (Au=b)
      REAL ga_mo                         ! the alpha Au matrix (Au=b)
      REAL gb_mo                         ! the beta  Au matrix (Au=b)
      REAL fa_ao                         ! the alpha perturbed fock
                                         ! matrix
      REAL fb_ao                         ! the beta  perturbed fock
                                         ! matrix
      REAL fa_mo                       
      REAL fb_mo                       
      REAL fxa_ao                        ! explicit derivative alpha
                                         ! fock matrix
      REAL fxb_ao                        ! explicit derivative beta
                                         ! fock matrix
      REAL fxa_mo                       
      REAL fxb_mo                       
      REAL hess             
      REAL smat

C  Work arrays
      REAL memory_fp(*)
      integer memory_int(*)
      REAL prpt(ngtypes,*)
      REAL prwt(ngtypes,*)
      integer bfn_radii_pt
      integer bfn_radii_ct
      integer nprt
c
c still to be allocated
c
      REAL arad(max_atom)
      integer first_bf(max_atom+1)
      integer active_bfn_list_pt
      integer active_bfn_indx_pt
      integer active_bfn_atms_pt
      integer active_wgh_atms_pt
      integer active_chf_pert_pt
      integer n_active_bfn
      integer n_active_bfn_atm
      integer n_active_wgh_atm
      integer n_active_chf_prt
      integer tau_pt
      integer xc_e_pt
      integer xc_v_pt, xc_dv_pt, xc_dt_pt
      integer xc_h_pt, xc_dh_pt
      integer wt_pt, wt2_pt, gwt_pt, g2wt_pt
      integer bfn_val_pt, bfng_val_pt, bfn_hess_pt, bfn_3rd_pt
      integer amo_val_pt, amo_grad_pt, amo_hess_pt
      integer bmo_val_pt, bmo_grad_pt, bmo_hess_pt
      integer gamma_pt
      integer drho_pt,  dgrho_pt, dgamma_pt
      integer ddrhoa_pt, ddgrhoa_pt
      integer ddrhob_pt, ddgrhob_pt
      integer apts_pt, awpt_pt, iprm_pt
cDEBUG
c     integer push_memory_estimate, pop_memory_estimate
c     integer igmem_get_estimate, igmem_max_estimate
c     integer icest, imest
c     integer imemcount
cDEBUG
c
C     Functions 
c
      logical opg_root, CD_gradcorr, CD_kinetic
      integer ipg_dlbtask, ipg_nodeid, ipg_nnodes
_IF(single)
      REAL sdot
      integer isamax
_ELSEIF(hp700)
      REAL `vec_$ddot'
      integer idamax
_ELSE
      REAL ddot
      integer idamax
_ENDIF
      integer null_memory
      integer allocate_memory
      integer allocate_memory2
      integer incr_memory
      integer incr_memory2
      integer memreq_pg_dgop, memreq_pg_igop

C Local variables
      integer idum
      integer latm,lrad,atmt,atom_num,i,j,k
      integer nradpt_num(max_gtype),mxang
      integer hi, lo, lang0
      integer idwght, idbfn, idmo, idden, idfun
      integer na_mo, nb_mo
      REAL rpt,rwt
      REAL atom_xce(max_atom),atom_den(max_atom,2)
      integer atom_pts(max_atom)
      REAL xc_e
      REAL del_psi, del_rho
      REAL rshell, rnear
      integer next, ltri
      REAL fac
      integer mxp
      integer npts
      integer nbatch, ntot, npack, npnew,ipnew
      integer iiwrk1
      integer iwrk1, iwrk2,  iwrk3,  iwrk4, iwrk5, iwrk6, iwrk7, iwrk8
      integer iwrk9, iwrk10, iwrk11, iwrk12
      integer nradmx,nradtot,itmx
      integer freq
      integer gridt
      integer irad, ig
      integer imax_active_atm
      integer imax_active_bfn
c
c local switches
c
      logical gradcorr_sw, kinetic_sw, gradkin_sw
      logical extout_sw
      logical gwt_avail_sw
      logical eval_mo_sw     ! true if MOs will be evaluated
      logical den_mo_sw      ! true if rho will be calculated from MOs
      logical lhs_mo_sw      ! true if lhs will be calculated from MOs
      logical rhs_mo_sw      ! true if rhs will be calculated from MOs
      logical dksm_mo_sw     ! true if the perturbed Kohn-Sham matrix 
                             ! will be calculated from MOs
      logical dksm_exp_mo_sw ! true if the explicit derivative of the
                             ! Kohn-Sham matrix will be calculated from
                             ! MOs
      logical hess_mo_sw     ! true if d2Exc will be calculated from MOs

      REAL cpu_prev
      REAL cpu1(3)
      integer ieshll, irshll, inmtyp(0:max_gtype)
c     integer atom_num_by_tag
chvd
c     integer ncurrrad(max_gtype)
      integer ntheta(max_gtype), nphi(max_gtype), nomega(max_gtype)
      integer nctht, ncphi, ncang, iradzone
chvd
c     
c     ============ temp timing analysis ==========
c     
      REAL t_now
      REAL t_start, t_then
      REAL t_becke, t_bas, t_den, t_pack
      REAL t_geom, t_xc, t_ks, t_grad
      LOGICAL time_cond

_IF(qsh)
      integer nbatch_qsh, npts_qsh
      REAL exc_qsh, den_qsh
      character*6 fname
      logical ocrap
_ENDIF

      logical odb
      character*4 fnm
      character*13 snm
c     
c     ==============================================

INCLUDE(../m4/common/parallel)
INCLUDE(../m4/common/parcntl)
c
      data fnm /'xc.m'/
      data snm /'memreq_exquad'/
c
C     
C     api entry point
c     
c
c define distinct timers for energy/gradient cases
c
c
c     initialise addresses
c
      bfn_val_pt  = null_memory()
      bfng_val_pt = null_memory()
      bfn_hess_pt = null_memory()
      bfn_3rd_pt  = null_memory()
      amo_val_pt  = null_memory()
      bmo_val_pt  = null_memory()
      amo_grad_pt = null_memory()
      bmo_grad_pt = null_memory()
      amo_hess_pt = null_memory()
      bmo_hess_pt = null_memory()
      tau_pt      = null_memory()
      xc_e_pt     = null_memory()
      xc_v_pt     = null_memory()
      xc_dv_pt    = null_memory()
      xc_dt_pt    = null_memory()
      xc_h_pt     = null_memory()
      xc_dh_pt    = null_memory()
      wt_pt       = null_memory()
      wt2_pt      = null_memory()
      gwt_pt      = null_memory()
      g2wt_pt     = null_memory()
      drho_pt     = null_memory()
      dgrho_pt    = null_memory()
      dgamma_pt   = null_memory()
      ddrhoa_pt   = null_memory()
      ddrhob_pt   = null_memory()
      ddgrhoa_pt  = null_memory()
      ddgrhob_pt  = null_memory()
      iiwrk1      = null_memory()
      iwrk1       = null_memory()
      iwrk2       = null_memory()
      iwrk3       = null_memory()
      iwrk4       = null_memory()
      iwrk5       = null_memory()
      iwrk6       = null_memory()
      iwrk7       = null_memory()
      iwrk8       = null_memory()
      iwrk9       = null_memory()
      iwrk10      = null_memory()
      iwrk11      = null_memory()
      iwrk12      = null_memory()
c
      gradcorr_sw = CD_gradcorr()
      kinetic_sw  = CD_kinetic()
      gradkin_sw = gradcorr_sw.or.kinetic_sw
      if (kinetic_sw) then
        if (dksm_exp_sw.or.lhs_sw.or.rhs_sw.or.dksm_sw.or.
     &      hess_sw) then
           call caserr("Meta-GGAs not supported yet")
        endif
      endif
      call find_num_grid_centres
C
C     Calculate inter atom distance array
C
_IFN(qmmm)
      call dijcalc
_ENDIF
c
c     build table of basis function radii
c     
      bfn_radii_ct = incr_memory(nao,'d')
      bfn_radii_pt = allocate_memory2(nao,'d',fnm,snm,
     &                                'bfn_radii')
      call calc_bfn_radii(memory_fp(bfn_radii_pt), ao_tag, nao, 
     &     ngridcentres, arad, first_bf, igrid, iout)
      call override_atom_radii(ngridcentres,igrid,arad)
c
      call calc_max_active_atm(screen_sw,weight_scheme,arad,
     &     ngridcentres,imax_active_atm)
      if (.not.geometric_pert_sw) then
         imax_active_atm = max(npert/3,imax_active_atm)
      endif
c
      call calc_max_active_bfn(screen_sw,arad,first_bf,ngridcentres,
     &     memory_fp(bfn_radii_pt),nao,
     &     imax_active_bfn)

C 
C     Build up angular and radial grid points
C 
      call npoints_by_accuracy(accuracy,ngtypes,igrid,nradpt_num,
     +                         ntheta,nphi,nomega)
c
c     This is used for null typed atoms (e.g. bqs)
c
      call build_radgrid(ngtypes,igrid,nradpt_num,
     +                   prpt,prwt,extwr_sw,iout)
c
c     Set the number of derivatives needed from the weights
c
      call set_weight_derivative_level(idwght,gradwght_sw,
     &                                 grad_sw,dksm_exp_sw,hess_sw)
c     
c     Scale accumulators as we will be summing at the end 
c     
      if(ks_sw)then
         ltri = ((BL_basis_size(1)+1)*BL_basis_size(1))/2
      endif
c
c     accumulators for stats - use float to avoid overflows
c
c     find maximum number of radial grid points
c
      nradtot = 0
      nradmx  = 0
      do latm=1,ngridcentres
         atom_num = ian(latm)
         atmt     = gtype_num(latm)
         if(atmt .gt. 0)then
            if(nradmx .lt. nradpt_num(atmt))then
               nradmx=nradpt_num(atmt)
               itmx=atmt
            endif
            nradtot = nradtot +  nradpt_num(atmt)
         endif
      enddo
      inmtyp(0)=0
      do latm = 1, ngtypes
         inmtyp(latm)=inmtyp(latm-1)+nradpt_num(latm)
      enddo
c
c     find maximum number of angular grid points
c
      mxang = 0
      do latm=1,ngridcentres
         gridt = gtype_num(latm)
         if (gridt.ne.0) then
            if (ang_grid_scheme(gridt,igrid).eq.AG_LEG) then
               do iradzone=1,radzones_num(gridt,igrid)
                  nctht=min(ntheta(gridt),
     +                  thetpt_radzn_num(iradzone,gridt,igrid))
                  ncphi=min(nphi(gridt),
     +                  phipt_radzn_num(iradzone,gridt,igrid))
                  mxang=max(mxang,ncphi*nctht)
               enddo
            endif
            if (ang_grid_scheme(gridt,igrid).eq.AG_LEB) then
               do iradzone=1,radzones_num(gridt,igrid)
                  mxang=max(mxang,
     &                  angpt_radzn_num(iradzone,gridt,igrid))
               enddo
            endif
         endif
      enddo
c     
c     == initialise dynamic load balancing
c     
      ieshll = incr_memory(inmtyp(ngtypes),'d')
      irshll = incr_memory(inmtyp(ngtypes),'d')

      if (kinetic_sw) then
         tau_pt = incr_memory2(mxp0*2,'d',fnm,snm,'tau')
      endif
c
c     Setup memory for screening tables
c
      active_bfn_list_pt = incr_memory2(imax_active_bfn,'i',
     +                     fnm,snm,'bfn_list')
      active_bfn_indx_pt = incr_memory2(imax_active_bfn,'i',
     +                     fnm,snm,'bfn_indx')
      active_wgh_atms_pt = incr_memory2(imax_active_atm,'i',
     +                     fnm,snm,'wgh_atms')
      active_bfn_atms_pt = incr_memory2(imax_active_atm,'i',
     +                     fnm,snm,'bfn_atms')
      if (lhs_sw.or.rhs_sw.or.dksm_sw) then
         active_chf_pert_pt = incr_memory2(3*imax_active_atm,'i',
     +                        fnm,snm,'chf_pert')
      endif
c
c     Setup memory for the weights
c
      wt_pt = incr_memory(mxp0,'d')
      wt2_pt = incr_memory(mxp0,'d')
      if (idwght.ge.1) then
         gwt_pt = incr_memory(mxp0*3*imax_active_atm,'d')
      endif
c
c     do latm=1,ngridcentres
c
c        Begin while loop over lrad
c
c        This while loop over radial shells also performs the
c        radial screening. This construct was needed to keep the basis
c        determined atom types and the grid types independent.
c
c10      if (gridt.ne.0.and.
c    +       lrad.le.nradpt_num(gridt).and.
c    +       (.not.screen_sw.or.
c    +        prpt(gridt,lrad).le.arad(latm)) ) then

               if (sort_points_sw.and.mxang.gt.mxp0) then
                  apts_pt = incr_memory2(3*mxang,'d',fnm,snm,
     +                              'apts_tmp')
                  awpt_pt = incr_memory2(mxang,'d',fnm,snm,
     +                              'awpt_tmp')
                  iprm_pt = incr_memory2(mxang,'i',fnm,snm,
     +                              'iprm_tmp')
c                       call lebedevlaikov
c                       call group_points
                  call decr_memory2(iprm_pt,'i',fnm,snm,
     +                                    'iprm_tmp')
                  call decr_memory2(awpt_pt,'d',fnm,snm,
     +                                    'awpt_tmp')
                  call decr_memory2(apts_pt,'d',fnm,snm,
     +                                    'apts_tmp')
               else
c                       call lebedevlaikov
               endif
c     
c              construct list of neighbouring atoms for use in weight 
c              computations
c     
c              call bld_active_wght_atm(max_active_atm,active_wgh_atms,
c    +              n_active_wgh_atm,latm,rpt,arad,rnear,ngridcentres,
c    +              screen_sw,weight_scheme,iout)
c     
c              construct list of neighbouring atoms for use in basis 
c              functions computations
c     
c              call bld_active_bfn_atm(max_active_atm,active_bfn_atms,
c    +              n_active_bfn_atm,latm,rpt,arad,ngridcentres,
c    +              screen_sw,iout)
c
c              construct list of perturbations for which the 
c              perturbing atom is close enough to the current batch
c              of grid points so that this batch will have significant
c              contributions
c
               if (lhs_sw.or.rhs_sw.or.dksm_sw) then
c                 call bld_active_prt_atm
               endif

c              do lang0=1,nang,mxp0
c                 lo = lang0
c                 hi = min( (lo+mxp0-1),nang)
c                 npts=hi-lo+1
c                 ntot = ntot + npts
                  mxp = min(mxp0,mxang)
c                 call bld_batch(mxp,natoms,latm,lo,hi,n_active_wgh_atm,
c    &                           active_wgh_atms,rpt,rwt,apts,atom_c,
c    &                           awpt,ra2_comp,ra2_val,memory_fp(wt_pt))
C    
C                 compute weights
c    
                  call memreq_calc_weights(igrid,mxp,npts,natoms,
     &                 ngridcentres,latm,0,imax_active_atm,
     &                 memory_fp,memory_int)
c     
c                 construct basis function values and derivatives
c     
                  if (hess_sw) then
                     idbfn=2
                     if(gradcorr_sw)idbfn=3
                  else if (grad_sw.or.dksm_exp_sw) then
                     idbfn=1
                     if(gradcorr_sw)idbfn=2
                  else
                     idbfn=0
                     if(gradcorr_sw)idbfn=1
                  endif
                  bfn_val_pt = incr_memory(mxp*nao,'d')
                  if (idbfn.ge.1) then
                     bfng_val_pt = incr_memory(mxp*nao*3,'d')
                  endif
                  if (idbfn.ge.2) then
                     bfn_hess_pt = incr_memory(mxp*nao*6,'d')
                  endif
                  if (idbfn.ge.3) then
                     bfn_3rd_pt = incr_memory(mxp*nao*10,'d')
                  endif

c                 call bas_val

c
c                 Plan how the various components will be evaluated
c
                  call route_batch(nvec,naocc,nbocc,nao,imax_active_bfn,
     &                 ao_in_sw,mo_in_sw,rks_sw,gradcorr_sw,
     &                 e_sw,ks_sw,grad_sw,dksm_exp_sw,
     &                 rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &                 eval_mo_sw,idmo,na_mo,nb_mo,den_mo_sw,
     &                 dksm_exp_mo_sw,rhs_mo_sw,lhs_mo_sw,dksm_mo_sw,
     &                 hess_mo_sw)
c
c                 Calculate the MOs
c
                  if (eval_mo_sw) then
                     amo_val_pt = incr_memory(mxp*na_mo,'d')
                     if (idmo.ge.1) then
                        amo_grad_pt = incr_memory(3*mxp*na_mo,'d')
                     endif
                     if (.not.rks_sw) then
                        bmo_val_pt = incr_memory(mxp*nb_mo,'d')
                        if (idmo.ge.1) then
                          bmo_grad_pt = incr_memory(3*mxp*nb_mo,'d')
                        endif
                     endif
                     if (screen_sw) then
c                       call calc_mo_val_scr
                     else
c                       call calc_mo_val
                     endif
                  endif
C
C                 Calculate rho and grad rho
C
                  idfun = 0
                  if (ks_sw.or.grad_sw) then
                     idfun = 1
                  endif
                  if (dksm_exp_sw.or.lhs_sw.or.rhs_sw.or.dksm_sw.or.
     &                hess_sw) then
                     idfun = 2
                  endif
                  xc_e_pt = incr_memory(mxp,'d')
                  if (idfun.ge.1) then
                     xc_v_pt = incr_memory(2*mxp,'d')
                     if (gradcorr_sw) then
                        xc_dv_pt = incr_memory(6*mxp,'d')
                     endif
                     if (kinetic_sw) then
                        xc_dt_pt = incr_memory2(2*mxp,'d',fnm,
     &                             snm,'xc_dt')
                     endif
                  endif
                  if (idfun.ge.2) then
                     xc_h_pt = incr_memory(3*mxp,'d')
                     if (gradcorr_sw) then
                        xc_dh_pt = incr_memory(12*mxp,'d')
                     endif
                  endif
c                 if(.not.screen_sw .or. rhomax.gt.rhotol(igrid) )then
C     
C                    Evaluate functional at grid points
C     
                     if (gradcorr_sw) then
                        gamma_pt = incr_memory(3*mxp,'d')
c                       call calc_gamma
                     endif
                     call memreq_xcfunc(idfun,rks_sw,
     &                    npts,mxp,
c    &                    rho,memory_fp(gamma_pt),
c    &                    memory_fp(xc_e_pt),
c    &                    memory_fp(xc_v_pt),memory_fp(xc_dv_pt),
c    &                    memory_fp(xc_h_pt),memory_fp(xc_dh_pt),
     &                    memory_fp)
                     if (gradcorr_sw) then
                        call decr_memory(gamma_pt,'d')
                     endif
C     
C                    compute weights
c     
c                    rshell - distance of this shell from the home atom
c                    since shells are not compacted yet, just take 1st 
c                    point
c     
                     if (idwght.gt.0) then
                        call memreq_calc_weights(igrid,mxp,npts,natoms,
     &                       ngridcentres,latm,idwght,imax_active_atm,
     &                       memory_fp,memory_int)
                     endif

c                    xc_e = 0.0d0
                     if(e_sw) then
c                      xc_e = ddot(npts,memory_fp(xc_e_pt),1,
c    +                                  memory_fp(wt_pt),1)
c           
c                      memory_fp(ieshll+inmtyp(gridt-1)+lrad-1) = xc_e +
c    +                    memory_fp(ieshll+inmtyp(gridt-1)+lrad-1)
                     endif

                     if (ks_sw) then
C     
C                       Add contribution to Kohn-Sham matrices
C     
                        if(screen_sw)then
                           iwrk1 = incr_memory(mxp*imax_active_bfn,'d')
                           iwrk2 = incr_memory(mxp*3,'d')
                           iwrk3 = incr_memory(mxp,'d')
c                          call kmaddcs_scr
                           call decr_memory(iwrk3,'d')
                           call decr_memory(iwrk2,'d')
                           call decr_memory(iwrk1,'d')
                        else
                           iwrk1 = incr_memory(nao,'d')
                           iwrk2 = incr_memory(mxp*3,'d')
c                          call kmaddcs
                           call decr_memory(iwrk2,'d')
                           call decr_memory(iwrk1,'d')
                        endif
                     endif
c
                     if (grad_sw) then
C
C                       Form the gradients for each atom centre
C
                        iwrk1 = incr_memory(3*nao,'d')
                        iwrk2 = incr_memory(nao,'d')
                        iwrk3 = incr_memory(mxp*3,'d')
                        if(screen_sw) then
c                          call xc_forces_scr
                        else
c                          call xc_forces
                        endif
                        call decr_memory(iwrk3,'d')
                        call decr_memory(iwrk2,'d')
                        call decr_memory(iwrk1,'d')
                     endif
c
                     if (rhs_sw) then
                        if (rhs_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp*3,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp*3,'d')
                                 iwrk5 = incr_memory(mxp*nvec,'d')
                                 iwrk6 = incr_memory(mxp*nvec,'d')
                                 iwrk7 = incr_memory(mxp,'d')
c                                call rks_rhs_mo_gga_scr
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp,'d')
                                 iwrk5 = incr_memory(mxp,'d')
c                                call rks_rhs_mo_scr
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp*3,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp*3,'d')
                                 iwrk5 = incr_memory(mxp*nvec,'d')
                                 iwrk6 = incr_memory(mxp*nvec,'d')
                                 iwrk7 = incr_memory(mxp,'d')
c                                call rks_rhs_mo_gga
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp,'d')
                                 iwrk5 = incr_memory(mxp,'d')
c                                call rks_rhs_mo
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           endif
                        else
                           drho_pt = incr_memory(mxp*2*npert,'d')
                           if (gradcorr_sw) then
                              dgrho_pt  = incr_memory(mxp*6*npert,'d')
                              dgamma_pt = incr_memory(mxp*3*npert,'d')
                           endif
                           if (screen_sw) then
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao_scr
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
c                             rhs: we can use lhs_dft_scr here
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao_scr
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           else
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
c                             rhs: we can use lhs_dft here
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           endif
                           if (gradcorr_sw) then
                              call decr_memory(dgamma_pt,'d')
                              call decr_memory(dgrho_pt,'d')
                           endif
                           call decr_memory(drho_pt,'d')
                        endif
                     endif
c
                     if (lhs_sw) then
                        if (lhs_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp*3,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp*3,'d')
                                 iwrk5 = incr_memory(mxp*nvec,'d')
                                 iwrk6 = incr_memory(mxp*nvec,'d')
                                 iwrk7 = incr_memory(mxp,'d')
c                                call rks_lhs_mo_gga_scr
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp,'d')
                                 iwrk5 = incr_memory(mxp,'d')
c                                call rks_lhs_mo_scr
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp*3,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp*3,'d')
                                 iwrk5 = incr_memory(mxp*nvec,'d')
                                 iwrk6 = incr_memory(mxp*nvec,'d')
                                 iwrk7 = incr_memory(mxp,'d')
c                                call rks_lhs_mo_gga
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1 = incr_memory(mxp,'d')
                                 iwrk2 = incr_memory(mxp,'d')
                                 iwrk3 = incr_memory(mxp,'d')
                                 iwrk4 = incr_memory(mxp,'d')
                                 iwrk5 = incr_memory(mxp,'d')
c                                call rks_lhs_mo
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           endif 
                        else
                           drho_pt = incr_memory(mxp*2*npert,'d')
                           if (gradcorr_sw) then
                              dgrho_pt  = incr_memory(mxp*6*npert,'d')
                              dgamma_pt = incr_memory(mxp*3*npert,'d')
                           endif
                           if (screen_sw) then
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao_scr
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao_scr
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           else
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           endif
                           if (gradcorr_sw) then
                              call decr_memory(dgamma_pt,'d')
                              call decr_memory(dgrho_pt,'d')
                           endif
                           call decr_memory(drho_pt,'d')
                        endif
                     endif
c
                     if (dksm_sw) then
                        if (dksm_mo_sw) then
                           drho_pt = incr_memory(mxp*2*npert,'d')
                           if (gradcorr_sw) then
                              dgrho_pt  = incr_memory(mxp*6*npert,'d')
                              dgamma_pt = incr_memory(mxp*3*npert,'d')
                           endif
                           if (screen_sw) then
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_dksm_mo_scr
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
c                             call dksm_dft_mo_scr
                           else
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_dksm_mo
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
c                             call dksm_dft_mo
                           endif
                           if (gradcorr_sw) then
                              call decr_memory(dgamma_pt,'d')
                              call decr_memory(dgrho_pt,'d')
                           endif
                           call decr_memory(drho_pt,'d')
                        else
                           drho_pt = incr_memory(mxp*2*npert,'d')
                           if (gradcorr_sw) then
                              dgrho_pt  = incr_memory(mxp*6*npert,'d')
                              dgamma_pt = incr_memory(mxp*3*npert,'d')
                           endif
                           if (screen_sw) then
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao_scr
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao_scr
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           else
                              iwrk1 = incr_memory(mxp,'d')
c                             call den_pert_ao
                              call decr_memory(iwrk1,'d')
                              if (gradcorr_sw) then
c                                call gamma_pert
                              endif
                              iwrk1 = incr_memory(mxp,'d')
                              iwrk2 = incr_memory(mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*nao*npert,'d')
                              iwrk4 = incr_memory(mxp*nao*nao,'d')
                              iwrk5 = incr_memory(mxp*nao*npert,'d')
                              iwrk6 = incr_memory(mxp*nao,'d')
                              iwrk7 = incr_memory(mxp*npert,'d')
                              iwrk8 = incr_memory(mxp*nao*npert,'d')
c                             call lhs_dft_ao
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           endif
                        endif
                     endif
c
                     if (dksm_exp_sw) then
                        if (dksm_exp_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1  = incr_memory(mxp*npert,'d')
                                 iwrk2  = incr_memory(mxp*npert*3,'d')
                                 iwrk3  = incr_memory(mxp*npert,'d')
                                 iwrk4  = incr_memory(mxp*npert*3,'d')
                                 iwrk5  = incr_memory(mxp*nvec,'d')
                                 iwrk6  = incr_memory(mxp,'d')
                                 iwrk7  = incr_memory(mxp,'d')
                                 iwrk8  = incr_memory(mxp,'d')
                                 iwrk9  = incr_memory(mxp,'d')
                                 iwrk10 = incr_memory(mxp,'d')
                                 iwrk11 = incr_memory(mxp*nvec,'d')
                                 if (gradwght_sw) then
c                                   call rks_dksm_exp_mo_gga_gwt_scr
                                 else
c                                   call rks_dksm_exp_mo_gga_scr
                                 endif
                                 call decr_memory(iwrk11,'d')
                                 call decr_memory(iwrk10,'d')
                                 call decr_memory(iwrk9,'d')
                                 call decr_memory(iwrk8,'d')
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1  = incr_memory(mxp*npert,'d')
                                 iwrk2  = incr_memory(mxp*npert,'d')
                                 iwrk3  = incr_memory(mxp*nvec,'d')
                                 if (gradwght_sw) then
c                                   call rks_dksm_exp_mo_gwt_scr
                                 else
c                                   call rks_dksm_exp_mo_scr
                                 endif
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1  = incr_memory(mxp*npert,'d')
                                 iwrk2  = incr_memory(mxp*npert*3,'d')
                                 iwrk3  = incr_memory(mxp*npert,'d')
                                 iwrk4  = incr_memory(mxp*npert*3,'d')
                                 iwrk5  = incr_memory(mxp*nvec,'d')
                                 iwrk6  = incr_memory(mxp,'d')
                                 iwrk7  = incr_memory(mxp,'d')
                                 iwrk8  = incr_memory(mxp,'d')
                                 iwrk9  = incr_memory(mxp,'d')
                                 iwrk10 = incr_memory(mxp,'d')
                                 iwrk11 = incr_memory(mxp*nvec,'d')
                                 if (gradwght_sw) then
c                                   call rks_dksm_exp_mo_gga_gwt
                                 else
c                                   call rks_dksm_exp_mo_gga
                                 endif
                                 call decr_memory(iwrk11,'d')
                                 call decr_memory(iwrk10,'d')
                                 call decr_memory(iwrk9,'d')
                                 call decr_memory(iwrk8,'d')
                                 call decr_memory(iwrk7,'d')
                                 call decr_memory(iwrk6,'d')
                                 call decr_memory(iwrk5,'d')
                                 call decr_memory(iwrk4,'d')
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              else
                                 iwrk1 = incr_memory(mxp*npert,'d')
                                 iwrk2 = incr_memory(mxp*npert,'d')
                                 iwrk3 = incr_memory(mxp*nvec,'d')
                                 if (gradwght_sw) then
c                                   call rks_dksm_exp_mo_gwt
                                 else
c                                   call rks_dksm_exp_mo
                                 endif
                                 call decr_memory(iwrk3,'d')
                                 call decr_memory(iwrk2,'d')
                                 call decr_memory(iwrk1,'d')
                              endif
                           endif
                        else
                           drho_pt = incr_memory(
     &                               mxp*6*imax_active_atm,'d')
                           if (gradcorr_sw) then
                              dgrho_pt = incr_memory(
     &                                   mxp*18*imax_active_atm,'d')
                           endif
                           if (screen_sw) then
c                             call den_pert_exp_ao_scr
c                             call dksm_exp_dft_ao_scr
                           else
c                             call den_pert_exp_ao
c                             call dksm_exp_dft_ao
                           endif
                           if (gradcorr_sw) then
                              call decr_memory(dgrho_pt,'d')
                           endif
                           call decr_memory(drho_pt,'d')
                        endif
                     endif
c
                     if (hess_sw) then
                        if (hess_mo_sw) then
                           if (screen_sw) then
                              iwrk1 = incr_memory(mxp*npert,'d')
                              iwrk2 = incr_memory(3*mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*npert,'d')
                              iwrk4 = incr_memory(9*mxp,'d')
                              iwrk5 = incr_memory(3*9*mxp,'d')
                              iwrk6 = incr_memory(3*mxp,'d')
                              iwrk7 = incr_memory(3*mxp,'d')
                              iwrk8 = incr_memory(6*mxp,'d')
                              iwrk9 = incr_memory(6*mxp,'d')
                              iwrk10= incr_memory(6*mxp,'d')
                              iwrk11= incr_memory(10*mxp,'d')
                              iwrk12= incr_memory(mxp,'d')
                              iiwrk1= incr_memory(
     &                                imax_active_atm+1,'i')
c                             call rks_hess_dft_mo_scr
                              call decr_memory(iiwrk1,'i')
                              call decr_memory(iwrk12,'d')
                              call decr_memory(iwrk11,'d')
                              call decr_memory(iwrk10,'d')
                              call decr_memory(iwrk9,'d')
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           else
                              iwrk1 = incr_memory(mxp*npert,'d')
                              iwrk2 = incr_memory(3*mxp*npert,'d')
                              iwrk3 = incr_memory(mxp*npert,'d')
                              iwrk4 = incr_memory(9*mxp,'d')
                              iwrk5 = incr_memory(3*9*mxp,'d')
                              iwrk6 = incr_memory(3*mxp,'d')
                              iwrk7 = incr_memory(3*mxp,'d')
                              iwrk8 = incr_memory(6*mxp,'d')
                              iwrk9 = incr_memory(6*mxp,'d')
                              iwrk10= incr_memory(6*mxp,'d')
                              iwrk11= incr_memory(10*mxp,'d')
                              iwrk12= incr_memory(mxp,'d')
c                             call rks_hess_dft_mo
                              call decr_memory(iwrk12,'d')
                              call decr_memory(iwrk11,'d')
                              call decr_memory(iwrk10,'d')
                              call decr_memory(iwrk9,'d')
                              call decr_memory(iwrk8,'d')
                              call decr_memory(iwrk7,'d')
                              call decr_memory(iwrk6,'d')
                              call decr_memory(iwrk5,'d')
                              call decr_memory(iwrk4,'d')
                              call decr_memory(iwrk3,'d')
                              call decr_memory(iwrk2,'d')
                              call decr_memory(iwrk1,'d')
                           endif
                        else
                           nprt = npert
                           if (screen_sw) nprt = 3*imax_active_atm
                           drho_pt  = incr_memory(mxp*2*npert,'d')
                           ddrhoa_pt = incr_memory(mxp*9,'d')
                           if (.not.rks_sw) then
                              ddrhob_pt = incr_memory(mxp*9,'d')
                           endif
                           if (gradcorr_sw) then
                              dgrho_pt  = incr_memory(mxp*6*npert,'d')
                              ddgrhoa_pt = incr_memory(mxp*3*9,'d')
                              if (.not.rks_sw) then
                                 ddgrhob_pt = incr_memory(mxp*3*9,'d')
                              endif
                           endif
                           if (screen_sw) then
c                             call den_pert_exp_ao_scr
                              iiwrk1 = incr_memory(
     &                                 imax_active_atm+1,'i')
c                             call hess_dft_ao_scr
                              call decr_memory(iiwrk1,'i')
                           else
c                             call den_pert_exp_ao
c                             call hess_dft_ao
                           endif
                           if (gradcorr_sw) then
                              if (.not.rks_sw) then
                                 call decr_memory(ddgrhob_pt,'d')
                              endif
                              call decr_memory(ddgrhoa_pt,'d')
                              call decr_memory(dgrho_pt,'d')
                           endif
                           if (.not.rks_sw) then
                              call decr_memory(ddrhob_pt,'d')
                           endif
                           call decr_memory(ddrhoa_pt,'d')
                           call decr_memory(drho_pt,'d')
                        endif
                     endif
C     *
C     *Perform various sums for output
C     *
                  if (idfun.ge.2) then
                     if (gradcorr_sw) then
                        call decr_memory(xc_dh_pt,'d')
                     endif
                     call decr_memory(xc_h_pt,'d')
                  endif
                  if (idfun.ge.1) then
                     if (kinetic_sw) then
                        call decr_memory2(xc_dt_pt,'d',fnm,snm,
     &                                    'xc_dt')
                     endif
                     if (gradcorr_sw) then
                        call decr_memory(xc_dv_pt,'d')
                     endif
                     call decr_memory(xc_v_pt,'d')
                  endif
                  call decr_memory(xc_e_pt,'d')

                  if (eval_mo_sw) then
                     if (.not.rks_sw) then
                        if (idmo.ge.1) then
                           call decr_memory(bmo_grad_pt,'d')
                        endif
                        call decr_memory(bmo_val_pt,'d')
                     endif
                     if (idmo.ge.1) then
                        call decr_memory(amo_grad_pt,'d')
                     endif
                     call decr_memory(amo_val_pt,'d')
                  endif
                  if (idbfn.ge.3) then
                     call decr_memory(bfn_3rd_pt,'d')
                  endif
                  if (idbfn.ge.2) then
                     call decr_memory(bfn_hess_pt,'d')
                  endif
                  if (idbfn.ge.1) then
                     call decr_memory(bfng_val_pt,'d')
                  endif
                  call decr_memory(bfn_val_pt,'d')
c
c              enddo 
c
c     
c              ps  end parallel section
c     
c           endif

c           lrad=lrad+1
c           goto 10
c        endif
c        end of while loop over lrad

c     enddo
c
c     Clearup memory for the weights
c
      if (idwght.ge.1) then
         call decr_memory(gwt_pt,'d')
      endif
      call decr_memory(wt2_pt,'d')
      call decr_memory(wt_pt,'d')
c
c     Clearup memory for screening tables
c
      if (lhs_sw.or.rhs_sw.or.dksm_sw) then
         call decr_memory2(active_chf_pert_pt,'i',fnm,snm,'chf_pert')
      endif
      call decr_memory2(active_bfn_atms_pt,'i',fnm,snm,'bfn_atms')
      call decr_memory2(active_wgh_atms_pt,'i',fnm,snm,'wgh_atms')
      call decr_memory2(active_bfn_indx_pt,'i',fnm,snm,'bfn_indx')
      call decr_memory2(active_bfn_list_pt,'i',fnm,snm,'bfn_list')
c
c     global sum of accumulated quantities
c
      idum = memreq_pg_dgop(ngridcentres,'+')
      idum = memreq_pg_dgop(ngridcentres,'+')
      idum = memreq_pg_dgop(ngridcentres,'+')
      idum = memreq_pg_igop(ngridcentres,'+')
      if(ks_sw)then
         idum = memreq_pg_dgop(ltri,'+')
         if(.not. rks_sw)idum = memreq_pg_dgop(ltri,'+')
      endif
      if(grad_sw)then
         idum = memreq_pg_dgop(natoms*3,'+')
      endif
      if(dksm_exp_sw)then
         if(ao_in_sw)then
            idum = memreq_pg_dgop(nao*(nao+1)/2*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nao*(nao+1)/2*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            idum = memreq_pg_dgop(nvec*(nvec+1)/2*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nvec*(nvec+1)/2*npert,'+')
            endif
         endif
      endif
      if(rhs_sw)then
         if(ao_in_sw)then
            idum = memreq_pg_dgop(nao*nao*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            idum = memreq_pg_dgop(naocc*(nvec-naocc)*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nbocc*(nvec-nbocc)*npert,'+')
            endif
         endif
      endif
      if(lhs_sw)then
         if(ao_in_sw)then
            idum = memreq_pg_dgop(nao*nao*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            idum = memreq_pg_dgop(naocc*(nvec-naocc)*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nbocc*(nvec-nbocc)*npert,'+')
            endif
         endif
      endif
      if(dksm_sw)then
         if(ao_in_sw)then
            idum = memreq_pg_dgop(nao*nao*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            idum = memreq_pg_dgop(nvec*(nvec+1)/2*npert,'+')
            if(.not.rks_sw)then
               idum = memreq_pg_dgop(nvec*(nvec+1)/2*npert,'+')
            endif
         endif
      endif
      if(hess_sw)then
         idum = memreq_pg_dgop((natoms*3)**2,'+')
      endif
c
      idum = memreq_pg_dgop(inmtyp(ngtypes),'+')
      idum = memreq_pg_dgop(inmtyp(ngtypes),'+')
c     
c     global sum of accumulated quantities
c     
      if (kinetic_sw) then
         call decr_memory2(tau_pt,'d',fnm,snm,'tau')
      endif
      call decr_memory(irshll,'d')
      call decr_memory(ieshll,'d')
      call free_memory2(bfn_radii_pt,'d',fnm,snm,'bfn_radii')
      call decr_memory(bfn_radii_ct,'d')

c     call exitc(0)

      return
      end
c
c---- the routines that do the real work -------------------------------
c     
      subroutine exquad(memory_fp,memory_int,
     &     igrid,
     &     ao_tag,nao,nvec,naocc,nbocc,chf_pert_atms,npert,
     &     geometric_pert_sw,
     &     apts,awpt,prpt,prwt,
     &     avec,bvec,adens,bdens,
     &     sa_ao,sb_ao,sa_mo,sb_mo,
     &     ua_ao,ub_ao,ua_mo,ub_mo,
     &     da_ao,db_ao,da_mo,db_mo,
     &     kma,kmb,grad,
     &     ba_ao,bb_ao,ba_mo,bb_mo,
     &     ga_ao,gb_ao,ga_mo,gb_mo,
     &     fa_ao,fb_ao,fa_mo,fb_mo,
     &     fxa_ao,fxb_ao,fxa_mo,fxb_mo,
     &     hess,
     &     ra2_val,ra2_comp,
     &     rho,grho,
     &     e_sw,ks_sw,grad_sw,
     &     dksm_exp_sw,rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &     mxp0,ao_in_sw,mo_in_sw,extwr_sw,accuracy,iout,smat)

      implicit none
c
c     Parameters  
c
INCLUDE(common/dft_parameters)
c
c     In variables  
c
INCLUDE(common/dft_api)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_xc)
INCLUDE(common/dft_dij)
      integer igrid             ! which grid to use
      integer ao_tag
      integer nao               ! the number of AOs
      integer nvec              ! the number of MOs
      integer naocc             ! the number of occupied alpha orbitals
      integer nbocc             ! the number of occupied beta  orbitals
      integer npert             ! the number of perturbations
      integer chf_pert_atms(npert) ! in the case of geometric 
                                   ! perturbations this lists the atom
                                   ! associated with a particular
                                   ! perturbation.
      logical geometric_pert_sw ! true if geometric perturbations are
                                ! considered, otherwise the source of
                                ! the perturbation is an external field.
      REAL avec(nao,nvec)       ! the alpha MOs
      REAL bvec(nao,nvec)       ! the beta  MOs
      REAL adens(nao*(nao+1)/2) ! the alpha density matrix in AO-basis
      REAL bdens(nao*(nao+1)/2) ! the beta  density matrix in AO-basis
      REAL sa_ao(nao,nao,npert) ! the derivative of the occupied-
                                ! occupied block of the overlap matrix
      REAL sb_ao(nao,nao,npert) ! the derivative of the occupied-
                                ! occupied block of the overlap matrix
      REAL sa_mo(nvec*(nvec+1)/2,npert)
      REAL sb_mo(nvec*(nvec+1)/2,npert)
      REAL ua_ao(nao,nao,npert) ! the alpha u matrix (Au=b)
      REAL ub_ao(nao,nao,npert) ! the beta  u matrix (Au=b)
      REAL ua_mo(naocc,nvec-naocc,npert) ! the alpha u matrix (Au=b)
      REAL ub_mo(nbocc,nvec-nbocc,npert) ! the beta  u matrix (Au=b)
      REAL da_ao(nao,nao,npert) ! the alpha perturbed density matrix 
      REAL db_ao(nao,nao,npert) ! the beta  perturbed density matrix 
      REAL da_mo(nvec*(nvec+1)/2,npert) 
      REAL db_mo(nvec*(nvec+1)/2,npert)
      REAL accuracy
      logical e_sw              ! calculate the energy
      logical ks_sw             ! calculate Kohn-Sham matrix elements
      logical grad_sw           ! calculate energy gradient contribution
      logical dksm_exp_sw       ! calculate the derivative of the Kohn-
                                ! Sham matrix with respect to the 
                                ! explicit dependencies on the nuclear
                                ! coordinates
      logical lhs_sw            ! calculate CHF "left-hand-side" 
                                ! contributions
      logical rhs_sw            ! calculate CHF "right-hand-side" 
                                ! contributions
      logical dksm_sw           ! calculate the "wavefunction" term to
                                ! the derivative of the Kohn-Sham matrix
      logical hess_sw           ! calculate energy 2nd derivative 
                                ! contribution
      logical ao_in_sw          ! true if adens and bdens present
      logical mo_in_sw          ! true if avec and bvec present
      logical extwr_sw
      integer mxp0
      integer iout
c
c     Out variables  
c
      REAL kma(*),kmb(*)
      REAL grad(3,*)
      REAL ba_ao(nao,nao,npert)          ! the alpha b vector (Au=b)
      REAL bb_ao(nao,nao,npert)          ! the beta  b vector (Au=b)
      REAL ba_mo(naocc,nvec-naocc,npert) ! the alpha b vector (Au=b)
      REAL bb_mo(nbocc,nvec-nbocc,npert) ! the beta  b vector (Au=b)
      REAL ga_ao(nao,nao,npert)          ! the alpha Au matrix (Au=b)
      REAL gb_ao(nao,nao,npert)          ! the beta  Au matrix (Au=b)
      REAL ga_mo(naocc,nvec-naocc,npert) ! the alpha Au matrix (Au=b)
      REAL gb_mo(nbocc,nvec-nbocc,npert) ! the beta  Au matrix (Au=b)
      REAL fa_ao(nao,nao,npert)          ! the alpha perturbed fock 
                                         ! matrix 
      REAL fb_ao(nao,nao,npert)          ! the beta  perturbed fock 
                                         ! matrix 
      REAL fa_mo(nvec*(nvec+1)/2,npert) 
      REAL fb_mo(nvec*(nvec+1)/2,npert)
      REAL fxa_ao(nao*(nao+1)/2,npert)   ! explicit derivative alpha
                                         ! fock matrix
      REAL fxb_ao(nao*(nao+1)/2,npert)   ! explicit derivative beta
                                         ! fock matrix
      REAL fxa_mo(nvec*(nvec+1)/2,npert)
      REAL fxb_mo(nvec*(nvec+1)/2,npert)
      REAL hess(npert,npert)
      REAL smat(*)

C  Work arrays
      REAL memory_fp(*)
      integer memory_int(*)
      REAL apts(3,*),awpt(*)
      REAL prpt(ngtypes,*)
      REAL prwt(ngtypes,*)
      REAL ra2_val(mxp0*natoms*2)
      REAL ra2_comp(mxp0*natoms*3)
      REAL rho(mxp0*2),grho(mxp0*2*3)
      integer bfn_radii_pt
      integer nprt
c
c still to be allocated
c
      REAL arad(max_atom)
      integer first_bf(max_atom+1)
      integer active_bfn_list_pt
      integer active_bfn_indx_pt
      integer active_bfn_atms_pt
      integer active_wgh_atms_pt
      integer active_chf_pert_pt
      integer n_active_bfn
      integer n_active_bfn_atm
      integer n_active_wgh_atm
      integer n_active_chf_prt
      integer tau_pt
      integer xc_e_pt
      integer xc_v_pt, xc_dv_pt, xc_dt_pt
      integer xc_h_pt, xc_dh_pt
      integer wt_pt, wt2_pt, gwt_pt, g2wt_pt
      integer bfn_val_pt, bfng_val_pt, bfn_hess_pt, bfn_3rd_pt
      integer amo_val_pt, amo_grad_pt, amo_hess_pt
      integer bmo_val_pt, bmo_grad_pt, bmo_hess_pt
      integer gamma_pt
      integer drho_pt,  dgrho_pt, dgamma_pt
      integer ddrhoa_pt, ddgrhoa_pt
      integer ddrhob_pt, ddgrhob_pt

      integer imemcount, imemusage, imemestimate
      integer push_memory_count,    pop_memory_count
      integer push_memory_estimate, pop_memory_estimate
cDEBUG
c     integer igmem_get_usage, igmem_max_usage
c     integer icuse, imuse
cDEBUG
c
C     Functions 
c
c     REAL SG1rad, lograd, srad
      logical opg_root, CD_gradcorr, CD_kinetic
      integer ipg_dlbtask, ipg_nodeid, ipg_nnodes
_IF(single)
      REAL sdot
      integer isamax
_ELSEIF(hp700)
      REAL `vec_$ddot'
      integer idamax
_ELSE
      REAL ddot
      integer idamax
_ENDIF
      integer null_memory
      integer allocate_memory
      integer allocate_memory2

C Local variables
      integer latm,lrad,atmt,atom_num,i,j,k
      integer nradpt_num(max_gtype),nang,mxang
      integer hi, lo, lang0
      integer idwght, idbfn, idmo, idden, idfun
      integer na_mo, nb_mo
      integer mxp, mxpmx
      REAL error, errtol
      REAL rpt,rwt
      REAL atom_xce(max_atom),atom_den(max_atom,2)
      integer atom_pts(max_atom)
      REAL xc_e
      REAL del_psi, del_rho, del_wgh
      REAL rshell, rnear
      integer next, ltri
      REAL fac
      integer npts
      integer nbatch, ntot, npack
      integer iiwrk1
      integer iwrk1, iwrk2,  iwrk3,  iwrk4, iwrk5, iwrk6, iwrk7, iwrk8
      integer iwrk9, iwrk10, iwrk11, iwrk12
      integer nradmx,nradtot,itmx
      integer freq
      integer gridt
      integer irad, ig, igmin
      REAL rhomax, rhoa, rhob
      integer imax_active_atm
      integer imax_active_bfn
      integer TP_TMP_EXQUAD, TP_TMP_EXQUAD_INTRO
      integer TP_TMP_EXQUAD_INTEG, TP_TMP_EXQUAD_DGOP
      integer iprm_pt, apts_pt, awpt_pt
c
c local switches
c
      logical gradcorr_sw, kinetic_sw, gradkin_sw
      logical extout_sw
      logical gwt_avail_sw
      logical eval_mo_sw     ! true if MOs will be evaluated
      logical den_mo_sw      ! true if rho will be calculated from MOs
      logical lhs_mo_sw      ! true if lhs will be calculated from MOs
      logical rhs_mo_sw      ! true if rhs will be calculated from MOs
      logical dksm_mo_sw     ! true if the perturbed Kohn-Sham matrix 
                             ! will be calculated from MOs
      logical dksm_exp_mo_sw ! true if the explicit derivative of the
                             ! Kohn-Sham matrix will be calculated from
                             ! MOs
      logical hess_mo_sw     ! true if d2Exc will be calculated from MOs
      logical accuracy_warned_sw
      save accuracy_warned_sw
      REAL accuracy_warn_tol
      REAL accuracy_fatal_tol
      parameter (accuracy_warn_tol=10.0d0)
      parameter (accuracy_fatal_tol=1200.0d0)

      REAL cpu_prev
      REAL cpu1(3)
      integer ieshll, irshll, inmtyp(0:max_gtype)
c     integer atom_num_by_tag
chvd
c     integer ncurrrad(max_gtype)
      integer ntheta(max_gtype), nphi(max_gtype), nomega(max_gtype)
      integer nctht, ncphi, ncang, iradzone
chvd
c     
c     ============ temp timing analysis ==========
c     
      REAL t_now
      REAL t_start, t_then
      REAL t_becke, t_bas, t_den, t_pack
      REAL t_geom, t_xc, t_ks, t_grad
      LOGICAL time_cond

_IF(qsh)
      integer nbatch_qsh, npts_qsh
      REAL exc_qsh, den_qsh
      character*6 fname
      logical ocrap
_ENDIF

      logical odb
      character*4 fnm
      character*6 snm
c     
c     ==============================================

INCLUDE(../m4/common/parallel)
INCLUDE(../m4/common/timeperiods)
INCLUDE(../m4/common/parcntl)
cDEBUG
c     integer batch_counter
c     logical screen_save_sw, big_mem_sw
c     common/gezeik/batch_counter,screen_save_sw,big_mem_sw
      REAL lograd, srad
cDEBUG
      data accuracy_warned_sw/.false./
c
      fnm = 'xc.m'
      snm = 'exquad'
      gradcorr_sw = CD_gradcorr()
      kinetic_sw = CD_kinetic()
      gradkin_sw = gradcorr_sw.or.kinetic_sw
      imemcount = push_memory_count()
c
      if (kinetic_sw) then
        if (dksm_exp_sw.or.lhs_sw.or.rhs_sw.or.dksm_sw.or.
     &      hess_sw) then
           call caserr("Meta-GGAs not supported yet")
        endif
      endif
c
c major debug mode
c
c     odb = .true.
      odb = .false.
cDEBUG
c     screen_save_sw = screen_sw
c     if (.not.lhs_sw) screen_sw = .false.
c     if (hess_sw) screen_sw = .false.
cDEBUG
C     
C     api entry point
c     
      alpha_den = 0.0d0
      beta_den  = 0.0d0
      totden    = 0.0d0
      totpts    = 0
      XC_energy = 0.0d0
c     
c     =======time========
c     call walltime(t_start)
      call gms_cputime(cpu1)
      t_start = cpu1(3)
      t_geom=0.0d0
      t_becke=0.0d0
      t_pack=0.0d0
      t_bas=0.0d0
      t_den=0.0d0
      t_xc=0.0d0
      t_ks=0.0d0
      t_grad=0.0d0
      freq = 99999999
      cpu_prev = -1.0d0
c     =======time========
c
c define distinct timers for energy/gradient cases
c
      if(hess_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADHES
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADHES_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADHES_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADHES_DGOP
      else if(dksm_exp_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADDKSX
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADDKSX_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADDKSX_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADDKSX_DGOP
      else if(lhs_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADLHS
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADLHS_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADLHS_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADLHS_DGOP
      else if(rhs_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADRHS
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADRHS_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADRHS_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADRHS_DGOP
      else if(dksm_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADDKS
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADDKS_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADDKS_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADDKS_DGOP
      else if(grad_sw)then
         TP_TMP_EXQUAD       = TP_DFT_EXQUADF
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUADF_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUADF_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUADF_DGOP
      else
         TP_TMP_EXQUAD       = TP_DFT_EXQUAD
         TP_TMP_EXQUAD_INTRO = TP_DFT_EXQUAD_INTRO
         TP_TMP_EXQUAD_INTEG = TP_DFT_EXQUAD_INTEG
         TP_TMP_EXQUAD_DGOP  = TP_DFT_EXQUAD_DGOP
      endif
c
c     initialise addresses
c
      bfn_val_pt  = null_memory()
      bfng_val_pt = null_memory()
      bfn_hess_pt = null_memory()
      bfn_3rd_pt  = null_memory()
      amo_val_pt  = null_memory()
      bmo_val_pt  = null_memory()
      amo_grad_pt = null_memory()
      bmo_grad_pt = null_memory()
      amo_hess_pt = null_memory()
      bmo_hess_pt = null_memory()
      tau_pt      = null_memory()
      xc_e_pt     = null_memory()
      xc_v_pt     = null_memory()
      xc_dv_pt    = null_memory()
      xc_dt_pt    = null_memory()
      xc_h_pt     = null_memory()
      xc_dh_pt    = null_memory()
      wt_pt       = null_memory()
      wt2_pt      = null_memory()
      gwt_pt      = null_memory()
      g2wt_pt     = null_memory()
      drho_pt     = null_memory()
      dgrho_pt    = null_memory()
      dgamma_pt   = null_memory()
      ddrhoa_pt   = null_memory()
      ddrhob_pt   = null_memory()
      ddgrhoa_pt  = null_memory()
      ddgrhob_pt  = null_memory()
      iiwrk1      = null_memory()
      iwrk1       = null_memory()
      iwrk2       = null_memory()
      iwrk3       = null_memory()
      iwrk4       = null_memory()
      iwrk5       = null_memory()
      iwrk6       = null_memory()
      iwrk7       = null_memory()
      iwrk8       = null_memory()
      iwrk9       = null_memory()
      iwrk10      = null_memory()
      iwrk11      = null_memory()
      iwrk12      = null_memory()

      if(opg_root() .and. print_sw(DEBUG_QUAD) )then

         write(iout,*)
         write(iout,'(a82)')
     &                '___ Debug Quadrature _________'//
     &                '______________________________'//
     &                '_____________________'
         write(iout,*)
         if (hess_sw) then
            write(iout,*)'start 2nd derivative quadrature'
         else if (dksm_exp_sw) then
            write(iout,*)'start explicit derivative KS matrix '//
     &                   'quadrature'
         else if (lhs_sw) then
            write(iout,*)'start left-hand-side quadrature'
         else if (rhs_sw) then
            write(iout,*)'start right-hand-side quadrature'
         else if (dksm_sw) then
            write(iout,*)'start perturbed KS matrix quadrature'
         else if (grad_sw) then
            write(iout,*)'start gradient quadrature'
         else
            write(iout,*)'start energy quadrature'
         endif
         write(iout,*)
         write(iout,*)'main quadrature parameters:'
         write(iout,*)'    max buffer size                 = ',mxp0
         if (screen_sw) then
            write(iout,*)'    tolerance on Pij in density     = ',
     &           dentol(igrid)
            write(iout,*)'    batch tolerance on rho          = ',
     &           rhotol(igrid)
            write(iout,*)'    tolerance on weight             = ',
     &           wghtol(igrid)
         else
            write(iout,*)'    no screening to be used'
         endif
         if (weight_scheme(igrid) .eq. WT_BECKE) then
            write(iout,*)'    use Becke atomic partition function'
         else if (weight_scheme(igrid) .eq. WT_BECKESCR) then
            write(iout,*)'    use Becke atomic partition function'
            write(iout,*)'        restricted to neighbouring atom list'
         else if (weight_scheme(igrid) .eq. WT_SSF) then
            write(iout,*)
     &   '    use Stratmann/Scuseria/Frisch atomic partition function'
         else if (weight_scheme(igrid) .eq. WT_SSFSCR) then
            write(iout,*)
     &   '    use Stratmann/Scuseria/Frisch atomic partition function'
            write(iout,*)'        restricted to neighbouring atom list'
         else if (weight_scheme(igrid) .eq. WT_MHL) then
            write(iout,*)
     &   '    use Murray/Handy/Laming atomic partition function'
         else if (weight_scheme(igrid) .eq. WT_MHLSCR) then
            write(iout,*)
     &   '    use Murray/Handy/Laming atomic partition function'
            write(iout,*)'        restricted to neighbouring atom list'
         else if (weight_scheme(igrid) .eq. WT_MHL4SSFSCR) then
            write(iout,'(a82)')
     &   '    use Murray/Handy/Laming 4 Stratmann/Scuseria/Frisch '//
     &       'atomic partition function'
            write(iout,*)'        restricted to neighbouring atom list'
         else if (weight_scheme(igrid) .eq. WT_MHL8SSFSCR) then
            write(iout,'(a82)')
     &   '    use Murray/Handy/Laming 8 Stratmann/Scuseria/Frisch '//
     &       'atomic partition function'
            write(iout,*)'        restricted to neighbouring atom list'
         endif
      endif

      call start_time_period(TP_TMP_EXQUAD)
      call start_time_period(TP_TMP_EXQUAD_INTRO)
c
      call find_num_grid_centres
C
C     Calculate inter atom distance array
C
_IFN(qmmm)
      call dijcalc
_ENDIF
c
c     build table of basis function radii
c     
      bfn_radii_pt = allocate_memory2(nao,'d',fnm,snm,
     &                                'bfn_radii')
      call calc_bfn_radii(memory_fp(bfn_radii_pt), ao_tag, nao, 
     &     ngridcentres, arad, first_bf, igrid, iout)
      call override_atom_radii(ngridcentres,igrid,arad)
c
      call calc_max_active_atm(screen_sw,weight_scheme,arad,
     &     ngridcentres,imax_active_atm)
      if (.not.geometric_pert_sw) then
         imax_active_atm = max(npert/3,imax_active_atm)
      endif
c
      call calc_max_active_bfn(screen_sw,arad,first_bf,ngridcentres,
     &     memory_fp(bfn_radii_pt),nao,
     &     imax_active_bfn)

      if(opg_root() .and. odb )then
         write(iout,*)'Basis Fn Radii'
         do i= 1,BL_basis_size(1)
            write(iout,'(i3,f8.4)')i, memory_fp(bfn_radii_pt+i-1)
         enddo
         write(iout,*)'Atom Radii'
         do i= 1,natoms
            write(iout,'(i3,f8.4)')i, arad(i)
         enddo
      endif
      extout_sw = extwr_sw
C 
C     Build up angular and radial grid points
C 
      call npoints_by_accuracy(accuracy,ngtypes,igrid,nradpt_num,
     +                         ntheta,nphi,nomega)
c
c     This is used for null typed atoms (e.g. bqs)
c
      call build_radgrid(ngtypes,igrid,nradpt_num,
     +                   prpt,prwt,extwr_sw,iout)
c
c     Set the number of derivatives needed from the weights
c
      call set_weight_derivative_level(idwght,gradwght_sw,
     &                                 grad_sw,dksm_exp_sw,hess_sw)
c     
c     Scale accumulators as we will be summing at the end 
c     
      call scale_accumulators(nao,nvec,naocc,nbocc,natoms,npert,
     &     rks_sw,ao_in_sw,mo_in_sw,ks_sw,grad_sw,dksm_exp_sw,
     &     lhs_sw,rhs_sw,dksm_sw,hess_sw,kma,kmb,grad,
     &     fxa_ao,fxb_ao,fxa_mo,fxb_mo,ba_ao,bb_ao,ba_mo,bb_mo,
     &     ga_ao,gb_ao,ga_mo,gb_mo,fa_ao,fb_ao,fa_mo,fb_mo,
     &     hess)
c
c     accumulators for stats - use float to avoid overflows
c
      del_psi = 0.0d0
      del_rho = 0.0d0
      del_wgh = 0.0d0
      nbatch = 0
      ntot = 0
      npack= 0
c
c     find maximum number of radial grid points
c
      nradtot = 0
      nradmx  = 0
      do latm=1,ngridcentres
         atom_num = ian(latm)
         atmt     = gtype_num(latm)
         if(atmt .gt. 0)then
            if(nradmx .lt. nradpt_num(atmt))then
               nradmx=nradpt_num(atmt)
               itmx=atmt
            endif
            nradtot = nradtot +  nradpt_num(atmt)
         endif
      enddo
      inmtyp(0)=0
      do latm = 1, ngtypes
         inmtyp(latm)=inmtyp(latm-1)+nradpt_num(latm)
      enddo
c
c     find maximum number of angular grid points
c
      mxang = 0
      do latm=1,ngridcentres
         gridt = gtype_num(latm)
         if (gridt.ne.0) then
            if (ang_grid_scheme(gridt,igrid).eq.AG_LEG) then
               do iradzone=1,radzones_num(gridt,igrid)
                  nctht=min(ntheta(gridt),
     +                  thetpt_radzn_num(iradzone,gridt,igrid))
                  ncphi=min(nphi(gridt),
     +                  phipt_radzn_num(iradzone,gridt,igrid))
                  mxang=max(mxang,ncphi*nctht)
               enddo
            endif
            if (ang_grid_scheme(gridt,igrid).eq.AG_LEB) then
               do iradzone=1,radzones_num(gridt,igrid)
                  mxang=max(mxang,
     &                  angpt_radzn_num(iradzone,gridt,igrid))
               enddo
            endif
         endif
      enddo

      if(opg_root().and.odb)then
         write(iout,*)'total radial shells tot,max:',nradtot,nradmx
         write(iout,*)'ntchnk',ntchnk
         write(iout,*)'chunk size:',nradtot/(ntchnk*ipg_nnodes())
      endif

c     
c     == initialise dynamic load balancing
c     
c     estimate total number of tasks
c     aim for about 40 calls to nxtval
c     set by ntchnk
c

c     write(iout,*)'xc:',nradtot,ntchnk,ipg_nnodes()

      call pg_dlbchunk(max(1,nradtot/(ntchnk*ipg_nnodes())),
     +                 print_sw(DEBUG_QUAD))
      call pg_dlbreset
      next = ipg_dlbtask()

      ieshll = allocate_memory2(inmtyp(ngtypes),'d',fnm,snm,
     +                          'ieshll')
      irshll = allocate_memory2(inmtyp(ngtypes),'d',fnm,snm,
     +                          'irshll')
      do lrad=1,inmtyp(ngtypes)
         memory_fp(ieshll+lrad-1) = 0.0d0
         memory_fp(irshll+lrad-1) = 0.0d0
      enddo

      if (kinetic_sw) then
         tau_pt = allocate_memory2(mxp0*2,'d',fnm,snm,'tau')
      endif

_IF(qsh)
      fname="qsh"
      call intwrt(fname,3,ipg_nodeid(),3,ocrap)
      open(unit=59,file=fname,form='formatted')
_ENDIF
c
      call end_time_period(TP_TMP_EXQUAD_INTRO)
      call start_time_period(TP_TMP_EXQUAD_INTEG)
c
c     Setup memory for screening tables
c
      active_bfn_list_pt = allocate_memory2(imax_active_bfn,'i',
     +                     fnm,snm,'bfn_list')
      active_bfn_indx_pt = allocate_memory2(imax_active_bfn,'i',
     +                     fnm,snm,'bfn_indx')
      active_wgh_atms_pt = allocate_memory2(imax_active_atm,'i',
     +                     fnm,snm,'wgh_atms')
      active_bfn_atms_pt = allocate_memory2(imax_active_atm,'i',
     +                     fnm,snm,'bfn_atms')
      if (lhs_sw.or.rhs_sw.or.dksm_sw) then
         active_chf_pert_pt = allocate_memory2(3*imax_active_atm,'i',
     +                        fnm,snm,'chf_pert')
      endif
c
c     Setup memory for the weights
c
      wt_pt = allocate_memory2(mxp0,'d',fnm,snm,'wt')
      wt2_pt = allocate_memory2(mxp0,'d',fnm,snm,'wt2')
      if (idwght.ge.1) then
         gwt_pt = allocate_memory2(mxp0*3*imax_active_atm,
     +                             'd',fnm,snm,'gwt')
      endif
c
      do latm=1,ngridcentres

         atom_den(latm,1)=0.0d0
         atom_den(latm,2)=0.0d0
         atom_xce(latm)=0.0d0
         atom_pts(latm)=0

         nang = -1
         atom_num = ian(latm)
         gridt = gtype_num(latm)

         lrad=1
         iradzone=1
c
c        Begin while loop over lrad
c
c        This while loop over radial shells also performs the
c        radial screening. This construct was needed to keep the basis
c        determined atom types and the grid types independent.
c
 10      if (gridt.ne.0.and.
     +       lrad.le.nradpt_num(gridt).and.
     +       (.not.screen_sw.or.
     +        prpt(gridt,lrad).le.arad(latm)) ) then

            icount_dlb = icount_dlb + 1

            if(icount_dlb . eq. next) then

               rpt=prpt(gridt,lrad)
               rwt=prwt(gridt,lrad)
c
c              Set up the angular grid as function of the radius
c
c              First we have to make sure that we know which zone
c              the current angular shell is in.
c
 20            if (iradzone.lt.radzones_num(gridt,igrid).and.
     +             bnd_radzn(iradzone,gridt,igrid).lt.rpt) then
                   iradzone=iradzone+1
                   goto 20
               endif
c
c              Second we need to find the required grid size and 
c              create the grids as appropriate.
c
               if (ang_grid_scheme(gridt,igrid).eq.AG_LEG) then
                  nctht = min(ntheta(gridt),
     +                        thetpt_radzn_num(iradzone,gridt,igrid))
                  ncphi = min(nphi(gridt),
     +                        phipt_radzn_num(iradzone,gridt,igrid))
                  if (nang.ne.ncphi*nctht) then
                     call glegend(nctht,ncphi,nang,apts,awpt)
                  endif
               else if (ang_grid_scheme(gridt,igrid).eq.AG_LEB) then
                  ncang = min(nomega(gridt),
     +                        angpt_radzn_num(iradzone,gridt,igrid))
                  if (nang.ne.ncang) then
                     if (sort_points_sw.and.mxang.gt.mxp0) then
                        apts_pt = allocate_memory2(3*mxang,'d',fnm,snm,
     +                              'apts_tmp')
                        awpt_pt = allocate_memory2(mxang,'d',fnm,snm,
     +                              'awpt_tmp')
                        iprm_pt = allocate_memory2(mxang,'i',fnm,snm,
     +                              'iprm_tmp')
                        if (ncang.gt.mxp0) then
                           call lebedevlaikov(ncang,memory_fp(apts_pt),
     +                                              memory_fp(awpt_pt))
                           call group_points(memory_fp(apts_pt),
     +                                       memory_fp(awpt_pt),
     +                                       memory_int(iprm_pt),
     +                                       apts,awpt,ncang)
                        else
                           call lebedevlaikov(ncang,apts,awpt)
                        endif
                        call free_memory2(iprm_pt,'i',fnm,snm,
     +                                    'iprm_tmp')
                        call free_memory2(awpt_pt,'d',fnm,snm,
     +                                    'awpt_tmp')
                        call free_memory2(apts_pt,'d',fnm,snm,
     +                                    'apts_tmp')
                     else
                        call lebedevlaikov(ncang,apts,awpt)
                     endif
                     nang=ncang
                  endif
               endif
c     
c              construct list of neighbouring atoms for use in weight 
c              computations
c     
               call bld_active_wght_atm(imax_active_atm,
     +              memory_int(active_wgh_atms_pt),
     +              n_active_wgh_atm,latm,rpt,arad,rnear,ngridcentres,
     +              screen_sw,weight_scheme,iout)
c     
c              construct list of neighbouring atoms for use in basis 
c              functions computations
c     
               call bld_active_bfn_atm(imax_active_atm,
     +              memory_int(active_bfn_atms_pt),
     +              n_active_bfn_atm,latm,rpt,arad,ngridcentres,
     +              screen_sw,iout)
c
c              construct list of perturbations for which the 
c              perturbing atom is close enough to the current batch
c              of grid points so that this batch will have significant
c              contributions
c
               if (lhs_sw.or.rhs_sw.or.dksm_sw) then
                  call bld_active_prt_atm(geometric_pert_sw,screen_sw,
     +                 iout,chf_pert_atms,npert,
     +                 arad,latm,natoms,
     +                 memory_int(active_chf_pert_pt),
     +                 imax_active_atm,n_active_chf_prt)
               endif

_IF(qsh)
c
cqsh           quadrature shell data for debugging
c                
               nbatch_qsh = 0
               npts_qsh   = 0
               exc_qsh    = 0.0d0
               den_qsh    = 0.0d0
_ENDIF
               do lang0=1,nang,mxp0
_IF(qsh)
                  nbatch_qsh = nbatch_qsh + 1
_ENDIF
                  lo    = lang0
                  hi    = min( (lo+mxp0-1),nang )
                  npts  = hi-lo+1
                  ntot  = ntot + npts
                  mxp   = npts
                  mxpmx = min(mxp0,mxang)
cDEBUG
c           if (hess_sw) then
c           write(*,*)'*** batch no:',latm*1000+lrad*10+lang0/mxp0+1
c           batch_counter = latm*1000+lrad*10+lang0/mxp0+1
c           endif
cDEBUG
c
c==== time===
                  time_cond = (opg_root() .and. 
     &                        (mod(nbatch,freq).eq.0
     &                        .or. print_sw(DEBUG_QUAD) ))
                  if (time_cond) then
                    call gms_cputime(cpu1)
                    t_then = cpu1(3)
                  endif
c==== time===
                  call bld_batch(mxp,natoms,ngridcentres,latm,lo,hi,
     &                           n_active_wgh_atm,
     &                           memory_int(active_wgh_atms_pt),
     &                           rpt,rwt,apts,atom_c,
     &                           awpt,ra2_comp,ra2_val,
     &                           memory_fp(wt2_pt))
c==== time===
                  if (time_cond) then
                    call gms_cputime(cpu1)
                    t_now = cpu1(3)
                    t_geom=t_geom+(t_now-t_then)
                    t_then=t_now
                  endif
c==== time===
C     
C                 compute weights
c     
c                 rshell - distance of this shell from the home atom
c                 since shells are not compacted yet, just take 1st 
c                 point
c     
                  call aclear_dp(memory_fp(wt_pt),npts,1.0d0)
                  call calc_weights(igrid,mxp,npts,natoms,ngridcentres,
     &                 latm,0,n_active_wgh_atm,
     &                 memory_int(active_wgh_atms_pt),rpt,
     &                 rnear,ra2_val,ra2_comp,arad,
     &                 memory_fp(wt_pt),memory_fp(gwt_pt),gwt_avail_sw,
     &                 memory_fp(xc_e_pt),hess,
     &                 memory_fp,memory_int)

c==== time===
                  if (time_cond) then
                    call gms_cputime(cpu1)
                    t_now = cpu1(3)
                    t_becke=t_becke+(t_now-t_then)
                    t_then=t_now
                  endif
c==== time===
c     
c                 now pack arrays over points to exclude points with 
c                 very small weights
c     
                  if(screen_sw)then
                     call screen_weight(npts,npack,ra2_val,ra2_comp,
     &                    memory_fp(wt_pt),memory_fp(wt2_pt),
     &                    wghtol(igrid),natoms,ngridcentres,
     &                    n_active_wgh_atm,
     &                    memory_int(active_wgh_atms_pt),mxp)

c==== time===
                     if (time_cond) then
                       call gms_cputime(cpu1)
                       t_now = cpu1(3)
                       t_pack=t_pack+(t_now-t_then)
                       t_then=t_now
                     endif
c==== time===
                  else
                     do i=0,npts-1
                        memory_fp(wt_pt+i)=memory_fp(wt_pt+i)
     &                                    *memory_fp(wt2_pt+i)
                     enddo
                  endif
c     

_IF(qsh)
                  npts_qsh = npts_qsh + npts
_ENDIF
                  if (npts.le.0) then
                     del_wgh = del_wgh + 1.0d0
                  else
c     
c                 construct basis function values and derivatives
c     
                  if (hess_sw) then
                     idbfn=2
                     if(gradkin_sw)idbfn=3
                  else if (grad_sw.or.dksm_exp_sw) then
                     idbfn=1
                     if(gradkin_sw)idbfn=2
                  else
                     idbfn=0
                     if(gradkin_sw)idbfn=1
                  endif
                  bfn_val_pt = allocate_memory2(mxpmx*nao,'d',fnm,
     +                            snm,'bfn_val')
                  if (idbfn.ge.1) then
                     bfng_val_pt = allocate_memory2(mxpmx*nao*3,'d',
     +                                fnm,snm,'bfng_val')
                  endif
                  if (idbfn.ge.2) then
                     bfn_hess_pt = allocate_memory2(mxpmx*nao*6,'d',
     +                                fnm,snm,'bfn_hess')
                  endif
                  if (idbfn.ge.3) then
                     bfn_3rd_pt = allocate_memory2(mxpmx*nao*10,'d',
     +                               fnm,snm,'bfn_3rd')
                  endif

                  call bas_val(ao_tag,nao,
     &                 ra2_comp,ra2_val,
     &                 memory_fp(bfn_val_pt),memory_fp(bfng_val_pt),
     &                 memory_fp(bfn_hess_pt),memory_fp(bfn_3rd_pt),
     &                 npts,mxp,screen_sw,first_bf,
     &                 memory_int(active_bfn_atms_pt),n_active_bfn_atm,
     &                 memory_int(active_bfn_list_pt),
     &                 memory_int(active_bfn_indx_pt),
     &                 imax_active_bfn,n_active_bfn,
     &                 psitol, del_psi, memory_fp(bfn_radii_pt), 
     &                 idbfn, odb, iout)

c                 if (lhs_sw.or.rhs_sw.or.dksm_sw) then
c                    call bld_lhs_atms(lhs_atms,npert/3,
c    &                                 active_bfn_atms,n_active_bfn_atm,
c    &                                 lhs_atm_list,n_lhs_atm)
c                 endif

                  if( odb )then
                     write(iout,*)'number of near functions',
     &                    n_active_bfn
                     write(iout,*)'fn map', 
     &                    (memory_int(active_bfn_list_pt+i),i=0
     &                    ,n_active_bfn-1)
                     do i=1,npts
                        write(iout,*)'bas vals'
                        write(iout,10101)(memory_fp(bfn_val_pt+
     &                        mxp*(i-1)+j),j=0,n_active_bfn-1)
                     enddo
10101                format(1x,5e20.14)
                  endif

c==== time===
                  if (time_cond) then
                    call gms_cputime(cpu1)
                    t_now = cpu1(3)
                    t_bas=t_bas+(t_now-t_then)
                    t_then=t_now
                  endif
c==== time===
c
c                 Plan how the various components will be evaluated
c
                  call route_batch(nvec,naocc,nbocc,nao,n_active_bfn,
     &                 ao_in_sw,mo_in_sw,rks_sw,gradkin_sw,
     &                 e_sw,ks_sw,grad_sw,dksm_exp_sw,
     &                 rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &                 eval_mo_sw,idmo,na_mo,nb_mo,den_mo_sw,
     &                 dksm_exp_mo_sw,rhs_mo_sw,lhs_mo_sw,dksm_mo_sw,
     &                 hess_mo_sw)
c
c                 Calculate the MOs
c
                  if (eval_mo_sw) then
                     amo_val_pt = allocate_memory2(mxpmx*na_mo,'d',fnm,
     &                            snm,'amo_val')
                     if (idmo.ge.1) then
                        amo_grad_pt = allocate_memory2(3*mxpmx*na_mo,
     &                                'd',fnm,snm,'amo_grad')
                     endif
                     if (.not.rks_sw) then
                        bmo_val_pt = allocate_memory2(mxpmx*nb_mo,'d',
     &                               fnm,snm,'bmo_val')
                        if (idmo.ge.1) then
                          bmo_grad_pt = allocate_memory2(3*mxpmx*nb_mo,
     &                                  'd',fnm,snm,'bmo_grad')
                        endif
                     endif
                     if (screen_sw) then
                        call calc_mo_val_scr(rks_sw,idmo,npts,mxp,
     &                       nao,nvec,na_mo,nb_mo,
     &                       memory_int(active_bfn_list_pt),
     &                       n_active_bfn,
     &                       avec,bvec,
     &                       memory_fp(bfn_val_pt),
     &                       memory_fp(bfng_val_pt),
     &                       memory_fp(bfn_hess_pt),
     &                       memory_fp(amo_val_pt),
     &                       memory_fp(amo_grad_pt),
     &                       memory_fp(amo_hess_pt),
     &                       memory_fp(bmo_val_pt),
     &                       memory_fp(bmo_grad_pt),
     &                       memory_fp(bmo_hess_pt))
                     else
                        call calc_mo_val(rks_sw,idmo,npts,mxp,nao,
     &                       nvec,na_mo,nb_mo,avec,bvec,
     &                       memory_fp(bfn_val_pt),
     &                       memory_fp(bfng_val_pt),
     &                       memory_fp(bfn_hess_pt),
     &                       memory_fp(amo_val_pt),
     &                       memory_fp(amo_grad_pt),
     &                       memory_fp(amo_hess_pt),
     &                       memory_fp(bmo_val_pt),
     &                       memory_fp(bmo_grad_pt),
     &                       memory_fp(bmo_hess_pt))
                     endif
                  endif
C
C                 Calculate rho and grad rho
C
                  idden = 0
                  if (gradcorr_sw) idden = 1

                  if (den_mo_sw) then
                     call den_val_mo(rks_sw,idden,npts,mxp,
     &                    na_mo,nb_mo,naocc,nbocc,
     &                    memory_fp(amo_val_pt),
     &                    memory_fp(amo_grad_pt),
     &                    memory_fp(bmo_val_pt),
     &                    memory_fp(bmo_grad_pt),
     &                    rho,grho,rhotol(igrid),screen_sw)
                  else if (screen_sw .and. n_active_bfn.ne.nao) then
                     call den_val_ao_scr(rks_sw,.false.,adens,bdens,
     &                    ao_tag,nao,
     &                    memory_fp(bfn_val_pt),memory_fp(bfng_val_pt),
     &                    memory_fp(bfn_hess_pt),rho,grho,
     &                    npts,mxp,screen_sw,
     &                    memory_int(active_bfn_list_pt),
     &                    n_active_bfn,dentol(igrid),
     &                    rhotol(igrid),idden)
                  else
                     call den_val_ao(rks_sw,.false.,adens,bdens,
     &                    ao_tag,nao,
     &                    memory_fp(bfn_val_pt),memory_fp(bfng_val_pt),
     &                    memory_fp(bfn_hess_pt),rho,grho,
     &                    npts,mxp,screen_sw,
     &                    memory_int(active_bfn_list_pt),
     &                    n_active_bfn,idden)
                  endif
c
c                 Calculate tau
c
                  if (kinetic_sw) then
                     if (den_mo_sw) then
                        call tau_val_mo(rks_sw,npts,mxp,na_mo,nb_mo,
     &                       naocc,nbocc,memory_fp(amo_grad_pt),
     &                       memory_fp(bmo_grad_pt),
     &                       memory_fp(tau_pt),rho,screen_sw)
                     else if (screen_sw .and. n_active_bfn.ne.nao) then
                        call tau_val_ao_scr(rks_sw,adens,bdens,
     &                       ao_tag,nao,memory_fp(bfng_val_pt),rho,
     &                       memory_fp(tau_pt),npts,mxp,screen_sw,
     &                       memory_int(active_bfn_list_pt),
     &                       n_active_bfn,dentol)
                     else
                        call tau_val_ao(rks_sw,adens,bdens,ao_tag,nao,
     &                       memory_fp(bfng_val_pt),rho,
     &                       memory_fp(tau_pt),
     &                       npts,mxp,screen_sw,
     &                       memory_int(active_bfn_list_pt),
     &                       n_active_bfn)
                     endif
                  endif

                  if( odb )then
                     write(iout,*)'density ',nbatch,npts
                     do i=1,npts
                        write(iout,1101)rho(i),(grho(i+(j-1)*2*mxp),
     &                                          j=1,3)
 1101                   format(1x,4f16.8)
                     enddo
                  endif
c
c==== time===
                  if (time_cond) then
                    call gms_cputime(cpu1)
                    t_now = cpu1(3)
                    t_den=t_den+(t_now-t_then)
                    t_then=t_now
                  endif
c==== time===
c     
c                 check there is some density here
c     
                  rhomax = rho(idamax(npts,rho,1))

                  idfun = 0
                  if (ks_sw.or.grad_sw) then
                     idfun = 1
                  endif
                  if (dksm_exp_sw.or.lhs_sw.or.rhs_sw.or.dksm_sw.or.
     &                hess_sw) then
                     idfun = 2
                  endif
                  xc_e_pt = allocate_memory2(mxpmx,'d',fnm,snm,
     &                      'xc_e')
                  if (idfun.ge.1) then
                     xc_v_pt = allocate_memory2(2*mxpmx,'d',fnm,
     &                         snm,'xc_v')
                     if (gradcorr_sw) then
                        xc_dv_pt = allocate_memory2(6*mxpmx,'d',fnm,
     &                             snm,'xc_dv')
                     endif
                     if (kinetic_sw) then
                        xc_dt_pt = allocate_memory2(2*mxpmx,'d',fnm,
     &                             snm,'xc_dt')
                     endif
                  endif
                  if (idfun.ge.2) then
                     xc_h_pt = allocate_memory2(3*mxpmx,'d',fnm,
     &                         snm,'xc_h')
                     if (gradcorr_sw) then
                        xc_dh_pt = allocate_memory2(12*mxpmx,'d',fnm,
     &                             snm,'xc_dh')
                     endif
                  endif
                  if(.not.screen_sw .or. rhomax.gt.rhotol(igrid) )then
C     
C                    Evaluate functional at grid points
C     
                     if (gradcorr_sw) then
                        gamma_pt = allocate_memory2(3*mxpmx,'d',fnm,
     &                             snm,'gamma')
                        call calc_gamma(mxp,npts,rks_sw,grho,
     &                       memory_fp(gamma_pt))
                     endif
                     call xcfunc(idfun,rks_sw,
     &                    npts,mxp,mxpmx,rho,memory_fp(gamma_pt),
     &                    memory_fp(tau_pt),
     &                    memory_fp(xc_e_pt),
     &                    memory_fp(xc_v_pt),memory_fp(xc_dv_pt),
     &                    memory_fp(xc_dt_pt),
     &                    memory_fp(xc_h_pt),memory_fp(xc_dh_pt),
     &                    memory_fp)
                     if (gradcorr_sw) then
                        call free_memory2(gamma_pt,'d',fnm,snm,
     &                                    'gamma')
                     endif

c==== time===
                     if (time_cond) then
                       call gms_cputime(cpu1)
                       t_now = cpu1(3)
                       t_xc=t_xc+(t_now-t_then)
                       t_then=t_now
                     endif
c==== time===
C     
C                    compute weights
c     
c                    rshell - distance of this shell from the home atom
c                    since shells are not compacted yet, just take 1st 
c                    point
c     
                     if (idwght.gt.0) then
                        call dcopy(npts,memory_fp(wt2_pt),1,
     &                                  memory_fp(wt_pt),1)
                        call calc_weights(igrid,mxp,npts,natoms,
     &                       ngridcentres,latm,idwght,n_active_wgh_atm,
     &                       memory_int(active_wgh_atms_pt),
     &                       rpt,rnear,ra2_val,
     &                       ra2_comp,arad,memory_fp(wt_pt),
     &                       memory_fp(gwt_pt),gwt_avail_sw,
     &                       memory_fp(xc_e_pt),hess,
     &                       memory_fp,memory_int)
                     endif
c==== time===
                     if (time_cond) then
                       call gms_cputime(cpu1)
                       t_now = cpu1(3)
                       t_becke=t_becke+(t_now-t_then)
                       t_then=t_now
                     endif
c==== time===
                     xc_e = 0.0d0
                     if(e_sw) then
                       xc_e = ddot(npts,memory_fp(xc_e_pt),1,
     +                                  memory_fp(wt_pt),1)
                       atom_xce(latm)=atom_xce(latm)+xc_e
            
c                      for debugging purposes:
                       memory_fp(ieshll+inmtyp(gridt-1)+lrad-1) = xc_e +
     +                    memory_fp(ieshll+inmtyp(gridt-1)+lrad-1)
                     endif

                     if (ks_sw) then
C     
C                       Add contribution to Kohn-Sham matrices
C     
                        if(screen_sw)then
                           iwrk1 = allocate_memory2(
     &                             mxpmx*imax_active_bfn,
     &                             'd',fnm,snm,'iwrk1')
                           iwrk2 = allocate_memory2(mxpmx*3,'d',
     &                                              fnm,snm,'iwrk2')
                           iwrk3 = allocate_memory2(mxpmx,'d',
     &                                              fnm,snm,'iwrk3')
                           call kmaddcs_scr(memory_fp(bfn_val_pt),
     &                          memory_fp(bfng_val_pt),memory_fp(wt_pt),
     &                          memory_fp(xc_v_pt),memory_fp(xc_dv_pt),
     &                          memory_fp(xc_dt_pt),
     &                          grho,gradcorr_sw,kinetic_sw,
     &                          kma,kmb,rks_sw,nao,npts,mxp,
     &                          memory_int(active_bfn_list_pt),
     &                          n_active_bfn,
     &                          memory_fp(iwrk1),memory_fp(iwrk2),
     &                          memory_fp(iwrk3),smat)
                           call free_memory2(iwrk3,'d',fnm,snm,'iwrk3')
                           call free_memory2(iwrk2,'d',fnm,snm,'iwrk2')
                           call free_memory2(iwrk1,'d',fnm,snm,'iwrk1')
                        else
                           iwrk1 = allocate_memory2(nao,'d',fnm,snm,
     &                                              'iwrk1')
                           iwrk2 = allocate_memory2(mxpmx*3,'d',fnm,snm,
     &                                              'iwrk2')
                           call kmaddcs(memory_fp(bfn_val_pt),
     &                          memory_fp(bfng_val_pt),memory_fp(wt_pt),
     &                          memory_fp(xc_v_pt),memory_fp(xc_dv_pt),
     &                          memory_fp(xc_dt_pt),
     &                          grho,gradcorr_sw,kinetic_sw,
     &                          kma,kmb,rks_sw,memory_fp(iwrk1),
     &                          memory_fp(iwrk2),nao,npts,mxp,smat)
                           call free_memory2(iwrk2,'d',fnm,snm,'iwrk2')
                           call free_memory2(iwrk1,'d',fnm,snm,'iwrk1')
                        endif
c==== time===
                        if (time_cond) then
                          call gms_cputime(cpu1)
                          t_now = cpu1(3)
                          t_ks=t_ks+(t_now-t_then)
                          t_then=t_now
                        endif
c==== time===
                     endif
c
                     if (grad_sw) then
C
C                       Form the gradients for each atom centre
C
                        iwrk1 = allocate_memory2(3*nao,'d',fnm,snm,
     &                                           'iwrk1') 
                        iwrk2 = allocate_memory2(nao,'d',fnm,snm,
     &                                           'iwrk2')
                        iwrk3 = allocate_memory2(mxpmx*3,'d',fnm,snm,
     &                                           'iwrk3')
                        if(screen_sw) then
                           call xc_forces_scr(ao_tag,nao,rks_sw,
     &                          gradcorr_sw,kinetic_sw,
     &                          adens,bdens,gradwght_sw,
     &                          gwt_avail_sw,
     &                          memory_fp(bfn_val_pt),
     &                          memory_fp(bfng_val_pt),
     &                          memory_fp(bfn_hess_pt),
     &                          memory_fp(wt_pt),memory_fp(gwt_pt),
     &                          memory_fp(xc_e_pt),memory_fp(xc_v_pt),
     &                          memory_fp(xc_dv_pt),memory_fp(xc_dt_pt),
     &                          grho,
     &                          memory_fp(iwrk1),memory_fp(iwrk2),
     &                          memory_fp(iwrk3),
     &                          memory_int(active_wgh_atms_pt), 
     &                          n_active_wgh_atm,
     &                          first_bf, 
     &                          memory_int(active_bfn_list_pt),
     &                          n_active_bfn,
     &                          memory_int(active_bfn_atms_pt), 
     &                          n_active_bfn_atm,
     &                          grad,npts,mxp,latm)
                        else
                           call xc_forces(ao_tag,nao,rks_sw,
     &                          gradcorr_sw,kinetic_sw,
     &                          adens,bdens,gradwght_sw,gwt_avail_sw,
     &                          memory_fp(bfn_val_pt),
     &                          memory_fp(bfng_val_pt),
     &                          memory_fp(bfn_hess_pt),
     &                          memory_fp(wt_pt),memory_fp(gwt_pt),
     &                          memory_fp(xc_e_pt),memory_fp(xc_v_pt),
     &                          memory_fp(xc_dv_pt),memory_fp(xc_dt_pt),
     &                          grho,
     &                          memory_fp(iwrk1),memory_fp(iwrk2),
     &                          memory_fp(iwrk3),
     &                          memory_int(active_wgh_atms_pt), 
     &                          n_active_wgh_atm,
     &                          grad,npts,mxp,latm)
                        endif
                        call free_memory2(iwrk3,'d',fnm,snm,'iwrk3')
                        call free_memory2(iwrk2,'d',fnm,snm,'iwrk2')
                        call free_memory2(iwrk1,'d',fnm,snm,'iwrk1')
c==== time===
                        if (time_cond) then
                          call gms_cputime(cpu1)
                          t_now = cpu1(3)
                          t_grad=t_grad+(t_now-t_then)
                          t_then=t_now
                        endif
c==== time===
                     endif
c
                     if (rhs_sw) then
                        if (rhs_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk5')
                                 iwrk6 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk6')
                                 iwrk7 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk7')
                                 call rks_rhs_mo_gga_scr(mxp,npts,npert,
     &                                nvec,naocc,
     &                                memory_int(active_chf_pert_pt),
     &                                n_active_chf_prt,
     &                                dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_dv_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(xc_dh_pt),sa_mo,
     &                                memory_fp(amo_val_pt),
     &                                memory_fp(amo_grad_pt),
     &                                grho,
     &                                memory_fp(iwrk1),memory_fp(iwrk2),
     &                                memory_fp(iwrk3),memory_fp(iwrk4),
     &                                memory_fp(iwrk5),memory_fp(iwrk6),
     &                                memory_fp(iwrk7),ba_mo)
                                 call free_memory2(iwrk7,'d',fnm,snm,
     &                                'iwrk7')
                                 call free_memory2(iwrk6,'d',fnm,snm,
     &                                'iwrk6')
                                 call free_memory2(iwrk5,'d',fnm,snm,
     &                                'iwrk5')
                                 call free_memory2(iwrk4,'d',fnm,snm,
     &                                'iwrk4')
                                 call free_memory2(iwrk3,'d',fnm,snm,
     &                                'iwrk3')
                                 call free_memory2(iwrk2,'d',fnm,snm,
     &                                'iwrk2')
                                 call free_memory2(iwrk1,'d',fnm,snm,
     &                                'iwrk1')
                              else
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk5')
                                 call rks_rhs_mo_scr(npts,nvec,naocc,
     &                                npert,
     &                                memory_int(active_chf_pert_pt),
     &                                n_active_chf_prt,dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(amo_val_pt),
     &                                sa_mo,ba_mo,memory_fp(iwrk1),
     &                                memory_fp(iwrk2),memory_fp(iwrk3),
     &                                memory_fp(iwrk4),memory_fp(iwrk5))
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk5')
                                 iwrk6 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk6')
                                 iwrk7 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk7')
                                 call rks_rhs_mo_gga(mxpmx,npts,npert,
     &                                nvec,naocc,dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_dv_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(xc_dh_pt),sa_mo,
     &                                memory_fp(amo_val_pt),
     &                                memory_fp(amo_grad_pt),
     &                                grho,
     &                                memory_fp(iwrk1),memory_fp(iwrk2),
     &                                memory_fp(iwrk3),memory_fp(iwrk4),
     &                                memory_fp(iwrk5),memory_fp(iwrk6),
     &                                memory_fp(iwrk7),ba_mo)
                                 call free_memory2(iwrk7,'d',
     &                                fnm,snm,'iwrk7')
                                 call free_memory2(iwrk6,'d',
     &                                fnm,snm,'iwrk6')
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              else
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk5')
                                 call rks_rhs_mo(npts,nvec,naocc,npert,
     &                                dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(amo_val_pt),
     &                                sa_mo,ba_mo,memory_fp(iwrk1),
     &                                memory_fp(iwrk2),memory_fp(iwrk3),
     &                                memory_fp(iwrk4),memory_fp(iwrk5))
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              endif
                           endif
                        else
                           drho_pt = allocate_memory2(mxpmx*2*npert,'d',
     &                               fnm,snm,'drho')
                           if (gradcorr_sw) then
                              dgrho_pt  = allocate_memory2(
     &                                    mxpmx*6*npert,
     &                                    'd',fnm,snm,'dgrho')
                              dgamma_pt = allocate_memory2(
     &                                    mxpmx*3*npert,
     &                                    'd',fnm,snm,'dgamma')
                           endif
                           if (screen_sw) then
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_ao_scr(rks_sw,sa_ao,sb_ao,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn, 
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),memory_fp(iwrk1),
     &                             npts,mxp,idden,dentol)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
c                             rhs: we can use lhs_dft_scr here
                              iwrk1 = allocate_memory2(mxpmx,'d',
     &                                fnm,snm,'iwrk1')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk2')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk3')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'iwrk4')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk5')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'iwrk6')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk7')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk8')
                              call lhs_dft_ao_scr(rks_sw,gradcorr_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npts,npert,
     &                             memory_fp(wt_pt),grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             ba_ao,bb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                             'iwrk8')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                             'iwrk7')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                             'iwrk6')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                             'iwrk5')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                             'iwrk4')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                             'iwrk3')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                             'iwrk2')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                           else
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_ao(rks_sw,sa_ao,sb_ao,nao,
     &                             npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(iwrk1),npts,mxp,idden)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
c                             rhs: we can use lhs_dft here
                              iwrk1 = allocate_memory2(mxpmx,'d',
     &                                fnm,snm,'iwrk1')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk2')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk3')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'iwrk4')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk5')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'iwrk6')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk7')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk8')
                              call lhs_dft_ao(rks_sw,gradcorr_sw,nao,
     &                             npts,npert,
     &                             memory_fp(wt_pt),grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             ba_ao,bb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                             'iwrk8')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                             'iwrk7')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                             'iwrk6')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                             'iwrk5')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                             'iwrk4')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                             'iwrk3')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                             'iwrk2')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                           endif
                           if (gradcorr_sw) then
                              call free_memory2(dgamma_pt,'d',fnm,snm,
     &                             'dgamma')
                              call free_memory2(dgrho_pt,'d',fnm,snm,
     &                             'dgrho')
                           endif
                           call free_memory2(drho_pt,'d',fnm,snm,
     &                             'drho')
                        endif
                     endif
c
                     if (lhs_sw) then
                        if (lhs_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk5')
                                 iwrk6 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk6')
                                 iwrk7 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk7')
                                 call rks_lhs_mo_gga_scr(mxp,npts,npert,
     &                                nvec,naocc,
     &                                memory_int(active_chf_pert_pt),
     &                                n_active_chf_prt,dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_dv_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(xc_dh_pt),ua_mo,
     &                                memory_fp(amo_val_pt),
     &                                memory_fp(amo_grad_pt),
     &                                grho,
     &                                memory_fp(iwrk1),memory_fp(iwrk2),
     &                                memory_fp(iwrk3),memory_fp(iwrk4),
     &                                memory_fp(iwrk5),memory_fp(iwrk6),
     &                                memory_fp(iwrk7),ga_mo)
                                 call free_memory2(iwrk7,'d',fnm,snm,
     &                                'iwrk7')
                                 call free_memory2(iwrk6,'d',fnm,snm,
     &                                'iwrk6')
                                 call free_memory2(iwrk5,'d',fnm,snm,
     &                                'iwrk5')
                                 call free_memory2(iwrk4,'d',fnm,snm,
     &                                'iwrk4')
                                 call free_memory2(iwrk3,'d',fnm,snm,
     &                                'iwrk3')
                                 call free_memory2(iwrk2,'d',fnm,snm,
     &                                'iwrk2')
                                 call free_memory2(iwrk1,'d',fnm,snm,
     &                                'iwrk1')
                              else
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk5')
                                 call rks_lhs_mo_scr(npts,nvec,naocc,
     &                                npert,
     &                                memory_int(active_chf_pert_pt),
     &                                n_active_chf_prt,dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(amo_val_pt),
     &                                ua_mo,ga_mo,memory_fp(iwrk1),
     &                                memory_fp(iwrk2),memory_fp(iwrk3),
     &                                memory_fp(iwrk4),memory_fp(iwrk5))
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx*3,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk5')
                                 iwrk6 = allocate_memory2(mxpmx*nvec,
     &                                   'd',fnm,snm,'iwrk6')
                                 iwrk7 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk7')
                                 call rks_lhs_mo_gga(mxp,npts,npert,
     &                                nvec,naocc,dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_dv_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(xc_dh_pt),ua_mo,
     &                                memory_fp(amo_val_pt),
     &                                memory_fp(amo_grad_pt),
     &                                grho,
     &                                memory_fp(iwrk1),memory_fp(iwrk2),
     &                                memory_fp(iwrk3),memory_fp(iwrk4),
     &                                memory_fp(iwrk5),memory_fp(iwrk6),
     &                                memory_fp(iwrk7),ga_mo)
                                 call free_memory2(iwrk7,'d',
     &                                fnm,snm,'iwrk7')
                                 call free_memory2(iwrk6,'d',
     &                                fnm,snm,'iwrk6')
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              else
                                 iwrk1 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk1')
                                 iwrk2 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk2')
                                 iwrk3 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk3')
                                 iwrk4 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk4')
                                 iwrk5 = allocate_memory2(mxpmx,'d',
     &                                   fnm,snm,'iwrk5')
                                 call rks_lhs_mo(npts,nvec,naocc,npert,
     &                                dentol,
     &                                memory_fp(wt_pt),
     &                                memory_fp(xc_h_pt),
     &                                memory_fp(amo_val_pt),
     &                                ua_mo,ga_mo,memory_fp(iwrk1),
     &                                memory_fp(iwrk2),memory_fp(iwrk3),
     &                                memory_fp(iwrk4),memory_fp(iwrk5))
                                 call free_memory2(iwrk5,'d',
     &                                fnm,snm,'iwrk5')
                                 call free_memory2(iwrk4,'d',
     &                                fnm,snm,'iwrk4')
                                 call free_memory2(iwrk3,'d',
     &                                fnm,snm,'iwrk3')
                                 call free_memory2(iwrk2,'d',
     &                                fnm,snm,'iwrk2')
                                 call free_memory2(iwrk1,'d',
     &                                fnm,snm,'iwrk1')
                              endif
                           endif 
                        else
                           drho_pt = allocate_memory2(mxpmx*2*npert,'d',
     &                               fnm,snm,'drho')
                           if (gradcorr_sw) then
                              dgrho_pt  = allocate_memory2(
     &                                    mxpmx*6*npert,
     &                                    'd',fnm,snm,'dgrho')
                              dgamma_pt = allocate_memory2(
     &                                    mxpmx*3*npert,
     &                                    'd',fnm,snm,'dgamma')
                           endif
                           if (screen_sw) then
                              iwrk1 = allocate_memory2(mxpmx,'d',
     &                                fnm,snm,'iwrk1')
                              call den_pert_ao_scr(rks_sw,ua_ao,ub_ao,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),memory_fp(iwrk1),
     &                             npts,mxp,idden,dentol)
                              call free_memory2(iwrk1,'d',
     &                             fnm,snm,'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              iwrk1 = allocate_memory2(mxpmx,'d',
     &                                fnm,snm,'iwrk1')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk2')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk3')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'iwrk4')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk5')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'iwrk6')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk7')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk8')
                              call lhs_dft_ao_scr(rks_sw,gradcorr_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npts,npert,memory_fp(wt_pt),
     &                             grho,memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             ga_ao,gb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                                         'iwrk8')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                                         'iwrk7')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                                         'iwrk6')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                                         'iwrk5')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                                         'iwrk4')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                                         'iwrk3')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                                         'iwrk2')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                                         'iwrk1')
                           else
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                                 snm,'iwrk1')
                              call den_pert_ao(rks_sw,ua_ao,ub_ao,
     &                             nao,npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(iwrk1),npts,mxp,idden)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                                         'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              iwrk1 = allocate_memory2(mxpmx,'d',
     &                                fnm,snm,'iwrk1')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk2')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk3')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'iwrk4')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk5')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'iwrk6')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'iwrk7')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'iwrk8')
                              call lhs_dft_ao(rks_sw,gradcorr_sw,nao,
     &                             npts,npert,
     &                             memory_fp(wt_pt),grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             ga_ao,gb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                                         'iwrk8')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                                         'iwrk7')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                                         'iwrk6')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                                         'iwrk5')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                                         'iwrk4')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                                         'iwrk3')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                                         'iwrk2')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                                         'iwrk1')
                           endif
                           if (gradcorr_sw) then
                              call free_memory2(dgamma_pt,'d',fnm,snm,
     &                                         'dgamma')
                              call free_memory2(dgrho_pt,'d',fnm,snm,
     &                                         'dgrho')
                           endif
                           call free_memory2(drho_pt,'d',fnm,snm,
     &                                      'drho')
                        endif
                     endif
c
                     if (dksm_sw) then
                        if (dksm_mo_sw) then
                           drho_pt = allocate_memory2(mxpmx*2*npert,'d',
     &                               fnm,snm,'drho_pt')
                           if (gradcorr_sw) then
                              dgrho_pt  = allocate_memory2(
     &                                    mxpmx*6*npert,'d',fnm,snm,
     &                                    'dgrho_pt')
                              dgamma_pt = allocate_memory2(
     &                                    mxpmx*3*npert,'d',fnm,snm,
     &                                    'dgamma_pt')
                           endif
                           if (screen_sw) then
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_dksm_mo_scr(rks_sw,idden,
     &                             npts,mxp,npert,nvec,naocc,nbocc,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             da_mo,db_mo,
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             memory_fp(bmo_val_pt),
     &                             memory_fp(bmo_grad_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(iwrk1),
     &                             dentol)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,
     &                                n_active_chf_prt,grho,
     &                                memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              call dksm_dft_mo_scr(rks_sw,gradcorr_sw,
     &                             npts,mxp,nvec,naocc,nbocc,npert,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             memory_fp(wt_pt),
     &                             grho,memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             memory_fp(bmo_val_pt),
     &                             memory_fp(bmo_grad_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),fa_mo,fb_mo)
                           else
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_dksm_mo(rks_sw,idden,npts,
     &                             mxp,npert,nvec,naocc,nbocc,
     &                             da_mo,db_mo,
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             memory_fp(bmo_val_pt),
     &                             memory_fp(bmo_grad_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),memory_fp(iwrk1),
     &                             0.0d0)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              call dksm_dft_mo(rks_sw,gradcorr_sw,
     &                             npts,mxp,nvec,naocc,nbocc,npert,
     &                             memory_fp(wt_pt),
     &                             grho,memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             memory_fp(bmo_val_pt),
     &                             memory_fp(bmo_grad_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),fa_mo,fb_mo)
                           endif
                           if (gradcorr_sw) then
                              call free_memory2(dgamma_pt,'d',fnm,snm,
     &                                          'dgamma_pt')
                              call free_memory2(dgrho_pt,'d',fnm,snm,
     &                                          'dgrho_pt')
                           endif
                           call free_memory2(drho_pt,'d',fnm,snm,
     &                                       'drho_pt')
                        else
                           drho_pt = allocate_memory2(mxpmx*2*npert,'d',
     &                               fnm,snm,'drho_pt')
                           if (gradcorr_sw) then
                              dgrho_pt  = allocate_memory2(
     &                                    mxpmx*6*npert,'d',fnm,snm,
     &                                    'dgrho_pt')
                              dgamma_pt = allocate_memory2(
     &                                    mxpmx*3*npert,'d',fnm,snm,
     &                                    'dgamma_pt')
                           endif
                           if (screen_sw) then
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_ao_scr(rks_sw,da_ao,db_ao,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),memory_fp(iwrk1),
     &                             npts,mxp,idden,dentol)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'wrrp')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'wrrpd')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'wrrpdn')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'drmn')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'dgrn')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'dvn')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'ddr')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'dhnd')
                              call lhs_dft_ao_scr(rks_sw,gradcorr_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             n_active_bfn,
     &                             memory_int(active_chf_pert_pt),
     &                             n_active_chf_prt,
     &                             nao,npts,npert,
     &                             memory_fp(wt_pt),grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             fa_ao,fb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                             'dhnd')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                             'ddr')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                             'dvn')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                             'dgrn')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                             'drmn')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                             'wrrpdn')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                             'wrrpd')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'wrrp')
                           else
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'iwrk1')
                              call den_pert_ao(rks_sw,da_ao,db_ao,
     &                             nao,npert,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(iwrk1),npts,mxp,idden)
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'iwrk1')
                              if (gradcorr_sw) then
                                 call gamma_pert(rks_sw,npts,npert,
     &                                grho,memory_fp(dgrho_pt),
     &                                memory_fp(dgamma_pt),mxp)
                              endif
                              iwrk1 = allocate_memory2(mxpmx,'d',fnm,
     &                                snm,'wrrp')
                              iwrk2 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'wrrpd')
                              iwrk3 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'wrrpdn')
                              iwrk4 = allocate_memory2(mxpmx*nao*nao,
     &                                'd',fnm,snm,'drmn')
                              iwrk5 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'dgrn')
                              iwrk6 = allocate_memory2(mxpmx*nao,'d',
     &                                fnm,snm,'dvn')
                              iwrk7 = allocate_memory2(mxpmx*npert,'d',
     &                                fnm,snm,'ddr')
                              iwrk8 = allocate_memory2(mxpmx*nao*npert,
     &                                'd',fnm,snm,'dhnd')
                              call lhs_dft_ao(rks_sw,gradcorr_sw,nao,
     &                             npts,npert,
     &                             memory_fp(wt_pt),
     &                             grho,memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(dgamma_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             fa_ao,fb_ao,mxp,
     &                             memory_fp(iwrk1),memory_fp(iwrk2),
     &                             memory_fp(iwrk3),memory_fp(iwrk4),
     &                             memory_fp(iwrk5),memory_fp(iwrk6),
     &                             memory_fp(iwrk7),memory_fp(iwrk8))
                              call free_memory2(iwrk8,'d',fnm,snm,
     &                             'dhnd')
                              call free_memory2(iwrk7,'d',fnm,snm,
     &                             'ddr')
                              call free_memory2(iwrk6,'d',fnm,snm,
     &                             'dvn')
                              call free_memory2(iwrk5,'d',fnm,snm,
     &                             'dgrn')
                              call free_memory2(iwrk4,'d',fnm,snm,
     &                             'drmn')
                              call free_memory2(iwrk3,'d',fnm,snm,
     &                             'wrrpdn')
                              call free_memory2(iwrk2,'d',fnm,snm,
     &                             'wrrpd')
                              call free_memory2(iwrk1,'d',fnm,snm,
     &                             'wrrp')
                           endif
                        endif
                     endif
c
                     if (dksm_exp_sw) then
                        if (dksm_exp_mo_sw) then
                           if (screen_sw) then
                              if (gradcorr_sw) then
                                 iwrk1  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk2  = allocate_memory(mxpmx*npert*3,
     &                                    'd')
                                 iwrk3  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk4  = allocate_memory(mxpmx*npert*3,
     &                                    'd')
                                 iwrk5  = allocate_memory(mxpmx*nvec,
     &                                    'd')
                                 iwrk6  = allocate_memory(mxpmx,'d')
                                 iwrk7  = allocate_memory(mxpmx,'d')
                                 iwrk8  = allocate_memory(mxpmx,'d')
                                 iwrk9  = allocate_memory(mxpmx,'d')
                                 iwrk10 = allocate_memory(mxpmx,'d')
                                 iwrk11 = allocate_memory(mxpmx*nvec,
     &                                    'd')
c                                iwrk12 = allocate_memory(
c    &                                    mxpmx*naocc*nvec,'d')
                                 if (gradwght_sw) then
                                    call rks_dksm_exp_mo_gga_gwt_scr(
     &                                   npts,nao,nvec,naocc,npert,
     &                                   mxp,
     &                                   memory_int(active_bfn_list_pt),
     &                                   memory_int(active_bfn_indx_pt),
     &                                   memory_int(active_bfn_atms_pt),
     &                                   n_active_bfn,n_active_bfn_atm,
     &                                   n_active_wgh_atm,
     &                                   memory_int(active_wgh_atms_pt),
     &                                   gwt_avail_sw,latm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(gwt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_dv_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(xc_dh_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(amo_grad_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   memory_fp(bfn_hess_pt),
     &                                   avec,grho,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   memory_fp(iwrk4),
     &                                   memory_fp(iwrk5),
     &                                   memory_fp(iwrk6),
     &                                   memory_fp(iwrk7),
     &                                   memory_fp(iwrk8),
     &                                   memory_fp(iwrk9),
     &                                   memory_fp(iwrk10),
     &                                   memory_fp(iwrk11),
     &                                   fxa_mo)
                                 else
                                    call rks_dksm_exp_mo_gga_scr(npts,
     &                                   nao,nvec,naocc,npert,mxp,
     &                                   memory_int(active_bfn_list_pt),
     &                                   memory_int(active_bfn_indx_pt),
     &                                   memory_int(active_bfn_atms_pt),
     &                                   n_active_bfn,n_active_bfn_atm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_dv_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(xc_dh_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(amo_grad_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   memory_fp(bfn_hess_pt),
     &                                   avec,grho,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   memory_fp(iwrk4),
     &                                   memory_fp(iwrk5),
     &                                   memory_fp(iwrk6),
     &                                   memory_fp(iwrk7),
     &                                   memory_fp(iwrk8),
     &                                   memory_fp(iwrk9),
     &                                   memory_fp(iwrk10),
     &                                   memory_fp(iwrk11),
     &                                   fxa_mo)
                                 endif
c                                call free_memory(iwrk12,'d')
                                 call free_memory(iwrk11,'d')
                                 call free_memory(iwrk10,'d')
                                 call free_memory(iwrk9,'d')
                                 call free_memory(iwrk8,'d')
                                 call free_memory(iwrk7,'d')
                                 call free_memory(iwrk6,'d')
                                 call free_memory(iwrk5,'d')
                                 call free_memory(iwrk4,'d')
                                 call free_memory(iwrk3,'d')
                                 call free_memory(iwrk2,'d')
                                 call free_memory(iwrk1,'d')
                              else
                                 iwrk1  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk2  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk3  = allocate_memory(mxpmx*nvec,
     &                                    'd')
                                 if (gradwght_sw) then
                                    call rks_dksm_exp_mo_gwt_scr(
     &                                   npts,nao,nvec,naocc,npert,
     &                                   mxp,
     &                                   memory_int(active_bfn_list_pt),
     &                                   memory_int(active_bfn_indx_pt),
     &                                   memory_int(active_bfn_atms_pt),
     &                                   n_active_bfn,n_active_bfn_atm,
     &                                   n_active_wgh_atm,
     &                                   memory_int(active_wgh_atms_pt),
     &                                   gwt_avail_sw,latm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(gwt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   avec,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   fxa_mo)
                                 else
                                    call rks_dksm_exp_mo_scr(npts,nao,
     &                                   nvec,naocc,npert,
     &                                   mxp,
     &                                   memory_int(active_bfn_list_pt),
     &                                   memory_int(active_bfn_indx_pt),
     &                                   memory_int(active_bfn_atms_pt),
     &                                   n_active_bfn,n_active_bfn_atm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   avec,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   fxa_mo)
                                 endif
                                 call free_memory(iwrk3,'d')
                                 call free_memory(iwrk2,'d')
                                 call free_memory(iwrk1,'d')
                              endif
                           else
                              if (gradcorr_sw) then
                                 iwrk1  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk2  = allocate_memory(mxpmx*npert*3,
     &                                    'd')
                                 iwrk3  = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk4  = allocate_memory(mxpmx*npert*3,
     &                                    'd')
                                 iwrk5  = allocate_memory(mxpmx*nvec,
     &                                    'd')
                                 iwrk6  = allocate_memory(mxpmx,'d')
                                 iwrk7  = allocate_memory(mxpmx,'d')
                                 iwrk8  = allocate_memory(mxpmx,'d')
                                 iwrk9  = allocate_memory(mxpmx,'d')
                                 iwrk10 = allocate_memory(mxpmx,'d')
                                 iwrk11 = allocate_memory(mxpmx*nvec,
     &                                    'd')
c                                iwrk12 = allocate_memory(
c    &                                    mxpmx*naocc*nvec,'d')
                                 if (gradwght_sw) then
                                    call rks_dksm_exp_mo_gga_gwt(
     &                                   npts,nao,nvec,naocc,npert,
     &                                   ngridcentres,mxp,gwt_avail_sw,
     &                                   latm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(gwt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_dv_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(xc_dh_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(amo_grad_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   memory_fp(bfn_hess_pt),
     &                                   avec,grho,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   memory_fp(iwrk4),
     &                                   memory_fp(iwrk5),
     &                                   memory_fp(iwrk6),
     &                                   memory_fp(iwrk7),
     &                                   memory_fp(iwrk8),
     &                                   memory_fp(iwrk9),
     &                                   memory_fp(iwrk10),
     &                                   memory_fp(iwrk11),
     &                                   fxa_mo)
                                 else
                                    call rks_dksm_exp_mo_gga(npts,nao,
     &                                   nvec,naocc,npert,ngridcentres,
     &                                   mxp,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_dv_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(xc_dh_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(amo_grad_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   memory_fp(bfn_hess_pt),
     &                                   avec,grho,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   memory_fp(iwrk4),
     &                                   memory_fp(iwrk5),
     &                                   memory_fp(iwrk6),
     &                                   memory_fp(iwrk7),
     &                                   memory_fp(iwrk8),
     &                                   memory_fp(iwrk9),
     &                                   memory_fp(iwrk10),
     &                                   memory_fp(iwrk11),
     &                                   fxa_mo)
                                 endif
c                                call free_memory(iwrk12,'d')
                                 call free_memory(iwrk11,'d')
                                 call free_memory(iwrk10,'d')
                                 call free_memory(iwrk9,'d')
                                 call free_memory(iwrk8,'d')
                                 call free_memory(iwrk7,'d')
                                 call free_memory(iwrk6,'d')
                                 call free_memory(iwrk5,'d')
                                 call free_memory(iwrk4,'d')
                                 call free_memory(iwrk3,'d')
                                 call free_memory(iwrk2,'d')
                                 call free_memory(iwrk1,'d')
                              else
                                 iwrk1 = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk2 = allocate_memory(mxpmx*npert,
     &                                    'd')
                                 iwrk3 = allocate_memory(mxpmx*nvec,
     &                                    'd')
                                 if (gradwght_sw) then
                                    call rks_dksm_exp_mo_gwt(npts,nao,
     &                                   nvec,naocc,npert,ngridcentres,
     &                                   mxp,gwt_avail_sw,latm,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(gwt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   avec,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   fxa_mo)
                                 else
                                    call rks_dksm_exp_mo(npts,nao,nvec,
     &                                   naocc,npert,ngridcentres,mxp,
     &                                   memory_fp(wt_pt),
     &                                   memory_fp(xc_v_pt),
     &                                   memory_fp(xc_h_pt),
     &                                   memory_fp(amo_val_pt),
     &                                   memory_fp(bfng_val_pt),
     &                                   avec,
     &                                   memory_fp(iwrk1),
     &                                   memory_fp(iwrk2),
     &                                   memory_fp(iwrk3),
     &                                   fxa_mo)
                                 endif
                                 call free_memory(iwrk3,'d')
                                 call free_memory(iwrk2,'d')
                                 call free_memory(iwrk1,'d')
                              endif
                           endif
                        else
                           drho_pt = allocate_memory(
     &                               mxpmx*6*n_active_bfn_atm,'d')
                           if (gradcorr_sw) then
                              dgrho_pt = allocate_memory(
     &                                   mxpmx*18*n_active_bfn_atm,'d')
                           endif
                           if (screen_sw) then
                              call den_pert_exp_ao_scr(rks_sw,
     &                             gradcorr_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             memory_int(active_bfn_indx_pt),
     &                             n_active_bfn,
     &                             n_active_bfn_atm,nao,npts,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),adens,bdens,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),mxp,dentol)
                              call dksm_exp_dft_ao_scr(rks_sw,
     &                             gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             memory_int(active_bfn_indx_pt),
     &                             memory_int(active_bfn_atms_pt),
     &                             n_active_bfn,n_active_bfn_atm,
     &                             memory_int(active_wgh_atms_pt),
     &                             n_active_wgh_atm,
     &                             nao,npts,npert,latm,
     &                             memory_fp(wt_pt),memory_fp(gwt_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),fxa_ao,fxb_ao,
     &                             mxp)
                           else
                              call den_pert_exp_ao(rks_sw,gradcorr_sw,
     &                             ao_tag,nao,npts,npert,ngridcentres,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),adens,bdens,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),mxp)
                              call dksm_exp_dft_ao(rks_sw,gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,nao,npts,
     &                             npert,ngridcentres,latm,
     &                             memory_fp(wt_pt),memory_fp(gwt_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),fxa_ao,fxb_ao,
     &                             mxp)
                           endif
                           if (gradcorr_sw) then
                              call free_memory(dgrho_pt,'d')
                           endif
                           call free_memory(drho_pt,'d')
                        endif
                     endif
c
                     if (hess_sw) then
                        if (hess_mo_sw) then
                           if (screen_sw) then
                              iwrk1 = allocate_memory(mxpmx*npert,'d')
                              iwrk2 = allocate_memory(3*mxpmx*npert,'d')
                              iwrk3 = allocate_memory(mxpmx*npert,'d')
                              iwrk4 = allocate_memory(9*mxpmx,'d')
                              iwrk5 = allocate_memory(3*9*mxpmx,'d')
                              iwrk6 = allocate_memory(3*mxpmx,'d')
                              iwrk7 = allocate_memory(3*mxpmx,'d')
                              iwrk8 = allocate_memory(6*mxpmx,'d')
                              iwrk9 = allocate_memory(6*mxpmx,'d')
                              iwrk10= allocate_memory(6*mxpmx,'d')
                              iwrk11= allocate_memory(10*mxpmx,'d')
                              iwrk12= allocate_memory(mxpmx,'d')
                              iiwrk1= allocate_memory(
     &                                n_active_bfn_atm+1,'i')
                              call rks_hess_dft_mo_scr(gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             memory_int(active_bfn_indx_pt),
     &                             memory_int(active_bfn_atms_pt),
     &                             n_active_bfn,n_active_bfn_atm,
     &                             memory_int(active_wgh_atms_pt),
     &                             n_active_wgh_atm,
     &                             npts,naocc,nao,ngridcentres,npert,
     &                             mxp,ao_tag,latm, 
     &                             memory_fp(wt_pt),
     &                             memory_fp(gwt_pt),
     &                             memory_fp(g2wt_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             memory_fp(bfn_3rd_pt),
     &                             memory_fp(xc_e_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             avec,
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             grho,
     &                             memory_fp(iwrk1),
     &                             memory_fp(iwrk2),
     &                             memory_fp(iwrk3),
     &                             memory_fp(iwrk4),
     &                             memory_fp(iwrk5),
     &                             memory_fp(iwrk6),
     &                             memory_fp(iwrk7),
     &                             memory_fp(iwrk8),
     &                             memory_fp(iwrk9),
     &                             memory_fp(iwrk10),
     &                             memory_fp(iwrk11),
     &                             memory_fp(iwrk12),
     &                             memory_int(iiwrk1),
     &                             hess)
                              call free_memory(iiwrk1,'i')
                              call free_memory(iwrk12,'d')
                              call free_memory(iwrk11,'d')
                              call free_memory(iwrk10,'d')
                              call free_memory(iwrk9,'d')
                              call free_memory(iwrk8,'d')
                              call free_memory(iwrk7,'d')
                              call free_memory(iwrk6,'d')
                              call free_memory(iwrk5,'d')
                              call free_memory(iwrk4,'d')
                              call free_memory(iwrk3,'d')
                              call free_memory(iwrk2,'d')
                              call free_memory(iwrk1,'d')
                           else
                              iwrk1 = allocate_memory(mxpmx*npert,'d')
                              iwrk2 = allocate_memory(3*mxpmx*npert,'d')
                              iwrk3 = allocate_memory(mxpmx*npert,'d')
                              iwrk4 = allocate_memory(9*mxpmx,'d')
                              iwrk5 = allocate_memory(3*9*mxpmx,'d')
                              iwrk6 = allocate_memory(3*mxpmx,'d')
                              iwrk7 = allocate_memory(3*mxpmx,'d')
                              iwrk8 = allocate_memory(6*mxpmx,'d')
                              iwrk9 = allocate_memory(6*mxpmx,'d')
                              iwrk10= allocate_memory(6*mxpmx,'d')
                              iwrk11= allocate_memory(10*mxpmx,'d')
                              iwrk12= allocate_memory(mxpmx,'d')
                              call rks_hess_dft_mo(gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,
     &                             npts,naocc,nao,ngridcentres,npert,
     &                             mxp,ao_tag,latm, 
     &                             memory_fp(wt_pt),
     &                             memory_fp(gwt_pt),
     &                             memory_fp(g2wt_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             memory_fp(bfn_3rd_pt),
     &                             memory_fp(xc_e_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             avec,
     &                             memory_fp(amo_val_pt),
     &                             memory_fp(amo_grad_pt),
     &                             grho,
     &                             memory_fp(iwrk1),
     &                             memory_fp(iwrk2),
     &                             memory_fp(iwrk3),
     &                             memory_fp(iwrk4),
     &                             memory_fp(iwrk5),
     &                             memory_fp(iwrk6),
     &                             memory_fp(iwrk7),
     &                             memory_fp(iwrk8),
     &                             memory_fp(iwrk9),
     &                             memory_fp(iwrk10),
     &                             memory_fp(iwrk11),
     &                             memory_fp(iwrk12),
     &                             hess)
                              call free_memory(iwrk12,'d')
                              call free_memory(iwrk11,'d')
                              call free_memory(iwrk10,'d')
                              call free_memory(iwrk9,'d')
                              call free_memory(iwrk8,'d')
                              call free_memory(iwrk7,'d')
                              call free_memory(iwrk6,'d')
                              call free_memory(iwrk5,'d')
                              call free_memory(iwrk4,'d')
                              call free_memory(iwrk3,'d')
                              call free_memory(iwrk2,'d')
                              call free_memory(iwrk1,'d')
                           endif
                        else
                           nprt = npert
                           if (screen_sw) nprt = 3*n_active_bfn_atm
                           drho_pt  = allocate_memory2(mxpmx*2*npert,
     &                                'd',fnm,snm,'drho')
                           ddrhoa_pt = allocate_memory2(mxpmx*9,'d',
     &                                 fnm,snm,'ddrhoa')
                           if (.not.rks_sw) then
                              ddrhob_pt = allocate_memory2(mxpmx*9,'d',
     &                                    fnm,snm,'ddrhob')
                           endif
                           if (gradcorr_sw) then
                              dgrho_pt  = allocate_memory2(
     &                                    mxpmx*6*npert,'d',fnm,
     &                                    snm,'dgrho')
                              ddgrhoa_pt = allocate_memory2(mxpmx*3*9,
     &                                     'd',fnm,snm,'ddgrhoa')
                              if (.not.rks_sw) then
                                 ddgrhob_pt = allocate_memory2(
     &                                        mxpmx*3*9,
     &                                        'd',fnm,snm,'ddgrhob')
                              endif
                           endif
                           if (screen_sw) then
                              call den_pert_exp_ao_scr(rks_sw,
     &                             gradcorr_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             memory_int(active_bfn_indx_pt),
     &                             n_active_bfn,
     &                             n_active_bfn_atm,
     &                             nao,npts,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             adens,bdens,memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),mxp,
     &                             dentol)
                              iiwrk1 = allocate_memory2(
     &                                 n_active_bfn_atm+1,'i',fnm,
     &                                 snm,'iiwrk1')
                              call hess_dft_ao_scr(rks_sw,gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,
     &                             memory_int(active_bfn_list_pt),
     &                             memory_int(active_bfn_indx_pt),
     &                             memory_int(active_bfn_atms_pt),
     &                             n_active_bfn,n_active_bfn_atm,
     &                             memory_int(active_wgh_atms_pt),
     &                             n_active_wgh_atm,
     &                             nao,npert,npts,
     &                             memory_fp(wt_pt),memory_fp(gwt_pt),
     &                             memory_fp(g2wt_pt),grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(ddrhoa_pt),
     &                             memory_fp(ddrhob_pt),
     &                             memory_fp(ddgrhoa_pt),
     &                             memory_fp(ddgrhob_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             memory_fp(bfn_3rd_pt),
     &                             adens,bdens,
     &                             memory_fp(xc_e_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             latm,hess,mxp,dentol,
     &                             memory_int(iiwrk1))
                              call free_memory2(iiwrk1,'i',fnm,
     &                                          snm,'iiwrk1')
                           else
                              call den_pert_exp_ao(rks_sw,gradcorr_sw,
     &                             ao_tag,nao,npts,npert,ngridcentres,
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             adens,bdens,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),mxp)
                              call hess_dft_ao(rks_sw,gradcorr_sw,
     &                             gradwght_sw,gwt_avail_sw,ao_tag,
     &                             nao,npert,ngridcentres,npts,
     &                             memory_fp(wt_pt),
     &                             memory_fp(gwt_pt),
     &                             memory_fp(g2wt_pt),
     &                             grho,
     &                             memory_fp(drho_pt),
     &                             memory_fp(dgrho_pt),
     &                             memory_fp(ddrhoa_pt),
     &                             memory_fp(ddrhob_pt),
     &                             memory_fp(ddgrhoa_pt),
     &                             memory_fp(ddgrhob_pt),
     &                             memory_fp(bfn_val_pt),
     &                             memory_fp(bfng_val_pt),
     &                             memory_fp(bfn_hess_pt),
     &                             memory_fp(bfn_3rd_pt),
     &                             adens,bdens,
     &                             memory_fp(xc_e_pt),
     &                             memory_fp(xc_v_pt),
     &                             memory_fp(xc_dv_pt),
     &                             memory_fp(xc_h_pt),
     &                             memory_fp(xc_dh_pt),
     &                             latm,hess,mxp)
                           endif
                           if (gradcorr_sw) then
                              if (.not.rks_sw) then
                                 call free_memory2(ddgrhob_pt,'d',fnm,
     &                                             snm,'ddgrhob')
                              endif
                              call free_memory2(ddgrhoa_pt,'d',fnm,
     &                                          snm,'ddgrhoa')
                              call free_memory2(dgrho_pt,'d',fnm,
     &                                          snm,'dgrho')
                           endif
                           if (.not.rks_sw) then
                              call free_memory2(ddrhob_pt,'d',fnm,
     &                                          snm,'ddrhob')
                           endif
                           call free_memory2(ddrhoa_pt,'d',fnm,
     &                                       snm,'ddrhoa')
                           call free_memory2(drho_pt,'d',fnm,
     &                                       snm,'drho')
                        endif
                     endif
C     *
C     *Perform various sums for output
C     *

                     if( odb )then
                        if(opg_root())then
                           write(iout,*)'atom',latm,'density',
     +                       ddot(npts,rho,1,memory_fp(wt_pt),1), 
     +                       'xc',xc_e
                        endif
                     endif
c
                     rhoa = ddot(npts,rho(1),1,memory_fp(wt_pt),1)
                     rhob = 0.0d0
                     atom_den(latm,1)=atom_den(latm,1) + rhoa
                     if (.not.rks_sw) then
                        rhob = ddot(npts,rho(1+mxp),1,
     &                              memory_fp(wt_pt),1)
                        atom_den(latm,2)=atom_den(latm,2) + rhob
                     endif
                     memory_fp(irshll+inmtyp(gridt-1)+lrad-1) = 
     +                  memory_fp(irshll+inmtyp(gridt-1)+lrad-1) +
     +                  rhoa + rhob
_IF(qsh)
                     exc_qsh = exc_qsh + xc_e
                     den_qsh = den_qsh + ddot(npts,rho(1,1),1,
     &                                        memory_fp(wt_pt),1)
_ENDIF
                  else
                     del_rho = del_rho + 1.0d0
                  endif

                  if (idfun.ge.2) then
                     if (gradcorr_sw) then
                        call free_memory2(xc_dh_pt,'d',fnm,snm,
     &                                    'xc_dh')
                     endif
                     call free_memory2(xc_h_pt,'d',fnm,snm,
     &                                 'xc_h')
                  endif
                  if (idfun.ge.1) then
                     if (kinetic_sw) then
                        call free_memory2(xc_dt_pt,'d',fnm,snm,
     &                                    'xc_dt')
                     endif
                     if (gradcorr_sw) then
                        call free_memory2(xc_dv_pt,'d',fnm,snm,
     &                                    'xc_dv')
                     endif
                     call free_memory2(xc_v_pt,'d',fnm,snm,
     &                                 'xc_v')
                  endif
                  call free_memory2(xc_e_pt,'d',fnm,snm,'xc_e')

                  if (eval_mo_sw) then
                     if (.not.rks_sw) then
                        if (idmo.ge.1) then
                           call free_memory2(bmo_grad_pt,'d',fnm,snm,
     &                                       'bmo_grad')
                        endif
                        call free_memory2(bmo_val_pt,'d',fnm,snm,
     &                                    'bmo_val')
                     endif
                     if (idmo.ge.1) then
                        call free_memory2(amo_grad_pt,'d',fnm,snm,
     &                                    'amo_grad')
                     endif
                     call free_memory2(amo_val_pt,'d',fnm,snm,
     &                                 'amo_val')
                  endif
                  if (idbfn.ge.3) then
                     call free_memory2(bfn_3rd_pt,'d',fnm,snm,'bfn_3rd')
                  endif
                  if (idbfn.ge.2) then
                     call free_memory2(bfn_hess_pt,'d',fnm,snm,
     &                                 'bfn_hess')
                  endif
                  if (idbfn.ge.1) then
                     call free_memory2(bfng_val_pt,'d',fnm,snm,
     &                                 'bfng_val')
                  endif
                  call free_memory2(bfn_val_pt,'d',fnm,snm,'bfn_val')
c
                  atom_pts(latm)  =atom_pts(latm)+npts
                  endif  ! end if on npts
c
                  nbatch = nbatch + 1

c==== time===
                  if(opg_root() .and. mod(nbatch,freq).eq.0)then
                     call gms_cputime(cpu1)
                     t_now = cpu1(3)
                     write(iout,*)
     &                  'detailed DFT Energy/KS matrix timings after ',
     &                  'batch:',nbatch
                     write(iout,101)'geom',t_geom
                     write(iout,101)'wght',t_becke
                     write(iout,101)'bas ',t_bas
                     write(iout,101)'dens',t_den
                     write(iout,101)'xc  ',t_xc
                     if(ks_sw)write(iout,101)'KS mat',t_ks
                     if(grad_sw)write(iout,101)'Grad',t_grad
                     write(iout,101)'other',
     &                  t_now-t_start-t_geom-t_becke-t_bas-t_den-t_xc-
     &                  t_ks-t_grad
                     
                     call gms_cputime(cpu1)
                     if(cpu_prev .gt. 0.0d0)then
                        write(iout,101)'batch CPU cost',
     &                       cpu1(1)-cpu_prev
                     endif
                     cpu_prev = cpu1(1)

                     if( odb )then
                        write(iout,*)'function after',nbatch
                        write(iout,1101)xc_e
c                       do i=1,npts
c                          write(iout,*)i, xc_vpt(i,1), xc_dvpt(i,1,1), 
c    &                          xc_dvpt(i,1,2), xc_dvpt(i,1,3)
c                       enddo
                     endif
                  endif
c
               enddo 
c
               if( odb )then
                  write(iout,*)'atom den/xc'
                  write(iout,1101)(atom_den(j,1),j=1,4)
                  write(iout,1101)(atom_xce(j),j=1,4)
               endif

_IF(qsh)
               write(59,1102)latm,lrad,nbatch_qsh,npts_qsh,
     &              exc_qsh,den_qsh,exc_qsh,den_qsh
 1102          format(1x,'Qsh ',i4,i6,i4,i4,f20.10,f20.10,1x,g22.16,
     &                1x,g22.16)
_ENDIF

c     
c              ps  end parallel section
c     
               next = ipg_dlbtask() 
            endif

            lrad=lrad+1
            goto 10
         endif
c        end of while loop over lrad

         if( odb )then
            write(iout,*) 'DENSITY:',atom_den(latm,1),atom_den(latm,2),
     &           ipg_nodeid()
         endif

         extout_sw = .false.
      enddo
c
c     Clearup memory for the weights
c
      if (idwght.ge.1) then
         call free_memory2(gwt_pt,'d',fnm,snm,'gwt')
      endif
      call free_memory2(wt2_pt,'d',fnm,snm,'wt2')
      call free_memory2(wt_pt,'d',fnm,snm,'wt')
c
c     Clearup memory for screening tables
c
      if (lhs_sw.or.rhs_sw.or.dksm_sw) then
         call free_memory2(active_chf_pert_pt,'i',fnm,snm,'chf_pert')
      endif
      call free_memory2(active_bfn_atms_pt,'i',fnm,snm,'bfn_atms')
      call free_memory2(active_wgh_atms_pt,'i',fnm,snm,'wgh_atms')
      call free_memory2(active_bfn_indx_pt,'i',fnm,snm,'bfn_indx')
      call free_memory2(active_bfn_list_pt,'i',fnm,snm,'bfn_list')

      call end_time_period(TP_TMP_EXQUAD_INTEG)
      call start_time_period(TP_TMP_EXQUAD_DGOP)

      call pg_dlbpush
c     
c     global sum of accumulated quantities
c     
      call pg_dgop(1001,atom_den(1,1),ngridcentres,'+')
      call pg_dgop(1002,atom_den(1,2),ngridcentres,'+')
      call pg_dgop(1003,atom_xce,ngridcentres,'+')
      call pg_igop(1004,atom_pts,ngridcentres,'+')
      if(ks_sw)then
         call pg_dgop(1005,kma,nao*(nao+1)/2,'+')
         if(.not. rks_sw)call pg_dgop(1006,kmb,nao*(nao+1)/2,'+')
      endif
      if(grad_sw)then
         call pg_dgop(1007,grad,natoms*3,'+')
      endif
      if(dksm_exp_sw)then
         if(ao_in_sw)then
            call pg_dgop(1008,fxa_ao,nao*(nao+1)/2*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1009,fxb_ao,nao*(nao+1)/2*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            call pg_dgop(1010,fxa_mo,nvec*(nvec+1)/2*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1011,fxb_mo,nvec*(nvec+1)/2*npert,'+')
            endif
         endif
      endif
      if(rhs_sw)then
         if(ao_in_sw)then
            call pg_dgop(1012,ba_ao,nao*nao*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1013,bb_ao,nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            call pg_dgop(1014,ba_mo,naocc*(nvec-naocc)*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1015,bb_mo,nbocc*(nvec-nbocc)*npert,'+')
            endif
         endif
      endif
      if(lhs_sw)then
         if(ao_in_sw)then
            call pg_dgop(1016,ga_ao,nao*nao*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1017,gb_ao,nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            call pg_dgop(1018,ga_mo,naocc*(nvec-naocc)*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1019,gb_mo,nbocc*(nvec-nbocc)*npert,'+')
            endif
         endif
      endif
      if(dksm_sw)then
         if(ao_in_sw)then
            call pg_dgop(1020,fa_ao,nao*nao*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1021,fb_ao,nao*nao*npert,'+')
            endif
         endif
         if(mo_in_sw)then
            call pg_dgop(1022,fa_mo,nvec*(nvec+1)/2*npert,'+')
            if(.not.rks_sw)then
               call pg_dgop(1023,fb_mo,nvec*(nvec+1)/2*npert,'+')
            endif
         endif
      endif
      if(hess_sw)then
         call pg_dgop(1024,hess,(natoms*3)**2,'+')
      endif

      call end_time_period(TP_TMP_EXQUAD_DGOP)
C     
C     api entry point
c     
      do latm=1,ngridcentres

         alpha_den = alpha_den + atom_den(latm,1)
         beta_den  = beta_den  + atom_den(latm,2)
         if(rks_sw)then
            totden = totden    + atom_den(latm,1)
         else
            totden = totden    + atom_den(latm,1) + 
     &                           atom_den(latm,2)
         endif
         totpts    = totpts    + atom_pts(latm)
         XC_energy = XC_energy + atom_xce(latm)
      enddo

      call end_time_period(TP_TMP_EXQUAD)

      call pg_dgop(1006,memory_fp(ieshll),inmtyp(ngtypes),'+')
      call pg_dgop(1006,memory_fp(irshll),inmtyp(ngtypes),'+')
      do ig = 1, ngtypes
         fac = 0.0d0
         do latm = 1, ngridcentres
            if (gtype_num(latm).eq.ig) fac = fac + 1.0d0
         enddo
         do irad = inmtyp(ig-1), inmtyp(ig)-1
            memory_fp(ieshll+irad) = memory_fp(ieshll+irad)/fac
            memory_fp(irshll+irad) = memory_fp(irshll+irad)/fac
         enddo
      enddo
cDEBUG
c     call CD_print_dftresults(.true.,.false.,iout)
c     screen_sw = screen_save_sw
cDEBUG
      if(opg_root()) then
         if (nelectrons.ne.0) then
            error  = (totDen-nelectrons)/nelectrons
         else
            error  = totDen
         endif
         errtol = warntol(1,igrid)
         igmin  = 1
         do ig = 2, ngtypes
            if (errtol.lt.warntol(ig,igrid)) then
               errtol = warntol(ig,igrid)
               igmin  = ig
            endif
         enddo
         if (( abs(error).gt.accuracy_warn_tol*errtol.and.
     &         .not.accuracy_warned_sw )
     &       .or.
     &       ( abs(error).gt.accuracy_fatal_tol*errtol.and.
     &         .not.ignore_accuracy_sw )) then
            accuracy_warned_sw = .true.
            write(iout,*)
            if (abs(error).gt.accuracy_fatal_tol*errtol) then
               write(iout,'(a101)')
     &              '*** ERROR: the achieved accuracy does not '//
     &                         'match the target accuracy '//
     &                         'for the current grid setting ***'
            else
               write(iout,'(a103)')
     &              '*** WARNING: the achieved accuracy does not '//
     &                           'match the target accuracy '//
     &                           'for the current grid setting ***'
            endif
            write(iout,*)
            write(iout,*)'*** Integrated Density    = ',totDen
            write(iout,*)'*** Number of electrons   = ',
     &                   dble(nelectrons)
            write(iout,*)'*** Relative error        = ',error
            if (gaccu_num(igmin,igrid).eq.GACC_LOW) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(LOW grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_LOWMEDIUM) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(LOW-MEDIUM grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_MEDIUM) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(MEDIUM grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_MEDIUMHIGH) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(MEDIUM-HIGH grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_HIGH) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(HIGH grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_VERYHIGH) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(VERYHIGH grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_REF) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(REF grid)'
            else if (gaccu_num(igmin,igrid).eq.GACC_SG1) then
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(SG1 grid)'
            else 
               write(iout,*)'*** Target relative error = ',errtol,
     &         '(unknown grid)'
            endif
            if (abs(XC_energy).gt.0.0d0) then
            write(iout,*)'*** XC energy             = ',XC_energy
            endif
            if (abs(error).gt.accuracy_fatal_tol*errtol) then
               if (ignore_accuracy_sw) then
                  write(iout,*)
                  write(iout,*)'*** Error ignored by user request ***'
               else
                  call caserr("Fatal error in accuracy of quadrature")
               endif
            endif
         endif
      endif

      if(opg_root() .and. print_sw(DEBUG_QUAD) )then

c==== time===

         call gms_cputime(cpu1)
         t_now = cpu1(3)
         write(iout,*)
         write(iout,*)'detailed DFT quadrature elapsed timings:'
         write(iout,101)'geom ' ,t_geom
         write(iout,101)'weight',t_becke
         write(iout,101)'pack ' ,t_pack
         write(iout,101)'bas  ' ,t_bas
         write(iout,101)'dens ' ,t_den
         write(iout,101)'xc  '  ,t_xc
         if(ks_sw)write(iout,101)'KS mat',t_ks
         if(grad_sw)write(iout,101)'Grad',t_grad
         write(iout,101)'other ',
     &t_now-t_start-t_geom-t_becke-t_pack-t_bas-t_den-t_xc-t_ks-t_grad
 101     format(a11,f12.5)

c==== time===
c
c...     main results
c
         write(iout,*)
         write(iout,*)'main DFT results:'
         if (rks_sw) then
            write(iout,*)'   integrated density ..... ',totDen
         else
            write(iout,*)'   alpha density .......... ',alpha_Den
            write(iout,*)'   beta density ........... ',beta_Den
            write(iout,*)'   total density .......... ',totDen
         endif
         write(iout,*)'   relative error ......... ',
     &                (totDen-nelectrons)/nelectrons
         if (e_sw) then
            write(iout,*)'   xc energy .............. ',XC_energy
         endif
         if (abs(J_energy).gt.0.0d0) then
            write(iout,*)'   fitted coulomb energy .. ',J_energy
         endif
         if (abs(lmult).gt.0.0d0) then
            write(iout,*)'   multipole coulomb energy ',lmult
         endif
c     
c...     screening stats
c     
         write(iout,*)
         write(iout,*)'screening statistics for node 0:'
         write(iout,*)'    no. batches of quadrature points     ',
     +        '            = ',nbatch
         write(iout,*)'    total number of quadrature points    ',
     +        '            = ',ntot
c        del_exp = del_exp / dble(nbatch)
         del_psi = del_psi / dble(nbatch)

         write(iout,*)'    average no. bfs per batch ignored on ',
     +        'basis value = ',nint(del_psi)
         write(iout,*)'    no. batches ignored based on ',
     +        'density value       = ',nint(del_rho)
         write(iout,*)'    no. batches ignored based on ',
     +        'weight value        = ',nint(del_wgh)
         write(iout,103)npack,nint(float(npack)/float(ntot)*100.0d0)
 103     format(5x,'number of points with zero wts  ',
     +          '                 = ',i8,' (',i3,
     +          ' percent)')

         write(iout,*)
         write(iout,*)'XC energy by atom'
         if (rks_sw) then
            write(iout,*)'                         Total'
            write(iout,*)'Atom       Energy        Density'
         else
            write(iout,*)'                         Total         ',
     +                   'Alpha         Beta'
            write(iout,*)'Atom       Energy        Density       ',
     +                   'Density       Density'
         endif
         do j=1,ngridcentres
            if (rks_sw) then
               write(iout,104)j,atom_xce(j),atom_den(j,1)
            else
               write(iout,105)j,atom_xce(j),atom_den(j,1)+atom_den(j,2),
     +                                      atom_den(j,1),atom_den(j,2)
            endif
         enddo
 104     format(1x,i4,2f14.7)
 105     format(1x,i4,4f14.7)

         write(iout,*)
         write(iout,*)'XC energy by radial shell'
         do j=1,ngtypes
            write(iout,*)'Grid type',j
            write(iout,*)'Shell   Energy      Density     Radius'
            do i=1,nradpt_num(j)
c              write(iout,102)i,eshell(i),prpt(itmx,i)
               write(iout,102)i,memory_fp(ieshll+inmtyp(j-1)+i-1),
     +                          memory_fp(irshll+inmtyp(j-1)+i-1),
     +                          prpt(j,i)
 102           format(1x,i3,3f12.7)
            enddo
         enddo
         write(iout,'(a82)')
     &                '______________________________'//
     &                '______________________________'//
     &                '______________________'
         write(iout,*)
      endif

      if (e_sw) then
         call auto_prune(igrid,nradpt_num,inmtyp,memory_fp(ieshll),prpt)
      endif

      if (kinetic_sw) then
         call free_memory2(tau_pt,'d',fnm,snm,'tau')
      endif

      call free_memory2(irshll,'d',fnm,snm,'irshll')
      call free_memory2(ieshll,'d',fnm,snm,'ieshll')
      call free_memory2(bfn_radii_pt,'d',fnm,snm,'bfn_radii')
_IF(qsh)
      close(unit=59)
_ENDIF

c     call exitc(0)
      imemusage = pop_memory_count(imemcount)
      imemcount = push_memory_estimate()
      call memreq_exquad(memory_fp,memory_int,
     &     igrid,
     &     ao_tag,nao,nvec,naocc,nbocc,npert,
     &     geometric_pert_sw,
     &     prpt,prwt,
     &     e_sw,ks_sw,grad_sw,
     &     dksm_exp_sw,rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &     mxp0,ao_in_sw,mo_in_sw,accuracy,iout)
      imemestimate = pop_memory_estimate(imemcount)
      if (opg_root().and.imemusage.ne.imemestimate) then
         if (
_IF(debug)
     &       imemusage.ne.imemestimate
_ELSE
     &       imemusage.lt.0.9d0*imemestimate.or.
     &       imemusage.gt.imemestimate
_ENDIF
     &   ) then
            write(iout,*)'*** estimated memory usage = ',imemestimate,
     &                   ' words'
            write(iout,*)'*** actual    memory usage = ',imemusage   ,
     &                   ' words'
            write(iout,*)'*** WARNING: the memory usage estimates for ',
     &                   'exquad seem to be incorrect'
         endif
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine scale_accumulators(nao,nvec,naocc,nbocc,natoms,npert,
     &           rks_sw,ao_in_sw,mo_in_sw,ks_sw,grad_sw,dksm_exp_sw,
     &           lhs_sw,rhs_sw,dksm_sw,hess_sw,kma,kmb,grad,
     &           fxa_ao,fxb_ao,fxa_mo,fxb_mo,ba_ao,bb_ao,ba_mo,bb_mo,
     &           ga_ao,gb_ao,ga_mo,gb_mo,fa_ao,fb_ao,fa_mo,fb_mo,
     &           hess)
      implicit none
c     
c     Scale accumulators as we will be summing at the end of the
c     parallel section.
c     
c     Inputs:
c
      integer nao,nvec,naocc,nbocc,natoms,npert
      logical rks_sw,ao_in_sw,mo_in_sw
      logical ks_sw,grad_sw,dksm_exp_sw,lhs_sw,rhs_sw,dksm_sw,hess_sw
c
c     Inputs/Outputs:
c
      REAL kma(nao*(nao+1)/2)
      REAL kmb(nao*(nao+1)/2)
      REAL grad(3,natoms)
      REAL ba_ao(nao,nao,npert)          ! the alpha b vector (Au=b)
      REAL bb_ao(nao,nao,npert)          ! the beta  b vector (Au=b)
      REAL ba_mo(naocc,nvec-naocc,npert) ! the alpha b vector (Au=b)
      REAL bb_mo(nbocc,nvec-nbocc,npert) ! the beta  b vector (Au=b)
      REAL ga_ao(nao,nao,npert)          ! the alpha Au matrix (Au=b)
      REAL gb_ao(nao,nao,npert)          ! the beta  Au matrix (Au=b)
      REAL ga_mo(naocc,nvec-naocc,npert) ! the alpha Au matrix (Au=b)
      REAL gb_mo(nbocc,nvec-nbocc,npert) ! the beta  Au matrix (Au=b)
      REAL fa_ao(nao,nao,npert)          ! the alpha perturbed fock
                                         ! matrix
      REAL fb_ao(nao,nao,npert)          ! the beta  perturbed fock
                                         ! matrix
      REAL fa_mo(nvec*(nvec+1)/2,npert)
      REAL fb_mo(nvec*(nvec+1)/2,npert)
      REAL fxa_ao(nao*(nao+1)/2,npert)
      REAL fxb_ao(nao*(nao+1)/2,npert)
      REAL fxa_mo(nvec*(nvec+1)/2,npert)
      REAL fxb_mo(nvec*(nvec+1)/2,npert)
      REAL hess(npert,npert)
c
c     Local data
c
      REAL fac
      integer nproc
      integer i, j, k, latm
c
c     Functions
c
      integer ipg_nnodes
c
      nproc = ipg_nnodes()
      if (nproc.eq.1) return
c
      fac = 1.0d0/dble(nproc)
c
      if (ks_sw) then
         do i=1,(nao+1)*nao/2
            kma(i) = kma(i) * fac
         enddo
         if (.not.rks_sw) then
            do i=1,(nao+1)*nao/2
               kmb(i) = kmb(i) * fac
            enddo
         endif
      endif
c
      if (grad_sw) then
         do latm = 1, natoms
            do i = 1, 3
               grad(i,latm) = grad(i,latm) * fac
            enddo
         enddo
      endif
c
      if(dksm_exp_sw)then
         if (ao_in_sw) then
            do i=1,npert
               do j=1,nao*(nao+1)/2
                  fxa_ao(j,i) = fxa_ao(j,i) * fac
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nao*(nao+1)/2
                     fxb_ao(j,i) = fxb_ao(j,i) * fac
                  enddo
               enddo
            endif
         endif
         if (mo_in_sw) then
            do i=1,npert
               do j=1,nvec*(nvec+1)/2
                  fxa_mo(j,i) = fxa_mo(j,i) * fac
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nvec*(nvec+1)/2
                     fxb_mo(j,i) = fxb_mo(j,i) * fac
                  enddo
               enddo
            endif
         endif
      endif
c
      if(rhs_sw)then
         if (ao_in_sw) then
            do i=1,npert
               do j=1,nao
                  do k=1,nao
                     ba_ao(k,j,i) = ba_ao(k,j,i) * fac
                  enddo
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nao
                     do k=1,nao
                        bb_ao(k,j,i) = bb_ao(k,j,i) * fac
                     enddo
                  enddo
               enddo
            endif
         endif
         if (mo_in_sw) then
            do i=1,npert
               do j=1,nvec-naocc
                  do k=1,naocc
                     ba_mo(k,j,i) = ba_mo(k,j,i) * fac
                  enddo
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nvec-nbocc
                     do k=1,nbocc
                        bb_mo(k,j,i) = bb_mo(k,j,i) * fac
                     enddo
                  enddo
               enddo
            endif
         endif
      endif
c
      if(lhs_sw)then
         if (ao_in_sw) then
            do i=1,npert
               do j=1,nao
                  do k=1,nao
                     ga_ao(k,j,i) = ga_ao(k,j,i) * fac
                  enddo
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nao
                     do k=1,nao
                        gb_ao(k,j,i) = gb_ao(k,j,i) * fac
                     enddo
                  enddo
               enddo
            endif
         endif
         if (mo_in_sw) then
            do i=1,npert
               do j=1,nvec-naocc
                  do k=1,naocc
                     ga_mo(k,j,i) = ga_mo(k,j,i) * fac
                  enddo
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nvec-nbocc
                     do k=1,nbocc
                        gb_mo(k,j,i) = gb_mo(k,j,i) * fac
                     enddo
                  enddo
               enddo
            endif
         endif
      endif
c
      if(dksm_sw)then
         if (ao_in_sw) then
            do i=1,npert
               do j=1,nao
                  do k=1,nao
                     fa_ao(k,j,i) = fa_ao(k,j,i) * fac
                  enddo
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nao
                     do k=1,nao
                        fb_ao(k,j,i) = fb_ao(k,j,i) * fac
                     enddo
                  enddo
               enddo
            endif
         endif
         if (mo_in_sw) then
            do i=1,npert
               do j=1,nvec*(nvec+1)/2
                  fa_mo(j,i) = fa_mo(j,i) * fac
               enddo
            enddo
            if (.not.rks_sw) then
               do i=1,npert
                  do j=1,nvec*(nvec+1)/2
                     fb_mo(j,i) = fb_mo(j,i) * fac
                  enddo
               enddo
            endif
         endif
      endif
c
      if(hess_sw)then
         do j=1,3*natoms
            do i=1,3*natoms
               hess(i,j) = hess(i,j) * fac
            enddo
         enddo
      endif
      end
c
c-----------------------------------------------------------------------
c
       subroutine group_points(p,w,iperm,po,wo,nang)
       implicit none
c
c      The grid points will be split into batches and all subsequent
c      processing executes on a per batch basis. For screening the
c      batches should be as compact as possible (i.e. cover a small
c      part of the surface of the unit sphere). To achieve this the
c      grid points are sorted:
c
c      - first into the quadrants:
c        1. y>0  and z>=0, or y=0 and z=0 and x>=0
c        2. y<=0 and z>0
c        3. y<0  and z<=0
c        4. y>=0 and z<0,  or y=0 and z=0 and x<0
c      - second within each quadrant according to x:
c        1. decreasing
c        2. increasing
c        3. decreasing
c        4. increasing
c
c      This ensures that any arbitrary set of consecutive points is 
c      relatively compact.
c
c      Inputs
c
       integer nang
       REAL p(3,nang),w(nang)
c
c      Outputs
c
       integer iperm(nang)
       REAL po(3,nang),wo(nang)
c
c      Local variables
c
       integer npoints(4)
       integer opoints(5)
       integer ipoints(4)
       REAL t
       integer i, j, it
       integer ipt
       integer ilvl
       integer iq
       integer ihome
       REAL fac
c
c      Code
c
       do i=1,4
          npoints(i)=0
       enddo
       do i=1,5
          opoints(i)=0
       enddo
c
c      First count the number of points in each quadrant
c
       do i=1,nang
          if (p(2,i).gt.0.0d0.and.p(3,i).ge.0.0d0) then
             npoints(1)=npoints(1)+1
          else if (p(2,i).le.0.0d0.and.p(3,i).gt.0.0d0) then
             npoints(2)=npoints(2)+1
          else if (p(2,i).lt.0.0d0.and.p(3,i).le.0.0d0) then
             npoints(3)=npoints(3)+1
          else if (p(2,i).ge.0.0d0.and.p(3,i).lt.0.0d0) then
             npoints(4)=npoints(4)+1
          else if (p(2,i).eq.0.0d0.and.p(3,i).eq.0.0d0) then
             if (p(1,i).ge.0.0d0) then
                npoints(1)=npoints(1)+1
             else
                npoints(4)=npoints(4)+1
             endif
          else
             write(*,*)"group_points confused: where does this point ",
     &                 "go?"
             write(*,*)p(1,i),p(2,i),p(3,i)
             call caserr("subroutine group_points failed in count")
          endif
       enddo
c
c      Set up offset for the quadrants
c
       do i=2,5
          opoints(i)=opoints(i-1)+npoints(i-1)
       enddo
       do i=1,4
          ipoints(i)=opoints(i)
       enddo
c
c      Sort the point into their appropriate quadrants
c
       do i=1,nang
          if (p(2,i).gt.0.0d0.and.p(3,i).ge.0.0d0) then
             ihome=1
          else if (p(2,i).le.0.0d0.and.p(3,i).gt.0.0d0) then
             ihome=2
          else if (p(2,i).lt.0.0d0.and.p(3,i).le.0.0d0) then
             ihome=3
          else if (p(2,i).ge.0.0d0.and.p(3,i).lt.0.0d0) then
             ihome=4
          else if (p(2,i).eq.0.0d0.and.p(3,i).eq.0.0d0) then
             if (p(1,i).ge.0.0d0) then
                ihome=1
             else
                ihome=4
             endif
          else
             write(*,*)"group_points confused: where does this point ",
     &                 "go?"
             write(*,*)p(1,i),p(2,i),p(3,i)
             call caserr("subroutine group_points failed")
          endif
          ipoints(ihome)=ipoints(ihome)+1
          iperm(ipoints(ihome))=i
       enddo
c
c      Sort the grid points in each quadrant according to x
c
       fac = 1.0d0
       do iq=1,4
          do i=opoints(iq)+1,opoints(iq+1)
             ipt=i
             do j=i+1,opoints(iq+1)
                if (fac*p(1,iperm(j)).gt.fac*p(1,iperm(ipt))) then
                   ipt=j
                endif
             enddo
             it=iperm(i)
             iperm(i)=iperm(ipt)
             iperm(ipt)=it
          enddo
          fac=-fac
       enddo
c
       do ipt=1,nang
          do i=1,3
             po(i,ipt)=p(i,iperm(ipt))
          enddo
          wo(ipt)=w(iperm(ipt))
       enddo
c
       end
c
c-----------------------------------------------------------------------
c
      subroutine route_batch(nvec,naocc,nbocc,nao,n_active_bfn,
     &                       ao_in_sw,mo_in_sw,rks_sw,gradcorr_sw,
     &                       e_sw,ks_sw,grad_sw,dksm_exp_sw,
     &                       rhs_sw,lhs_sw,dksm_sw,hess_sw,
     &                       eval_mo_sw,idmo,na_mo,nb_mo,den_mo_sw,
     &                       dksm_exp_mo_sw,rhs_mo_sw,lhs_mo_sw,
     &                       dksm_mo_sw,hess_mo_sw)
      implicit none
c
c     Decides which quantities should be evaluated from MOs if any.
c
c     This subroutine analyses the cost of calculating the various 
c     quantities from MOs and AOs and selects the cheapest option 
c     available. It sets the various options appropriately to evaluate 
c     at least the required quantities.
c
c     Input:
c
      integer nvec         ! the total number of MOs
      integer naocc        ! the number of occupied alpha MOs
      integer nbocc        ! the number of occupied beta  MOs
      integer nao          ! the total number of basis functions
      integer n_active_bfn ! the number of active basis functions for 
                           ! this batch
      logical ao_in_sw     ! true if adenm and bdenm available
      logical mo_in_sw     ! true if avec  and bvec  available
      logical rks_sw       ! true if this is a closed shell calculation
      logical gradcorr_sw  ! true if GGA functionals are used
      logical e_sw         ! true if energy to be evaluated
      logical ks_sw        ! true if Kohn-Sham matrix to be evaluated
      logical grad_sw      ! true if forces to be evaluated
      logical dksm_exp_sw  ! true if explicit dKS matrix to be evaluated
      logical rhs_sw       ! true if right-hand-side to be evaluated
      logical lhs_sw       ! true if left-hand-side to be evaluated
      logical dksm_sw      ! true if perturbed KS matrix to be evaluated
      logical hess_sw      ! true if hessian matrix to be evaluated
c
c     Input/Output:
c
      logical eval_mo_sw ! true if MOs should be evaluated
      integer idmo       ! the number of derivatives of the MOs required
      integer na_mo      ! the number of alpha MOs required
      integer nb_mo      ! the number of beta  MOs required
c
c     Output:
c
      logical den_mo_sw      ! true if density be calculated from MOs
      logical dksm_exp_mo_sw ! true if explicit derivative Kohn-Sham
                             ! matrix be calculated from MOs
      logical rhs_mo_sw      ! true if RHS be calculated from MOs
      logical lhs_mo_sw      ! true if LHS be calculated from MOs
      logical dksm_mo_sw     ! true if perturbed Kohn-Sham matrix be 
                             ! calculated from MOs
      logical hess_mo_sw     ! true if hessian be calculated from MOs
c
c     Local:
c                            ! Note: all costs are in FLOPs/grid point
c
      integer cost_occ_amo   ! the cost of evaluating the occupied alpha
      integer cost_occ_bmo   ! and beta MOs respectively
      integer cost_occ_admo  ! the cost of evaluating the gradient of
      integer cost_occ_bdmo  ! the occupied alpha and beta MOs 
                             ! respectively
      integer cost_occ_ad2mo ! the cost of evaluating the hessian of
      integer cost_occ_bd2mo ! the occupied alpha and beta MOs 
                             ! respectively
c
      integer cost_den_amo   ! cost of evaluating rho from MOs
      integer cost_den_bmo   ! cost of evaluating rho from MOs
      integer cost_den_admo  ! cost of evaluating grad rho from MOs
      integer cost_den_bdmo  ! cost of evaluating grad rho from MOs
      integer cost_den_aao   ! cost of evaluating rho from AOs
      integer cost_den_bao   ! cost of evaluating rho from AOs
      integer cost_den_adao  ! cost of evaluating grad rho from AOs
      integer cost_den_bdao  ! cost of evaluating grad rho from AOs
c
c     Code:
c
      if (ao_in_sw.and..not.mo_in_sw) then
c
c        No vectors available therefore MOs can not be evaluated.
c
         den_mo_sw      = .false.
         eval_mo_sw     = .false.
         dksm_exp_mo_sw = .false.
         rhs_mo_sw      = .false.
         lhs_mo_sw      = .false.
         dksm_mo_sw     = .false.
         hess_mo_sw     = .false.
c
      else if (.not.ao_in_sw.and.mo_in_sw) then
c
c        Only vectors available therefore MOs have to be evaluated.
c
         den_mo_sw      = .true.
         eval_mo_sw     = .true.
         dksm_exp_mo_sw = dksm_exp_sw
         rhs_mo_sw      = rhs_sw
         lhs_mo_sw      = lhs_sw
         dksm_mo_sw     = dksm_sw
         hess_mo_sw     = hess_sw
c
      else if (ao_in_sw.and.mo_in_sw) then
c
c        Both vectors and density matrices available so we have to 
c        choose
c
         cost_occ_amo   = 2*naocc*n_active_bfn
         cost_occ_bmo   = 2*nbocc*n_active_bfn
         cost_occ_admo  = 3*cost_occ_amo
         cost_occ_bdmo  = 3*cost_occ_bmo
         cost_occ_ad2mo = 6*cost_occ_amo
         cost_occ_bd2mo = 6*cost_occ_bmo
c
         cost_den_amo   = 2*naocc
         cost_den_bmo   = 2*nbocc
         cost_den_admo  = 3*cost_den_amo
         cost_den_bdmo  = 3*cost_den_bmo
c
         cost_den_aao   = 2*n_active_bfn*n_active_bfn+2*n_active_bfn
         cost_den_bao   = 2*n_active_bfn*n_active_bfn+2*n_active_bfn
         cost_den_adao  = 6*n_active_bfn
         cost_den_bdao  = 6*n_active_bfn
c
c        First establish whether we want to calculate the MOs at all
c        I.e. there must be at least 1 component for which the cost
c        of evaluating the MOs plus the component from the MOs is 
c        cheaper than evaluating that component from the AOs
c      
         if (rks_sw) then
            if (gradcorr_sw) then
               eval_mo_sw = 
     &         (cost_occ_amo+cost_occ_admo+cost_den_amo+cost_den_admo
     &         .lt.cost_den_aao+cost_den_adao)
            else
               eval_mo_sw = 
     &         (cost_occ_amo+cost_den_amo
     &         .lt.cost_den_aao)
            endif
         else
            if (gradcorr_sw) then
               eval_mo_sw = 
     &         (cost_occ_amo+cost_occ_admo+cost_den_amo+cost_den_admo
     &         +cost_occ_bmo+cost_occ_bdmo+cost_den_bmo+cost_den_bdmo
     &         .lt.cost_den_aao+cost_den_adao
     &            +cost_den_bao+cost_den_bdao)
            else
               eval_mo_sw = 
     &         (cost_occ_amo+cost_den_amo
     &         +cost_occ_bmo+cost_den_bmo
     &         .lt.cost_den_aao
     &            +cost_den_bao)
            endif
         endif
c
c        If we are going to calculate the MOs anyway then decide per
c        component whether it is best to calculate them from the MOs
c
         if (eval_mo_sw) then
            if (rks_sw) then
               if (gradcorr_sw) then
                  den_mo_sw = (cost_den_amo+cost_den_admo
     &                         .lt.cost_den_aao+cost_den_adao)
               else
                  den_mo_sw = (cost_den_amo.lt.cost_den_aao)
               endif
            else
               if (gradcorr_sw) then
                  den_mo_sw = (cost_den_amo+cost_den_admo
     &                        +cost_den_bmo+cost_den_bdmo
     &                         .lt.cost_den_aao+cost_den_adao
     &                            +cost_den_bao+cost_den_bdao)
               else
                  den_mo_sw = (cost_den_amo+cost_den_bmo
     &                         .lt.cost_den_aao+cost_den_bao)
               endif
            endif
         else
            den_mo_sw = .false.
         endif
c
c        For now just switch to AO basis (still need to sort out the
c        cost functions).
c
         dksm_exp_mo_sw = .false.
         rhs_mo_sw      = .false.
         lhs_mo_sw      = .false.
         dksm_mo_sw     = .false.
         hess_mo_sw     = .false.
c
      else
         call caserr('route_batch: No MOs and no AOs???')
      endif
c
c     What MOs (if any) are we going to need?
c
      na_mo = 0
      nb_mo = 0
      idmo  = 0
      if (eval_mo_sw.and.gradcorr_sw) idmo = 1
      if (den_mo_sw) then
         if (rks_sw) then
            na_mo = naocc
            nb_mo = 0
         else
            na_mo = naocc
            nb_mo = nbocc
         endif
      endif
      if (dksm_exp_sw.or.lhs_mo_sw.or.rhs_mo_sw.or.dksm_mo_sw) then
         if (rks_sw) then
            na_mo = nvec
            nb_mo = 0
         else
            na_mo = nvec
            nb_mo = nvec
         endif
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_lhs_atms(lhs_atms,natm,
     &                        active_bfn_atms,n_active_atm,
     &                        lhs_atm_list,n_lhs_atm)
      implicit none
c
c     As not all left-hand-sides have to present we have to assure that
c     the screened version only loops over the left-hand-sides that
c     are present and for which the contributions of the current batch
c     of grid points will be significant.
c
c     Given that lhs_atms are the atoms for which we have the
c     left-hand-sides, and active_bfn_atms are the perturbation offsets
c     for which significant contributions are expected, we now need
c     to compile lhs_atm_list which includes only the offsets which
c     relate to both lhs_atms and active_bfn_atms.
c
      integer natm, lhs_atms(natm)
      integer n_active_atm, active_bfn_atms(n_active_atm)
      integer n_lhs_atm, lhs_atm_list(*)
c
      integer ia, il
c
      n_lhs_atm = 0
      do 10 ia=1,n_active_atm
         do il=1,natm
            if (active_bfn_atms(ia).eq.3*(lhs_atms(il)-1)) then
               n_lhs_atm=n_lhs_atm+1
c              note: the left-hand-sides are a compressed set
               lhs_atm_list(n_lhs_atm)=3*(il-1)
               goto 10
            endif
         enddo
 10   continue
      end
c
c-----------------------------------------------------------------------
c
      subroutine grid_init_gen1
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_xc)
c
c     Initialises the grid parameters according to the generation 1
c     settings.
c
c     Some parameters have been explicitly set in the input phase,
c     some have been defined more implicitly. This routine propagates
c     the generic parameters to undefined specific parameters. 
c
c     This process is performed in 3 steps. 
c     1) Foreach grid, foreach parameter: if the parameter is undefined
c        copy the value of the generic grid (grid 0).
c     2) Foreach grid g: if a parameter is undefined set it in 
c        accordance with the grid specified by value gaccu_num(g).
c     3) Foreach grid g: if ang_prune_scheme(g) is either AP_SG1 or
c        AP_MHL then set the radial zone stuff to the appropriate 
c        values.
c
      REAL srad, sg1rad
      external srad, sg1rad
c
      integer ig, i, it
      integer rad_low(7),      ang_low(7)
      integer rad_medium(7),   ang_medium(7)
      integer rad_high(7),     ang_high(7)
      integer rad_veryhigh(7), ang_veryhigh(7)
      integer rad_ref(7),      ang_ref(7)
c
      integer row_by_atomnum, atom_num_by_grid, row, atmnum
c
      data rad_low     / 15, 20, 25, 30, 35, 40, 45/
      data rad_medium  / 25, 35, 45, 55, 65, 75, 85/
      data rad_high    / 35, 50, 65, 80, 95,110,125/
      data rad_veryhigh/ 45, 65, 85,105,125,145,165/
      data rad_ref     /300,300,300,300,300,300,300/
c
      data ang_low     / 194, 194, 194, 194, 194, 194, 194/
      data ang_medium  / 302, 302, 302, 302, 302, 302, 302/
      data ang_high    / 434, 434, 434, 434, 434, 434, 434/
      data ang_veryhigh/ 590, 590, 590, 590, 590, 590, 590/
      data ang_ref     /2354,2354,2354,2354,2354,2354,2354/
c
      rad_scale_scheme = SC_GAM1
c
c     loop over all grids
c
      do it=1,max_grids
c
c        Sort out the weighting scheme
c
         if (weight_scheme(it).eq.DFT_UNDEF) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               weight_scheme(it)=WT_MHL4SSFSCR
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_VERYHIGH) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_REF) then
               weight_scheme(it)=WT_MHL
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               weight_scheme(it)=WT_BECKE
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
c
c        Sort out the screening stuff
c
         if (rhotol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               rhotol(it)=1.0d-8
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               rhotol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               rhotol(it)=1.0d-12
            else if (gaccu_num(0,it).eq.GACC_VERYHIGH) then
               rhotol(it)=1.0d-14
            else if (gaccu_num(0,it).eq.GACC_REF) then
               rhotol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               rhotol(it)=1.0d-10
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
         if (dentol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_VERYHIGH) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_REF) then
               dentol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               dentol(it)=0.0d0
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
         if (wghtol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_VERYHIGH) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_REF) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               wghtol(it)=1.0d-20
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
c
c        Phase 1)
c
         do ig=1,ngtypes
            if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
               psitol(ig,it)=psitol(0,it)
            endif
            if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
               warntol(ig,it)=warntol(0,it)
            endif
            if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
               radm_num(ig,it)=radm_num(0,it)
            endif
            if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
               grid_scale(ig,it)=grid_scale(0,it)
            endif
            if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               grid_atom_radius(ig,it)=grid_atom_radius(0,it)
            endif
            if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               weight_atom_radius(ig,it)=weight_atom_radius(0,it)
            endif
            if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               prune_atom_radius(ig,it)=prune_atom_radius(0,it)
            endif
            if (screen_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               screen_atom_radius(ig,it)=screen_atom_radius(0,it)
            endif
            do i = 1, maxradzn-1
               if (bnd_radzn(i,ig,it).eq.dble(DFT_UNDEF)) then
                  bnd_radzn(i,ig,it)=bnd_radzn(i,0,it)
               endif
            enddo
            do i = 1, maxradzn
               if (angpt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  angpt_radzn_num(i,ig,it)=angpt_radzn_num(i,0,it)
               endif
               if (phipt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  phipt_radzn_num(i,ig,it)=phipt_radzn_num(i,0,it)
               endif
               if (thetpt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  thetpt_radzn_num(i,ig,it)=thetpt_radzn_num(i,0,it)
               endif
            enddo
            if (radzones_num(ig,it).eq.DFT_UNDEF) then
               radzones_num(ig,it)=radzones_num(0,it)
            endif
            if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
               ang_prune_scheme(ig,it)=ang_prune_scheme(0,it)
            endif
            if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
               rad_grid_scheme(ig,it)=rad_grid_scheme(0,it)
            endif
            if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
               ang_grid_scheme(ig,it)=ang_grid_scheme(0,it)
            endif
            if (radpt_num(ig,it).eq.DFT_UNDEF) then
               radpt_num(ig,it)=radpt_num(0,it)
            endif
            if (gaccu_num(ig,it).eq.DFT_UNDEF) then
               gaccu_num(ig,it)=gaccu_num(0,it)
            endif
         enddo
c
c        Phase 2)
c
         do ig=1,ngtypes
c
c           Default grid for non-atomic centres is that of carbon
c           See message in master.m
c
            atmnum = atom_num_by_grid(ig)
            if(atmnum.eq.0)atmnum = 6
            row = row_by_atomnum(atmnum)
c
            if (radpt_num(ig,it).eq.DFT_UNDEF) then
               radpt_num(ig,it)=radnpt_row(row)
            endif
            if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF) then
               angpt_radzn_num(1,ig,it)=angnpt_row(row)
            endif
c
            if (gaccu_num(ig,it).eq.GACC_LOW) then
c
c              Low accuracy grid (rel.err. 1.0d-4)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-6
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-4
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_low(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_low(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_MEDIUM) then
c
c              Medium accuracy grid (rel.err. 1.0d-6)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-7
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-6
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_medium(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_medium(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_HIGH) then
c
c              High accuracy grid (rel.err. 1.0d-8)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-13
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-8
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_high(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_high(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_VERYHIGH) then
c
c              Veryhigh accuracy grid (rel.err. 1.0d-10)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-16
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-10
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_veryhigh(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_veryhigh(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_REF) then
c
c              Reference accuracy grid
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-20
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-12
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_ref(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_ref(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_SG1) then
c
c              SG1 (Standard grid #1)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=0.0d0
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-4
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=sg1rad(atmnum,0,.false.)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_SG1
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_SG1
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=302
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=50
               endif
            else
               call caserr
     +              ('CD_defaults is messed up or never called!!!')
            endif
         enddo
c
c        Phase 3)
c
         do ig = 1, ngtypes
            if (ang_prune_scheme(ig,it).eq.AP_MHL) then
               call mhl_prune(ig,it)
            else if (ang_prune_scheme(ig,it).eq.AP_AUTO) then
c
c              initially use MHL pruning, the grid parameter will later
c              be modified in auto_prune
c
               call mhl_prune(ig,it)
            else if (ang_prune_scheme(ig,it).eq.AP_SG1) then
               call sg1_prune(ig,it,1)
            else if (ang_prune_scheme(ig,it).eq.AP_SG1a) then
               call sg1_prune(ig,it,2)
            else if (ang_prune_scheme(ig,it).ne.AP_RADZONE) then
               call caserr('grid_init: error angular pruning scheme')
            endif
         enddo
c
      enddo ! it
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine grid_init_gen2
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_xc)
c
c     Initialises the grid parameters according to the generation 2
c     settings.
c
c     Some parameters have been explicitly set in the input phase,
c     some have been defined more implicitly. This routine propagates
c     the generic parameters to undefined specific parameters. 
c
c     This process is performed in 3 steps. 
c     1) Foreach grid, foreach parameter: if the parameter is undefined
c        copy the value of the generic grid (grid 0).
c     2) Foreach grid g: if a parameter is undefined set it in 
c        accordance with the grid specified by value gaccu_num(g).
c     3) Foreach grid g: if ang_prune_scheme(g) is either AP_SG1 or
c        AP_MHL then set the radial zone stuff to the appropriate 
c        values.
c
      REAL srad, sg1rad, mkrad, lograd
      external srad, sg1rad, mkrad, lograd
c
      integer ig, i, it
      integer rad_low(7),      ang_low(7)
      integer rad_lowmed(7),   ang_lowmed(7)
      integer rad_medium(7),   ang_medium(7)
      integer rad_medhigh(7),  ang_medhigh(7)
      integer rad_high(7),     ang_high(7)
      integer rad_ref(7),      ang_ref(7)
c
      integer row_by_atomnum, atom_num_by_grid, row, atmnum
c
      data rad_low     / 20, 40, 35, 40, 40, 40, 35/
c     data rad_lowmed  / 49, 52, 58, 59, 59, 67, 58/
      data rad_lowmed  / 38, 49, 51, 54, 54, 58, 51/
      data rad_medium  / 60, 60, 70, 70, 70, 80, 70/
c     data rad_medhigh /145,107,131,146,177,147, 95/
      data rad_medhigh /110, 89,108,116,133,122, 87/
      data rad_high    /180,130,160,180,220,180,110/
      data rad_ref     /300,300,300,300,300,300,300/
c
      data ang_low     / 194, 266, 302, 302, 302, 302, 266/
c     data ang_lowmed  / 230, 434, 434, 434, 590, 434, 350/
      data ang_lowmed  / 230, 434, 434, 434, 434, 434, 350/
      data ang_medium  / 266, 590, 590, 590, 770, 590, 434/
c     data ang_medhigh / 590, 974,1202,1202,1202,1202, 974/
      data ang_medhigh / 434, 770, 974, 974,1202, 974, 770/
      data ang_high    / 770,1202,1730,1454,1730,1454,1454/
      data ang_ref     /2354,2354,2354,2354,2354,2354,2354/
c
c     loop over all grids
c
      do it=1,max_grids
c
c        Sort out the weighting scheme
c
         if (weight_scheme(it).eq.DFT_UNDEF) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               weight_scheme(it)=WT_MHL4SSFSCR
            else if (gaccu_num(0,it).eq.GACC_LOWMEDIUM) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_MEDIUMHIGH) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               weight_scheme(it)=WT_MHL8SSFSCR
            else if (gaccu_num(0,it).eq.GACC_REF) then
               weight_scheme(it)=WT_MHL
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               weight_scheme(it)=WT_BECKE
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
c
c        Sort out the screening stuff
c
         if (rhotol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               rhotol(it)=1.0d-8
            else if (gaccu_num(0,it).eq.GACC_LOWMEDIUM) then
               rhotol(it)=1.0d-9
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               rhotol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_MEDIUMHIGH) then
               rhotol(it)=1.0d-11
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               rhotol(it)=1.0d-12
            else if (gaccu_num(0,it).eq.GACC_REF) then
               rhotol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               rhotol(it)=1.0d-10
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
         if (dentol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_LOWMEDIUM) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_MEDIUMHIGH) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               dentol(it)=1.0d-10
            else if (gaccu_num(0,it).eq.GACC_REF) then
               dentol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               dentol(it)=0.0d0
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
         if (wghtol(it).eq.dble(DFT_UNDEF)) then
            if (gaccu_num(0,it).eq.GACC_LOW) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_LOWMEDIUM) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_MEDIUM) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_MEDIUMHIGH) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_HIGH) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_REF) then
               wghtol(it)=1.0d-20
            else if (gaccu_num(0,it).eq.GACC_SG1) then
               wghtol(it)=1.0d-20
            else
               call caserr('CD_init DFT defaults messed up!')
            endif
         endif
c
c        Phase 1)
c
         do ig=1,ngtypes
            if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
               psitol(ig,it)=psitol(0,it)
            endif
            if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
               warntol(ig,it)=warntol(0,it)
            endif
            if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
               radm_num(ig,it)=radm_num(0,it)
            endif
            if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
               grid_scale(ig,it)=grid_scale(0,it)
            endif
            if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               grid_atom_radius(ig,it)=grid_atom_radius(0,it)
            endif
            if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               weight_atom_radius(ig,it)=weight_atom_radius(0,it)
            endif
            if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               prune_atom_radius(ig,it)=prune_atom_radius(0,it)
            endif
            if (screen_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
               screen_atom_radius(ig,it)=screen_atom_radius(0,it)
            endif
            do i = 1, maxradzn-1
               if (bnd_radzn(i,ig,it).eq.dble(DFT_UNDEF)) then
                  bnd_radzn(i,ig,it)=bnd_radzn(i,0,it)
               endif
            enddo
            do i = 1, maxradzn
               if (angpt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  angpt_radzn_num(i,ig,it)=angpt_radzn_num(i,0,it)
               endif
               if (phipt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  phipt_radzn_num(i,ig,it)=phipt_radzn_num(i,0,it)
               endif
               if (thetpt_radzn_num(i,ig,it).eq.DFT_UNDEF) then
                  thetpt_radzn_num(i,ig,it)=thetpt_radzn_num(i,0,it)
               endif
            enddo
            if (radzones_num(ig,it).eq.DFT_UNDEF) then
               radzones_num(ig,it)=radzones_num(0,it)
            endif
            if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
               ang_prune_scheme(ig,it)=ang_prune_scheme(0,it)
            endif
            if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
               rad_grid_scheme(ig,it)=rad_grid_scheme(0,it)
            endif
            if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
               ang_grid_scheme(ig,it)=ang_grid_scheme(0,it)
            endif
            if (radpt_num(ig,it).eq.DFT_UNDEF) then
               radpt_num(ig,it)=radpt_num(0,it)
            endif
            if (gaccu_num(ig,it).eq.DFT_UNDEF) then
               gaccu_num(ig,it)=gaccu_num(0,it)
            endif
         enddo
c
c        Phase 2)
c
         do ig=1,ngtypes
c
c           Default grid for non-atomic centres is that of carbon
c           See message in master.m
c
            atmnum = atom_num_by_grid(ig)
            if(atmnum.eq.0)atmnum = 6
            row = row_by_atomnum(atmnum)
c
            if (radpt_num(ig,it).eq.DFT_UNDEF) then
               radpt_num(ig,it)=radnpt_row(row)
            endif
            if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF) then
               angpt_radzn_num(1,ig,it)=angnpt_row(row)
            endif
c
            if (gaccu_num(ig,it).eq.GACC_LOW) then
c
c              Low accuracy grid (rel.err. 1.0d-4)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-6
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-4
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_low(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_low(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_LOWMEDIUM) then
c
c              Low-Medium accuracy grid (rel.err. 1.0d-4)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-7
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-4
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_lowmed(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_lowmed(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_MEDIUM) then
c
c              Medium accuracy grid (rel.err. 1.0d-6)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-7
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-6
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_medium(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_medium(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_MEDIUMHIGH) then
c
c              Medium accuracy grid (rel.err. 1.0d-6)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-10
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-6
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_medhigh(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_medhigh(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_HIGH) then
c
c              High accuracy grid (rel.err. 1.0d-8)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-13
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-8
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_high(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_high(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_REF) then
c
c              Reference accuracy grid
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=1.0d-20
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-12
               endif
               if (radm_num(ig,it).eq.dble(DFT_UNDEF)) then
                  radm_num(ig,it)=3.0d0
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  if (rad_scale_scheme.eq.SC_MK) then
                     grid_atom_radius(ig,it)=mkrad(atmnum,0,.false.)
                  else if (rad_scale_scheme.eq.SC_GAM1) then
                     grid_atom_radius(ig,it)=3.3d0*srad(atmnum)
                  else if (rad_scale_scheme.eq.SC_GAM2) then
                     grid_atom_radius(ig,it)=lograd(atmnum,0,.false.)
                  endif
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_MHL
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_MK
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=ang_ref(row)
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=rad_ref(row)
               endif
            else if (gaccu_num(ig,it).eq.GACC_SG1) then
c
c              SG1 (Standard grid #1)
c
               if (psitol(ig,it).eq.dble(DFT_UNDEF)) then
                  psitol(ig,it)=0.0d0
               endif
               if (warntol(ig,it).eq.dble(DFT_UNDEF)) then
                  warntol(ig,it)=1.0d-4
               endif
               if (grid_scale(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_scale(ig,it)=1.0d0
               endif
               if (grid_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  grid_atom_radius(ig,it)=sg1rad(atmnum,0,.false.)
               endif
               if (weight_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  weight_atom_radius(ig,it)=srad(atmnum)
               endif
               if (prune_atom_radius(ig,it).eq.dble(DFT_UNDEF)) then
                  prune_atom_radius(ig,it)=srad(atmnum)
               endif
               if (ang_prune_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_prune_scheme(ig,it)=AP_SG1
               endif
               if (rad_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  rad_grid_scheme(ig,it)=RG_SG1
               endif
               if (ang_grid_scheme(ig,it).eq.DFT_UNDEF) then
                  ang_grid_scheme(ig,it)=AG_LEB 
               endif
               if (angpt_radzn_num(1,ig,it).eq.DFT_UNDEF.and.
     +             radzones_num(ig,it).eq.DFT_UNDEF) then
                  radzones_num(ig,it)=1
                  angpt_radzn_num(1,ig,it)=302
               endif
               if (radpt_num(ig,it).eq.DFT_UNDEF) then
                  radpt_num(ig,it)=50
               endif
            else
               call caserr
     +              ('CD_defaults is messed up or never called!!!')
            endif
         enddo
c
c        Phase 3)
c
         do ig = 1, ngtypes
            if (ang_prune_scheme(ig,it).eq.AP_MHL) then
               call mhl_prune(ig,it)
            else if (ang_prune_scheme(ig,it).eq.AP_AUTO) then
c
c              initially use MHL pruning, the grid parameter will later
c              be modified in auto_prune
c
               call mhl_prune(ig,it)
            else if (ang_prune_scheme(ig,it).eq.AP_SG1) then
               call sg1_prune(ig,it,1)
            else if (ang_prune_scheme(ig,it).eq.AP_SG1a) then
               call sg1_prune(ig,it,2)
            else if (ang_prune_scheme(ig,it).ne.AP_RADZONE) then
               call caserr('grid_init: error angular pruning scheme')
            endif
         enddo
c
      enddo ! it
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine grid_sane_check(tag,iwr)
      implicit none
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_order_info)
INCLUDE(common/dft_xc)
c in
      integer tag,iwr

      integer lcent
      integer centre
      logical opg_root
c
c...  Check for some abnormalities
c
      do lcent=1, natoms
         centre = atom_tag(tag,lcent)
         if (gtype_num(lcent).eq.0) then
            if (centre.ne.0) then
               call caserr('atom with basis but no grid assigned')
            endif
         endif
         if (centre.eq.0) then
            if (gtype_num(lcent).ne.0) then
               if (opg_root()) then
                  write(iwr,600)lcent
               endif
            endif
         endif
      enddo
 600  format('WARNING: centre ',i4,' has a grid but no basis assigned')
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_batch(mxp,natm,ngridcent,latm,lo,hi,num_near,
     &                     atom_list,rpt,rwt,apts,atom_c,awpt,
     &                     ra2_comp,ra2_val,wt)
      implicit none
c
c     Description:
c
c     Collect a number of grid points and put them in a single batch.
c     All quantities are calculated where needed only.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
c
c     Input:
c
      integer mxp       ! maximum number of points in batch
      integer natm      ! the number of atoms in the system
      integer ngridcent ! the number of grid centres in the system
      integer latm      ! the current atom
      integer lo        ! the index of the first point
      integer hi        ! the index of the last point
      integer num_near
      integer atom_list(num_near)
      REAL rpt     ! the radius of the current shell
      REAL rwt     ! the weight of the radial point
      REAL apts(3,1202)
      REAL atom_c(max_atom,3)
      REAL awpt(1202)
c
c     Output:
c
      REAL ra2_comp(mxp,natm,3) ! difference between point and atom
      REAL ra2_val(mxp,natm,2)  ! distance between point and atom
      REAL wt(mxp)              ! the weight of the point
c
c     Local:
c
      integer ipt, iatm, lang, i
      REAL px, py, pz
cDEBUG
      integer iunk(2)
      REAL junk
      equivalence(junk,iunk)
      iunk(1) = -1
      iunk(2) = -1
c     call aclear_dp(ra2_comp(1,1,1),mxp*ngridcent,junk)
c     call aclear_dp(ra2_comp(1,1,2),mxp*ngridcent,junk)
c     call aclear_dp(ra2_comp(1,1,3),mxp*ngridcent,junk)
c     call aclear_dp(ra2_val(1,1,1),mxp*ngridcent,junk)
c     call aclear_dp(ra2_val(1,1,2),mxp*ngridcent,junk)
cDEBUG
c
c     Code:
c
      call aclear_dp(ra2_val(1,1,2),mxp*ngridcent,1.0d10)
      ipt=1
      do lang = lo, hi
c
c...     coords of pt
c
         px=apts(1,lang)*rpt+atom_c(latm,1)
         py=apts(2,lang)*rpt+atom_c(latm,2)
         pz=apts(3,lang)*rpt+atom_c(latm,3)
c
         wt(ipt)=awpt(lang)*rwt
c
c...     distances and atomic vectors
c
         do i=1,num_near
            iatm = atom_list(i)
c
            ra2_comp(ipt,iatm,1)=px - atom_c(iatm,1)
            ra2_comp(ipt,iatm,2)=py - atom_c(iatm,2)
            ra2_comp(ipt,iatm,3)=pz - atom_c(iatm,3)
            ra2_val(ipt,iatm,1) =
     &           ra2_comp(ipt,iatm,1)*ra2_comp(ipt,iatm,1)
     &         + ra2_comp(ipt,iatm,2)*ra2_comp(ipt,iatm,2)
     &         + ra2_comp(ipt,iatm,3)*ra2_comp(ipt,iatm,3)
c
            ra2_val(ipt,iatm,2)=sqrt(ra2_val(ipt,iatm,1))
         enddo
c
         ipt = ipt + 1
      enddo
      end
c
c-----------------------------------------------------------------------
c
      subroutine find_num_grid_centres
      implicit none
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
      integer i, ao_tag
      parameter(ao_tag=1)
c
c     In the quadrature we loop over the number of centres in a few
c     innerloops. This can be expensive, e.g. if there are a lot of
c     point charges. What we really want is to loop over all centres
c     that have a grid on them. Currently these centres are only those
c     centres which have AO-basis functions or have a grid type not 
c     equal 0. 
c
c     In the case BQ centres and other centres are mixed it gets 
c     complicated to keep a list of just the relevant atoms. So we set
c     ngridcentres to the last centre with a grid.
c
c     All this should be changed once we have the option to specify
c     a grid per centre.
c
      ngridcentres = 0
      do i = 1, natoms
         if (BL_get_atom_type(ao_tag,i).ne.0.or.
     +       gtype_num(i).ne.0) then
            ngridcentres = i
         endif
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_active_wght_atm(max_near_atms,near_atom_list,
     +                        num_near_atoms,latm,
     +                        ri,arad,rnear,natoms1,screen_sw,weights,
     +                        iwr)
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
      integer max_near_atms
      integer num_near_atoms, natoms1, near_atom_list(natoms1), latm
      integer iwr
      REAL ri, arad(natoms), rnear
      logical screen_sw
      integer weights
c
c     This subroutine constructs the list of neighbouring atoms.
c     The list is used in the weighting schemes.
c     While we are at it the distance to the closest neighbour is
c     registered also.
c
c     Although the current atom is stored on the same array it must
c     not be considered as part of the list.
c
c     num_near_atoms: On return the number of atoms on the near atom 
c                     list.
c     near_atom_list: On return the list of near atoms and appended to
c                     the list the label of the current atom.
c     latm          : The current atom.
c     ri            : The radius of the current angular shell.
c     arad          : The atomic radii.
c     rnear         : On return the distance to the closest atom
c                     neigbouring latm.
c     natoms1       : The number of atoms in the molecule.
c     screen_sw     : True if screening turned on.
c
c...  Local variables
c
      integer jatm
      REAL rtest, radl, radj, a
      logical use_screen_sw
c
c...  Code
c
      num_near_atoms=0
      rnear = 1.0d10
      use_screen_sw = screen_sw.and..not.(
     +                weights.eq.WT_BECKE.or.
     +                weights.eq.WT_MHL.or.
     +                weights.eq.WT_SSF)
      if (use_screen_sw) then
c
c       Screening turned on so record the nearest atoms only
c
        do jatm = 1,natoms1
          if (gtype_num(jatm).ne.0) then
            rtest = abs( dij(latm,jatm) - ri) 
            if (jatm .ne. latm) then
              if (arad(jatm) .gt. rtest) then
                if (num_near_atoms.ge.max_near_atms) then
                  write(iwr,*)'BLD_ACTIVE_WGHT_ATM:'
                  write(iwr,*)'num_near_atoms = ',num_near_atoms
                  write(iwr,*)'max_near_atms  = ',max_near_atms 
                  write(iwr,*)'Change max_near_atoms in ',
     +                        'common/dft_parameters and recompile'
                  call caserr(
     +             "actual number of near atoms exceeds code parameter")
                endif
                num_near_atoms = num_near_atoms+1
                near_atom_list(num_near_atoms) = jatm
              endif
              rnear = min(rnear,dij(latm,jatm))
cDEBUG
c           else
c             num_near_atoms = num_near_atoms+1
c             near_atom_list(num_near_atoms) = latm
cDEBUG
            endif
          endif
        enddo
c
c       last entry (not included in num_near_atoms) is home atom
c
        num_near_atoms = num_near_atoms+1
        near_atom_list(num_near_atoms) = latm
c
      else
c
c       Screening turned off so record all atoms
c
        if (natoms1.gt.max_near_atms) then
          write(iwr,*)'BLD_ACTIVE_WGHT_ATM:'
          write(iwr,*)'num_near_atoms = ',num_near_atoms
          write(iwr,*)'max_near_atms  = ',max_near_atms
          write(iwr,*)'Change max_near_atoms in ',
     +                'common/dft_parameters and recompile'
          call caserr("actual number of atoms exceeds code parameter")
        endif
        do jatm = 1, natoms1
          near_atom_list(jatm) = jatm
          if (gtype_num(jatm).ne.0.and.jatm.ne.latm) then
             rnear = min(rnear,dij(latm,jatm))
          endif
        enddo
        num_near_atoms = natoms1
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_active_bfn_atm(max_near_atms,near_atom_list,
     +                        num_near_atoms,latm,
     +                        ri,arad,natoms1,screen_sw,
     +                        iwr)
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
      integer max_near_atms
      integer num_near_atoms, natoms1, near_atom_list(natoms1), latm
      integer iwr
      REAL ri, arad(natoms), rnear
      logical screen_sw
      integer weights
c
c     This subroutine constructs the list of neighbouring atoms.
c     The list is used in the calculation of the basis functions.
c
c     num_near_atoms: On return the number of atoms on the near atom 
c                     list.
c     near_atom_list: On return the list of near atoms and appended to
c                     the list the label of the current atom.
c     latm          : The current atom.
c     ri            : The radius of the current angular shell.
c     arad          : The atomic radii.
c     natoms1       : The number of atoms in the molecule.
c     screen_sw     : True if screening turned on.
c
c...  Local variables
c
      integer jatm, ao_tag
      REAL rtest, radl, radj, a
c
c...  Functions
c
      integer BL_get_atom_type
      external BL_get_atom_type
c
c...  Code
c
      ao_tag = 1
      num_near_atoms=0
      if (screen_sw) then
c
c       Screening turned on so record the nearest atoms only
c
        do jatm = 1,natoms1
          if (BL_get_atom_type(ao_tag,jatm).ne.0) then
            rtest = abs( dij(latm,jatm) - ri) 
            if (arad(jatm) .gt. rtest) then
              if (num_near_atoms.ge.max_near_atms) then
                write(iwr,*)'BLD_ACTIVE_BFN_ATM:'
                write(iwr,*)'num_near_atoms = ',num_near_atoms
                write(iwr,*)'max_near_atms  = ',max_near_atms 
                write(iwr,*)'Change max_near_atoms in ',
     +                      'common/dft_parameters and recompile'
                call caserr(
     +           "actual number of near atoms exceeds code parameter")
              endif
              num_near_atoms = num_near_atoms+1
              near_atom_list(num_near_atoms) = jatm
            endif
          endif
        enddo
c
      else
c
c       Screening turned off so record all atoms
c
        if (natoms1.gt.max_near_atms) then
          write(iwr,*)'BLD_ACTIVE_BFN_ATM:'
          write(iwr,*)'num_near_atoms = ',num_near_atoms
          write(iwr,*)'max_near_atms  = ',max_near_atms
          write(iwr,*)'Change max_near_atoms in ',
     +                'common/dft_parameters and recompile'
          call caserr("actual number of atoms exceeds code parameter")
        endif
        do jatm = 1, natoms1
          near_atom_list(jatm) = jatm
        enddo
        num_near_atoms = natoms1
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine calc_max_active_atm(screen_sw,weights,arad,natoms,
     +           mx_active_atm)
      implicit none
c
c     Description:
c
c     This subroutine calculates the maximum number of active atoms
c     using the atom radii. An atom will be deemed active if their 
c     relative distance is less than the sum of their radii.
c
c     Input:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
c
      logical screen_sw            ! true if screening is to be used.
      integer natoms               ! the number of atoms
      REAL    arad(natoms)         ! the radii of the atoms
      integer weights              ! the weighting scheme
c
c     Output:
c
      integer mx_active_atm        ! the maximum number of active atoms
                                   ! at any time throughout the 
                                   ! quadrature
c
c     Local:
c
      integer i,j, num_active_atm
      logical use_screen_sw
      REAL disttol
      parameter(disttol=1.0d0)
c
c     Code:
c
      use_screen_sw = screen_sw.and..not.(
     +                weights.eq.WT_BECKE.or.
     +                weights.eq.WT_MHL.or.
     +                weights.eq.WT_SSF)
      if (use_screen_sw.and.1.eq.0) then
c...  disabled because mx_active_atm was wrong in hessian calc  (jvl,2010)
         mx_active_atm = 0
         do i=1,natoms
            if (gtype_num(i).ne.0) then
               num_active_atm = 0
               do j=1,natoms
                  if (gtype_num(j).ne.0.and.
     &                dij(i,j).le.disttol*(arad(i)+arad(j))) then
                     num_active_atm = num_active_atm + 1
                  endif
               enddo
               mx_active_atm = max(mx_active_atm,num_active_atm)
            endif
         enddo
      else
         mx_active_atm = natoms
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine calc_max_active_bfn(screen_sw,arad,first_bf,natoms,
     +           brad,nao,
     +           mx_active_bfn)
      implicit none
c
c     Description:
c
c     This subroutine calculates the maximum number of active basis
c     functions using the atom radii and basis function radii. 
c     For each atom A all its own basis functions will be deemed active
c     as well as all basis functions on all other atoms B for which 
c     the relative distance between the atoms A and B is less than 
c     the atom radius of A plus the basis function radius on B.
c
c     Input:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dij)
c
      logical screen_sw            ! true if screening is to be used.
      integer natoms               ! the number of atoms
      REAL    arad(natoms)         ! the radii of the atoms
      integer first_bf(natoms+1)   ! the first basis function on an
                                   ! atom
      integer nao                  ! the number of basis functions
      REAL    brad(nao)            ! the radii of the basis functions
c
c     Output:
c
      integer mx_active_bfn        ! the maximum number of active basis
                                   ! functions at any time throughout 
                                   ! the quadrature
c
c     Local:
c
      integer i,j,k, num_active_bfn
      REAL disttol
      parameter(disttol=1.0d0)
c
c     Code:
c
      if (screen_sw) then
         mx_active_bfn = 0
         do i=1,natoms
            num_active_bfn = 0
            do j=1,natoms
               do k=first_bf(j),first_bf(j+1)-1
                  if (dij(i,j).le.disttol*(arad(i)+brad(k))) then
                     num_active_bfn = num_active_bfn + 1
                  endif
               enddo
            enddo
            mx_active_bfn = max(mx_active_bfn,num_active_bfn)
         enddo
      else
         mx_active_bfn = nao
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_active_prt_atm(geometric_pert_sw,screen_sw,iwr,
     +              chf_pert_atms,npert,
     +              arad,latm,natoms,
     +              active_chf_pert,max_active_atm,n_active_chf_prt)
      implicit none
c
c     Description:
c
c     This subroutine constructs a list that contains a subset of all
c     perturbations passed to the CHF equations. In the case of external
c     field perturbations the subset will equal the set of all 
c     perturbations. In case of geometric perturbations the subset
c     will contain only those perturbations that are associated with
c     atoms that appear on the active basis function atoms list.
c
c     Input:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dij)
c
      logical geometric_pert_sw    ! true if we are considering
                                   ! geometric perturbations
      logical screen_sw            ! true if screening is to be used.
      integer iwr                  ! unit of standard out
      integer npert                ! the number of perturbations
      integer chf_pert_atms(npert) ! lists for each perturbation the
                                   ! associated atom
      integer natoms               ! the number of atoms
      REAL    arad(natoms)         ! the radii of the atoms
      integer latm                 ! the current grid atom
      integer max_active_atm
c
c     Output:
c
      integer active_chf_pert(3*max_active_atm) 
      integer n_active_chf_prt     ! the number of perturbations to be
                                   ! considered
c
c     Local:
c
      integer i,j
      logical otest(max_atom)
      REAL disttol
      parameter(disttol=1.25d0)
c
c     Code:
c
      if (geometric_pert_sw.and.screen_sw) then
         do i=1,natoms
            otest(i)=dij(i,latm).le.disttol*(arad(i)+arad(latm))
         enddo
         n_active_chf_prt = 0
         do 10 j = 1, npert
            do i = 1, natoms
               if (otest(i).and.i.eq.chf_pert_atms(j)) then
                  if (n_active_chf_prt.ge.3*max_active_atm) then
                     write(iwr,*)'BLD_ACTIVE_PRT_ATM:'
                     write(iwr,*)'n_active_chf_prt = ',n_active_chf_prt
                     write(iwr,*)'3*max_active_atm = ',3*max_active_atm
cjvl                     write(iwr,*)'Change max_active_atm in ',
cjvl     +                           'common/dft_parameters and recompile'
                     call caserr("actual number of perturbations "//
     +                           "exceeds code parameter")
                  endif
                  n_active_chf_prt = n_active_chf_prt + 1
                  active_chf_pert(n_active_chf_prt) = j
                  goto 10
               endif
            enddo
 10      continue
      else
         if (npert.gt.3*max_active_atm) then
            write(iwr,*)'BLD_ACTIVE_PRT_ATM:'
            write(iwr,*)'npert            = ',npert
            write(iwr,*)'3*max_active_atm = ',3*max_active_atm
cjvl            write(iwr,*)'Change max_active_atm in ',
cjvl     +                  'common/dft_parameters and recompile'
            call caserr(
     +       "actual number of perturbations exceeds code parameter")
         endif
         do i = 1, npert
            active_chf_pert(i) = i
         enddo
         n_active_chf_prt = npert
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine bld_active_prt_atm_old(geometric_pert_sw,screen_sw,iwr,
     +              max_active_atm,
     +              chf_pert_atms,npert,
     +              active_bfn_atms,n_active_bfn_atm,
     +              active_chf_pert,n_active_chf_prt)
      implicit none
c
c     Description:
c
c     This subroutine constructs a list that contains a subset of all
c     perturbations passed to the CHF equations. In the case of external
c     field perturbations the subset will equal the set of all 
c     perturbations. In case of geometric perturbations the subset
c     will contain only those perturbations that are associated with
c     atoms that appear on the active basis function atoms list.
c
c     Input:
c
      logical geometric_pert_sw    ! true if we are considering
                                   ! geometric perturbations
      logical screen_sw            ! true if screening is to be used.
      integer iwr                  ! unit of standard out
      integer max_active_atm       ! the maximum number of active atoms
      integer npert                ! the number of perturbations
      integer chf_pert_atms(npert) ! lists for each perturbation the
                                   ! associated atom
      integer n_active_bfn_atm     ! the number of atoms with active
                                   ! (ie, non-zero) basis functions
      integer active_bfn_atms(n_active_bfn_atm) ! list of active atoms
c
c     Output:
c
      integer active_chf_pert(3*max_active_atm) 
      integer n_active_chf_prt     ! the number of perturbations to be
                                   ! considered
c
c     Local:
c
      integer i,j
c
c     Code:
c
      if (geometric_pert_sw.and.screen_sw) then
         n_active_chf_prt = 0
         do 10 j = 1, npert
            do i = 1, n_active_bfn_atm
               if (active_bfn_atms(i).eq.chf_pert_atms(j)) then
                  if (n_active_chf_prt.ge.3*max_active_atm) then
                     write(iwr,*)'BLD_ACTIVE_PRT_ATM:'
                     write(iwr,*)'n_active_chf_prt = ',n_active_chf_prt
                     write(iwr,*)'3*max_active_atm = ',3*max_active_atm
cjvl                     write(iwr,*)'Change max_active_atm in ',
cjvl     +                           'common/dft_parameters and recompile'
                     call caserr("actual number of perturbations "//
     +                           "exceeds code parameter")
                  endif
                  n_active_chf_prt = n_active_chf_prt + 1
                  active_chf_pert(n_active_chf_prt) = j
                  goto 10
               endif
            enddo
 10      continue
      else
         if (npert.gt.3*max_active_atm) then
            write(iwr,*)'BLD_ACTIVE_PRT_ATM:'
            write(iwr,*)'npert            = ',npert
            write(iwr,*)'3*max_active_atm = ',3*max_active_atm
cjvl            write(iwr,*)'Change max_active_atm in ',
cjvl     +                  'common/dft_parameters and recompile'
            call caserr(
     +       "actual number of perturbations exceeds code parameter")
         endif
         do i = 1, npert
            active_chf_pert(i) = i
         enddo
         n_active_chf_prt = npert
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine bas_val(tag,nao,
     &     ra2_comp,ra2_val,
     &     bfn_val,bfng_val,bfn_hess,bfn_3rd,
     &     npts,mxp,screen_sw,first_bf,
     &     active_bfn_atms,n_active_bfn_atm,
     &     active_bfn_list,active_bfn_indx,max_active_bfn,n_active_bfn,
     &     psitol, del_psi, bfn_radii, ider, odb, iout)

C     ******************************************************************
C     * Description:					
C     * Calculate value of basis functions at at point		
c     *
c     * Computation of basis values and derivatives
c     * this version assumes we only want significant basis values
c     *
c     * screening now applied using basis function radii
c     *
C     ******************************************************************
      implicit none
C     *****************************************************************
C     *Declarations						
C     *								
C     *Parameters
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basder)
INCLUDE(common/dft_physical_constants)
c
C     *In variables
c				
INCLUDE(common/dft_numbers)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_order_info)
      integer mxp
      integer first_bf(max_atom+1)
      integer n_active_bfn_atm
      integer active_bfn_atms(n_active_bfn_atm)
      REAL ra2_comp(mxp,natoms,3)
      REAL ra2_val(mxp,natoms,2)
      integer tag, nao
      integer ider, iout
c
C     *Out variables
c						
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      REAL bfn_3rd(mxp,nao,10)
      integer npts
      logical screen_sw
c
c     active_bfn_list: contains for every active basis function its
c                      number in the total basis.
c     active_bfn_indx: contains for every active basis function the
c                      position of its atom in the active basis function
c                      atom list (i.e. active_bfn_atms).
      integer max_active_bfn
      integer active_bfn_list(max_active_bfn)
      integer active_bfn_indx(max_active_bfn)
      integer n_active_bfn
      REAL psitol
      REAL bfn_radii(*)
c
c     *In/Out variables
c
      REAL del_psi
c
C     *Local variables
c      
      integer lcent,lshl,lprm,lhyb
      integer icent
      integer num_cent
      integer bpos,bfn,istart,icount,iloc,ilpos
      integer nshells,centre,nprm,ploc,ploc_save,ll,lh
      REAL tol
      REAL xd(500),yd(500),zd(500),ra2(500),ra1(500)
      REAL xx,yy,zz,xy,xz,yz
_IF(single)
      integer isamin, isamax
_ELSE
      integer idamin, idamax 
_ENDIF
      integer i
      REAL rmin
      REAL bfn_rmax

      REAL alp,cc,expo,minalp,minexp,alpexp
      integer px,py,pz,dxx,dyy,dzz,dxy,dxz,dyz
      integer ipt

      REAL g,dg,ddg,dgu,ddgu
      REAL gx,gy,gz
      REAL gxr,gyr,gzr
      REAL dgx,dgy,dgz
      REAL dgxx,dgyy,dgzz,dgxy,dgxz,dgyz
      REAL dgxxr,dgyyr,dgzzr,dgxyr,dgxzr,dgyzr
      REAL dgxxx,dgxxy,dgxxz,dgxyy,dgyyy
      REAL dgyyz,dgxzz,dgyzz,dgzzz,dgxyz
      REAL ddgx,ddgy,ddgz
      REAL ddgxxx,ddgxxy,ddgxxz,ddgxyy,ddgxyz,ddgxzz,ddgyyy
      REAL ddgyyz,ddgyzz,ddgzzz
      REAL ddgxxxx,ddgxxxy,ddgxxxz,ddgxxyy,ddgxxyz,ddgxxzz
      REAL ddgxyyy,ddgxyyz,ddgxyzz,ddgxzzz,ddgyyyy,ddgyyyz
      REAL ddgyyzz,ddgyzzz,ddgzzzz

      REAL g3
      REAL alp2, alp4, alp6, alp8, alp10, alp12, alp24, alp36
      REAL x0d1, x0d2, x0d3
      REAL x1d0, x1d1, x1d2, x1d3
      REAL x2d0, x2d1, x2d2, x2d3
      REAL y0d1, y0d2, y0d3
      REAL y1d0, y1d1, y1d2, y1d3
      REAL y2d0, y2d1, y2d2, y2d3
      REAL z0d1, z0d2, z0d3
      REAL z1d0, z1d1, z1d2, z1d3
      REAL z2d0, z2d1, z2d2, z2d3

      integer nftab(5,5)

      integer k,l,m, ilm
      integer kq(25), lq(25), mq(25)

      REAL small
      parameter (small=1.0d-20)
      REAL x,y,z,e

      REAL ffac(25)

      logical odb

C     *End declarations							      *
C     ******************************************************************

      data nftab/ 1, 0, 0, 0, 0,
     p            4, 3, 0, 0, 0,
     d           10, 9, 6, 0, 0,
     f           20,19,16,10, 0,
     g           35,34,31,25,15/
      data kq /3,0,0,2,2,1,0,1,0,1,
     g         4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data lq /0,3,0,1,0,2,2,0,1,1,
     g         0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data mq /0,0,3,0,1,0,1,2,2,1,
     g         0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/

      ffac(1)=1.0d0
      ffac(2)=1.0d0
      ffac(3)=1.0d0

      ffac(4)=sqrt5
      ffac(5)=sqrt5
      ffac(6)=sqrt5
      ffac(7)=sqrt5
      ffac(8)=sqrt5
      ffac(9)=sqrt5

      ffac(10)=sqrt5*sqrt3

      ffac(11) = 1.0d0
      ffac(12) = 1.0d0
      ffac(13) = 1.0d0

      ffac(14) = sqrt7
      ffac(15) = sqrt7
      ffac(16) = sqrt7
      ffac(17) = sqrt7
      ffac(18) = sqrt7
      ffac(19) = sqrt7

      ffac(20) = sqrt7*sqrt5/sqrt3
      ffac(21) = sqrt7*sqrt5/sqrt3
      ffac(22) = sqrt7*sqrt5/sqrt3

      ffac(23) = sqrt7*sqrt5
      ffac(24) = sqrt7*sqrt5
      ffac(25) = sqrt7*sqrt5
C     
C     Clean arrays
C     
      call aclear_dp(bfn_val,mxp*totbfn(tag),0.0d0)
      if(ider.ge.1)call aclear_dp(bfng_val,3*mxp*totbfn(tag),0.0d0)
      if(ider.ge.2)call aclear_dp(bfn_hess,6*mxp*totbfn(tag),0.0d0)
      if(ider.ge.3)call aclear_dp(bfn_3rd,10*mxp*totbfn(tag),0.0d0)


      n_active_bfn = 0

      tol   = -36.9d0    ! dlog(1.0d-16)
c     num_cent=ngridcentres
      num_cent=n_active_bfn_atm
c     
c  general version, with no bfn value tests
c     
c  index of next free location in table of selected functios
      istart = 1                

c  index of next location in complete table of functions
      bfn    = 1               

c     do lcent=1,num_cent
      do icent=1,num_cent
         lcent=active_bfn_atms(icent)
c     
c        store vectors for points wrt to the centre
c     
         centre = BL_get_atom_type(tag,lcent)
         if (centre.eq.0) goto 100
         bfn = first_bf(lcent)
         nshells = num_shl(tag,centre)
         do ipt=1,npts
            ra2(ipt)    = ra2_val(ipt,lcent,1)
            ra1(ipt)    = ra2_val(ipt,lcent,2)
            xd(ipt)     = ra2_comp(ipt,lcent,1)
            yd(ipt)     = ra2_comp(ipt,lcent,2)
            zd(ipt)     = ra2_comp(ipt,lcent,3)
         enddo
c
c        shortest contact from nucleus to this batch
c
         rmin = ra1(idamin(npts,ra1,1)) 

         ploc   = 0
         do lshl=1,nshells
            nprm   = nprim(tag,centre,lshl)
            lh     = hybrid(tag,centre,lshl)
            ll     = angmom(tag,centre,lshl)
            bpos=nftab(lh,ll)   ! size of shell
c
c           most diffuse fn in this shell
c
            bfn_rmax = bfn_radii(bfn-1+idamax(bpos,bfn_radii(bfn),1)) 

            if(screen_sw .and. rmin .gt. bfn_rmax) then
c
c              we can reject all points for the whole shell 
c     
               ploc = ploc + nprm
               del_psi = del_psi + dble(bpos)

c              write(6,*)'reject shell',lcent,lshl,rmin, bfn_rmax

            else

               minalp = alpha(tag,centre,ploc+
     &                        idamin(nprm,alpha(tag,centre,ploc+1),1))

               ploc_save = ploc

               do ipt=1,npts 

                  ploc=ploc_save
                  minexp = minalp*ra2(ipt)
c
c                 if this point is too distant the whole shell will be
c                 skipped
c
                  if( .not. screen_sw .or.
     &                (ra1(ipt) .le. bfn_rmax)  ) then 

                     xx      = xd(ipt)*xd(ipt)
                     xy      = xd(ipt)*yd(ipt)
                     xz      = xd(ipt)*zd(ipt)
                     yy      = yd(ipt)*yd(ipt)
                     yz      = yd(ipt)*zd(ipt)
                     zz      = zd(ipt)*zd(ipt)

C                    loop over primitives in shell

                     do 200 lprm=1,nprm

                        ploc    = ploc + 1
                        alp     = alpha(tag,centre,ploc)
                        alpexp  = -alp*ra2(ipt)
c                       if (minexp-alpexp.lt.tol) goto 200
                        expo    = exp(alpexp)
                        dgu     = -2.0d0*alp
                        ddgu    = dgu*dgu
C     
C                       loop over hybrid angular momentum numbers 
c                       e.g. sp shells
c     
                        iloc = istart
                        do lhyb=lh,ll
C     
C                          Basis function values
C     
                           cc      = cont_coeff(tag,centre,ploc,lhyb)
                           g       = expo*cc
                           dg      = dgu*g
                           ddg     = ddgu*g
c
                           if(lhyb.eq.1) then
c
c                             S functions
c

                              bfn_val(ipt,iloc) = bfn_val(ipt,iloc) + g

                              if(ider.ge.1)then
                                 bfng_val(ipt,iloc,1) 
     &                           = bfng_val(ipt,iloc,1) + dg*xd(ipt)
                                 bfng_val(ipt,iloc,2) 
     &                           = bfng_val(ipt,iloc,2) + dg*yd(ipt)
                                 bfng_val(ipt,iloc,3) 
     &                           = bfng_val(ipt,iloc,3) + dg*zd(ipt)
                              endif

                              if(ider.ge.2)then
                                 ddgx = ddg*xd(ipt)
                                 ddgy = ddg*yd(ipt)
                                 ddgz = ddg*zd(ipt)
                                 bfn_hess(ipt,iloc,hxx) 
     &                           = bfn_hess(ipt,iloc,hxx) 
     &                             + dg + xd(ipt)*ddgx
                                 bfn_hess(ipt,iloc,hyy) 
     &                           = bfn_hess(ipt,iloc,hyy) 
     &                             + dg + yd(ipt)*ddgy
                                 bfn_hess(ipt,iloc,hzz) 
     &                           = bfn_hess(ipt,iloc,hzz) 
     &                             + dg + zd(ipt)*ddgz
                                 bfn_hess(ipt,iloc,hxy) 
     &                           = bfn_hess(ipt,iloc,hxy) 
     &                             + xd(ipt)*ddgy
                                 bfn_hess(ipt,iloc,hxz) 
     &                           = bfn_hess(ipt,iloc,hxz) 
     &                             + xd(ipt)*ddgz
                                 bfn_hess(ipt,iloc,hyz) 
     &                           = bfn_hess(ipt,iloc,hyz) 
     &                             + yd(ipt)*ddgz
                              endif

                              if (ider.ge.3) then
                                 alp2 = -2.0d0*alp
                                 alp4 = alp2*alp2
                                 alp8 = alp4*alp2
                                 alp12 = 3.0d0*alp4
                                 bfn_3rd(ipt,iloc,txxx)
     &                           = bfn_3rd(ipt,iloc,txxx)
     &                             + (alp12+alp8*xx)*xd(ipt)*g
                                 bfn_3rd(ipt,iloc,tyyy)
     &                           = bfn_3rd(ipt,iloc,tyyy)
     &                             + (alp12+alp8*yy)*yd(ipt)*g
                                 bfn_3rd(ipt,iloc,tzzz)
     &                           = bfn_3rd(ipt,iloc,tzzz)
     &                             + (alp12+alp8*zz)*zd(ipt)*g
                                 bfn_3rd(ipt,iloc,txxy)
     &                           = bfn_3rd(ipt,iloc,txxy)
     &                             + (alp4+alp8*xx)*yd(ipt)*g
                                 bfn_3rd(ipt,iloc,txxz)
     &                           = bfn_3rd(ipt,iloc,txxz)
     &                             + (alp4+alp8*xx)*zd(ipt)*g
                                 bfn_3rd(ipt,iloc,txyy)
     &                           = bfn_3rd(ipt,iloc,txyy)
     &                             + (alp4+alp8*yy)*xd(ipt)*g
                                 bfn_3rd(ipt,iloc,tyyz)
     &                           = bfn_3rd(ipt,iloc,tyyz)
     &                             + (alp4+alp8*yy)*zd(ipt)*g
                                 bfn_3rd(ipt,iloc,txzz)
     &                           = bfn_3rd(ipt,iloc,txzz)
     &                             + (alp4+alp8*zz)*xd(ipt)*g
                                 bfn_3rd(ipt,iloc,tyzz)
     &                           = bfn_3rd(ipt,iloc,tyzz)
     &                             + (alp4+alp8*zz)*yd(ipt)*g
                                 bfn_3rd(ipt,iloc,txyz)
     &                           = bfn_3rd(ipt,iloc,txyz)
     &                             + alp8*xy*zd(ipt)*g
                              endif
                              bpos=1
c
                           else if(lhyb.eq.2) then
c
c                             P functions
c
                              px=iloc
                              py=iloc+1
                              pz=iloc+2
                              bfn_val(ipt,px) 
     &                        = bfn_val(ipt,px) + g*xd(ipt)
                              bfn_val(ipt,py) 
     &                        = bfn_val(ipt,py) + g*yd(ipt)
                              bfn_val(ipt,pz) 
     &                        = bfn_val(ipt,pz) + g*zd(ipt)

                              if(ider.ge.1)then
                                 bfng_val(ipt,px,1) 
     &                           = bfng_val(ipt,px,1) + g + xx*dg
                                 bfng_val(ipt,px,2) 
     &                           = bfng_val(ipt,px,2) +     xy*dg
                                 bfng_val(ipt,px,3) 
     &                           = bfng_val(ipt,px,3) +     xz*dg
                                 bfng_val(ipt,py,1) 
     &                           = bfng_val(ipt,py,1) +     xy*dg
                                 bfng_val(ipt,py,2) 
     &                           = bfng_val(ipt,py,2) + g + yy*dg
                                 bfng_val(ipt,py,3) 
     &                           = bfng_val(ipt,py,3) +     yz*dg
                                 bfng_val(ipt,pz,1) 
     &                           = bfng_val(ipt,pz,1) +     xz*dg
                                 bfng_val(ipt,pz,2) 
     &                           = bfng_val(ipt,pz,2) +     yz*dg
                                 bfng_val(ipt,pz,3) 
     &                           = bfng_val(ipt,pz,3) + g + zz*dg
                              endif

                              if(ider.ge.2)then
c     
                                 dgx  = dg*xd(ipt)
                                 dgy  = dg*yd(ipt)
                                 dgz  = dg*zd(ipt)
                                 ddgx = ddg*xd(ipt)
                                 ddgy = ddg*yd(ipt)
                                 ddgz = ddg*zd(ipt)
                                 ddgxxx = ddgx*xx
                                 ddgxxy = ddgx*xy
                                 ddgxxz = ddgx*xz
                                 ddgxyy = ddgx*yy
                                 ddgxyz = ddgx*yz
                                 ddgxzz = ddgx*zz
                                 ddgyyy = ddgy*yy
                                 ddgyyz = ddgy*yz
                                 ddgyzz = ddgy*zz
                                 ddgzzz = ddgz*zz
C     
                                 bfn_hess(ipt,px,hxx) 
     &                           = bfn_hess(ipt,px,hxx) 
     &                             + 3.d0*dgx + ddgxxx
                                 bfn_hess(ipt,px,hyy) 
     &                           = bfn_hess(ipt,px,hyy) 
     &                             +      dgx + ddgxyy
                                 bfn_hess(ipt,px,hzz) 
     &                           = bfn_hess(ipt,px,hzz) 
     &                             +      dgx + ddgxzz
                                 bfn_hess(ipt,px,hxy) 
     &                           = bfn_hess(ipt,px,hxy) 
     &                             +      dgy + ddgxxy
                                 bfn_hess(ipt,px,hxz) 
     &                           = bfn_hess(ipt,px,hxz) 
     &                             +      dgz + ddgxxz
                                 bfn_hess(ipt,px,hyz) 
     &                           = bfn_hess(ipt,px,hyz) 
     &                             +            ddgxyz
C     
                                 bfn_hess(ipt,py,hxx) 
     &                           = bfn_hess(ipt,py,hxx) 
     &                             +      dgy + ddgxxy
                                 bfn_hess(ipt,py,hyy) 
     &                           = bfn_hess(ipt,py,hyy) 
     &                             + 3.d0*dgy + ddgyyy
                                 bfn_hess(ipt,py,hzz) 
     &                           = bfn_hess(ipt,py,hzz) 
     &                             +      dgy + ddgyzz
                                 bfn_hess(ipt,py,hxy) 
     &                           = bfn_hess(ipt,py,hxy) 
     &                             +      dgx + ddgxyy
                                 bfn_hess(ipt,py,hxz) 
     &                           = bfn_hess(ipt,py,hxz) 
     &                             +            ddgxyz
                                 bfn_hess(ipt,py,hyz) 
     &                           = bfn_hess(ipt,py,hyz) 
     &                             +      dgz + ddgyyz
C     
                                 bfn_hess(ipt,pz,hxx) 
     &                           = bfn_hess(ipt,pz,hxx) 
     &                             +      dgz + ddgxxz
                                 bfn_hess(ipt,pz,hyy) 
     &                           = bfn_hess(ipt,pz,hyy) 
     &                             +      dgz + ddgyyz
                                 bfn_hess(ipt,pz,hzz) 
     &                           = bfn_hess(ipt,pz,hzz) 
     &                             + 3.d0*dgz + ddgzzz
                                 bfn_hess(ipt,pz,hxy) 
     &                           = bfn_hess(ipt,pz,hxy) 
     &                             +            ddgxyz
                                 bfn_hess(ipt,pz,hxz) 
     &                           = bfn_hess(ipt,pz,hxz) 
     &                             +      dgx + ddgxzz
                                 bfn_hess(ipt,pz,hyz) 
     &                           = bfn_hess(ipt,pz,hyz) 
     &                             +      dgy + ddgyzz
                              endif 

                              if (ider.ge.3) then
                                 alp2 = -2.0d0*alp
                                 alp4 = alp2*alp2
                                 alp6 = -6.0d0*alp
                                 alp8 = alp4*alp2
                                 alp12 = 3.0d0*alp4
                                 alp24 = 6.0d0*alp4
                                 bfn_3rd(ipt,px,txxx) 
     &                           = bfn_3rd(ipt,px,txxx) 
     &                             + (alp6+alp24*xx+alp8*xx*xx)*g
                                 bfn_3rd(ipt,px,tyyy) 
     &                           = bfn_3rd(ipt,px,tyyy) 
     &                             + (alp12+alp8*yy)*xy*g
                                 bfn_3rd(ipt,px,tzzz) 
     &                           = bfn_3rd(ipt,px,tzzz) 
     &                             + (alp12+alp8*zz)*xz*g
                                 bfn_3rd(ipt,px,txxy) 
     &                           = bfn_3rd(ipt,px,txxy) 
     &                             + (alp12+alp8*xx)*xy*g
                                 bfn_3rd(ipt,px,txxz) 
     &                           = bfn_3rd(ipt,px,txxz) 
     &                             + (alp12+alp8*xx)*xz*g
                                 bfn_3rd(ipt,px,txyy) 
     &                           = bfn_3rd(ipt,px,txyy) 
     &                             + (alp2+alp4*(xx+yy)+alp8*xx*yy)*g
                                 bfn_3rd(ipt,px,txzz) 
     &                           = bfn_3rd(ipt,px,txzz) 
     &                             + (alp2+alp4*(xx+zz)+alp8*xx*zz)*g
                                 bfn_3rd(ipt,px,txyz) 
     &                           = bfn_3rd(ipt,px,txyz) 
     &                             + (alp4+alp8*xx)*yz*g
                                 bfn_3rd(ipt,px,tyyz) 
     &                           = bfn_3rd(ipt,px,tyyz) 
     &                             + (alp4+alp8*yy)*xz*g
                                 bfn_3rd(ipt,px,tyzz) 
     &                           = bfn_3rd(ipt,px,tyzz) 
     &                             + (alp4+alp8*zz)*xy*g
c
                                 bfn_3rd(ipt,py,txxx) 
     &                           = bfn_3rd(ipt,py,txxx) 
     &                             + (alp12+alp8*xx)*xy*g
                                 bfn_3rd(ipt,py,tyyy) 
     &                           = bfn_3rd(ipt,py,tyyy) 
     &                             + (alp6+alp24*yy+alp8*yy*yy)*g
                                 bfn_3rd(ipt,py,tzzz) 
     &                           = bfn_3rd(ipt,py,tzzz) 
     &                             + (alp12+alp8*zz)*yz*g
                                 bfn_3rd(ipt,py,txyy) 
     &                           = bfn_3rd(ipt,py,txyy) 
     &                             + (alp12+alp8*yy)*xy*g
                                 bfn_3rd(ipt,py,tyyz) 
     &                           = bfn_3rd(ipt,py,tyyz) 
     &                             + (alp12+alp8*yy)*yz*g
                                 bfn_3rd(ipt,py,txxy) 
     &                           = bfn_3rd(ipt,py,txxy) 
     &                             + (alp2+alp4*(yy+xx)+alp8*yy*xx)*g
                                 bfn_3rd(ipt,py,tyzz) 
     &                           = bfn_3rd(ipt,py,tyzz) 
     &                             + (alp2+alp4*(yy+zz)+alp8*yy*zz)*g
                                 bfn_3rd(ipt,py,txyz) 
     &                           = bfn_3rd(ipt,py,txyz) 
     &                             + (alp4+alp8*yy)*xz*g
                                 bfn_3rd(ipt,py,txxz) 
     &                           = bfn_3rd(ipt,py,txxz) 
     &                             + (alp4+alp8*xx)*yz*g
                                 bfn_3rd(ipt,py,txzz) 
     &                           = bfn_3rd(ipt,py,txzz) 
     &                             + (alp4+alp8*zz)*xy*g
c
                                 bfn_3rd(ipt,pz,txxx) 
     &                           = bfn_3rd(ipt,pz,txxx) 
     &                             + (alp12+alp8*xx)*xz*g
                                 bfn_3rd(ipt,pz,tyyy) 
     &                           = bfn_3rd(ipt,pz,tyyy) 
     &                             + (alp12+alp8*yy)*yz*g
                                 bfn_3rd(ipt,pz,tzzz) 
     &                           = bfn_3rd(ipt,pz,tzzz) 
     &                             + (alp6+alp24*zz+alp8*zz*zz)*g
                                 bfn_3rd(ipt,pz,txzz) 
     &                           = bfn_3rd(ipt,pz,txzz) 
     &                             + (alp12+alp8*zz)*xz*g
                                 bfn_3rd(ipt,pz,tyzz) 
     &                           = bfn_3rd(ipt,pz,tyzz) 
     &                             + (alp12+alp8*zz)*yz*g
                                 bfn_3rd(ipt,pz,txxz) 
     &                           = bfn_3rd(ipt,pz,txxz) 
     &                             + (alp2+alp4*(zz+xx)+alp8*zz*xx)*g
                                 bfn_3rd(ipt,pz,tyyz) 
     &                           = bfn_3rd(ipt,pz,tyyz) 
     &                             + (alp2+alp4*(zz+yy)+alp8*zz*yy)*g
                                 bfn_3rd(ipt,pz,txyz) 
     &                           = bfn_3rd(ipt,pz,txyz) 
     &                             + (alp4+alp8*zz)*xy*g
                                 bfn_3rd(ipt,pz,txxy) 
     &                           = bfn_3rd(ipt,pz,txxy) 
     &                             + (alp4+alp8*xx)*yz*g
                                 bfn_3rd(ipt,pz,txyy) 
     &                           = bfn_3rd(ipt,pz,txyy) 
     &                             + (alp4+alp8*yy)*xz*g
                              endif
                              bpos = 3
c
                           else if(lhyb.eq.3) then
c     
c                             D functions
c
                              dxx=iloc
                              dyy=iloc+1
                              dzz=iloc+2
                              dxy=iloc+3
                              dxz=iloc+4
                              dyz=iloc+5

                              bfn_val(ipt,dxx) = bfn_val(ipt,dxx) + g*xx
                              bfn_val(ipt,dyy) = bfn_val(ipt,dyy) + g*yy
                              bfn_val(ipt,dzz) = bfn_val(ipt,dzz) + g*zz
                              bfn_val(ipt,dxy) = bfn_val(ipt,dxy) 
     &                                         + g*xy*sqrt3
                              bfn_val(ipt,dxz) = bfn_val(ipt,dxz) 
     &                                         + g*xz*sqrt3
                              bfn_val(ipt,dyz) = bfn_val(ipt,dyz) 
     &                                         + g*yz*sqrt3

                              if(ider.ge.1)then

                                 gx = g*xd(ipt)
                                 gy = g*yd(ipt)
                                 gz = g*zd(ipt)
                                 gxr = gx*sqrt3
                                 gyr = gy*sqrt3
                                 gzr = gz*sqrt3
                                 dgxx  = dg*xx
                                 dgyy  = dg*yy
                                 dgzz  = dg*zz
                                 dgxy  = dg*xy
                                 dgxz  = dg*xz
                                 dgyz  = dg*yz

                                 dgxxr = dgxx*sqrt3
                                 dgyyr = dgyy*sqrt3
                                 dgzzr = dgzz*sqrt3
                                 dgxyr = dgxy*sqrt3
                                 dgxzr = dgxz*sqrt3
                                 dgyzr = dgyz*sqrt3

                                 dgxxx = dgxx*xd(ipt)
                                 dgxxy = dgxx*yd(ipt)
                                 dgxxz = dgxx*zd(ipt)
                                 dgxyy = dgxy*yd(ipt)
                                 dgyyy = dgyy*yd(ipt)
                                 dgyyz = dgyy*zd(ipt)
                                 dgxzz = dgxz*zd(ipt)
                                 dgyzz = dgyz*zd(ipt)
                                 dgzzz = dgzz*zd(ipt)
                                 dgxyz = dgxy*zd(ipt)


                                 bfng_val(ipt,dxx,1) 
     &                           = bfng_val(ipt,dxx,1) + 2.d0*gx + dgxxx
                                 bfng_val(ipt,dxx,2) 
     &                           = bfng_val(ipt,dxx,2) +           dgxxy
                                 bfng_val(ipt,dxx,3) 
     &                           = bfng_val(ipt,dxx,3) +           dgxxz
                                 bfng_val(ipt,dyy,1)
     &                           = bfng_val(ipt,dyy,1) +           dgxyy
                                 bfng_val(ipt,dyy,2) 
     &                           = bfng_val(ipt,dyy,2) + 2.d0*gy + dgyyy
                                 bfng_val(ipt,dyy,3) 
     &                           = bfng_val(ipt,dyy,3) +           dgyyz
                                 bfng_val(ipt,dzz,1) 
     &                           = bfng_val(ipt,dzz,1) +           dgxzz
                                 bfng_val(ipt,dzz,2) 
     &                           = bfng_val(ipt,dzz,2) +           dgyzz
                                 bfng_val(ipt,dzz,3) 
     &                           = bfng_val(ipt,dzz,3) + 2.d0*gz + dgzzz
                                 bfng_val(ipt,dxy,1) 
     &                           = bfng_val(ipt,dxy,1) 
     &                             +     gyr + dgxxy*sqrt3
                                 bfng_val(ipt,dxy,2) 
     &                           = bfng_val(ipt,dxy,2) 
     &                             +     gxr + dgxyy*sqrt3
                                 bfng_val(ipt,dxy,3) 
     &                           = bfng_val(ipt,dxy,3) 
     &                             +           dgxyz*sqrt3
                                 bfng_val(ipt,dxz,1) 
     &                           = bfng_val(ipt,dxz,1) 
     &                             +     gzr + dgxxz*sqrt3
                                 bfng_val(ipt,dxz,2) 
     &                           = bfng_val(ipt,dxz,2) 
     &                             +           dgxyz*sqrt3
                                 bfng_val(ipt,dxz,3) 
     &                           = bfng_val(ipt,dxz,3) 
     &                             +     gxr + dgxzz*sqrt3
                                 bfng_val(ipt,dyz,1) 
     &                           = bfng_val(ipt,dyz,1) 
     &                             +           dgxyz*sqrt3
                                 bfng_val(ipt,dyz,2) 
     &                           = bfng_val(ipt,dyz,2) 
     &                             +     gzr + dgyyz*sqrt3
                                 bfng_val(ipt,dyz,3) 
     &                           = bfng_val(ipt,dyz,3) 
     &                             +     gyr + dgyzz*sqrt3

                              endif

                              if(ider.ge.2)then

                                 ddgxxxx = ddg*xx*xx
                                 ddgxxxy = ddg*xx*xy
                                 ddgxxxz = ddg*xx*xz
                                 ddgxxyy = ddg*xx*yy
                                 ddgxxyz = ddg*xx*yz
                                 ddgxxzz = ddg*xx*zz
                                 ddgxyyy = ddg*xy*yy
                                 ddgxyyz = ddg*xy*yz
                                 ddgxyzz = ddg*xy*zz
                                 ddgxzzz = ddg*xz*zz
                                 ddgyyyy = ddg*yy*yy
                                 ddgyyyz = ddg*yy*yz
                                 ddgyyzz = ddg*yy*zz
                                 ddgyzzz = ddg*yz*zz
                                 ddgzzzz = ddg*zz*zz
                                 bfn_hess(ipt,dxx,hxx) 
     &                           = bfn_hess(ipt,dxx,hxx)
     &                             +g+g+5.d0*dgxx + ddgxxxx
                                 bfn_hess(ipt,dxx,hyy) 
     &                           = bfn_hess(ipt,dxx,hyy)
     &                             +        dgxx + ddgxxyy
                                 bfn_hess(ipt,dxx,hzz) 
     &                           = bfn_hess(ipt,dxx,hzz)
     &                             +        dgxx + ddgxxzz
                                 bfn_hess(ipt,dxx,hxy) 
     &                           = bfn_hess(ipt,dxx,hxy)
     &                             +   dgxy+dgxy + ddgxxxy
                                 bfn_hess(ipt,dxx,hyz) 
     &                           = bfn_hess(ipt,dxx,hyz)
     &                             +               ddgxxyz
                                 bfn_hess(ipt,dxx,hxz) 
     &                           = bfn_hess(ipt,dxx,hxz)
     &                             +   dgxz+dgxz + ddgxxxz
C     
                                 bfn_hess(ipt,dyy,hxx) 
     &                           = bfn_hess(ipt,dyy,hxx)
     &                             +        dgyy + ddgxxyy
                                 bfn_hess(ipt,dyy,hyy) 
     &                           = bfn_hess(ipt,dyy,hyy)
     &                             +g+g+5.d0*dgyy + ddgyyyy
                                 bfn_hess(ipt,dyy,hzz) 
     &                           = bfn_hess(ipt,dyy,hzz)
     &                             +        dgyy + ddgyyzz
                                 bfn_hess(ipt,dyy,hxy) 
     &                           = bfn_hess(ipt,dyy,hxy)
     &                             +   dgxy+dgxy + ddgxyyy
                                 bfn_hess(ipt,dyy,hyz) 
     &                           = bfn_hess(ipt,dyy,hyz)
     &                             +   dgyz+dgyz + ddgyyyz
                                 bfn_hess(ipt,dyy,hxz) 
     &                           = bfn_hess(ipt,dyy,hxz)
     &                             +               ddgxyyz
c     
                                 bfn_hess(ipt,dzz,hxx) 
     &                           = bfn_hess(ipt,dzz,hxx)
     &                             +        dgzz + ddgxxzz
                                 bfn_hess(ipt,dzz,hyy) 
     &                           = bfn_hess(ipt,dzz,hyy)
     &                             +        dgzz + ddgyyzz
                                 bfn_hess(ipt,dzz,hzz) 
     &                           = bfn_hess(ipt,dzz,hzz)
     &                             +g+g+5.d0*dgzz + ddgzzzz
                                 bfn_hess(ipt,dzz,hxy) 
     &                           = bfn_hess(ipt,dzz,hxy)
     &                             +               ddgxyzz
                                 bfn_hess(ipt,dzz,hyz) 
     &                           = bfn_hess(ipt,dzz,hyz)
     &                             +   dgyz+dgyz + ddgyzzz
                                 bfn_hess(ipt,dzz,hxz) 
     &                           = bfn_hess(ipt,dzz,hxz)
     &                             +   dgxz+dgxz + ddgxzzz

                                 bfn_hess(ipt,dxy,hxx) 
     &                           = bfn_hess(ipt,dxy,hxx)
     &                             +(  3.d0*dgxy  + ddgxxxy)*sqrt3
                                 bfn_hess(ipt,dxy,hyy) 
     &                           = bfn_hess(ipt,dxy,hyy)
     &                             +(  3.d0*dgxy  + ddgxyyy)*sqrt3
                                 bfn_hess(ipt,dxy,hzz) 
     &                           = bfn_hess(ipt,dxy,hzz)
     &                             +(       dgxy  + ddgxyzz)*sqrt3
                                 bfn_hess(ipt,dxy,hxy) 
     &                           = bfn_hess(ipt,dxy,hxy)
     &                             +(g+dgxx+dgyy  + ddgxxyy)*sqrt3
                                 bfn_hess(ipt,dxy,hyz) 
     &                           = bfn_hess(ipt,dxy,hyz)
     &                             +(       dgxz  + ddgxyyz)*sqrt3    
                                 bfn_hess(ipt,dxy,hxz) 
     &                           = bfn_hess(ipt,dxy,hxz)
     &                             +(       dgyz  + ddgxxyz)*sqrt3

                                 bfn_hess(ipt,dxz,hxx) 
     &                           = bfn_hess(ipt,dxz,hxx)
     &                             +( 3.d0*dgxz + ddgxxxz)*sqrt3
                                 bfn_hess(ipt,dxz,hyy) 
     &                           = bfn_hess(ipt,dxz,hyy)
     &                             +(  dgxz     + ddgxyyz)*sqrt3
                                 bfn_hess(ipt,dxz,hzz) 
     &                           = bfn_hess(ipt,dxz,hzz)
     &                             +( 3.d0*dgxz + ddgxzzz)*sqrt3
                                 bfn_hess(ipt,dxz,hxy) 
     &                           = bfn_hess(ipt,dxz,hxy)
     &                             +(  dgyz     + ddgxxyz)*sqrt3
                                 bfn_hess(ipt,dxz,hyz) 
     &                           = bfn_hess(ipt,dxz,hyz)
     &                             +(  dgxy     + ddgxyzz)*sqrt3
                                 bfn_hess(ipt,dxz,hxz) 
     &                           = bfn_hess(ipt,dxz,hxz)
     &                             +(g+dgxx+dgzz+ ddgxxzz)*sqrt3
C     
                                 bfn_hess(ipt,dyz,hxx) 
     &                           = bfn_hess(ipt,dyz,hxx)
     &                             +(  dgyz     + ddgxxyz)*sqrt3     
                                 bfn_hess(ipt,dyz,hyy) 
     &                           = bfn_hess(ipt,dyz,hyy)
     &                             +( 3.d0*dgyz + ddgyyyz)*sqrt3
                                 bfn_hess(ipt,dyz,hzz) 
     &                           = bfn_hess(ipt,dyz,hzz)
     &                             +( 3.d0*dgyz + ddgyzzz)*sqrt3
                                 bfn_hess(ipt,dyz,hxy) 
     &                           = bfn_hess(ipt,dyz,hxy)
     &                             +(  dgxz     + ddgxyyz)*sqrt3
                                 bfn_hess(ipt,dyz,hyz) 
     &                           = bfn_hess(ipt,dyz,hyz)
     &                             +(g+dgyy+dgzz+ ddgyyzz)*sqrt3
                                 bfn_hess(ipt,dyz,hxz) 
     &                           = bfn_hess(ipt,dyz,hxz)
     &                             +(  dgxy     + ddgxyzz)*sqrt3

                              endif

                              if (ider.ge.3) then
                                 g3   = g*sqrt3
                                 alp2 = -2.0d0*alp
                                 alp4 = alp2*alp2
                                 alp6 = -6.0d0*alp
                                 alp8 = alp4*alp2
                                 alp10 = -10.0d0*alp
                                 alp12 = 3.0d0*alp4
                                 alp24 = 6.0d0*alp4
                                 alp36 = alp12+alp24
                                 x0d1 = alp2*xd(ipt)
                                 x0d2 = alp2+alp4*xx
                                 x0d3 = (alp12+alp8*xx)*xd(ipt)
                                 x1d0 = xd(ipt)
                                 x1d1 = 1.0d0+alp2*xx
                                 x1d2 = (alp6+alp4*xx)*xd(ipt)
                                 x1d3 = alp6+(alp24+alp8*xx)*xx
                                 x2d0 = xx
                                 x2d1 = (2.0d0+alp2*xx)*xd(ipt)
                                 x2d2 = 2.0d0+xx*(alp10+alp4*xx)
                                 x2d3 = xd(ipt)
     &                                * (-24*alp+xx*(alp36+alp8*xx))
                                 y0d1 = alp2*yd(ipt)
                                 y0d2 = alp2+alp4*yy
                                 y0d3 = (alp12+alp8*yy)*xd(ipt)
                                 y1d0 = yd(ipt)
                                 y1d1 = 1.0d0+alp2*yy
                                 y1d2 = (alp6+alp4*yy)*yd(ipt)
                                 y1d3 = alp6+(alp24+alp8*yy)*yy
                                 y2d0 = yy
                                 y2d1 = (2.0d0+alp2*yy)*yd(ipt)
                                 y2d2 = 2.0d0+yy*(alp10+alp4*yy)
                                 y2d3 = yd(ipt)
     &                                * (-24*alp+yy*(alp36+alp8*yy))
                                 z0d1 = alp2*zd(ipt)
                                 z0d2 = alp2+alp4*zz
                                 z0d3 = (alp12+alp8*zz)*zd(ipt)
                                 z1d0 = zd(ipt)
                                 z1d1 = 1.0d0+alp2*zz
                                 z1d2 = (alp6+alp4*zz)*zd(ipt)
                                 z1d3 = alp6+(alp24+alp8*zz)*zz
                                 z2d0 = zz
                                 z2d1 = (2.0d0+alp2*zz)*zd(ipt)
                                 z2d2 = 2.0d0+zz*(alp10+alp4*zz)
                                 z2d3 = zd(ipt)
     &                                * (-24*alp+zz*(alp36+alp8*zz))
c
                                 bfn_3rd(ipt,dxx,txxx)
     &                           = bfn_3rd(ipt,dxx,txxx)
     &                             + x2d3*g
                                 bfn_3rd(ipt,dxx,txxy)
     &                           = bfn_3rd(ipt,dxx,txxy)
     &                             + x2d2*y0d1*g
                                 bfn_3rd(ipt,dxx,txxz)
     &                           = bfn_3rd(ipt,dxx,txxz)
     &                             + x2d2*z0d1*g
                                 bfn_3rd(ipt,dxx,txyy)
     &                           = bfn_3rd(ipt,dxx,txyy)
     &                             + x2d1*y0d2*g
                                 bfn_3rd(ipt,dxx,txyz)
     &                           = bfn_3rd(ipt,dxx,txyz)
     &                             + x2d1*y0d1*z0d1*g
                                 bfn_3rd(ipt,dxx,txzz)
     &                           = bfn_3rd(ipt,dxx,txzz)
     &                             + x2d1*z0d2*g
                                 bfn_3rd(ipt,dxx,tyyy)
     &                           = bfn_3rd(ipt,dxx,tyyy)
     &                             + x2d0*y0d3*g
                                 bfn_3rd(ipt,dxx,tyyz)
     &                           = bfn_3rd(ipt,dxx,tyyz)
     &                             + x2d0*y0d2*z0d1*g
                                 bfn_3rd(ipt,dxx,tyzz)
     &                           = bfn_3rd(ipt,dxx,tyzz)
     &                             + x2d0*y0d1*z0d2*g
                                 bfn_3rd(ipt,dxx,tzzz)
     &                           = bfn_3rd(ipt,dxx,tzzz)
     &                             + x2d0*z0d3*g
c
                                 bfn_3rd(ipt,dyy,txxx)
     &                           = bfn_3rd(ipt,dyy,txxx)
     &                             + x0d3*y2d0*g
                                 bfn_3rd(ipt,dyy,txxy)
     &                           = bfn_3rd(ipt,dyy,txxy)
     &                             + x0d2*y2d1*g
                                 bfn_3rd(ipt,dyy,txxz)
     &                           = bfn_3rd(ipt,dyy,txxz)
     &                             + x0d2*y2d0*z0d1*g
                                 bfn_3rd(ipt,dyy,txyy)
     &                           = bfn_3rd(ipt,dyy,txyy)
     &                             + x0d1*y2d2*g
                                 bfn_3rd(ipt,dyy,txyz)
     &                           = bfn_3rd(ipt,dyy,txyz)
     &                             + x0d1*y2d1*z0d1*g
                                 bfn_3rd(ipt,dyy,txzz)
     &                           = bfn_3rd(ipt,dyy,txzz)
     &                             + x0d1*y2d0*z0d2*g
                                 bfn_3rd(ipt,dyy,tyyy)
     &                           = bfn_3rd(ipt,dyy,tyyy)
     &                             + y2d3*g
                                 bfn_3rd(ipt,dyy,tyyz)
     &                           = bfn_3rd(ipt,dyy,tyyz)
     &                             + y2d2*z0d1*g
                                 bfn_3rd(ipt,dyy,tyzz)
     &                           = bfn_3rd(ipt,dyy,tyzz)
     &                             + y2d1*z0d2*g
                                 bfn_3rd(ipt,dyy,tzzz)
     &                           = bfn_3rd(ipt,dyy,tzzz)
     &                             + y2d0*z0d3*g
c
                                 bfn_3rd(ipt,dzz,txxx)
     &                           = bfn_3rd(ipt,dzz,txxx)
     &                             + x0d3*z2d0*g
                                 bfn_3rd(ipt,dzz,txxy)
     &                           = bfn_3rd(ipt,dzz,txxy)
     &                             + x0d2*y0d1*z2d0*g
                                 bfn_3rd(ipt,dzz,txxz)
     &                           = bfn_3rd(ipt,dzz,txxz)
     &                             + x0d2*z2d1*g
                                 bfn_3rd(ipt,dzz,txyy)
     &                           = bfn_3rd(ipt,dzz,txyy)
     &                             + x0d1*y0d2*z2d0*g
                                 bfn_3rd(ipt,dzz,txyz)
     &                           = bfn_3rd(ipt,dzz,txyz)
     &                             + x0d1*y0d1*z2d1*g
                                 bfn_3rd(ipt,dzz,txzz)
     &                           = bfn_3rd(ipt,dzz,txzz)
     &                             + x0d1*z2d2*g
                                 bfn_3rd(ipt,dzz,tyyy)
     &                           = bfn_3rd(ipt,dzz,tyyy)
     &                             + y0d3*z2d0*g
                                 bfn_3rd(ipt,dzz,tyyz)
     &                           = bfn_3rd(ipt,dzz,tyyz)
     &                             + y0d2*z2d1*g
                                 bfn_3rd(ipt,dzz,tyzz)
     &                           = bfn_3rd(ipt,dzz,tyzz)
     &                             + y0d1*z2d2*g
                                 bfn_3rd(ipt,dzz,tzzz)
     &                           = bfn_3rd(ipt,dzz,tzzz)
     &                             + z2d3*g
c
                                 bfn_3rd(ipt,dxy,txxx)
     &                           = bfn_3rd(ipt,dxy,txxx)
     &                             + x1d3*y1d0*g3
                                 bfn_3rd(ipt,dxy,txxy)
     &                           = bfn_3rd(ipt,dxy,txxy)
     &                             + x1d2*y1d1*g3
                                 bfn_3rd(ipt,dxy,txxz)
     &                           = bfn_3rd(ipt,dxy,txxz)
     &                             + x1d2*y1d0*z0d1*g3
                                 bfn_3rd(ipt,dxy,txyy)
     &                           = bfn_3rd(ipt,dxy,txyy)
     &                             + x1d1*y1d2*g3
                                 bfn_3rd(ipt,dxy,txyz)
     &                           = bfn_3rd(ipt,dxy,txyz)
     &                             + x1d1*y1d1*z0d1*g3
                                 bfn_3rd(ipt,dxy,txzz)
     &                           = bfn_3rd(ipt,dxy,txzz)
     &                             + x1d1*y1d0*z0d2*g3
                                 bfn_3rd(ipt,dxy,tyyy)
     &                           = bfn_3rd(ipt,dxy,tyyy)
     &                             + x1d0*y1d3*g3
                                 bfn_3rd(ipt,dxy,tyyz)
     &                           = bfn_3rd(ipt,dxy,tyyz)
     &                             + x1d0*y1d2*z0d1*g3
                                 bfn_3rd(ipt,dxy,tyzz)
     &                           = bfn_3rd(ipt,dxy,tyzz)
     &                             + x1d0*y1d1*z0d2*g3
                                 bfn_3rd(ipt,dxy,tzzz)
     &                           = bfn_3rd(ipt,dxy,tzzz)
     &                             + x1d0*y1d0*z0d3*g3
c
                                 bfn_3rd(ipt,dxz,txxx)
     &                           = bfn_3rd(ipt,dxz,txxx)
     &                             + x1d3*z1d0*g3
                                 bfn_3rd(ipt,dxz,txxy)
     &                           = bfn_3rd(ipt,dxz,txxy)
     &                             + x1d2*y0d1*z1d0*g3
                                 bfn_3rd(ipt,dxz,txxz)
     &                           = bfn_3rd(ipt,dxz,txxz)
     &                             + x1d2*z1d1*g3
                                 bfn_3rd(ipt,dxz,txyy)
     &                           = bfn_3rd(ipt,dxz,txyy)
     &                             + x1d1*y0d2*z1d0*g3
                                 bfn_3rd(ipt,dxz,txyz)
     &                           = bfn_3rd(ipt,dxz,txyz)
     &                             + x1d1*y0d1*z1d1*g3
                                 bfn_3rd(ipt,dxz,txzz)
     &                           = bfn_3rd(ipt,dxz,txzz)
     &                             + x1d1*z1d2*g3
                                 bfn_3rd(ipt,dxz,tyyy)
     &                           = bfn_3rd(ipt,dxz,tyyy)
     &                             + x1d0*y0d3*z1d0*g3
                                 bfn_3rd(ipt,dxz,tyyz)
     &                           = bfn_3rd(ipt,dxz,tyyz)
     &                             + x1d0*y0d2*z1d1*g3
                                 bfn_3rd(ipt,dxz,tyzz)
     &                           = bfn_3rd(ipt,dxz,tyzz)
     &                             + x1d0*y0d1*z1d2*g3
                                 bfn_3rd(ipt,dxz,tzzz)
     &                           = bfn_3rd(ipt,dxz,tzzz)
     &                             + x1d0*z1d3*g3
c
                                 bfn_3rd(ipt,dyz,txxx)
     &                           = bfn_3rd(ipt,dyz,txxx)
     &                             + x0d3*y1d0*z1d0*g3
                                 bfn_3rd(ipt,dyz,txxy)
     &                           = bfn_3rd(ipt,dyz,txxy)
     &                             + x0d2*y1d1*z1d0*g3
                                 bfn_3rd(ipt,dyz,txxz)
     &                           = bfn_3rd(ipt,dyz,txxz)
     &                             + x0d2*y1d0*z1d1*g3
                                 bfn_3rd(ipt,dyz,txyy)
     &                           = bfn_3rd(ipt,dyz,txyy)
     &                             + x0d1*y1d2*z1d0*g3
                                 bfn_3rd(ipt,dyz,txyz)
     &                           = bfn_3rd(ipt,dyz,txyz)
     &                             + x0d1*y1d1*z1d1*g3
                                 bfn_3rd(ipt,dyz,txzz)
     &                           = bfn_3rd(ipt,dyz,txzz)
     &                             + x0d1*y1d0*z1d2*g3
                                 bfn_3rd(ipt,dyz,tyyy)
     &                           = bfn_3rd(ipt,dyz,tyyy)
     &                             + y1d3*z1d0*g3
                                 bfn_3rd(ipt,dyz,tyyz)
     &                           = bfn_3rd(ipt,dyz,tyyz)
     &                             + y1d2*z1d1*g3
                                 bfn_3rd(ipt,dyz,tyzz)
     &                           = bfn_3rd(ipt,dyz,tyzz)
     &                             + y1d1*z1d2*g3
                                 bfn_3rd(ipt,dyz,tzzz)
     &                           = bfn_3rd(ipt,dyz,tzzz)
     &                             + y1d0*z1d3*g3
                              endif
                              bpos=6

                           elseif(lhyb .eq. 4 .or.
     1                            lhyb .eq. 5)then
c
c                             F and G functions
c
                              x=xd(ipt)
                              y=yd(ipt)
                              z=zd(ipt)
                              
                              if (x.eq.0d0) x=small
                              if (y.eq.0d0) y=small
                              if (z.eq.0d0) z=small

                              if (lhyb.eq.4) ilpos = 0
                              if (lhyb.eq.5) ilpos = 10

                              bpos = nftab(lhyb,lhyb)

                              do ilm  = 1,bpos
                                 
                                 k = kq(ilpos+ilm)
                                 l = lq(ilpos+ilm)
                                 m = mq(ilpos+ilm)

                                 e = g * ffac(ilpos+ilm)

                                 bfn_val(ipt,iloc+ilm-1) = 
     &                                bfn_val(ipt,iloc+ilm-1) + e * x**k
     $                                * y**l * z**m 
                                 
                                 if(ider.ge.1)then

                                    bfng_val(ipt,iloc+ilm-1,1) =
     $                                   bfng_val(ipt,iloc+ilm-1,1) +e
     $                                   * y**l * z**m *(k*x**(k-1) - 2
     $                                   *alp*x**(k+1)) 

                                    bfng_val(ipt,iloc+ilm-1,2) =
     $                                   bfng_val(ipt,iloc+ilm-1,2) +e
     $                                   * x**k * z**m *(l*y**(l-1) - 2
     $                                   *alp*y**(l+1))

                                    bfng_val(ipt,iloc+ilm-1,3) =
     $                                   bfng_val(ipt,iloc+ilm-1,3) +e
     $                                   * x**k * y**l *(m*z**(m-1) - 2
     $                                   *alp*z**(m+1))

                                 endif

                                 if(ider.ge.2)then

                                    bfn_hess(ipt,iloc+ilm-1,hxx) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hxx) +e
     $                                   * y**l * z**m * (k*(k-1)*x**(k
     $                                   -2)-2*alp*(2*k+1)*x**k+4*alp**2
     $                                   *x**(k+2))

                                    bfn_hess(ipt,iloc+ilm-1,hyy) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hyy) +e
     $                                   * x**k * z**m * (l*(l-1)*y**(l
     $                                   -2)-2*alp*(2*l+1)*y**l+4*alp**2
     $                                   *y**(l+2))

                                    bfn_hess(ipt,iloc+ilm-1,hzz) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hzz) +e
     $                                   * x**k * y**l * (m*(m-1)*z**(m
     $                                   -2)-2*alp*(2*m+1)*z**m+4*alp**2
     $                                   *z**(m+2))

                                    bfn_hess(ipt,iloc+ilm-1,hxy) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hxy) +e
     $                                   * z**m * (k*x**(k-1) - 2*alp*x
     $                                   **(k+1)) * (l*y**(l-1) - 2*alp
     $                                   *y**(l+1)) 

                                    bfn_hess(ipt,iloc+ilm-1,hxz) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hxz) +e
     $                                   * y**l * (k*x**(k-1) - 2*alp*x
     $                                   **(k+1)) * (m*z**(m-1) - 2*alp
     $                                   *z**(m+1)) 

                                    bfn_hess(ipt,iloc+ilm-1,hyz) =
     $                                   bfn_hess(ipt,iloc+ilm-1,hyz) +e
     $                                   * x**k * (l*y**(l-1) - 2*alp*y
     $                                   **(l+1)) * (m*z**(m-1) - 2*alp
     $                                   *z**(m+1))
                                    
                                 endif

                                 if (ider.ge.3) then
                                    bfn_3rd(ipt,iloc+ilm-1,txxx) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txxx) + e
     &                                 * y**l * z**m * (
     &                                 k*(k-1)*(k-2)*x**(k-3)
     &                                 -6*alp*k**2*x**(k-1)
     &                                 +12*alp**2*(k+1)*x**(k+1)
     &                                 -8*alp**3*x**(k+3))
                                    bfn_3rd(ipt,iloc+ilm-1,tyyy) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,tyyy) + e
     &                                 * x**k * z**m * (
     &                                 l*(l-1)*(l-2)*y**(l-3)
     &                                 -6*alp*l**2*y**(l-1)
     &                                 +12*alp**2*(l+1)*y**(l+1)
     &                                 -8*alp**3*y**(l+3))
                                    bfn_3rd(ipt,iloc+ilm-1,tzzz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,tzzz) + e
     &                                 * x**k * y**l * (
     &                                 m*(m-1)*(m-2)*z**(m-3)
     &                                 -6*alp*m**2*z**(m-1)
     &                                 +12*alp**2*(m+1)*z**(m+1)
     &                                 -8*alp**3*z**(m+3))
                                    bfn_3rd(ipt,iloc+ilm-1,txxy) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txxy) + e
     &                                 * z**m *
     &                                 (k*(k-1)*x**(k-2)
     &                                  -2*alp*(2*k+1)*x**k
     &                                  +4*alp**2*x**(k+2))*
     &                                 (l*y**(l-1)-2*alp*y**(l+1))
                                    bfn_3rd(ipt,iloc+ilm-1,txxz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txxz) + e
     &                                 * y**l *
     &                                 (k*(k-1)*x**(k-2)
     &                                  -2*alp*(2*k+1)*x**k
     &                                  +4*alp**2*x**(k+2))*
     &                                 (m*z**(m-1)-2*alp*z**(m+1))
                                    bfn_3rd(ipt,iloc+ilm-1,txyy) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txyy) + e
     &                                 * z**m *
     &                                 (l*(l-1)*y**(l-2)
     &                                  -2*alp*(2*l+1)*y**l
     &                                  +4*alp**2*y**(l+2))*
     &                                 (k*x**(k-1)-2*alp*x**(k+1))
                                    bfn_3rd(ipt,iloc+ilm-1,tyyz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,tyyz) + e
     &                                 * x**k *
     &                                 (l*(l-1)*y**(l-2)
     &                                  -2*alp*(2*l+1)*y**l
     &                                  +4*alp**2*y**(l+2))*
     &                                 (m*z**(m-1)-2*alp*z**(m+1))
                                    bfn_3rd(ipt,iloc+ilm-1,txzz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txzz) + e
     &                                 * y**l *
     &                                 (m*(m-1)*z**(m-2)
     &                                  -2*alp*(2*m+1)*z**m
     &                                  +4*alp**2*z**(m+2))*
     &                                 (k*x**(k-1)-2*alp*x**(k+1))
                                    bfn_3rd(ipt,iloc+ilm-1,tyzz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,tyzz) + e
     &                                 * x**k *
     &                                 (m*(m-1)*z**(m-2)
     &                                  -2*alp*(2*m+1)*z**m
     &                                  +4*alp**2*z**(m+2))*
     &                                 (l*y**(l-1)-2*alp*y**(l+1))
                                    bfn_3rd(ipt,iloc+ilm-1,txyz) = 
     &                                 bfn_3rd(ipt,iloc+ilm-1,txyz) + e*
     &                                 (k*x**(k-1)-2*alp*x**(k+1))*
     &                                 (l*y**(l-1)-2*alp*y**(l+1))*
     &                                 (m*z**(m-1)-2*alp*z**(m+1))
                                 endif
                                 
                              enddo
                              
                           else
                              call caserr('unimplemented shell type')
                           endif ! shell type
c
c                          update destination for next component
c
                           iloc=iloc + bpos 
                        enddo   ! lhyb
 200                 continue   ! lprm
                  else          ! while shell deleted for this point
                     del_psi = del_psi + float(bpos) / float(npts)
                     ploc = ploc + nprm
                  endif
               enddo            ! ipt
c     
c              store addresses of bfns we have just computed
c     
               icount=nftab(lh,ll)
               if (icount+istart-1.gt.max_active_bfn) then
                  write(iout,*)'n_active_bfn  =',icount+istart-1
                  write(iout,*)'max_active_bfn=',max_active_bfn
                  write(iout,*)'Change max_active_bfn in ',
     +                         'common/dft_parameters and recompile'
                  call caserr(
     +              "No. active basis functions exceeds code parameter")
               endif
               do i=1,icount
                  active_bfn_list(istart+i-1) = bfn+i-1
c                 active_bfn_indx(istart+i-1) = lcent
                  active_bfn_indx(istart+i-1) = icent
               enddo
               istart = istart + icount
            endif               ! test on shell
            bfn=bfn+nftab(lh,ll) ! counter for index in complete list
         enddo                  ! lshl
 100     continue
      enddo                     ! lcent
      n_active_bfn=istart-1
cDEBUG
c     if (ider.ge.1) then
c        do i=1,3
c           do iloc=1,n_active_bfn
c              do ipt=1,npts
c                 if (dabs(bfng_val(ipt,iloc,i)).lt.1.0d-2*psitol)
c    &                bfng_val(ipt,iloc,i)=0.0d0
c              enddo
c           enddo
c        enddo
c     endif
c     if (ider.ge.2) then
c        do i=1,6
c           do iloc=1,n_active_bfn
c              do ipt=1,npts
c                 if (dabs(bfn_hess(ipt,iloc,i)).lt.1.0d-2*psitol)
c    &                bfn_hess(ipt,iloc,i)=0.0d0
c              enddo
c           enddo
c        enddo
c     endif
c     if (ider.ge.3) then
c        do i=1,10
c           do iloc=1,n_active_bfn
c              do ipt=1,npts
c                 if (dabs(bfn_3rd(ipt,iloc,i)).lt.1.0d-2*psitol)
c    &                bfn_3rd(ipt,iloc,i)=0.0d0
c              enddo
c           enddo
c        enddo
c     endif
c     write(*,*)'*** bas_val ',bfng_val(155,9,3)
c     write(*,*)'*** bas_val2',loc(bfng_val(155,9,3))
cDEBUG
      end
c
c-----------------------------------------------------------------------
c
      subroutine calc_mo_val(rks_sw,ider,npts,mxp,nao,
     &     nvec,navec,nbvec,
     &     avec,bvec,
     &     bfn_val,bfng_val,bfn_hess,
     &     amo_val,amo_grad,amo_hess,
     &     bmo_val,bmo_grad,bmo_hess)
      implicit none
c
c     Constructs the MO values, gradients and hessians from the AO 
c     quantities.
c
c     Inputs:
c
      logical rks_sw 
      integer npts    ! the number of grid points
      integer mxp     ! the maximum number of grid points
      integer nao     ! the number of AOs
      integer nvec    ! the number of vectors
      integer navec   ! the number of alpha-vectors
      integer nbvec   ! the number of beta-vectors
      REAL avec(nao,navec) ! the alpha-vectors
      REAL bvec(nao,nbvec) ! the beta-vectors
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      integer ider
c
c     Outputs:
c
      REAL amo_val(npts,navec)
      REAL amo_grad(npts,navec,3)
      REAL amo_hess(npts,navec,6)
      REAL bmo_val(npts,nbvec)
      REAL bmo_grad(npts,nbvec,3)
      REAL bmo_hess(npts,nbvec,6)
c
c     Local variables
c
      integer ivec
      integer iao
      integer ipt
      integer i
c
c     Code:
c
      do ivec=1,navec
         do ipt=1,npts
            amo_val(ipt,ivec) = 0.0d0
         enddo
         do iao=1,nao
            do ipt=1,npts
               amo_val(ipt,ivec)=amo_val(ipt,ivec)
     &                          +bfn_val(ipt,iao)*avec(iao,ivec)
            enddo
         enddo
      enddo
      if (ider.ge.1) then
         do i=1,3
            do ivec=1,navec
               do ipt=1,npts
                  amo_grad(ipt,ivec,i) = 0.0d0
               enddo
               do iao=1,nao
                  do ipt=1,npts
                     amo_grad(ipt,ivec,i)=amo_grad(ipt,ivec,i)
     &               +bfng_val(ipt,iao,i)*avec(iao,ivec)
                  enddo
               enddo
            enddo
         enddo
      endif
      if (ider.ge.2) then
         do i=1,6
            do ivec=1,navec
               do ipt=1,npts
                  amo_hess(ipt,ivec,i) = 0.0d0
               enddo
               do iao=1,nao
                  do ipt=1,npts
                     amo_hess(ipt,ivec,i)=amo_hess(ipt,ivec,i)
     &               +bfn_hess(ipt,iao,i)*avec(iao,ivec)
                  enddo
               enddo
            enddo
         enddo
      endif
      if (.not.rks_sw) then
         do ivec=1,nbvec
            do ipt=1,npts
               bmo_val(ipt,ivec) = 0.0d0
            enddo
            do iao=1,nao
               do ipt=1,npts
                  bmo_val(ipt,ivec)=bmo_val(ipt,ivec)
     &                             +bfn_val(ipt,iao)*bvec(iao,ivec)
               enddo
            enddo
         enddo
         if (ider.ge.1) then
            do i=1,3
               do ivec=1,nbvec
                  do ipt=1,npts
                     bmo_grad(ipt,ivec,i) = 0.0d0
                  enddo
                  do iao=1,nao
                     do ipt=1,npts
                        bmo_grad(ipt,ivec,i)=bmo_grad(ipt,ivec,i)
     &                  +bfng_val(ipt,iao,i)*bvec(iao,ivec)
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (ider.ge.2) then
            do i=1,6
               do ivec=1,nbvec
                  do ipt=1,npts
                     bmo_hess(ipt,ivec,i) = 0.0d0
                  enddo
                  do iao=1,nao
                     do ipt=1,npts
                        bmo_hess(ipt,ivec,i)=bmo_hess(ipt,ivec,i)
     &                  +bfn_hess(ipt,iao,i)*bvec(iao,ivec)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine calc_mo_val_scr(rks_sw,ider,npts,mxp,nao,
     &     nvec,navec,nbvec,
     &     active_bfn_list,n_active_bfn,
     &     avec,bvec,
     &     bfn_val,bfng_val,bfn_hess,
     &     amo_val,amo_grad,amo_hess,
     &     bmo_val,bmo_grad,bmo_hess)
      implicit none
c
c     Constructs the MO values, gradients and hessians from the AO 
c     quantities.
c
c     Inputs:
c
      logical rks_sw 
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer npts    ! the number of grid points
      integer mxp     ! the maximum number of grid points
      integer nao     ! the number of AOs
      integer nvec    ! the number of vectors
      integer navec   ! the number of alpha-vectors
      integer nbvec   ! the number of beta-vectors
      REAL avec(nao,nvec) ! the alpha-vectors
      REAL bvec(nao,nvec) ! the beta-vectors
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      integer ider
c
c     Outputs:
c
      REAL amo_val(npts,navec)
      REAL amo_grad(npts,navec,3)
      REAL amo_hess(npts,navec,6)
      REAL bmo_val(npts,nbvec)
      REAL bmo_grad(npts,nbvec,3)
      REAL bmo_hess(npts,nbvec,6)
c
c     Local variables
c
      integer ivec
      integer iao
      integer ipt
      integer i
c
c     Code:
c
      do ivec=1,navec
         do ipt=1,npts
            amo_val(ipt,ivec) = 0.0d0
         enddo
         do iao=1,n_active_bfn
            do ipt=1,npts
               amo_val(ipt,ivec)=amo_val(ipt,ivec)
     &                          +bfn_val(ipt,iao)*
     &                           avec(active_bfn_list(iao),ivec)
            enddo
         enddo
      enddo
      if (ider.ge.1) then
         do i=1,3
            do ivec=1,navec
               do ipt=1,npts
                  amo_grad(ipt,ivec,i) = 0.0d0
               enddo
               do iao=1,n_active_bfn
                  do ipt=1,npts
cDEBUG
c                 write(*,*)'*** ',i,ivec,iao,ipt,active_bfn_list(iao)
c                 write(*,*)'*** ',amo_grad(ipt,ivec,i)
c                 write(*,*)'*** ',bfn_val(ipt,iao)
c                 write(*,*)'*** ',bfng_val(ipt,iao,i)
c                 write(*,*)'*** ',loc(bfng_val(ipt,iao,i))
c                 write(*,*)'*** ',avec(active_bfn_list(iao),ivec)
cDEBUG
                     amo_grad(ipt,ivec,i)=amo_grad(ipt,ivec,i)
     &               +bfng_val(ipt,iao,i)*
     &                avec(active_bfn_list(iao),ivec)
                  enddo
               enddo
            enddo
         enddo
      endif
      if (ider.ge.2) then
         do i=1,6
            do ivec=1,navec
               do ipt=1,npts
                  amo_hess(ipt,ivec,i) = 0.0d0
               enddo
               do iao=1,n_active_bfn
                  do ipt=1,npts
                     amo_hess(ipt,ivec,i)=amo_hess(ipt,ivec,i)
     &               +bfn_hess(ipt,iao,i)*
     &                avec(active_bfn_list(iao),ivec)
                  enddo
               enddo
            enddo
         enddo
      endif
      if (.not.rks_sw) then
         do ivec=1,nbvec
            do ipt=1,npts
               bmo_val(ipt,ivec) = 0.0d0
            enddo
            do iao=1,n_active_bfn
               do ipt=1,npts
                  bmo_val(ipt,ivec)=bmo_val(ipt,ivec)
     &                             +bfn_val(ipt,iao)*
     &                              bvec(active_bfn_list(iao),ivec)
               enddo
            enddo
         enddo
         if (ider.ge.1) then
            do i=1,3
               do ivec=1,nbvec
                  do ipt=1,npts
                     bmo_grad(ipt,ivec,i) = 0.0d0
                  enddo
                  do iao=1,n_active_bfn
                     do ipt=1,npts
                        bmo_grad(ipt,ivec,i)=bmo_grad(ipt,ivec,i)
     &                  +bfng_val(ipt,iao,i)*
     &                   bvec(active_bfn_list(iao),ivec)
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (ider.ge.2) then
            do i=1,6
               do ivec=1,nbvec
                  do ipt=1,npts
                     bmo_hess(ipt,ivec,i) = 0.0d0
                  enddo
                  do iao=1,n_active_bfn
                     do ipt=1,npts
                        bmo_hess(ipt,ivec,i)=bmo_hess(ipt,ivec,i)
     &                  +bfn_hess(ipt,iao,i)*
     &                   bvec(active_bfn_list(iao),ivec)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_val_mo(rks_sw,ider,npts,mxp,navec,nbvec,
     &                      naocc,nbocc,amo_val,amo_grad,bmo_val,
     &                      bmo_grad,rho,grho,rhotol,screen_sw)
      implicit none
c
c     Compute the density and the gradient of the density from the
c     MO values.
c
c     Inputs
c
      logical rks_sw
      integer ider
      integer npts
      integer mxp
      integer navec
      integer nbvec
      integer naocc
      integer nbocc
      REAL amo_val(npts,navec)
      REAL amo_grad(npts,navec,3)
      REAL bmo_val(npts,nbvec)
      REAL bmo_grad(npts,nbvec,3)
      REAL rhotol
      logical screen_sw
c
c     Outputs
c
      REAL rho(mxp,2)
      REAL grho(mxp,2,3)
c
c     Local variables
c
      integer ipt
      integer ivec
      integer i
c
c     Code
c
      call aclear_dp(rho,mxp*2,0.0d0)
      if (ider.ge.1) then
         call aclear_dp(grho,mxp*2*3,0.0d0)
      endif
c
      do ivec=1,naocc
         do ipt=1,npts
            rho(ipt,1)=rho(ipt,1)+amo_val(ipt,ivec)*amo_val(ipt,ivec)
         enddo
      enddo
      if (ider.ge.1) then
         do i=1,3
            do ivec=1,naocc
               do ipt=1,npts
                  grho(ipt,1,i)=grho(ipt,1,i)
     &            +amo_val(ipt,ivec)*amo_grad(ipt,ivec,i)
               enddo
            enddo
            do ipt=1,npts
               grho(ipt,1,i)=2*grho(ipt,1,i)
            enddo
         enddo
      endif
      if (rks_sw) then
         do ipt=1,npts
            rho(ipt,1)=2*rho(ipt,1)
         enddo
         if (ider.ge.1) then
            do i=1,3
               do ipt=1,npts
                  grho(ipt,1,i)=2*grho(ipt,1,i)
               enddo
            enddo
         endif
      else
         do ivec=1,nbocc
            do ipt=1,npts
               rho(ipt,2)=rho(ipt,2)+bmo_val(ipt,ivec)*bmo_val(ipt,ivec)
            enddo
         enddo
         if (ider.ge.1) then
            do i=1,3
               do ivec=1,nbocc
                  do ipt=1,npts
                     grho(ipt,2,i)=grho(ipt,2,i)
     &               +bmo_val(ipt,ivec)*bmo_grad(ipt,ivec,i)
                  enddo
               enddo
               do ipt=1,npts
                  grho(ipt,2,i)=2*grho(ipt,2,i)
               enddo
            enddo
         endif
      endif
c
      if (screen_sw) then
         do ipt=1,npts
            if (rho(ipt,1).lt.rhotol) rho(ipt,1)=0.0d0
         enddo
         if (.not.rks_sw) then
            do ipt=1,npts
               if (rho(ipt,2).lt.rhotol) rho(ipt,2)=0.0d0
            enddo
         endif
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_val_ao(rks_sw,hessian_sw,adens,bdens,
     &     tag,nao,
     &     bfn_val,bfng_val,bfn_hess,
     &     rho,grho,
     &     npts,mxp,screen_sw,
     &     active_bfn_list,n_active_bfn,ider)
C     ******************************************************************
C     *Description:						       *
C     *Calculate density and grad density at a point                   *
C     ******************************************************************
      implicit none

C     *Parameters	
INCLUDE(common/dft_parameters)

C In variables	
      integer npts, mxp
      logical rks_sw,hessian_sw
      REAL adens(*),bdens(*)
      REAL bfn_val(mxp,*)
      integer tag, nao
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      
      logical screen_sw
      integer active_bfn_list(*), n_active_bfn
      integer ider
C  Out variables
      REAL rho(mxp,2),grho(mxp,2,3)

C  Local variables
      integer ialpha, ibeta
      integer lbasi,i,ioff
      REAL t(500)
      REAL fac,eps_prim

      integer ipt, ixyz
      character*1 xn, xt
      data xt,xn/'t','n'/

      if(npts.gt.500) call caserr("npts exceeds size of t in den_val")

      ialpha = 1
      ibeta  = 2

      call aclear_dp(rho,2*mxp,0.0d0)
      if(ider.gt.0)call aclear_dp(grho,3*2*mxp,0.0d0)

      eps_prim=1.0d-15
      fac=2.0d0
C     
C  Calculate rho and grad rho
c
      ioff = 1
      do lbasi=1,nao
c     
c  Adjust off-diagonal density for matrix mult
c
         do i=0,lbasi-2
            adens(ioff+i) = adens(ioff+i) * 2.0d0
         enddo

         if(.not. rks_sw)then
            do i=0,lbasi-2
               bdens(ioff+i) = bdens(ioff+i) * 2.0d0
            enddo
         endif
c
c multiply a row of the triangulated density by the 
c basis function values, store in t
c
         call dgemv(xn,npts,lbasi,1.0d0,bfn_val,mxp
     &        ,adens(ioff),1,0.0d0,t,1)

         do ipt = 1,npts
            rho(ipt,ialpha)=rho(ipt,ialpha) + 
     &           bfn_val(ipt,lbasi)*t(ipt)
         enddo

         if(ider.gt.0)then
            do ipt = 1,npts
               grho(ipt,ialpha,1)=grho(ipt,ialpha,1) + 
     &              bfng_val(ipt,lbasi,1)*t(ipt)
               grho(ipt,ialpha,2)=grho(ipt,ialpha,2) + 
     &              bfng_val(ipt,lbasi,2)*t(ipt)
               grho(ipt,ialpha,3)=grho(ipt,ialpha,3) + 
     &              bfng_val(ipt,lbasi,3)*t(ipt)
            enddo
         endif

         if(.not. rks_sw)then
            call dgemv(xn,npts,lbasi,1.0d0,bfn_val,mxp
     &           ,bdens(ioff),1,0.0d0,t,1)

            do ipt = 1,npts
               rho(ipt,ibeta)=rho(ipt,ibeta) + 
     &              bfn_val(ipt,lbasi)*t(ipt)
            enddo

            if(ider.gt.0)then
               do ipt = 1,npts
                  grho(ipt,ibeta,1)=grho(ipt,ibeta,1) + 
     &                 bfng_val(ipt,lbasi,1)*t(ipt)
                  grho(ipt,ibeta,2)=grho(ipt,ibeta,2) + 
     &                 bfng_val(ipt,lbasi,2)*t(ipt)
                  grho(ipt,ibeta,3)=grho(ipt,ibeta,3) + 
     &                 bfng_val(ipt,lbasi,3)*t(ipt)
               enddo
               
            endif
         endif

         if(ider.gt.0)then
            do ixyz=1,3
               call dgemv(xn,npts,lbasi
     &              ,1.0d0,bfng_val(1,1,ixyz)
     &              ,mxp,adens(ioff),1,0.0d0,t,1)
               do ipt = 1,npts
                  grho(ipt,ialpha,ixyz)=grho(ipt,ialpha,ixyz) + 
     &                 bfn_val(ipt,lbasi)*t(ipt)
               enddo
            enddo
            if(.not. rks_sw)then
               do ixyz=1,3
                  call dgemv(xn,npts,lbasi
     &                 ,1.0d0,bfng_val(1,1,ixyz)
     &                 ,mxp,bdens(ioff),1,0.0d0,t,1)
                  do ipt = 1,npts
                     grho(ipt,ibeta,ixyz)=grho(ipt,ibeta,ixyz) + 
     &                    bfn_val(ipt,lbasi)*t(ipt)
                  enddo
               enddo
            endif
         endif

         do i=0,lbasi-2
            adens(ioff+i) = adens(ioff+i) * 0.5d0
         enddo
         if(.not. rks_sw)then
            do i=0,lbasi-2
               bdens(ioff+i) = bdens(ioff+i) * 0.5d0
            enddo
         endif
         ioff = ioff + lbasi

      enddo
c
c...  Protection against round-off errors with small densities...
c
      do ipt = 1, npts
         rho(ipt,ialpha) = max(0.0d0,rho(ipt,ialpha))
         rho(ipt,ibeta)  = max(0.0d0,rho(ipt,ibeta))
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
c density and gradient for screened basis function list, also 
c incorporates test on density elements
c
      subroutine den_val_ao_scr(rks_sw,hessian_sw,adens,bdens,
     &     tag,nao,bfn_val,bfng_val,bfn_hess,
     &     rho,grho,npts,mxp,screen_sw,
     &     active_bfn_list,n_active_bfn,dentol,rhotol,ider)

C     ******************************************************************
C     *Description:						       *
C     *Calculate density and grad density at a point                   *
C     ******************************************************************
      implicit none
C     ******************************************************************
C     *Declarations				
C     *Parameters				
INCLUDE(common/dft_parameters)
C     *In variables				
      integer npts, mxp
      logical rks_sw,hessian_sw
      REAL adens(*),bdens(*)
      REAL bfn_val(mxp,*)
      integer tag, nao
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)

      logical screen_sw
      integer active_bfn_list(*), n_active_bfn
      REAL dentol
      REAL rhotol
      integer ider

C     *Out variables
      REAL rho(mxp,2),grho(mxp,2,3)
C     *Local variables
      integer ialpha,ibeta
      integer lbasi,lbasj
      integer iki, ibi
      REAL t(500)
      REAL fac,eps_prim
      REAL ptmp

      integer ipt
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)

C     *End declarations                                                *
C     ******************************************************************
      ialpha = 1
      ibeta  = 2

c     call aclear_dp(rho,2*mxp,0.0d0)
c     if(ider.gt.0)call aclear_dp(grho,3*2*mxp,0.0d0)
      call aclear_dp(rho(1,1),npts,0.0d0)
      if(ider.gt.0) then
         call aclear_dp(grho(1,1,1),npts,0.0d0)
         call aclear_dp(grho(1,1,2),npts,0.0d0)
         call aclear_dp(grho(1,1,3),npts,0.0d0)
      endif

      eps_prim=1.0d-15
      fac=2.0d0
C     
C     Calculate rho and grad rho

      do lbasi=1,n_active_bfn

         call aclear_dp(t,npts,0.0d0)
         ibi=active_bfn_list(lbasi)
         iki=iky(ibi)
         do lbasj=1,lbasi
            ptmp = 2.0d0*adens(iki + active_bfn_list(lbasj))
            if(abs(ptmp) .gt. dentol)
     +call daxpy(npts,ptmp,bfn_val(1,lbasj),1,t,1)
         enddo
         do lbasj=lbasi+1,n_active_bfn
            ptmp = 2.0d0*adens(iky(active_bfn_list(lbasj)) + ibi)
            if(abs(ptmp) .gt. dentol)
     +call daxpy(npts,ptmp,bfn_val(1,lbasj),1,t,1)
         enddo

         do ipt = 1,npts
            rho(ipt,ialpha)=rho(ipt,ialpha) + 
     &           bfn_val(ipt,lbasi)*t(ipt)
         enddo
         if(ider.gt.0)then
            do ipt = 1,npts
               grho(ipt,ialpha,1)=grho(ipt,ialpha,1) + 
     &              bfng_val(ipt,lbasi,1)*t(ipt)
            enddo
            do ipt = 1,npts
               grho(ipt,ialpha,2)=grho(ipt,ialpha,2) + 
     &              bfng_val(ipt,lbasi,2)*t(ipt)
            enddo
            do ipt = 1,npts
               grho(ipt,ialpha,3)=grho(ipt,ialpha,3) + 
     &              bfng_val(ipt,lbasi,3)*t(ipt)
            enddo
         endif
      enddo

      if( .not. rks_sw) then
         call aclear_dp(rho(1,2),npts,0.0d0)
         if(ider.gt.0) then
            call aclear_dp(grho(1,2,1),npts,0.0d0)
            call aclear_dp(grho(1,2,2),npts,0.0d0)
            call aclear_dp(grho(1,2,3),npts,0.0d0)
         endif
         do lbasi=1,n_active_bfn
            call aclear_dp(t,npts,0.0d0)
            ibi=active_bfn_list(lbasi)
            iki=iky(ibi)
            do lbasj=1,lbasi
               ptmp = 2.0d0*bdens(iki + active_bfn_list(lbasj))
               if(abs(ptmp) .gt. dentol)
     +      call daxpy(npts,ptmp,bfn_val(1,lbasj),1,t,1)
            enddo
            do lbasj=lbasi+1,n_active_bfn
               ptmp = 2.0d0*bdens(iky(active_bfn_list(lbasj)) + ibi)
               if(abs(ptmp) .gt. dentol)
     +      call daxpy(npts,ptmp,bfn_val(1,lbasj),1,t,1)
            enddo

            do ipt = 1,npts
               rho(ipt,ibeta)=rho(ipt,ibeta) + 
     &              bfn_val(ipt,lbasi)*t(ipt)
            enddo
            if(ider.gt.0)then
               do ipt = 1,npts
                  grho(ipt,ibeta,1)=grho(ipt,ibeta,1) + 
     &                 bfng_val(ipt,lbasi,1)*t(ipt)
               enddo
               do ipt = 1,npts
                  grho(ipt,ibeta,2)=grho(ipt,ibeta,2) + 
     &                 bfng_val(ipt,lbasi,2)*t(ipt)
               enddo
               do ipt = 1,npts
                  grho(ipt,ibeta,3)=grho(ipt,ibeta,3) + 
     &                 bfng_val(ipt,lbasi,3)*t(ipt)
               enddo
            endif
         enddo
      endif
         
      call dscal(npts,0.5d0,rho(1,1),1)
      do ipt = 1, npts
         if (rho(ipt,1).lt.rhotol) rho(ipt,1) = 0.0d0
      enddo
      if(.not.rks_sw) then
         call dscal(npts,0.5d0,rho(1,2),1)
         do ipt = 1, npts
            if (rho(ipt,2).lt.rhotol) rho(ipt,2) = 0.0d0
         enddo
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tau_val_mo(rks_sw,npts,mxp,navec,nbvec,
     &                      naocc,nbocc,amo_grad,
     &                      bmo_grad,tau,rho,screen_sw)
      implicit none
c
c     Compute the kinetic energy density from the MO values.
c
c     Inputs
c
      logical rks_sw
      integer npts
      integer mxp
      integer navec
      integer nbvec
      integer naocc
      integer nbocc
      REAL amo_grad(npts,navec,3)
      REAL bmo_grad(npts,nbvec,3)
      REAL rho(mxp,2)
      logical screen_sw
c
c     Outputs
c
      REAL tau(mxp,2)
c
c     Local variables
c
      integer ipt
      integer ivec
      integer i
c
c     Code
c
      call aclear_dp(tau,mxp*2,0.0d0)
c
      do ivec=1,naocc
         do ipt=1,npts
            do i=1,3
               tau(ipt,1)=tau(ipt,1)
     &         +amo_grad(ipt,ivec,i)*amo_grad(ipt,ivec,i)
            enddo
         enddo
      enddo
      if (rks_sw) then
         do ipt=1,npts
            tau(ipt,1)=2*tau(ipt,1)
         enddo
      else
         do ivec=1,nbocc
            do ipt=1,npts
               do i=1,3
                  tau(ipt,2)=tau(ipt,2)
     &            +bmo_grad(ipt,ivec,i)*bmo_grad(ipt,ivec,i)
               enddo
            enddo
         enddo
      endif
c
      if (screen_sw) then
         do ipt=1,npts
            if (rho(ipt,1).le.0.0d0) tau(ipt,1)=0.0d0
         enddo
         if (.not.rks_sw) then
            do ipt=1,npts
               if (rho(ipt,2).lt.0.0d0) tau(ipt,2)=0.0d0
            enddo
         endif
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine tau_val_ao(rks_sw,adens,bdens,
     &     tag,nao,
     &     bfng_val,
     &     rho,tau,
     &     npts,mxp,screen_sw,
     &     active_bfn_list,n_active_bfn)
C     ******************************************************************
C     *Description:						       *
C     *Calculate kinetic energy density at a point                     *
C     ******************************************************************
      implicit none

C     *Parameters	
INCLUDE(common/dft_parameters)

C In variables	
      integer npts, mxp
      logical rks_sw
      integer tag, nao
      integer n_active_bfn, active_bfn_list(n_active_bfn)
      REAL adens(nao*(nao+1)/2),bdens(nao*(nao+1)/2)
      REAL rho(mxp,2)
      REAL bfng_val(mxp,nao,3)
      
      logical screen_sw
C  Out variables
      REAL tau(mxp,2)

C  Local variables
      integer ialpha, ibeta
      integer lbasi,i,ioff
      REAL t(500)
      REAL fac,eps_prim

      integer ipt, ixyz
      character*1 xn, xt
      data xt,xn/'t','n'/

      if(npts.gt.500) call caserr("npts exceeds size of t in den_val")

      ialpha = 1
      ibeta  = 2

      call aclear_dp(tau,2*mxp,0.0d0)

      eps_prim=1.0d-15
      fac=2.0d0
C     
C  Calculate rho and grad rho
c
      ioff = 1
      do lbasi=1,nao
c     
c  Adjust off-diagonal density for matrix mult
c
         do i=0,lbasi-2
            adens(ioff+i) = adens(ioff+i) * 2.0d0
         enddo

         if(.not. rks_sw)then
            do i=0,lbasi-2
               bdens(ioff+i) = bdens(ioff+i) * 2.0d0
            enddo
         endif

         do ixyz=1,3
            call dgemv(xn,npts,lbasi
     &           ,1.0d0,bfng_val(1,1,ixyz)
     &           ,mxp,adens(ioff),1,0.0d0,t,1)
            do ipt = 1,npts
               tau(ipt,ialpha)=tau(ipt,ialpha) + 
     &            bfng_val(ipt,lbasi,ixyz)*t(ipt)
            enddo
         enddo
         if(.not. rks_sw)then
            do ixyz=1,3
               call dgemv(xn,npts,lbasi
     &              ,1.0d0,bfng_val(1,1,ixyz)
     &              ,mxp,bdens(ioff),1,0.0d0,t,1)
               do ipt = 1,npts
                  tau(ipt,ibeta)=tau(ipt,ibeta) + 
     &               bfng_val(ipt,lbasi,ixyz)*t(ipt)
               enddo
            enddo
         endif

         do i=0,lbasi-2
            adens(ioff+i) = adens(ioff+i) * 0.5d0
         enddo
         if(.not. rks_sw)then
            do i=0,lbasi-2
               bdens(ioff+i) = bdens(ioff+i) * 0.5d0
            enddo
         endif
         ioff = ioff + lbasi

      enddo
c
c...  Protection against round-off errors with small densities...
c
      do ipt = 1, npts
         if (rho(ipt,ialpha).le.0.0d0) tau(ipt,ialpha) = 0.0d0
         if (rho(ipt,ibeta) .le.0.0d0) tau(ipt,ibeta)  = 0.0d0
cDEBUG
c     write(*,'("tau =",2f16.6)')tau(ipt,1),tau(ipt,2)
cDEBUG
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
c density and gradient for screened basis function list, also 
c incorporates test on density elements
c
      subroutine tau_val_ao_scr(rks_sw,adens,bdens,
     &     tag,nao,bfng_val,
     &     rho,tau,npts,mxp,screen_sw,
     &     active_bfn_list,n_active_bfn,dentol)

C     ******************************************************************
C     *Description:						       *
C     *Calculate kinetic energy density at a point                     *
C     ******************************************************************
      implicit none
C     ******************************************************************
C     *Declarations				
C     *Parameters				
INCLUDE(common/dft_parameters)
C     *In variables				
      integer npts, mxp
      logical rks_sw
      integer tag, nao
      integer n_active_bfn, active_bfn_list(n_active_bfn)
      REAL adens(nao*(nao+1)/2),bdens(nao*(nao+1)/2)
      REAL bfng_val(mxp,nao,3)
      REAL rho(mxp,2)

      logical screen_sw
      REAL dentol

C     *Out variables
      REAL tau(mxp,2)
C     *Local variables
      integer ialpha,ibeta
      integer lbasi,lbasj
      integer iki, ibi
      REAL t(200,3)
      REAL fac,eps_prim
      REAL ptmp

      integer ipt,i
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)

C     *End declarations                                                *
C     ******************************************************************
      ialpha = 1
      ibeta  = 2
      if(npts.gt.200) 
     &   call caserr("npts exceeds size of t in tau_val_ao_scr")

      call aclear_dp(tau(1,1),npts,0.0d0)

      eps_prim=1.0d-15
      fac=2.0d0
C     
C     Calculate rho and grad rho

      do lbasi=1,n_active_bfn

         call aclear_dp(t,3*200,0.0d0)
         ibi=active_bfn_list(lbasi)
         iki=iky(ibi)
         do lbasj=1,lbasi
            ptmp = adens(iki + active_bfn_list(lbasj))
            if(abs(ptmp) .gt. 0.5d0*dentol) then
               do i=1,3
                  call daxpy(npts,ptmp,bfng_val(1,lbasj,i),1,t(1,i),1)
               enddo
            endif
         enddo
         do lbasj=lbasi+1,n_active_bfn
            ptmp = adens(iky(active_bfn_list(lbasj)) + ibi)
            if(abs(ptmp) .gt. 0.5d0*dentol) then
               do i=1,3
                  call daxpy(npts,ptmp,bfng_val(1,lbasj,i),1,t(1,i),1)
               enddo
            endif
         enddo

         do ipt = 1,npts
            do i=1,3
               tau(ipt,ialpha)=tau(ipt,ialpha) + 
     &              bfng_val(ipt,lbasi,i)*t(ipt,i)
            enddo
         enddo
      enddo

      if( .not. rks_sw) then
         call aclear_dp(tau(1,ibeta),npts,0.0d0)
         do lbasi=1,n_active_bfn
            call aclear_dp(t,3*200,0.0d0)
            ibi=active_bfn_list(lbasi)
            iki=iky(ibi)
            do lbasj=1,lbasi
               ptmp = bdens(iki + active_bfn_list(lbasj))
               if(abs(ptmp) .gt. 0.5d0*dentol) then
                  do i=1,3
                    call daxpy(npts,ptmp,bfng_val(1,lbasj,i),1,t(1,i),1)
                  enddo
               endif
            enddo
            do lbasj=lbasi+1,n_active_bfn
               ptmp = bdens(iky(active_bfn_list(lbasj)) + ibi)
               if(abs(ptmp) .gt. 0.5d0*dentol) then
                  do i=1,3
                    call daxpy(npts,ptmp,bfng_val(1,lbasj,i),1,t(1,i),1)
                  enddo
               endif
            enddo

            do ipt = 1,npts
               do i=1,3
                  tau(ipt,ibeta)=tau(ipt,ibeta) + 
     &                 bfng_val(ipt,lbasi,i)*t(ipt,i)
               enddo
            enddo
         enddo
      endif
         
      do ipt = 1, npts
         if (rho(ipt,1).le.0.0d0) tau(ipt,1) = 0.0d0
      enddo
      if(.not.rks_sw) then
         do ipt = 1, npts
            if (rho(ipt,2).le.0.0d0) tau(ipt,2) = 0.0d0
         enddo
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine calc_gamma(mxp,npts,rks_sw,grho,gamma)
      implicit none
c
c     Computes (grad rhoa,grad rhoa), (grad rhoa, grad rhob) and
c     (grad rhob, grad rhob).
c
c     Parameters:
c
INCLUDE(common/dft_dfder)
c
c     Inputs:
c
      integer mxp
      integer npts
      logical rks_sw
      REAL grho(mxp,2,3)
c
c     Outputs
c
      REAL gamma(mxp,3)
c
c     Local:
c
      integer i
c
      if (rks_sw) then
         do i = 1, npts
            gamma(i,igaa) = grho(i,1,1)*grho(i,1,1)+
     &                      grho(i,1,2)*grho(i,1,2)+
     &                      grho(i,1,3)*grho(i,1,3)
         enddo
      else
         do i = 1, npts
            gamma(i,igaa) = grho(i,1,1)*grho(i,1,1)+
     &                      grho(i,1,2)*grho(i,1,2)+
     &                      grho(i,1,3)*grho(i,1,3)
            gamma(i,igab) = grho(i,1,1)*grho(i,2,1)+
     &                      grho(i,1,2)*grho(i,2,2)+
     &                      grho(i,1,3)*grho(i,2,3)
            gamma(i,igbb) = grho(i,2,1)*grho(i,2,1)+
     &                      grho(i,2,2)*grho(i,2,2)+
     &                      grho(i,2,3)*grho(i,2,3)
         enddo
      endif
      end
c
c-----------------------------------------------------------------------
c
c  routine to tabulate the extent of each basis function
c  Also stores the largest radius for functions on each atom
c
      subroutine calc_bfn_radii(tab,tag,nao,natm,arad,first_bf,igrid,
     &                          iwr)
      implicit none
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_order_info)
INCLUDE(common/dft_xc)
c in
      integer tag,nao,natm,igrid,iwr
c out
      REAL tab(nao)
      REAL arad(natm)
      integer first_bf(*)
c local
      logical owarned
      data owarned/.false./
      save owarned

      REAL r, rtest, delr, rmin
      REAL gg, cc, expo, alp, s, p, d, f, g
      integer lcent, count, ploc, lprm, lhyb, lshl
      integer nprm, nshells, ll, lh, centre
      logical otest, opg_root
c
c...  First check for some abnormalities
c
c     do lcent=1, natm
c        centre = atom_tag(tag,lcent)
c        if (gtype_num(lcent).eq.0) then
c           if (centre.ne.0) then
c              call caserr('atom with basis but no grid assigned')
c           endif
c        endif
c        if (centre.eq.0) then
c           if (gtype_num(lcent).ne.0) then
c              gtype_num(lcent)=0
c              if (.not.owarned.and.opg_root()) then
c                 write(iwr,600)lcent
c                 write(iwr,610)lcent
c                 owarned = .true.
c              endif
c           endif
c        endif
c     enddo
c600  format('WARNING: centre ',i4,' has a grid but no basis assigned')
c610  format('WARNING: grid on centre ',i4,' has been turned off')
c
c smallest r to consider
      rmin = 0.5d0
c increment in r
      delr = 0.1d0

      r = rmin

      otest = .true.

      do lcent=1, natoms
         arad(lcent) = 0.0d0
      enddo

      do count =1, totbfn(tag)
         tab(count)=rmin
      enddo

      first_bf(natm+1) = totbfn(tag) + 1

      do while (otest) 

         if (r.gt.10000.0d0) call caserr('calc_bfn_radii: bad r')

         otest = .false.

         count = 1
         
         do lcent=1, natm


            first_bf(lcent) = count

            centre = BL_get_atom_type(tag,lcent)
            if(centre.gt.0)then

               nshells = num_shl(tag,centre)
               ploc   = 0
               do lshl=1,nshells

                  rtest = 0.0d0

                  nprm   = nprim(tag,centre,lshl)
                  lh     = hybrid(tag,centre,lshl)
                  ll     = angmom(tag,centre,lshl)

                  s = 0.0d0
                  p = 0.0d0
                  d = 0.0d0
                  f = 0.0d0
                  g = 0.0d0

                  do lprm=1,nprm

                     ploc    = ploc + 1
                     alp     = alpha(tag,centre,ploc)
                     expo    = exp(-alp*r*r)

                     do lhyb=lh,ll

                        cc      = cont_coeff(tag,centre,ploc,lhyb)
                        gg      = expo*cc

                        if(lhyb.eq.1) then
                           s = s  + gg
                        else if(lhyb.eq.2) then
                           p = p  + gg*r
                        else if(lhyb.eq.3) then
                           d = d  + gg*r*r
                        else if(lhyb.eq.4) then
                           f = f  + gg*r*r*r
                        else if(lhyb.eq.5) then
                           g = g  + gg*r*r*r*r
                        endif
                     enddo
                  enddo
                  
                  do lhyb=lh,ll
                     if(lhyb.eq.1) then
                        if(dabs(s) .gt. psitol(gtype_num(lcent),igrid))
     &                  then
                           tab(count)=r
                           rtest = r
                           otest=.true.
                        endif
                        count=count+1
                     else if(lhyb.eq.2) then
                        if(dabs(p) .gt. psitol(gtype_num(lcent),igrid))
     &                  then
                           call aclear_dp(tab(count),3,r)
                           rtest = r
                           otest=.true.
                        endif
                        count=count+3
                     else if(lhyb.eq.3) then
                        if(dabs(d) .gt. psitol(gtype_num(lcent),igrid))
     &                  then
                           call aclear_dp(tab(count),6,r)
                           rtest = r
                           otest=.true.
                        endif
                        count=count+6
                     else if(lhyb.eq.4) then
                        if(dabs(f) .gt. psitol(gtype_num(lcent),igrid))
     &                  then
                           call aclear_dp(tab(count),10,r)
                           rtest = r
                           otest=.true.
                        endif
                        count=count+10
                     else if(lhyb.eq.5) then
                        if(dabs(g) .gt. psitol(gtype_num(lcent),igrid))
     &                  then
                           call aclear_dp(tab(count),15,r)
                           rtest = r
                           otest=.true.
                        endif
                        count=count+15
                     endif
                  enddo

                  arad(lcent) = max(arad(lcent),rtest)
                  
               enddo

            endif
            
         enddo

         r = r + delr

         if( r .lt. 60.0d0)otest = .true.

      enddo

      if (.not.screen_sw) then
         do lcent=1,natoms
            arad(lcent)=1.0d6
         enddo
         do count=1,totbfn(tag)
            tab(count)=1.0d6
         enddo
      endif
      end
c
c-----------------------------------------------------------------------
c
c     Override the calculated atomic radii in arad with user defined
c     radii stored in screen_atom_radius, if any.
c
      subroutine override_atom_radii(ngridcentres,igrid,arad)
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
      integer ngridcentres, igrid
      REAL arad(ngridcentres)
      integer i
c
      do i = 1, ngridcentres
         if (screen_atom_radius(gtype_num(i),igrid).ne.
     &       dble(DFT_UNDEF)) then
            arad(i) = screen_atom_radius(gtype_num(i),igrid)
         endif
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine glegend(nthet,nphi,nang,apts,awpt)
C **********************************************************************
C *Description:				Version: 1.0                   *
C *Gauss-Legendre angular quadrature                                   *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *Parameters                                                          *
INCLUDE(common/dft_parameters) 
INCLUDE(common/dft_physical_constants)
C *In variables                                                        *
      integer nthet,nphi,nang
      logical extwr_sw
C *Out variables                                                       *
      REAL apts(3,*),awpt(*)
C *Local variables

      integer maxgleg, maxgleg2
      parameter (maxgleg=80)
      parameter (maxgleg2=(maxgleg+1)/2)

      REAL pfact(2,maxgleg),tarr1(maxgleg2),tarr2(maxgleg2)
      REAL cosphi(maxgleg), sinphi(maxgleg)
      REAL tol,thet,x,p,p1,dp,dp1,q,dq,d,dphi,pw,dxy,phi
      integer nhpoints,nhpointe,i,k,l
C *Function
      logical opg_root
C *End declarations                                                    *
C **********************************************************************
      tol=1.0d-14
C *
C *      
      if(nthet .gt. maxgleg)call caserr('ntheta too big in glegend')

      do 10 i=2,nthet
        pfact(1,i)=2.0d0-(1.0d0/i)
        pfact(2,i)=1.0d0-(1.0d0/i)
10    continue
C
c...  Beware nthet can be an ODD-number !!!
c
      nhpointe=nthet/2
      nhpoints=(nthet+1)/2
      do 20 k=1,nhpoints
         thet = pi*(4.0d0*k-1.0d0)/(4.0d0*nthet+2.0d0)
         x    = cos(thet+1.0d0/(8.0d0*(nthet**2)*tan(thet)))
 35      p    = x
         p1   = 1.0d0
         dp   = 1.0d0
         dp1  = 0.0d0
         do 30 i=2,nthet
            q   = pfact(1,i)*x*p-pfact(2,i)*p1
            dq  = pfact(1,i)*(x*dp+p)-pfact(2,i)*dp1
            p1  = p
            p   = q
            dp1 = dp
            dp  = dq
 30      continue
         d=p/dp
         x=x-d
         if (abs(d).gt.tol) goto 35
         tarr1(k)=x
         tarr2(k)=2.0d0/((1.0d0-x**2)*(dp**2))
c
c         write(6,*)tarr1(k),tarr2(k)
c
 20   continue
C
      dphi=pi/nphi
      do l = 1, (nphi+1)/2
         phi=(2*l-1)*dphi
         cosphi(l) = cos(phi)
         sinphi(l) = sin(phi)
      enddo
c
      do 40 i=1,nhpoints
         k=(i-1)*nphi
         pw=tarr2(i)*2.0d0*dphi
         dxy=sqrt(1.0d0-tarr1(i)**2)
         do 50 l=1,(nphi+1)/2
            awpt(k+l)=pw
            awpt(k+nphi/2+l)=pw
c           phi=(2*l-1)*dphi
            apts(1,k+l)=dxy*cosphi(l)
            apts(2,k+l)=dxy*sinphi(l)
            apts(3,k+l)=tarr1(i)
            apts(1,k+nphi+1-l)=apts(1,k+l)
            apts(2,k+nphi+1-l)=-apts(2,k+l)
            apts(3,k+nphi+1-l)=tarr1(i)
 50      continue
 40   continue
C
      k=nhpoints*nphi
      do 60 l=1,nhpointe*nphi
        apts(1,k+l) =  apts(1,l)
        apts(2,k+l) =  apts(2,l)
        apts(3,k+l) = -apts(3,l)
        awpt(k+l)   =  awpt(l)
60    continue
c
      nang=nthet*nphi
      return
      end


      subroutine build_radgrid(ngtyp,igrid,nradpt_num,
     +                         prpt,prwt,extwr_sw,iout)
c **********************************************************************
c *                                                                    *
c * This subroutine builds the radial grids for all atom types.        *
c *                                                                    *
c **********************************************************************
      implicit none
c
c     Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_xc)
c
c     In variables
c
      integer ngtyp,igrid
      integer nradpt_num(ngtyp)
      logical extwr_sw
      integer iout
c
c     Out variables
c
      REAL prpt(ngtyp,*)
      REAL prwt(ngtyp,*)
c
c     Local variables
c
      integer i, j, atom_num
      REAL alph, alph3
c
c     Functions
c
      integer atom_num_by_grid
      REAL lograd, srad
c
c  This is used for null typed atoms (e.g. bqs)
c
      do i = 1, ngtyp
         if (rad_grid_scheme(i,igrid).eq.RG_SG1) then
            call premac(ngtyp,i,nradpt_num(i),prpt,prwt) 
         else if (rad_grid_scheme(i,igrid).eq.RG_EML) then
            call premac(ngtyp,i,nradpt_num(i),prpt,prwt) 
         elseif (rad_grid_scheme(i,igrid).eq.RG_MK) then
            call prelog(ngtyp,i,nradpt_num(i),radm_num(i,igrid),
     +                prpt,prwt) 
         endif
         atom_num = atom_num_by_grid(i)
         alph = grid_scale(i,igrid)*grid_atom_radius(i,igrid)
         alph3 = alph**3
         do j = 1, nradpt_num(i)
            prpt(i,j) = prpt(i,j)*alph
            prwt(i,j) = prwt(i,j)*alph3
         enddo
      enddo
c
      end


      subroutine premac(natyp,in,ntrad,prpt,prwt)
C **********************************************************************
C *Description:				Version: 2.1                   *
C * Radial quadrature using Euler-Maclaurin formula (6) and (7) [1]    *
c * with m = 2. The use of qr = i/(n+1), i=1,n requires a weight       *
C * adjustment for the 2 outer integration points.                     *
c *                                                                    *
c * No Points are removed, i.e. no screening.                          *
c *                                                                    *
c * [1] C.W. Murray, N.C. Handy, G.J. Laming                           *
c *     Quadrature schemes for integrals of density functional theory  *
c *     Molecular Physics, Vol. 78 (1993) pages 997-1014               *
C **********************************************************************
      implicit none
C *********************************************************************
C *Declarations                                                       
C *                                                                    
C *Parameters                                                         
INCLUDE(common/dft_parameters)
C *                                                                   
C *In variables                                                       
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_mol_info)
INCLUDE(../m4/common/errcodes)
      integer natyp, in
      integer ntrad
C *Out variables                                                     
      REAL prpt(natyp,*)
      REAL prwt(natyp,*)
C *Local variables
      integer latype,lrpt,nt,iat,irow

      nt=ntrad+1
      do lrpt=1,ntrad
         prpt(in,lrpt)=
     &        (dble(lrpt)/dble(nt-lrpt))**2
         prwt(in,lrpt)=
     &        (dble(nt)*2.0d0*dble(lrpt)**5)/(dble(nt-lrpt)**7)
      enddo
c
c...  Formally the following weight adjustments are correct only if
c...  the number of points exceeds 1. However a radial grid with
c...  only 1 point will yield such poor results that multiplying
c...  the weight by 2.25 in stead of 2 is not worth the extra
c...  logic to solve this.
c
      prwt(in,1)     = 1.5d0*prwt(in,1)
      prwt(in,ntrad) = 1.5d0*prwt(in,ntrad)
c
      return
      end

      subroutine prelog(natyp,in,ntrad,radm,prpt,prwt)
c **********************************************************************
c *   This subroutine computes the Log transformed radial quadrature   *
c *   grid due to Mura and Knowles [1] for all atom types.             *
c *   No screening is applied to the quadrature points.                *
c *                                                                    *
c *   [1] Michael E. Mura, Peter J. Knowles                            *
c *       Improved radial grids for quadrature in molecular density-   *
c *       functional calculations                                      *
c *       Journal of Chemical Physics, Vol. 104, (1996)                *
c *       pages 9848-9858                                              *
c **********************************************************************
      implicit none
C *********************************************************************
C *Declarations                                                       
C *                                                                    
C *Parameters                                                         
INCLUDE(common/dft_parameters)
C *                                                                   
C *In variables                                                       
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_mol_info)
INCLUDE(../m4/common/errcodes)
      integer natyp
      integer in
      integer ntrad
      REAL radm
C *Out variables                                                     
      REAL prpt(natyp,*)
      REAL prwt(natyp,*)
C *Local variables                                                    
      integer lrpt, iat
      REAL rpt, n2, ln2, n2m, i2p1m1, i2p1m, m
      REAL x, r, dr
c
      r(x,m)=-log(1.0d0-x**m)
      dr(x,m)=m*x**(m-1.0d0)/(1.0d0-x**m)
c
      m     = radm
      do lrpt=0,ntrad-1
         x = dble(2*lrpt+1)/dble(2*ntrad)
         rpt = r(x,m)
         prpt(in,lrpt+1)=rpt
         prwt(in,lrpt+1)=rpt*rpt*dr(x,m)/dble(ntrad)
      enddo 
c
      return
      end


      subroutine print_angangpt(gtyp,igrid,iout)
      implicit none
      integer gtyp,igrid,iout
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_parameters)
ccccINCLUDE(common/dft_basis)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
      integer i,j,k, atom_num_by_grid
      character*4 yelm
      REAL r
      write(iout,600)
chvd  r = srad(atom_num_by_grid(gtyp))
      r = 1.0d0
      write(iout,610)(angpt_radzn_num(k,gtyp,igrid),
     1   bnd_radzn(k,gtyp,igrid)*r,k=1,radzones_num(gtyp,igrid)-1),
     1   angpt_radzn_num(radzones_num(gtyp,igrid),gtyp,igrid)
      write(iout,620)prune_atom_radius(gtyp,igrid)
      write(iout,630)weight_atom_radius(gtyp,igrid)
      if(screen_atom_radius(gtyp,igrid).ne.dble(DFT_UNDEF)) then
         write(iout,640)screen_atom_radius(gtyp,igrid)
      endif
600   format(1x,'No ang. grid points for interval ',
     1          '(boundaries are fractions of atom radius):')
610   format(100(3x,10(i5,f6.3),/))
620   format(1x,
     1'Atom size for radius dependent angular grid pruning   = ',
     !f15.9,' Bohr')
630   format(1x,
     1'Atom size for weighting scheme atomic size adjustments= ',
     1f15.9,' Bohr')
640   format(1x,
     1'Atom size for radial grid screening                   = ',
     1f15.9,' Bohr')
      end


      subroutine mhl_prune(gtyp,it)
      implicit none
      integer gtyp,it
c
c     This subroutine sorts out the radial zones and their number
c     of grid points for a specified grid in accordance with the
c     pruning scheme suggested by Murray, Handy and Laming [1].
c
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
c
      integer mxl
      integer mxang, mxtht, mxphi, iat, i, j, ileb
      integer mnang, mntht, mnphi
      integer ao_tag
      parameter(ao_tag=1)
      REAL ratm, rphi, rtht
c
      integer atom_num_by_grid
      REAL srad
c
      REAL ktht
      parameter (ktht=5.0d0)
      integer maxleb
      parameter (maxleb=32)
      integer leba(maxleb)
      data leba/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,
     1          590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,
     1          3890,4334,4802,5294,5810/
c
      iat  = atom_num_by_grid(gtyp)
      ratm = prune_atom_radius(gtyp,it)
      mxl  = 0
      do i = 1, natoms
         if (gtype_num(i).eq.gtyp) then
            mxl = max(mxl,BL_maxang_on_atom(ao_tag,i))
         endif
      enddo
c
      if (ang_grid_scheme(gtyp,it).eq.AG_LEB) then
c
c        Sort out Lebedev grid
c
         if (mxl.le.2) then
            mnang = 6
            ileb  = 1
         else if (mxl.le.3) then
            mnang = 14
            ileb  = 2
         else if (mxl.le.4) then
            mnang = 26
            ileb  = 3
         else if (mxl.le.5) then
            mnang = 38
            ileb  = 4
         else
            mnang = 50
            ileb  = 5
         endif
c
c        Pick up the maximum number of radial grid points.
c        I.e. the size of the outer most angular shell.
c
         mxang = angpt_radzn_num(radzones_num(gtyp,it),gtyp,it)
         angpt_radzn_num(1,gtyp,it) = mnang
         radzones_num(gtyp,it)      = 1
         do i = ileb+1, maxleb
            if (leba(i).le.mxang) then
               if (radzones_num(gtyp,it).ge.maxradzn) then
                  call caserr('mhl_prune: too many radial zones!!!')
               endif
               bnd_radzn(radzones_num(gtyp,it),gtyp,it)=
     +            dsqrt(dble(leba(i-1))/mxang)*ratm/ktht
               radzones_num(gtyp,it)=radzones_num(gtyp,it)+1
               angpt_radzn_num(radzones_num(gtyp,it),gtyp,it)=leba(i)
            endif
         enddo
c
      else if (ang_grid_scheme(gtyp,it).eq.AG_LEG) then
c
c        Sort out Gauss-Legendre grid
c
         if (mxl.le.2) then
            mntht = 3
            mnphi = 6
         else 
            mntht = 2*mxl-1
            mnphi = 2*mntht
         endif
         mxtht = thetpt_radzn_num(1,gtyp,it)
         mxphi = phipt_radzn_num(1,gtyp,it)
         angpt_radzn_num(1,gtyp,it)  = mnphi*mntht
         thetpt_radzn_num(1,gtyp,it) = mntht
         phipt_radzn_num(1,gtyp,it)  = mnphi
         radzones_num(gtyp,it)       = 1
 10      if (mnphi.lt.mxphi.or.mntht.lt.mxtht) then
            if (radzones_num(gtyp,it).ge.maxradzn) then
               call caserr('mhl_prune: too many radial zones!!!')
            endif
            rphi = mnphi*ratm/(mxphi*ktht)
            rtht = mntht*ratm/(mxtht*ktht)
            if (mnphi.ge.mxphi.or.rtht.le.rphi) then
               bnd_radzn(radzones_num(gtyp,it),gtyp,it)=rtht
               radzones_num(gtyp,it)=radzones_num(gtyp,it)+1
               mntht=mntht+1
               thetpt_radzn_num(radzones_num(gtyp,it),gtyp,it)=mntht
               phipt_radzn_num(radzones_num(gtyp,it),gtyp,it) =mnphi
               angpt_radzn_num(radzones_num(gtyp,it),gtyp,it) =mnphi*
     &                                                         mntht
            else if (mntht.ge.mxtht.or.rtht.gt.rphi) then
               bnd_radzn(radzones_num(gtyp,it),gtyp,it)=rphi
               radzones_num(gtyp,it)=radzones_num(gtyp,it)+1
               mnphi=mnphi+1
               thetpt_radzn_num(radzones_num(gtyp,it),gtyp,it)=mntht
               phipt_radzn_num(radzones_num(gtyp,it),gtyp,it) =mnphi
               angpt_radzn_num(radzones_num(gtyp,it),gtyp,it) =mnphi*
     &                                                         mntht
            endif
            goto 10
         endif
      else
         call caserr('mhl_prune')
      endif
c
      end


      subroutine sg1_prune(gtyp,it,iwhich)
      implicit none
      integer gtyp, it, iwhich
c
c     This subroutine sorts out the radial zones for the specified
c     grid type in accordance with the SG1 grid specification [1].
c     However, once upon a time there was a modified specification
c     of this grid that was used as the default grid. The iwhich
c     variable can be used to select which one will be used. 
c     iwhich = 1, gives the real SG1 grids [1],
c     iwhich = 2, gives the modified version.
c     The modified version is only supported so that ancient 
c     calculations can be reproduced...
c
c     [1] P.M.W. Gill, B.G. Johnson, J.A. Pople
c         "A standard grid for density functional calculations"
c         Chem.Phys.Lett. Vol 209, No. 5,6, 1993, pp. 506-512.
c
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
c
      integer row_by_atomnum, atom_num_by_grid
      REAL ratm
c
      integer iat, irow, i, j
c
      integer acc(2,5)
      REAL alpha(3,4)
      data (alpha(1,i),i=1,4)/0.25d0, 0.5d0, 1.0d0, 4.50/
      data (alpha(2,i),i=1,4)/0.1667d0, 0.5d0, 0.9d0, 3.5d0/
      data (alpha(3,i),i=1,4)/0.1d0, 0.4d0, 0.8d0, 2.5d0/
      data (acc(1,i),i=1,5)/6, 38, 86, 194, 86/
c     second set should give more accurate results (but is not conform
c     the specification).
      data (acc(2,i),i=1,5)/38, 50, 110, 194, 110/
c
      if (iwhich.lt.1.or.iwhich.gt.2) then
         call caserr('sg1_prune: illegal value for iwhich')
      endif
      iat = atom_num_by_grid(gtyp)
      ratm = prune_atom_radius(gtyp,it)
      irow = min(3,row_by_atomnum(iat))
c
      angpt_radzn_num(1,gtyp,it) = acc(iwhich,1)
      do i = 1, 4
         angpt_radzn_num(1+i,gtyp,it) = acc(iwhich,1+i)
         bnd_radzn(i,gtyp,it)         = alpha(irow,i)*ratm
      enddo
      radzones_num(gtyp,it) = 5
c
      end


      subroutine print_atmgrid(igrid,iout)
      implicit none
      integer igrid,iout
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_xc)
      integer i,j,k,gtyp
      integer labels(max_atom), nlab, nlab1, nlab2
      character*4 yelm
      character*16 atxt,btxt
      REAL srad,r,SG1rad
c     Loop over all elements
c     element 0 is a BQ centre
      do i = 0, 118
         do 30 j = 1, ngridcentres
            gtyp = gtype_num(j)
            if (ian(j).eq.i.and.gtyp.gt.0) then
               do k = 1, j-1
                  if (gtyp.eq.gtype_num(k)) goto 30
               enddo
               call ztoname(i,yelm)
               r = srad(i)
               nlab = 0
               do k = j, ngridcentres
                  if (gtyp.eq.gtype_num(k)) then
                     nlab=nlab+1
                     labels(nlab)=k
                  endif
               enddo
               write(iout,600)yelm
               write(iout,605)
               nlab1=1
               nlab2=min(nlab,10)
 10            if (nlab1.le.nlab) then
                  write(iout,610)(labels(k),k=nlab1,nlab2)
                  nlab1=nlab2+1
                  nlab2=min(nlab,nlab2+10)
                  goto 10
               endif
               if (screen_sw) then
                  if (rad_grid_scheme(gtyp,igrid).eq.RG_SG1) then
                     write(iout,620)'Standard Grid 1',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid),
     +                              psitol(gtyp,igrid)
                  else if (rad_grid_scheme(gtyp,igrid).eq.RG_EML) then
                     write(iout,620)'Euler-MacLaurin',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid),
     +                              psitol(gtyp,igrid)
                  else if (rad_grid_scheme(gtyp,igrid).eq.RG_MK) then
                     write(iout,630)'Logarithmic    ',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid),
     +                              radm_num(gtyp,igrid),
     +                              psitol(gtyp,igrid)
                     if (rad_scale_scheme.eq.SC_MK) then
                        write(iout,*)'Mura-Knowles scale factors'
                     else if (rad_scale_scheme.eq.SC_GAM1) then
                        write(iout,*)'3.3*Bragg-Slater radius scale ',
     +                               'factors'
                     else if (rad_scale_scheme.eq.SC_GAM2) then
                        write(iout,*)'optimised scale factors'
                     else 
                        write(iout,*)'huh? rad_scale_scheme = ',
     +                               rad_scale_scheme
                     endif
                  else
                     call caserr('Print_atmgrid: Invalid radial grid')
                  endif
               else
                  if (rad_grid_scheme(gtyp,igrid).eq.RG_SG1) then
                     write(iout,625)'Standard Grid 1',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid)
                  else if (rad_grid_scheme(gtyp,igrid).eq.RG_EML) then
                     write(iout,625)'Euler-MacLaurin',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid)
                  else if (rad_grid_scheme(gtyp,igrid).eq.RG_MK) then
                     write(iout,635)'Logarithmic    ',
     +                              grid_scale(gtyp,igrid),
     +                              grid_atom_radius(gtyp,igrid),
     +                              radpt_num(gtyp,igrid),
     +                              radm_num(gtyp,igrid)
                  else
                     call caserr('Print_atmgrid: Invalid radial grid')
                  endif
               endif
               if (ang_grid_scheme(gtyp,igrid).eq.AG_LEG) then
                  atxt = 'Gauss-Legendre'
               else if (ang_grid_scheme(gtyp,igrid).eq.AG_LEB) then
                  atxt = 'Lebedev'
               else
                  call caserr('Print_atmgrid: Invalid angular grid')
               endif
               if (ang_prune_scheme(gtyp,igrid).eq.AP_SG1) then
                  btxt = 'SG1'
                  write(iout,650)atxt,btxt
                  call print_angangpt(gtyp,igrid,iout)
               else if (ang_prune_scheme(gtyp,igrid).eq.AP_SG1a) then
                  btxt = 'modSG1'
                  write(iout,650)atxt,btxt
                  call print_angangpt(gtyp,igrid,iout)
c              else if (ang_prune_scheme(gtyp,igrid).eq.AP_NONE) then
c                 btxt = 'no'
c                 write(iout,640)atxt,btxt,angupt_num(gtyp)
               else if (ang_prune_scheme(gtyp,igrid).eq.AP_MHL) then
                  btxt = 'MHL'
                  write(iout,650)atxt,btxt
                  call print_angangpt(gtyp,igrid,iout)
               else if (ang_prune_scheme(gtyp,igrid).eq.AP_AUTO) then
                  btxt = 'auto'
                  write(iout,650)atxt,btxt
                  call print_angangpt(gtyp,igrid,iout)
               else if (ang_prune_scheme(gtyp,igrid).eq.AP_RADZONE) then
                  btxt = 'manual'
                  write(iout,650)atxt,btxt
                  call print_angangpt(gtyp,igrid,iout)
               else
                  call caserr('Print_atmgrid: Invalid pruning scheme')
               endif
            endif
30       continue
      enddo
600   format(1x,a3)
605   format(1x,'Atom numbers:')
610   format(3x,10i4)
620   format(1x,'Radial grid : ',a14,' scale=',f8.4,' alpha=',f8.4,
     +          ' npts=',i4,' psitol=',e7.1)
625   format(1x,'Radial grid : ',a14,' scale=',f8.4,' alpha=',f8.4,
     +          ' npts=',i4)
630   format(1x,'Radial grid : ',a14,' scale=',f8.4,' alpha=',f8.4,
     +          ' npts=',i4,' m=',f8.4,' psitol=',e7.1)
635   format(1x,'Radial grid : ',a14,' scale=',f8.4,' alpha=',f8.4,
     +          ' npts=',i4,' m=',f8.4)
640   format(1x,'Angular grid: ',a14,' using ',a6,' pruning, npts=',i6)
650   format(1x,'Angular grid: ',a14,' using ',a6,' pruning')
      end

_IF(qmmm)
      REAL function dij(i,j)
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
      integer i,j
      REAL xt,yt,zt
      xt=atom_c(i,1)-atom_c(j,1)
      yt=atom_c(i,2)-atom_c(j,2)
      zt=atom_c(i,3)-atom_c(j,3)
      dij=sqrt(xt*xt+yt*yt+zt*zt)
      return
      end 
_ELSE
      subroutine dijcalc
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
      integer i,j
      REAL xt,yt,zt
      do i=1,ngridcentres
        do j=1,i
          xt=atom_c(i,1)-atom_c(j,1)
          yt=atom_c(i,2)-atom_c(j,2)
          zt=atom_c(i,3)-atom_c(j,3)
          dij(i,j)=sqrt(xt*xt+yt*yt+zt*zt)
          dij(j,i)=dij(i,j)
        enddo
      enddo
      return
      end 
_ENDIF

      function SG1rad(natomic,iout,extwr_sw)
C **********************************************************************
C *Description:				Version: 1.0                   *
C *Radii sizes for Euler-Maclaurin radial integration taken from SG1   *
C *	grid of Pople et. al Chem. Phys. Lett. 209 (1993) 506.         *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *In variables                                                        *
      integer natomic
      integer iout
      logical extwr_sw
C *Out function                                                        *
      REAL SG1rad
C *Local variables                                                     *
      REAL SG1radb(18)
c *function
      logical opg_root
C *End declarations                                                    *
C **********************************************************************
      data SG1radb/1.0000d0, 0.5882d0, 3.0769d0, 2.0513d0, 1.5385d0,
     +             1.2308d0, 
     +             1.0256d0, 0.8791d0, 0.7692d0, 0.6838d0, 4.0909d0, 
     +             3.1579d0, 
     +             2.5714d0, 2.1687d0, 1.8750d0, 1.6514d0, 1.4754d0, 
     +             1.3333d0/
      if(natomic .eq. 0)then
c
c Apply Carbon grid for BQ centres
c
         SG1rad=SG1radb(6)
      else if(natomic .gt. 18)then
        SG1rad = 5.0d0
        if(opg_root().and.extwr_sw)
     +       write(iout,*)
     + '************** WARNING arbitrary SG1 radius of 5 used for z=',
     +       natomic
      else
        SG1rad=SG1radb(natomic)
      endif

      return
      end 

      REAL function mkrad(natomic,iout,extwr_sw)
c **********************************************************************
c *   Description:                                                     *
c *   Radii sizes for logarithmic radial integration as recommended    *
c *   by Mura and Knowles:                                             *
C *   [1] Mura, Knowles, J.Chem.Phys.,Vol.104 (1996) pages 9849        *
c **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *In variables                                                        *
      integer natomic
      integer iout
      logical extwr_sw, opg_root
C *End declarations                                                    *
C **********************************************************************
      mkrad = -1.0d0
      if (natomic.lt.0) then
         call caserr("mkrad failed: negative atomic number")
      else if (natomic.eq.0) then
c        BQ center: use Carbon parameters
         mkrad = 5.0d0
      else if (natomic.le.2) then
c        H-He
         mkrad = 5.0d0
      else if (natomic.le.4) then
c        Li-Be
         mkrad = 7.0d0
      else if (natomic.le.10) then
c        B-Ne
         mkrad = 5.0d0
      else if (natomic.le.12) then
c        Na-Mg
         mkrad = 7.0d0
      else if (natomic.le.18) then
c        Al-Ar
         mkrad = 5.0d0
      else if (natomic.le.20) then
c        K-Ca
         mkrad = 7.0d0
      else if (natomic.le.36) then
c        Sc-Kr
         mkrad = 5.0d0
      else if (natomic.le.38) then
c        Rb-Sr
         mkrad = 7.0d0
      else if (natomic.le.54) then
c        Y-Xe
         mkrad = 5.0d0
      else if (natomic.le.56) then
c        Cs-Ba
         mkrad = 7.0d0
      else if (natomic.le.86) then
c        La-Rn
         mkrad = 5.0d0
      else if (natomic.le.88) then
c        Fr-Ra
         mkrad = 7.0d0
      else if (natomic.le.118) then
c        Ac-???
         mkrad = 5.0d0
      else
         call caserr("mkrad failed: natomic exceeds 118")
      endif
      if (natomic.gt.30) then
         if(opg_root().and.extwr_sw)
     +   write(iout,*)'************** WARNING unpublished radius ',
     +                'used for z=',natomic
      endif
      end

      REAL function lograd(natomic,iout,extwr_sw)
C **********************************************************************
C *   Description:				Version: 1.0           *
C *   Radii sizes for Logarithmic radial integration due to Mura et al.*
C *   [1] Mura, Knowles, J.Chem.Phys.,Vol.104 (1996) pages 9849        *
c *                                                                    *
c *   Radii established from atomic tests: see RESULT_alpha            *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *In variables                                                        *
      integer natomic
      integer iout
      logical extwr_sw
C *Local variables                                                     *
      REAL logradb(118)
c *function
      logical opg_root
C *End declarations                                                    *
C **********************************************************************
      data logradb/
c      H       He
     + 9.46d0, 5.46d0, 
c      Li              B                                       Ne
     +10.92d0, 9.05d0,6.290d0, 5.63d0, 4.84d0, 4.25d0, 3.92d0, 3.63d0,
c      Na              Al                                      Ar
     + 9.25d0, 7.54d0, 6.42d0, 5.75d0, 5.50d0, 5.04d0, 4.71d0, 4.42d0,
c      K       Ca
     + 8.80d0, 7.88d0, 
c      Sc                              Mn
     + 8.00d0, 7.34d0, 7.25d0, 6.79d0, 6.54d0,
c      Fe                              Zn
     + 6.42d0, 5.88d0, 5.75d0, 5.54d0, 5.34d0,
c      Ga                                      Kr
     + 6.29d0, 5.84d0, 5.42d0, 5.25d0, 5.21d0, 4.92d0,
c      Rb      Sr
     + 8.80d0, 7.88d0,
c      Y                               Tc
     + 8.00d0, 7.34d0, 7.25d0, 6.79d0, 6.54d0,
c      Ru                              Cd
     + 6.42d0, 5.88d0, 5.75d0, 5.54d0, 5.34d0,
c      In                                      Xe
     + 6.29d0, 5.84d0, 5.42d0, 5.25d0, 5.21d0, 4.92d0,
c      Cs      Ba      La
     + 8.80d0, 7.88d0, 8.00d0,
c      Ce                                              Gd
     + 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0,
c      Tb                                              Lu
     + 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0,
c      Hf                      Re
     + 7.34d0, 7.25d0, 6.79d0, 6.54d0, 
c      Os                              Hg
     + 6.42d0, 5.88d0, 5.75d0, 5.54d0, 5.34d0,
c      Tl                                      Rn
     + 6.29d0, 5.84d0, 5.42d0, 5.25d0, 5.21d0, 4.92d0,
c      Fr      Ra      Ac
     + 8.80d0, 7.88d0, 8.00d0,
c      Th                                              Cm
     + 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0,
c      Bk                                              Lr
     + 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0, 5.00d0,
c      Unq                     Uns
     + 7.34d0, 7.25d0, 6.79d0, 6.54d0,
c      Uno                             112
     + 6.42d0, 5.88d0, 5.75d0, 5.54d0, 5.34d0,
c      113                                     118
     + 6.29d0, 5.84d0, 5.42d0, 5.25d0, 5.21d0, 4.92d0/
      lograd = 1.0d0
      if(natomic.eq.0)then
c
c Carbon grid for BQ centres
c
         lograd = logradb(6)
      else if (0.lt.natomic.and.natomic.le.36) then
         lograd = logradb(natomic)
      else if (0.lt.natomic) then
         lograd = logradb(natomic)
         if(opg_root().and.extwr_sw)
     +  write(iout,*)'************** WARNING arbitrary log radius of ',
     +               logradb(natomic),' used for z=',natomic
      else
         call caserr('negative atomic number in lograd !!!')
      endif
      return
      end 

      function srad(natomic)
C **********************************************************************
C * Description:                        Version: 1.0                   *
C *                                                                    *
C * Returns the Bragg-Slater radius [1] of an element.                 *
C * For the noble gasses no Bragg-Slater radii exist and the Clementi  *
C * radii [2] are returned.                                            *
C * For the elements Pb and At the Clementi radii [2] were used also.  *
C * For the element Fr the Bragg-Slater radius [1] of Cs was used.     *
C * For the hydrogen atom the Bragg-Slater radius (0.25 Angstrom) is   *
C * very small, and instead of using the value of 0.35 Angstrom        *
C * proposed by Becke [3] we choose to use 1 Bohr in accordance with   *
C * the observation that the radius of an atom is approximately equal  *
C * to the radius at which the radial distribution function of the     *
C * valence orbital has its maximum value [1].                         *
C * For the elements Cm to Lr we took the values of the elements Gd to *
C * Lu [1] minus the difference between the values of Eu and Am.       *
C *                                                                    *
C * [1] "Atomic Radii in Crystals", J.C. Slater                        *
C *     J. Chem. Phys. Vol. 41, No. 10, 1964, pages 3199-3204          *
C *                                                                    *
C * [2] "Atomic Screening Constants from SCF functions                 *
C *      II Atoms with 37 to 86 electrons                              *
C *     E. Clementi, D. L. Raimondi, W. P. Reinhardt                   *
C *     J. Chem. Phys. Vol. 47, No 4, 1967, pages 1300-1307            *
C *                                                                    *
C * [3] "A multicenter numerical integration scheme for polyatomic     *
C *      molecules", A.D. Becke                                        *
C *     J. Chem. Phys. Vol. 88 No 4, 1988, pages 2547-2553             *
C *     Appendix                                                       *
C *                                                                    *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *Parameters                                                          *
INCLUDE(common/dft_physical_constants)
C *In variables                                                        *
      integer natomic
C *Out variables                                                       *
      REAL srad
      REAL slatrad(103)
C *Local variables                                                     *
      integer i
C *End declarations                                                    *
C **********************************************************************
      data(slatrad(i),i=1,103)
C      H       He                      B
     */0.53d0, 0.31d0, 1.45d0, 1.05d0, 0.85d0, 
c    */0.53d0, 0.31d0, 0.598916d0, 1.05d0, 0.85d0, 
C      C                               Ne
     * 0.70d0, 0.65d0, 0.60d0, 0.50d0, 0.38d0, 
c    * 0.70d0, 0.65d0, 0.60d0, 0.952788d0, 0.38d0, 
C      Na                              P
     * 1.80d0, 1.50d0, 1.25d0, 1.10d0, 1.00d0, 
C      S               Ar              Ca
     * 1.00d0, 1.00d0, 0.71d0, 2.20d0, 1.80d0, 
C      Sc                              Mn
     * 1.60d0, 1.40d0, 1.35d0, 1.40d0, 1.40d0, 
C      Fe                              Zn
     * 1.40d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, 
C      Ga                              Br
     * 1.30d0, 1.25d0, 1.15d0, 1.15d0, 1.15d0,
C      Kr                              Zr
     * 0.88d0, 2.35d0, 2.00d0, 1.80d0, 1.55d0, 
C      Nb                              Rh
     * 1.45d0, 1.45d0, 1.35d0, 1.30d0, 1.35d0, 
C      Pd                              Sn
     * 1.40d0, 1.60d0, 1.55d0, 1.55d0, 1.45d0, 
C      Sb                      Xe      Cs
     * 1.45d0, 1.40d0, 1.40d0, 1.08d0, 2.60d0, 
C      Ba                              Nd
     * 2.15d0, 1.95d0, 1.85d0, 1.85d0, 1.85d0, 
C      Pm                              Tb
     * 1.85d0, 1.85d0, 1.85d0, 1.80d0, 1.75d0, 
C      Dy                              Yb
     * 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,
C      Lu                              Re
     * 1.75d0, 1.55d0, 1.45d0, 1.35d0, 1.35d0, 
C      Os                              Hg
     * 1.30d0, 1.35d0, 1.35d0, 1.35d0, 1.50d0, 
C      Tl                              At
     * 1.90d0, 1.54d0, 1.60d0, 1.90d0, 1.27d0, 
C      Rn                              Th
     * 1.20d0, 2.60d0, 2.15d0, 1.95d0, 1.80d0, 
C      Pa                              Am
     * 1.80d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 
C      Cm                              Fm
     * 1.70d0, 1.65d0, 1.65d0, 1.65d0, 1.65d0, 
C      Md              Lr
     * 1.65d0, 1.65d0, 1.65d0/
c
      if (natomic .lt. 0 .or. natomic .gt. 103) then
         call caserr("srad: illegal atomic number, must be in 0-103")
      else if (natomic .eq. 0) then
c        This is a BQ centre
c        Use Carbon radius as a default
         srad = slatrad(6)/atob
      else
         srad=slatrad(natomic)/atob
      endif
      return
      end


c
C     ******************************************************************
c     *
c     * kmaddcs
C     * Add contibution to Kohn-Sham matrix                            *
c     *
C     ******************************************************************
c
      subroutine kmaddcs(bfn_val,bfng_val,wt,xc_vpt,xc_dvpt,xc_dtpt,
     &     grho,gradcorr_sw,kinetic_sw,kma,kmb,rks_sw,om_scr,om_scr1,
     &     nao,npts,mxp,sma)

      implicit none
c
c     Description:
c
c     Calculate the contributions to the Kohn-Sham matrix in AO basis
c     and without screening.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dfder)
c
c     Inputs:
c
      logical gradcorr_sw,kinetic_sw,rks_sw
      integer nao,npts,mxp
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL wt(npts)
      REAL xc_vpt(mxp,2)
      REAL xc_dvpt(mxp,3) 
      REAL xc_dtpt(mxp,2) 
      REAL grho(mxp,2,3)
c
c     Input/Output:
c
      REAL kma(nao*(nao+1)/2), kmb(nao*(nao+1)/2)
c
c     Output:
c
      REAL sma(nao*(nao+1)/2)
c
c     Workspace:
c
      REAL om_scr(nao), om_scr1(mxp,3)
c
c     Local:
c
      integer i,j,n
      integer ipt
      REAL T
c
c     Code:
c
_IF(debug_S)
      n=0 
      do i=1,nao
         do j=1,i
            n=n+1
            do ipt=1,npts
               sma(n)= sma(n) + 
     &              wt(ipt)*bfn_val(ipt,i)*bfn_val(ipt,j)
            enddo
         enddo
      enddo
_ENDIF
      if(.not.gradcorr_sw) then
c
c        Local density functional
c
         n=0 
         do i=1,nao
            do j=1,i
               n=n+1
               do ipt=1,npts
                  kma(n)= kma(n) + wt(ipt)*xc_vpt(ipt,ira)*
     &                             bfn_val(ipt,i)*bfn_val(ipt,j)
               enddo
            enddo
         enddo
         if (.not.rks_sw) then
            n=0 
            do i=1,nao
               do j=1,i
                  n=n+1
                  do ipt=1,npts
                     kmb(n)= kmb(n) + wt(ipt)*xc_vpt(ipt,irb)*
     &                                bfn_val(ipt,i)*bfn_val(ipt,j)
                  enddo
               enddo
            enddo
         endif
c
      else
c
c        Gradient corrected functionals
c
         if (rks_sw) then
            do ipt = 1, npts
               om_scr1(ipt,1)=0.5d0*grho(ipt,1,1)*xc_dvpt(ipt,igaa)
               om_scr1(ipt,2)=0.5d0*grho(ipt,1,2)*xc_dvpt(ipt,igaa)
               om_scr1(ipt,3)=0.5d0*grho(ipt,1,3)*xc_dvpt(ipt,igaa)
            enddo
         else
            do ipt = 1, npts
               om_scr1(ipt,1)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,1))
               om_scr1(ipt,2)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,2))
               om_scr1(ipt,3)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,3))
            enddo
         endif
c
         do ipt = 1, npts
            do i=1,nao
               om_scr(i)=bfng_val(ipt,i,1)*om_scr1(ipt,1) +
     &                   bfng_val(ipt,i,2)*om_scr1(ipt,2) +
     &                   bfng_val(ipt,i,3)*om_scr1(ipt,3)
               om_scr(i)=wt(ipt)*om_scr(i)
            enddo
            n=0
            do i=1,nao
               T=wt(ipt)*xc_vpt(ipt,ira)*bfn_val(ipt,i)+om_scr(i)
               do j=1,i
                  n=n+1
                  kma(n)=kma(n)+ T*bfn_val(ipt,j) + 
     &                 bfn_val(ipt,i)*om_scr(j)
               enddo 
            enddo
         enddo
c
         if(.not.rks_sw)then
            do ipt = 1, npts
               om_scr1(ipt,1)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,1))
               om_scr1(ipt,2)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,2))
               om_scr1(ipt,3)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,3))
            enddo
c
            do ipt = 1, npts
               n=0
               do i=1,nao
                  om_scr(i)=bfng_val(ipt,i,1)*om_scr1(ipt,1) +
     &                      bfng_val(ipt,i,2)*om_scr1(ipt,2) +
     &                      bfng_val(ipt,i,3)*om_scr1(ipt,3)
                  om_scr(i)=wt(ipt)*om_scr(i)
               enddo
               do i=1,nao
                  T=wt(ipt)*xc_vpt(ipt,irb)*bfn_val(ipt,i)+om_scr(i)
                  do j=1,i
                     n=n+1
                     kmb(n)= kmb(n)+T*bfn_val(ipt,j) + 
     &                    bfn_val(ipt,i)*om_scr(j)
                  enddo 
               enddo
            enddo
         endif
      endif
      if (kinetic_sw) then
         n=0 
         do i=1,nao
            do j=1,i
               n=n+1
               do ipt=1,npts
                  kma(n)= kma(n)+wt(ipt)*xc_dtpt(ipt,ita)*
     &                          (bfng_val(ipt,i,1)*bfng_val(ipt,j,1)+
     &                           bfng_val(ipt,i,2)*bfng_val(ipt,j,2)+
     &                           bfng_val(ipt,i,3)*bfng_val(ipt,j,3))
               enddo
            enddo
         enddo
         if (.not.rks_sw) then
            n=0 
            do i=1,nao
               do j=1,i
                  n=n+1
                  do ipt=1,npts
                     kmb(n)= kmb(n)+wt(ipt)*xc_dtpt(ipt,itb)*
     &                       (bfng_val(ipt,i,1)*bfng_val(ipt,j,1)+
     &                        bfng_val(ipt,i,2)*bfng_val(ipt,j,2)+
     &                        bfng_val(ipt,i,3)*bfng_val(ipt,j,3))
                  enddo
               enddo
            enddo
         endif
      endif
      return
      end
c
C     ****************************************************************
c     *                                                              *
C     * Description:				Version: 1.0         *
C     * Add contibution to Kohn-Sham matrix, based on screening      *
c     *                                                              *
C     ****************************************************************
c

      subroutine kmaddcs_scr(bfn_val,bfng_val,wt,xc_vpt,xc_dvpt,
     &     xc_dtpt,grho,gradcorr_sw,kinetic_sw,kma,kmb,rks_sw,nao,
     &     npts,mxp,active_bfn_list,n_active_bfn,om_scr,om_scr1,T,sma)

      implicit none
c
c     Description:
c
c     Calculate the contributions to the Kohn-Sham matrix in AO basis
c     and with screening.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dfder)
c
c     Inputs:
c
      logical gradcorr_sw, kinetic_sw, rks_sw
      integer npts,mxp,nao
      integer n_active_bfn, active_bfn_list(n_active_bfn)
      REAL bfn_val(mxp,n_active_bfn)
      REAL bfng_val(mxp,nao,3)
      REAL wt(npts)
      REAL xc_vpt(mxp,2)
      REAL xc_dvpt(mxp,3) 
      REAL xc_dtpt(mxp,2) 
      REAL grho(mxp,2,3)
c
c     Input/Output:
c
      REAL kma(nao*(nao+1)/2), kmb(nao*(nao+1)/2)
c
c     Output:
c
      REAL sma(nao*(nao+1)/2)
c
c     Workspaces:
c
      REAL om_scr(mxp,n_active_bfn), om_scr1(mxp,3), T(mxp)
c
c     Local:
c
      integer i,j,n,iki
      integer ipt
_IF(single)
      REAL sdot
_ELSEIF(hp700)
      REAL `vec_$ddot'
_ELSE
      REAL ddot
_ENDIF
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)
c
c     Code:
c
_IF(debug_S)
      do i=1,n_active_bfn
         do j=1,i
            n=iky(active_bfn_list(i))+active_bfn_list(j)
            do ipt=1,npts
               sma(n)= sma(n)+
     &              wt(ipt)*bfn_val(ipt,i)*bfn_val(ipt,j)
            enddo
         enddo
      enddo
_ENDIF
      if(.not.gradcorr_sw) then
c
c        Local density functional
c
         do i=1,n_active_bfn
            do j=1,i
               n=iky(active_bfn_list(i))+active_bfn_list(j)
               do ipt=1,npts
                  kma(n)= kma(n)+wt(ipt)*xc_vpt(ipt,ira)*
     &                           bfn_val(ipt,i)*bfn_val(ipt,j)
               enddo
            enddo
         enddo
         if (.not.rks_sw) then
            do i=1,n_active_bfn
               do j=1,i
                  n=iky(active_bfn_list(i))+active_bfn_list(j)
                  do ipt=1,npts
                     kmb(n)= kmb(n)+wt(ipt)*xc_vpt(ipt,irb)*
     &                              bfn_val(ipt,i)*bfn_val(ipt,j)
                  enddo
               enddo
            enddo
         endif
c
      else
c
c        Gradient corrected density functional
c
         if (rks_sw) then
            do ipt = 1, npts
               om_scr1(ipt,1)=0.5d0*grho(ipt,1,1)*xc_dvpt(ipt,igaa)
               om_scr1(ipt,2)=0.5d0*grho(ipt,1,2)*xc_dvpt(ipt,igaa)
               om_scr1(ipt,3)=0.5d0*grho(ipt,1,3)*xc_dvpt(ipt,igaa)
            enddo
         else
            do ipt = 1, npts
               om_scr1(ipt,1)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,1))
               om_scr1(ipt,2)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,2))
               om_scr1(ipt,3)=(2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,2,3))
            enddo
         endif
c
         do i=1,n_active_bfn
            do ipt = 1,npts
               om_scr(ipt,i)=bfng_val(ipt,i,1)*om_scr1(ipt,1) +
     &                       bfng_val(ipt,i,2)*om_scr1(ipt,2) +
     &                       bfng_val(ipt,i,3)*om_scr1(ipt,3)
               om_scr(ipt,i)=wt(ipt)*om_scr(ipt,i)
            enddo
         enddo
         do i=1,n_active_bfn
            do ipt=1,npts
               T(ipt) =wt(ipt)*xc_vpt(ipt,ira)*bfn_val(ipt,i)+
     &                 om_scr(ipt,i)
            enddo
            iki=iky(active_bfn_list(i))
            do j=1,i
               n=iki+active_bfn_list(j)
               kma(n)= kma(n)+ddot(npts,t,1,bfn_val(1,j),1) + 
     &                        ddot(npts,bfn_val(1,i),1,om_scr(1,j),1)
            enddo
         enddo

         if (.not.rks_sw)then
            do ipt = 1, npts
               om_scr1(ipt,1)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,1))
               om_scr1(ipt,2)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,2))
               om_scr1(ipt,3)=(2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                       +       xc_dvpt(ipt,igab)*grho(ipt,1,3))
            enddo
            do i=1,n_active_bfn
               do ipt = 1,npts
                  om_scr(ipt,i)=bfng_val(ipt,i,1)*om_scr1(ipt,1) +
     &                          bfng_val(ipt,i,2)*om_scr1(ipt,2) +
     &                          bfng_val(ipt,i,3)*om_scr1(ipt,3)
                  om_scr(ipt,i)=wt(ipt)*om_scr(ipt,i)
               enddo
            enddo
            do i=1,n_active_bfn
               do ipt=1,npts
                  T(ipt) =wt(ipt)*xc_vpt(ipt,irb)*bfn_val(ipt,i)+
     &                    om_scr(ipt,i)
               enddo
               iki=iky(active_bfn_list(i))
               do j=1,i
                  n=iki+active_bfn_list(j)
c                 do ipt=1,npts
c                    kmb(n)= kmb(n)+t(ipt)*bfn_val(ipt,j)+ 
c    &                              bfn_val(ipt,i)*om_scr(ipt,j)
c                 enddo
                  kmb(n)= kmb(n)+ddot(npts,t,1,bfn_val(1,j),1)+ 
     &                           ddot(npts,bfn_val(1,i),1,om_scr(1,j),1)
               enddo
            enddo
         endif
      endif
      if (kinetic_sw) then
         do i=1,n_active_bfn
            do j=1,i
               n=iky(active_bfn_list(i))+active_bfn_list(j)
               do ipt=1,npts
                  kma(n)= kma(n)+wt(ipt)*xc_dtpt(ipt,ira)*
     &                          (bfng_val(ipt,i,1)*bfng_val(ipt,j,1)+
     &                           bfng_val(ipt,i,2)*bfng_val(ipt,j,2)+
     &                           bfng_val(ipt,i,3)*bfng_val(ipt,j,3))
               enddo
            enddo
         enddo
         if (.not.rks_sw) then
            do i=1,n_active_bfn
               do j=1,i
                  n=iky(active_bfn_list(i))+active_bfn_list(j)
                  do ipt=1,npts
                     kmb(n)= kmb(n)+wt(ipt)*xc_dtpt(ipt,irb)*
     &                             (bfng_val(ipt,i,1)*bfng_val(ipt,j,1)+
     &                              bfng_val(ipt,i,2)*bfng_val(ipt,j,2)+
     &                              bfng_val(ipt,i,3)*bfng_val(ipt,j,3))
                  enddo
               enddo
            enddo
         endif
      endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine xc_forces_scr(ao_tag,nao,rks_sw,gradcorr_sw,kinetic_sw,
     &     adens,bdens,gwt_sw,gwt_avail_sw,
     &     bfn_val,bfng_val,bfn_hess,
     &     wt,gwt,xc_ept,xc_vpt,xc_dvpt,xc_dtpt,grho,
     &     scr2,scr3,scr4,
     &     near_atom_list, num_near_atoms,
     &     first_bf, active_bfn_list,n_active_bfn,
     &     active_bfn_atms,n_active_bfn_atm,
     &     grad,npts,mxp,iatom)
      implicit none
c
c     Description:
c
c     Calculate the first derivative of Exc with respect to the 
c     nuclear coordinates using screening.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dfder)
c
c     Inputs:
c
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_basis_api)
cINCLUDE(common/dft_module_comm)
      integer mxp, npts
      integer ao_tag, nao
      integer num_near_atoms
      integer near_atom_list(num_near_atoms)
      logical rks_sw,gradcorr_sw,kinetic_sw
      REAL adens(nao*(nao+1)/2), bdens(nao*(nao+1)/2)
      logical gwt_sw, gwt_avail_sw
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      REAL wt(mxp),gwt(3,mxp,num_near_atoms)
      REAL xc_ept(mxp),xc_vpt(mxp,2),xc_dvpt(mxp,3),xc_dtpt(mxp,2)
      REAL grho(mxp,2,3)
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer n_active_bfn_atm
      integer active_bfn_atms(n_active_bfn_atm)
      integer first_bf(*)
      integer iatom
c
c     Input/Output:
c
      REAL grad(3,*)
c
c     Workspace:
c
      REAL scr2(3,nao),scr3(nao),scr4(npts,3)
c
c     Local:
c
      integer alph,beta,spin,nspin
      integer latm,nu,mu,ixf,ixl
      integer atmt,bstart,bend
      REAL scr1(3),tmp
      REAL grad1(3)
      integer ipt
      integer nchk
      integer near
      character*1 xn, xt
      data xt,xn/'t','n'/
c
c     Code:
c
      alph = 1
      beta  = 2
c
      nspin = 1
      if(.not. rks_sw) nspin = 2
c
c     main loop over centres for which the current quad shell 
c     contributes to the gradient - limited to those that are close 
c     enough for the basis functions on the centre (mu loop) to be 
c     significant for points in this quadrature shell
c
c     add the contribution from the gradients of the weights
c
      if (gwt_sw.and.gwt_avail_sw) then
         do near=1,num_near_atoms
            latm = near_atom_list(near)
            do ipt = 1, npts
               grad(1,latm) = grad(1,latm) 
     +                      + gwt(1,ipt,near)*xc_ept(ipt)
               grad(2,latm) = grad(2,latm) 
     +                      + gwt(2,ipt,near)*xc_ept(ipt)
               grad(3,latm) = grad(3,latm) 
     +                      + gwt(3,ipt,near)*xc_ept(ipt)
            enddo
         enddo
      endif
c
c     NB screening is done twice as active_bfn_list will have no 
c     entries on a distant atom, but we save the effort of working 
c     this out).
c
      if (gradcorr_sw) then
         do spin = 1, nspin
            if (nspin.eq.1) then
               do ipt = 1, npts
                  scr4(ipt,1) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
                  scr4(ipt,2) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
                  scr4(ipt,3) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
               enddo
            else if (spin.eq.1) then
               do ipt = 1, npts
                  scr4(ipt,1) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,1))
                  scr4(ipt,2) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,2))
                  scr4(ipt,3) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,3))
               enddo
            else
               do ipt = 1, npts
                  scr4(ipt,1) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,1))
                  scr4(ipt,2) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,2))
                  scr4(ipt,3) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,3))
               enddo
            endif
            nchk  = 0
            bend  = 0
            do 100 near=1,n_active_bfn_atm
               latm = active_bfn_atms(near)
               atmt = BL_get_atom_type(ao_tag,latm)
               if (atmt.eq.0) goto 100

               bstart=first_bf(latm)
               bend = first_bf(latm+1) - 1
               call aclear_dp(grad1,3,0.0d0)
c
c              ixf, ixl are the indices into active_bfn_list for the
c              functions on this atom
c
               call atombf(active_bfn_list,n_active_bfn,bstart,bend,
     &                     ixf,ixl)
               if (latm.eq.iatom.and.gwt_sw) then
                  nchk = nchk + ixl - ixf + 1
                  goto 100
               endif
c
c              loop over basis functions on the current atom
c
               do mu=ixf,ixl          
c
c                 count basis functions treated
c
                  nchk = nchk + 1     
c
c                 construct a row of density matrix
c
                  if(spin .eq. 1)then
                     call gather_mu_scr(adens,nao,
     &                    active_bfn_list,n_active_bfn,mu,scr3) 
                  else
                     call gather_mu_scr(bdens,nao,
     &                    active_bfn_list,n_active_bfn,mu,scr3) 
                  endif

                  do ipt = 1,npts
                     scr1(1) = bfng_val(ipt,mu,1)*xc_vpt(ipt,spin)
                     scr1(2) = bfng_val(ipt,mu,2)*xc_vpt(ipt,spin)
                     scr1(3) = bfng_val(ipt,mu,3)*xc_vpt(ipt,spin)
c
                     scr1(1) = scr1(1) + 
     &                         scr4(ipt,1)*bfn_hess(ipt,mu,1)+
     &                         scr4(ipt,2)*bfn_hess(ipt,mu,4)+
     &                         scr4(ipt,3)*bfn_hess(ipt,mu,5)
                     scr1(2) = scr1(2) + 
     &                         scr4(ipt,1)*bfn_hess(ipt,mu,4)+
     &                         scr4(ipt,2)*bfn_hess(ipt,mu,2)+
     &                         scr4(ipt,3)*bfn_hess(ipt,mu,6)
                     scr1(3) = scr1(3) + 
     &                         scr4(ipt,1)*bfn_hess(ipt,mu,5)+
     &                         scr4(ipt,2)*bfn_hess(ipt,mu,6)+
     &                         scr4(ipt,3)*bfn_hess(ipt,mu,3)
c 
                     scr1(1) = scr1(1)*wt(ipt)
                     scr1(2) = scr1(2)*wt(ipt)
                     scr1(3) = scr1(3)*wt(ipt)
c     
                     do nu=1,n_active_bfn
                        scr2(1,nu)=scr1(1)*bfn_val(ipt,nu)
                        scr2(2,nu)=scr1(2)*bfn_val(ipt,nu)
                        scr2(3,nu)=scr1(3)*bfn_val(ipt,nu)
                     enddo
c     
                     do nu=1,n_active_bfn
                        tmp=scr4(ipt,1)*bfng_val(ipt,nu,1)+
     &                      scr4(ipt,2)*bfng_val(ipt,nu,2)+
     &                      scr4(ipt,3)*bfng_val(ipt,nu,3)
                        tmp=tmp*wt(ipt)
                        scr2(1,nu) = scr2(1,nu) + bfng_val(ipt,mu,1)*tmp
                        scr2(2,nu) = scr2(2,nu) + bfng_val(ipt,mu,2)*tmp
                        scr2(3,nu) = scr2(3,nu) + bfng_val(ipt,mu,3)*tmp
                     enddo
c
                     if (kinetic_sw) then
                        do nu=1,n_active_bfn
                           scr2(1,nu) = scr2(1,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,1)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,4)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,5))
                           scr2(2,nu) = scr2(2,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,4)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,2)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,6))
                           scr2(3,nu) = scr2(3,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,5)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,6)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,3))
                        enddo
                     endif
c     
                     call dgemv(xn,3,n_active_bfn,-2.0d0
     +                    ,scr2,3,scr3,1,1.0d0,grad1,1)
                  enddo               ! loop over quadrature points
               enddo
               grad(1,latm) = grad(1,latm)+grad1(1)
               grad(2,latm) = grad(2,latm)+grad1(2)
               grad(3,latm) = grad(3,latm)+grad1(3)
               if(gwt_sw) then
c
c...              If we use the gradients of the weights then we must 
c...              also take the gradients of the quadrature points into
c...              account. We do that using the translational invariance
c...              of the energy.
c
                  grad(1,iatom) = grad(1,iatom)-grad1(1)
                  grad(2,iatom) = grad(2,iatom)-grad1(2)
                  grad(3,iatom) = grad(3,iatom)-grad1(3)
               endif
 100        continue
         enddo
      else
         do spin = 1, nspin
            nchk  = 0
            bend  = 0
            do 200 near=1,n_active_bfn_atm
               latm = active_bfn_atms(near)
               atmt = BL_get_atom_type(ao_tag,latm)
               if (atmt.eq.0) goto 200
c
               bstart=first_bf(latm)
               bend = first_bf(latm+1) - 1
               call aclear_dp(grad1,3,0.0d0)
c
c              ixf, ixl are the indices into active_bfn_list for the
c              functions on this atom
c
               call atombf(active_bfn_list,n_active_bfn,bstart,bend,
     &                     ixf,ixl)
               if (latm.eq.iatom.and.gwt_sw) then
                  nchk = nchk + ixl - ixf + 1
                  goto 200
               endif
c
c              loop over basis functions on the current atom
c
               do mu=ixf,ixl          
c
c                 count basis functions treated
c
                  nchk = nchk + 1     
c
c                 construct a row of density matrix
c
                  if(spin .eq. 1)then
                     call gather_mu_scr(adens,nao,
     &                    active_bfn_list,n_active_bfn,mu,scr3) 
                  else
                     call gather_mu_scr(bdens,nao,
     &                    active_bfn_list,n_active_bfn,mu,scr3) 
                  endif
c
                  do ipt = 1,npts
                     scr1(1) = bfng_val(ipt,mu,1)*xc_vpt(ipt,spin)
                     scr1(2) = bfng_val(ipt,mu,2)*xc_vpt(ipt,spin)
                     scr1(3) = bfng_val(ipt,mu,3)*xc_vpt(ipt,spin)
c 
                     scr1(1) = scr1(1)*wt(ipt)
                     scr1(2) = scr1(2)*wt(ipt)
                     scr1(3) = scr1(3)*wt(ipt)
c     
                     do nu=1,n_active_bfn
                        scr2(1,nu)=scr1(1)*bfn_val(ipt,nu)
                        scr2(2,nu)=scr1(2)*bfn_val(ipt,nu)
                        scr2(3,nu)=scr1(3)*bfn_val(ipt,nu)
                     enddo
c
                     if (kinetic_sw) then
                        do nu=1,n_active_bfn
                           scr2(1,nu) = scr2(1,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,1)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,4)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,5))
                           scr2(2,nu) = scr2(2,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,4)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,2)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,6))
                           scr2(3,nu) = scr2(3,nu) +
     &                         wt(ipt)*xc_dtpt(ipt,spin)*
     &                         (bfng_val(ipt,nu,1)*bfn_hess(ipt,mu,5)+
     &                          bfng_val(ipt,nu,2)*bfn_hess(ipt,mu,6)+
     &                          bfng_val(ipt,nu,3)*bfn_hess(ipt,mu,3))
                        enddo
                     endif
c     
                     call dgemv(xn,3,n_active_bfn,-2.0d0
     +                    ,scr2,3,scr3,1,1.0d0,grad1,1)
                  enddo               ! loop over quadrature points
               enddo
               grad(1,latm) = grad(1,latm)+grad1(1)
               grad(2,latm) = grad(2,latm)+grad1(2)
               grad(3,latm) = grad(3,latm)+grad1(3)
               if(gwt_sw) then
c
c...              If we use the gradients of the weights then we must
c...              also take the gradients of the quadrature points into
c...              account. We do that using the translational invariance
c...              of the energy.
c
                  grad(1,iatom) = grad(1,iatom)-grad1(1)
                  grad(2,iatom) = grad(2,iatom)-grad1(2)
                  grad(3,iatom) = grad(3,iatom)-grad1(3)
               endif
 200        continue
         enddo
      endif
c
c     DEBUG PRINT
c
      if(nchk .ne. n_active_bfn)then
         write(6,*)'********** screening error**********',
     &        nchk,n_active_bfn
         write(6,*)'near atoms',(active_bfn_atms(latm),
     &        latm=1,n_active_bfn_atm)

         write(6,*)'near funcs',(active_bfn_list(nu),nu=1,n_active_bfn)
         do near=1,n_active_bfn_atm
            latm = active_bfn_atms(near)
            bstart=first_bf(latm)
            bend = first_bf(latm+1) - 1
            call atombf(active_bfn_list,n_active_bfn,
     &           bstart,bend,ixf,ixl)
            write(6,*)'bf indices for atom=',latm,':',ixf,ixl,
     &           active_bfn_list(ixf),active_bfn_list(ixl)
         enddo
      endif
c
      return
      end

      subroutine atombf(active_bfn_list,n_active_bfn,
     &     bstart,bend,ixf,ixl)
      implicit none
c in
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer bstart, bend
c out
      integer ixf, ixl
c local
      integer i
      ixf=999999999
      ixl=-1
      do i=1,n_active_bfn
         if(active_bfn_list(i) .ge. bstart)ixf=min(ixf,i)
         if(active_bfn_list(i) .le. bend)ixl=max(ixl,i)
      enddo
      end

      subroutine gather_mu_scr(in,n,active_bfn_list,n_active_bfn,i,out)
      implicit none
c in
      REAL in(*)
      integer n, i
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/mapper)
c out
      REAL out(*)
c local
      integer j, ii, iki

      ii = active_bfn_list(i)
      iki = iky(ii)

      do j = 1,i
         out(j) = in(iki + active_bfn_list(j))
      enddo

      do j = i+1,n_active_bfn
         out(j) = in(iky(active_bfn_list(j)) + ii)
      enddo

      end



C **********************************************************************
c
C xc_forces:  Calculate the derivative of Exc
c 
c  this version does not use screening
c
C **********************************************************************
      subroutine xc_forces(ao_tag,nao,rks_sw,gradcorr_sw,kinetic_sw,
     &                     adens,bdens,gwt_sw,gwt_avail_sw,
     &                     bfn_val,bfng_val,bfn_hess,wt,gwt,
     &                     xc_ept,xc_vpt,xc_dvpt,xc_dtpt,grho,
     &                     scr2,scr3,scr4,
     &                     near_atom_list, num_near_atoms,
     &                     grad,npts,mxp,iatom)
      implicit none
c
c     Description:
c
c     Calculate the first derivative of Exc with respect to the nuclear
c     coordinates without screening.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_dfder)
c
c     Inputs:
c
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_basis)
cINCLUDE(common/dft_module_comm)
      integer mxp, npts
      integer ao_tag, nao
      logical rks_sw,gradcorr_sw,kinetic_sw
      REAL adens(nao*(nao+1)/2), bdens(nao*(nao+1)/2)
      logical gwt_sw,gwt_avail_sw
      REAL bfn_val(mxp,nao)
      REAL bfng_val(mxp,nao,3)
      REAL bfn_hess(mxp,nao,6)
      REAL wt(mxp),gwt(3,mxp,*)
      REAL grho(mxp,2,3)
      REAL xc_ept(mxp),xc_vpt(mxp,2),xc_dvpt(mxp,3),xc_dtpt(mxp,2)
      integer iatom
      integer num_near_atoms, near_atom_list(num_near_atoms)
c
c     Input/Output:
c
      REAL grad(3,natoms)
c
c     Workspace:
c
      REAL scr2(3,nao),scr3(nao),scr4(npts,3)
c
c     Local:
c
      integer alph,beta,spin, nspin
      integer latm,mu,lbas2,lbfn
      integer atmt,bstart,bend,near
      REAL scr1(3),tmp
      REAL grad1(3)
      integer ipt
      character*1 xn, xt
      data xt,xn/'t','n'/
c
C End declarations	
C *********************************************************************
      alph = 1
      beta = 2

      nspin = 1
      if(.not.rks_sw) nspin = 2
c
c     add the contribution from the gradients of the weights
c
      if (gwt_sw.and.gwt_avail_sw) then
         do near=1,num_near_atoms
            latm = near_atom_list(near)
            do ipt = 1, npts
               grad(1,latm) = grad(1,latm)
     +                      + gwt(1,ipt,near)*xc_ept(ipt)
               grad(2,latm) = grad(2,latm)
     +                      + gwt(2,ipt,near)*xc_ept(ipt)
               grad(3,latm) = grad(3,latm)
     +                      + gwt(3,ipt,near)*xc_ept(ipt)
            enddo
         enddo
      endif
c
      do spin = 1, nspin
c
         if (gradcorr_sw) then
            if (nspin.eq.1) then
               do ipt = 1, npts
                  scr4(ipt,1) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
                  scr4(ipt,2) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
                  scr4(ipt,3) = 0.5d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
               enddo
            else if (spin.eq.1) then
               do ipt = 1, npts
                  scr4(ipt,1) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,1))
                  scr4(ipt,2) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,2))
                  scr4(ipt,3) = (2.0d0*xc_dvpt(ipt,igaa)*grho(ipt,1,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,2,3))
               enddo
            else
               do ipt = 1, npts
                  scr4(ipt,1) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,1)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,1))
                  scr4(ipt,2) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,2)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,2))
                  scr4(ipt,3) = (2.0d0*xc_dvpt(ipt,igbb)*grho(ipt,2,3)
     &                        +        xc_dvpt(ipt,igab)*grho(ipt,1,3))
               enddo
            endif
         endif
c
         bend = 0
         do 100 latm=1,ngridcentres
            atmt = atom_tag(ao_tag,latm)
            if(atmt .eq. 0) goto 100
            bstart=bend+1
            bend=bend+Abfn(ao_tag,atmt)
            if (latm.eq.iatom.and.gwt_sw) goto 100
            call aclear_dp(grad1,3,0.0d0)
            do mu=bstart,bend

               if(spin .eq. 1)then
                  call gather_mu(adens,mu,nao,scr3)
               else
                  call gather_mu(bdens,mu,nao,scr3)
               endif

               do ipt = 1,npts

                  scr1(1) = bfng_val(ipt,mu,1)*xc_vpt(ipt,spin)
                  scr1(2) = bfng_val(ipt,mu,2)*xc_vpt(ipt,spin)
                  scr1(3) = bfng_val(ipt,mu,3)*xc_vpt(ipt,spin)
c     
                  if(gradcorr_sw) then
                     scr1(1) = scr1(1) 
     &                    + scr4(ipt,1)*bfn_hess(ipt,mu,1)
     &                    + scr4(ipt,2)*bfn_hess(ipt,mu,4)
     &                    + scr4(ipt,3)*bfn_hess(ipt,mu,5)
                     scr1(2) = scr1(2) 
     &                    + scr4(ipt,1)*bfn_hess(ipt,mu,4)
     &                    + scr4(ipt,2)*bfn_hess(ipt,mu,2)
     &                    + scr4(ipt,3)*bfn_hess(ipt,mu,6)
                     scr1(3) = scr1(3) 
     &                    + scr4(ipt,1)*bfn_hess(ipt,mu,5)
     &                    + scr4(ipt,2)*bfn_hess(ipt,mu,6)
     &                    + scr4(ipt,3)*bfn_hess(ipt,mu,3)
                  endif
c     
                  scr1(1) = scr1(1)*wt(ipt)
                  scr1(2) = scr1(2)*wt(ipt)
                  scr1(3) = scr1(3)*wt(ipt)
c     
                  do lbfn=1,nao
                     scr2(1,lbfn)=scr1(1)*bfn_val(ipt,lbfn)
                     scr2(2,lbfn)=scr1(2)*bfn_val(ipt,lbfn)
                     scr2(3,lbfn)=scr1(3)*bfn_val(ipt,lbfn)
                  enddo
c     
                  if (gradcorr_sw) then
                     do lbas2=1,nao
                        tmp=scr4(ipt,1)*bfng_val(ipt,lbas2,1)+
     &                      scr4(ipt,2)*bfng_val(ipt,lbas2,2)+
     &                      scr4(ipt,3)*bfng_val(ipt,lbas2,3)
                        tmp=tmp*wt(ipt)
                        scr2(1,lbas2) = scr2(1,lbas2)
     &                                + bfng_val(ipt,mu,1)*tmp
                        scr2(2,lbas2) = scr2(2,lbas2)
     &                                + bfng_val(ipt,mu,2)*tmp
                        scr2(3,lbas2) = scr2(3,lbas2)
     &                                + bfng_val(ipt,mu,3)*tmp
                     enddo
                  endif
c
                  if (kinetic_sw) then
                     do lbas2=1,nao
                        scr2(1,lbas2) = scr2(1,lbas2) + 
     &                      wt(ipt)*xc_dtpt(ipt,spin)*
     &                      (bfng_val(ipt,lbas2,1)*bfn_hess(ipt,mu,1)+
     &                       bfng_val(ipt,lbas2,2)*bfn_hess(ipt,mu,4)+
     &                       bfng_val(ipt,lbas2,3)*bfn_hess(ipt,mu,5))
                        scr2(2,lbas2) = scr2(2,lbas2) + 
     &                      wt(ipt)*xc_dtpt(ipt,spin)*
     &                      (bfng_val(ipt,lbas2,1)*bfn_hess(ipt,mu,4)+
     &                       bfng_val(ipt,lbas2,2)*bfn_hess(ipt,mu,2)+
     &                       bfng_val(ipt,lbas2,3)*bfn_hess(ipt,mu,6))
                        scr2(3,lbas2) = scr2(3,lbas2) + 
     &                      wt(ipt)*xc_dtpt(ipt,spin)*
     &                      (bfng_val(ipt,lbas2,1)*bfn_hess(ipt,mu,5)+
     &                       bfng_val(ipt,lbas2,2)*bfn_hess(ipt,mu,6)+
     &                       bfng_val(ipt,lbas2,3)*bfn_hess(ipt,mu,3))
                     enddo
                  endif
c     
                  call dgemv(xn,3,nao,-2.0d0
     +                 ,scr2,3,scr3,1,1.0d0,grad1,1)
c     
               enddo               ! loop over quadrature points
            enddo
c     
            grad(1,latm) = grad(1,latm)+grad1(1)
            grad(2,latm) = grad(2,latm)+grad1(2)
            grad(3,latm) = grad(3,latm)+grad1(3)
            if(gwt_sw) then
c     
c...           If we use the gradients of the weights then we must also 
c...           take the gradients of the quadrature points into account.
c...           We do that using the translational invariance of the 
c...           energy.
c     
               grad(1,iatom) = grad(1,iatom)-grad1(1)
               grad(2,iatom) = grad(2,iatom)-grad1(2)
               grad(3,iatom) = grad(3,iatom)-grad1(3)
            endif
 100     continue                  ! loop over atoms
      enddo                        ! loop over spin
      return
      end

      subroutine gather_mu(dens,mu,nu_num,scr3)
      implicit none
      REAL dens(*)
      integer mu,nu_num
      REAL scr3(*)
      integer lbas,maxa,mina,numu_pos
      do lbas=1,nu_num
        maxa=max(mu,lbas)
        mina=min(mu,lbas)
        numu_pos=(((maxa*maxa)-maxa)/2)+mina
        scr3(lbas)=dens(numu_pos)
      enddo
      return
      end

      integer function atom_num_by_tag(atm_tag)
      implicit none
      integer atm_tag
c
c     Returns the atomic number Z of an atom specified by its tag.
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_basis_api)
INCLUDE(../m4/common/errcodes)
      integer i, ao_tag
      logical opg_root
      call caserr('I thought this function was obsolete...')
      atom_num_by_tag = -1
      ao_tag = 1
      if(atm_tag .eq. 0)then
c
c null tag - can't get atom number
c
         return
      endif
      do i = 1, ngridcentres
         if (BL_get_atom_type(ao_tag,i).eq.atm_tag) then
            atom_num_by_tag = ian(i)
            return
         endif
      enddo
      if(opg_root())write(6,*)'ATOM_NUM_BY_TAG atm_tag = ',atm_tag

      call gamerr('ATOM_NUM_BY_TAG atm_tag not found',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_SYNC, ERR_NO_SYS)

      return
      end


      integer function atom_num_by_grid(grid_tag)
      implicit none
      integer grid_tag
c
c     Returns the atomic number Z of an atom specified by its grid 
c     identifier.
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(../m4/common/errcodes)
      integer i
      logical opg_root
      atom_num_by_grid = -1
      if(grid_tag .eq. 0)then
c
c null tag - can't get atom number
c
         return
      else if (grid_tag.lt.0.or.grid_tag.gt.ngtypes) then
         write(6,*)'Grid_tag        = ',grid_tag
         write(6,*)'Number of grids = ',ngtypes
         call caserr('Grid_tag out of range')
      endif
      do i = 1, natoms
         if (gtype_num(i).eq.grid_tag) then
            atom_num_by_grid = ian(i)
            return
         endif
      enddo
      if(opg_root())write(6,*)'ATOM_NUM_BY_GRID grid_tag = ',grid_tag

      call gamerr('ATOM_NUM_BY_GRID grid_tag not found',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_SYNC, ERR_NO_SYS)

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine screen_weight(npts,npack,ra2_val,ra2_comp,wt,wt2,
     &     tolerance,natoms,ngridcentres,
     &     n_active_wgh_atm,active_wgh_atms,mxp)
      implicit none
INCLUDE(common/dft_parameters)
c
c...  Input
c
      integer natoms
      integer ngridcentres
      integer mxp
      integer n_active_wgh_atm
      integer active_wgh_atms(n_active_wgh_atm) ! the active atoms
      REAL tolerance
c
c...  Input/Output
c
      integer npts   ! the number of points in this batch
      integer npack  ! the number of skipped points overall
      REAL ra2_val(mxp,natoms,2)
      REAL ra2_comp(mxp,natoms,3)
      REAL wt(mxp)   ! On input : the partition function
                     ! On output: the non-zero weights
      REAL wt2(mxp)  ! copy of the original Lebedev or Legendre weights
c
c...  Local variables
c
      integer ipt, npnew, i, iatm
      logical omove
      REAL test, weight
c
      npnew = 0
      omove = .false.
      do ipt = 1, npts
c
c...     Absolute value needed because Lebedev weights may be
c...     negative
c
         test   = wt(ipt)
         weight = test*wt2(ipt)
         test   = max(test,abs(weight))
         if (test .gt. tolerance) then
            npnew = npnew + 1
            wt(npnew) = weight
            if (omove) then
               wt2(npnew) = wt2(ipt)
               do i = 1, n_active_wgh_atm
                  iatm = active_wgh_atms(i)
                  ra2_val(npnew,iatm,1)  = ra2_val(ipt,iatm,1)
                  ra2_val(npnew,iatm,2)  = ra2_val(ipt,iatm,2)
                  ra2_comp(npnew,iatm,1) = ra2_comp(ipt,iatm,1)
                  ra2_comp(npnew,iatm,2) = ra2_comp(ipt,iatm,2)
                  ra2_comp(npnew,iatm,3) = ra2_comp(ipt,iatm,3)
               enddo
            endif
         else
            npack = npack + 1
            omove = .true.
         endif
      enddo
c
      npts = npnew
      end
c
c-----------------------------------------------------------------------
c
      subroutine npoints_by_accuracy(grad,natyp,igrid,nrad,
     &           ntheta,nphi,nomega)
      implicit none
INCLUDE(common/dft_parameters)
      REAL grad
      integer natyp,igrid
      integer nrad(natyp), ntheta(natyp), nphi(natyp), nomega(natyp)
      integer mtheta, mphi, mang
INCLUDE(common/dft_module_comm)
c
c     Computes the number of grid points needed to obtain appropriate
c     integration precision in a situation where the SCF proces has
c     converged to within a orbital gradient <grad>. The returned number
c     of grid points never exceeds the limits set in the input
c
      integer i, j
      REAL qlimit
      parameter (qlimit = 1.0d-4)
      REAL scale
c
      if (grad.eq.0.0d0.or..not.conv_prune_sw) then
         do i = 1, ngtypes
            nrad(i)   = radpt_num(i,igrid)
            ntheta(i) = 0
            nphi(i)   = 0
            nomega(i) = 0
            do j = 1, radzones_num(i,igrid)
               ntheta(i) = max(ntheta(i),thetpt_radzn_num(j,i,igrid))
               nphi(i)   = max(nphi(i)  ,phipt_radzn_num(j,i,igrid))
               nomega(i) = max(nomega(i),angpt_radzn_num(j,i,igrid))
            enddo
         enddo
      else
c
c...     Compute the scale factor for the number of grid points. 
c...     The number of grid points should scale linearly with the 
c...     accuracy in the energy. 
c
c...     Assuming:
c...     - that the accuracy in the energy is quadratic with the size of
c...       the SCF gradient, 
c...     - that the quadrature accuracy is linear in the molecular grid 
c...       size, 
c...     - that the molecular grid size is linear in the atomic grid 
c...       sizes,
c...     - that the atomic grid size is cubic in the number of points
c...       per polar coordinate
c...     - and that the same scaling parameter is used for the number of
c...       grid point for all atomic coordinates,
c...     then the scale factor should be related to the "accuracy" by
c...     a power 2/3. The factor <qlimit> guarantees that the maximum
c...     grid size is used if the SCF gradient falls below <qlimit>.
c
         scale  = abs(qlimit/grad)**(2.0d0/3.0d0)
c        scale  = (1.0d0+2.2d-3)-scale
         do i = 1, ngtypes
            nrad(i)   = min(max(int(scale*radpt_num(i,igrid)),20),
     +                          radpt_num(i,igrid))
            mtheta    = 0
            mphi      = 0
            mang      = 0
            do j = 1, radzones_num(i,igrid)
               mtheta = max(mtheta,thetpt_radzn_num(j,i,igrid))
               mphi   = max(mphi  ,phipt_radzn_num(j,i,igrid))
               mang   = max(mang  ,angpt_radzn_num(j,i,igrid))
            enddo
            ntheta(i) = min(max(int(scale*mtheta),10),mtheta)
            nphi(i)   = min(max(int(scale*mphi),20),mphi)
            nomega(i) = min(max(int(scale*scale*mang),110),mang)
         enddo
c        scale  = -log10(min(1.0d0,max(1.0d-10,abs(accuracy))))/4.5d0
c        nrad   = 20+int(scale*max(0,radpt_num-20))
c        ntheta = 10+int(scale*max(0,thetpt_num-10))
c        nphi   = 20+int(scale*max(0,phipt_num-20))
c        nomega = 200+int(scale*max(0,angupt_num-200))
c
c        nrad   = min(radpt_num,nrad)
c        ntheta = min(thetpt_num,ntheta)
c        nphi   = min(phipt_num,nphi)
c        nomega = min(angupt_num,nomega)
      endif
c
      end


      integer function row_by_atomnum(Z)
      implicit none
c
c     This function returns the row-number in the periodic table of 
c     the elements given the nuclear charge Z.
c
      integer Z
      if (Z.le.0) then
         call caserr("row_by_atomnum: Atomic number must be at least 1")
      else if (Z.le.2) then
         row_by_atomnum = 1
      else if (Z.le.10) then
         row_by_atomnum = 2
      else if (Z.le.18) then
         row_by_atomnum = 3
      else if (Z.le.36) then
         row_by_atomnum = 4
      else if (Z.le.54) then
         row_by_atomnum = 5
      else if (Z.le.86) then
         row_by_atomnum = 6
      else if (Z.le.118) then
         row_by_atomnum = 7
      else
         call caserr("row_by_atomnum:Atomic number must be at most 118")
      endif
      return
      end
c
      subroutine auto_prune(igrid,nradpt_num,inmtyp,xc_e,rad)
      implicit none
c
c     This subroutine selects the angular grids as a function of
c     the radius on basis of the total energy accumulated at
c     in a particular shell. 
c
c     The shell with the largest amount of energy is taken to
c     require 100% of the angular grid size. Other shells get
c     less angular grid points based on the fraction of the
c     energy they bring.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
      integer maxleb, maxpt
      parameter (maxleb=32)
      integer leba(maxleb)
      data leba/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,
     1          590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,
     1          3890,4334,4802,5294,5810/
      parameter(maxpt=5810)
      integer ao_tag
      parameter(ao_tag=1)
c
c     In/Out variables:
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
c
c     In variables:
c
      integer igrid
      integer nradpt_num(ngtypes)
      integer inmtyp(0:max_gtype)
      REAL xc_e(*)
      REAL rad(ngtypes,*)
c
c     Local variables:
c
      integer ig, it, i, ileb, ileb2, iang, iend, irad
      integer maxang1, mxl, minleb
      integer imaxe
      REAL maxe
      REAL rend
c
c     Functions:
c
      integer BL_maxang_on_atom
      external BL_maxang_on_atom
c
c     Code:
c
      iauto_cnt = iauto_cnt+1
      if (iauto_cnt.ne.1.and.
     +    iauto_cnt.ne.2.and.
     +    iauto_cnt.ne.7) return
c
c     Find current maximum number of angular grid points
c
      maxang1 = 0
      do ig = 1, ngtypes
         do i = 1, radzones_num(ig,igrid)
            maxang1 = max(maxang1,angpt_radzn_num(i,ig,igrid))
         enddo
      enddo
c
c     Adapt the energy per shell to the function we need
c     The multiplication with (1+r)r**2 is purely emperical
c
      do ig = 1, ngtypes
         do irad = 2, nradpt_num(ig)
            xc_e(inmtyp(ig-1)+irad) =
     +         xc_e(inmtyp(ig-1)+irad)*
     +         rad(ig,irad)**2
     +         *(1.0d0+rad(ig,irad))
         enddo
      enddo
c
c     Adjust the grid parameters 
c
      do ig = 1, ngtypes
         if (ang_prune_scheme(ig,igrid) .eq. AP_AUTO) then
c
c           find current maximum energy in a shell
c
            maxe  = 0.0d0
            imaxe = 0
            do i = inmtyp(ig-1)+1, inmtyp(ig)
               if (abs(xc_e(i)).gt.maxe) then
                  maxe  = abs(xc_e(i))
                  imaxe = i
               endif
            enddo
c
c           find maximum angular momentum
c
            mxl = 0
            do i = 1, ngridcentres
               if (gtype_num(i).eq.ig) then
                  mxl = max(mxl,BL_maxang_on_atom(ao_tag,i))
               endif
            enddo
            if (mxl.le.2) then
               minleb  = 1
            else if (mxl.le.3) then
               minleb  = 2
            else if (mxl.le.4) then
               minleb  = 3
            else if (mxl.le.5) then
               minleb  = 4
            else
               minleb  = 5
            endif
            radzones_num(ig,igrid) = 1
            rend   = 0.0d0
            iend   = 0
            iang   = int(sqrt(abs(xc_e(inmtyp(ig-1)+iend+1))/maxe)
     &                   *maxang1+0.5d0)
            do i = maxleb, minleb, -1
               if (leba(i).ge.iang) ileb = leba(i)
            enddo
            ileb2 = ileb
 10         if (iend.lt.nradpt_num(ig).and.ileb2.eq.ileb) then
               iend = iend+1
               iang = int(sqrt(abs(xc_e(inmtyp(ig-1)+iend+1))/maxe)
     &                    *maxang1+0.5d0)
               do i = maxleb, minleb, -1
                  if (leba(i).ge.iang) ileb2 = leba(i)
               enddo
               goto 10
            endif
            if (iend.lt.nradpt_num(ig)) then
c
c              this radial zone is NOT the outer most radial zone
c
               rend = 0.5d0*(rad(ig,iend)+rad(ig,iend+1))
               if (radzones_num(ig,igrid).ge.maxradzn) then
                  call caserr('auto_prune: too many radial zones!!!')
               endif
               bnd_radzn(radzones_num(ig,igrid),ig,igrid) = rend
               angpt_radzn_num(radzones_num(ig,igrid),ig,igrid) = ileb
               radzones_num(ig,igrid) = radzones_num(ig,igrid) + 1
               ileb = ileb2
               goto 10
            else
c
c              this radial zone is the outer most radial zone
c              therefore the boundary is implied to be infinity
c
               angpt_radzn_num(radzones_num(ig,igrid),ig,igrid) = ileb
            endif
         endif
      enddo
      return
      end
c     
      subroutine ver_dft_xc(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/xc.m,v $
     +     "/
      data revision /
     +     "$Revision: 6206 $"
     +      /
      data date /
     +     "$Date: 2010-10-31 18:05:53 +0100 (Sun, 31 Oct 2010) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
